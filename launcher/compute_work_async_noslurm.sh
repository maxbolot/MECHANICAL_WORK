#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

# Serial runner for systems without Slurm.
# Dates are read from two list files, each mapped to its own source root.

SIMULATION=${SIMULATION:-control}

case "$SIMULATION" in
    control)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp
        LIST_FILE_PART1=${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_control_part1.txt}
        LIST_FILE_PART2=${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_control_part2.txt}
        TARGET_DIR_COMPUTE=${TARGET_DIR_COMPUTE:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720}
        TARGET_DIR_HISTOGRAMS=${TARGET_DIR_HISTOGRAMS:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms}
        ;;
    warming)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        LIST_FILE_PART1=${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt}
        LIST_FILE_PART2=${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt}
        TARGET_DIR_COMPUTE=${TARGET_DIR_COMPUTE:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_PLUS_4K_CO2_1270ppmv}
        TARGET_DIR_HISTOGRAMS=${TARGET_DIR_HISTOGRAMS:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_PLUS_4K_CO2_1270ppmv}
        ;;
    *)
        echo "Error: unsupported SIMULATION='$SIMULATION'. Use 'control' or 'warming'." >&2
        exit 1
        ;;
esac

module purge || true
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

# Hardcoded OpenMP setting for serial execution.
export OMP_NUM_THREADS=2

LOG_DIR=$PROJECT_ROOT/logs
NAMELIST_DIR=$PROJECT_ROOT/output/namelists

# Prefer local/bin build output, fallback to legacy bin path.
COMPUTE_WORK_BIN=$PROJECT_ROOT/local/bin/compute_work_async
if [[ ! -x "$COMPUTE_WORK_BIN" ]]; then
    COMPUTE_WORK_BIN=$PROJECT_ROOT/bin/compute_work_async
fi

if [[ ! -f "$LIST_FILE_PART1" ]]; then
    echo "Error: list file not found: $LIST_FILE_PART1" >&2
    exit 1
fi

if [[ ! -f "$LIST_FILE_PART2" ]]; then
    echo "Error: list file not found: $LIST_FILE_PART2" >&2
    exit 1
fi

if [[ ! -d "$SOURCE_ROOT_PART1" ]]; then
    echo "Error: source root directory not found: $SOURCE_ROOT_PART1" >&2
    exit 1
fi

if [[ ! -d "$SOURCE_ROOT_PART2" ]]; then
    echo "Error: source root directory not found: $SOURCE_ROOT_PART2" >&2
    exit 1
fi

if [[ ! -x "$COMPUTE_WORK_BIN" ]]; then
    echo "Error: compute_work binary not found or not executable." >&2
    echo "Tried: ./local/bin/compute_work and ./bin/compute_work" >&2
    exit 1
fi

mkdir -p "$TARGET_DIR_COMPUTE" "$TARGET_DIR_HISTOGRAMS" "$LOG_DIR" "$NAMELIST_DIR"

RUN_LOG="$LOG_DIR/compute_work_noslurm_$(date +%Y%m%d_%H%M%S).log"
echo "[$(date +%F\ %T)] Starting serial run for simulation: $SIMULATION" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] List file part1: $LIST_FILE_PART1" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] List file part2: $LIST_FILE_PART2" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Source root part1: $SOURCE_ROOT_PART1" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Source root part2: $SOURCE_ROOT_PART2" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Using binary: $COMPUTE_WORK_BIN" | tee -a "$RUN_LOG"

FAILED_DATES=()

process_list() {
    local list_file="$1"
    local source_root="$2"
    local date
    local source_dir
    local config_file
    local date_log
    local rc

    while IFS= read -r raw_line; do
        date=$(echo "$raw_line" | tr -d '[:space:]')

        [[ -z "$date" ]] && continue
        [[ "$date" == \#* ]] && continue

        source_dir="$source_root/$date"
        if [[ ! -d "$source_dir" ]]; then
            echo "Warning: source directory missing for date $date in $source_root, skipping" | tee -a "$RUN_LOG"
            FAILED_DATES+=("$date (missing source directory)")
            continue
        fi

        config_file="$NAMELIST_DIR/config_${date}.nml"
        if ! cat > "$config_file" << EOF
&config
    path_dz         = '$source_dir/DZ_C3072_1440x720.fre.nc',
    path_temp       = '$source_dir/temp_coarse_C3072_1440x720.fre.nc',
    path_omega      = '$source_dir/ptend_coarse_C3072_1440x720.fre.nc',
    path_qv         = '$source_dir/sphum_coarse_C3072_1440x720.fre.nc',
    path_qw         = '$source_dir/liq_wat_coarse_C3072_1440x720.fre.nc',
    path_qr         = '$source_dir/rainwat_coarse_C3072_1440x720.fre.nc',
    path_qi         = '$source_dir/ice_wat_coarse_C3072_1440x720.fre.nc',
    path_qs         = '$source_dir/snowwat_coarse_C3072_1440x720.fre.nc',
    path_qg         = '$source_dir/graupel_coarse_C3072_1440x720.fre.nc',
    path_omt        = '$source_dir/omT_coarse_C3072_1440x720.fre.nc',
    path_omqv       = '$source_dir/omqv_coarse_C3072_1440x720.fre.nc',
    path_omqw       = '$source_dir/omql_coarse_C3072_1440x720.fre.nc',
    path_omqr       = '$source_dir/omqr_coarse_C3072_1440x720.fre.nc',
    path_omqi       = '$source_dir/omqi_coarse_C3072_1440x720.fre.nc',
    path_omqs       = '$source_dir/omqs_coarse_C3072_1440x720.fre.nc',
    path_omqg       = '$source_dir/omqg_coarse_C3072_1440x720.fre.nc',
    path_pr         = '$source_dir/PRATEsfc_coarse_C3072_1440x720.fre.nc',
    path_work_out   = '$TARGET_DIR_COMPUTE/work_$date.nc',
    path_hist_out   = '$TARGET_DIR_HISTOGRAMS/hist_$date.nc',
/
EOF
        then
            echo "[$(date +%F\ %T)] Failed to create config for date $date, skipping" | tee -a "$RUN_LOG"
            FAILED_DATES+=("$date (failed to create namelist)")
            continue
        fi

        date_log="$LOG_DIR/compute_work_${date}.log"
        echo "[$(date +%F\ %T)] Running compute_work for date $date (log: $date_log)" | tee -a "$RUN_LOG"
        if "$COMPUTE_WORK_BIN" "$config_file" > "$date_log" 2>&1; then
            echo "[$(date +%F\ %T)] Finished date $date" | tee -a "$RUN_LOG"
        else
            rc=$?
            echo "[$(date +%F\ %T)] Failed date $date (exit code: $rc). See $date_log" | tee -a "$RUN_LOG"
            FAILED_DATES+=("$date (compute_work exit $rc)")
        fi

        if ! rm -f "$config_file"; then
            echo "[$(date +%F\ %T)] Warning: failed to remove temporary config $config_file" | tee -a "$RUN_LOG"
        fi
    done < "$list_file"
}

process_list "$LIST_FILE_PART1" "$SOURCE_ROOT_PART1"
process_list "$LIST_FILE_PART2" "$SOURCE_ROOT_PART2"

if (( ${#FAILED_DATES[@]} > 0 )); then
    echo "[$(date +%F\ %T)] All dates processed serially with ${#FAILED_DATES[@]} failure(s)." | tee -a "$RUN_LOG"
    echo "[$(date +%F\ %T)] Failed dates summary:" | tee -a "$RUN_LOG"
    for failed_entry in "${FAILED_DATES[@]}"; do
        echo "  - $failed_entry" | tee -a "$RUN_LOG"
    done
    exit 1
fi

echo "[$(date +%F\ %T)] All dates processed serially with no failures." | tee -a "$RUN_LOG"
