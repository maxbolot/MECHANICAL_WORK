#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

# Serial runner (no Slurm) for thresholded async work/lift computation.

SIMULATION=${SIMULATION:-control}

case "$SIMULATION" in
    control)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp
        LIST_FILE_PART1=${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_control_part1.txt}
        LIST_FILE_PART2=${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_control_part2.txt}
        TARGET_DIR_COMPUTE=${TARGET_DIR_COMPUTE:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_prate_thresholded}
        THRESHOLD_FILE=${THRESHOLD_FILE:-$PROJECT_ROOT/output/thresholds/thresholds_control.txt}
        ;;
    warming)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        LIST_FILE_PART1=${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt}
        LIST_FILE_PART2=${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt}
        TARGET_DIR_COMPUTE=${TARGET_DIR_COMPUTE:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_PLUS_4K_CO2_1270ppmv_prate_thresholded}
        THRESHOLD_FILE=${THRESHOLD_FILE:-$PROJECT_ROOT/output/thresholds/thresholds_warming.txt}
        ;;
    *)
        echo "Error: unsupported SIMULATION='$SIMULATION'. Use 'control' or 'warming'." >&2
        exit 1
        ;;
esac

LOG_DIR=${LOG_DIR:-$PROJECT_ROOT/logs}
NAMELIST_DIR=$PROJECT_ROOT/output/namelists

module purge || true
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-3}

COMPUTE_BIN=${COMPUTE_BIN:-$PROJECT_ROOT/local/bin/compute_work_async_prate_threshold}
if [[ ! -x "$COMPUTE_BIN" && "$COMPUTE_BIN" == "$PROJECT_ROOT/local/bin/compute_work_async_prate_threshold" ]]; then
    COMPUTE_BIN=$PROJECT_ROOT/bin/compute_work_async_prate_threshold
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

if [[ ! -f "$THRESHOLD_FILE" ]]; then
    echo "Error: threshold file not found: $THRESHOLD_FILE" >&2
    exit 1
fi

if [[ ! -x "$COMPUTE_BIN" ]]; then
    echo "Error: thresholded compute binary not found or not executable." >&2
    echo "Tried: $PROJECT_ROOT/local/bin/compute_work_async_prate_threshold and $PROJECT_ROOT/bin/compute_work_async_prate_threshold" >&2
    exit 1
fi

mkdir -p "$TARGET_DIR_COMPUTE" "$LOG_DIR" "$NAMELIST_DIR"

RUN_LOG="$LOG_DIR/compute_work_async_prate_threshold_noslurm_$(date +%Y%m%d_%H%M%S).log"
echo "[$(date +%F\ %T)] Starting serial run for simulation: $SIMULATION" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] List file part1: $LIST_FILE_PART1" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] List file part2: $LIST_FILE_PART2" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Using binary: $COMPUTE_BIN" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Threshold file: $THRESHOLD_FILE" | tee -a "$RUN_LOG"

process_list() {
    local list_file="$1"
    local source_root="$2"

    while IFS= read -r raw_line; do
        date=$(echo "$raw_line" | tr -d '[:space:]')

        [[ -z "$date" ]] && continue
        [[ "$date" == \#* ]] && continue

        source_dir="$source_root/$date"
        if [[ ! -d "$source_dir" ]]; then
            echo "Warning: source directory missing for date $date in $source_root, skipping" | tee -a "$RUN_LOG"
            continue
        fi

        config_file="$NAMELIST_DIR/config_prate_work_${date}.nml"
        cat > "$config_file" << EOF
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
    path_work_out   = '$TARGET_DIR_COMPUTE/work_prate_threshold_${date}.nc',
    path_thresholds = '$THRESHOLD_FILE',
/
EOF

        date_log="$LOG_DIR/compute_work_async_prate_threshold_${date}.log"
        echo "[$(date +%F\ %T)] Running thresholded compute for date $date (log: $date_log)" | tee -a "$RUN_LOG"
        "$COMPUTE_BIN" "$config_file" > "$date_log" 2>&1
        rm -f "$config_file"
        echo "[$(date +%F\ %T)] Finished date $date" | tee -a "$RUN_LOG"
    done < "$list_file"
}

process_list "$LIST_FILE_PART1" "$SOURCE_ROOT_PART1"
process_list "$LIST_FILE_PART2" "$SOURCE_ROOT_PART2"

echo "[$(date +%F\ %T)] All dates processed serially." | tee -a "$RUN_LOG"
