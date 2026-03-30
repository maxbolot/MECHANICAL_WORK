#!/usr/bin/env bash
set -euo pipefail

# Serial runner for systems without Slurm.
# One date (yyyymmdd) per line in list file.
LIST_FILE=${1:-list.txt}

# Hardcoded OpenMP setting for serial execution.
export OMP_NUM_THREADS=1

# Paths
SOURCE_ROOT=/scratch/gpfs/mbolot/data/20191020.00Z.C3072.L79x2_pire
TARGET_DIR_COMPUTE=/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720
TARGET_DIR_HISTOGRAMS=/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms
LOG_DIR=./logs

# Prefer local/bin build output, fallback to legacy bin path.
COMPUTE_WORK_BIN=./local/bin/compute_work
if [[ ! -x "$COMPUTE_WORK_BIN" ]]; then
    COMPUTE_WORK_BIN=./bin/compute_work
fi

if [[ ! -f "$LIST_FILE" ]]; then
    echo "Error: list file not found: $LIST_FILE" >&2
    exit 1
fi

if [[ ! -x "$COMPUTE_WORK_BIN" ]]; then
    echo "Error: compute_work binary not found or not executable." >&2
    echo "Tried: ./local/bin/compute_work and ./bin/compute_work" >&2
    exit 1
fi

mkdir -p "$TARGET_DIR_COMPUTE" "$TARGET_DIR_HISTOGRAMS" "$LOG_DIR"

RUN_LOG="$LOG_DIR/compute_work_fork_$(date +%Y%m%d_%H%M%S).log"
echo "[$(date +%F\ %T)] Starting serial run with list: $LIST_FILE" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Using binary: $COMPUTE_WORK_BIN" | tee -a "$RUN_LOG"

while IFS= read -r raw_line; do
    date=$(echo "$raw_line" | tr -d '[:space:]')

    # Skip empty lines and comments.
    [[ -z "$date" ]] && continue
    [[ "$date" == \#* ]] && continue

    source_dir="$SOURCE_ROOT/$date"
    if [[ ! -d "$source_dir" ]]; then
        echo "Warning: source directory missing for date $date, skipping" >&2
        continue
    fi

    config_file="config_${date}.nml"
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
    path_work_out   = '$TARGET_DIR_COMPUTE/work_$date.nc',
    path_hist_out   = '$TARGET_DIR_HISTOGRAMS/hist_$date.nc',
/
EOF

    date_log="$LOG_DIR/compute_work_${date}.log"
    echo "[$(date +%F\ %T)] Running compute_work for date $date (log: $date_log)" | tee -a "$RUN_LOG"
    "$COMPUTE_WORK_BIN" "$config_file" > "$date_log" 2>&1
    echo "[$(date +%F\ %T)] Finished date $date" | tee -a "$RUN_LOG"
    rm -f "$config_file"
done < "$LIST_FILE"

echo "[$(date +%F\ %T)] All dates processed serially." | tee -a "$RUN_LOG"
