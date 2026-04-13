#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

# Minimal one-date test launcher for compute_work_async.
# Source root and date are intentionally hardcoded for quick validation.
SOURCE_ROOT=/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history
DATE=2020010300

TARGET_DIR_COMPUTE=${TARGET_DIR_COMPUTE:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_test}
TARGET_DIR_HISTOGRAMS=${TARGET_DIR_HISTOGRAMS:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_test}
LOG_DIR=${LOG_DIR:-$PROJECT_ROOT/logs}
NAMELIST_DIR=$PROJECT_ROOT/output/namelists

module purge || true
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-3}

COMPUTE_WORK_BIN=${COMPUTE_WORK_BIN:-$PROJECT_ROOT/local/bin/compute_work_async}
if [[ ! -x "$COMPUTE_WORK_BIN" && "$COMPUTE_WORK_BIN" == "$PROJECT_ROOT/local/bin/compute_work_async" ]]; then
    COMPUTE_WORK_BIN=$PROJECT_ROOT/bin/compute_work_async
fi

if [[ ! -d "$SOURCE_ROOT" ]]; then
    echo "Error: source root directory not found: $SOURCE_ROOT" >&2
    exit 1
fi

if [[ ! -x "$COMPUTE_WORK_BIN" ]]; then
    echo "Error: compute_work_async binary not found or not executable." >&2
    echo "Tried: $PROJECT_ROOT/local/bin/compute_work_async and $PROJECT_ROOT/bin/compute_work_async" >&2
    exit 1
fi

source_dir="$SOURCE_ROOT/$DATE"
if [[ ! -d "$source_dir" ]]; then
    echo "Error: source directory not found: $source_dir" >&2
    exit 1
fi

mkdir -p "$TARGET_DIR_COMPUTE" "$TARGET_DIR_HISTOGRAMS" "$LOG_DIR" "$NAMELIST_DIR"

RUN_LOG="$LOG_DIR/compute_work_async_test_${DATE}_$(date +%Y%m%d_%H%M%S).log"
CONFIG_FILE="$NAMELIST_DIR/config_test_${DATE}_$$.nml"

cat > "$CONFIG_FILE" << EOF
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
    path_work_out   = '$TARGET_DIR_COMPUTE/work_${DATE}.nc',
    path_hist_out   = '$TARGET_DIR_HISTOGRAMS/hist_${DATE}.nc',
/
EOF

echo "[$(date +%F\ %T)] Starting compute_work_async test run" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Date: $DATE" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Source root: $SOURCE_ROOT" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Binary: $COMPUTE_WORK_BIN" | tee -a "$RUN_LOG"

"$COMPUTE_WORK_BIN" "$CONFIG_FILE" 2>&1 | tee -a "$RUN_LOG"

echo "[$(date +%F\ %T)] Done. Outputs written to $TARGET_DIR_COMPUTE and $TARGET_DIR_HISTOGRAMS" | tee -a "$RUN_LOG"