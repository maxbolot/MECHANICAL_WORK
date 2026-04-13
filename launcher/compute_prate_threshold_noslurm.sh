#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

# Serial launcher (no Slurm) for compute_prate_thresholds.
# The Fortran program reads two history roots and two date lists.

# Simulation selector:
#   control (default)
#   warming
# Manual override example:
#   SIMULATION=warming ./compute_prate_threshold_noslurm.sh
SIMULATION=${SIMULATION:-control}

case "$SIMULATION" in
    control)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp
        LIST_FILE_PART1=$PROJECT_ROOT/launcher/list/list_control_part1.txt
        LIST_FILE_PART2=$PROJECT_ROOT/launcher/list/list_control_part2.txt
        OUTPUT_FILE=$PROJECT_ROOT/output/thresholds/thresholds_control.txt
        ;;
    warming)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        LIST_FILE_PART1=$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt
        LIST_FILE_PART2=$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt
        OUTPUT_FILE=$PROJECT_ROOT/output/thresholds/thresholds_warming.txt
        ;;
    *)
        echo "Error: unsupported SIMULATION='$SIMULATION'. Use 'control' or 'warming'." >&2
        exit 1
        ;;
esac

# Fixed output location is selected by SIMULATION (no separate output override).
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
NAMELIST_DIR=$PROJECT_ROOT/output/namelists

LOG_DIR=${LOG_DIR:-$PROJECT_ROOT/logs}
mkdir -p "$LOG_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$NAMELIST_DIR"

module purge || true
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

# Prefer local/bin build output, fallback to legacy bin path.
COMPUTE_PRATE_THRESH_BIN=${COMPUTE_PRATE_THRESH_BIN:-$PROJECT_ROOT/local/bin/compute_prate_thresholds}
if [[ ! -x "$COMPUTE_PRATE_THRESH_BIN" && "$COMPUTE_PRATE_THRESH_BIN" == "$PROJECT_ROOT/local/bin/compute_prate_thresholds" ]]; then
    COMPUTE_PRATE_THRESH_BIN=$PROJECT_ROOT/bin/compute_prate_thresholds
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

if [[ ! -x "$COMPUTE_PRATE_THRESH_BIN" ]]; then
    echo "Error: compute_prate_thresholds binary not found or not executable." >&2
    echo "Tried: $PROJECT_ROOT/local/bin/compute_prate_thresholds and $PROJECT_ROOT/bin/compute_prate_thresholds" >&2
    exit 1
fi

CONFIG_FILE="$NAMELIST_DIR/config_prate_threshold_$$.nml"
RUN_LOG="$LOG_DIR/compute_prate_threshold_noslurm_$(date +%Y%m%d_%H%M%S).log"

cat > "$CONFIG_FILE" << EOF
&config
    history_root_part1   = '$SOURCE_ROOT_PART1',
    history_root_part2   = '$SOURCE_ROOT_PART2',
    date_list_file_part1 = '$LIST_FILE_PART1',
    date_list_file_part2 = '$LIST_FILE_PART2',
    output_file          = '$OUTPUT_FILE',
/
EOF

echo "[$(date +%F\ %T)] Starting compute_prate_thresholds" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Simulation: $SIMULATION" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Binary: $COMPUTE_PRATE_THRESH_BIN" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Source root part1: $SOURCE_ROOT_PART1" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Source root part2: $SOURCE_ROOT_PART2" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] List file part1: $LIST_FILE_PART1" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] List file part2: $LIST_FILE_PART2" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Output file: $OUTPUT_FILE" | tee -a "$RUN_LOG"

"$COMPUTE_PRATE_THRESH_BIN" "$CONFIG_FILE" 2>&1 | tee -a "$RUN_LOG"
rm -f "$CONFIG_FILE"

echo "[$(date +%F\ %T)] Done. Thresholds written to: $OUTPUT_FILE" | tee -a "$RUN_LOG"
