#!/usr/bin/env bash
# sync_hpc_data.sh - Synchronize data and libraries between Princeton HPC and local Thinkpad
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Source local environment defaults if available
if [[ -f "$PROJECT_ROOT/env.sh" ]]; then
    # shellcheck disable=SC1091
    source "$PROJECT_ROOT/env.sh"
fi

HPC_HOST="${HPC_HOST:-stellar-vis1.princeton.edu}"
HPC_USER="${HPC_USER:-mbolot}"
SSH_TARGET="${HPC_USER}@${HPC_HOST}"

HPC_PROCESSED_DIR="${HPC_PROCESSED_DIR:-/scratch/gpfs/mbolot/results/GLOBALFV3/}"
HPC_RAW_DIR="${HPC_RAW_DIR:-/scratch/cimes/GLOBALFV3/}"
HPC_LOCAL_DIR="${HPC_LOCAL_DIR:-~mbolot/local/}"

LOCAL_PROCESSED_DIR="${GLOBALFV3_DATA_ROOT:-/media/max/data/climate_work/processed}"
LOCAL_RAW_DIR="${GLOBALFV3_RAW_ROOT:-/media/max/data/climate_work/raw}"
LOCAL_LIBS_DIR="${EXT_LIB_DIR:-$HOME/local}"

PROFILE="analysis"
DRY_RUN=""
SIM_FILTER="all"
DATE_OVERRIDE=""
VARS_PROFILE="all"
EXTRA_RSYNC_ARGS=()

# Exact 17 coarse 1440x720 variables consumed by compute_work Fortran binary.
# Ingesting only these files saves ~99% bandwidth and local storage by avoiding
# multi-terabyte native cubed-sphere C3072 tiles (*.tile*.nc) and restart files.
RAW_WORKFLOW_FILES=(
    "DZ_C3072_1440x720.fre.nc"
    "temp_coarse_C3072_1440x720.fre.nc"
    "ptend_coarse_C3072_1440x720.fre.nc"
    "sphum_coarse_C3072_1440x720.fre.nc"
    "liq_wat_coarse_C3072_1440x720.fre.nc"
    "rainwat_coarse_C3072_1440x720.fre.nc"
    "ice_wat_coarse_C3072_1440x720.fre.nc"
    "snowwat_coarse_C3072_1440x720.fre.nc"
    "graupel_coarse_C3072_1440x720.fre.nc"
    "omT_coarse_C3072_1440x720.fre.nc"
    "omqv_coarse_C3072_1440x720.fre.nc"
    "omql_coarse_C3072_1440x720.fre.nc"
    "omqr_coarse_C3072_1440x720.fre.nc"
    "omqi_coarse_C3072_1440x720.fre.nc"
    "omqs_coarse_C3072_1440x720.fre.nc"
    "omqg_coarse_C3072_1440x720.fre.nc"
    "PRATEsfc_coarse_C3072_1440x720.fre.nc"
)

RAW_PRATE_FILES=(
    "PRATEsfc_coarse_C3072_1440x720.fre.nc"
)

usage() {
    cat << EOF
Usage: $(basename "$0") [options]

Options:
  --profile <name>         Sync profile:
                             analysis    : Sync complete processed 360x180 datasets, histograms, and thresholds (Default)
                             sample-raw  : Sync only the 17 workflow variables for benchmark sample dates from /scratch/cimes
                             vendor-libs : Sync ~/local libraries (libinterp, libnc, libhist) from HPC
                             all         : Run analysis, sample-raw, and vendor-libs
  --sim <sim>              Simulation for sample-raw: control, warming, or all (Default: all)
  --date <YYYYMMDDHH>      Sync a single specific date instead of reading list files
  --vars <mode>            Variables to sync for sample-raw:
                             all         : All 17 coarse 1440x720 variables required by compute_work (Default)
                             prate       : Only PRATEsfc coarse field (for compute_prate_thresholds)
                             <f1,f2,...> : Custom comma-separated list of filenames
  --dry-run                Show files that would be transferred without downloading
  --host <hpc_host>        HPC Hostname (Default: $HPC_HOST)
  --user <hpc_user>        HPC Username (Default: $HPC_USER)
  --bwlimit <KB/s>         Bandwidth limit in KB/s (e.g. --bwlimit 10000 for ~10MB/s)
  -h, --help               Show this help message

Examples:
  $(basename "$0") --profile analysis --dry-run
  $(basename "$0") --profile sample-raw --dry-run
  $(basename "$0") --profile sample-raw --sim control
  $(basename "$0") --profile sample-raw --sim warming --date 2020080500
  $(basename "$0") --profile sample-raw --vars prate
  $(basename "$0") --profile vendor-libs
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --sim)
            SIM_FILTER="$2"
            shift 2
            ;;
        --date)
            DATE_OVERRIDE="$2"
            shift 2
            ;;
        --vars)
            VARS_PROFILE="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN="--dry-run"
            shift
            ;;
        --host)
            HPC_HOST="$2"
            SSH_TARGET="${HPC_USER}@${HPC_HOST}"
            shift 2
            ;;
        --user)
            HPC_USER="$2"
            SSH_TARGET="${HPC_USER}@${HPC_HOST}"
            shift 2
            ;;
        --bwlimit)
            EXTRA_RSYNC_ARGS+=("--bwlimit=$2")
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown argument: $1" >&2
            usage
            ;;
    esac
done

RSYNC_BASE_FLAGS=(-avP --human-readable)
if [[ -n "$DRY_RUN" ]]; then
    RSYNC_BASE_FLAGS+=("$DRY_RUN")
    echo "=== RUNNING IN DRY-RUN MODE (No files will be transferred) ==="
fi

echo "Connecting to: $SSH_TARGET"

sync_vendor_libs() {
    echo -e "\n--- [Profile: vendor-libs] Syncing ~/local libraries from HPC ---"
    local dest_dir="$HOME/local"
    mkdir -p "$dest_dir" || true

    rsync "${RSYNC_BASE_FLAGS[@]}" "${EXTRA_RSYNC_ARGS[@]}" \
        "${SSH_TARGET}:local/" "${dest_dir}/"

    echo "Vendor libraries synced to: $dest_dir"
}

sync_analysis() {
    echo -e "\n--- [Profile: analysis] Syncing complete processed dataset (360x180 + histograms + thresholds) ---"
    mkdir -p "$LOCAL_PROCESSED_DIR" || true

    # Synchronize processed data directory while excluding large 1440x720 intermediate directories
    rsync "${RSYNC_BASE_FLAGS[@]}" "${EXTRA_RSYNC_ARGS[@]}" \
        --exclude="*1440x720*/" \
        --include="*/" \
        --include="*.nc" \
        --include="*.txt" \
        --exclude="*" \
        "${SSH_TARGET}:${HPC_PROCESSED_DIR}" "${LOCAL_PROCESSED_DIR}/"

    # Also sync thresholds directory if it exists under output/
    if [[ -d "$PROJECT_ROOT/output/thresholds" ]] || [[ -n "$DRY_RUN" ]]; then
        echo "Syncing threshold files..."
        mkdir -p "$PROJECT_ROOT/output/thresholds" || true
        rsync "${RSYNC_BASE_FLAGS[@]}" "${EXTRA_RSYNC_ARGS[@]}" \
            "${SSH_TARGET}:/home/${HPC_USER}/PROJECTS/MECHANICAL_WORK/output/thresholds/" \
            "$PROJECT_ROOT/output/thresholds/" || true
    fi

    echo "Processed analysis dataset synced to: $LOCAL_PROCESSED_DIR"
}

sync_sample_raw() {
    echo -e "\n--- [Profile: sample-raw] Syncing workflow raw variables from /scratch/cimes ---"

    local control_list="$PROJECT_ROOT/launcher/list/list_sample_control.txt"
    local warming_list="$PROJECT_ROOT/launcher/list/list_sample_warming.txt"

    local control_remote_base="/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history"
    local control_local_base="$LOCAL_RAW_DIR/20191020.00Z.C3072.L79x2_pire/history"
    local warming_remote_base="/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    local warming_local_base="$LOCAL_RAW_DIR/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"

    # Determine files to transfer
    local -a target_files=()
    case "$VARS_PROFILE" in
        all)
            target_files=("${RAW_WORKFLOW_FILES[@]}")
            ;;
        prate|precip)
            target_files=("${RAW_PRATE_FILES[@]}")
            ;;
        *)
            IFS=',' read -ra target_files <<< "$VARS_PROFILE"
            ;;
    esac

    echo "Selected ${#target_files[@]} variable file(s) per date (excluding multi-terabyte native grids):"
    for f in "${target_files[@]}"; do
        echo "  - $f"
    done

    # Build rsync include rules specifically for these files, then exclude everything else
    local -a rsync_var_filters=()
    for f in "${target_files[@]}"; do
        rsync_var_filters+=(--include="$f")
    done
    rsync_var_filters+=(--exclude="*")

    sync_date_dir() {
        local sim_name="$1"
        local date="$2"
        local remote_base="$3"
        local local_base="$4"

        echo "  -> [$sim_name] Syncing date: $date"
        mkdir -p "${local_base}/${date}" 2>/dev/null || true
        rsync "${RSYNC_BASE_FLAGS[@]}" "${EXTRA_RSYNC_ARGS[@]}" \
            "${rsync_var_filters[@]}" \
            "${SSH_TARGET}:${remote_base}/${date}/" \
            "${local_base}/${date}/"
    }

    if [[ -n "$DATE_OVERRIDE" ]]; then
        echo "Single date override specified: $DATE_OVERRIDE (sim: $SIM_FILTER)"
        if [[ "$SIM_FILTER" == "control" || "$SIM_FILTER" == "all" ]]; then
            sync_date_dir "Control" "$DATE_OVERRIDE" "$control_remote_base" "$control_local_base"
        fi
        if [[ "$SIM_FILTER" == "warming" || "$SIM_FILTER" == "all" ]]; then
            sync_date_dir "Warming" "$DATE_OVERRIDE" "$warming_remote_base" "$warming_local_base"
        fi
    else
        # 1. Sync Control benchmark dates
        if [[ "$SIM_FILTER" == "control" || "$SIM_FILTER" == "all" ]]; then
            if [[ -f "$control_list" ]]; then
                echo "Syncing Control sample dates from: $control_list"
                while IFS= read -r date_line || [[ -n "$date_line" ]]; do
                    local date
                    date=$(echo "$date_line" | tr -d '[:space:]')
                    [[ -z "$date" || "$date" == \#* ]] && continue
                    sync_date_dir "Control" "$date" "$control_remote_base" "$control_local_base"
                done < "$control_list"
            else
                echo "Warning: Control list file not found: $control_list" >&2
            fi
        fi

        # 2. Sync Warming benchmark dates
        if [[ "$SIM_FILTER" == "warming" || "$SIM_FILTER" == "all" ]]; then
            if [[ -f "$warming_list" ]]; then
                echo "Syncing Warming sample dates from: $warming_list"
                while IFS= read -r date_line || [[ -n "$date_line" ]]; do
                    local date
                    date=$(echo "$date_line" | tr -d '[:space:]')
                    [[ -z "$date" || "$date" == \#* ]] && continue
                    sync_date_dir "Warming" "$date" "$warming_remote_base" "$warming_local_base"
                done < "$warming_list"
            else
                echo "Warning: Warming list file not found: $warming_list" >&2
            fi
        fi
    fi

    echo "Sample raw variables synced to: $LOCAL_RAW_DIR"
}

case "$PROFILE" in
    analysis)
        sync_analysis
        ;;
    sample-raw)
        sync_sample_raw
        ;;
    vendor-libs)
        sync_vendor_libs
        ;;
    all)
        sync_vendor_libs
        sync_analysis
        sync_sample_raw
        ;;
    *)
        echo "Error: Unknown profile '$PROFILE'. Use: analysis, sample-raw, vendor-libs, or all." >&2
        exit 1
        ;;
esac

echo -e "\n[$(date +%T)] Synchronization completed successfully."
