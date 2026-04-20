#!/bin/bash
#SBATCH -A geoclim
#SBATCH --time=22:00:00
#SBATCH --qos=cimes-short
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mbolot@princeton.edu

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}
SCRIPT_PATH="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"

NTASKS=${NTASKS:-1}
CPUS_PER_TASK=${CPUS_PER_TASK:-2}
# remapcon + in-memory operator chaining can hold multiple full-resolution slices.
# A conservative default is 8G per CPU; increase if jobs OOM on your partition.
MEM_PER_CPU=${MEM_PER_CPU:-8G}
LOG_DIR=${LOG_DIR:-$PROJECT_ROOT/logs}

# Keep CDO internal OpenMP threads aligned with Slurm CPU placement.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-$CPUS_PER_TASK}}
CDO_P_THREADS=${CDO_P_THREADS:-$OMP_NUM_THREADS}

SOURCE_ROOT=${SOURCE_ROOT:-/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history}
LIST_FILE=${LIST_FILE:-$PROJECT_ROOT/launcher/list/list_control_11520x5760.txt}
TARGET_ROOT=${TARGET_ROOT:-/scratch/gpfs/mbolot/data/20191020.00Z.C3072.L79x2_pire.ea_to_0.25}
# group1: fields 1-12, group2: fields 13-17, all: submit both groups in one array.
RUN_GROUP=${RUN_GROUP:-all}

GRID_FILE=${GRID_FILE:-$PROJECT_ROOT/output/grids/ea_25km_grid.txt}
WEIGHTS_DIR=${WEIGHTS_DIR:-$PROJECT_ROOT/output/grids}
WEIGHTS_NATIVE_TO_EA=${WEIGHTS_NATIVE_TO_EA:-$WEIGHTS_DIR/weights_native_to_ea.nc}
WEIGHTS_EA_TO_025=${WEIGHTS_EA_TO_025:-$WEIGHTS_DIR/weights_ea_to_0.25.nc}

# Native high-resolution defaults (11520x5760).
NATIVE_DZ_FILE=${NATIVE_DZ_FILE:-delz_C3072_11520x5760.fre.nc}
NATIVE_OMEGA_FILE=${NATIVE_OMEGA_FILE:-omega_C3072_11520x5760.fre.nc}
NATIVE_QV_FILE=${NATIVE_QV_FILE:-specific_humidity_C3072_11520x5760.fre.nc}
NATIVE_QW_FILE=${NATIVE_QW_FILE:-clw_C3072_11520x5760.fre.nc}
NATIVE_QR_FILE=${NATIVE_QR_FILE:-rainwat_C3072_11520x5760.fre.nc}
NATIVE_QI_FILE=${NATIVE_QI_FILE:-cli_C3072_11520x5760.fre.nc}
NATIVE_QS_FILE=${NATIVE_QS_FILE:-snowwat_C3072_11520x5760.fre.nc}
NATIVE_QG_FILE=${NATIVE_QG_FILE:-graupel_C3072_11520x5760.fre.nc}
NATIVE_TEMP_FILE=${NATIVE_TEMP_FILE:-ta_C3072_11520x5760.fre.nc}
NATIVE_PR_FILE=${NATIVE_PR_FILE:-pr_C3072_11520x5760.fre.nc}
CDO_MODULE_LOAD=${CDO_MODULE_LOAD:-intel/2021.1.2 hdf5/intel-2021.1/1.10.6 netcdf/intel-2021.1/hdf5-1.10.6/4.7.4 cdo/netcdf-4.7.4/hdf5-1.10.6/2.0.1}

require_cmd() {
    local cmd=$1
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Error: required command '$cmd' not found in PATH" >&2
        exit 1
    fi
}

ensure_cdo_environment() {
    if command -v cdo >/dev/null 2>&1; then
        return
    fi

    if ! command -v module >/dev/null 2>&1; then
        for module_init in /usr/share/Modules/init/bash /etc/profile.d/modules.sh; do
            if [[ -f "$module_init" ]]; then
                # shellcheck source=/dev/null
                source "$module_init"
                break
            fi
        done
    fi

    if command -v module >/dev/null 2>&1; then
        # Load the exact compiler/netcdf/cdo stack required on stellar-vis1.
        module load $CDO_MODULE_LOAD
    fi

    if ! command -v cdo >/dev/null 2>&1; then
        echo "Error: required command 'cdo' not found in PATH. The command should be module load $CDO_MODULE_LOAD" >&2
        exit 1
    fi
}

trimmed_dates() {
    sed -E '/^[[:space:]]*($|#)/d; s/^[[:space:]]+//; s/[[:space:]]+$//' "$LIST_FILE"
}

detect_var_name() {
    local file=$1
    cdo -s -P "$CDO_P_THREADS" showname "$file" | awk '{print $1}'
}

cdo_run() {
    # Centralize CDO parallel flags so every operator uses the same thread count.
    cdo -L -P "$CDO_P_THREADS" "$@"
}

remap_mean_to_fortran() {
    local native_file=$1
    local reference_file=$2
    local out_file=$3

    local src_var dst_var
    src_var=$(detect_var_name "$native_file")
    dst_var=$(detect_var_name "$reference_file")

    if [[ -z "$src_var" || -z "$dst_var" ]]; then
        echo "Error: failed to detect variable names for $out_file" >&2
        exit 1
    fi

    if [[ "$src_var" != "$dst_var" ]]; then
        # Preserve Fortran-facing variable names from legacy coarse files.
        cdo_run chname,"$src_var","$dst_var" -remap,r1440x720,"$WEIGHTS_EA_TO_025" -remap,"$GRID_FILE","$WEIGHTS_NATIVE_TO_EA" "$native_file" "$out_file"
    else
        cdo_run remap,r1440x720,"$WEIGHTS_EA_TO_025" -remap,"$GRID_FILE","$WEIGHTS_NATIVE_TO_EA" "$native_file" "$out_file"
    fi
}

remap_covariance_to_fortran() {
    local native_w=$1
    local native_x=$2
    local reference_cov=$3
    local out_file=$4

    local tmp_ea_wx tmp_ea_wx_prod src_var dst_var
    tmp_ea_wx=$(mktemp -p "$TMP_WORKDIR" ea_wx_XXXXXX.nc)
    tmp_ea_wx_prod=$(mktemp -p "$TMP_WORKDIR" ea_wx_prod_XXXXXX.nc)

    # First term: E[w*x] on EA grid.
    cdo_run remap,"$GRID_FILE","$WEIGHTS_NATIVE_TO_EA" -mul "$native_w" "$native_x" "$tmp_ea_wx"
    # Second term: E[w]E[x] on EA grid.
    cdo_run mul -remap,"$GRID_FILE","$WEIGHTS_NATIVE_TO_EA" "$native_w" -remap,"$GRID_FILE","$WEIGHTS_NATIVE_TO_EA" "$native_x" "$tmp_ea_wx_prod"

    src_var=$(detect_var_name "$native_w")
    dst_var=$(detect_var_name "$reference_cov")

    if [[ -z "$src_var" || -z "$dst_var" ]]; then
        echo "Error: failed to detect variable names for covariance output $out_file" >&2
        exit 1
    fi

    if [[ "$src_var" != "$dst_var" ]]; then
        # Covariance = E[w*x] - E[w]E[x], then remap to Fortran grid.
        cdo_run chname,"$src_var","$dst_var" -remap,r1440x720,"$WEIGHTS_EA_TO_025" -sub "$tmp_ea_wx" "$tmp_ea_wx_prod" "$out_file"
    else
        cdo_run remap,r1440x720,"$WEIGHTS_EA_TO_025" -sub "$tmp_ea_wx" "$tmp_ea_wx_prod" "$out_file"
    fi

    rm -f "$tmp_ea_wx" "$tmp_ea_wx_prod"
}

maybe_copy_or_fail() {
    local src_file=$1
    local dst_file=$2
    if [[ ! -f "$src_file" ]]; then
        echo "Error: required fallback file missing: $src_file" >&2
        exit 1
    fi
    cp -f "$src_file" "$dst_file"
}

produce_mean_field() {
    local source_dir=$1
    local target_dir=$2
    local out_name=$3
    local native_name=$4

    local native_file="$source_dir/$native_name"
    local reference_file="$source_dir/$out_name"
    local out_file="$target_dir/$out_name"

    if [[ -f "$native_file" ]]; then
        remap_mean_to_fortran "$native_file" "$reference_file" "$out_file"
    else
        # If a native variable is unavailable, keep pipeline moving with legacy coarse input.
        echo "Warning: native file missing for $out_name ($native_name). Copying existing coarse file." >&2
        maybe_copy_or_fail "$reference_file" "$out_file"
    fi
}

produce_cov_field() {
    local source_dir=$1
    local target_dir=$2
    local out_name=$3
    local native_x_name=$4

    local native_w_file="$source_dir/$NATIVE_OMEGA_FILE"
    local native_x_file="$source_dir/$native_x_name"
    local reference_cov="$source_dir/$out_name"
    local out_file="$target_dir/$out_name"

    if [[ -f "$native_w_file" && -f "$native_x_file" ]]; then
        remap_covariance_to_fortran "$native_w_file" "$native_x_file" "$reference_cov" "$out_file"
    else
        # Covariance fallback mirrors current operational behavior when native pairs are missing.
        echo "Warning: native files missing for $out_name ($NATIVE_OMEGA_FILE, $native_x_name). Copying existing coarse file." >&2
        maybe_copy_or_fail "$reference_cov" "$out_file"
    fi
}

produce_dz_field() {
    local source_dir=$1
    local target_dir=$2

    local native_file="$source_dir/$NATIVE_DZ_FILE"
    local reference_file="$source_dir/DZ_C3072_1440x720.fre.nc"
    local out_file="$target_dir/DZ_C3072_1440x720.fre.nc"
    local src_var dst_var

    if [[ ! -f "$native_file" ]]; then
        echo "Warning: native DZ file missing ($NATIVE_DZ_FILE). Copying existing coarse DZ file." >&2
        maybe_copy_or_fail "$reference_file" "$out_file"
        return
    fi

    src_var=$(detect_var_name "$native_file")
    dst_var=$(detect_var_name "$reference_file")
    if [[ -z "$src_var" || -z "$dst_var" ]]; then
        echo "Error: failed to detect variable names for DZ output" >&2
        exit 1
    fi

    if [[ "$src_var" != "$dst_var" ]]; then
        # Native delz is opposite-sign relative to legacy DZ; apply -1 to match Fortran input convention.
        cdo_run chname,"$src_var","$dst_var" -mulc,-1 -remap,r1440x720,"$WEIGHTS_EA_TO_025" -remap,"$GRID_FILE","$WEIGHTS_NATIVE_TO_EA" "$native_file" "$out_file"
    else
        cdo_run mulc,-1 -remap,r1440x720,"$WEIGHTS_EA_TO_025" -remap,"$GRID_FILE","$WEIGHTS_NATIVE_TO_EA" "$native_file" "$out_file"
    fi
}

ensure_weights_from_template() {
    local template_native=$1

    mkdir -p "$WEIGHTS_DIR"

    if [[ ! -f "$WEIGHTS_NATIVE_TO_EA" ]]; then
        echo "Generating weights: native -> EA" >&2
        cdo_run gencon,"$GRID_FILE" "$template_native" "$WEIGHTS_NATIVE_TO_EA"
    fi

    if [[ ! -f "$WEIGHTS_EA_TO_025" ]]; then
        echo "Generating weights: EA -> r1440x720" >&2
        # Build a lightweight EA-grid template (2D only) to avoid remapping a huge native file.
        local tmp_r_template tmp_ea_template
        tmp_r_template=$(mktemp -p "$TMP_WORKDIR" r_template_XXXXXX.nc)
        tmp_ea_template=$(mktemp -p "$TMP_WORKDIR" ea_template_XXXXXX.nc)
        cdo_run -f nc topo,r1440x720 "$tmp_r_template"
        cdo_run setgrid,"$GRID_FILE" "$tmp_r_template" "$tmp_ea_template"
        cdo_run gencon,r1440x720 "$tmp_ea_template" "$WEIGHTS_EA_TO_025"
        rm -f "$tmp_r_template"
        rm -f "$tmp_ea_template"
    fi
}

find_template_native_file() {
    local date candidate
    mapfile -t template_dates < <(trimmed_dates)
    for date in "${template_dates[@]}"; do
        [[ -z "$date" ]] && continue
        candidate="$SOURCE_ROOT/$date/$NATIVE_OMEGA_FILE"
        if [[ -f "$candidate" ]]; then
            echo "$candidate"
            return 0
        fi
    done
    return 1
}

run_group1_fields() {
    local date=$1
    local source_dir=$2
    local target_dir=$3

    echo "[$date][group1] Generating DZ_C3072_1440x720.fre.nc"
    produce_dz_field "$source_dir" "$target_dir"

    echo "[$date][group1] Generating ptend_coarse_C3072_1440x720.fre.nc"
    produce_mean_field "$source_dir" "$target_dir" "ptend_coarse_C3072_1440x720.fre.nc" "$NATIVE_OMEGA_FILE"
    echo "[$date][group1] Generating sphum_coarse_C3072_1440x720.fre.nc"
    produce_mean_field "$source_dir" "$target_dir" "sphum_coarse_C3072_1440x720.fre.nc" "$NATIVE_QV_FILE"
    echo "[$date][group1] Generating liq_wat_coarse_C3072_1440x720.fre.nc"
    produce_mean_field "$source_dir" "$target_dir" "liq_wat_coarse_C3072_1440x720.fre.nc" "$NATIVE_QW_FILE"
    echo "[$date][group1] Generating rainwat_coarse_C3072_1440x720.fre.nc"
    produce_mean_field "$source_dir" "$target_dir" "rainwat_coarse_C3072_1440x720.fre.nc" "$NATIVE_QR_FILE"
    echo "[$date][group1] Generating ice_wat_coarse_C3072_1440x720.fre.nc"
    produce_mean_field "$source_dir" "$target_dir" "ice_wat_coarse_C3072_1440x720.fre.nc" "$NATIVE_QI_FILE"
    echo "[$date][group1] Generating snowwat_coarse_C3072_1440x720.fre.nc"
    produce_mean_field "$source_dir" "$target_dir" "snowwat_coarse_C3072_1440x720.fre.nc" "$NATIVE_QS_FILE"
    echo "[$date][group1] Generating graupel_coarse_C3072_1440x720.fre.nc"
    produce_mean_field "$source_dir" "$target_dir" "graupel_coarse_C3072_1440x720.fre.nc" "$NATIVE_QG_FILE"
    echo "[$date][group1] Generating temp_coarse_C3072_1440x720.fre.nc"
    produce_mean_field "$source_dir" "$target_dir" "temp_coarse_C3072_1440x720.fre.nc" "$NATIVE_TEMP_FILE"
    echo "[$date][group1] Generating PRATEsfc_coarse_C3072_1440x720.fre.nc"
    produce_mean_field "$source_dir" "$target_dir" "PRATEsfc_coarse_C3072_1440x720.fre.nc" "$NATIVE_PR_FILE"

    echo "[$date][group1] Generating omT_coarse_C3072_1440x720.fre.nc"
    produce_cov_field "$source_dir" "$target_dir" "omT_coarse_C3072_1440x720.fre.nc" "$NATIVE_TEMP_FILE"
    echo "[$date][group1] Generating omqv_coarse_C3072_1440x720.fre.nc"
    produce_cov_field "$source_dir" "$target_dir" "omqv_coarse_C3072_1440x720.fre.nc" "$NATIVE_QV_FILE"
}

run_group2_fields() {
    local date=$1
    local source_dir=$2
    local target_dir=$3

    echo "[$date][group2] Generating omql_coarse_C3072_1440x720.fre.nc"
    produce_cov_field "$source_dir" "$target_dir" "omql_coarse_C3072_1440x720.fre.nc" "$NATIVE_QW_FILE"
    echo "[$date][group2] Generating omqr_coarse_C3072_1440x720.fre.nc"
    produce_cov_field "$source_dir" "$target_dir" "omqr_coarse_C3072_1440x720.fre.nc" "$NATIVE_QR_FILE"
    echo "[$date][group2] Generating omqi_coarse_C3072_1440x720.fre.nc"
    produce_cov_field "$source_dir" "$target_dir" "omqi_coarse_C3072_1440x720.fre.nc" "$NATIVE_QI_FILE"
    echo "[$date][group2] Generating omqs_coarse_C3072_1440x720.fre.nc"
    produce_cov_field "$source_dir" "$target_dir" "omqs_coarse_C3072_1440x720.fre.nc" "$NATIVE_QS_FILE"
    echo "[$date][group2] Generating omqg_coarse_C3072_1440x720.fre.nc"
    produce_cov_field "$source_dir" "$target_dir" "omqg_coarse_C3072_1440x720.fre.nc" "$NATIVE_QG_FILE"
}

process_single_date() {
    local date=$1
    local group=$2
    local source_dir="$SOURCE_ROOT/$date"
    local target_dir="$TARGET_ROOT/$date"

    if [[ ! -d "$source_dir" ]]; then
        echo "Error: source directory not found: $source_dir" >&2
        exit 1
    fi

    mkdir -p "$target_dir"

    # Each task processes only one group to keep per-task runtime under queue limits.
    case "$group" in
        group1)
            run_group1_fields "$date" "$source_dir" "$target_dir"
            ;;
        group2)
            run_group2_fields "$date" "$source_dir" "$target_dir"
            ;;
        *)
            echo "Error: unsupported RUN_GROUP='$group'. Use 'group1', 'group2', or 'all'." >&2
            exit 1
            ;;
    esac
}

ensure_cdo_environment
require_cmd cdo
require_cmd sed

if [[ ! -f "$LIST_FILE" ]]; then
    echo "Error: list file not found: $LIST_FILE" >&2
    exit 1
fi

if [[ ! -f "$GRID_FILE" ]]; then
    echo "Error: EA grid file not found: $GRID_FILE" >&2
    exit 1
fi

case "$RUN_GROUP" in
    group1|group2|all)
        ;;
    *)
        echo "Error: unsupported RUN_GROUP='$RUN_GROUP'. Use 'group1', 'group2', or 'all'." >&2
        exit 1
        ;;
esac

TMP_WORKDIR=${TMP_WORKDIR:-/tmp/ea_to_0.25_${USER}_$$}
mkdir -p "$TMP_WORKDIR"
trap 'rm -rf "$TMP_WORKDIR"' EXIT

# Launch mode: submit array from list.
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    mkdir -p "$LOG_DIR"

    mapfile -t launch_dates < <(trimmed_dates)
    n_dates=${#launch_dates[@]}
    if [[ "$n_dates" -lt 1 ]]; then
        echo "Error: date list is empty: $LIST_FILE" >&2
        exit 1
    fi

    case "$RUN_GROUP" in
        group1|group2)
            # One task per date when only one group is requested.
            n_tasks=$n_dates
            ;;
        all)
            # all = two tasks per date: first pass group1, second pass group2.
            n_tasks=$(( 2 * n_dates ))
            ;;
    esac

    template_native=$(find_template_native_file || true)
    if [[ -z "${template_native:-}" ]]; then
        echo "Error: unable to find a native template file '$NATIVE_OMEGA_FILE' from list dates." >&2
        exit 1
    fi

    ensure_weights_from_template "$template_native"

    echo "Submitting ea_to_0.25 job array with $n_tasks tasks (RUN_GROUP=$RUN_GROUP, dates=$n_dates)"
    echo "Resources per array task: ntasks=$NTASKS, cpus-per-task=$CPUS_PER_TASK, mem-per-cpu=$MEM_PER_CPU"
    echo "CDO threading: -P $CDO_P_THREADS (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

    sbatch --array=1-"$n_tasks" --ntasks="$NTASKS" --cpus-per-task="$CPUS_PER_TASK" --mem-per-cpu="$MEM_PER_CPU" \
           --output="$LOG_DIR/ea_to_0.25_%A_%a.out" --error="$LOG_DIR/ea_to_0.25_%A_%a.err" \
           --export=ALL,PROJECT_ROOT="$PROJECT_ROOT",SOURCE_ROOT="$SOURCE_ROOT",LIST_FILE="$LIST_FILE",TARGET_ROOT="$TARGET_ROOT",GRID_FILE="$GRID_FILE",WEIGHTS_DIR="$WEIGHTS_DIR",WEIGHTS_NATIVE_TO_EA="$WEIGHTS_NATIVE_TO_EA",WEIGHTS_EA_TO_025="$WEIGHTS_EA_TO_025",NTASKS="$NTASKS",CPUS_PER_TASK="$CPUS_PER_TASK",MEM_PER_CPU="$MEM_PER_CPU",LOG_DIR="$LOG_DIR",RUN_GROUP="$RUN_GROUP" "$SCRIPT_PATH"
    exit $?
fi

mapfile -t dates < <(trimmed_dates)
idx=${SLURM_ARRAY_TASK_ID}
n_dates=${#dates[@]}
if (( n_dates < 1 )); then
    echo "Error: no dates found in task mode." >&2
    exit 1
fi

group_to_run=$RUN_GROUP
date_idx=$idx

if [[ "$RUN_GROUP" == "all" ]]; then
    # In all-mode, first half of array indices map to group1 and second half to group2.
    total_tasks=$(( 2 * n_dates ))
    if (( idx < 1 || idx > total_tasks )); then
        echo "Error: array index out of bounds: $idx (n=$total_tasks)" >&2
        exit 1
    fi
    if (( idx <= n_dates )); then
        group_to_run=group1
        date_idx=$idx
    else
        group_to_run=group2
        date_idx=$(( idx - n_dates ))
    fi
else
    if (( idx < 1 || idx > n_dates )); then
        echo "Error: array index out of bounds: $idx (n=$n_dates)" >&2
        exit 1
    fi
fi

if (( date_idx < 1 || date_idx > n_dates )); then
    echo "Error: resolved date index out of bounds: $date_idx (n=$n_dates)" >&2
    exit 1
fi

if [[ ! -f "$WEIGHTS_NATIVE_TO_EA" || ! -f "$WEIGHTS_EA_TO_025" ]]; then
    template_native="$SOURCE_ROOT/${dates[0]}/$NATIVE_OMEGA_FILE"
    if [[ ! -f "$template_native" ]]; then
        echo "Error: weights are missing and template native file not found: $template_native" >&2
        exit 1
    fi
    ensure_weights_from_template "$template_native"
fi

date=${dates[date_idx-1]}

if [[ "$RUN_GROUP" == "all" ]]; then
    echo "Processing date $date (array task $idx/$((2 * n_dates)), group=$group_to_run)"
else
    echo "Processing date $date (array task $idx/$n_dates, group=$group_to_run)"
fi

process_single_date "$date" "$group_to_run"
echo "Finished date $date (group=$group_to_run)"
