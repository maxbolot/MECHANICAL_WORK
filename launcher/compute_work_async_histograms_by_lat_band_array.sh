#!/bin/bash
#SBATCH -A geoclim
#SBATCH --time=04:00:00
#SBATCH --qos=cimes-short
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mbolot@princeton.edu

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}
SCRIPT_PATH="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"

NTASKS=${NTASKS:-1}
CPUS_PER_TASK=${CPUS_PER_TASK:-2}
MEM_PER_CPU=${MEM_PER_CPU:-16G}
RUN_ID=${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}
LOG_DIR=${LOG_DIR:-$PROJECT_ROOT/logs/slurm}
MANIFEST_DIR=${MANIFEST_DIR:-$PROJECT_ROOT/logs/manifests}
MANIFEST_TOOL=${MANIFEST_TOOL:-$PROJECT_ROOT/python/workflow_manifest.py}
PERIOD_TASK_DIR=${PERIOD_TASK_DIR:-}
PERIOD_KEYS_FILE=${PERIOD_KEYS_FILE:-}

SUBMIT_ENV_FILE=${SUBMIT_ENV_FILE:-}
if [[ -n "$SUBMIT_ENV_FILE" && -f "$SUBMIT_ENV_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$SUBMIT_ENV_FILE"
fi

SIMULATION=${SIMULATION:-control_monthly_by_lat_band}

# Select input roots, date lists, and output target by experiment.
case "$SIMULATION" in
    control_monthly_by_lat_band)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp
        LIST_FILE_PART1=${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_control_part1.txt}
        LIST_FILE_PART2=${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_control_part2.txt}
        TARGET_DIR_HIST_BASE=${TARGET_DIR_HIST_BASE:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_monthly_by_lat_band}
        ;;
    warming_monthly_by_lat_band)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        LIST_FILE_PART1=${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt}
        LIST_FILE_PART2=${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt}
        TARGET_DIR_HIST_BASE=${TARGET_DIR_HIST_BASE:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_monthly_by_lat_band_PLUS_4K_CO2_1270ppmv}
        ;;
    *)
        echo "Error: unsupported SIMULATION='$SIMULATION'. Use 'control_monthly_by_lat_band' or 'warming_monthly_by_lat_band'." >&2
        exit 1
        ;;
esac

TARGET_DIR_HIST=${TARGET_DIR_HIST:-$TARGET_DIR_HIST_BASE/$RUN_ID}

NLAT_BANDS=${NLAT_BANDS:-18}
USE_CUSTOM_LAT_BAND_BOUNDS=${USE_CUSTOM_LAT_BAND_BOUNDS:-false}
LAT_BAND_BOUNDS=${LAT_BAND_BOUNDS:-}
AGGREGATION_MODE=${AGGREGATION_MODE:-monthly}

case "$AGGREGATION_MODE" in
    monthly)
        OUTPUT_FILE_PREFIX=${OUTPUT_FILE_PREFIX:-hist_monthly_by_lat_band}
        ;;
    daily)
        OUTPUT_FILE_PREFIX=${OUTPUT_FILE_PREFIX:-hist_daily_by_lat_band}
        ;;
    *)
        echo "Error: unsupported AGGREGATION_MODE='$AGGREGATION_MODE'. Use 'daily' or 'monthly'." >&2
        exit 1
        ;;
esac

period_key_from_date() {
    local date_key="$1"
    case "$AGGREGATION_MODE" in
        monthly)
            echo "${date_key:0:6}"
            ;;
        daily)
            echo "${date_key:0:8}"
            ;;
        *)
            echo ""
            ;;
    esac
}

# Submission mode: compute the combined task count and submit the array job.
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    mkdir -p "$LOG_DIR" "$MANIFEST_DIR"

    if [[ ! -f "$LIST_FILE_PART1" ]]; then
        echo "Error: list file not found: $LIST_FILE_PART1" >&2
        exit 1
    fi

    if [[ ! -f "$LIST_FILE_PART2" ]]; then
        echo "Error: list file not found: $LIST_FILE_PART2" >&2
        exit 1
    fi

    # Persist submission settings to a shell-readable file so array tasks see
    # the exact values even when Slurm export parsing is unfriendly to commas.
    submit_stamp=$(date +%Y%m%d_%H%M%S)
    PERIOD_TASK_DIR="$LOG_DIR/compute_work_async_histograms_by_lat_band_tasks_${submit_stamp}_$$"
    PERIOD_KEYS_FILE="$PERIOD_TASK_DIR/period_keys.txt"
    combined_file="$PERIOD_TASK_DIR/dates_with_roots.txt"
    sorted_file="$PERIOD_TASK_DIR/dates_with_roots_sorted.txt"
    mkdir -p "$PERIOD_TASK_DIR"
    : > "$combined_file"

    while IFS= read -r raw_line; do
        date=$(echo "$raw_line" | tr -d '[:space:]')
        [[ -z "$date" ]] && continue
        [[ "$date" == \#* ]] && continue
        echo "$date|$SOURCE_ROOT_PART1" >> "$combined_file"
    done < "$LIST_FILE_PART1"

    while IFS= read -r raw_line; do
        date=$(echo "$raw_line" | tr -d '[:space:]')
        [[ -z "$date" ]] && continue
        [[ "$date" == \#* ]] && continue
        echo "$date|$SOURCE_ROOT_PART2" >> "$combined_file"
    done < "$LIST_FILE_PART2"

    sort -t'|' -k1,1 "$combined_file" > "$sorted_file"
    : > "$PERIOD_KEYS_FILE"

    current_key=""
    while IFS='|' read -r date source_root; do
        period_key=$(period_key_from_date "$date")
        if [[ -z "$period_key" ]]; then
            echo "Error: failed to derive period key for date '$date'" >&2
            exit 1
        fi
        if [[ "$period_key" != "$current_key" ]]; then
            echo "$period_key" >> "$PERIOD_KEYS_FILE"
            current_key="$period_key"
        fi
        echo "$date|$source_root" >> "$PERIOD_TASK_DIR/${period_key}.lst"
    done < "$sorted_file"

    n_tasks=$(wc -l < "$PERIOD_KEYS_FILE")
    if [[ "$n_tasks" -lt 1 ]]; then
        echo "Error: no period tasks were created from list files" >&2
        exit 1
    fi

    SUBMIT_ENV_FILE="$LOG_DIR/compute_work_async_histograms_by_lat_band_${submit_stamp}_$$.env"
    {
        printf 'PROJECT_ROOT=%q\n' "$PROJECT_ROOT"
        printf 'SIMULATION=%q\n' "$SIMULATION"
        printf 'LIST_FILE_PART1=%q\n' "$LIST_FILE_PART1"
        printf 'LIST_FILE_PART2=%q\n' "$LIST_FILE_PART2"
        printf 'NTASKS=%q\n' "$NTASKS"
        printf 'CPUS_PER_TASK=%q\n' "$CPUS_PER_TASK"
        printf 'MEM_PER_CPU=%q\n' "$MEM_PER_CPU"
        printf 'RUN_ID=%q\n' "$RUN_ID"
        printf 'LOG_DIR=%q\n' "$LOG_DIR"
        printf 'MANIFEST_DIR=%q\n' "$MANIFEST_DIR"
        printf 'MANIFEST_TOOL=%q\n' "$MANIFEST_TOOL"
        printf 'TARGET_DIR_HIST=%q\n' "$TARGET_DIR_HIST"
        printf 'TARGET_DIR_HIST_BASE=%q\n' "$TARGET_DIR_HIST_BASE"
        printf 'NLAT_BANDS=%q\n' "$NLAT_BANDS"
        printf 'USE_CUSTOM_LAT_BAND_BOUNDS=%q\n' "$USE_CUSTOM_LAT_BAND_BOUNDS"
        printf 'LAT_BAND_BOUNDS=%q\n' "$LAT_BAND_BOUNDS"
        printf 'AGGREGATION_MODE=%q\n' "$AGGREGATION_MODE"
        printf 'OUTPUT_FILE_PREFIX=%q\n' "$OUTPUT_FILE_PREFIX"
        printf 'PERIOD_TASK_DIR=%q\n' "$PERIOD_TASK_DIR"
        printf 'PERIOD_KEYS_FILE=%q\n' "$PERIOD_KEYS_FILE"
    } > "$SUBMIT_ENV_FILE"
    export SUBMIT_ENV_FILE

        manifest_path="$MANIFEST_DIR/histograms/manifest_histograms_${RUN_ID}.json"
        cat << EOF | python3 "$MANIFEST_TOOL" --output "$manifest_path" --workflow-type histograms --run-id "$RUN_ID" --project-root "$PROJECT_ROOT" --launcher "launcher/compute_work_async_histograms_by_lat_band_array.sh" --mode "slurm_array" --root "$TARGET_DIR_HIST_BASE"
{
    "simulation": "$SIMULATION",
    "list_files": {
        "part1": "$LIST_FILE_PART1",
        "part2": "$LIST_FILE_PART2"
    },
    "source_roots": {
        "part1": "$SOURCE_ROOT_PART1",
        "part2": "$SOURCE_ROOT_PART2"
    },
    "resources": {
        "ntasks": "$NTASKS",
        "cpus_per_task": "$CPUS_PER_TASK",
        "mem_per_cpu": "$MEM_PER_CPU"
    },
    "log_dir": "$LOG_DIR",
    "output": {
        "hist_base": "$TARGET_DIR_HIST_BASE",
        "hist_dir": "$TARGET_DIR_HIST"
    },
    "aggregation_mode": "$AGGREGATION_MODE",
    "output_file_prefix": "$OUTPUT_FILE_PREFIX",
    "lat_banding": {
        "enabled": true,
        "nlat_bands": "$NLAT_BANDS",
        "use_custom_bounds": $([[ "$USE_CUSTOM_LAT_BAND_BOUNDS" == "true" ]] && echo true || echo false),
        "custom_bounds": "$LAT_BAND_BOUNDS"
    }
}
EOF

    echo "Submitting histogram-by-lat-band array with $n_tasks period tasks for simulation $SIMULATION (aggregation_mode=$AGGREGATION_MODE)"
    echo "Resources per array task: ntasks=$NTASKS, cpus-per-task=$CPUS_PER_TASK, mem-per-cpu=$MEM_PER_CPU"
    sbatch --array=1-"$n_tasks" --ntasks="$NTASKS" --cpus-per-task="$CPUS_PER_TASK" --mem-per-cpu="$MEM_PER_CPU" \
        --output="$LOG_DIR/compute_work_async_histograms_by_lat_band_%A_%a.out" --error="$LOG_DIR/compute_work_async_histograms_by_lat_band_%A_%a.err" \
        --export=ALL "$SCRIPT_PATH"
    exit $?
fi

if [[ ! -f "$PERIOD_KEYS_FILE" || ! -d "$PERIOD_TASK_DIR" ]]; then
    echo "Error: period task metadata not found in task mode." >&2
    exit 1
fi

period_key=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PERIOD_KEYS_FILE")
if [[ -z "$period_key" ]]; then
    echo "Error: no period key found at task index $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi
period_list="$PERIOD_TASK_DIR/${period_key}.lst"
if [[ ! -f "$period_list" ]]; then
    echo "Error: period list file not found: $period_list" >&2
    exit 1
fi

module purge || true
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-$CPUS_PER_TASK}

COMPUTE_BIN=${COMPUTE_BIN:-$PROJECT_ROOT/local/bin/compute_work_async_histograms_by_lat_band}
if [[ ! -x "$COMPUTE_BIN" && "$COMPUTE_BIN" == "$PROJECT_ROOT/local/bin/compute_work_async_histograms_by_lat_band" ]]; then
    COMPUTE_BIN=$PROJECT_ROOT/bin/compute_work_async_histograms_by_lat_band
fi
if [[ ! -x "$COMPUTE_BIN" ]]; then
    echo "Error: histogram binary not found or not executable." >&2
    echo "Tried: $PROJECT_ROOT/local/bin/compute_work_async_histograms_by_lat_band and $PROJECT_ROOT/bin/compute_work_async_histograms_by_lat_band" >&2
    exit 1
fi

namelist_dir=$PROJECT_ROOT/output/namelists
[[ -d $TARGET_DIR_HIST ]] || mkdir -p $TARGET_DIR_HIST
[[ -d $namelist_dir ]] || mkdir -p $namelist_dir

first_entry=$(head -n 1 "$period_list")
if [[ -z "$first_entry" ]]; then
    echo "Error: period list '$period_list' is empty" >&2
    exit 1
fi
IFS='|' read -r date source_root <<< "$first_entry"
source_dir=$source_root/$date
if [[ ! -d "$source_dir" ]]; then
    echo "Error: source directory not found: $source_dir" >&2
    exit 1
fi

config_file="$namelist_dir/config_hist_by_lat_band_${period_key}_${SLURM_ARRAY_TASK_ID:-0}_$$.nml"
{
    # Build a task-local namelist so each array element is self-contained.
    cat << EOF
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
    path_hist_out   = '$TARGET_DIR_HIST/${OUTPUT_FILE_PREFIX}_${period_key}.nc',
    nlat_bands      = $NLAT_BANDS,
    aggregation_mode = '$AGGREGATION_MODE',
    date_list_file  = '$period_list',
    target_period_key = '$period_key',
EOF

    if [[ "$USE_CUSTOM_LAT_BAND_BOUNDS" == "true" ]]; then
        IFS=',' read -r -a bounds <<< "$LAT_BAND_BOUNDS"
        echo "    use_custom_lat_band_bounds = .true.,"
        echo "    lat_band_bounds_count = ${#bounds[@]},"
        printf '    lat_band_bounds = '
        for ((i=0; i<${#bounds[@]}; i++)); do
            if [[ $i -gt 0 ]]; then
                printf ', '
            fi
            printf '%s' "${bounds[$i]}"
        done
        printf ',\n'
    else
        echo "    use_custom_lat_band_bounds = .false.,"
        echo "    lat_band_bounds_count = 0,"
    fi

    cat << EOF
/
EOF
} > "$config_file"

"$COMPUTE_BIN" "$config_file"
rm -f "$config_file"
