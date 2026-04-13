#!/bin/bash
#SBATCH -A geoclim
#SBATCH --time=02:00:00
#SBATCH --qos=cimes-short
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mbolot@princeton.edu

set -euo pipefail

# In submission mode BASH_SOURCE[0] is the real script path; in task mode Slurm
# copies the script to /var/spool/slurmd, so PROJECT_ROOT is wrong if recomputed.
# The sbatch --export line pins the correct value; only derive it when not already set.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}
SCRIPT_PATH="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"

NTASKS=${NTASKS:-1}
CPUS_PER_TASK=${CPUS_PER_TASK:-1}
MEM_PER_CPU=${MEM_PER_CPU:-3G}
LOG_DIR=${LOG_DIR:-$PROJECT_ROOT/logs}
# Number of compute_work instances launched concurrently inside one array task.
INNER_PARALLEL=${INNER_PARALLEL:-1}

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

if ! [[ "$INNER_PARALLEL" =~ ^[0-9]+$ ]] || [[ "$INNER_PARALLEL" -lt 1 ]]; then
    echo "Error: INNER_PARALLEL must be a positive integer (got: $INNER_PARALLEL)" >&2
    exit 1
fi

# Launch mode (no SLURM_ARRAY_TASK_ID): read date lists and submit this script as a job array.
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    mkdir -p "$LOG_DIR"

    if [[ ! -f "$LIST_FILE_PART1" ]]; then
        echo "Error: list file not found: $LIST_FILE_PART1" >&2
        exit 1
    fi

    if [[ ! -f "$LIST_FILE_PART2" ]]; then
        echo "Error: list file not found: $LIST_FILE_PART2" >&2
        exit 1
    fi

    n_part1=$(sed -E '/^[[:space:]]*($|#)/d' "$LIST_FILE_PART1" | wc -l)
    n_part2=$(sed -E '/^[[:space:]]*($|#)/d' "$LIST_FILE_PART2" | wc -l)
    n_dates=$(( n_part1 + n_part2 ))
    if [[ "$n_dates" -lt 1 ]]; then
        echo "Error: both list files are empty" >&2
        exit 1
    fi

    # Each array task handles up to INNER_PARALLEL dates.
    n_tasks=$(( (n_dates + INNER_PARALLEL - 1) / INNER_PARALLEL ))

    echo "Submitting job array with $n_tasks tasks for simulation $SIMULATION"
    echo "Dates to process: $n_dates"
    echo "Inner parallel runs per array task: $INNER_PARALLEL"
    echo "Resources per array task: ntasks=$NTASKS, cpus-per-task=$CPUS_PER_TASK, mem-per-cpu=$MEM_PER_CPU"
    sbatch --array=1-"$n_tasks" --ntasks="$NTASKS" --cpus-per-task="$CPUS_PER_TASK" --mem-per-cpu="$MEM_PER_CPU" \
            --output="$LOG_DIR/compute_work_%A_%a.out" --error="$LOG_DIR/compute_work_%A_%a.err" \
            --export=ALL,PROJECT_ROOT="$PROJECT_ROOT",SIMULATION="$SIMULATION",LIST_FILE_PART1="$LIST_FILE_PART1",LIST_FILE_PART2="$LIST_FILE_PART2",NTASKS="$NTASKS",CPUS_PER_TASK="$CPUS_PER_TASK",MEM_PER_CPU="$MEM_PER_CPU",LOG_DIR="$LOG_DIR",INNER_PARALLEL="$INNER_PARALLEL",TARGET_DIR_COMPUTE="$TARGET_DIR_COMPUTE",TARGET_DIR_HISTOGRAMS="$TARGET_DIR_HISTOGRAMS" "$SCRIPT_PATH"
    exit $?
fi

# Task mode: map array index -> a chunk of dates (one yyyymmdd per non-empty, non-comment line).
if [[ ! -f "$LIST_FILE_PART1" || ! -f "$LIST_FILE_PART2" ]]; then
    echo "Error: list files not found in task mode." >&2
    exit 1
fi
mapfile -t dates_part1 < <(sed -E '/^[[:space:]]*($|#)/d; s/^[[:space:]]+//; s/[[:space:]]+$//' "$LIST_FILE_PART1")
mapfile -t dates_part2 < <(sed -E '/^[[:space:]]*($|#)/d; s/^[[:space:]]+//; s/[[:space:]]+$//' "$LIST_FILE_PART2")
entries=()
for date in "${dates_part1[@]}"; do
    entries+=("$SOURCE_ROOT_PART1|$date")
done
for date in "${dates_part2[@]}"; do
    entries+=("$SOURCE_ROOT_PART2|$date")
done
if [[ "${#entries[@]}" -eq 0 ]]; then
    echo "Error: no valid dates found in either list file" >&2
    exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID is not set in task mode" >&2
    exit 1
fi

start_idx=$(( (SLURM_ARRAY_TASK_ID - 1) * INNER_PARALLEL ))
end_idx=$(( start_idx + INNER_PARALLEL - 1 ))
if [[ "$start_idx" -ge "${#entries[@]}" ]]; then
    echo "Info: no dates assigned for array task $SLURM_ARRAY_TASK_ID" >&2
    exit 0
fi
if [[ "$end_idx" -ge "${#entries[@]}" ]]; then
    end_idx=$(( ${#entries[@]} - 1 ))
fi

# Print task-level summary so failures can be diagnosed at a glance in the log.
n_total="${#entries[@]}"
chunk_size=$(( end_idx - start_idx + 1 ))
echo "=== Task ${SLURM_ARRAY_TASK_ID} / dates ${start_idx}-${end_idx} (1-based: $((start_idx+1))-$((end_idx+1))) ==="
echo "    Total dates in list : $n_total"
echo "    Dates for this task : $chunk_size"
echo "    Assigned entries    : ${entries[*]:$start_idx:$chunk_size}"
echo "    INNER_PARALLEL      : $INNER_PARALLEL"
echo "    OMP threads/instance: (to be set after module load)"

module purge || true
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

# Split CPU threads across inner jobs (rounded up so each instance gets at least one).
threads_per_instance=$(( (${SLURM_CPUS_PER_TASK:-$CPUS_PER_TASK} + INNER_PARALLEL - 1) / INNER_PARALLEL ))
if [[ "$threads_per_instance" -lt 1 ]]; then
    threads_per_instance=1
fi
export OMP_NUM_THREADS=$threads_per_instance
echo "    OMP_NUM_THREADS     : $OMP_NUM_THREADS"

target_dir_compute=$TARGET_DIR_COMPUTE
target_dir_histograms=$TARGET_DIR_HISTOGRAMS

# Prefer local/bin build output, fallback to legacy bin path.
COMPUTE_WORK_BIN=${COMPUTE_WORK_BIN:-$PROJECT_ROOT/local/bin/compute_work_async}
if [[ ! -x "$COMPUTE_WORK_BIN" && "$COMPUTE_WORK_BIN" == "$PROJECT_ROOT/local/bin/compute_work_async" ]]; then
    COMPUTE_WORK_BIN=$PROJECT_ROOT/bin/compute_work_async
fi
if [[ ! -x "$COMPUTE_WORK_BIN" ]]; then
    echo "Error: compute_work binary not found or not executable." >&2
    echo "Tried: $PROJECT_ROOT/local/bin/compute_work_async" >&2
    echo "  and: $PROJECT_ROOT/bin/compute_work_async" >&2
    echo "PROJECT_ROOT resolved to: $PROJECT_ROOT" >&2
    exit 1
fi

namelist_dir=$PROJECT_ROOT/output/namelists
[[ -d $target_dir_compute ]] || mkdir -p $target_dir_compute
[[ -d $target_dir_histograms ]] || mkdir -p $target_dir_histograms
[[ -d $LOG_DIR ]] || mkdir -p $LOG_DIR
[[ -d $namelist_dir ]] || mkdir -p $namelist_dir

run_one_date() {
    local date="$1"
    local source_root="$2"
    local source_dir="$source_root/$date"
    local config_file="$namelist_dir/config_${date}_${SLURM_ARRAY_TASK_ID}_$$.nml"
    local date_log="$LOG_DIR/compute_work_${SLURM_ARRAY_JOB_ID:-job}_${SLURM_ARRAY_TASK_ID}_${date}.log"

    if [[ ! -d "$source_dir" ]]; then
        echo "Error: source directory not found: $source_dir" >&2
        return 1
    fi

    # Isolated namelist per date avoids write races between background jobs.
    cat << EOF > "$config_file"
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
    path_work_out   = '$target_dir_compute/work_$date.nc',
    path_hist_out   = '$target_dir_histograms/hist_$date.nc',
/
EOF

    # Keep per-date logs to make failures easy to locate.
    "$COMPUTE_WORK_BIN" "$config_file" > "$date_log" 2>&1
    rm -f "$config_file"
}

pids=()
for ((i=start_idx; i<=end_idx; i++)); do
    # Launch assigned dates concurrently inside this single array task.
    IFS='|' read -r source_root date <<< "${entries[i]}"
    run_one_date "$date" "$source_root" &
    pids+=("$!")
done

rc=0
for pid in "${pids[@]}"; do
    # Preserve failure if any inner process exits non-zero.
    if ! wait "$pid"; then
        rc=1
    fi
done

exit "$rc"
