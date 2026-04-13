#!/bin/bash
#SBATCH -A geoclim
#SBATCH --time=04:00:00
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

# Launch mode (no SLURM_ARRAY_TASK_ID): read two date lists and submit this script as a job array.
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
    n_tasks=$(( n_part1 + n_part2 ))
    if [[ "$n_tasks" -lt 1 ]]; then
        echo "Error: both list files are empty" >&2
        exit 1
    fi

    echo "Submitting job array with $n_tasks tasks for simulation $SIMULATION"
    echo "Resources per array task: ntasks=$NTASKS, cpus-per-task=$CPUS_PER_TASK, mem-per-cpu=$MEM_PER_CPU"
    sbatch --array=1-"$n_tasks" --ntasks="$NTASKS" --cpus-per-task="$CPUS_PER_TASK" --mem-per-cpu="$MEM_PER_CPU" \
            --output="$LOG_DIR/compute_work_%A_%a.out" --error="$LOG_DIR/compute_work_%A_%a.err" \
            --export=ALL,PROJECT_ROOT="$PROJECT_ROOT",SIMULATION="$SIMULATION",LIST_FILE_PART1="$LIST_FILE_PART1",LIST_FILE_PART2="$LIST_FILE_PART2",NTASKS="$NTASKS",CPUS_PER_TASK="$CPUS_PER_TASK",MEM_PER_CPU="$MEM_PER_CPU",LOG_DIR="$LOG_DIR",TARGET_DIR_COMPUTE="$TARGET_DIR_COMPUTE",TARGET_DIR_HISTOGRAMS="$TARGET_DIR_HISTOGRAMS" "$SCRIPT_PATH"
    exit $?
fi

# Task mode: map array index onto the concatenation of part1 and part2 lists.
if [[ ! -f "$LIST_FILE_PART1" || ! -f "$LIST_FILE_PART2" ]]; then
    echo "Error: list files not found in task mode." >&2
    exit 1
fi

mapfile -t dates_part1 < <(sed -E '/^[[:space:]]*($|#)/d; s/^[[:space:]]+//; s/[[:space:]]+$//' "$LIST_FILE_PART1")
mapfile -t dates_part2 < <(sed -E '/^[[:space:]]*($|#)/d; s/^[[:space:]]+//; s/[[:space:]]+$//' "$LIST_FILE_PART2")
n_part1=${#dates_part1[@]}

if (( SLURM_ARRAY_TASK_ID <= n_part1 )); then
    date=${dates_part1[SLURM_ARRAY_TASK_ID-1]}
    source_root=$SOURCE_ROOT_PART1
else
    idx_part2=$(( SLURM_ARRAY_TASK_ID - n_part1 ))
    if (( idx_part2 < 1 || idx_part2 > ${#dates_part2[@]} )); then
        echo "Error: no date found at combined line $SLURM_ARRAY_TASK_ID" >&2
        exit 1
    fi
    date=${dates_part2[idx_part2-1]}
    source_root=$SOURCE_ROOT_PART2
fi

module purge || true
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-$CPUS_PER_TASK}
source_dir=$source_root/$date
target_dir_compute=$TARGET_DIR_COMPUTE
target_dir_histograms=$TARGET_DIR_HISTOGRAMS

if [[ ! -d "$source_dir" ]]; then
    echo "Error: source directory not found: $source_dir" >&2
    exit 1
fi

# Prefer local/bin build output, fallback to legacy bin path.
COMPUTE_WORK_BIN=${COMPUTE_WORK_BIN:-$PROJECT_ROOT/local/bin/compute_work_async}
if [[ ! -x "$COMPUTE_WORK_BIN" && "$COMPUTE_WORK_BIN" == "$PROJECT_ROOT/local/bin/compute_work_async" ]]; then
    COMPUTE_WORK_BIN=$PROJECT_ROOT/bin/compute_work_async
fi
if [[ ! -x "$COMPUTE_WORK_BIN" ]]; then
    echo "Error: compute_work binary not found or not executable." >&2
    echo "Tried: ./local/bin/compute_work_async and ./bin/compute_work_async" >&2
    exit 1
fi

namelist_dir=$PROJECT_ROOT/output/namelists
[[ -d $target_dir_compute ]] || mkdir -p $target_dir_compute
[[ -d $target_dir_histograms ]] || mkdir -p $target_dir_histograms
[[ -d $namelist_dir ]] || mkdir -p $namelist_dir

config_file="$namelist_dir/config_${date}_${SLURM_ARRAY_TASK_ID:-0}_$$.nml"
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

"$COMPUTE_WORK_BIN" "$config_file"
rm -f "$config_file"
