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
LOG_DIR=${LOG_DIR:-$PROJECT_ROOT/logs}

# Launch mode (no SLURM_ARRAY_TASK_ID): read date list and submit this script as a job array.
LIST_FILE=${LIST_FILE:-${1:-$SCRIPT_DIR/list.txt}}
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    mkdir -p "$LOG_DIR"

    if [[ ! -f "$LIST_FILE" ]]; then
        echo "Error: list file not found: $LIST_FILE" >&2
        exit 1
    fi

    n_tasks=$(wc -l < "$LIST_FILE")
    if [[ "$n_tasks" -lt 1 ]]; then
        echo "Error: list file is empty: $LIST_FILE" >&2
        exit 1
    fi

    echo "Submitting job array with $n_tasks tasks using dates from $LIST_FILE"
    echo "Resources per array task: ntasks=$NTASKS, cpus-per-task=$CPUS_PER_TASK"
    sbatch --array=1-"$n_tasks" --ntasks="$NTASKS" --cpus-per-task="$CPUS_PER_TASK" \
            --output="$LOG_DIR/compute_work_%A_%a.out" --error="$LOG_DIR/compute_work_%A_%a.err" \
            --export=ALL,PROJECT_ROOT="$PROJECT_ROOT",LIST_FILE="$LIST_FILE",NTASKS="$NTASKS",CPUS_PER_TASK="$CPUS_PER_TASK",LOG_DIR="$LOG_DIR" "$SCRIPT_PATH"
    exit $?
fi

# Task mode: map array index -> date (expects one yyyymmdd date per line).
date=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST_FILE" | tr -d '[:space:]')
if [[ -z "$date" ]]; then
    echo "Error: no date found at line $SLURM_ARRAY_TASK_ID in $LIST_FILE" >&2
    exit 1
fi

module purge || true
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-$CPUS_PER_TASK}
# source_dir=/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history/$date
# target_dir_compute=/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720
# target_dir_histograms=/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms
source_dir=/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp/$date
target_dir_compute=/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_PLUS_4K_CO2_1270ppmv
target_dir_histograms=/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_PLUS_4K_CO2_1270ppmv

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

[[ -d $target_dir_compute ]] || mkdir -p $target_dir_compute
[[ -d $target_dir_histograms ]] || mkdir -p $target_dir_histograms

config_file="$PROJECT_ROOT/config_${date}_${SLURM_ARRAY_TASK_ID:-0}_$$.nml"
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
