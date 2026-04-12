#!/usr/bin/env bash
set -euo pipefail

# Ensure module command is available in non-interactive shells, then load tools.
if ! type module >/dev/null 2>&1; then
  if [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  elif [[ -f /usr/share/Modules/init/bash ]]; then
    # shellcheck disable=SC1091
    source /usr/share/Modules/init/bash
  fi
fi

if ! type module >/dev/null 2>&1; then
  echo "Error: environment modules are unavailable; cannot run module load." >&2
  exit 1
fi

module purge
module load intel/2021.1.2 hdf5/intel-2021.1/1.10.6 netcdf/intel-2021.1/hdf5-1.10.6/4.7.4 cdo/netcdf-4.7.4/hdf5-1.10.6/2.0.1 nco/netcdf-4.7.4/hdf5-1.10.6/5.0.3

SRC_DIR="${SRC_DIR:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720}"
MAX_DIFF_DAYS="${MAX_DIFF_DAYS:-5}"
CUTOFF_DATE="${CUTOFF_DATE:-2021-05-27 00:00:00}"
DRY_RUN="${DRY_RUN:-0}"

for cmd in cdo date awk sort; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Error: required command '$cmd' is not in PATH." >&2
    exit 1
  fi
done

if [[ ! -d "$SRC_DIR" ]]; then
  echo "Error: source directory does not exist: $SRC_DIR" >&2
  exit 1
fi

if ! [[ "$MAX_DIFF_DAYS" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Error: MAX_DIFF_DAYS must be numeric (got: $MAX_DIFF_DAYS)" >&2
  exit 1
fi

if [[ "$DRY_RUN" != "0" && "$DRY_RUN" != "1" ]]; then
  echo "Error: DRY_RUN must be 0 or 1 (got: $DRY_RUN)" >&2
  exit 1
fi

cutoff_epoch=$(date -u -d "$CUTOFF_DATE" +%s)
threshold_sec=$(awk -v d="$MAX_DIFF_DAYS" 'BEGIN{printf "%.0f", d*86400.0}')

shopt -s nullglob
inputs=()
for f in "$SRC_DIR"/work_*.nc; do
  b="$(basename "$f")"
  if [[ "$b" =~ ^work_[0-9]{10}\.nc$ ]]; then
    inputs+=("$f")
  fi
done
shopt -u nullglob

if [[ ${#inputs[@]} -eq 0 ]]; then
  echo "Error: no files matching work_YYYYMMDDHH.nc found in $SRC_DIR" >&2
  exit 1
fi

IFS=$'\n' sorted_inputs=($(printf '%s\n' "${inputs[@]}" | sort))
unset IFS

total=0
repaired=0
skipped=0

for in_file in "${sorted_inputs[@]}"; do
  total=$((total + 1))
  base="$(basename "$in_file")"

  if [[ ! "$base" =~ ^work_([0-9]{10})\.nc$ ]]; then
    skipped=$((skipped + 1))
    echo "SKIP   $base  (name does not match expected pattern)"
    continue
  fi

  label="${BASH_REMATCH[1]}"
  start_iso="${label:0:4}-${label:4:2}-${label:6:2} ${label:8:2}:00:00"
  start_epoch=$(date -u -d "$start_iso" +%s)

  # A priori cadence: 5-day before 2021-05-27, 1-day on/after 2021-05-27.
  if (( start_epoch < cutoff_epoch )); then
    cadence_days=5
    midpoint_sec=$((5 * 86400 / 2))
  else
    cadence_days=1
    midpoint_sec=$((1 * 86400 / 2))
  fi

  expected_epoch=$((start_epoch + midpoint_sec))
  expected_iso=$(date -u -d "@$expected_epoch" '+%Y-%m-%dT%H:%M:%S')
  expected_set_date=$(date -u -d "@$expected_epoch" '+%Y-%m-%d')
  expected_set_time=$(date -u -d "@$expected_epoch" '+%H:%M:%S')

  # Read first timestamp from file time axis.
  actual_raw=$(cdo -s showtimestamp "$in_file" 2>/dev/null | awk '{print $1}')
  if [[ -z "$actual_raw" ]]; then
    skipped=$((skipped + 1))
    echo "SKIP   $base  (unable to read timestamp with cdo showtimestamp)"
    continue
  fi

  actual_iso=${actual_raw/Z/}
  actual_iso=${actual_iso/T/ }
  actual_epoch=$(date -u -d "$actual_iso" +%s)

  diff_sec=$(( actual_epoch - expected_epoch ))
  if (( diff_sec < 0 )); then
    diff_sec=$(( -diff_sec ))
  fi
  diff_days=$(awk -v s="$diff_sec" 'BEGIN{printf "%.6f", s/86400.0}')

  if (( diff_sec > threshold_sec )); then
    out_file="${in_file%.nc}.taxis_repaired.nc"
    if [[ "$DRY_RUN" == "1" ]]; then
      echo "REPAIR $base  actual=$actual_raw expected=$expected_iso diff_days=$diff_days -> $(basename "$out_file")"
    else
      cdo -O -L -settaxis,"$expected_set_date","$expected_set_time","${cadence_days}day" "$in_file" "$out_file"
      echo "REPAIR $base  actual=$actual_raw expected=$expected_iso diff_days=$diff_days -> $(basename "$out_file")"
    fi
    repaired=$((repaired + 1))
  else
    echo "OK     $base  actual=$actual_raw expected=$expected_iso diff_days=$diff_days"
  fi
done

echo
echo "Summary: total=$total repaired=$repaired skipped=$skipped threshold_days=$MAX_DIFF_DAYS dry_run=$DRY_RUN"
