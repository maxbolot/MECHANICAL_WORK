#!/usr/bin/env bash
set -euo pipefail

# Check for time-axis drift in per-date work files.
#
# For each canonical file named work_YYYYMMDDHH.nc:
#   1) parse the filename start timestamp,
#   2) read the actual first timestamp from the NetCDF time axis,
#   3) if the drift from filename start exceeds MAX_DIFF_DAYS, report DRIFT
#      in the output and optionally remove the source file when
#      REMOVE_SOURCE_FILES=1.
#
# The script never modifies source files unless REMOVE_SOURCE_FILES=1.

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

# Match the same simulation naming used in launchers and post-processing scripts.
SIMULATION="${SIMULATION:-control}"

case "$SIMULATION" in
  control)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720"
    file_prefix="work_"
    ;;
  warming)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_PLUS_4K_CO2_1270ppmv"
    file_prefix="work_"
    ;;
  control_prate_thresholded)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_prate_thresholded"
    file_prefix="work_prate_threshold_"
    ;;
  warming_prate_thresholded)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_PLUS_4K_CO2_1270ppmv_prate_thresholded"
    file_prefix="work_prate_threshold_"
    ;;
  *)
    echo "Error: unsupported SIMULATION='$SIMULATION'." >&2
    echo "Use one of: control, warming, control_prate_thresholded, warming_prate_thresholded" >&2
    exit 1
    ;;
esac

# SRC_DIR can still be explicitly overridden for ad-hoc runs.
SRC_DIR="${SRC_DIR:-$default_src_dir}"
# Maximum tolerated drift between expected and observed first timestamp.
MAX_DIFF_DAYS="${MAX_DIFF_DAYS:-5}"
# DRY_RUN=1 prints planned deletions without removing files.
DRY_RUN="${DRY_RUN:-0}"
# REMOVE_SOURCE_FILES=1 deletes source files with drift > MAX_DIFF_DAYS.
REMOVE_SOURCE_FILES="${REMOVE_SOURCE_FILES:-0}"

# Hard-fail early if required runtime tools are missing.
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

# MAX_DIFF_DAYS accepts integer or decimal days.
if ! [[ "$MAX_DIFF_DAYS" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Error: MAX_DIFF_DAYS must be numeric (got: $MAX_DIFF_DAYS)" >&2
  exit 1
fi

# Keep DRY_RUN values strict for predictable behavior.
if [[ "$DRY_RUN" != "0" && "$DRY_RUN" != "1" ]]; then
  echo "Error: DRY_RUN must be 0 or 1 (got: $DRY_RUN)" >&2
  exit 1
fi

if [[ "$REMOVE_SOURCE_FILES" != "0" && "$REMOVE_SOURCE_FILES" != "1" ]]; then
  # Keep mode values strict to avoid accidental deletions from typos.
  echo "Error: REMOVE_SOURCE_FILES must be 0 or 1 (got: $REMOVE_SOURCE_FILES)" >&2
  exit 1
fi

# Convert day threshold to seconds for integer arithmetic comparisons.
threshold_sec=$(awk -v d="$MAX_DIFF_DAYS" 'BEGIN{printf "%.0f", d*86400.0}')

# Collect only canonical per-date files to avoid reprocessing concatenated products.
# The file_prefix differs between simulations (e.g. work_ vs work_prate_threshold_).
shopt -s nullglob
inputs=()
for f in "$SRC_DIR"/${file_prefix}[0-9]*.nc; do
  b="$(basename "$f")"
  if [[ "$b" =~ ^${file_prefix}([0-9]{10})\.nc$ ]]; then
    inputs+=("$f")
  fi
done
shopt -u nullglob

if [[ ${#inputs[@]} -eq 0 ]]; then
  echo "Error: no files matching ${file_prefix}YYYYMMDDHH.nc found in $SRC_DIR" >&2
  exit 1
fi

# Sort by embedded date so reports are chronological and reproducible.
IFS=$'\n' sorted_inputs=($(printf '%s\n' "${inputs[@]}" | sort))
unset IFS

# Summary counters printed at the end.
total=0
drifted=0
skipped=0
removed=0

for in_file in "${sorted_inputs[@]}"; do
  total=$((total + 1))
  base="$(basename "$in_file")"

  if [[ ! "$base" =~ ^${file_prefix}([0-9]{10})\.nc$ ]]; then
    skipped=$((skipped + 1))
    echo "SKIP   $base  (name does not match expected pattern)"
    continue
  fi

  label="${BASH_REMATCH[1]}"
  # Derive nominal file start time from the filename stamp: ${file_prefix}YYYYMMDDHH.nc
  start_iso="${label:0:4}-${label:4:2}-${label:6:2} ${label:8:2}:00:00"
  start_epoch=$(date -u -d "$start_iso" +%s)
  start_iso_fmt=$(date -u -d "@$start_epoch" '+%Y-%m-%dT%H:%M:%S')

  # Read first timestamp from file time axis.
  # showtimestamp output is tokenized; the first token is used as reference.
  actual_raw=$(cdo -s showtimestamp "$in_file" 2>/dev/null | awk '{print $1}')
  if [[ -z "$actual_raw" ]]; then
    skipped=$((skipped + 1))
    echo "SKIP   $base  (unable to read timestamp with cdo showtimestamp)"
    continue
  fi

  # Normalize timestamp format for GNU date parsing.
  actual_iso=${actual_raw/Z/}
  actual_iso=${actual_iso/T/ }
  actual_epoch=$(date -u -d "$actual_iso" +%s)

  # Compare observed first timestamp against the filename start.
  diff_sec=$(( actual_epoch - start_epoch ))
  if (( diff_sec < 0 )); then
    diff_sec=$(( -diff_sec ))
  fi
  diff_days=$(awk -v s="$diff_sec" 'BEGIN{printf "%.6f", s/86400.0}')

  if (( diff_sec > threshold_sec )); then
    drifted=$((drifted + 1))
    if [[ "$REMOVE_SOURCE_FILES" == "1" ]]; then
      if [[ "$DRY_RUN" == "1" ]]; then
        echo "DRIFT  $base  actual=$actual_raw filename_start=$start_iso_fmt diff_days=$diff_days  [would remove]"
      else
        rm -f -- "$in_file"
        echo "DRIFT  $base  actual=$actual_raw filename_start=$start_iso_fmt diff_days=$diff_days  [removed]"
        removed=$((removed + 1))
      fi
    else
      echo "DRIFT  $base  actual=$actual_raw filename_start=$start_iso_fmt diff_days=$diff_days"
    fi
  else
    echo "OK     $base  actual=$actual_raw filename_start=$start_iso_fmt diff_days=$diff_days"
  fi
done

echo
echo "Summary: total=$total drifted=$drifted removed=$removed skipped=$skipped threshold_days=$MAX_DIFF_DAYS dry_run=$DRY_RUN remove_source_files=$REMOVE_SOURCE_FILES"
