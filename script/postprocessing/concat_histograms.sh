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

# Shared filesystems (e.g., GPFS/Lustre) can trigger HDF5 backend write errors
# unless file locking is disabled. Allow user override when needed.
export HDF5_USE_FILE_LOCKING="${HDF5_USE_FILE_LOCKING:-FALSE}"

# Path resolution (run-aware):
#   Source:
#     1) SRC_DIR (if set) wins.
#     2) Else RUN_ID -> SRC_DIR_BASE/RUN_ID.
#     3) Else legacy flat files in SRC_DIR_BASE -> use SRC_DIR_BASE.
#     4) Else auto-pick latest run subdirectory in SRC_DIR_BASE.
#   Output:
#     1) OUT_DIR (if set) wins.
#     2) Else if RUN_ID is known -> OUT_DIR_BASE/RUN_ID.
#     3) Else -> OUT_DIR_BASE.
#
# Set SRC_DIR and OUT_DIR explicitly to force legacy flat behavior.

SIMULATION="${SIMULATION:-control}"

case "$SIMULATION" in
  control)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms"
    default_file_prefix="hist_"
    default_output_prefix="hist"
    default_date_digits=10
    ;;
  warming)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_PLUS_4K_CO2_1270ppmv"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_PLUS_4K_CO2_1270ppmv"
    default_file_prefix="hist_"
    default_output_prefix="hist"
    default_date_digits=10
    ;;
  control_monthly_by_lat_band)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_monthly_by_lat_band"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_monthly_by_lat_band"
    default_file_prefix="hist_monthly_by_lat_band_"
    default_output_prefix="hist_monthly_by_lat_band"
    default_date_digits=6
    ;;
  warming_monthly_by_lat_band)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_monthly_by_lat_band_PLUS_4K_CO2_1270ppmv"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_monthly_by_lat_band_PLUS_4K_CO2_1270ppmv"
    default_file_prefix="hist_monthly_by_lat_band_"
    default_output_prefix="hist_monthly_by_lat_band"
    default_date_digits=6
    ;;
  *)
    echo "Error: unsupported SIMULATION='$SIMULATION'. Use 'control', 'warming', 'control_monthly_by_lat_band', or 'warming_monthly_by_lat_band'." >&2
    exit 1
    ;;
esac

# Source folder containing histogram files.
# Can be overridden via environment variable SRC_DIR.
RUN_ID="${RUN_ID:-}"
SRC_DIR_BASE="${SRC_DIR_BASE:-$default_src_dir}"
if [[ -z "${SRC_DIR:-}" ]]; then
  if [[ -n "$RUN_ID" ]]; then
    SRC_DIR="$SRC_DIR_BASE/$RUN_ID"
  elif compgen -G "$SRC_DIR_BASE/${default_file_prefix}*.nc" > /dev/null; then
    # Backward compatibility for legacy flat output layout.
    SRC_DIR="$SRC_DIR_BASE"
  else
    latest_run="$(find "$SRC_DIR_BASE" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null | sort | tail -n1)"
    if [[ -n "$latest_run" ]]; then
      SRC_DIR="$SRC_DIR_BASE/$latest_run"
      RUN_ID="$latest_run"
    else
      SRC_DIR="$SRC_DIR_BASE"
    fi
  fi
else
  SRC_DIR="$SRC_DIR"
fi

# Prefix used for default filename matching.
FILE_PREFIX="${FILE_PREFIX:-$default_file_prefix}"
# Glob and regex can be overridden for custom naming conventions.
FILE_GLOB="${FILE_GLOB:-${FILE_PREFIX}*.nc}"
# Keep default regex assignment outside ${var:-...} because unescaped '}' in
# regex quantifiers (e.g., {10}) can terminate parameter expansion early.
if [[ -z "${FILE_DATE_REGEX:-}" ]]; then
  FILE_DATE_REGEX="^${FILE_PREFIX}([0-9]{${default_date_digits}})\\.nc$"
fi

# Destination folder for concatenated output.
# Can be overridden via environment variable OUT_DIR.
OUT_DIR_BASE="${OUT_DIR_BASE:-$default_out_dir}"
if [[ -z "${OUT_DIR:-}" ]]; then
  if [[ -n "$RUN_ID" ]]; then
    OUT_DIR="$OUT_DIR_BASE/$RUN_ID"
  else
    OUT_DIR="$OUT_DIR_BASE"
  fi
else
  OUT_DIR="$OUT_DIR"
fi

# Output filename prefix for histogram outputs.
OUTPUT_PREFIX="${OUTPUT_PREFIX:-$default_output_prefix}"
# NCO format flag. Default to NetCDF-4 to allow unlimited dimensions in
# non-leading positions for multi-dimensional histogram variables.
NCO_FORMAT_FLAG="${NCO_FORMAT_FLAG:--4}"
# Reorder temporary files to put time first before concatenation. This avoids
# ncrcat write failures with some NCO/HDF5 builds when record dims trail.
REORDER_TIME_FIRST="${REORDER_TIME_FIRST:-1}"
# Optional final-file reordering: when enabled and time is the trailing
# dimension in timed variables, reverse to this explicit order.
# With REVERSED_DIM_ORDER='time,nbin_work,nbin_pr,nlat', 4D variables become
# (time,nbin_work,nbin_pr,nlat) and 3D variables become (time,nbin_pr,nlat).
REVERSE_DIMS_IF_TIME_LAST="${REVERSE_DIMS_IF_TIME_LAST:-0}"
REVERSED_DIM_ORDER="${REVERSED_DIM_ORDER:-time,nbin_work,nbin_pr,nlat}"

if ! command -v ncrcat >/dev/null 2>&1; then
  echo "Error: required command 'ncrcat' is not in PATH." >&2
  exit 1
fi

if ! command -v ncks >/dev/null 2>&1; then
  echo "Error: required command 'ncks' is not in PATH." >&2
  exit 1
fi

if ! command -v ncdump >/dev/null 2>&1; then
  echo "Error: required command 'ncdump' is not in PATH." >&2
  exit 1
fi

if ! command -v ncpdq >/dev/null 2>&1; then
  echo "Error: required command 'ncpdq' is not in PATH." >&2
  exit 1
fi

if [[ ! -d "$SRC_DIR" ]]; then
  echo "Error: source directory does not exist: $SRC_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

shopt -s nullglob
inputs=()
candidates=("$SRC_DIR"/$FILE_GLOB)
for f in "${candidates[@]}"; do
  b="$(basename "$f")"
  # Keep only raw single-date histogram files, excluding prior concatenated outputs.
  if [[ "$b" =~ $FILE_DATE_REGEX ]]; then
    inputs+=("$f")
  fi
done
shopt -u nullglob

if [[ ${#inputs[@]} -eq 0 ]]; then
  if [[ ${#candidates[@]} -eq 0 ]]; then
    echo "Error: no input files found for glob '$FILE_GLOB' in $SRC_DIR" >&2
  else
    echo "Error: FILE_GLOB matched ${#candidates[@]} files, but FILE_DATE_REGEX matched none." >&2
    echo "  FILE_GLOB='$FILE_GLOB'" >&2
    echo "  FILE_DATE_REGEX='$FILE_DATE_REGEX'" >&2
  fi
  exit 1
fi

# Sort by date label in filename to preserve chronological order.
IFS=$'\n' sorted_inputs=($(printf '%s\n' "${inputs[@]}" | sort))
unset IFS

first_base="$(basename "${sorted_inputs[0]}")"
last_base="$(basename "${sorted_inputs[${#sorted_inputs[@]}-1]}")"

if [[ "$first_base" =~ $FILE_DATE_REGEX ]]; then
  date1="${BASH_REMATCH[1]}"
else
  echo "Error: cannot parse start date from filename: $first_base" >&2
  exit 1
fi

if [[ "$last_base" =~ $FILE_DATE_REGEX ]]; then
  date2="${BASH_REMATCH[1]}"
else
  echo "Error: cannot parse end date from filename: $last_base" >&2
  exit 1
fi

final_file="$OUT_DIR/${OUTPUT_PREFIX}_${date1}_${date2}.nc"

echo "Concatenating ${#sorted_inputs[@]} files -> $(basename "$final_file")"
echo "Using HDF5_USE_FILE_LOCKING=$HDF5_USE_FILE_LOCKING"

# ncrcat requires at least one record variable (unlimited dimension).
# If input files have only fixed dimensions (e.g., time=1 fixed),
# convert time to a record dimension in temporary copies.
if ncdump -h "${sorted_inputs[0]}" | grep -qE '^\s*time\s*=\s*UNLIMITED\s*;'; then
  ncrcat "$NCO_FORMAT_FLAG" "${sorted_inputs[@]}" "$final_file"
else
  tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/hist_concat.XXXXXX")"
  cleanup() {
    rm -rf "$tmp_dir"
  }
  trap cleanup EXIT

  tmp_inputs=()
  for f in "${sorted_inputs[@]}"; do
    tf_rec="$tmp_dir/rec_$(basename "$f")"
    tf="$tmp_dir/$(basename "$f")"
    ncks "$NCO_FORMAT_FLAG" --mk_rec_dmn time "$f" "$tf_rec"
    if [[ "$REORDER_TIME_FIRST" == "1" ]]; then
      ncpdq "$NCO_FORMAT_FLAG" -a time "$tf_rec" "$tf"
    else
      mv "$tf_rec" "$tf"
    fi
    tmp_inputs+=("$tf")
  done

  ncrcat "$NCO_FORMAT_FLAG" "${tmp_inputs[@]}" "$final_file"
fi

if [[ "$REVERSE_DIMS_IF_TIME_LAST" == "1" ]]; then
  # Detect variables where time is the last dimension (excluding time(time)).
  if ncdump -h "$final_file" | grep -Eq '^[[:space:]]*[[:alnum:]_]+[[:space:]]+[[:alnum:]_]+\([^)]*,time\)[[:space:]]*;'; then
    echo "Detected variables with trailing time dimension; reordering final output to '$REVERSED_DIM_ORDER'"
    final_tmp="$final_file.reordered.tmp"
    ncpdq "$NCO_FORMAT_FLAG" -a "$REVERSED_DIM_ORDER" "$final_file" "$final_tmp"
    mv -f "$final_tmp" "$final_file"
  fi
fi

echo "Done: $final_file"
