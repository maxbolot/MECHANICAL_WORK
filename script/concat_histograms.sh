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

# Source folder containing histogram files named like hist_YYYYMMDDHH.nc
# SRC_DIR="${SRC_DIR:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms}"
SRC_DIR="${SRC_DIR:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_PLUS_4K_CO2_1270ppmv}"

# Destination folder for concatenated output.
# OUT_DIR="${OUT_DIR:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms}"
OUT_DIR="${OUT_DIR:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_PLUS_4K_CO2_1270ppmv}"

# Output filename prefix for histogram outputs.
OUTPUT_PREFIX="${OUTPUT_PREFIX:-hist}"

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

if [[ ! -d "$SRC_DIR" ]]; then
  echo "Error: source directory does not exist: $SRC_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

shopt -s nullglob
inputs=()
for f in "$SRC_DIR"/hist_*.nc; do
  b="$(basename "$f")"
  # Keep only raw single-date histogram files, excluding prior concatenated outputs.
  if [[ "$b" =~ ^hist_[0-9]{10}\.nc$ ]]; then
    inputs+=("$f")
  fi
done
shopt -u nullglob

if [[ ${#inputs[@]} -eq 0 ]]; then
  echo "Error: no input files matching hist_YYYYMMDDHH.nc found in $SRC_DIR" >&2
  exit 1
fi

# Sort by date label in filename to preserve chronological order.
IFS=$'\n' sorted_inputs=($(printf '%s\n' "${inputs[@]}" | sort))
unset IFS

first_base="$(basename "${sorted_inputs[0]}")"
last_base="$(basename "${sorted_inputs[${#sorted_inputs[@]}-1]}")"

if [[ "$first_base" =~ ^hist_([0-9]{10})\.nc$ ]]; then
  date1="${BASH_REMATCH[1]}"
else
  echo "Error: cannot parse start date from filename: $first_base" >&2
  exit 1
fi

if [[ "$last_base" =~ ^hist_([0-9]{10})\.nc$ ]]; then
  date2="${BASH_REMATCH[1]}"
else
  echo "Error: cannot parse end date from filename: $last_base" >&2
  exit 1
fi

final_file="$OUT_DIR/${OUTPUT_PREFIX}_${date1}_${date2}.nc"

echo "Concatenating ${#sorted_inputs[@]} files -> $(basename "$final_file")"

# ncrcat requires at least one record variable (unlimited dimension).
# If input files have only fixed dimensions (e.g., time=1 fixed),
# convert time to a record dimension in temporary copies.
if ncdump -h "${sorted_inputs[0]}" | grep -qE '^\s*time\s*=\s*UNLIMITED\s*;'; then
  ncrcat "${sorted_inputs[@]}" "$final_file"
else
  tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/hist_concat.XXXXXX")"
  cleanup() {
    rm -rf "$tmp_dir"
  }
  trap cleanup EXIT

  tmp_inputs=()
  for f in "${sorted_inputs[@]}"; do
    tf="$tmp_dir/$(basename "$f")"
    ncks --mk_rec_dmn time "$f" "$tf"
    tmp_inputs+=("$tf")
  done

  ncrcat "${tmp_inputs[@]}" "$final_file"
fi

echo "Done: $final_file"
