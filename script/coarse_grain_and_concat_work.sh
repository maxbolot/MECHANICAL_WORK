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

SIMULATION="${SIMULATION:-control}"

case "$SIMULATION" in
  control)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180"
    file_prefix="work_"
    ;;
  warming)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_PLUS_4K_CO2_1270ppmv"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv"
    file_prefix="work_"
    ;;
  control_prate_thresholded)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_prate_thresholded"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_prate_thresholded"
    file_prefix="work_prate_threshold_"
    ;;
  warming_prate_thresholded)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_PLUS_4K_CO2_1270ppmv_prate_thresholded"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded"
    file_prefix="work_prate_threshold_"
    ;;
  *)
    echo "Error: unsupported SIMULATION='$SIMULATION'." >&2
    echo "Use one of: control, warming, control_prate_thresholded, warming_prate_thresholded" >&2
    exit 1
    ;;
esac

# Source folder containing input files named like ${file_prefix}YYYYMMDDHH.nc.
# Can be overridden via environment variable SRC_DIR.
SRC_DIR="${SRC_DIR:-$default_src_dir}"

# Destination folder for intermediate and final outputs.
# Can be overridden via environment variable OUT_DIR.
OUT_DIR="${OUT_DIR:-$default_out_dir}"

# CDO remapping operator. For coarse-graining, remapcon is usually appropriate.
REMAP_METHOD="${REMAP_METHOD:-remapcon}"
TARGET_GRID="${TARGET_GRID:-r360x180}"

# Temporary workspace (can be overridden if needed).
TMP_DIR="${TMP_DIR:-$OUT_DIR/tmp_remap_$$}"

for cmd in cdo ncrcat; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Error: required command '$cmd' is not in PATH." >&2
    exit 1
  fi
done

if [[ ! -d "$SRC_DIR" ]]; then
  echo "Error: source directory does not exist: $SRC_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
mkdir -p "$TMP_DIR"

cleanup() {
  rm -rf "$TMP_DIR"
}
trap cleanup EXIT

shopt -s nullglob
inputs=()
for f in "$SRC_DIR"/${file_prefix}[0-9]*.nc; do
  b="$(basename "$f")"
  # Keep only canonical single-date files. Repaired files are selected as
  # aliases below to avoid duplicate entries in the timeline.
  if [[ "$b" =~ ^${file_prefix}[0-9]{10}\.nc$ ]]; then
    inputs+=("$f")
  fi
done
shopt -u nullglob

if [[ ${#inputs[@]} -eq 0 ]]; then
  echo "Error: no input files matching ${file_prefix}YYYYMMDDHH.nc found in $SRC_DIR" >&2
  exit 1
fi

# Sort by date label in filename to preserve chronological order.
IFS=$'\n' sorted_inputs=($(printf '%s\n' "${inputs[@]}" | sort))
unset IFS

first_base="$(basename "${sorted_inputs[0]}")"
last_base="$(basename "${sorted_inputs[${#sorted_inputs[@]}-1]}")"

if [[ "$first_base" =~ ^${file_prefix}([0-9]{10})\.nc$ ]]; then
  date1="${BASH_REMATCH[1]}"
else
  echo "Error: cannot parse start date from filename: $first_base" >&2
  exit 1
fi

if [[ "$last_base" =~ ^${file_prefix}([0-9]{10})\.nc$ ]]; then
  date2="${BASH_REMATCH[1]}"
else
  echo "Error: cannot parse end date from filename: $last_base" >&2
  exit 1
fi

selected_inputs=()
for in_file in "${sorted_inputs[@]}"; do
  repaired_file="${in_file%.nc}.taxis_repaired.nc"
  if [[ -f "$repaired_file" ]]; then
    selected_inputs+=("$repaired_file")
    echo "Using repaired time axis: $(basename "$repaired_file") (instead of $(basename "$in_file"))"
  else
    selected_inputs+=("$in_file")
  fi
done

remapped_files=()
for i in "${!selected_inputs[@]}"; do
  in_file="${selected_inputs[$i]}"
  base_name="$(basename "${sorted_inputs[$i]}")"
  out_file="$TMP_DIR/$base_name"

  echo "Remapping $base_name -> $(basename "$out_file")"
  cdo -L "${REMAP_METHOD},${TARGET_GRID}" "$in_file" "$out_file"
  remapped_files+=("$out_file")
done

final_file="$OUT_DIR/work_${date1}_${date2}.nc"

echo "Concatenating ${#remapped_files[@]} files -> $(basename "$final_file")"
ncrcat "${remapped_files[@]}" "$final_file"

echo "Done: $final_file"