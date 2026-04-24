#!/usr/bin/env bash
set -euo pipefail

# Usage examples:
#   SIMULATION=control ./coarse_grain_and_concat_precip.sh
#   SIMULATION=warming ./coarse_grain_and_concat_precip.sh
#   SIMULATION=control_prate_thresholded ./coarse_grain_and_concat_precip.sh
#   SIMULATION=warming_prate_thresholded ./coarse_grain_and_concat_precip.sh
#   SIMULATION=control_prate_thresholded_by_lat_band ./coarse_grain_and_concat_precip.sh
#   SIMULATION=warming_prate_thresholded_by_lat_band ./coarse_grain_and_concat_precip.sh

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

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

SIMULATION="${SIMULATION:-control}"
INPUT_MODE=""

case "$SIMULATION" in
  control)
    SOURCE_ROOT_PART1="/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history"
    SOURCE_ROOT_PART2="/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp"
    LIST_FILE_PART1="${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_control_part1.txt}"
    LIST_FILE_PART2="${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_control_part2.txt}"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180"
    INPUT_MODE="list"
    ;;
  warming)
    SOURCE_ROOT_PART1="/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    SOURCE_ROOT_PART2="/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    LIST_FILE_PART1="${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt}"
    LIST_FILE_PART2="${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt}"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv"
    INPUT_MODE="list"
    ;;
  control_prate_thresholded)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_prate_thresholded"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_prate_thresholded"
    file_prefix="work_prate_threshold_"
    INPUT_MODE="directory"
    ;;
  warming_prate_thresholded)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_PLUS_4K_CO2_1270ppmv_prate_thresholded"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded"
    file_prefix="work_prate_threshold_"
    INPUT_MODE="directory"
    ;;
  control_prate_thresholded_by_lat_band)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_prate_thresholded_by_lat_band"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_prate_thresholded_by_lat_band"
    file_prefix="work_prate_threshold_by_lat_band_"
    INPUT_MODE="directory"
    ;;
  warming_prate_thresholded_by_lat_band)
    default_src_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_1440x720_PLUS_4K_CO2_1270ppmv_prate_thresholded_by_lat_band"
    default_out_dir="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded_by_lat_band"
    file_prefix="work_prate_threshold_by_lat_band_"
    INPUT_MODE="directory"
    ;;
  *)
    echo "Error: unsupported SIMULATION='$SIMULATION'." >&2
    echo "Use one of: control, warming, control_prate_thresholded, warming_prate_thresholded, control_prate_thresholded_by_lat_band, warming_prate_thresholded_by_lat_band" >&2
    exit 1
    ;;
esac

# Source folder containing input files named like ${file_prefix}YYYYMMDDHH.nc.
# Used for INPUT_MODE=directory and can be overridden via environment variable SRC_DIR.
SRC_DIR="${SRC_DIR:-${default_src_dir:-}}"

# Destination folder for intermediate and final outputs.
# Can be overridden via environment variable OUT_DIR.
OUT_DIR="${OUT_DIR:-$default_out_dir}"

# Input file naming/variable defaults by mode.
if [[ "$INPUT_MODE" == "list" ]]; then
  INPUT_FILENAME="${INPUT_FILENAME:-PRATEsfc_coarse_C3072_1440x720.fre.nc}"
  PRECIP_VAR="${PRECIP_VAR:-PRATEsfc_coarse}"
else
  INPUT_FILENAME="${INPUT_FILENAME:-}"
  PRECIP_VAR="${PRECIP_VAR:-precip}"
fi

# CDO remapping operator. For coarse-graining, remapcon is usually appropriate.
REMAP_METHOD="${REMAP_METHOD:-remapcon}"
TARGET_GRID="${TARGET_GRID:-r360x180}"

# Temporary workspace (can be overridden if needed).
TMP_DIR="${TMP_DIR:-$OUT_DIR/tmp_remap_$$}"

for cmd in cdo ncrcat ncks sed; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Error: required command '$cmd' is not in PATH." >&2
    exit 1
  fi
done

if [[ "$INPUT_MODE" == "directory" ]]; then
  if [[ ! -d "$SRC_DIR" ]]; then
    echo "Error: source directory does not exist: $SRC_DIR" >&2
    exit 1
  fi
else
  if [[ ! -f "$LIST_FILE_PART1" ]]; then
    echo "Error: list file not found: $LIST_FILE_PART1" >&2
    exit 1
  fi
  if [[ ! -f "$LIST_FILE_PART2" ]]; then
    echo "Error: list file not found: $LIST_FILE_PART2" >&2
    exit 1
  fi
fi

mkdir -p "$OUT_DIR"
mkdir -p "$TMP_DIR"

cleanup() {
  rm -rf "$TMP_DIR"
}
trap cleanup EXIT

selected_inputs=()
sorted_inputs=()

if [[ "$INPUT_MODE" == "directory" ]]; then
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

  for in_file in "${sorted_inputs[@]}"; do
    repaired_file="${in_file%.nc}.taxis_repaired.nc"
    if [[ -f "$repaired_file" ]]; then
      selected_inputs+=("$repaired_file")
      echo "Using repaired time axis: $(basename "$repaired_file") (instead of $(basename "$in_file"))"
    else
      selected_inputs+=("$in_file")
    fi
  done
else
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

  date1="${entries[0]#*|}"
  date2="${entries[${#entries[@]}-1]#*|}"

  for i in "${!entries[@]}"; do
    source_root="${entries[$i]%|*}"
    date="${entries[$i]#*|}"
    in_file="$source_root/$date/$INPUT_FILENAME"

    if [[ ! -d "$source_root/$date" ]]; then
      echo "Error: source directory does not exist: $source_root/$date" >&2
      exit 1
    fi
    if [[ ! -f "$in_file" ]]; then
      echo "Error: input file not found: $in_file" >&2
      exit 1
    fi

    repaired_file="${in_file%.nc}.taxis_repaired.nc"
    if [[ -f "$repaired_file" ]]; then
      selected_inputs+=("$repaired_file")
      sorted_inputs+=("$in_file")
      echo "Using repaired time axis: $(basename "$repaired_file") (instead of $(basename "$in_file"))"
    else
      selected_inputs+=("$in_file")
      sorted_inputs+=("$in_file")
    fi
  done
fi

if [[ ${#selected_inputs[@]} -eq 0 ]]; then
  echo "Error: no valid selected inputs found" >&2
  exit 1
fi

remapped_files=()
for i in "${!selected_inputs[@]}"; do
  in_file="${selected_inputs[$i]}"
  if [[ "$INPUT_MODE" == "directory" ]]; then
    base_name="$(basename "${sorted_inputs[$i]}")"
  else
    date_label="${in_file%/*}"
    date_label="${date_label##*/}"
    base_name="${date_label}.nc"
  fi

  # Step 1: extract precipitation variable only.
  extracted_file="$TMP_DIR/extracted_$base_name"
  echo "Extracting ${PRECIP_VAR} from $base_name"
  if ! ncks -O -v "$PRECIP_VAR" "$in_file" "$extracted_file"; then
    echo "Error: failed to extract variable '$PRECIP_VAR' from $in_file" >&2
    exit 1
  fi

  # Unify interface across all modes: output variable should be named 'precip'.
  if [[ "$PRECIP_VAR" != "precip" ]]; then
    if ! ncrename -O -v "${PRECIP_VAR},precip" "$extracted_file"; then
      echo "Error: failed to rename variable '${PRECIP_VAR}' to 'precip' in $extracted_file" >&2
      exit 1
    fi
  fi

  # Step 2: remap extracted precip to target grid.
  out_file="$TMP_DIR/remapped_$base_name"
  echo "Remapping $(basename "$extracted_file") -> $(basename "$out_file")"
  cdo -L "${REMAP_METHOD},${TARGET_GRID}" "$extracted_file" "$out_file"

  remapped_files+=("$out_file")
done

final_file="$OUT_DIR/precip_${date1}_${date2}.nc"

echo "Concatenating ${#remapped_files[@]} files -> $(basename "$final_file")"
ncrcat "${remapped_files[@]}" "$final_file"

echo "Done: $final_file"
