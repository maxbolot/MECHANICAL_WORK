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

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

SIMULATION="${SIMULATION:-control}"

case "$SIMULATION" in
  control)
    SOURCE_ROOT_PART1="/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history"
    SOURCE_ROOT_PART2="/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp"
    LIST_FILE_PART1="${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_control_part1.txt}"
    LIST_FILE_PART2="${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_control_part2.txt}"
    DEFAULT_OUT_DIR="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180"
    ;;
  warming)
    SOURCE_ROOT_PART1="/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    SOURCE_ROOT_PART2="/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    LIST_FILE_PART1="${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt}"
    LIST_FILE_PART2="${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt}"
    DEFAULT_OUT_DIR="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv"
    ;;
  *)
    echo "Error: unsupported SIMULATION='$SIMULATION'." >&2
    echo "Use one of: control, warming" >&2
    exit 1
    ;;
esac

INPUT_FILENAME="${INPUT_FILENAME:-PRATEsfc_coarse_C3072_1440x720.fre.nc}"
OUT_DIR="${OUT_DIR:-$DEFAULT_OUT_DIR}"
REMAP_METHOD="${REMAP_METHOD:-remapcon}"
TARGET_GRID="${TARGET_GRID:-r360x180}"
SOURCE_GRID="${SOURCE_GRID:-r1440x720}"
FORCE_SETGRID="${FORCE_SETGRID:-1}"
TMP_DIR="${TMP_DIR:-$OUT_DIR/tmp_remap_$$}"

for cmd in cdo ncrcat; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Error: required command '$cmd' is not in PATH." >&2
    exit 1
  fi
done

if [[ ! -f "$LIST_FILE_PART1" ]]; then
  echo "Error: list file not found: $LIST_FILE_PART1" >&2
  exit 1
fi

if [[ ! -f "$LIST_FILE_PART2" ]]; then
  echo "Error: list file not found: $LIST_FILE_PART2" >&2
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

date1="${entries[0]#*|}"
date2="${entries[${#entries[@]}-1]#*|}"

mkdir -p "$OUT_DIR"
mkdir -p "$TMP_DIR"

cleanup() {
  rm -rf "$TMP_DIR"
}
trap cleanup EXIT

remapped_files=()
for i in "${!entries[@]}"; do
  source_root="${entries[$i]%|*}"
  date="${entries[$i]#*|}"
  source_dir="$source_root/$date"
  in_file="$source_dir/$INPUT_FILENAME"

  if [[ ! -d "$source_dir" ]]; then
    echo "Error: source directory does not exist: $source_dir" >&2
    exit 1
  fi

  if [[ ! -f "$in_file" ]]; then
    echo "Error: input file not found: $in_file" >&2
    exit 1
  fi

  out_file="$TMP_DIR/$(printf '%06d' "$i")_${date}.nc"
  echo "Remapping ${date}/${INPUT_FILENAME} -> $(basename "$out_file")"
  if [[ "$FORCE_SETGRID" == "1" ]]; then
    if ! cdo -L "${REMAP_METHOD},${TARGET_GRID}" -setgrid,"${SOURCE_GRID}" "$in_file" "$out_file"; then
      echo "setgrid remap failed for $in_file; retrying direct remap without setgrid"
      rm -f "$out_file"
      cdo -L "${REMAP_METHOD},${TARGET_GRID}" "$in_file" "$out_file"
    fi
  else
    if ! cdo -L "${REMAP_METHOD},${TARGET_GRID}" "$in_file" "$out_file"; then
      echo "Direct remap failed for $in_file; retrying with SOURCE_GRID=$SOURCE_GRID"
      rm -f "$out_file"
      cdo -L "${REMAP_METHOD},${TARGET_GRID}" -setgrid,"${SOURCE_GRID}" "$in_file" "$out_file"
    fi
  fi
  remapped_files+=("$out_file")
done

final_file="$OUT_DIR/precip_${date1}_${date2}.nc"

echo "Concatenating ${#remapped_files[@]} files -> $(basename "$final_file")"
ncrcat "${remapped_files[@]}" "$final_file"

# Daily aggregation: one-pass step to enforce daily means on concatenated output.
# IMPORTANT: Use -daymean instead of dayavg to preserve first-of-day timestamp
# (consistent with Fortran compute_work_async.f90 which uses floor(time) grouping).
DAILY_AGGREGATE="${DAILY_AGGREGATE:-1}"
if [[ "$DAILY_AGGREGATE" == "1" ]]; then
  echo "Computing daily aggregates from concatenated output (preserving first-of-day timestamps)"
  temp_daily="$final_file.daily.tmp"
  # -daymean groups by calendar day and selects first timestep of each day.
  # This aligns with the Fortran program behavior.
  cdo -L -daymean "$final_file" "$temp_daily"
  mv "$temp_daily" "$final_file"
  echo "Output is now daily-aggregated with first-of-day timestamps: $final_file"
else
  echo "Daily aggregation skipped (DAILY_AGGREGATE=0): $final_file"
fi

echo "Done: $final_file"
