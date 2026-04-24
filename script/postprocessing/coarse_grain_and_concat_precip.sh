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
PRATE_THRESHOLDED=0
LAT_BAND_THRESHOLDED=0
THRESHOLD_FILE=""

case "$SIMULATION" in
  control)
    SOURCE_ROOT_PART1="/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history"
    SOURCE_ROOT_PART2="/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp"
    LIST_FILE_PART1="${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_control_part1.txt}"
    LIST_FILE_PART2="${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_control_part2.txt}"
    DEFAULT_OUT_DIR="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180"
    ;;
  control_prate_thresholded)
    SOURCE_ROOT_PART1="/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history"
    SOURCE_ROOT_PART2="/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp"
    LIST_FILE_PART1="${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_control_part1.txt}"
    LIST_FILE_PART2="${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_control_part2.txt}"
    DEFAULT_OUT_DIR="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_prate_thresholded"
    PRATE_THRESHOLDED=1
    THRESHOLD_FILE="${THRESHOLD_FILE:-$PROJECT_ROOT/output/thresholds/thresholds_control.txt}"
    ;;
  control_prate_thresholded_by_lat_band)
    SOURCE_ROOT_PART1="/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history"
    SOURCE_ROOT_PART2="/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp"
    LIST_FILE_PART1="${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_control_part1.txt}"
    LIST_FILE_PART2="${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_control_part2.txt}"
    DEFAULT_OUT_DIR="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_prate_thresholded_by_lat_band"
    PRATE_THRESHOLDED=1
    LAT_BAND_THRESHOLDED=1
    THRESHOLD_FILE="${THRESHOLD_FILE:-$PROJECT_ROOT/output/thresholds/thresholds_control_by_lat_band.txt}"
    ;;
  warming)
    SOURCE_ROOT_PART1="/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    SOURCE_ROOT_PART2="/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    LIST_FILE_PART1="${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt}"
    LIST_FILE_PART2="${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt}"
    DEFAULT_OUT_DIR="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv"
    ;;
  warming_prate_thresholded)
    SOURCE_ROOT_PART1="/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    SOURCE_ROOT_PART2="/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    LIST_FILE_PART1="${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt}"
    LIST_FILE_PART2="${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt}"
    DEFAULT_OUT_DIR="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded"
    PRATE_THRESHOLDED=1
    THRESHOLD_FILE="${THRESHOLD_FILE:-$PROJECT_ROOT/output/thresholds/thresholds_warming.txt}"
    ;;
  warming_prate_thresholded_by_lat_band)
    SOURCE_ROOT_PART1="/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    SOURCE_ROOT_PART2="/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp"
    LIST_FILE_PART1="${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt}"
    LIST_FILE_PART2="${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt}"
    DEFAULT_OUT_DIR="/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded_by_lat_band"
    PRATE_THRESHOLDED=1
    LAT_BAND_THRESHOLDED=1
    THRESHOLD_FILE="${THRESHOLD_FILE:-$PROJECT_ROOT/output/thresholds/thresholds_warming_by_lat_band.txt}"
    ;;
  *)
    echo "Error: unsupported SIMULATION='$SIMULATION'." >&2
    echo "Use one of: control, warming, control_prate_thresholded, warming_prate_thresholded, control_prate_thresholded_by_lat_band, warming_prate_thresholded_by_lat_band" >&2
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

if [[ "$PRATE_THRESHOLDED" == "1" && ! -f "$THRESHOLD_FILE" ]]; then
  echo "Error: threshold file not found for SIMULATION='$SIMULATION': $THRESHOLD_FILE" >&2
  exit 1
fi

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

# For thresholded mode: threshold + daily-aggregate BEFORE remap/concat
# (matching Fortran workflow: compute_work_async_prate_threshold.f90 outputs daily aggregates with percentile dim,
#  then coarse_grain_and_concat_work.sh remaps and concats those daily files).
# For non-thresholded: remap each file, then concat, then daily-aggregate.

if [[ "$PRATE_THRESHOLDED" == "1" ]]; then
  # Thresholded workflow (mirrors Fortran compute_work_async_prate_threshold.f90):
  #   per input file: mask all thresholds at once → daily aggregate → remap
  #   then: ncrcat per-threshold files → ncecat into percentile dimension
  # Never concatenate raw native-resolution files; process each file individually
  # to avoid OOM.
  if ! command -v ncap2 >/dev/null 2>&1; then
    echo "Error: ncap2 is required for PRATE thresholded output." >&2
    exit 1
  fi
  if ! command -v ncks >/dev/null 2>&1; then
    echo "Error: ncks is required for PRATE thresholded output." >&2
    exit 1
  fi

  # Parse thresholds first.
  # Non-lat-band file format: percentile threshold
  # Lat-band file format: # lat_band p ... ; then rows: band_id t1 t2 ... tN
  pvals_arr=()
  tvals_arr=()
  declare -A tvals_2d=()

  if [[ "$LAT_BAND_THRESHOLDED" == "1" ]]; then
    mapfile -t pvals_arr < <(awk '
      /^#/ {
        for (i=1; i<=NF; i++) {
          if ($i == "p" && (i+1) <= NF) {
            print $(i+1)
          } else if ($i ~ /^p[0-9]/) {
            v=$i
            sub(/^p/, "", v)
            print v
          }
        }
        if (length($0) > 0) exit
      }
    ' "$THRESHOLD_FILE")

    if [[ "${#pvals_arr[@]}" -eq 0 ]]; then
      mapfile -t pvals_arr < <(awk 'NF && $1 !~ /^#/ {
        for (i=2; i<=NF; i++) print i-1;
        exit
      }' "$THRESHOLD_FILE")
    fi

    nthreshold=${#pvals_arr[@]}
    if [[ "$nthreshold" -eq 0 ]]; then
      echo "Error: failed to parse percentile columns from lat-band threshold file: $THRESHOLD_FILE" >&2
      exit 1
    fi

    while read -r band ip tval; do
      tvals_2d["$band,$ip"]="$tval"
    done < <(awk -v n="${nthreshold}" 'NF && $1 !~ /^#/ {
      band=$1
      for (i=1; i<=n; i++) {
        print band, i, $(i+1)
      }
    }' "$THRESHOLD_FILE")

    missing=0
    for ((ilat=1; ilat<=18; ilat++)); do
      for ((ip=1; ip<=nthreshold; ip++)); do
        if [[ -z "${tvals_2d[$ilat,$ip]:-}" ]]; then
          echo "Error: missing threshold value for band=$ilat percentile_index=$ip in $THRESHOLD_FILE" >&2
          missing=1
          break
        fi
      done
      [[ "$missing" -eq 1 ]] && break
    done
    if [[ "$missing" -eq 1 ]]; then
      exit 1
    fi
  else
    while read -r pval tval; do
      pvals_arr+=("$pval")
      tvals_arr+=("$tval")
    done < <(awk 'NF && $1 !~ /^#/ {print $1, $2}' "$THRESHOLD_FILE")

    nthreshold=${#pvals_arr[@]}
    if [[ "$nthreshold" -eq 0 ]]; then
      echo "Error: no valid threshold rows found in $THRESHOLD_FILE" >&2
      exit 1
    fi
  fi

  # After daily sums are computed on native grid, convert p* to conditional means
  # by dividing by c* at each grid point/day and writing FillValue-like sentinel
  # where no threshold exceedance occurred.
  daily_cond_expr=""
  for ((ip=0; ip<nthreshold; ip++)); do
    daily_cond_expr+="where(c${ip}>0){p${ip}=p${ip}/c${ip};} elsewhere {p${ip}=-9999.0;};"
  done

  echo "Found $nthreshold thresholds from $THRESHOLD_FILE"

  # Detect precipitation variable from the first input file
  first_source_root="${entries[0]%|*}"
  first_date="${entries[0]#*|}"
  first_file="$first_source_root/$first_date/$INPUT_FILENAME"
  if [[ ! -f "$first_file" ]]; then
    echo "Error: first input file not found: $first_file" >&2
    exit 1
  fi
  ncks_header="$TMP_DIR/ncks_header.txt"
  ncks -m "$first_file" > "$ncks_header" 2>&1 || true
  precip_var=$(awk '
    /^[[:space:]]*(float|double|short|int|byte|char)[[:space:]]+[A-Za-z_][A-Za-z0-9_]*\(/ {
      v=$2; sub(/\(.*/, "", v); gsub(/;/, "", v);
      if (v != "time" && v != "lon" && v != "lat" && v != "percentile" && v != "prate_threshold") {
        print v; exit;
      }
    }
  ' "$ncks_header")
  if [[ -z "$precip_var" ]]; then
    echo "Error: unable to detect precipitation variable in $first_file" >&2
    cat "$ncks_header" >&2
    exit 1
  fi
  echo "Detected precipitation variable: $precip_var"

  # Build ncap2 expression: compute one masked precip variable and one event-count
  # indicator per threshold in a single pass.
  # For lat-band mode, thresholds depend on latitude band.
  ncap_expr=""
  if [[ "$LAT_BAND_THRESHOLDED" == "1" ]]; then
    lat_var="grid_yt_coarse"
    if ! grep -q "[[:space:]]${lat_var}(" "$ncks_header"; then
      if grep -q "[[:space:]]lat(" "$ncks_header"; then
        lat_var="lat"
      else
        echo "Error: unable to detect latitude coordinate variable (tried grid_yt_coarse and lat) in $first_file" >&2
        exit 1
      fi
    fi
    for ((ip=0; ip<nthreshold; ip++)); do
      ncap_expr+="p${ip}=0.0*${precip_var};"
      ncap_expr+="c${ip}=0.0*${precip_var};"
      for ((ilat=1; ilat<=18; ilat++)); do
        south=$(( -90 + (ilat - 1) * 10 ))
        north=$(( south + 10 ))
        tval="${tvals_2d[$ilat,$((ip+1))]}"
        if [[ "$ilat" -lt 18 ]]; then
          ncap_expr+="p${ip}=where(((${lat_var}>=${south})&&(${lat_var}<${north})),${precip_var}*(${precip_var}>=${tval}),p${ip});"
          ncap_expr+="c${ip}=where(((${lat_var}>=${south})&&(${lat_var}<${north})),(${precip_var}>=${tval}),c${ip});"
        else
          ncap_expr+="p${ip}=where(((${lat_var}>=${south})&&(${lat_var}<=${north})),${precip_var}*(${precip_var}>=${tval}),p${ip});"
          ncap_expr+="c${ip}=where(((${lat_var}>=${south})&&(${lat_var}<=${north})),(${precip_var}>=${tval}),c${ip});"
        fi
      done
    done
  else
    for ((ip=0; ip<nthreshold; ip++)); do
      ncap_expr+="p${ip}=${precip_var}*(${precip_var}>=${tvals_arr[$ip]});"
      ncap_expr+="c${ip}=(${precip_var}>=${tvals_arr[$ip]});"
    done
  fi

  # Initialise per-threshold arrays to collect remapped per-date files
  for ((ip=0; ip<nthreshold; ip++)); do
    declare -a "thresh_files_${ip}=()"
  done

  # Process each input file: mask → daily aggregate → remap
  nentries=${#entries[@]}
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

    echo "[$((i+1))/$nentries] Processing ${date}/${INPUT_FILENAME}"

    # Step 1: apply all threshold masks in one ncap2 pass → multi-variable file
    masked="$TMP_DIR/masked_$(printf '%06d' "$i").nc"
    if ! ncap2 -O -s "$ncap_expr" "$in_file" "$masked" 2>&1; then
      echo "Error: ncap2 masking failed for $in_file" >&2
      exit 1
    fi
    # Drop the original precip variable; keep only p0..p(N-1)
    if ! ncks -O -x -v "${precip_var}" "$masked" "$masked" 2>&1; then
      echo "Error: ncks drop of ${precip_var} failed for $in_file" >&2
      exit 1
    fi

    # Step 2: daily aggregate on native grid.
    # Use daily sums for both masked precip (p*) and event counts (c*), then
    # compute conditional means later as p_sum / c_sum.
    daily="$TMP_DIR/daily_$(printf '%06d' "$i").nc"
    if ! cdo -s -L -daysum "$masked" "$daily" 2>&1; then
      echo "Error: cdo -daysum failed for $in_file" >&2
      exit 1
    fi
    rm -f "$masked"

    if ! ncap2 -O -s "$daily_cond_expr" "$daily" "$daily" 2>&1; then
      echo "Error: ncap2 conditional-mean step failed for $in_file" >&2
      exit 1
    fi

    # Step 3: remap to target grid
    remapped="$TMP_DIR/remapped_$(printf '%06d' "$i").nc"
    if [[ "$FORCE_SETGRID" == "1" ]]; then
      if ! cdo -s -L "${REMAP_METHOD},${TARGET_GRID}" -setgrid,"${SOURCE_GRID}" \
          "$daily" "$remapped" 2>&1; then
        echo "setgrid remap failed for $date; retrying without setgrid"
        rm -f "$remapped"
        if ! cdo -s -L "${REMAP_METHOD},${TARGET_GRID}" "$daily" "$remapped" 2>&1; then
          echo "Error: both remap attempts failed for $date" >&2
          exit 1
        fi
      fi
    else
      if ! cdo -s -L "${REMAP_METHOD},${TARGET_GRID}" "$daily" "$remapped" 2>&1; then
        echo "Direct remap failed for $date; retrying with setgrid"
        rm -f "$remapped"
        if ! cdo -s -L "${REMAP_METHOD},${TARGET_GRID}" -setgrid,"${SOURCE_GRID}" \
            "$daily" "$remapped" 2>&1; then
          echo "Error: both remap attempts failed for $date" >&2
          exit 1
        fi
      fi
    fi
    rm -f "$daily"

    # Record this remapped file for each threshold
    for ((ip=0; ip<nthreshold; ip++)); do
      ref="thresh_files_${ip}[@]"
      eval "thresh_files_${ip}+=('$remapped')"
    done
  done

  # For each threshold: extract p/c variables from all remapped files and concat.
  # p* is already a daily conditional mean from the native-grid daily step.
  per_thresh_concat=()
  for ((ip=0; ip<nthreshold; ip++)); do
    pval="${pvals_arr[$ip]}"
    tval="${tvals_arr[$ip]}"
    echo "Threshold $((ip+1))/$nthreshold: extracting p${ip}, c${ip} (percentile=$pval threshold=$tval)"

    ref="thresh_files_${ip}[@]"
    mapfile -t tfiles < <(printf '%s\n' "${!ref}" | sort -u)

    # Extract this threshold's masked-sum and count variables from every remapped file, then concat
    extracted_files=()
    for f in "${tfiles[@]}"; do
      ext="$TMP_DIR/ext_p${ip}_$(basename "$f")"
      if ! ncks -O -v "p${ip},c${ip}" "$f" "$ext" 2>&1; then
        echo "Error: ncks extraction of p${ip},c${ip} failed for $f" >&2
        exit 1
      fi
      extracted_files+=("$ext")
    done

    concat_thresh="$TMP_DIR/concat_thresh_${ip}.nc"
    if ! ncrcat -O "${extracted_files[@]}" "$concat_thresh" 2>&1; then
      echo "Error: ncrcat failed for threshold $ip" >&2
      exit 1
    fi
    # Clean up extracted files
    rm -f "${extracted_files[@]}"

    if ! ncrename -O -v "p${ip},precip" -v "c${ip},event_count" "$concat_thresh" 2>&1; then
      echo "Error: ncrename failed for threshold $ip" >&2
      exit 1
    fi

    if ! ncks -O -v "precip,event_count,time,lon,lat" "$concat_thresh" "$concat_thresh" 2>&1; then
      echo "Error: ncks keep of precip/event_count failed for threshold $ip" >&2
      exit 1
    fi

    per_thresh_concat+=("$concat_thresh")
  done

  # Clean up all remapped per-date files (shared across thresholds)
  for i in "${!entries[@]}"; do
    rm -f "$TMP_DIR/remapped_$(printf '%06d' "$i").nc"
  done

  # Stack per-threshold files along a new 'percentile' dimension
  final_file="$OUT_DIR/precip_${date1}_${date2}.nc"
  stacked="$TMP_DIR/precip_stacked.nc"
  echo "Stacking $nthreshold per-threshold files along percentile dimension (ncecat)"
  if ! ncecat -O "${per_thresh_concat[@]}" "$stacked" 2>&1; then
    echo "Error: ncecat failed to stack per-threshold files" >&2
    exit 1
  fi
  rm -f "${per_thresh_concat[@]}"

  # ncecat names the new outer dimension 'record'; rename to 'percentile'
  ncrename -O -d record,percentile "$stacked"

  # Assign percentile coordinate values
  pvals_ncap="{$(IFS=','; echo "${pvals_arr[*]}")}"
  ncap2 -O -s "percentile[percentile]=${pvals_ncap};" "$stacked" "$final_file"
  rm -f "$stacked"

  ncatted -O \
    -a long_name,percentile,o,c,"precipitation percentile threshold" \
    -a units,percentile,o,c,"1" \
    -a long_name,precip,o,c,"precipitation rate filtered by precipitation threshold" \
    -a units,precip,o,c,"kg m-2 s-1" \
    -a coordinates,precip,o,c,"percentile lon lat" \
    -a time_stat,precip,o,c,"daily_conditional_mean" \
    -a time_aggregation,precip,o,c,"mean over threshold exceedance events within each day" \
    -a _FillValue,precip,o,d,-9999.0 \
    -a long_name,event_count,o,c,"number of events exceeding precipitation threshold" \
    -a units,event_count,o,c,"count" \
    -a coordinates,event_count,o,c,"percentile lon lat" \
    -a time_stat,event_count,o,c,"daily_count" \
    -a time_aggregation,event_count,o,c,"sum of threshold exceedance events over native timesteps within each day" \
    -a threshold_file,global,o,c,"$THRESHOLD_FILE" \
    -a daily_stat,global,o,c,"mixed_by_variable" \
    -a daily_aggregation,global,o,c,"mean over native timesteps within each day" \
    "$final_file"

  if [[ "$LAT_BAND_THRESHOLDED" == "1" ]]; then
    for ((ip=1; ip<=nthreshold; ip++)); do
      pval="${pvals_arr[$((ip-1))]}"
      ptxt=$(awk -v p="$pval" 'BEGIN {printf "%.4f", 100.0*p}')
      packed=$(echo "$ptxt" | sed -E 's/0+$//; s/\.$//; s/\.//g')
      if [[ -z "$packed" ]]; then
        packed=$((ip))
      fi
      for ((ilat=1; ilat<=18; ilat++)); do
        tval="${tvals_2d[$ilat,$ip]}"
        ncatted -O -a "band$(printf '%02d' "$ilat")_p${packed}_threshold,global,o,d,${tval}" "$final_file"
      done
    done
  else
    while read -r pval tval; do
      ptxt=$(awk -v p="$pval" 'BEGIN {printf "%.4f", 100.0*p}')
      packed=$(echo "$ptxt" | sed -E 's/0+$//; s/\.$//; s/\.//g')
      ncatted -O -a "p${packed}_threshold,global,o,d,${tval}" "$final_file"
    done < <(awk 'NF && $1 !~ /^#/ {print $1, $2}' "$THRESHOLD_FILE")
  fi

  echo "Thresholded percentile output written: $final_file"

else
  # Non-thresholded workflow: remap each file → concat → daily aggregate
  echo "Non-thresholded workflow: remapping each file individually"
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
    echo "[$((i+1))/${#entries[@]}] Remapping ${date}/${INPUT_FILENAME} -> $(basename "$out_file")"
    if [[ "$FORCE_SETGRID" == "1" ]]; then
      if ! cdo -L "${REMAP_METHOD},${TARGET_GRID}" -setgrid,"${SOURCE_GRID}" "$in_file" "$out_file" 2>&1; then
        echo "setgrid remap failed for $in_file; retrying direct remap without setgrid"
        rm -f "$out_file"
        if ! cdo -L "${REMAP_METHOD},${TARGET_GRID}" "$in_file" "$out_file" 2>&1; then
          echo "Error: both remap attempts failed for $in_file" >&2
          exit 1
        fi
      fi
    else
      if ! cdo -L "${REMAP_METHOD},${TARGET_GRID}" "$in_file" "$out_file" 2>&1; then
        echo "Direct remap failed for $in_file; retrying with SOURCE_GRID=$SOURCE_GRID"
        rm -f "$out_file"
        if ! cdo -L "${REMAP_METHOD},${TARGET_GRID}" -setgrid,"${SOURCE_GRID}" "$in_file" "$out_file" 2>&1; then
          echo "Error: both remap attempts failed for $in_file" >&2
          exit 1
        fi
      fi
    fi
    remapped_files+=("$out_file")
  done

  final_file="$OUT_DIR/precip_${date1}_${date2}.nc"
  concat_file="$TMP_DIR/precip_concat_${date1}_${date2}.nc"

  echo "Concatenating ${#remapped_files[@]} remapped files"
  if ! ncrcat "${remapped_files[@]}" "$concat_file" 2>&1; then
    echo "Error: ncrcat failed to concatenate remapped files" >&2
    exit 1
  fi
  echo "Concatenation complete: $concat_file"

  # Daily aggregation: one-pass step to enforce daily means on concatenated output.
  # IMPORTANT: Use -daymean instead of dayavg to preserve first-of-day timestamp
  # (consistent with Fortran compute_work_async.f90 which uses floor(time) grouping).
  DAILY_AGGREGATE="${DAILY_AGGREGATE:-1}"
  if [[ "$DAILY_AGGREGATE" == "1" ]]; then
    echo "Computing daily aggregates from concatenated output (preserving first-of-day timestamps)"
    temp_daily="$TMP_DIR/precip_daily_${date1}_${date2}.nc"
    if ! cdo -L -daymean "$concat_file" "$temp_daily" 2>&1; then
      echo "Error: cdo -daymean failed" >&2
      exit 1
    fi
    cp "$temp_daily" "$final_file"
    echo "Output is daily-aggregated with first-of-day timestamps: $final_file"
  else
    echo "Daily aggregation skipped (DAILY_AGGREGATE=0)"
    cp "$concat_file" "$final_file"
  fi
fi

