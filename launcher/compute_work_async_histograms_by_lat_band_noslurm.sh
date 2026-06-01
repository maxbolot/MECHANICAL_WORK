#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
RUN_ID=${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}

SIMULATION=${SIMULATION:-control_monthly_by_lat_band}

# Select input roots, date lists, and output target by experiment.
case "$SIMULATION" in
    control_monthly_by_lat_band)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp
        LIST_FILE_PART1=${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_control_part1.txt}
        LIST_FILE_PART2=${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_control_part2.txt}
        TARGET_DIR_HIST_BASE=${TARGET_DIR_HIST_BASE:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_monthly_by_lat_band}
        ;;
    warming_monthly_by_lat_band)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        LIST_FILE_PART1=${LIST_FILE_PART1:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt}
        LIST_FILE_PART2=${LIST_FILE_PART2:-$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt}
        TARGET_DIR_HIST_BASE=${TARGET_DIR_HIST_BASE:-/scratch/gpfs/mbolot/results/GLOBALFV3/work_histograms_monthly_by_lat_band_PLUS_4K_CO2_1270ppmv}
        ;;
    *)
        echo "Error: unsupported SIMULATION='$SIMULATION'. Use 'control_monthly_by_lat_band' or 'warming_monthly_by_lat_band'." >&2
        exit 1
        ;;
esac

TARGET_DIR_HIST=${TARGET_DIR_HIST:-$TARGET_DIR_HIST_BASE/$RUN_ID}

# Optional latitude-band controls
NLAT_BANDS=${NLAT_BANDS:-18}
USE_CUSTOM_LAT_BAND_BOUNDS=${USE_CUSTOM_LAT_BAND_BOUNDS:-false}
# Comma-separated numeric boundaries, only used when USE_CUSTOM_LAT_BAND_BOUNDS=true
# Example: LAT_BAND_BOUNDS="-90,-60,-30,0,30,60,90"
LAT_BAND_BOUNDS=${LAT_BAND_BOUNDS:-}
AGGREGATION_MODE=${AGGREGATION_MODE:-monthly}
GROUP_BY_PERIOD=${GROUP_BY_PERIOD:-0}

case "$AGGREGATION_MODE" in
    monthly)
        OUTPUT_FILE_PREFIX=${OUTPUT_FILE_PREFIX:-hist_monthly_by_lat_band}
        ;;
    daily)
        OUTPUT_FILE_PREFIX=${OUTPUT_FILE_PREFIX:-hist_daily_by_lat_band}
        ;;
    *)
        echo "Error: unsupported AGGREGATION_MODE='$AGGREGATION_MODE'. Use 'daily' or 'monthly'." >&2
        exit 1
        ;;
esac

period_key_from_date() {
    local date_key="$1"
    case "$AGGREGATION_MODE" in
        monthly)
            echo "${date_key:0:6}"
            ;;
        daily)
            echo "${date_key:0:8}"
            ;;
        *)
            echo ""
            ;;
    esac
}

module purge || true
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-3}

LOG_DIR=${LOG_DIR:-$PROJECT_ROOT/logs/slurm}
MANIFEST_DIR=${MANIFEST_DIR:-$PROJECT_ROOT/logs/manifests}
MANIFEST_TOOL=${MANIFEST_TOOL:-$PROJECT_ROOT/python/workflow_manifest.py}
NAMELIST_DIR=$PROJECT_ROOT/output/namelists

COMPUTE_BIN=${COMPUTE_BIN:-$PROJECT_ROOT/local/bin/compute_work_async_histograms_by_lat_band}
if [[ ! -x "$COMPUTE_BIN" && "$COMPUTE_BIN" == "$PROJECT_ROOT/local/bin/compute_work_async_histograms_by_lat_band" ]]; then
    COMPUTE_BIN=$PROJECT_ROOT/bin/compute_work_async_histograms_by_lat_band
fi

for f in "$LIST_FILE_PART1" "$LIST_FILE_PART2"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: list file not found: $f" >&2
        exit 1
    fi
done

for d in "$SOURCE_ROOT_PART1" "$SOURCE_ROOT_PART2"; do
    if [[ ! -d "$d" ]]; then
        echo "Error: source root directory not found: $d" >&2
        exit 1
    fi
done

if [[ ! -x "$COMPUTE_BIN" ]]; then
    echo "Error: histogram binary not found or not executable." >&2
    echo "Tried: $PROJECT_ROOT/local/bin/compute_work_async_histograms_by_lat_band and $PROJECT_ROOT/bin/compute_work_async_histograms_by_lat_band" >&2
    exit 1
fi

mkdir -p "$TARGET_DIR_HIST" "$LOG_DIR" "$MANIFEST_DIR" "$NAMELIST_DIR"

manifest_path="$MANIFEST_DIR/histograms/manifest_histograms_${RUN_ID}.json"
cat << EOF | python3 "$MANIFEST_TOOL" --output "$manifest_path" --workflow-type histograms --run-id "$RUN_ID" --project-root "$PROJECT_ROOT" --launcher "launcher/compute_work_async_histograms_by_lat_band_noslurm.sh" --mode "noslurm" --root "$TARGET_DIR_HIST_BASE"
{
    "simulation": "$SIMULATION",
    "list_files": {
        "part1": "$LIST_FILE_PART1",
        "part2": "$LIST_FILE_PART2"
    },
    "source_roots": {
        "part1": "$SOURCE_ROOT_PART1",
        "part2": "$SOURCE_ROOT_PART2"
    },
    "resources": {
        "omp_num_threads": "$OMP_NUM_THREADS"
    },
    "log_dir": "$LOG_DIR",
    "output": {
        "hist_base": "$TARGET_DIR_HIST_BASE",
        "hist_dir": "$TARGET_DIR_HIST"
    },
    "aggregation_mode": "$AGGREGATION_MODE",
    "output_file_prefix": "$OUTPUT_FILE_PREFIX",
    "group_by_period": "$GROUP_BY_PERIOD",
    "lat_banding": {
        "enabled": true,
        "nlat_bands": "$NLAT_BANDS",
        "use_custom_bounds": $([[ "$USE_CUSTOM_LAT_BAND_BOUNDS" == "true" ]] && echo true || echo false),
        "custom_bounds": "$LAT_BAND_BOUNDS"
    }
}
EOF

RUN_LOG="$LOG_DIR/compute_work_async_histograms_by_lat_band_noslurm_$(date +%Y%m%d_%H%M%S).log"
echo "[$(date +%F\ %T)] Starting serial run for simulation: $SIMULATION" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Using binary: $COMPUTE_BIN" | tee -a "$RUN_LOG"

make_namelist() {
    # Emit a per-date namelist so the Fortran executable can run standalone.
    local source_dir="$1"
    local out_hist="$2"
    local config_file="$3"

    {
        cat << EOF
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
    path_hist_out   = '$out_hist',
    nlat_bands      = $NLAT_BANDS,
    aggregation_mode = '$AGGREGATION_MODE',
EOF

        if [[ "$USE_CUSTOM_LAT_BAND_BOUNDS" == "true" ]]; then
            IFS=',' read -r -a bounds <<< "$LAT_BAND_BOUNDS"
            echo "    use_custom_lat_band_bounds = .true.,"
            echo "    lat_band_bounds_count = ${#bounds[@]},"
            printf '    lat_band_bounds = '
            for ((i=0; i<${#bounds[@]}; i++)); do
                if [[ $i -gt 0 ]]; then
                    printf ', '
                fi
                printf '%s' "${bounds[$i]}"
            done
            printf ',\n'
        else
            echo "    use_custom_lat_band_bounds = .false.,"
            echo "    lat_band_bounds_count = 0,"
        fi

        cat << EOF
/
EOF
    } > "$config_file"
}

process_list() {
    # Iterate one date per line (blank/comment lines ignored) and run serially.
    local list_file="$1"
    local source_root="$2"
    local date source_dir config_file date_log out_hist

    while IFS= read -r raw_line; do
        date=$(echo "$raw_line" | tr -d '[:space:]')
        [[ -z "$date" ]] && continue
        [[ "$date" == \#* ]] && continue

        source_dir="$source_root/$date"
        if [[ ! -d "$source_dir" ]]; then
            echo "Warning: source directory missing for date $date in $source_root, skipping" | tee -a "$RUN_LOG"
            continue
        fi

        out_hist="$TARGET_DIR_HIST/${OUTPUT_FILE_PREFIX}_${date}.nc"
        config_file="$NAMELIST_DIR/config_histograms_by_lat_band_${date}.nml"
        make_namelist "$source_dir" "$out_hist" "$config_file"

        date_log="$LOG_DIR/compute_work_async_histograms_by_lat_band_${date}.log"
        echo "[$(date +%F\ %T)] Running histogram-only compute for date $date (log: $date_log)" | tee -a "$RUN_LOG"
        "$COMPUTE_BIN" "$config_file" > "$date_log" 2>&1
        rm -f "$config_file"
        echo "[$(date +%F\ %T)] Finished date $date" | tee -a "$RUN_LOG"
    done < "$list_file"
}

if [[ "$GROUP_BY_PERIOD" == "1" ]]; then
    period_task_dir="$LOG_DIR/compute_work_async_histograms_by_lat_band_period_tasks_$(date +%Y%m%d_%H%M%S)_$$"
    period_keys_file="$period_task_dir/period_keys.txt"
    combined_file="$period_task_dir/dates_with_roots.txt"
    sorted_file="$period_task_dir/dates_with_roots_sorted.txt"

    mkdir -p "$period_task_dir"
    : > "$combined_file"

    while IFS= read -r raw_line; do
        date=$(echo "$raw_line" | tr -d '[:space:]')
        [[ -z "$date" ]] && continue
        [[ "$date" == \#* ]] && continue
        echo "$date|$SOURCE_ROOT_PART1" >> "$combined_file"
    done < "$LIST_FILE_PART1"

    while IFS= read -r raw_line; do
        date=$(echo "$raw_line" | tr -d '[:space:]')
        [[ -z "$date" ]] && continue
        [[ "$date" == \#* ]] && continue
        echo "$date|$SOURCE_ROOT_PART2" >> "$combined_file"
    done < "$LIST_FILE_PART2"

    sort -t'|' -k1,1 "$combined_file" > "$sorted_file"
    : > "$period_keys_file"

    current_key=""
    while IFS='|' read -r date source_root; do
        period_key=$(period_key_from_date "$date")
        if [[ -z "$period_key" ]]; then
            echo "Error: failed to derive period key for date '$date'" >&2
            exit 1
        fi
        if [[ "$period_key" != "$current_key" ]]; then
            echo "$period_key" >> "$period_keys_file"
            current_key="$period_key"
        fi
        echo "$date|$source_root" >> "$period_task_dir/${period_key}.lst"
    done < "$sorted_file"

    while IFS= read -r period_key; do
        [[ -z "$period_key" ]] && continue
        echo "[$(date +%F\ %T)] Processing period $period_key" | tee -a "$RUN_LOG"
        period_list="$period_task_dir/${period_key}.lst"
        first_entry=$(head -n 1 "$period_list")
        if [[ -z "$first_entry" ]]; then
            echo "Warning: period list is empty for period $period_key, skipping" | tee -a "$RUN_LOG"
            continue
        fi
        IFS='|' read -r date source_root <<< "$first_entry"
        source_dir="$source_root/$date"
        if [[ ! -d "$source_dir" ]]; then
            echo "Warning: source directory missing for representative date $date in $source_root, skipping period $period_key" | tee -a "$RUN_LOG"
            continue
        fi

        out_hist="$TARGET_DIR_HIST/${OUTPUT_FILE_PREFIX}_${period_key}.nc"
        config_file="$NAMELIST_DIR/config_histograms_by_lat_band_${period_key}.nml"
        {
            cat << EOF
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
    path_hist_out   = '$out_hist',
    nlat_bands      = $NLAT_BANDS,
    aggregation_mode = '$AGGREGATION_MODE',
    date_list_file  = '$period_list',
    target_period_key = '$period_key',
EOF

            if [[ "$USE_CUSTOM_LAT_BAND_BOUNDS" == "true" ]]; then
                IFS=',' read -r -a bounds <<< "$LAT_BAND_BOUNDS"
                echo "    use_custom_lat_band_bounds = .true.,"
                echo "    lat_band_bounds_count = ${#bounds[@]},"
                printf '    lat_band_bounds = '
                for ((i=0; i<${#bounds[@]}; i++)); do
                    if [[ $i -gt 0 ]]; then
                        printf ', '
                    fi
                    printf '%s' "${bounds[$i]}"
                done
                printf ',\n'
            else
                echo "    use_custom_lat_band_bounds = .false.,"
                echo "    lat_band_bounds_count = 0,"
            fi

            cat << EOF
/
EOF
        } > "$config_file"

        period_log="$LOG_DIR/compute_work_async_histograms_by_lat_band_${period_key}.log"
        echo "[$(date +%F\ %T)] Running one-pass histogram compute for period $period_key (log: $period_log)" | tee -a "$RUN_LOG"
        "$COMPUTE_BIN" "$config_file" > "$period_log" 2>&1
        rm -f "$config_file"
        echo "[$(date +%F\ %T)] Finished period $period_key" | tee -a "$RUN_LOG"
    done < "$period_keys_file"

    echo "[$(date +%F\ %T)] All grouped periods processed serially." | tee -a "$RUN_LOG"
    exit 0
fi

process_list "$LIST_FILE_PART1" "$SOURCE_ROOT_PART1"
process_list "$LIST_FILE_PART2" "$SOURCE_ROOT_PART2"

echo "[$(date +%F\ %T)] All dates processed serially." | tee -a "$RUN_LOG"
