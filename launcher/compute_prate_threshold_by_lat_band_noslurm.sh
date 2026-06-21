#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
RUN_ID=${RUN_ID:-""}
RUN_ID_ROOT=${RUN_ID_ROOT:-""}
RUN_ID_TIMESTAMP=$(date -u +%Y%m%dT%H%M%SZ)
if [[ -n "$RUN_ID" ]]; then
    :
elif [[ -n "$RUN_ID_ROOT" ]]; then
    RUN_ID="$RUN_ID_ROOT"_"$RUN_ID_TIMESTAMP"
else
    RUN_ID="$RUN_ID_TIMESTAMP"
fi
NLAT_BANDS=${NLAT_BANDS:-18}
USE_CUSTOM_LAT_BAND_BOUNDS=${USE_CUSTOM_LAT_BAND_BOUNDS:-false}
LAT_BAND_BOUNDS=${LAT_BAND_BOUNDS:-}
USE_CUSTOM_LAT_BAND_SEGMENTS=${USE_CUSTOM_LAT_BAND_SEGMENTS:-false}
# Segment spec format: "south:north[,south:north...];south:north[...]" (one ';'-separated band spec per band).
# Example for 2 bands: "-60:-30,30:60;-30:30"
LAT_BAND_SEGMENTS=${LAT_BAND_SEGMENTS:-}

# Serial launcher (no Slurm) for compute_prate_thresholds_by_lat_band.
# The Fortran program reads two history roots and two date lists.

# Simulation selector:
#   control (default)
#   warming
# Manual override example:
#   SIMULATION=warming ./compute_prate_threshold_by_lat_band_noslurm.sh
SIMULATION=${SIMULATION:-control}

case "$SIMULATION" in
    control)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp
        LIST_FILE_PART1=$PROJECT_ROOT/launcher/list/list_control_part1.txt
        LIST_FILE_PART2=$PROJECT_ROOT/launcher/list/list_control_part2.txt
        OUTPUT_FILE_BASENAME=thresholds_control_by_lat_band.txt
        ;;
    warming)
        SOURCE_ROOT_PART1=/scratch/cimes/GLOBALFV3/stellar_run/processed/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        SOURCE_ROOT_PART2=/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv/pp
        LIST_FILE_PART1=$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part1.txt
        LIST_FILE_PART2=$PROJECT_ROOT/launcher/list/list_PLUS_4K_CO2_1270ppmv_part2.txt
        OUTPUT_FILE_BASENAME=thresholds_warming_by_lat_band.txt
        ;;
    *)
        echo "Error: unsupported SIMULATION='$SIMULATION'. Use 'control' or 'warming'." >&2
        exit 1
        ;;
esac

if [[ "$USE_CUSTOM_LAT_BAND_BOUNDS" == "true" && "$USE_CUSTOM_LAT_BAND_SEGMENTS" == "true" ]]; then
    echo "Error: only one of USE_CUSTOM_LAT_BAND_BOUNDS or USE_CUSTOM_LAT_BAND_SEGMENTS may be true." >&2
    exit 1
fi

if [[ "$USE_CUSTOM_LAT_BAND_SEGMENTS" == "true" && -z "$LAT_BAND_SEGMENTS" ]]; then
    echo "Error: LAT_BAND_SEGMENTS is required when USE_CUSTOM_LAT_BAND_SEGMENTS=true." >&2
    exit 1
fi

OUTPUT_DIR_BASE=${OUTPUT_DIR_BASE:-$PROJECT_ROOT/output/thresholds}
OUTPUT_FILE=${OUTPUT_FILE:-$OUTPUT_DIR_BASE/$RUN_ID/$OUTPUT_FILE_BASENAME}

# Fixed output location is selected by SIMULATION (no separate output override).
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
NAMELIST_DIR=$PROJECT_ROOT/output/namelists

LOG_DIR=${LOG_DIR:-$PROJECT_ROOT/logs/slurm}
MANIFEST_DIR=${MANIFEST_DIR:-$PROJECT_ROOT/logs/manifests}
MANIFEST_TOOL=${MANIFEST_TOOL:-$PROJECT_ROOT/python/workflow_manifest.py}
mkdir -p "$LOG_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$MANIFEST_DIR"
mkdir -p "$NAMELIST_DIR"

manifest_path="$MANIFEST_DIR/thresholds/manifest_thresholds_${RUN_ID}.json"
cat << EOF | python3 "$MANIFEST_TOOL" --output "$manifest_path" --workflow-type thresholds --run-id "$RUN_ID" --project-root "$PROJECT_ROOT" --launcher "launcher/compute_prate_threshold_by_lat_band_noslurm.sh" --mode "noslurm" --root "$OUTPUT_DIR_BASE"
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
    "log_dir": "$LOG_DIR",
    "output": {
        "thresholds_base": "$OUTPUT_DIR_BASE",
        "thresholds_file": "$OUTPUT_FILE"
    },
    "lat_banding": {
        "enabled": true,
        "nlat_bands": "$NLAT_BANDS",
        "use_custom_bounds": $([[ "$USE_CUSTOM_LAT_BAND_BOUNDS" == "true" ]] && echo true || echo false),
        "custom_bounds": "$LAT_BAND_BOUNDS",
        "use_custom_segments": $([[ "$USE_CUSTOM_LAT_BAND_SEGMENTS" == "true" ]] && echo true || echo false),
        "custom_segments": "$LAT_BAND_SEGMENTS"
    }
}
EOF

module purge || true
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

# Prefer local/bin build output, fallback to legacy bin path.
COMPUTE_PRATE_THRESH_BIN=${COMPUTE_PRATE_THRESH_BIN:-$PROJECT_ROOT/local/bin/compute_prate_thresholds_by_lat_band}
if [[ ! -x "$COMPUTE_PRATE_THRESH_BIN" && "$COMPUTE_PRATE_THRESH_BIN" == "$PROJECT_ROOT/local/bin/compute_prate_thresholds_by_lat_band" ]]; then
    COMPUTE_PRATE_THRESH_BIN=$PROJECT_ROOT/bin/compute_prate_thresholds_by_lat_band
fi

if [[ ! -f "$LIST_FILE_PART1" ]]; then
    echo "Error: list file not found: $LIST_FILE_PART1" >&2
    exit 1
fi

if [[ ! -f "$LIST_FILE_PART2" ]]; then
    echo "Error: list file not found: $LIST_FILE_PART2" >&2
    exit 1
fi

if [[ ! -d "$SOURCE_ROOT_PART1" ]]; then
    echo "Error: source root directory not found: $SOURCE_ROOT_PART1" >&2
    exit 1
fi

if [[ ! -d "$SOURCE_ROOT_PART2" ]]; then
    echo "Error: source root directory not found: $SOURCE_ROOT_PART2" >&2
    exit 1
fi

if [[ ! -x "$COMPUTE_PRATE_THRESH_BIN" ]]; then
    echo "Error: compute_prate_thresholds_by_lat_band binary not found or not executable." >&2
    echo "Tried: $PROJECT_ROOT/local/bin/compute_prate_thresholds_by_lat_band and $PROJECT_ROOT/bin/compute_prate_thresholds_by_lat_band" >&2
    exit 1
fi

CONFIG_FILE="$NAMELIST_DIR/config_prate_threshold_by_lat_band_$$.nml"
RUN_LOG="$LOG_DIR/compute_prate_threshold_by_lat_band_noslurm_$(date +%Y%m%d_%H%M%S).log"

emit_lat_banding_namelist() {
    if [[ "$USE_CUSTOM_LAT_BAND_BOUNDS" == "true" ]]; then
        IFS=',' read -r -a bounds <<< "$LAT_BAND_BOUNDS"
        echo "    use_custom_lat_band_bounds = .true.,"
        echo "    use_custom_lat_band_segments = .false.,"
        echo "    lat_band_bounds_count = ${#bounds[@]},"
        printf '    lat_band_bounds = '
        for ((i=0; i<${#bounds[@]}; i++)); do
            if [[ $i -gt 0 ]]; then
                printf ', '
            fi
            printf '%s' "${bounds[$i]}"
        done
        printf ',\n'
        return
    fi

    if [[ "$USE_CUSTOM_LAT_BAND_SEGMENTS" == "true" ]]; then
        IFS=';' read -r -a band_specs <<< "$LAT_BAND_SEGMENTS"
        if [[ ${#band_specs[@]} -ne $NLAT_BANDS ]]; then
            echo "Error: LAT_BAND_SEGMENTS must provide exactly $NLAT_BANDS ';'-separated band specs." >&2
            exit 1
        fi

        echo "    use_custom_lat_band_bounds = .false.,"
        echo "    lat_band_bounds_count = 0,"
        echo "    use_custom_lat_band_segments = .true.,"

        for ((iband=0; iband<${#band_specs[@]}; iband++)); do
            band_clean=$(echo "${band_specs[$iband]}" | tr -d '[:space:]')
            if [[ -z "$band_clean" ]]; then
                echo "Error: empty segment spec for band $((iband + 1))." >&2
                exit 1
            fi
            IFS=',' read -r -a segments <<< "$band_clean"
            echo "    lat_band_segment_count($((iband + 1))) = ${#segments[@]},"

            for ((iseg=0; iseg<${#segments[@]}; iseg++)); do
                seg="${segments[$iseg]}"
                if [[ "$seg" != *:* ]]; then
                    echo "Error: invalid segment '$seg' in band $((iband + 1)); expected south:north." >&2
                    exit 1
                fi
                south=${seg%%:*}
                north=${seg##*:}
                if [[ -z "$south" || -z "$north" ]]; then
                    echo "Error: invalid segment '$seg' in band $((iband + 1)); south/north must be non-empty." >&2
                    exit 1
                fi
                echo "    lat_band_segment_south($((iseg + 1)),$((iband + 1))) = $south,"
                echo "    lat_band_segment_north($((iseg + 1)),$((iband + 1))) = $north,"
            done
        done
        return
    fi

    echo "    use_custom_lat_band_bounds = .false.,"
    echo "    use_custom_lat_band_segments = .false.,"
    echo "    lat_band_bounds_count = 0,"
}

cat > "$CONFIG_FILE" << EOF
&config
    history_root_part1   = '$SOURCE_ROOT_PART1',
    history_root_part2   = '$SOURCE_ROOT_PART2',
    date_list_file_part1 = '$LIST_FILE_PART1',
    date_list_file_part2 = '$LIST_FILE_PART2',
    output_file          = '$OUTPUT_FILE',
    nlat_bands           = $NLAT_BANDS,
EOF

emit_lat_banding_namelist >> "$CONFIG_FILE"

cat >> "$CONFIG_FILE" << EOF
/
EOF

echo "[$(date +%F\ %T)] Starting compute_prate_thresholds_by_lat_band" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Simulation: $SIMULATION" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Binary: $COMPUTE_PRATE_THRESH_BIN" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Source root part1: $SOURCE_ROOT_PART1" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Source root part2: $SOURCE_ROOT_PART2" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] List file part1: $LIST_FILE_PART1" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] List file part2: $LIST_FILE_PART2" | tee -a "$RUN_LOG"
echo "[$(date +%F\ %T)] Output file: $OUTPUT_FILE" | tee -a "$RUN_LOG"

"$COMPUTE_PRATE_THRESH_BIN" "$CONFIG_FILE" 2>&1 | tee -a "$RUN_LOG"
rm -f "$CONFIG_FILE"

echo "[$(date +%F\ %T)] Done. Thresholds written to: $OUTPUT_FILE" | tee -a "$RUN_LOG"
