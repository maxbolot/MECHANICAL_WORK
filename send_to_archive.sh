#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT="$SCRIPT_DIR"
ARCHIVE_ROOT="$PROJECT_ROOT/archive"
DRY_RUN="${DRY_RUN:-0}"
ARCHIVE_BATCH_MODE="${ARCHIVE_BATCH_MODE:-timestamp}"
ARCHIVE_BATCH_ID="${ARCHIVE_BATCH_ID:-}"

usage() {
  cat <<'EOF'
Usage:
  DRY_RUN=1 ./send_to_archive.sh <file1> [file2 ...]
  ./send_to_archive.sh <file1> [file2 ...]
  ARCHIVE_BATCH_MODE=flat ./send_to_archive.sh <file1> [file2 ...]

Description:
  Moves one or more files into $PROJECT_ROOT/archive while preserving each
  file's project-relative layout.

  Logging:
    Writes a per-batch manifest at DEST_ROOT/MANIFEST.txt.
    Appends all archive actions to a global index at archive/_index.tsv.

  Batch modes:
    ARCHIVE_BATCH_MODE=timestamp (default)
      archive/<YYYYmmdd_HHMMSS>/<project-relative-path>
    ARCHIVE_BATCH_MODE=flat
      archive/<project-relative-path>

  Optional:
    ARCHIVE_BATCH_ID=<id>
      Overrides auto-generated timestamp when ARCHIVE_BATCH_MODE=timestamp.

Example:
  ./send_to_archive.sh local/src/program.f90

  Moves:
    $PROJECT_ROOT/local/src/program.f90
  To:
    $PROJECT_ROOT/archive/local/src/program.f90
EOF
}

is_truthy() {
  case "$1" in
    1|true|TRUE|yes|YES|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

if [[ "$#" -lt 1 ]]; then
  usage
  exit 1
fi

case "$ARCHIVE_BATCH_MODE" in
  timestamp|flat) ;;
  *)
    echo "Error: unsupported ARCHIVE_BATCH_MODE='$ARCHIVE_BATCH_MODE' (use 'timestamp' or 'flat')." >&2
    exit 1
    ;;
esac

if [[ "$ARCHIVE_BATCH_MODE" == "timestamp" ]]; then
  if [[ -z "$ARCHIVE_BATCH_ID" ]]; then
    ARCHIVE_BATCH_ID=$(date +%Y%m%d_%H%M%S)
  fi
  DEST_ROOT="$ARCHIVE_ROOT/$ARCHIVE_BATCH_ID"
else
  DEST_ROOT="$ARCHIVE_ROOT"
fi

MANIFEST_FILE="$DEST_ROOT/MANIFEST.txt"
GLOBAL_INDEX_FILE="$ARCHIVE_ROOT/_index.tsv"

if is_truthy "$DRY_RUN"; then
  dry_run_enabled=1
else
  dry_run_enabled=0
fi

if [[ "$dry_run_enabled" -eq 1 ]]; then
  echo "DRY_RUN enabled: no files will be moved."
  echo "Archive mode: $ARCHIVE_BATCH_MODE"
  if [[ "$ARCHIVE_BATCH_MODE" == "timestamp" ]]; then
    echo "Archive batch id: $ARCHIVE_BATCH_ID"
  fi
  if [[ ! -d "$DEST_ROOT" ]]; then
    echo "Would create directory: $DEST_ROOT"
  fi
  if [[ ! -f "$MANIFEST_FILE" ]]; then
    echo "Would create manifest: $MANIFEST_FILE"
  fi
  if [[ ! -f "$GLOBAL_INDEX_FILE" ]]; then
    echo "Would create global index: $GLOBAL_INDEX_FILE"
  fi
else
  mkdir -p "$ARCHIVE_ROOT"
  mkdir -p "$DEST_ROOT"

  if [[ ! -f "$MANIFEST_FILE" ]]; then
    printf 'archived_at_utc\tbatch_mode\tbatch_id\tsource_rel\tdestination_rel\n' > "$MANIFEST_FILE"
  fi

  if [[ ! -f "$GLOBAL_INDEX_FILE" ]]; then
    printf 'archived_at_utc\tbatch_mode\tbatch_id\tsource_rel\tdestination_rel\n' > "$GLOBAL_INDEX_FILE"
  fi
fi

for input_path in "$@"; do
  if [[ ! -e "$input_path" ]]; then
    echo "Error: file does not exist: $input_path" >&2
    exit 1
  fi

  src_abs=$(realpath "$input_path")

  if [[ -d "$src_abs" ]]; then
    echo "Error: directories are not supported, expected files only: $input_path" >&2
    exit 1
  fi

  case "$src_abs" in
    "$PROJECT_ROOT"/*) ;;
    *)
      echo "Error: path is outside project root: $input_path" >&2
      exit 1
      ;;
  esac

  rel_path="${src_abs#"$PROJECT_ROOT"/}"

  if [[ "$rel_path" == archive/* ]]; then
    echo "Skipping already archived file: $rel_path"
    continue
  fi

  dst_abs="$DEST_ROOT/$rel_path"
  dst_dir=$(dirname "$dst_abs")
  dst_rel="${dst_abs#"$PROJECT_ROOT"/}"
  archived_at_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  batch_id_for_log="${ARCHIVE_BATCH_ID:-flat}"

  if [[ -e "$dst_abs" ]]; then
    echo "Error: destination already exists: $dst_abs" >&2
    exit 1
  fi

  if [[ "$dry_run_enabled" -eq 1 ]]; then
    echo "Would create directory: $dst_dir"
    echo "Would move: $src_abs -> $dst_abs"
    echo "Would append manifest row: ${archived_at_utc}\t${ARCHIVE_BATCH_MODE}\t${batch_id_for_log}\t${rel_path}\t${dst_rel}"
    echo "Would append global index row: ${archived_at_utc}\t${ARCHIVE_BATCH_MODE}\t${batch_id_for_log}\t${rel_path}\t${dst_rel}"
  else
    mkdir -p "$dst_dir"
    mv "$src_abs" "$dst_abs"
    printf '%s\t%s\t%s\t%s\t%s\n' "$archived_at_utc" "$ARCHIVE_BATCH_MODE" "$batch_id_for_log" "$rel_path" "$dst_rel" >> "$MANIFEST_FILE"
    printf '%s\t%s\t%s\t%s\t%s\n' "$archived_at_utc" "$ARCHIVE_BATCH_MODE" "$batch_id_for_log" "$rel_path" "$dst_rel" >> "$GLOBAL_INDEX_FILE"
    echo "Moved: $src_abs -> $dst_abs"
  fi
done
