#!/usr/bin/env python3
"""Write workflow run manifests with a stable envelope and caller payload.

Usage example:
  cat <<'JSON' | python3 workflow_manifest.py \
    --output logs/manifests/manifest_work_20260601T010203Z.json \
    --workflow-type work \
    --run-id 20260601T010203Z \
    --project-root /path/to/repo \
    --launcher launcher/compute_work_async_array.sh \
    --mode slurm_array \
    --root /scratch/results/work
  {"simulation": "control", "thresholding_enabled": false}
  JSON
"""

import argparse
import datetime as dt
import json
import os
import pathlib
import sys
import tempfile
from typing import Any, Dict


def _now_utc_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def _read_payload_from_stdin() -> Dict[str, Any]:
    if sys.stdin.isatty():
        return {}
    text = sys.stdin.read().strip()
    if not text:
        return {}
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"Invalid JSON payload on stdin: {exc}") from exc
    if not isinstance(parsed, dict):
        raise SystemExit("JSON payload on stdin must be an object")
    return parsed


def _atomic_write_json(path: pathlib.Path, obj: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2, sort_keys=True)
            fh.write("\n")
        os.replace(tmp_name, path)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)


def main() -> int:
    parser = argparse.ArgumentParser(description="Write workflow run manifest JSON")
    parser.add_argument("--output", required=True, help="Manifest output file path")
    parser.add_argument(
        "--workflow-type",
        required=True,
        choices=["work", "histograms", "thresholds", "percentile-thresholds"],
        help="Workflow family",
    )
    parser.add_argument("--run-id", required=True, help="Run identifier")
    parser.add_argument("--project-root", required=True, help="Project root path")
    parser.add_argument("--launcher", required=True, help="Launcher script path/name")
    parser.add_argument("--mode", required=True, help="Execution mode, e.g. slurm_array/noslurm")
    parser.add_argument("--root", default="", help="Primary output root/base path")

    args = parser.parse_args()

    payload = _read_payload_from_stdin()
    manifest: Dict[str, Any] = {
        "schema_version": 1,
        "created_utc": _now_utc_iso(),
        "workflow_type": args.workflow_type,
        "project_root": args.project_root,
        "run_id": args.run_id,
        "launcher": args.launcher,
        "mode": args.mode,
        "root": args.root,
    }
    manifest.update(payload)

    _atomic_write_json(pathlib.Path(args.output), manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
