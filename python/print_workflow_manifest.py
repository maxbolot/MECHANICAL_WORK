#!/usr/bin/env python3
"""Pretty-print a workflow manifest JSON by workflow type and run id.

Examples:
  python3 python/print_workflow_manifest.py --workflow-type work --run-id 20260601T010203Z
  python3 python/print_workflow_manifest.py --workflow-type histograms --run-id 20260601T010203Z
"""

import argparse
import json
import pathlib
import sys
from typing import Dict


def _normalize_workflow_type(value: str) -> str:
    normalized = value.strip().lower()
    aliases: Dict[str, str] = {
        "work": "work",
        "hist": "histograms",
        "histogram": "histograms",
        "histograms": "histograms",
        "threshold": "thresholds",
        "thresholds": "thresholds",
        "percentile-threshold": "percentile-thresholds",
        "percentile-thresholds": "percentile-thresholds",
    }
    if normalized not in aliases:
        valid = ", ".join(sorted(set(aliases.values())))
        raise SystemExit(
            f"Unsupported workflow type '{value}'. Supported values map to: {valid}"
        )
    return aliases[normalized]


def _default_manifest_dir() -> pathlib.Path:
    # Script lives in <repo>/python, so parent of parent is repo root.
    project_root = pathlib.Path(__file__).resolve().parent.parent
    return project_root / "logs" / "manifests"


def _resolve_manifest_path(workflow_type: str, manifest_dir: pathlib.Path) -> pathlib.Path:
    filename = f"manifests_{workflow_type}.json"
    return manifest_dir / filename


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Pretty-print a workflow manifest by workflow type and run id"
    )
    parser.add_argument(
        "--workflow-type",
        required=True,
        help="Workflow type (work, histograms, thresholds; aliases accepted)",
    )
    parser.add_argument("--run-id", required=True, help="Run identifier")
    parser.add_argument(
        "--manifest-dir",
        default=str(_default_manifest_dir()),
        help="Directory containing manifest JSON files (default: <repo>/logs/manifests)",
    )

    args = parser.parse_args()

    workflow_type = _normalize_workflow_type(args.workflow_type)
    manifest_dir = pathlib.Path(args.manifest_dir).expanduser().resolve()
    manifest_path = _resolve_manifest_path(workflow_type, manifest_dir)

    if not manifest_path.exists():
        raise SystemExit(
            f"Concatenated manifest file not found: {manifest_path}\n"
            "Run python/concat_workflow_manifests.py first to build grouped manifest files."
        )

    try:
        with manifest_path.open("r", encoding="utf-8") as fh:
            payload = json.load(fh)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"Invalid JSON in manifest file '{manifest_path}': {exc}") from exc

    if not isinstance(payload, dict):
        raise SystemExit(f"Expected object payload in '{manifest_path}'")

    manifests = payload.get("manifests")
    if not isinstance(manifests, list):
        raise SystemExit(f"Expected key 'manifests' as a list in '{manifest_path}'")

    matches = [entry for entry in manifests if isinstance(entry, dict) and entry.get("run_id") == args.run_id]
    if not matches:
        raise SystemExit(
            f"Run id '{args.run_id}' not found in '{manifest_path}'.\n"
            "Run python/concat_workflow_manifests.py if new per-run manifests were recently created."
        )

    # If duplicates exist, print the most recently created entry.
    matches.sort(key=lambda m: str(m.get("created_utc", "")))
    payload = matches[-1]

    json.dump(payload, sys.stdout, indent=2, sort_keys=False)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
