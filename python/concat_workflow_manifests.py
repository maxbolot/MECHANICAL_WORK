#!/usr/bin/env python3
"""Concatenate per-run workflow manifests into grouped JSON files.

Writes one output file per workflow under logs/manifests:
  manifests_work.json
  manifests_histograms.json
  manifests_thresholds.json
  manifests_percentile-thresholds.json
"""

import argparse
import json
import pathlib
from typing import Dict, List


WORKFLOW_TYPES = [
    "work",
    "histograms",
    "thresholds",
    "percentile-thresholds",
]


def _default_manifest_dir() -> pathlib.Path:
    project_root = pathlib.Path(__file__).resolve().parent.parent
    return project_root / "logs" / "manifests"


def _load_json(path: pathlib.Path) -> Dict:
    with path.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)
    if not isinstance(payload, dict):
        raise SystemExit(f"Manifest must be a JSON object: {path}")
    return payload


def _sort_key(item: Dict) -> str:
    created = str(item.get("created_utc", ""))
    run_id = str(item.get("run_id", ""))
    return f"{created}|{run_id}"


def _write_atomic(path: pathlib.Path, payload: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=False)
        fh.write("\n")
    tmp_path.replace(path)


def _concat_workflow(manifest_dir: pathlib.Path, workflow_type: str) -> pathlib.Path:
    per_run_dir = manifest_dir / workflow_type
    pattern = f"manifest_{workflow_type}_*.json"

    manifests: List[Dict] = []
    if per_run_dir.exists():
        for file_path in sorted(per_run_dir.glob(pattern)):
            manifests.append(_load_json(file_path))

    manifests.sort(key=_sort_key)

    output_payload = {
        "workflow_type": workflow_type,
        "count": len(manifests),
        "manifests": manifests,
    }
    output_path = manifest_dir / f"manifests_{workflow_type}.json"
    _write_atomic(output_path, output_payload)
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Concatenate per-run manifest JSON files by workflow type"
    )
    parser.add_argument(
        "--manifest-dir",
        default=str(_default_manifest_dir()),
        help="Base manifest directory (default: <repo>/logs/manifests)",
    )
    args = parser.parse_args()

    manifest_dir = pathlib.Path(args.manifest_dir).expanduser().resolve()
    manifest_dir.mkdir(parents=True, exist_ok=True)

    for workflow_type in WORKFLOW_TYPES:
        output_path = _concat_workflow(manifest_dir, workflow_type)
        print(f"wrote {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
