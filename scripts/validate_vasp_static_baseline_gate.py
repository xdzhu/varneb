"""Block VASP VCNEB production unless its fixed-endpoint static baseline matches."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


def _load(path: str, label: str) -> dict[str, Any]:
    source = Path(path).resolve()
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read {label} {source}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{label} {source} is not a JSON object")
    return payload


def _fingerprints(directory: str | Path) -> dict[str, str]:
    root = Path(directory).resolve()
    result = {}
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        path = root / name
        if not path.is_file():
            raise ValueError(f"missing required VASP input: {path}")
        result[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def validate(summary: dict[str, Any], endpoint_gate: dict[str, Any], *, initial_dir: str | Path, git_revision: str) -> None:
    """Verify a completed VASP static SCF used the exact planned input basis."""

    if summary.get("status") != "completed" or summary.get("execution_mode") != "fixed_initial_endpoint_static_scf":
        raise ValueError("VASP static baseline is not a completed fixed-initial-endpoint SCF")
    if summary.get("n_images") != 7 or summary.get("n_interior_images") != 5 or summary.get("evaluated_image_index") != 0:
        raise ValueError("VASP static baseline is not the required 7-total-image endpoint-00 layout")
    if not git_revision or git_revision == "remote-sync-unknown":
        raise ValueError("production requires an explicit VCNEB_GIT_REVISION, not remote-sync-unknown")
    if summary.get("git_revision") != git_revision:
        raise ValueError("VASP static baseline code revision differs from planned production revision")
    endpoints = endpoint_gate.get("endpoints")
    initial = endpoints.get("initial") if isinstance(endpoints, dict) else None
    if endpoint_gate.get("matches") is not True or not isinstance(initial, dict) or initial.get("matches") is not True:
        raise ValueError("endpoint identity gate did not accept the initial endpoint")
    baseline_initial = ((summary.get("endpoint_structures") or {}).get("initial") or {}).get("sha256")
    if baseline_initial != initial.get("candidate_sha256"):
        raise ValueError("VASP static baseline initial endpoint differs from the accepted VASP path endpoint")
    recorded = summary.get("licensed_input_fingerprints")
    if not isinstance(recorded, dict):
        raise ValueError("VASP static baseline lacks licensed-input fingerprints")
    for name, digest in _fingerprints(initial_dir).items():
        item = recorded.get(name)
        if not isinstance(item, dict) or item.get("sha256") != digest:
            raise ValueError(f"VASP {name} differs from the static baseline")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, help="completed vasp_static_summary.json")
    parser.add_argument("--endpoint-gate", required=True, help="passing endpoint_identity_gate.json")
    parser.add_argument("--initial-dir", required=True, help="production VASP initial input directory")
    parser.add_argument("--git-revision", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        validate(
            _load(args.summary, "VASP static baseline"),
            _load(args.endpoint_gate, "endpoint identity gate"),
            initial_dir=args.initial_dir,
            git_revision=args.git_revision,
        )
    except ValueError as exc:
        print(f"[BLOCKED] {exc}", file=sys.stderr)
        return 2
    print(f"[OK] VASP static baseline gate passed: {Path(args.summary).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
