"""Block QE VCNEB production unless its static cutoff audit matches the run."""

from __future__ import annotations

import argparse
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


def _as_pair(value: Any, label: str) -> tuple[float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f"{label} must be a two-value cutoff pair")
    return float(value[0]), float(value[1])


def validate(
    audit: dict[str, Any],
    endpoint_gate: dict[str, Any],
    *,
    ecutwfc: float,
    ecutrho: float,
    kpts: tuple[int, int, int],
    scf_thr: float,
    git_revision: str,
) -> None:
    """Validate that an accepted static series is exactly the planned run basis."""

    if audit.get("status") != "ok" or int(audit.get("n_static_points", 0)) < 3:
        raise ValueError("QE static convergence audit is not an accepted three-or-more-point series")
    if not git_revision or git_revision == "remote-sync-unknown":
        raise ValueError("production requires an explicit VCNEB_GIT_REVISION, not remote-sync-unknown")
    conditions = audit.get("fixed_conditions")
    if not isinstance(conditions, dict):
        raise ValueError("QE static audit lacks fixed_conditions")
    if conditions.get("git_revision") != git_revision:
        raise ValueError("QE static audit code revision differs from planned production revision")
    parameters = conditions.get("calculator_parameters_except_cutoffs")
    if not isinstance(parameters, dict):
        raise ValueError("QE static audit lacks calculator parameters")
    if tuple(parameters.get("kpts") or ()) != tuple(kpts):
        raise ValueError("QE static audit k-point mesh differs from planned production mesh")
    electrons = ((parameters.get("input_data") or {}).get("electrons") or {})
    if float(electrons.get("conv_thr", "nan")) != float(scf_thr):
        raise ValueError("QE static audit SCF threshold differs from planned production threshold")
    accepted = audit.get("accepted_final_pair")
    if not isinstance(accepted, dict):
        raise ValueError("QE static audit lacks accepted_final_pair")
    highest = _as_pair(accepted.get("higher_cutoffs_Ry"), "accepted highest cutoff")
    if highest != (float(ecutwfc), float(ecutrho)):
        raise ValueError("QE static audit highest cutoff differs from planned production cutoff")
    endpoints = endpoint_gate.get("endpoints")
    initial = endpoints.get("initial") if isinstance(endpoints, dict) else None
    candidate_initial_sha = initial.get("candidate_sha256") if isinstance(initial, dict) else None
    if endpoint_gate.get("matches") is not True or not isinstance(initial, dict) or initial.get("matches") is not True:
        raise ValueError("endpoint identity gate did not accept the initial endpoint")
    if conditions.get("initial_endpoint_sha256") != candidate_initial_sha:
        raise ValueError("QE static audit initial endpoint differs from the accepted QE path endpoint")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, help="passing qe_static_convergence_audit.json")
    parser.add_argument("--endpoint-gate", required=True, help="passing endpoint_identity_gate.json")
    parser.add_argument("--ecutwfc", required=True, type=float)
    parser.add_argument("--ecutrho", required=True, type=float)
    parser.add_argument("--kpts", required=True, type=int, nargs=3)
    parser.add_argument("--scf-thr", required=True, type=float)
    parser.add_argument("--git-revision", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        validate(
            _load(args.audit, "QE static convergence audit"),
            _load(args.endpoint_gate, "endpoint identity gate"),
            ecutwfc=args.ecutwfc,
            ecutrho=args.ecutrho,
            kpts=tuple(args.kpts),
            scf_thr=args.scf_thr,
            git_revision=args.git_revision,
        )
    except ValueError as exc:
        print(f"[BLOCKED] {exc}", file=sys.stderr)
        return 2
    print(f"[OK] QE static convergence gate passed: {Path(args.audit).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
