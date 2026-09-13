"""Audit an ABACUS-native endpoint relaxation summary without rerunning DFT."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def audit(summary_path: Path, *, force_thr: float | None = None, stress_thr: float | None = None,
          natoms: int | None = None, composition: str | None = None) -> dict:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if summary.get("returncode") != 0:
        issues.append(f"returncode={summary.get('returncode')!r}")
    if summary.get("converged") is not True:
        issues.append("native summary is not converged")
    final_natoms = summary.get("final_natoms", summary.get("natoms"))
    if natoms is not None and final_natoms != natoms:
        issues.append(f"final_natoms={final_natoms!r} != expected {natoms}")
    if composition is not None:
        expected = dict(item.split("=", 1) for item in composition.split(",") if "=" in item)
        got = summary.get("final_composition", summary.get("composition", {}))
        got = {str(key): int(value) for key, value in got.items()}
        expected = {str(key): int(value) for key, value in expected.items()}
        if got != expected:
            issues.append(f"composition={got!r} != expected {expected!r}")
    force = summary.get("last_max_atom_force_eV_per_A")
    stress = summary.get("last_max_stress_kbar")
    if force is None:
        issues.append("missing last_max_atom_force_eV_per_A")
    elif force_thr is not None and float(force) > force_thr:
        issues.append(f"max_atom_force={force} > {force_thr}")
    if stress is None:
        issues.append("missing last_max_stress_kbar")
    elif stress_thr is not None and float(stress) > stress_thr:
        issues.append(f"max_stress={stress} > {stress_thr}")
    return {
        "status": "ok" if not issues else "failed",
        "summary": str(summary_path.resolve()),
        "issues": issues,
        "returncode": summary.get("returncode"),
        "converged": summary.get("converged"),
        "final_natoms": final_natoms,
        "final_composition": summary.get("final_composition", summary.get("composition")),
        "last_max_atom_force_eV_per_A": force,
        "last_max_stress_kbar": stress,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, help="native_relax_summary.json")
    parser.add_argument("--force-thr", type=float, default=0.02)
    parser.add_argument("--stress-thr", type=float, default=0.1)
    parser.add_argument("--natoms", type=int, default=None)
    parser.add_argument("--composition", default=None, help="e.g. Hf=4,O=8")
    args = parser.parse_args()
    report = audit(args.summary, force_thr=args.force_thr, stress_thr=args.stress_thr,
                   natoms=args.natoms, composition=args.composition)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    raise SystemExit(0 if report["status"] == "ok" else 1)


if __name__ == "__main__":
    main()
