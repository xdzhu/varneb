"""Validate and atomically promote converged HfO2 endpoint CONTCAR files."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import tempfile


def _composition_from_vasp(path: Path) -> dict[str, int] | None:
    """Read a VASP5 element/count header when an older summary lacks it."""

    try:
        lines = path.read_text(encoding="utf-8").splitlines()
        if len(lines) < 7:
            return None
        symbols = lines[5].split()
        counts = [int(value) for value in lines[6].split()]
        if len(symbols) != len(counts):
            return None
        return dict(sorted(zip(symbols, counts)))
    except (OSError, ValueError, UnicodeError):
        return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t-workdir", type=Path, required=True)
    parser.add_argument("--po-workdir", type=Path, required=True)
    parser.add_argument("--target-root", type=Path, required=True)
    parser.add_argument("--force-thr", type=float, default=0.02)
    parser.add_argument("--stress-thr-kbar", type=float, default=0.1)
    parser.add_argument("--natoms", type=int, default=12)
    return parser.parse_args()


def _check_endpoint(path: Path, *, force_thr: float, stress_thr: float, natoms: int) -> dict:
    summary_path = path / "relax_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    issues: list[str] = []
    # The ASE endpoint driver writes a summary only after the Python process
    # reaches the final bookkeeping block; older summaries predate the
    # explicit returncode field and are therefore successful by construction.
    if summary.get("returncode", 0) != 0:
        issues.append(f"returncode={summary.get('returncode')!r}")
    if summary.get("converged") is not True:
        issues.append("converged flag is not true")
    if int(summary.get("natoms", -1)) != natoms:
        issues.append(f"natoms={summary.get('natoms')!r}, expected {natoms}")
    final_natoms = summary.get("final_natoms", summary.get("natoms"))
    if int(final_natoms) != natoms:
        issues.append(f"final_natoms={final_natoms!r}, expected {natoms}")
    force = summary.get("max_generalized_force_eV_per_A")
    if force is None:
        force = summary.get("max_force_eV_per_A")
    if force is None or float(force) >= force_thr:
        issues.append(f"force={force!r} is not below {force_thr}")
    stress = summary.get("max_abs_stress_kbar")
    if stress is None:
        stress = summary.get("last_max_stress_kbar")
    if stress is None or float(stress) >= stress_thr:
        issues.append(f"stress={stress!r} is not below {stress_thr} kbar")
    contcar = path / "CONTCAR"
    if not contcar.exists():
        issues.append(f"missing {contcar}")
    if issues:
        raise RuntimeError(f"endpoint {path} failed promotion gate: {'; '.join(issues)}")
    composition = summary.get("final_composition", summary.get("composition"))
    if composition is None:
        composition = _composition_from_vasp(contcar)
    return {
        "workdir": str(path),
        "summary": str(summary_path),
        "contcar": str(contcar),
        "force_eV_per_A": float(force),
        "stress_kbar": float(stress),
        "natoms": natoms,
        "composition": composition,
    }


def _atomic_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=target.parent)
    os.close(fd)
    try:
        shutil.copy2(source, temporary)
        os.replace(temporary, target)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def promote(t_workdir: Path, po_workdir: Path, target_root: Path, *, force_thr: float, stress_thr: float, natoms: int) -> dict:
    t = _check_endpoint(t_workdir.resolve(), force_thr=force_thr, stress_thr=stress_thr, natoms=natoms)
    po = _check_endpoint(po_workdir.resolve(), force_thr=force_thr, stress_thr=stress_thr, natoms=natoms)
    if t["composition"] is not None and po["composition"] is not None:
        if t["composition"] != po["composition"]:
            raise RuntimeError(
                f"endpoint compositions differ: T={t['composition']!r}, PO={po['composition']!r}"
            )
    _atomic_copy(Path(t["contcar"]), target_root / "relaxed_T" / "CONTCAR")
    _atomic_copy(Path(po["contcar"]), target_root / "relaxed_PO" / "CONTCAR")
    report = {
        "status": "promoted",
        "force_threshold_eV_per_A": force_thr,
        "stress_threshold_kbar": stress_thr,
        "endpoints": {"T": t, "PO": po},
    }
    report_path = target_root / "endpoint_promotion_gate.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> int:
    args = parse_args()
    try:
        report = promote(
            args.t_workdir,
            args.po_workdir,
            args.target_root,
            force_thr=args.force_thr,
            stress_thr=args.stress_thr_kbar,
            natoms=args.natoms,
        )
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"[ERROR] {exc}")
        return 1
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
