"""Audit four GaN atomic-transverse statics with the path cell held fixed."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_600eV_atomic_tube_canary import ANCHORS
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work: Path, trajectory: Path, feasibility_path: Path,
          hessian_npz: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_fixed_cell_atomic_transverse_canary"
            or manifest.get("anchors") != list(ANCHORS)
            or manifest.get("q_atom_A") != [-0.2, 0.2]
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or len(manifest.get("cases", [])) != 4
            or manifest["source_sha256"].get("trajectory") != sha256(trajectory)
            or manifest["source_sha256"].get("atomic_feasibility") != sha256(feasibility_path)
            or manifest["source_sha256"].get("hessian_npz") != sha256(hessian_npz)
            or manifest["source_sha256"].get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_atomic_tube_canary.py"))
            or feasibility.get("status")
            != "GaN_600eV_atomic_transverse_tube_geometry_only_no_DFT"
            or feasibility["source_sha256"].get("trajectory") != sha256(trajectory)):
        raise ValueError("atomic-tube geometry or 600-eV electronic provenance changed")
    frames = read(trajectory, index=":")
    pressure = 45.7 * GPa
    by_anchor: dict[int, dict[str, dict]] = {index: {} for index in ANCHORS}
    cases = []
    failures = []
    for record in manifest["cases"]:
        directory = work / "cases" / record["name"]
        if any(record.get("input_sha256", {}).get(name) != digest
               or sha256(directory / name) != digest
               for name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"{record['name']} changed electronic inputs")
        stdout_path = directory / "vasp.stdout"
        stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
        if "Inconsistent Bravais lattice types found" in stdout:
            failures.append({
                "case": record["name"],
                "failure": "VASP_direct_reciprocal_Bravais_classification_conflict",
                "vasp_stdout_sha256": sha256(stdout_path),
                "input_sha256": record["input_sha256"],
            })
            continue
        atoms, detail = completed_case(
            directory, {**record, "POSCAR_sha256": record["input_sha256"]["POSCAR"]}
        )
        index = int(record["image_index"])
        if not np.allclose(atoms.cell.array, frames[index].cell.array, atol=2e-5, rtol=0):
            raise ValueError(f"atomic-only cell changed in VASP at {record['name']}")
        q = float(record["q_atom_A"])
        side = "plus" if q > 0 else "minus"
        h0 = float(frames[index].get_potential_energy() + pressure * frames[index].get_volume())
        h = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
        detail.update({
            "image_index": index,
            "q_atom_A": q,
            "enthalpy_eV_per_cell": h,
            "delta_enthalpy_meV_per_GaN": float((h - h0) * 500),
            "vasp_stdout_sha256": sha256(stdout_path),
        })
        if side in by_anchor[index]:
            raise ValueError("duplicate atomic-canary side")
        by_anchor[index][side] = detail
        cases.append(detail)
    if len(cases) + len(failures) != 4:
        raise ValueError("atomic-canary outcomes incomplete")
    return {
        "status": (
            "GaN_600eV_atomic_transverse_canary_all_raw_audited"
            if not failures else "GaN_600eV_atomic_transverse_canary_partial_failure_raw_audited"
        ),
        "claim_limit": "Two path anchors at q_atom=+/-0.2 A; cells fixed to archived VCNEB images; not a global surface",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "anchors": list(ANCHORS),
        "n_completed_static_cases": len(cases),
        "n_failed_cases": len(failures),
        "cases": cases,
        "failures": failures,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "trajectory": sha256(trajectory),
            "atomic_feasibility": sha256(feasibility_path),
            "hessian_npz": sha256(hessian_npz),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "Atomic-only transverse mode differs from the joint atom-strain transverse mode.",
            "Two anchors do not support smooth full-path interpolation or a unique-MEP claim.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "trajectory", "feasibility", "hessian-npz", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.work, args.trajectory, args.feasibility, args.hessian_npz)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"],
                      "completed": report["n_completed_static_cases"],
                      "failed": report["n_failed_cases"]}))


if __name__ == "__main__":
    main()
