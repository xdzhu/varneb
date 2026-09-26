"""Audit same-setting static VASP E+P*V for a relaxed GaN basin."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_ts_basin_followups import residual_stress_kbar
from scripts.audit_gan_ts_basin_pilot import ga_n_coordination
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256


def audit(work_root: Path, relax_root: Path, relax_audit_path: Path,
          structural_match_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    relax_audit = json.loads(relax_audit_path.read_text(encoding="utf-8"))
    match = json.loads(structural_match_path.read_text(encoding="utf-8"))
    if (manifest.get("purpose")
            != "GaN_1000eV_audited_basin_static_E_plus_PV_not_TS_certificate"
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("pressure_for_postprocessing_GPa") != 45.7
            or manifest.get("endpoint") != match.get("endpoint")
            or manifest.get("source_sha256", {}).get("relax_manifest")
            != sha256(relax_root / "manifest.json")
            or manifest.get("source_sha256", {}).get("relax_audit")
            != sha256(relax_audit_path)
            or manifest.get("source_sha256", {}).get("structural_match")
            != sha256(structural_match_path)
            or manifest.get("source_sha256", {}).get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_ts_basin_static.py"))):
        raise ValueError("static basin provenance changed")
    relaxed = relax_root / "case"
    if (not relax_audit["force_pass_0p02eV_per_A"]
            or not relax_audit["stress_pass_1kbar"]
            or sha256(relaxed / "CONTCAR") != manifest["source_sha256"]["relaxed_CONTCAR"]):
        raise ValueError("static basin source is not force/stress-audited")
    case_dir = work_root / "case"
    if (any(sha256(case_dir / filename) != expected
            for filename, expected in manifest["input_sha256"].items())
            or b"PSTRESS" in (case_dir / "INCAR").read_bytes()
            or b"\r" in (case_dir / "INCAR").read_bytes()
            or not same_geometry(read(case_dir / "POSCAR", format="vasp"),
                                 read(relaxed / "CONTCAR", format="vasp"), tolerance=2e-5)):
        raise ValueError("static VASP input contract or geometry changed")
    outcar = case_dir / "OUTCAR"
    raw = outcar.read_text(encoding="utf-8", errors="replace")
    if ("General timing and accounting informations" not in raw
            or "aborting loop because EDIFF is reached" not in raw
            or "free  energy   TOTEN" not in raw):
        raise ValueError("static VASP raw output incomplete")
    atoms = read(outcar)
    if (len(atoms) != 4 or not same_geometry(atoms,
                                             read(case_dir / "POSCAR", format="vasp"),
                                             tolerance=2e-5)
            or not np.isfinite(atoms.get_forces()).all()
            or not np.isfinite(atoms.get_stress(voigt=False)).all()):
        raise ValueError("static VASP geometry, forces or stress invalid")
    E0 = float(atoms.get_potential_energy())
    free_E = float(atoms.get_potential_energy(force_consistent=True))
    PV = float(45.7 * GPa * atoms.get_volume())
    H = float(E0 + PV)
    fmax = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    stress_residual = residual_stress_kbar(atoms.get_stress(voigt=False), 45.7)
    H_delta = H - relax_audit["last_enthalpy_eV_per_cell"]
    closure = abs(H_delta) <= 0.002
    result = {
        "status": ("GaN_basin_static_1000eV_raw_audited_force_stress_and_H_match"
                   if closure and fmax <= 0.02 and stress_residual <= 1.0 else
                   "GaN_basin_static_1000eV_raw_audited_review_required"),
        "endpoint_candidate": manifest["endpoint"],
        "E0_eV_per_cell": E0,
        "free_E_eV_per_cell": free_E,
        "PV_eV_per_cell": PV,
        "H_E0_plus_PV_eV_per_cell": H,
        "H_minus_relax_H_eV_per_cell": H_delta,
        "H_agrees_with_relax_within_2meV": closure,
        "volume_A3": float(atoms.get_volume()),
        "maximum_atomic_force_eV_per_A": fmax,
        "maximum_stress_residual_from_45p7GPa_kbar": stress_residual,
        "coordination_GaN_2p4A": ga_n_coordination(atoms),
        "source_sha256": {"manifest": sha256(manifest_path),
                          "relax_audit": sha256(relax_audit_path),
                          "structural_match": sha256(structural_match_path),
                          "OUTCAR": sha256(outcar),
                          "auditor": sha256(Path(__file__))},
        "limitations": [
            "This static energy is at 1000 eV and cannot be subtracted from original 600 eV path energies.",
            "B4/B1 basin checks require both branches; one static does not establish a transition-state connection.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "H_eV": H,
                      "H_difference_eV": H_delta, "fmax": fmax,
                      "stress_residual_kbar": stress_residual}))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work-root", "relax-root", "relax-audit", "structural-match", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.work_root, args.relax_root, args.relax_audit,
          args.structural_match, args.output)


if __name__ == "__main__":
    main()
