"""Audit native VASP pressure refinement of the GaN B1 basin candidate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_ts_basin_followups import residual_stress_kbar
from scripts.audit_gan_ts_basin_pilot import energy_triples, ga_n_coordination
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256


def audit(work_root: Path, source_root: Path, source_audit_path: Path,
          structural_match_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    source_audit = json.loads(source_audit_path.read_text(encoding="utf-8"))
    if (manifest.get("purpose")
            != "GaN_1000eV_B1_candidate_pressure_refinement_not_TS_certificate"
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("PSTRESS_kbar") != 457.0
            or manifest.get("source_sha256", {}).get("source_manifest")
            != sha256(source_root / "manifest.json")
            or manifest.get("source_sha256", {}).get("source_audit")
            != sha256(source_audit_path)
            or manifest.get("source_sha256", {}).get("structural_match")
            != sha256(structural_match_path)
            or manifest.get("source_sha256", {}).get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_ts_b1_pressure_refine.py"))):
        raise ValueError("B1 pressure-refinement provenance changed")
    source = source_root / "cases/negative_mode_minus"
    if (source_audit.get("case_kind") != "continued_minus"
            or sha256(source / "OUTCAR") != manifest["source_sha256"]["source_OUTCAR"]
            or sha256(source / "CONTCAR") != manifest["source_sha256"]["source_CONTCAR"]):
        raise ValueError("source minus basin geometry/output changed")
    case_dir = work_root / "case"
    if (any(sha256(case_dir / filename) != expected
            for filename, expected in manifest["input_sha256"].items())
            or b"\r" in (case_dir / "INCAR").read_bytes()
            or not same_geometry(read(case_dir / "POSCAR", format="vasp"),
                                 read(source / "CONTCAR", format="vasp"), tolerance=2e-5)):
        raise ValueError("B1 refinement input changed or starts from wrong cell")
    outcar = case_dir / "OUTCAR"
    raw = outcar.read_text(encoding="utf-8", errors="replace")
    if ("General timing and accounting informations" not in raw
            or "aborting loop because EDIFF is reached" not in raw
            or "PSTRESS=  457.0" not in raw):
        raise ValueError("B1 refinement VASP output incomplete or pressure changed")
    triples = energy_triples(raw)
    if not 1 <= len(triples) <= 30:
        raise ValueError("unexpected ionic energy record count")
    atoms = read(outcar)
    contcar = case_dir / "CONTCAR"
    candidate = read(contcar, format="vasp")
    fmax = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    stress = residual_stress_kbar(atoms.get_stress(voigt=False), 45.7)
    pressure = 45.7 * GPa
    if (len(atoms) != 4 or not np.isfinite(atoms.get_forces()).all()
            or not np.isfinite(atoms.get_stress(voigt=False)).all()
            or not np.isclose(triples[-1][2], pressure * atoms.get_volume(),
                              atol=2e-4, rtol=0)):
        raise ValueError("last evaluated force, stress, volume or P*V invalid")
    force_pass = fmax <= 0.02
    stress_pass = stress <= 1.0
    result = {
        "status": ("GaN_B1_pressure_refine_raw_audited_force_stress_pass_not_TS_certificate"
                   if force_pass and stress_pass else
                   "GaN_B1_pressure_refine_raw_audited_force_or_stress_gate_failed"),
        "n_ionic_energy_records": len(triples),
        "initial_enthalpy_minus_source_final_eV_per_cell": (
            triples[0][1] - source_audit["last_enthalpy_eV_per_cell"]
        ),
        "last_free_energy_eV_per_cell": triples[-1][0],
        "last_PV_eV_per_cell": triples[-1][2],
        "last_enthalpy_eV_per_cell": triples[-1][1],
        "last_enthalpy_minus_first_eV_per_cell": triples[-1][1] - triples[0][1],
        "last_volume_A3": float(atoms.get_volume()),
        "last_max_atomic_force_eV_per_A": fmax,
        "last_max_stress_residual_from_45p7GPa_kbar": stress,
        "last_coordination_GaN_2p4A": ga_n_coordination(atoms),
        "contcar_is_last_evaluated_geometry": same_geometry(atoms, candidate, tolerance=2e-5),
        "vasp_reports_ionic_convergence": (
            "reached required accuracy - stopping structural energy minimisation" in raw
        ),
        "force_pass_0p02eV_per_A": force_pass,
        "stress_pass_1kbar": stress_pass,
        "source_sha256": {"manifest": sha256(manifest_path),
                          "source_audit": sha256(source_audit_path),
                          "structural_match": sha256(structural_match_path),
                          "OUTCAR": sha256(outcar), "CONTCAR": sha256(contcar),
                          "auditor": sha256(Path(__file__))},
        "limitations": [
            "Passing force and stress gates confirms a B1-like basin candidate, not both TS basin links.",
            "Structural identity and same-setting static E/force/stress checks remain separate gates.",
            "The VASP enthalpy field already includes P*V; do not add pressure work twice.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "H_eV": triples[-1][1],
                      "fmax": fmax, "stress_residual_kbar": stress}))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work-root", "source-root", "source-audit", "structural-match", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.work_root, args.source_root, args.source_audit,
          args.structural_match, args.output)


if __name__ == "__main__":
    main()
