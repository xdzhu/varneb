"""Audit raw VASP GaN basin continuations at 45.7 GPa.

Supports an uninterrupted plus branch from the original seed and a minus
branch restarted from the last *evaluated* ten-step CONTCAR. Even a converged
branch is only a basin candidate until the endpoint identity is checked.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_ts_basin_pilot import energy_triples, ga_n_coordination
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256


PURPOSES = {
    "uninterrupted_plus": "GaN_plus_basin_100step_uninterrupted_from_audited_seed_not_endpoint_certificate",
    "continued_minus": "GaN_1000eV_basin_continuation_from_audited_pilot_not_phase_certificate",
}


def residual_stress_kbar(stress_eV_per_A3: np.ndarray, pressure_GPa: float) -> float:
    stress = np.asarray(stress_eV_per_A3, dtype=float)
    if stress.shape != (3, 3) or not np.isfinite(stress).all():
        raise ValueError("invalid VASP stress tensor")
    residual_GPa = stress / GPa + pressure_GPa * np.eye(3)
    return float(np.max(np.abs(residual_GPa)) * 10.0)


def audit(case_kind: str, work_root: Path, pilot_root: Path,
          pilot_audit_path: Path, failed_root: Path | None,
          output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    if case_kind not in PURPOSES:
        raise ValueError(f"unknown case kind: {case_kind}")
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    pilot_manifest_path = pilot_root / "manifest.json"
    pilot_audit = json.loads(pilot_audit_path.read_text(encoding="utf-8"))
    if (manifest.get("purpose") != PURPOSES[case_kind]
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or pilot_audit.get("status")
            != "GaN_native_basin_10step_pilots_raw_audited_not_phase_or_endpoint_certificates"
            or pilot_audit.get("source_sha256", {}).get("manifest")
            != sha256(pilot_manifest_path)):
        raise ValueError("source manifest/audit contract differs")
    case_name = "negative_mode_plus" if case_kind == "uninterrupted_plus" \
        else "negative_mode_minus"
    case_dir = work_root / ("case" if case_kind == "uninterrupted_plus"
                            else f"cases/{case_name}")
    case_record = (manifest if case_kind == "uninterrupted_plus" else next(
        item for item in manifest["cases"] if item["name"] == case_name
    ))
    if any(sha256(case_dir / filename) != expected
           for filename, expected in case_record["input_sha256"].items()):
        raise ValueError("VASP input hash differs from staged manifest")
    if b"\r" in (case_dir / "INCAR").read_bytes():
        raise ValueError("INCAR contains CRLF")
    pilot_record = next(item for item in pilot_audit["cases"]
                        if item["case"] == case_name)
    pilot_dir = pilot_root / "cases" / case_name
    pilot_outcar = pilot_dir / "OUTCAR"
    if sha256(pilot_outcar) != pilot_record["outcar_sha256"]:
        raise ValueError("pilot OUTCAR differs from raw audit")
    pilot_last = read(pilot_outcar)
    source_geometry = read(case_dir / "POSCAR", format="vasp")
    if case_kind == "continued_minus":
        if (manifest.get("source_sha256", {}).get("pilot_audit")
                != sha256(pilot_audit_path)
                or case_record["source_CONTCAR_sha256"]
                != sha256(pilot_dir / "CONTCAR")
                or not same_geometry(source_geometry, pilot_last, tolerance=2e-5)):
            raise ValueError("minus continuation source geometry/provenance differs")
    else:
        if (manifest.get("source_sha256", {}).get("pilot_audit")
                != sha256(pilot_audit_path)
                or failed_root is None
                or manifest.get("source_sha256", {}).get("failed_restart_OUTCAR")
                != sha256(failed_root / "cases/negative_mode_plus/OUTCAR")
                or not same_geometry(source_geometry,
                                     read(pilot_dir / "POSCAR", format="vasp"),
                                     tolerance=2e-5)):
            raise ValueError("uninterrupted plus seed/provenance differs")
    outcar = case_dir / "OUTCAR"
    raw = outcar.read_text(encoding="utf-8", errors="replace")
    if ("General timing and accounting informations" not in raw
            or "aborting loop because EDIFF is reached" not in raw
            or "PSTRESS=  457.0" not in raw):
        raise ValueError("VASP run lacks terminal accounting, SCF or target pressure")
    triples = energy_triples(raw)
    if not triples or len(triples) > manifest.get("ionic_settings", {}).get("NSW", 100):
        # The uninterrupted manifest declares its cap in changed text, not an
        # ionic_settings object; 100 is its deliberately fixed upper bound.
        raise ValueError("ionic step count exceeds declared cap")
    pressure_GPa = 45.7
    last_atoms = read(outcar)
    if (len(last_atoms) != 4 or not np.isfinite(last_atoms.get_forces()).all()
            or not np.isfinite(last_atoms.get_stress(voigt=False)).all()
            or not np.isclose(triples[-1][2], pressure_GPa * GPa
                              * last_atoms.get_volume(), atol=2e-4, rtol=0)):
        raise ValueError("last evaluated geometry/force/stress/PV invalid")
    pilot_triples = energy_triples(pilot_outcar.read_text(encoding="utf-8", errors="replace"))
    if len(pilot_triples) != 10:
        raise ValueError("pilot did not have ten evaluated geometries")
    if case_kind == "uninterrupted_plus":
        if len(triples) < 10:
            raise ValueError("long run stopped before reproducing the pilot")
        first_ten_max_H_difference = max(abs(triples[i][1] - pilot_triples[i][1])
                                         for i in range(10))
        initial_H_difference = triples[0][1] - pilot_triples[0][1]
    else:
        first_ten_max_H_difference = None
        initial_H_difference = triples[0][1] - pilot_triples[-1][1]
    contcar = case_dir / "CONTCAR"
    next_atoms = read(contcar, format="vasp")
    fmax = float(np.max(np.linalg.norm(last_atoms.get_forces(), axis=1)))
    stress_residual = residual_stress_kbar(last_atoms.get_stress(voigt=False),
                                           pressure_GPa)
    result = {
        "status": "GaN_basin_followup_raw_audited_not_endpoint_phase_certificate",
        "case_kind": case_kind, "case": case_name,
        "n_ionic_energy_records": len(triples),
        "nsw_cap": 100 if case_kind == "uninterrupted_plus" else 40,
        "initial_H_difference_from_audited_pilot_eV": initial_H_difference,
        "first_ten_max_H_difference_from_pilot_eV": first_ten_max_H_difference,
        "last_free_energy_eV_per_cell": triples[-1][0],
        "last_PV_eV_per_cell": triples[-1][2],
        "last_enthalpy_eV_per_cell": triples[-1][1],
        "last_enthalpy_minus_first_eV_per_cell": triples[-1][1] - triples[0][1],
        "lowest_evaluated_enthalpy_eV_per_cell": min(item[1] for item in triples),
        "last_volume_A3": float(last_atoms.get_volume()),
        "last_max_atomic_force_eV_per_A": fmax,
        "last_max_stress_residual_from_45p7GPa_kbar": stress_residual,
        "last_coordination_GaN_2p4A": ga_n_coordination(last_atoms),
        "contcar_is_last_evaluated_geometry": same_geometry(
            last_atoms, next_atoms, tolerance=2e-5),
        "vasp_reports_ionic_convergence": (
            "reached required accuracy - stopping structural energy minimisation" in raw
        ),
        "force_pass_0p02eV_per_A": fmax <= 0.02,
        "stress_pass_1kbar": stress_residual <= 1.0,
        "source_sha256": {"manifest": sha256(manifest_path),
                          "pilot_manifest": sha256(pilot_manifest_path),
                          "pilot_audit": sha256(pilot_audit_path),
                          "pilot_OUTCAR": sha256(pilot_outcar),
                          "OUTCAR": sha256(outcar),
                          "CONTCAR": sha256(contcar),
                          "auditor": sha256(Path(__file__))},
        "limitations": [
            "A 2.4 Å coordination count is a topology screen, not a phase certificate.",
            "Raw E, H and P*V are per four-atom cell; never add P*V to the VASP H field twice.",
            "If CONTCAR is unevaluated, assign the last H only to the last OUTCAR geometry.",
            "B4/B1 phase assignment requires structural comparison and final comparable static checks.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"],
                      "n_records": len(triples),
                      "H_eV": triples[-1][1], "fmax": fmax,
                      "stress_residual_kbar": stress_residual}))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-kind", choices=sorted(PURPOSES), required=True)
    for name in ("work-root", "pilot-root", "pilot-audit", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--failed-root", type=Path)
    args = parser.parse_args()
    audit(args.case_kind, args.work_root, args.pilot_root, args.pilot_audit,
          args.failed_root, args.output)


if __name__ == "__main__":
    main()
