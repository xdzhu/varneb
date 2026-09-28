"""Audit both 600-eV GaN uninterrupted basin tests from complete VASP output.

The first ten evaluated enthalpies must reproduce the corresponding signed
pilot before any continuation is used as a continuous trajectory.  Phase
identity still requires structure, force, stress and static-energy checks.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_ts_basin_followups import residual_stress_kbar
from scripts.audit_gan_ts_basin_pilot import energy_triples, ga_n_coordination
from scripts.compare_gan_basin_endpoint import relative_displacement_metrics
from scripts.prepare_gan_600eV_basin_restart_canary import CASES, sha256
from scripts.prepare_gan_600eV_basin_uninterrupted import PURPOSE
from scripts.prepare_gan_600eV_ts_hessian import same_geometry


def audit(work: Path, pilot: Path, pilot_audit_path: Path,
          chain_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    pilot_audit = json.loads(pilot_audit_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != PURPOSE
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("PSTRESS_kbar") != 457.0
            or manifest.get("encut_eV") != 600
            or manifest.get("n_cases") != 2
            or tuple(item["name"] for item in manifest["cases"]) != CASES
            or manifest.get("source_sha256", {}).get("pilot_manifest")
            != sha256(pilot / "manifest.json")
            or manifest.get("source_sha256", {}).get("pilot_audit")
            != sha256(pilot_audit_path)
            or manifest.get("source_sha256", {}).get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_basin_uninterrupted.py"))
            or pilot_audit.get("source_sha256", {}).get("manifest")
            != sha256(pilot / "manifest.json")):
        raise ValueError("uninterrupted basin source/provenance differs")
    reference = read(chain_path, index=":")
    if (len(reference) != 29
            or ga_n_coordination(reference[0]) != [4, 4]
            or ga_n_coordination(reference[-1]) != [6, 6]):
        raise ValueError("source B4/B1 endpoint identities differ")
    pilot_results = {item["case"]: item for item in pilot_audit["cases"]}
    cases = []
    for record in manifest["cases"]:
        name = record["name"]
        directory = work / "cases" / name
        inputs = record["input_sha256"]
        original = pilot / "cases" / name
        if (json.loads((directory / "sha256.inputs.json").read_text(encoding="utf-8")) != inputs
                or any(sha256(directory / key) != digest for key, digest in inputs.items())
                or sha256(original / "OUTCAR") != record["source_pilot_OUTCAR_sha256"]
                or record["source_pilot_input_sha256"]
                != pilot_results[name]["input_sha256"]
                or not same_geometry(read(directory / "POSCAR", format="vasp"),
                                     read(original / "POSCAR", format="vasp"))):
            raise ValueError(f"uninterrupted input or original pilot differs: {name}")
        outcar = directory / "OUTCAR"
        raw = outcar.read_text(encoding="utf-8", errors="replace")
        if ("General timing and accounting informations" not in raw
                or "PSTRESS=  457.0" not in raw):
            raise ValueError(f"uninterrupted output incomplete or pressure wrong: {name}")
        triples = energy_triples(raw)
        if (len(triples) < 10 or len(triples) > 100
                or raw.count("aborting loop because EDIFF is reached") != len(triples)):
            raise ValueError(f"ionic evaluations or SCF completion invalid: {name}")
        pilot_triples = energy_triples((original / "OUTCAR").read_text(
            encoding="utf-8", errors="replace"
        ))
        if len(pilot_triples) != 10:
            raise ValueError("original signed pilot does not have ten evaluations")
        first_ten_H_differences = [float(triples[i][1] - pilot_triples[i][1])
                                   for i in range(10)]
        final = read(outcar)
        if (len(final) != 4
                or final.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
                or not np.isfinite(final.get_forces()).all()
                or not np.isfinite(final.get_stress(voigt=False)).all()
                or not np.isclose(triples[-1][2], 45.7 * GPa * final.get_volume(),
                                  atol=2e-4, rtol=0)):
            raise ValueError(f"final force/stress/enthalpy closure invalid: {name}")
        next_geometry = read(directory / "CONTCAR", format="vasp")
        fmax = float(np.linalg.norm(final.get_forces(), axis=1).max())
        stress = residual_stress_kbar(final.get_stress(voigt=False), 45.7)
        wavecar = directory / "WAVECAR"
        cases.append({
            "case": name,
            "n_ionic_energy_records": len(triples),
            "first_ten_H_differences_from_pilot_eV_per_cell": first_ten_H_differences,
            "first_ten_max_abs_H_difference_eV_per_cell": max(
                abs(value) for value in first_ten_H_differences
            ),
            "first_ten_reproduce_within_1meV_per_cell": max(
                abs(value) for value in first_ten_H_differences
            ) <= 0.001,
            "enthalpy_eV_per_cell_by_ionic_step": [item[1] for item in triples],
            "lowest_evaluated_enthalpy_eV_per_cell": min(item[1] for item in triples),
            "last_evaluated_enthalpy_eV_per_cell": triples[-1][1],
            "last_evaluated_volume_A3": float(final.get_volume()),
            "last_evaluated_max_atomic_force_eV_per_A": fmax,
            "last_evaluated_max_stress_residual_kbar": stress,
            "last_evaluated_GaN_coordination_2p4A": ga_n_coordination(final),
            "mapped_endpoint_distances_A": {
                phase: relative_displacement_metrics(final, endpoint)
                for phase, endpoint in (("B4", reference[0]), ("B1", reference[-1]))
            },
            "contcar_is_last_evaluated_geometry": same_geometry(final, next_geometry),
            "vasp_reports_ionic_convergence":
                "reached required accuracy - stopping structural energy minimisation" in raw,
            "force_pass_0p02eV_per_A": fmax <= 0.02,
            "stress_pass_2kbar": stress <= 2.0,
            "wavecar_archived_locally": wavecar.exists(),
            "outcar_sha256": sha256(outcar),
            "contcar_sha256": sha256(directory / "CONTCAR"),
            "wavecar_size_bytes": wavecar.stat().st_size if wavecar.exists() else 0,
        })
    summary = {
        "status": "GaN_600eV_uninterrupted_basin_tests_raw_audited_not_phase_certificates",
        "pressure_GPa": 45.7,
        "encut_eV": 600,
        "n_cases": 2,
        "all_first_ten_reproduce_within_1meV_per_cell": all(
            case["first_ten_reproduce_within_1meV_per_cell"] for case in cases
        ),
        "cases": cases,
        "source_sha256": {"manifest": sha256(manifest_path),
                          "pilot_audit": sha256(pilot_audit_path),
                          "reference_chain": sha256(chain_path),
                          "auditor": sha256(Path(__file__))},
        "limitations": [
            "The 1-meV first-ten gate is a reproducibility screen, not a basis-set convergence study.",
            "Coordination and endpoint distances are structural screens, not full phase certificates.",
            "True B4/B1 basin connection requires converged forces/stress and comparable final statics.",
            "The endpoint stress acceptance threshold is 2 kbar; it is not a DFT input parameter.",
            "WAVECAR size refers only to an optional local archive, not to the remote Slurm run.",
            "A VASP-reported stop does not alone certify a variable-cell transition state.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "pilot", "pilot-audit", "chain", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    summary = audit(args.work, args.pilot, args.pilot_audit,
                    args.chain, args.output)
    print(json.dumps({"status": summary["status"],
                      "all_first_ten_reproduce_within_1meV_per_cell":
                      summary["all_first_ten_reproduce_within_1meV_per_cell"]}))


if __name__ == "__main__":
    main()
