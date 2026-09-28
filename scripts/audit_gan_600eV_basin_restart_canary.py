"""Audit same-600-eV GaN restart canaries before any longer continuation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_ts_basin_pilot import energy_triples
from scripts.prepare_gan_600eV_basin_restart_canary import CASES, PURPOSE, sha256
from scripts.prepare_gan_600eV_ts_hessian import same_geometry


def audit(work: Path, pilot: Path, pilot_audit_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    pilot_audit = json.loads(pilot_audit_path.read_text(encoding="utf-8"))
    if (manifest.get("purpose") != PURPOSE
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("PSTRESS_kbar") != 457.0
            or manifest.get("encut_eV") != 600
            or manifest.get("n_cases") != 2
            or manifest.get("source_sha256", {}).get("pilot_manifest")
            != sha256(pilot / "manifest.json")
            or manifest.get("source_sha256", {}).get("pilot_audit")
            != sha256(pilot_audit_path)
            or manifest.get("source_sha256", {}).get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_basin_restart_canary.py"))
            or tuple(case["name"] for case in manifest["cases"]) != CASES):
        raise ValueError("restart canary manifest/provenance changed")
    audited = {item["case"]: item for item in pilot_audit["cases"]}
    cases = []
    for record in manifest["cases"]:
        name = record["name"]
        directory = work / "cases" / name
        inputs = record["input_sha256"]
        if (json.loads((directory / "sha256.inputs.json").read_text(encoding="utf-8")) != inputs
                or any(sha256(directory / key) != digest for key, digest in inputs.items())
                or record["source_OUTCAR_sha256"] != audited[name]["outcar_sha256"]
                or record["source_CONTCAR_sha256"] != audited[name]["contcar_sha256"]
                or sha256(pilot / "cases" / name / "OUTCAR") != record["source_OUTCAR_sha256"]
                or sha256(pilot / "cases" / name / "CONTCAR") != record["source_CONTCAR_sha256"]):
            raise ValueError(f"restart canary input/pilot provenance changed: {name}")
        outcar = directory / "OUTCAR"
        raw = outcar.read_text(encoding="utf-8", errors="replace")
        if ("General timing and accounting informations" not in raw
                or "PSTRESS=  457.0" not in raw
                or raw.count("aborting loop because EDIFF is reached") != 1):
            raise ValueError(f"restart canary output incomplete: {name}")
        triples = energy_triples(raw)
        if len(triples) != 1:
            raise ValueError(f"restart canary has wrong ionic step count: {name}")
        atoms = read(outcar)
        seed = read(directory / "POSCAR", format="vasp")
        if (not same_geometry(seed, atoms)
                or len(atoms) != 4
                or atoms.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
                or not np.isfinite(atoms.get_forces()).all()
                or not np.isfinite(atoms.get_stress(voigt=False)).all()
                or not np.isclose(triples[0][2], 45.7 * GPa * atoms.get_volume(),
                                  atol=2e-4, rtol=0)):
            raise ValueError(f"restart canary geometry/forces/stress/PV failed: {name}")
        delta = float(triples[0][1] - record["source_last_enthalpy_eV_per_cell"])
        cases.append({
            "case": name,
            "n_ionic_energy_records": 1,
            "pressure_GPa": 45.7,
            "enthalpy_eV_per_cell": triples[0][1],
            "enthalpy_difference_from_pilot_last_eV_per_cell": delta,
            "reproduces_pilot_enthalpy_within_1meV_per_cell": abs(delta) <= 0.001,
            "initial_geometry_matches_last_evaluated_pilot": True,
            "outcar_sha256": sha256(outcar),
            "geometry_preflight": record["geometry_preflight"],
        })
    result = {
        "status": "GaN_600eV_basin_restart_canaries_raw_audited_not_future_cell_guarantee",
        "n_cases": 2,
        "pressure_GPa": 45.7,
        "cases": cases,
        "all_reproduce_within_1meV_per_cell": all(
            case["reproduces_pilot_enthalpy_within_1meV_per_cell"] for case in cases
        ),
        "source_sha256": {"manifest": sha256(manifest_path),
                          "pilot_audit": sha256(pilot_audit_path),
                          "auditor": sha256(Path(__file__))},
        "limitations": [
            "A successful initial restart cannot predict future Bravais classification during cell relaxation.",
            "This is a same-geometry numerical and parser check, not a basin or TS certificate.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--pilot-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.work, args.pilot, args.pilot_audit, args.output)
    print(json.dumps({"status": result["status"],
                      "all_reproduce_within_1meV_per_cell":
                      result["all_reproduce_within_1meV_per_cell"]}))


if __name__ == "__main__":
    main()
