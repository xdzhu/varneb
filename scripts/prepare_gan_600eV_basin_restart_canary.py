"""Stage a one-step, same-contract VASP restart test from audited GaN pilots.

The negative branch ended in an empirically Bravais-sensitive cell.  These
canaries test the *actual restart POSCAR* before a longer native relaxation;
they do not certify future cells or either basin.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from ase.io import read

from scripts.prepare_gan_600eV_ts_basin_pilot import CASES, sha256
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


PURPOSE = "GaN_45p7_600eV_native_VASP_basin_restart_one_step_canary"
PILOT_PURPOSE = "GaN_45p7_600eV_signed_native_VASP_basin_10step_pilot"
AUDIT_STATUS = "GaN_600eV_basin_10step_pilots_raw_audited_not_phase_certificates"


def canary_incar(source: bytes) -> bytes:
    if (b"\r" in source or source.count(b" NSW = 10\n") != 1
            or source.count(b" ENCUT = 600.000000\n") != 1
            or source.count(b" PSTRESS = 457.0\n") != 1
            or source.count(b" ISIF = 3\n") != 1
            or source.count(b" IBRION = 2\n") != 1):
        raise ValueError("pilot INCAR is not the audited 600-eV pressure contract")
    return source.replace(b" NSW = 10\n", b" NSW = 1\n")


def prepare(pilot_root: Path, audit_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = pilot_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != PILOT_PURPOSE
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("PSTRESS_kbar") != 457.0
            or manifest.get("encut_eV") != 600
            or audit.get("status") != AUDIT_STATUS
            or audit.get("pressure_GPa") != 45.7
            or audit.get("PSTRESS_kbar") != 457.0
            or audit.get("source_sha256", {}).get("manifest") != sha256(manifest_path)):
        raise ValueError("source pilot pressure/provenance audit differs")
    records = {item["name"]: item for item in manifest["cases"]}
    audited = {item["case"]: item for item in audit["cases"]}
    if set(records) != set(CASES) or set(audited) != set(CASES):
        raise ValueError("both signed pilot branches are required")
    checked = []
    for name in CASES:
        source = pilot_root / "cases" / name
        item, result = records[name], audited[name]
        hashes = item["input_sha256"]
        if (result["n_ionic_energy_records"] != 10
                or not result["contcar_is_last_evaluated_geometry"]
                or result["input_sha256"] != hashes
                or any(sha256(source / key) != digest for key, digest in hashes.items())
                or sha256(source / "OUTCAR") != result["outcar_sha256"]
                or sha256(source / "CONTCAR") != result["contcar_sha256"]
                or hashes["KPOINTS"] != PRODUCTION_INPUT_SHA256["KPOINTS"]
                or hashes["POTCAR"] != PRODUCTION_INPUT_SHA256["POTCAR"]):
            raise ValueError(f"unaudited restart source: {name}")
        atoms = read(source / "CONTCAR", format="vasp")
        if (not same_geometry(atoms, read(source / "OUTCAR"))
                or atoms.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]):
            raise ValueError(f"CONTCAR was not the evaluated structure: {name}")
        geometry = geometry_preflight(atoms)
        if (not 25.0 < geometry["volume_A3"] < 50.0
                or geometry["minimum_distance_A"] < 1.65):
            raise ValueError(f"unsafe restart geometry: {name}")
        checked.append((name, source, result, geometry))
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    cases = []
    for name, source, result, geometry in checked:
        destination = output / "cases" / name
        destination.mkdir()
        shutil.copy2(source / "CONTCAR", destination / "POSCAR")
        for filename in ("KPOINTS", "POTCAR"):
            shutil.copy2(source / filename, destination / filename)
        (destination / "INCAR").write_bytes(canary_incar((source / "INCAR").read_bytes()))
        hashes = {key: sha256(destination / key)
                  for key in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
        (destination / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        cases.append({
            "name": name,
            "source_last_enthalpy_eV_per_cell": result["last_evaluated_enthalpy_eV_per_cell"],
            "source_OUTCAR_sha256": result["outcar_sha256"],
            "source_CONTCAR_sha256": result["contcar_sha256"],
            "input_sha256": hashes,
            "geometry_preflight": geometry,
        })
    summary = {
        "status": "inputs_finalized_no_DFT",
        "purpose": PURPOSE,
        "pressure_GPa": 45.7,
        "PSTRESS_kbar": 457.0,
        "encut_eV": 600,
        "ionic_settings": {"IBRION": 2, "ISIF": 3, "NSW": 1,
                           "EDIFFG_eV_per_A": -0.02, "POTIM": 0.25},
        "n_cases": 2,
        "cases": cases,
        "source_sha256": {"pilot_manifest": sha256(manifest_path),
                          "pilot_audit": sha256(audit_path),
                          "preparer": sha256(Path(__file__))},
        "limitations": [
            "One successful restart step does not guarantee all later cells avoid VASP's Bravais classifier.",
            "No one-step output is a converged phase or TS certificate.",
        ],
    }
    (output / "manifest.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-root", type=Path, required=True)
    parser.add_argument("--pilot-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = prepare(args.pilot_root, args.pilot_audit, args.output)
    print(json.dumps({"status": summary["status"], "n_cases": summary["n_cases"]}))


if __name__ == "__main__":
    main()
