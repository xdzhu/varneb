"""Stage isolated, pressure-consistent GaN basin continuations from audited pilots.

The input is not a phase certificate. No original VCNEB or ten-step file is
modified; each branch starts from its last *evaluated* pilot CONTCAR.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil

import numpy as np
from ase.io import read

from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256


CASES = {"negative_mode_plus", "negative_mode_minus"}


def continuation_incar(raw: bytes) -> bytes:
    if b"\r" in raw or not raw.endswith(b"\n"):
        raise ValueError("source VASP INCAR must be LF-only")
    text = raw.decode("ascii")
    expected = {
        "IBRION": "2", "ISIF": "3", "NSW": "10", "PSTRESS": "457.0",
        "EDIFFG": "-0.02", "POTIM": "0.25", "ENCUT": "1000.000000",
        "EDIFF": "1.00e-07", "ISYM": "-1", "SYMPREC": "1.00e-04",
    }
    found: dict[str, str] = {}
    lines = text.splitlines()
    for index, line in enumerate(lines):
        match = re.fullmatch(r"\s*([A-Za-z_]+)\s*=\s*(.*?)\s*", line)
        if not match:
            continue
        key, value = match.group(1).upper(), match.group(2)
        if key in found:
            raise ValueError(f"duplicate INCAR key: {key}")
        found[key] = value
        if key == "NSW":
            lines[index] = " NSW = 40"
    if any(found.get(key) != value for key, value in expected.items()):
        raise ValueError("pilot physical or relaxation setting changed")
    return ("\n".join(lines) + "\n").encode("ascii")


def geometry_metrics(atoms) -> dict:
    cell = np.asarray(atoms.cell.array, dtype=float)
    distances = atoms.get_all_distances(mic=True)
    minimum = min(float(distances[i, j]) for i in range(len(atoms))
                  for j in range(i + 1, len(atoms)))
    metrics = {
        "volume_A3": float(atoms.get_volume()),
        "cell_condition": float(np.linalg.cond(cell)),
        "minimum_pair_distance_A": minimum,
    }
    if (len(atoms) != 4 or atoms.get_chemical_symbols().count("Ga") != 2
            or atoms.get_chemical_symbols().count("N") != 2
            or not all(np.isfinite(value) for value in metrics.values())
            or not 25.0 < metrics["volume_A3"] < 50.0
            or metrics["cell_condition"] >= 10.0
            or minimum < 1.65):
        raise ValueError(f"unsafe GaN continuation geometry: {metrics}")
    return metrics


def prepare(pilot_root: Path, pilot_audit_path: Path, work_root: Path) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    pilot_manifest_path = pilot_root / "manifest.json"
    pilot_manifest = json.loads(pilot_manifest_path.read_text(encoding="utf-8"))
    audit = json.loads(pilot_audit_path.read_text(encoding="utf-8"))
    if (pilot_manifest.get("purpose")
            != "GaN_1000eV_signed_basin_native_VASP_10step_pilot_not_endpoint_certificate"
            or audit.get("status")
            != "GaN_native_basin_10step_pilots_raw_audited_not_phase_or_endpoint_certificates"
            or audit.get("source_sha256", {}).get("manifest") != sha256(pilot_manifest_path)
            or audit.get("pressure_GPa") != 45.7
            or audit.get("PSTRESS_kbar") != 457.0):
        raise ValueError("source pilot audit or pressure contract differs")
    audited = {record["case"]: record for record in audit["cases"]}
    if set(audited) != CASES:
        raise ValueError("both pilot branches must be audited")
    cases = []
    for name in sorted(CASES):
        source_dir = pilot_root / "cases" / name
        record = audited[name]
        if (record["n_ionic_energy_records"] != 10
                or not record["contcar_is_last_evaluated_geometry"]
                or sha256(source_dir / "OUTCAR") != record["outcar_sha256"]
                or sha256(source_dir / "CONTCAR") != record["contcar_sha256"]):
            raise ValueError(f"pilot branch is not an audited evaluated geometry: {name}")
        atoms = read(source_dir / "CONTCAR", format="vasp")
        if not same_geometry(atoms, read(source_dir / "OUTCAR"), tolerance=2e-5):
            raise ValueError(f"CONTCAR differs from last evaluated OUTCAR: {name}")
        metrics = geometry_metrics(atoms)
        destination = work_root / "cases" / name
        destination.mkdir(parents=True, exist_ok=False)
        shutil.copy2(source_dir / "CONTCAR", destination / "POSCAR")
        for filename in ("KPOINTS", "POTCAR"):
            shutil.copy2(source_dir / filename, destination / filename)
        (destination / "INCAR").write_bytes(
            continuation_incar((source_dir / "INCAR").read_bytes())
        )
        hashes = {filename: sha256(destination / filename)
                  for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
        (destination / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
        )
        cases.append({
            "name": name, "input_sha256": hashes,
            "source_OUTCAR_sha256": record["outcar_sha256"],
            "source_CONTCAR_sha256": record["contcar_sha256"],
            "source_last_enthalpy_eV_per_cell": record["last_evaluated_enthalpy_eV"],
            "source_last_max_atomic_force_eV_per_A": record[
                "last_evaluated_max_atomic_force_eV_per_A"],
            "initial_geometry_preflight": metrics,
        })
    result = {
        "purpose": "GaN_1000eV_basin_continuation_from_audited_pilot_not_phase_certificate",
        "status": "inputs_finalized_no_DFT", "pressure_GPa": 45.7,
        "PSTRESS_kbar": 457.0, "encut_eV": 1000,
        "ionic_settings": {"IBRION": 2, "ISIF": 3, "NSW": 40,
                           "EDIFFG_eV_per_A": -0.02, "POTIM": 0.25},
        "cases": cases,
        "source_sha256": {"pilot_manifest": sha256(pilot_manifest_path),
                          "pilot_audit": sha256(pilot_audit_path),
                          "preparer": sha256(Path(__file__))},
        "limitations": [
            "Forty further ionic steps are a bounded continuation, not a guaranteed local minimum.",
            "Do not assign B4/B1 phases until geometry, force, stress and static enthalpy are audited.",
        ],
    }
    (work_root / "manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-root", type=Path, required=True)
    parser.add_argument("--pilot-audit", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.pilot_root, args.pilot_audit, args.work_root)
    print(json.dumps({"status": result["status"], "n_cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
