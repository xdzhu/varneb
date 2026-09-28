"""Compare GaN's original FFT-grid restart with its automatic-grid restart."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_ts_basin_pilot import energy_triples
from scripts.prepare_gan_600eV_fft_restart_canary import PURPOSE, sha256
from scripts.prepare_gan_600eV_ts_hessian import same_geometry


def audit(work: Path, automatic: Path, automatic_audit_path: Path,
          pilot: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    automatic_audit = json.loads(automatic_audit_path.read_text(encoding="utf-8"))
    if (manifest.get("purpose") != PURPOSE
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("encut_eV") != 600
            or manifest.get("source_sha256", {}).get("automatic_restart_manifest")
            != sha256(automatic / "manifest.json")
            or manifest.get("source_sha256", {}).get("automatic_restart_audit")
            != sha256(automatic_audit_path)
            or manifest.get("source_sha256", {}).get("pilot_OUTCAR")
            != sha256(pilot / "cases/grid_um_vz/OUTCAR")
            or manifest.get("source_sha256", {}).get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_fft_restart_canary.py"))):
        raise ValueError("FFT-grid canary provenance differs")
    case = work / "case"
    if any(sha256(case / key) != digest
           for key, digest in manifest["input_sha256"].items()):
        raise ValueError("FFT-grid canary input differs")
    raw = (case / "OUTCAR").read_text(encoding="utf-8", errors="replace")
    if ("General timing and accounting informations" not in raw
            or "PSTRESS=  457.0" not in raw
            or raw.count("aborting loop because EDIFF is reached") != 1
            or "dimension x,y,z NGX =    24 NGY =   24 NGZ =   40" not in raw
            or "dimension x,y,z NGXF=    48 NGYF=   48 NGZF=   80" not in raw):
        raise ValueError("locked FFT-grid canary did not complete as declared")
    values = energy_triples(raw)
    atoms = read(case / "OUTCAR")
    if (len(values) != 1
            or not same_geometry(read(case / "POSCAR", format="vasp"), atoms)
            or not np.isfinite(atoms.get_forces()).all()
            or not np.isfinite(atoms.get_stress(voigt=False)).all()
            or not np.isclose(values[0][2], 45.7 * GPa * atoms.get_volume(),
                              atol=2e-4, rtol=0)):
        raise ValueError("locked-grid enthalpy/geometry/force/stress invalid")
    negative = next(item for item in automatic_audit["cases"]
                    if item["case"] == "grid_um_vz")
    if (sha256(automatic / "cases/grid_um_vz/OUTCAR") != negative["outcar_sha256"]
            or manifest["source_automatic_restart_enthalpy_eV_per_cell"]
            != negative["enthalpy_eV_per_cell"]):
        raise ValueError("automatic-grid comparator differs")
    pilot_H = manifest["source_pilot_last_enthalpy_eV_per_cell"]
    auto_H = manifest["source_automatic_restart_enthalpy_eV_per_cell"]
    result = {
        "status": "GaN_600eV_original_FFT_grid_one_step_raw_audited",
        "pressure_GPa": 45.7,
        "encut_eV": 600,
        "n_ionic_energy_records": 1,
        "pilot_last_H_eV_per_cell": pilot_H,
        "automatic_restart_H_eV_per_cell": auto_H,
        "original_grid_restart_H_eV_per_cell": values[0][1],
        "automatic_minus_pilot_meV_per_cell": 1000 * (auto_H - pilot_H),
        "original_grid_minus_pilot_meV_per_cell": 1000 * (values[0][1] - pilot_H),
        "original_grid_reproduces_pilot_within_1meV_per_cell":
            abs(values[0][1] - pilot_H) <= 0.001,
        "original_grid_coarse": [24, 24, 40],
        "original_grid_fine": [48, 48, 80],
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "OUTCAR": sha256(case / "OUTCAR"),
            "automatic_restart_audit": sha256(automatic_audit_path),
            "pilot_OUTCAR": sha256(pilot / "cases/grid_um_vz/OUTCAR"),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "This isolates one numerical-grid restart difference; no future cell is guaranteed grid-safe.",
            "A one-step static comparison is not a basin, phase, or TS certificate.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "automatic", "automatic-audit", "pilot", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.work, args.automatic, args.automatic_audit,
                   args.pilot, args.output)
    print(json.dumps({"status": result["status"],
                      "original_grid_minus_pilot_meV_per_cell":
                      result["original_grid_minus_pilot_meV_per_cell"]}))


if __name__ == "__main__":
    main()
