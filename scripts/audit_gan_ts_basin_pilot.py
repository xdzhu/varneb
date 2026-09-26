"""Audit native VASP GaN basin pilots without confusing E, H and next POSCAR.

At a finite NSW cap, CONTCAR may be an unevaluated next geometry. This audit
reports that distinction and never assigns a B4/B1 phase prematurely.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list
from ase.units import GPa

from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256


NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
FREE = re.compile(rf"free\s+energy\s+TOTEN\s*=\s*({NUMBER})\s*eV", re.I)
ENTHALPY = re.compile(rf"enthalpy\s+is\s+TOTEN\s*=\s*({NUMBER})\s*eV\s+P\s*V\s*=\s*({NUMBER})", re.I)


def energy_triples(outcar_text: str) -> list[tuple[float, float, float]]:
    triples = []
    latest_free = None
    for line in outcar_text.splitlines():
        free_match = FREE.search(line)
        if free_match:
            # OUTCAR also prints intermediate electronic-iteration free energies.
            # The ionic enthalpy line pairs only with the last preceding one.
            latest_free = float(free_match.group(1).replace("D", "E"))
        enthalpy_match = ENTHALPY.search(line)
        if enthalpy_match:
            if latest_free is None:
                raise ValueError("ionic enthalpy has no preceding free energy")
            h = float(enthalpy_match.group(1).replace("D", "E"))
            pv = float(enthalpy_match.group(2).replace("D", "E"))
            triples.append((latest_free, h, pv))
            latest_free = None
    if not triples:
        raise ValueError("incomplete paired VASP free-energy / enthalpy records")
    if any(not np.all(np.isfinite(triple)) or abs(triple[0] + triple[2] - triple[1]) > 5e-7
           for triple in triples):
        raise ValueError("VASP reported E + PV and H do not close")
    return triples


def ga_n_coordination(atoms, cutoff_A: float = 2.4) -> list[int]:
    """Count periodic Ga--N pairs around each Ga; not a phase certificate."""

    symbols = atoms.get_chemical_symbols()
    if symbols.count("Ga") != 2 or symbols.count("N") != 2:
        raise ValueError("GaN pilot requires two Ga and two N in the mapped cell")
    i, j = neighbor_list("ij", atoms, cutoff_A)
    return [int(sum(symbols[neighbor] == "N" for owner, neighbor in zip(i, j)
                    if owner == ga))
            for ga, symbol in enumerate(symbols) if symbol == "Ga"]


def audit(work_root: Path, signed_root: Path, signed_audit: Path,
          path_trajectory: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("purpose") != "GaN_1000eV_signed_basin_native_VASP_10step_pilot_not_endpoint_certificate"
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("PSTRESS_kbar") != 457.0
            or manifest.get("ionic_settings", {}).get("NSW") != 10
            or manifest.get("source_sha256", {}).get("signed_manifest") != sha256(signed_root / "manifest.json")
            or manifest.get("source_sha256", {}).get("signed_audit") != sha256(signed_audit)
            or manifest.get("source_sha256", {}).get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_ts_basin_pilot.py"))):
        raise ValueError("pilot source, pressure or input contract differs")
    pressure = float(manifest["pressure_GPa"]) * GPa
    chain = read(path_trajectory, index="-29:")
    if len(chain) != 29:
        raise ValueError("the source VCNEB chain is not 29 images")
    endpoint_coordination = {
        "B4_image00": ga_n_coordination(chain[0]),
        "B1_image28": ga_n_coordination(chain[-1]),
    }
    if endpoint_coordination != {"B4_image00": [4, 4], "B1_image28": [6, 6]}:
        raise ValueError("the declared 2.4 Å topology screen does not separate source endpoints")
    results = []
    for record in manifest["cases"]:
        directory = work_root / "cases" / record["name"]
        source = signed_root / "cases" / record["name"]
        hashes = json.loads((directory / "sha256.inputs.json").read_text(encoding="utf-8"))
        if (hashes != record["input_sha256"]
                or any(sha256(directory / filename) != expected
                       for filename, expected in hashes.items())
                or sha256(source / "OUTCAR") != record["source_OUTCAR_sha256"]
                or b"\r" in (directory / "INCAR").read_bytes()):
            raise ValueError(f"source or LF-only VASP input changed: {directory}")
        outcar = directory / "OUTCAR"
        text = outcar.read_text(encoding="utf-8", errors="replace")
        if ("General timing and accounting informations" not in text
                or "aborting loop because EDIFF is reached" not in text
                or "PSTRESS=  457.0" not in text):
            raise ValueError(f"pilot VASP run incomplete or pressure differs: {outcar}")
        triples = energy_triples(text)
        last_atoms = read(outcar)
        if (len(last_atoms) != 4 or not np.isfinite(last_atoms.get_forces()).all()
                or not np.isfinite(last_atoms.get_stress(voigt=False)).all()
                or not np.isclose(triples[-1][2], pressure * last_atoms.get_volume(), atol=2e-4, rtol=0)):
            raise ValueError(f"final VASP volume, force, stress or pressure work differs: {outcar}")
        contcar = directory / "CONTCAR"
        if not contcar.is_file():
            raise FileNotFoundError(contcar)
        next_geometry = read(contcar, format="vasp")
        results.append({
            "case": record["name"], "sign": record["sign"],
            "n_ionic_energy_records": len(triples),
            "initial_free_energy_eV": triples[0][0],
            "initial_enthalpy_eV": triples[0][1],
            "initial_enthalpy_minus_audited_signed_static_eV": (
                triples[0][1] - record["source_static_enthalpy_eV_per_cell"]
            ),
            "last_evaluated_free_energy_eV": triples[-1][0],
            "last_evaluated_PV_eV": triples[-1][2],
            "last_evaluated_enthalpy_eV": triples[-1][1],
            "enthalpy_drop_from_first_evaluated_step_eV": triples[-1][1] - triples[0][1],
            "lowest_evaluated_enthalpy_eV": min(item[1] for item in triples),
            "last_evaluated_volume_A3": float(last_atoms.get_volume()),
            "last_evaluated_max_atomic_force_eV_per_A": float(
                np.max(np.linalg.norm(last_atoms.get_forces(), axis=1))
            ),
            "contcar_is_last_evaluated_geometry": same_geometry(
                last_atoms, next_geometry, tolerance=2e-5,
            ),
            "contcar_volume_A3": float(next_geometry.get_volume()),
            "last_evaluated_GaN_coordination_2p4A": ga_n_coordination(last_atoms),
            "contcar_GaN_coordination_2p4A": ga_n_coordination(next_geometry),
            "input_sha256": hashes,
            "outcar_sha256": sha256(outcar), "contcar_sha256": sha256(contcar),
        })
    result = {
        "status": "GaN_native_basin_10step_pilots_raw_audited_not_phase_or_endpoint_certificates",
        "n_cases": len(results), "pressure_GPa": float(manifest["pressure_GPa"]),
        "PSTRESS_kbar": float(manifest["PSTRESS_kbar"]),
        "source_endpoint_GaN_coordination_2p4A": endpoint_coordination,
        "cases": results,
        "source_sha256": {
            "manifest": sha256(manifest_path), "signed_audit": sha256(signed_audit),
            "source_trajectory": sha256(path_trajectory),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "NSW=10 can stop before force/stress convergence; no phase identity is assigned.",
            "CONTCAR may be the next unevaluated geometry and must not be given the last OUTCAR energy.",
            "The VASP enthalpy field already includes P*V; use either it or free E plus P*V once.",
            "A 2.4 Å Ga--N coordination count is a topology screen calibrated on this path, not a converged phase label.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"],
                      "n_records": [item["n_ionic_energy_records"] for item in results],
                      "last_enthalpies": [item["last_evaluated_enthalpy_eV"] for item in results]}))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work-root", "signed-root", "signed-audit", "path-trajectory", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.work_root, args.signed_root, args.signed_audit,
          args.path_trajectory, args.output)


if __name__ == "__main__":
    main()
