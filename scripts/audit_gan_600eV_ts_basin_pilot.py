"""Audit the 600-eV GaN signed native-VASP basin pilots after completion.

Only raw electronic/ionic output and input identity are certified here.  A
ten-step cap, coordination count, or changed cell is not an endpoint or TS
certificate; unevaluated CONTCAR geometry is never assigned an OUTCAR energy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_ts_basin_pilot import energy_triples, ga_n_coordination
from scripts.compare_gan_basin_endpoint import relative_displacement_metrics
from scripts.prepare_gan_600eV_ts_basin_pilot import CASES, sha256
from scripts.prepare_gan_600eV_ts_hessian import same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


EXPECTED_CHAIN_SHA256 = "952298c1830b293690b8fc4722649147cba236e50e149d07dae45ff75fc2bb7e"
ARCHIVED_INPUTS = (Path(__file__).resolve().parents[1]
                   / "benchmarks/numerical_integrity/gan_600eV_ts_basin_pilot_inputs_20260928.json")


def audit(work: Path, chain_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    archived_manifest = json.loads(ARCHIVED_INPUTS.read_text(encoding="utf-8"))
    if (manifest != archived_manifest
            or manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_signed_native_VASP_basin_10step_pilot"
            or manifest.get("encut_eV") != 600
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("PSTRESS_kbar") != 457.0
            or manifest.get("n_cases") != 2
            or manifest.get("source_sha256", {}).get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_ts_basin_pilot.py"))
            or sha256(chain_path) != EXPECTED_CHAIN_SHA256):
        raise ValueError("600-eV basin pilot source or final chain changed")
    cases = manifest["cases"]
    if tuple(record["name"] for record in cases) != CASES:
        raise ValueError("signed basin case order changed")
    reference = read(chain_path, index=":")
    if len(reference) != 29:
        raise ValueError("original GaN VASP chain is not 29 total images")
    endpoint_coordination = {
        "B4": ga_n_coordination(reference[0]),
        "B1": ga_n_coordination(reference[-1]),
    }
    if endpoint_coordination != {"B4": [4, 4], "B1": [6, 6]}:
        raise ValueError("source endpoint topology reference differs")
    pressure = 45.7 * GPa
    results = []
    for record in cases:
        directory = work / "cases" / record["name"]
        hashes = json.loads((directory / "sha256.inputs.json").read_text(encoding="utf-8"))
        if (hashes != record["input_sha256"]
                or any(sha256(directory / name) != digest for name, digest in hashes.items())
                or hashes["KPOINTS"] != PRODUCTION_INPUT_SHA256["KPOINTS"]
                or hashes["POTCAR"] != PRODUCTION_INPUT_SHA256["POTCAR"]
                or b"\r" in (directory / "INCAR").read_bytes()):
            raise ValueError(f"VASP pilot input hash or electronic contract changed: {record['name']}")
        outcar = directory / "OUTCAR"
        raw = outcar.read_text(encoding="utf-8", errors="replace")
        if ("General timing and accounting informations" not in raw
                or "PSTRESS=  457.0" not in raw):
            raise ValueError(f"VASP pilot raw output is incomplete or pressure differs: {outcar}")
        triples = energy_triples(raw)
        if raw.count("aborting loop because EDIFF is reached") != len(triples):
            raise ValueError(f"not every ionic step reached the electronic SCF criterion: {outcar}")
        final = read(outcar)
        if (len(final) != 4
                or final.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
                or not np.isfinite(final.get_forces()).all()
                or not np.isfinite(final.get_stress(voigt=False)).all()
                or not np.isclose(triples[-1][2], pressure * final.get_volume(), atol=2e-4, rtol=0)):
            raise ValueError(f"VASP pilot final E/PV, forces, stress or geometry differ: {outcar}")
        contcar = directory / "CONTCAR"
        next_geometry = read(contcar, format="vasp")
        seed = read(directory / "POSCAR", format="vasp")
        results.append({
            "case": record["name"],
            "q_u_A": record["q_u_A"],
            "n_ionic_energy_records": len(triples),
            "first_evaluated_enthalpy_eV_per_cell": triples[0][1],
            "last_evaluated_enthalpy_eV_per_cell": triples[-1][1],
            "enthalpy_change_from_first_evaluated_eV_per_cell": triples[-1][1] - triples[0][1],
            "last_evaluated_volume_A3": float(final.get_volume()),
            "last_evaluated_max_atomic_force_eV_per_A": float(
                np.max(np.linalg.norm(final.get_forces(), axis=1))
            ),
            "last_evaluated_GaN_coordination_2p4A": ga_n_coordination(final),
            "contcar_GaN_coordination_2p4A": ga_n_coordination(next_geometry),
            "contcar_is_last_evaluated_geometry": same_geometry(final, next_geometry),
            "last_evaluated_volume_difference_from_B4_A3": float(final.get_volume() - reference[0].get_volume()),
            "last_evaluated_volume_difference_from_B1_A3": float(final.get_volume() - reference[-1].get_volume()),
            "mapped_endpoint_distances_A": {
                phase: {
                    "seed": relative_displacement_metrics(seed, endpoint),
                    "last_evaluated": relative_displacement_metrics(final, endpoint),
                }
                for phase, endpoint in (("B4", reference[0]), ("B1", reference[-1]))
            },
            "input_sha256": hashes,
            "outcar_sha256": sha256(outcar),
            "contcar_sha256": sha256(contcar),
        })
    result = {
        "status": "GaN_600eV_basin_10step_pilots_raw_audited_not_phase_certificates",
        "n_cases": 2,
        "pressure_GPa": 45.7,
        "PSTRESS_kbar": 457.0,
        "n_images_source_chain": 29,
        "source_endpoint_GaN_coordination_2p4A": endpoint_coordination,
        "cases": results,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "archived_input_manifest_normalized_copy": sha256(ARCHIVED_INPUTS),
            "600eV_final_chain": sha256(chain_path),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "Ten native VASP ionic steps do not demonstrate a basin connection or phase convergence.",
            "The 2.4-A Ga-N coordination count is a topology screen, not a phase proof.",
            "CONTCAR can be an unevaluated next geometry after a finite NSW cap.",
            "VASP printed enthalpy already contains PSTRESS*V; no second PV term was added.",
            "The 600-eV energy/stress derivative mismatch remains a strict-TS certification limit.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--chain", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.work, args.chain, args.output)
    print(json.dumps({"status": result["status"], "n_cases": result["n_cases"],
                      "enthalpy_changes": [record["enthalpy_change_from_first_evaluated_eV_per_cell"]
                                           for record in result["cases"]]}))


if __name__ == "__main__":
    main()
