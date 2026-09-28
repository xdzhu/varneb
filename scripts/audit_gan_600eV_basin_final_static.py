"""Audit a same-600-eV static check of a completed native GaN basin geometry."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_ts_basin_pilot import ga_n_coordination
from scripts.audit_gan_ts_basin_followups import residual_stress_kbar
from scripts.compare_gan_basin_endpoint import relative_displacement_metrics
from scripts.prepare_gan_600eV_basin_final_static import PURPOSE
from scripts.prepare_gan_600eV_ts_hessian import same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256, sha256


def audit(work: Path, relax_case: Path, reference_chain: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    case = work / "case"
    inputs = json.loads((case / "sha256.inputs.json").read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != PURPOSE
            or manifest.get("case") != relax_case.name
            or manifest.get("pressure_for_postprocessing_GPa") != 45.7
            or manifest.get("electronic_contract", {}).get("static_input_sha256")
            != PRODUCTION_INPUT_SHA256
            or manifest.get("source_sha256", {}).get("native_OUTCAR")
            != sha256(relax_case / "OUTCAR")
            or manifest.get("source_sha256", {}).get("native_CONTCAR")
            != sha256(relax_case / "CONTCAR")
            or inputs != manifest.get("input_sha256")
            or any(sha256(case / key) != digest for key, digest in inputs.items())):
        raise ValueError("static inputs or native relaxation provenance differ")
    raw = (case / "OUTCAR").read_text(encoding="utf-8", errors="replace")
    if ("General timing and accounting informations" not in raw
            or raw.count("aborting loop because EDIFF is reached") != 1):
        raise ValueError("static VASP output incomplete or electronic SCF not converged")
    static = read(case / "OUTCAR")
    original = read(case / "POSCAR", format="vasp")
    native = read(relax_case / "OUTCAR")
    reference = read(reference_chain, index=":")
    if (len(reference) != 29 or ga_n_coordination(reference[0]) != [4, 4]
            or ga_n_coordination(reference[-1]) != [6, 6]
            or static.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
            or not same_geometry(static, original)
            or not same_geometry(static, native)
            or not np.isfinite(static.get_forces()).all()
            or not np.isfinite(static.get_stress(voigt=False)).all()):
        raise ValueError("static geometry/force/stress or phase reference invalid")
    pressure = 45.7 * GPa
    H_static = float(static.get_potential_energy() + pressure * static.get_volume())
    H_native = float(native.get_potential_energy() + pressure * native.get_volume())
    H_B4 = float(reference[0].get_potential_energy() + pressure * reference[0].get_volume())
    H_B1 = float(reference[-1].get_potential_energy() + pressure * reference[-1].get_volume())
    result = {
        "status": "GaN_600eV_relaxed_geometry_static_raw_audited_not_TS_certificate",
        "case": manifest["case"], "pressure_GPa": 45.7, "ENCUT_eV": 600,
        "static_enthalpy_eV_per_cell": H_static,
        "native_relax_enthalpy_eV_per_cell": H_native,
        "static_minus_native_enthalpy_meV_per_cell": 1000.0 * (H_static - H_native),
        "static_minus_reference_B4_enthalpy_meV_per_cell": 1000.0 * (H_static - H_B4),
        "static_minus_reference_B1_enthalpy_meV_per_cell": 1000.0 * (H_static - H_B1),
        "static_max_atomic_force_eV_per_A": float(np.linalg.norm(static.get_forces(), axis=1).max()),
        "static_stress_residual_kbar": residual_stress_kbar(static.get_stress(voigt=False), 45.7),
        "static_force_pass_0p02eV_per_A": float(np.linalg.norm(static.get_forces(), axis=1).max()) <= 0.02,
        "static_stress_pass_2kbar": residual_stress_kbar(static.get_stress(voigt=False), 45.7) <= 2.0,
        "static_GaN_coordination_2p4A": ga_n_coordination(static),
        "mapped_endpoint_distances_A": {
            phase: relative_displacement_metrics(static, endpoint)
            for phase, endpoint in (("B4", reference[0]), ("B1", reference[-1]))
        },
        "source_sha256": {"manifest": sha256(manifest_path),
                          "static_OUTCAR": sha256(case / "OUTCAR"),
                          "native_OUTCAR": sha256(relax_case / "OUTCAR"),
                          "reference_chain": sha256(reference_chain),
                          "auditor": sha256(Path(__file__))},
        "limitations": [
            "A static calculation at one relaxed geometry does not prove a unique phase basin.",
            "Native variable-cell and fresh static energies may differ from basis history; compare only static enthalpies on one electronic contract.",
            "The endpoint stress acceptance threshold is 2 kbar; it is not a DFT input parameter.",
            "A B4/B1 endpoint match and low force/stress would still not certify an index-one TS.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "relax-case", "reference-chain", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.work, args.relax_case, args.reference_chain, args.output)
    print(json.dumps({key: result[key] for key in (
        "case", "static_minus_native_enthalpy_meV_per_cell",
        "static_max_atomic_force_eV_per_A", "static_stress_residual_kbar",
        "static_GaN_coordination_2p4A")}))


if __name__ == "__main__":
    main()
