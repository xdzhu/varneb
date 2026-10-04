"""Audit both completed GaN fixed-center halves against their raw VASP outputs.

Run on the calculation host. The report records hashes, never POTCAR contents.
The shared center remains a candidate rather than a certified saddle.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.units import GPa

from vcneb.vasp import vasp_input_fingerprints


PRESSURE_GPA = 45.7
N_GAN = 2


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _same_geometry(a, b, tolerance: float = 3e-5) -> bool:
    if a.get_chemical_symbols() != b.get_chemical_symbols():
        return False
    if np.max(np.abs(a.cell.array - b.cell.array)) > tolerance:
        return False
    delta = a.get_scaled_positions(wrap=False) - b.get_scaled_positions(wrap=False)
    delta -= np.rint(delta)
    return bool(np.max(np.linalg.norm(delta @ b.cell.array, axis=1)) < tolerance)


def _segment(root: Path, name: str, manifest: dict) -> tuple[dict, list]:
    data = manifest["segments"][name]
    run = root / name / "run"
    summary_path = run / "vcneb_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    audit = json.loads((run / "vcneb_audit.json").read_text(encoding="utf-8"))
    n = data["total_images"]
    if (summary["status"] != "completed" or not summary["converged"]
            or summary["n_images"] != n
            or summary["final_max_generalized_force_eV_per_A"] > 0.10
            or summary["endpoint_evaluation_policy"] != "fixed_cached_once"
            or audit["status"] != "ok"):
        raise ValueError(f"{name}: manager convergence or geometry audit failed")
    with Trajectory(str(run / "vcneb.traj"), "r") as trajectory:
        if len(trajectory) < n or len(trajectory) % n:
            raise ValueError(f"{name}: incomplete final chain")
        chain = [trajectory[i] for i in range(len(trajectory) - n, len(trajectory))]
    if len(chain) != n:
        raise ValueError(f"{name}: incomplete final chain")
    if [im.get_chemical_symbols() for im in chain] != [chain[0].get_chemical_symbols()] * n:
        raise ValueError(f"{name}: atom order varies")
    params = summary["calculator_parameters"]
    if (params["encut"] != 600 or params["kpts"] != [8, 8, 6]
            or params["isym"] != -1 or params["symprec"] != 1e-4
            or params["ediff"] != 1e-7 or params["setups"] != {"Ga": "_d", "N": ""}
            or summary["optimizer"] != "FIRE" or summary["climbing_image_requested"]):
        raise ValueError(f"{name}: electronic or NEB contract changed")
    reference_potcar = summary["licensed_input_fingerprints"]["POTCAR"]["sha256"]
    generated_contract = None
    records = []
    for index, image in enumerate(chain):
        e = float(image.get_potential_energy())
        f = np.asarray(image.get_forces(), dtype=float)
        s = np.asarray(image.get_stress(voigt=True), dtype=float)
        h = e + PRESSURE_GPA * GPa * image.get_volume()
        if (not np.isfinite([e, h]).all() or not np.isfinite(f).all()
                or not np.isfinite(s).all()
                or abs(h - summary["image_enthalpies_eV"][index]) > 3e-4):
            raise ValueError(f"{name} image {index}: trajectory properties incomplete")
        record = {"segment": name, "index": index, "energy_eV_cell": e,
                  "enthalpy_eV_cell": h, "volume_A3": float(image.get_volume()),
                  "max_force_eV_A": float(np.linalg.norm(f, axis=1).max())}
        if 0 < index < n - 1:
            image_dir = run / f"{index:02d}"
            raw_outcar = image_dir / "OUTCAR"
            text = raw_outcar.read_text(encoding="utf-8", errors="replace")
            if ("aborting loop because EDIFF is reached" not in text
                    or "General timing and accounting informations for this job" not in text):
                raise ValueError(f"{name} image {index}: incomplete SCF/OUTCAR")
            raw = read(raw_outcar)
            if (not _same_geometry(raw, image)
                    or abs(raw.get_potential_energy() - e) > 3e-4
                    or np.max(np.abs(raw.get_forces() - f)) > 3e-4
                    or np.max(np.abs(raw.get_stress(voigt=True) - s)) > 2e-5):
                raise ValueError(f"{name} image {index}: raw/trajectory mismatch")
            inputs = vasp_input_fingerprints(image_dir)
            contract = json.loads((image_dir / "vasp_input_contract.json").read_text(encoding="utf-8"))
            actual_hashes = {key: value["sha256"] for key, value in inputs.items()}
            if (actual_hashes != {key: contract["input_sha256"][key] for key in actual_hashes}
                    or actual_hashes["POTCAR"] != reference_potcar
                    or contract["effective_symmetry"] != {"isym": -1, "symprec": 1e-4}):
                raise ValueError(f"{name} image {index}: input hash or symmetry contract changed")
            effective = (actual_hashes["INCAR"], actual_hashes["KPOINTS"],
                         actual_hashes["POTCAR"], contract["parameter_sha256"])
            if generated_contract is None:
                generated_contract = effective
            elif effective != generated_contract:
                raise ValueError(f"{name} image {index}: generated static inputs varied")
            record["raw_OUTCAR_sha256"] = sha256(raw_outcar)
        records.append(record)
    report = {
        "job_segment": name,
        "n_images": n,
        "n_raw_interior_OUTCARs_audited": n - 2,
        "final_max_generalized_force_eV_per_A": summary["final_max_generalized_force_eV_per_A"],
        "summary_sha256": sha256(summary_path),
        "trajectory_sha256": sha256(run / "vcneb.traj"),
        "manager_audit_sha256": sha256(run / "vcneb_audit.json"),
        "shared_center_structure_sha256": (summary["endpoint_structures"]["final" if name == "left" else "initial"]["sha256"]),
        "generated_static_input_sha256": dict(zip(
            ("INCAR", "KPOINTS", "POTCAR", "effective_parameters"), generated_contract,
        )),
        "records": records,
    }
    return report, chain


def audit(root: Path) -> dict:
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest["pressure_GPa"] != PRESSURE_GPA
            or manifest["status"] != "two_segment_preflight_passed_no_DFT"):
        raise ValueError("invalid preparation manifest")
    left, left_chain = _segment(root, "left", manifest)
    right, right_chain = _segment(root, "right", manifest)
    if (left["shared_center_structure_sha256"] != right["shared_center_structure_sha256"]
            or not _same_geometry(left_chain[-1], right_chain[0])
            or abs(left["records"][-1]["enthalpy_eV_cell"]
                   - right["records"][0]["enthalpy_eV_cell"]) > 1e-8):
        raise ValueError("shared center is inconsistent")
    if left["generated_static_input_sha256"] != right["generated_static_input_sha256"]:
        raise ValueError("left and right generated electronic inputs differ")
    stitched = left["records"] + right["records"][1:]
    peak = max(range(len(stitched)), key=lambda i: stitched[i]["enthalpy_eV_cell"])
    if peak != len(left["records"]) - 1:
        raise ValueError("interior image overtakes shared center")
    h = np.array([r["enthalpy_eV_cell"] for r in stitched])
    return {
        "status": "raw_per_image_audit_passed",
        "claim_limit": "fixed center is an index-one saddle candidate, not a certified stationary TS",
        "pressure_GPa": PRESSURE_GPA,
        "n_total_stitched_images": len(stitched),
        "n_total_raw_interior_OUTCARs_audited": left["n_raw_interior_OUTCARs_audited"] + right["n_raw_interior_OUTCARs_audited"],
        "peak_stitched_index": peak,
        "forward_barrier_meV_per_GaN": float(1000 * (h[peak] - h[0]) / N_GAN),
        "reverse_barrier_meV_per_GaN": float(1000 * (h[peak] - h[-1]) / N_GAN),
        "manifest_sha256": sha256(manifest_path),
        "segments": {"left": left, "right": right},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.case)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "n_total_raw_interior_OUTCARs_audited", "forward_barrier_meV_per_GaN",
        "reverse_barrier_meV_per_GaN")}))


if __name__ == "__main__":
    main()
