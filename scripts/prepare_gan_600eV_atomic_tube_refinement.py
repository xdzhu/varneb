"""Targeted same-600-eV GaN atomic-tube anchors and transverse holdouts.

The four new path anchors address measured LOO failures at images 8 and 20.
The two signed half-step pairs independently test local q interpolation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.analyze_gan_600eV_atomic_tube_feasibility import (
    CELL_SCALE_A,
    atomic_normals,
)
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


NEW_ANCHORS = (6, 7, 19, 21)
HALFSTEP_HOLDOUTS = (8, 20)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(trajectory: Path, atomic_feasibility_path: Path,
            joint_feasibility_path: Path, transport_path: Path,
            first_audit_path: Path, hessian_npz: Path,
            source_static: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    atomic_feasibility = json.loads(atomic_feasibility_path.read_text(encoding="utf-8"))
    joint_feasibility = json.loads(joint_feasibility_path.read_text(encoding="utf-8"))
    transport = json.loads(transport_path.read_text(encoding="utf-8"))
    first = json.loads(first_audit_path.read_text(encoding="utf-8"))
    worst = {(item["image_index"], item["side"]): item["absolute_error_meV_per_GaN"]
             for item in first.get("leave_one_anchor_out_excess_enthalpy_errors", [])}
    if (atomic_feasibility.get("status")
            != "GaN_600eV_atomic_transverse_tube_geometry_only_no_DFT"
            or joint_feasibility.get("status")
            != "GaN_600eV_path_adapted_joint_normal_tube_geometry_only"
            or transport.get("status")
            != "GaN_600eV_atomic_transverse_transport_offline_comparison"
            or first.get("status") != "GaN_600eV_central_atomic_transverse_tube_raw_audited"
            or first.get("central_contour_gate_pass") is not False
            or first.get("maximum_leave_one_out_error_meV_per_GaN", 0) <= 1.0
            or worst.get((8, "plus"), 0) <= 1.0
            or worst.get((20, "minus"), 0) <= 1.0
            or set(NEW_ANCHORS) & set(first["anchors"])
            or not set(HALFSTEP_HOLDOUTS) <= set(first["anchors"])
            or atomic_feasibility["source_sha256"]["trajectory"] != sha256(trajectory)
            or atomic_feasibility["source_sha256"]["hessian_npz"] != sha256(hessian_npz)
            or joint_feasibility["source_sha256"]["trajectory"] != sha256(trajectory)
            or first["source_sha256"]["trajectory"] != sha256(trajectory)
            or first["source_sha256"]["atomic_feasibility"]
            != sha256(atomic_feasibility_path)
            or first["source_sha256"]["joint_feasibility"]
            != sha256(joint_feasibility_path)
            or first["source_sha256"]["transport_comparison"]
            != sha256(transport_path)
            or any(sha256(source_static / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("targeted refinement rationale or same-contract sources changed")
    frames = read(trajectory, index=":")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    normals, metrics = atomic_normals(frames, seed)
    if len(frames) != 29 or any(not np.isclose(value, atomic_feasibility["transport"][key], atol=1e-10, rtol=0)
                                for key, value in metrics.items()):
        raise ValueError("atomic normal changed")
    requested = [(index, 0.05, "interpolation_anchor") for index in NEW_ANCHORS]
    requested += [(index, 0.025, "transverse_holdout") for index in HALFSTEP_HOLDOUTS]
    proposed = []
    for index, amplitude, role in requested:
        chart = JointCurvatureCoordinates(frames[index], CELL_SCALE_A)
        for sign, side in ((-1, "minus"), (1, "plus")):
            atoms = chart.displaced(np.r_[sign * amplitude * normals[index], np.zeros(6)])
            geometry = geometry_preflight(atoms)
            if not np.allclose(atoms.cell.array, frames[index].cell.array,
                               atol=1e-10, rtol=0):
                raise ValueError(f"atomic-only chart changed cell at image {index} {side}")
            proposed.append((index, sign, side, amplitude, role, atoms, geometry))
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    cases = []
    for index, sign, side, amplitude, role, atoms, geometry in proposed:
        suffix = "" if role == "interpolation_anchor" else "_half"
        name = f"image_{index:02d}{suffix}_q_{side}"
        target = output / "cases" / name
        target.mkdir()
        write(target / "POSCAR", atoms, format="vasp", direct=True, sort=False)
        if not same_geometry(atoms, read(target / "POSCAR", format="vasp")):
            raise ValueError(f"POSCAR roundtrip failed at {name}")
        for filename in PRODUCTION_INPUT_SHA256:
            shutil.copy2(source_static / filename, target / filename)
        hashes = {filename: sha256(target / filename)
                  for filename in ("POSCAR", *PRODUCTION_INPUT_SHA256)}
        if any(hashes[filename] != digest
               for filename, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"electronic contract changed at {name}")
        (target / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        cases.append({
            "name": name,
            "role": role,
            "image_index": index,
            "q_atom_A": sign * amplitude,
            "arc_fraction_s": joint_feasibility["arc_fraction_s"][index],
            "path_center_enthalpy_eV_per_cell":
                joint_feasibility["path_enthalpy_eV_per_cell"][index],
            "cell_unchanged_from_audited_path": True,
            **geometry,
            "input_sha256": hashes,
        })
    manifest = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_atomic_only_central_tube_targeted_refinement",
        "claim_limit": "Four new path anchors and two signed q-halfstep holdouts; central segment only",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "new_anchors": list(NEW_ANCHORS),
        "transverse_holdout_anchors": list(HALFSTEP_HOLDOUTS),
        "anchor_q_A": [-0.05, 0.05],
        "holdout_q_A": [-0.025, 0.025],
        "predeclared_LOO_gate_meV_per_GaN": 1.0,
        "predeclared_halfstep_curvature_relative_gate": 0.05,
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "atomic_feasibility": sha256(atomic_feasibility_path),
            "joint_feasibility": sha256(joint_feasibility_path),
            "transport_comparison": sha256(transport_path),
            "first_audit": sha256(first_audit_path),
            "hessian_npz": sha256(hessian_npz),
            "source_static_OUTCAR": sha256(source_static / "OUTCAR"),
            "preparer": sha256(Path(__file__)),
        },
        "cases": cases,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "atomic-feasibility", "joint-feasibility",
                 "transport", "first-audit", "hessian-npz", "source-static", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.trajectory, args.atomic_feasibility,
                     args.joint_feasibility, args.transport, args.first_audit,
                     args.hessian_npz, args.source_static, args.output)
    print(json.dumps({"status": result["status"], "cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
