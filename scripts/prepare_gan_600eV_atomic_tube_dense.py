"""Stage only missing points of the GaN 600-eV central atomic tube.

The centerline is the audited 29-image variable-cell VCNEB chain.  Each
transverse displacement changes atomic positions only and keeps that image's
cell fixed.  The resulting chart covers images 5--22, not both endpoint
neighborhoods, where the transported normal loses identity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.analyze_gan_600eV_atomic_tube_feasibility import CELL_SCALE_A, atomic_normals
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


CENTRAL_IMAGES = tuple(range(5, 23))
TRANSVERSE_Q_A = (-0.05, -0.025, 0.0, 0.025, 0.05)
PROSPECTIVE_HOLDOUT_IMAGES = (10, 16)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def missing_points(first_audit: dict, refinement_audit: dict) -> list[tuple[int, float]]:
    """Leave all previously completed ±q cases out of the new DFT array."""
    prior = {
        (int(case["image_index"]), round(float(case["q_atom_A"]), 6))
        for case in first_audit["cases"] + refinement_audit["new_cases"]
    }
    if (len(prior) != 28
            or set(first_audit["anchors"]) != {5, 8, 11, 14, 17, 18, 20, 22}
            or set(refinement_audit["anchors"]) != {5, 6, 7, 8, 11, 14, 17, 18, 19, 20, 21, 22}
            or any(index not in CENTRAL_IMAGES or q not in TRANSVERSE_Q_A for index, q in prior)):
        raise ValueError("prior central-tube point set changed; refuse duplicate submission")
    requested = [
        (index, q)
        for index in CENTRAL_IMAGES
        for q in TRANSVERSE_Q_A
        if q != 0.0 and (index, q) not in prior
    ]
    if len(requested) != 44:
        raise ValueError("expected exactly 44 new central-grid statics")
    return requested


def prepare(
    trajectory: Path,
    atomic_feasibility_path: Path,
    joint_feasibility_path: Path,
    transport_path: Path,
    first_audit_path: Path,
    refinement_audit_path: Path,
    hessian_npz: Path,
    source_static: Path,
    output: Path,
) -> dict:
    if output.exists():
        raise FileExistsError(output)
    atomic = json.loads(atomic_feasibility_path.read_text(encoding="utf-8"))
    joint = json.loads(joint_feasibility_path.read_text(encoding="utf-8"))
    transport = json.loads(transport_path.read_text(encoding="utf-8"))
    first = json.loads(first_audit_path.read_text(encoding="utf-8"))
    refinement = json.loads(refinement_audit_path.read_text(encoding="utf-8"))
    source_hashes = {
        "trajectory": sha256(trajectory),
        "atomic_feasibility": sha256(atomic_feasibility_path),
        "joint_feasibility": sha256(joint_feasibility_path),
        "transport_comparison": sha256(transport_path),
        "hessian_npz": sha256(hessian_npz),
    }
    if (atomic.get("status") != "GaN_600eV_atomic_transverse_tube_geometry_only_no_DFT"
            or joint.get("status") != "GaN_600eV_path_adapted_joint_normal_tube_geometry_only"
            or transport.get("status") != "GaN_600eV_atomic_transverse_transport_offline_comparison"
            or first.get("status") != "GaN_600eV_central_atomic_transverse_tube_raw_audited"
            or refinement.get("status") != "GaN_600eV_central_atomic_tube_targeted_refinement_raw_audited"
            or refinement.get("central_contour_gate_pass") is not False
            or refinement.get("predeclared_LOO_gate_meV_per_GaN") != 1.0
            or transport["fixed_seed_projection"]["minimum_central_5_to_22_overlap"] <= 0.95
            or any(first["source_sha256"].get(key) != digest
                   or refinement["source_sha256"].get(key) != digest
                   for key, digest in source_hashes.items())
            or any(sha256(source_static / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("audited central chart or 600-eV production contract changed")
    requested = missing_points(first, refinement)
    frames = read(trajectory, index=":")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    normals, metrics = atomic_normals(frames, seed)
    if (len(frames) != 29
            or any(not np.isclose(value, atomic["transport"][key], atol=1e-10, rtol=0)
                   for key, value in metrics.items())):
        raise ValueError("atomic normal transport changed")
    proposed = []
    for index, q in requested:
        chart = JointCurvatureCoordinates(frames[index], CELL_SCALE_A)
        atoms = chart.displaced(np.r_[q * normals[index], np.zeros(6)])
        geometry = geometry_preflight(atoms)
        if not np.allclose(atoms.cell.array, frames[index].cell.array, atol=1e-10, rtol=0):
            raise ValueError(f"atomic transverse displacement changed cell at image {index}")
        proposed.append((index, q, atoms, geometry))
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    cases = []
    for index, q, atoms, geometry in proposed:
        side = "minus" if q < 0 else "plus"
        width = "025" if abs(q) == 0.025 else "050"
        name = f"image_{index:02d}_q{width}_{side}"
        case_dir = output / "cases" / name
        case_dir.mkdir()
        write(case_dir / "POSCAR", atoms, format="vasp", direct=True, sort=False)
        if not same_geometry(atoms, read(case_dir / "POSCAR", format="vasp")):
            raise ValueError(f"POSCAR roundtrip failed at {name}")
        for filename in PRODUCTION_INPUT_SHA256:
            shutil.copy2(source_static / filename, case_dir / filename)
        hashes = {
            filename: sha256(case_dir / filename)
            for filename in ("POSCAR", *PRODUCTION_INPUT_SHA256)
        }
        if any(hashes[filename] != digest
               for filename, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"600-eV input contract changed at {name}")
        (case_dir / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        cases.append({
            "name": name,
            "image_index": index,
            "q_atom_A": q,
            "role": "prospective_s_holdout" if index in PROSPECTIVE_HOLDOUT_IMAGES
                    else "dense_grid_training",
            "arc_fraction_s": joint["arc_fraction_s"][index],
            "path_center_enthalpy_eV_per_cell": joint["path_enthalpy_eV_per_cell"][index],
            "cell_unchanged_from_audited_path": True,
            **geometry,
            "input_sha256": hashes,
        })
    manifest = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_atomic_only_central_dense_2D_grid",
        "claim_limit": "Images 5-22 only; frozen atomic transverse coordinate, path s variable-cell",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "central_image_indices": list(CENTRAL_IMAGES),
        "q_atom_A": list(TRANSVERSE_Q_A),
        "prospective_s_holdout_images": list(PROSPECTIVE_HOLDOUT_IMAGES),
        "predeclared_s_holdout_gate_meV_per_GaN": 1.0,
        "n_reused_offpath_statics": 28,
        "n_reused_path_centers": 18,
        "n_new_statics": len(cases),
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "source_sha256": {
            **source_hashes,
            "first_audit": sha256(first_audit_path),
            "refinement_audit": sha256(refinement_audit_path),
            "source_static_OUTCAR": sha256(source_static / "OUTCAR"),
            "preparer": sha256(Path(__file__)),
        },
        "cases": cases,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "atomic-feasibility", "joint-feasibility", "transport",
                 "first-audit", "refinement-audit", "hessian-npz", "source-static", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(
        args.trajectory, args.atomic_feasibility, args.joint_feasibility,
        args.transport, args.first_audit, args.refinement_audit,
        args.hessian_npz, args.source_static, args.output,
    )
    print(json.dumps({"status": result["status"], "n_new_statics": len(result["cases"])}))


if __name__ == "__main__":
    main()
