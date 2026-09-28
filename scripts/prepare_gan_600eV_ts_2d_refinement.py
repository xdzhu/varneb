"""Stage the 12 missing static points of the GaN local 5x5 joint cut.

Reuse the original 600-eV center and twelve audited pilot points. This only
prepares frozen atom--strain structures; it neither changes VASP parameters
nor asserts a certified transition state or a whole-path surface.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coordinate_key(q_u: float, q_v: float) -> tuple[float, float]:
    return round(q_u, 8), round(q_v, 8)


def missing_coordinates(prior_coords: set[tuple[float, float]]) -> tuple[tuple[float, ...], tuple[float, ...], list[tuple[int, int, float, float]]]:
    q_u_values = (-0.02, -0.01, 0.0, 0.01, 0.02)
    q_v_values = (-0.0125, -0.00625, 0.0, 0.00625, 0.0125)
    full = {coordinate_key(q_u, q_v) for q_u in q_u_values for q_v in q_v_values}
    if len(prior_coords) != 13 or not prior_coords <= full:
        raise ValueError("previous center and twelve statics do not belong to the 5x5 grid")
    missing = [(i, j, q_u, q_v)
               for j, q_v in enumerate(q_v_values)
               for i, q_u in enumerate(q_u_values)
               if coordinate_key(q_u, q_v) not in prior_coords]
    if len(missing) != 12:
        raise ValueError(f"expected twelve new joint coordinates, got {len(missing)}")
    return q_u_values, q_v_values, missing


def prepare(center_work: Path, hessian_root: Path, prior_manifest_path: Path,
            pilot_audit_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    old = json.loads(prior_manifest_path.read_text(encoding="utf-8"))
    pilot = json.loads(pilot_audit_path.read_text(encoding="utf-8"))
    hessian_audit_path = hessian_root / "audit.json"
    hessian_npz = hessian_root / "joint_hessian.npz"
    center_outcar = center_work / "OUTCAR"
    if (old.get("status") != "inputs_finalized_no_DFT"
            or old.get("purpose") != "GaN_45p7_600eV_local_joint_mode_2D_frozen_pilot"
            or len(old.get("cases", [])) != 12
            or old.get("pressure_GPa") != 45.7
            or old.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or old.get("source_sha256", {}).get("center_OUTCAR") != sha256(center_outcar)
            or old.get("source_sha256", {}).get("hessian_audit") != sha256(hessian_audit_path)
            or old.get("source_sha256", {}).get("hessian_npz") != sha256(hessian_npz)
            or pilot.get("status") != "GaN_600eV_local_2D_frozen_enthalpy_pilot_raw_audited"
            or pilot.get("source_sha256", {}).get("manifest") != sha256(prior_manifest_path)
            or pilot.get("source_sha256", {}).get("center_OUTCAR") != sha256(center_outcar)
            or pilot.get("source_sha256", {}).get("hessian_audit") != sha256(hessian_audit_path)
            or pilot.get("source_sha256", {}).get("hessian_npz") != sha256(hessian_npz)
            or any(sha256(center_work / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("old 600-eV local cut or electronic-input provenance changed")
    center = read(center_outcar)
    chart = JointCurvatureCoordinates(center, float(old["axes"]["cell_scale_A"]))
    with np.load(hessian_npz, allow_pickle=False) as archive:
        eigenvalues = np.asarray(archive["eigenvalues"], dtype=float)
        modes = np.asarray(archive["eigenvectors"], dtype=float)
    if (eigenvalues.shape != (15,) or modes.shape != (18, 15)
            or not np.allclose(modes.T @ modes, np.eye(15), atol=1e-8)
            or not np.allclose(eigenvalues[:2], old["axes"]["eigenvalues_eV_per_A2"], atol=1e-8)
            or not eigenvalues[0] < 0 < eigenvalues[1]):
        raise ValueError("local joint-mode basis differs from audited pilot")
    prior_coords = {coordinate_key(float(item["q_u_A"]), float(item["q_v_A"]))
                    for item in old["cases"]}
    prior_coords.add((0.0, 0.0))
    q_u_values, q_v_values, missing = missing_coordinates(prior_coords)
    staged = []
    for i, j, q_u, q_v in missing:
        atoms = chart.displaced(q_u * modes[:, 0] + q_v * modes[:, 1])
        geometry = geometry_preflight(atoms)
        if geometry["empirical_near_symmetry_warning"]:
            raise ValueError(f"Bravais risk at joint coordinate ({i}, {j})")
        staged.append((i, j, q_u, q_v, atoms, geometry))
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    records = []
    for i, j, q_u, q_v, atoms, geometry in staged:
        name = f"grid5_u{i}_v{j}"
        case = output / "cases" / name
        case.mkdir()
        write(case / "POSCAR", atoms, format="vasp", direct=True, sort=False)
        if not same_geometry(atoms, read(case / "POSCAR", format="vasp")):
            raise ValueError(f"POSCAR roundtrip failed at {name}")
        for input_name in PRODUCTION_INPUT_SHA256:
            shutil.copy2(center_work / input_name, case / input_name)
        hashes = {input_name: sha256(case / input_name)
                  for input_name in ("POSCAR", *PRODUCTION_INPUT_SHA256)}
        if any(hashes[name] != digest
               for name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"electronic contract changed at {name}")
        (case / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        records.append({
            "name": name, "grid_index": [i, j], "q_u_A": q_u, "q_v_A": q_v,
            **geometry, "input_sha256": hashes,
        })
    result = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_local_joint_mode_5x5_frozen_refinement",
        "pressure_GPa": 45.7, "formula_units_per_cell": 2,
        "q_u_grid_A": q_u_values, "q_v_grid_A": q_v_values,
        "n_prior_coordinates": 13, "n_new_coordinates": 12,
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "axes": old["axes"],
        "source_sha256": {
            "center_OUTCAR": sha256(center_outcar),
            "hessian_audit": sha256(hessian_audit_path),
            "hessian_npz": sha256(hessian_npz),
            "prior_manifest": sha256(prior_manifest_path),
            "pilot_audit": sha256(pilot_audit_path),
            "preparer": sha256(Path(__file__)),
        },
        "cases": records,
        "limitations": [
            "Frozen local enthalpy cut only; neither a conditional PES nor a whole-path surface.",
            "Twelve new off-axis statics do not repair the existing energy-stress derivative mismatch.",
        ],
    }
    (output / "manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--center-work", type=Path, required=True)
    parser.add_argument("--hessian-root", type=Path, required=True)
    parser.add_argument("--prior-manifest", type=Path, required=True)
    parser.add_argument("--pilot-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.center_work, args.hessian_root, args.prior_manifest,
                     args.pilot_audit, args.output)
    print(json.dumps({"status": result["status"], "n_new_cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
