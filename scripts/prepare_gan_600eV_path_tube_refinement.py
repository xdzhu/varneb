"""Targeted same-contract GaN path-tube anchors after measured LOO failure."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.analyze_gan_600eV_path_tube_feasibility import transported_normals
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


REFINEMENT_ANCHORS = (1, 6, 8, 9, 16, 18, 19, 21, 22)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(trajectory: Path, feasibility_path: Path, first_audit_path: Path,
            transport_path: Path, screen_path: Path, hessian_npz: Path, source_static: Path,
            output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    first = json.loads(first_audit_path.read_text(encoding="utf-8"))
    transport = json.loads(transport_path.read_text(encoding="utf-8"))
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    if (feasibility.get("status") != "GaN_600eV_path_adapted_joint_normal_tube_geometry_only"
            or feasibility["source_sha256"]["trajectory"] != sha256(trajectory)
            or feasibility["source_sha256"]["hessian_npz"] != sha256(hessian_npz)
            or first.get("status") != "GaN_600eV_path_adapted_frozen_transverse_tube_raw_audited"
            or first["source_sha256"]["feasibility"] != sha256(feasibility_path)
            or first["maximum_leave_one_out_error_meV_per_GaN"] <= 0.10
            or transport.get("status") != "GaN_600eV_path_normal_transport_offline_comparison"
            or transport["source_sha256"]["trajectory"] != sha256(trajectory)
            or transport["source_sha256"]["hessian_npz"] != sha256(hessian_npz)
            or transport["recursive_parallel_transport"]["max_adjacent_slope_jump_eV_per_A"]
            <= transport["fixed_seed_projection"]["max_adjacent_slope_jump_eV_per_A"]
            or screen.get("status") != "GaN_600eV_offpath_refinement_geometry_only_no_DFT"
            or screen["source_sha256"]["trajectory"] != sha256(trajectory)
            or screen["source_sha256"]["hessian_npz"] != sha256(hessian_npz)
            or screen["source_sha256"]["feasibility"] != sha256(feasibility_path)
            or any(sha256(source_static / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("targeted refinement rationale or 600-eV source changed")
    if set(REFINEMENT_ANCHORS) & set(first["anchors"]):
        raise ValueError("refinement repeats an already computed path anchor")
    frames = read(trajectory, index=":")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        seed = np.asarray(archive["eigenvectors"][:, 1], dtype=float)
    scale = 3.3982714330050063
    normals, calculated_transport = transported_normals(frames, seed, scale)
    for key, value in calculated_transport.items():
        if not np.isclose(value, feasibility["transport"][key], atol=1e-10, rtol=0):
            raise ValueError("path-tube normal differs from the first batch")
    safe = {row["image_index"] for row in screen["rows"]
            if row["amplitude_A"] == 0.015 and row["both_sides_pass_empirical_screen"]}
    if not set(REFINEMENT_ANCHORS) <= safe:
        raise ValueError("a proposed anchor did not pass the recorded geometry screen")
    # Validate the whole proposed batch before creating any output directory.
    proposed = []
    for index in REFINEMENT_ANCHORS:
        chart = JointCurvatureCoordinates(frames[index], scale)
        for sign, side in ((-1, "minus"), (1, "plus")):
            atoms = chart.displaced(sign * 0.015 * normals[index])
            geometry = geometry_preflight(atoms)
            if geometry["empirical_near_symmetry_warning"]:
                raise ValueError(f"near-symmetry risk at refinement image {index} {side}")
            proposed.append((index, sign, side, atoms, geometry))
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    cases = []
    for index, sign, side, atoms, geometry in proposed:
        name = f"image_{index:02d}_q_{side}"
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
            "name": name, "image_index": index,
            "q_perp_A": sign * 0.015,
            "arc_fraction_s": feasibility["arc_fraction_s"][index],
            "path_center_enthalpy_eV_per_cell": feasibility["path_enthalpy_eV_per_cell"][index],
            **geometry, "input_sha256": hashes,
        })
    result = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_path_tube_error_targeted_refinement",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "anchors": list(REFINEMENT_ANCHORS),
        "q_perp_A": [-0.015, 0.015],
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "feasibility": sha256(feasibility_path),
            "first_tube_audit": sha256(first_audit_path),
            "transport_comparison": sha256(transport_path),
            "geometry_screen": sha256(screen_path),
            "hessian_npz": sha256(hessian_npz),
            "source_static_OUTCAR": sha256(source_static / "OUTCAR"),
            "preparer": sha256(Path(__file__)),
        },
        "cases": cases,
        "claim_limit": "Frozen path-tube refinement; no conditional relaxation or unique-MEP claim",
    }
    (output / "manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "feasibility", "first-audit", "transport", "screen",
                 "hessian-npz", "source-static", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.trajectory, args.feasibility, args.first_audit,
                     args.transport, args.screen, args.hessian_npz,
                     args.source_static, args.output)
    print(json.dumps({"status": result["status"], "n_new_cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
