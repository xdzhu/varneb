"""Stage only 56 unmeasured points in the local GaN 9x9 joint-mode cut.

The existing 25 DFT coordinates, production VASP inputs, Hessian axes, and
45.7-GPa enthalpy reference are immutable inputs. No DFT runs in this script.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.prepare_gan_600eV_ts_2d_refinement import coordinate_key, sha256
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


def missing_9x9(prior: set[tuple[float, float]]) -> tuple[list[float], list[float], list[tuple[int, int, float, float]]]:
    u = [round(float(x), 8) for x in np.linspace(-0.02, 0.02, 9)]
    v = [round(float(x), 8) for x in np.linspace(-0.0125, 0.0125, 9)]
    full = {coordinate_key(x, y) for x in u for y in v}
    if len(prior) != 25 or not prior <= full:
        raise ValueError("the audited 5x5 grid is not a subset of the 9x9 grid")
    missing = [(i, j, x, y) for j, y in enumerate(v) for i, x in enumerate(u)
               if coordinate_key(x, y) not in prior]
    if len(missing) != 56:
        raise ValueError("the 9x9 increment must have exactly 56 fresh nodes")
    return u, v, missing


def prepare(center_work: Path, hessian_root: Path, pilot_manifest_path: Path,
            pilot_path: Path, refined5_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    pilot = json.loads(pilot_path.read_text(encoding="utf-8"))
    pilot_manifest = json.loads(pilot_manifest_path.read_text(encoding="utf-8"))
    refined5 = json.loads(refined5_path.read_text(encoding="utf-8"))
    hessian_audit = hessian_root / "audit.json"
    hessian_npz = hessian_root / "joint_hessian.npz"
    center_outcar = center_work / "OUTCAR"
    if (pilot.get("status") != "GaN_600eV_local_2D_frozen_enthalpy_pilot_raw_audited"
            or refined5.get("status") != "GaN_600eV_local_joint_5x5_refinement_raw_audited"
            or refined5.get("source_sha256", {}).get("original_pilot_audit") != sha256(pilot_path)
            or pilot.get("source_sha256", {}).get("manifest") != sha256(pilot_manifest_path)
            or pilot_manifest.get("status") != "inputs_finalized_no_DFT"
            or pilot.get("source_sha256", {}).get("center_OUTCAR") != sha256(center_outcar)
            or pilot.get("source_sha256", {}).get("hessian_audit") != sha256(hessian_audit)
            or pilot.get("source_sha256", {}).get("hessian_npz") != sha256(hessian_npz)
            or any(sha256(center_work / key) != digest
                   for key, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("9x9 staging cannot verify original 600-eV input and 5x5 audit")
    prior = {coordinate_key(float(item["q_u_A"]), float(item["q_v_A"]))
             for item in [*pilot["cases"], *refined5["cases"]]}
    prior.add((0.0, 0.0))
    u, v, missing = missing_9x9(prior)
    center = read(center_outcar)
    chart = JointCurvatureCoordinates(center, float(pilot_manifest["axes"]["cell_scale_A"]))
    with np.load(hessian_npz, allow_pickle=False) as archive:
        eig = np.asarray(archive["eigenvalues"], dtype=float)
        modes = np.asarray(archive["eigenvectors"], dtype=float)
    if (eig.shape != (15,) or modes.shape != (18, 15)
            or not eig[0] < 0 < eig[1]
            or not np.allclose(modes.T @ modes, np.eye(15), atol=1e-8)
            or not np.allclose(eig[:2], pilot["hessian_eigenvalues_uv_eV_per_A2"], atol=1e-8)):
        raise ValueError("joint axes differ from the audited 5x5 cut")
    staged = []
    for i, j, q_u, q_v in missing:
        atoms = chart.displaced(q_u * modes[:, 0] + q_v * modes[:, 1])
        geometry = geometry_preflight(atoms)
        if geometry["empirical_near_symmetry_warning"]:
            raise ValueError(f"Bravais-risk geometry at ({i}, {j})")
        staged.append((i, j, q_u, q_v, atoms, geometry))
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    records = []
    for i, j, q_u, q_v, atoms, geometry in staged:
        name = f"grid9_u{i}_v{j}"
        case = output / "cases" / name
        case.mkdir()
        write(case / "POSCAR", atoms, format="vasp", direct=True, sort=False)
        if not same_geometry(atoms, read(case / "POSCAR", format="vasp")):
            raise ValueError(f"POSCAR roundtrip changed geometry: {name}")
        for key in PRODUCTION_INPUT_SHA256:
            shutil.copy2(center_work / key, case / key)
        hashes = {key: sha256(case / key) for key in ("POSCAR", *PRODUCTION_INPUT_SHA256)}
        if any(hashes[key] != digest for key, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"VASP input contract changed at {name}")
        (case / "sha256.inputs.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")
        records.append({"name": name, "grid_index": [i, j], "q_u_A": q_u,
                        "q_v_A": q_v, **geometry, "input_sha256": hashes})
    result = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_local_joint_mode_9x9_frozen_refinement",
        "pressure_GPa": 45.7, "formula_units_per_cell": 2,
        "q_u_grid_A": u, "q_v_grid_A": v,
        "n_prior_coordinates": 25, "n_new_coordinates": 56,
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "axes": pilot_manifest["axes"],
        "source_sha256": {
            "center_OUTCAR": sha256(center_outcar),
            "hessian_audit": sha256(hessian_audit),
            "hessian_npz": sha256(hessian_npz),
            "pilot_audit": sha256(pilot_path),
            "pilot_manifest": sha256(pilot_manifest_path),
            "refinement5_audit": sha256(refined5_path),
            "preparer": sha256(Path(__file__)),
        },
        "cases": records,
        "claim_limit": "Frozen local joint enthalpy cut only; no full-path or certified TS claim.",
    }
    (output / "manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("center-work", "hessian-root", "pilot-manifest", "pilot-audit", "refinement5-audit", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.center_work, args.hessian_root, args.pilot_manifest, args.pilot_audit,
                     args.refinement5_audit, args.output)
    print(json.dumps({"status": result["status"], "n_new_cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
