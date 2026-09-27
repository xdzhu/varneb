"""Stage GaN 45.7-GPa path-adapted transverse statics at the 600-eV contract.

Nine audited VCNEB anchors provide q=0 values. This stages ±0.015-A
transverse displacements at those anchors; it does not relax orthogonal
coordinates or claim a global two-phonon plane.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.analyze_gan_600eV_path_tube_feasibility import (
    ANCHORS, transported_normals,
)
from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(trajectory: Path, feasibility_path: Path, hessian_root: Path,
            source_static: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    hessian_npz = hessian_root / "joint_hessian.npz"
    hessian_audit = hessian_root / "audit.json"
    if (feasibility.get("status") != "GaN_600eV_path_adapted_joint_normal_tube_geometry_only"
            or feasibility.get("largest_preflighted_amplitude_A") != 0.015
            or feasibility.get("anchors") != list(ANCHORS)
            or feasibility.get("source_sha256", {}).get("trajectory") != sha256(trajectory)
            or feasibility.get("source_sha256", {}).get("hessian_npz") != sha256(hessian_npz)
            or feasibility.get("source_sha256", {}).get("hessian_audit") != sha256(hessian_audit)
            or feasibility.get("source_sha256", {}).get("analyzer")
            != sha256(Path(__file__).with_name("analyze_gan_600eV_path_tube_feasibility.py"))
            or any(sha256(source_static / name) != digest
                   for name, digest in PRODUCTION_INPUT_SHA256.items())):
        raise ValueError("600-eV path-tube feasibility or input contract changed")
    frames = read(trajectory, index=":")
    with np.load(hessian_npz, allow_pickle=False) as archive:
        modes = np.asarray(archive["eigenvectors"], dtype=float)
    scale = 3.3982714330050063
    normals, transport = transported_normals(frames, modes[:, 1], scale)
    for key, value in transport.items():
        if not np.isclose(value, feasibility["transport"][key], rtol=0, atol=1e-10):
            raise ValueError("transported transverse mode changed from preflight")
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    records = []
    for index in ANCHORS:
        chart = JointCurvatureCoordinates(frames[index], scale)
        for sign, side in ((-1, "minus"), (1, "plus")):
            displaced = chart.displaced(sign * 0.015 * normals[index])
            geometry = geometry_preflight(displaced)
            if geometry["empirical_near_symmetry_warning"]:
                raise ValueError(f"near-symmetry risk for image {index} {side}")
            name = f"image_{index:02d}_q_{side}"
            case = output / "cases" / name
            case.mkdir()
            write(case / "POSCAR", displaced, format="vasp", direct=True, sort=False)
            if not same_geometry(displaced, read(case / "POSCAR", format="vasp")):
                raise ValueError(f"off-path POSCAR changed at {name}")
            for filename in PRODUCTION_INPUT_SHA256:
                shutil.copy2(source_static / filename, case / filename)
            hashes = {filename: sha256(case / filename)
                      for filename in ("POSCAR", *PRODUCTION_INPUT_SHA256)}
            if any(hashes[filename] != digest
                   for filename, digest in PRODUCTION_INPUT_SHA256.items()):
                raise ValueError(f"VASP electronic contract changed at {name}")
            (case / "sha256.inputs.json").write_text(
                json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
            )
            records.append({
                "name": name,
                "image_index": index, "q_perp_A": sign * 0.015,
                "arc_fraction_s": feasibility["arc_fraction_s"][index],
                "path_center_enthalpy_eV_per_cell": feasibility["path_enthalpy_eV_per_cell"][index],
                **geometry, "input_sha256": hashes,
            })
    result = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_path_adapted_transverse_enthalpy_tube_pilot",
        "pressure_GPa": 45.7, "formula_units_per_cell": 2,
        "n_path_images": 29, "anchors": list(ANCHORS),
        "q_perp_A": [-0.015, 0.0, 0.015],
        "coordinate_definition": (
            "s = VCNEB joint arc fraction; q_perp = transported 600-eV candidate "
            "stable joint direction projected normal to each local path tangent"
        ),
        "q0_path_values_reused": True,
        "input_contract_sha256": PRODUCTION_INPUT_SHA256,
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "feasibility": sha256(feasibility_path),
            "hessian_npz": sha256(hessian_npz),
            "source_static_OUTCAR": sha256(source_static / "OUTCAR"),
            "preparer": sha256(Path(__file__)),
        },
        "cases": records,
        "limitations": [
            "Frozen off-path statics, not conditional minimization at fixed coordinates.",
            "Nine anchor locations are a pilot; interpolation requires holdout/error checks.",
            "The transverse joint normal is transported, not a fixed phonon eigenmode.",
        ],
    }
    (output / "manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "feasibility", "hessian-root", "source-static", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.trajectory, args.feasibility, args.hessian_root,
                     args.source_static, args.output)
    print(json.dumps({"status": report["status"], "n_offpath_cases": len(report["cases"])}))


if __name__ == "__main__":
    main()
