"""The GaN candidate Hessian probes must stay in the original four-atom cell."""

import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.audit_gan_joint_curvature import tangent_in_joint_coordinates
from scripts.prepare_gan_joint_curvature import prepare, same_geometry
from vcneb.joint_curvature import JointCurvatureCoordinates


ROOT = Path(__file__).resolve().parents[1]


def test_prepare_complete_joint_probe_set_without_supercell(tmp_path) -> None:
    work = tmp_path / "joint"
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "status": "completed", "converged": True, "n_images": 29,
        "endpoint_evaluation_policy": "fixed_cached_once",
        "path_diagnostics": {"cell_scale_A": 3.3982714330050063},
    }), encoding="utf-8")
    manifest = prepare(
        ROOT / "outputs/gan_b4_b1_gamma_1x1x1_20260926/gan_vasp_tetragonal_final_29_images.traj",
        summary,
        work, 0.01,
    )
    assert manifest["status"] == "geometry_staged_no_DFT"
    assert manifest["supercell_matrix"] == np.eye(3, dtype=int).tolist()
    assert len(manifest["cases"]) == 36
    assert {case["n_atoms"] for case in manifest["cases"]} == {4}
    assert all(case["minimum_distance_A"] > 1.4 for case in manifest["cases"])
    center = read(work / "center_POSCAR", format="vasp")
    wrapped = center.copy()
    wrapped.wrap()
    assert same_geometry(center, wrapped)
    for axis in range(18):
        plus = read(work / f"cases/axis{axis:02d}_plus/POSCAR", format="vasp")
        minus = read(work / f"cases/axis{axis:02d}_minus/POSCAR", format="vasp")
        assert plus.get_chemical_symbols() == center.get_chemical_symbols()
        assert minus.get_chemical_symbols() == center.get_chemical_symbols()
        assert plus.get_volume() > 0 and minus.get_volume() > 0
        if axis < 12:
            np.testing.assert_allclose(plus.cell.array, center.cell.array, atol=1e-12)
            np.testing.assert_allclose(minus.cell.array, center.cell.array, atol=1e-12)
        else:
            np.testing.assert_allclose(
                plus.cell.array + minus.cell.array, 2 * center.cell.array, atol=1e-10,
            )


def test_archived_neighbor_tangent_is_well_defined() -> None:
    chain = read(
        ROOT / "outputs/gan_b4_b1_gamma_1x1x1_20260926/gan_vasp_tetragonal_final_29_images.traj",
        index=":",
    )
    coordinates = JointCurvatureCoordinates(chain[15], cell_scale_A=3.3982714330050063)
    tangent, rotation = tangent_in_joint_coordinates(coordinates, chain[14], chain[16])
    np.testing.assert_allclose(np.linalg.norm(tangent), 1.0, atol=1e-12)
    assert rotation < 1e-4
