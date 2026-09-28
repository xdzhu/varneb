"""No-DFT selection checks for the frozen BTO and local GaN refinements."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import csv

import pytest
import numpy as np

from examples.prepare_bto_transverse_soft_missing22 import stage
from scripts.prepare_gan_600eV_ts_2d_refinement import coordinate_key, missing_coordinates
from scripts.prepare_gan_600eV_ts_2d_dense9 import missing_9x9
from scripts.plot_gan_600eV_local_joint_cut import measured_interpolant, source_rows


ROOT = Path(__file__).resolve().parents[1]


def _sha256(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_bto_59_to_81_reuses_only_identical_structures(tmp_path) -> None:
    dense_dir = tmp_path / "dense"
    dense_dir.mkdir()
    points = []
    for j in range(9):
        for i in range(9):
            name = f"grid-q1-{i:02d}-q2-{j:02d}"
            relative = f"{name}/POSCAR"
            path = dense_dir / relative
            path.parent.mkdir()
            path.write_text(f"geometry {i} {j}\n", encoding="ascii")
            points.append({
                "name": name, "structure": relative,
                "structure_sha256": _sha256(path),
                "grid_index_q1_q2": [i, j], "q1": float(i), "q2": float(j),
            })
    dense = {
        "kind": "inert_frozen_cubic_cell_transverse_soft_mode_points_not_T_to_C_barrier",
        "n_points": 81, "n_atoms": 5, "points": points,
        "axis_labels": ["Q1", "Q2"], "axis_definition": "soft/soft",
        "axis_mode_weights": [[0.0, 0.0]] * 15,
        "boundary_condition": "frozen", "amplitude_unit": "sqrt(amu)*angstrom",
        "reference_id": "reference", "input_sha256": {"reference": "hash"},
    }
    dense_path = dense_dir / "manifest.json"
    dense_path.write_text(json.dumps(dense), encoding="utf-8")
    old = {
        "kind": "BTO_transverse_soft_frozen_C_cell_59_real_DFT_points_not_conditional_PES_or_MEP",
        "status": "eighteen_independent_adaptive_edge_holdouts_scored",
        "n_real_DFT_points": 59,
        "source_sha256": {"dense_manifest": _sha256(dense_path)},
        "axis_labels": dense["axis_labels"],
        "axis_mode_weights": dense["axis_mode_weights"],
        "boundary_condition": dense["boundary_condition"],
        "samples": [{
            "dense_i_q1": point["grid_index_q1_q2"][0],
            "dense_j_q2": point["grid_index_q1_q2"][1],
            "q_parallel_sqrt_amu_A": point["q1"],
            "q_transverse_sqrt_amu_A": point["q2"],
            "source_structure_sha256": point["structure_sha256"],
        } for point in points[:59]],
    }
    old_path = tmp_path / "old.json"
    old_path.write_text(json.dumps(old), encoding="utf-8")
    output = tmp_path / "new"
    result = stage(dense_dir, old_path, output)
    assert result["n_points"] == 22
    assert {tuple(point["grid_index_q1_q2"]) for point in result["points"]} == {
        tuple(point["grid_index_q1_q2"]) for point in points[59:]
    }
    assert all(_sha256(output / point["structure"]) == point["structure_sha256"]
               for point in result["points"])
    with pytest.raises(FileExistsError):
        stage(dense_dir, old_path, output)


def test_gan_13_to_25_joint_grid_is_disjoint_and_complete() -> None:
    prior = {coordinate_key(q_u, q_v)
             for q_u in (-0.02, 0.0, 0.02)
             for q_v in (-0.0125, 0.0, 0.0125)}
    prior.update({(-0.01, 0.0), (0.01, 0.0), (0.0, -0.00625), (0.0, 0.00625)})
    u, v, missing = missing_coordinates(prior)
    fresh = {coordinate_key(q_u, q_v) for _, _, q_u, q_v in missing}
    assert len(fresh) == 12 and not fresh & prior
    assert fresh | prior == {coordinate_key(q_u, q_v) for q_u in u for q_v in v}
    with pytest.raises(ValueError):
        missing_coordinates(prior | {(0.03, 0.0)})


def test_gan_dft_contour_uses_complete_measured_grid() -> None:
    u = np.linspace(-0.02, 0.02, 5)
    v = np.linspace(-0.0125, 0.0125, 5)
    u_mesh, v_mesh = np.meshgrid(u, v)
    rows = [{"q_u_A": float(x), "q_v_A": float(y),
             "measured_delta_H_meV_per_GaN": float(-x*x + y*y)}
            for x, y in zip(u_mesh.ravel(), v_mesh.ravel())]
    interpolated = measured_interpolant(rows, u_mesh, v_mesh)
    assert interpolated == pytest.approx(-u_mesh*u_mesh + v_mesh*v_mesh, abs=1e-12)
    with pytest.raises(ValueError):
        measured_interpolant(rows[:-1], u_mesh, v_mesh)
    u9, v9 = np.meshgrid(np.linspace(-0.02, 0.02, 9),
                         np.linspace(-0.0125, 0.0125, 9))
    dense_rows = [{"q_u_A": float(x), "q_v_A": float(y),
                   "measured_delta_H_meV_per_GaN": float(-x*x + y*y)}
                  for x, y in zip(u9.ravel(), v9.ravel())]
    assert measured_interpolant(dense_rows, u9, v9) == pytest.approx(
        -u9*u9 + v9*v9, abs=1e-12)


def test_gan_5x5_to_9x9_selects_only_unmeasured_nodes() -> None:
    prior = {coordinate_key(float(x), float(y))
             for x in np.linspace(-0.02, 0.02, 5)
             for y in np.linspace(-0.0125, 0.0125, 5)}
    u, v, missing = missing_9x9(prior)
    fresh = {coordinate_key(x, y) for _, _, x, y in missing}
    assert len(fresh) == 56 and not fresh & prior
    assert fresh | prior == {coordinate_key(x, y) for x in u for y in v}
    with pytest.raises(ValueError):
        missing_9x9(prior | {(0.03, 0.0)})


def test_bto_81_point_figure_tracks_raw_audit() -> None:
    audit_path = (ROOT / "benchmarks/numerical_integrity/"
                  "bto_transverse_soft_frozen81_20260928.json")
    qa_path = (ROOT / "paper/VARNEB_CPC/figures/"
               "bto_frozen_soft_mode_81_20260928_qa.json")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    qa = json.loads(qa_path.read_text(encoding="utf-8"))
    assert audit["status"] == "22_new_nodes_raw_energy_force_stress_audited"
    assert audit["n_total_DFT_points"] == 81
    assert qa["source_sha256"]["analysis81"] == _sha256(audit_path)
    assert qa["interpolation"]["n_measured_DFT_samples"] == 81
    assert qa["interpolation"]["independent_22_missing_node_max_abs_error_meV_per_BTO"] < 0.5


def test_gan_81_row_join_reuses_25_audited_points() -> None:
    pilot_path = ROOT / "benchmarks/numerical_integrity/gan_600eV_ts_2d_pilot_20260928.json"
    refined5_path = ROOT / "benchmarks/numerical_integrity/gan_600eV_ts_2d_5x5_refinement_20260928.json"
    pilot = json.loads(pilot_path.read_text(encoding="utf-8"))
    refined5 = json.loads(refined5_path.read_text(encoding="utf-8"))
    u, v, missing = missing_9x9({coordinate_key(float(x), float(y))
                                for x in np.linspace(-0.02, 0.02, 5)
                                for y in np.linspace(-0.0125, 0.0125, 5)})
    fake_dense = {
        "status": "GaN_600eV_local_joint_9x9_refinement_raw_audited",
        "n_total_DFT_points": 81,
        "source_sha256": {"pilot_audit": _sha256(pilot_path),
                          "refinement5_audit": _sha256(refined5_path)},
        "cases": [{"case": f"synthetic_{i}_{j}", "q_u_A": x, "q_v_A": y,
                   "delta_enthalpy_meV_per_GaN": float(-x*x + y*y),
                   "outcar_sha256": "synthetic-test-only"}
                  for i, j, x, y in missing],
    }
    rows = source_rows(pilot, _sha256(pilot_path), refined5,
                       _sha256(refined5_path), fake_dense, "synthetic-test-only")
    assert len(rows) == 81
    assert len({(row["q_u_A"], row["q_v_A"]) for row in rows}) == 81
    assert {float(row["q_u_A"]) for row in rows} == set(u)
    assert {float(row["q_v_A"]) for row in rows} == set(v)


def test_gan_81_point_figure_tracks_raw_audit() -> None:
    audit_path = (ROOT / "benchmarks/numerical_integrity/"
                  "gan_600eV_ts_2d_9x9_refinement_20260928.json")
    prefix = (ROOT / "paper/VARNEB_CPC/figures/"
              "gan_600eV_local_joint_dft81_v2_20260928")
    qa = json.loads(prefix.with_name(prefix.name + "_qa.json").read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    with prefix.with_name(prefix.name + "_source_data.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 81 and audit["n_total_DFT_points"] == 81
    assert audit["status"] == "GaN_600eV_local_joint_9x9_refinement_raw_audited"
    assert qa["source_sha256"]["dense9_raw_audit"] == _sha256(audit_path)
    assert qa["TS_certified"] is False and qa["whole_path_2D_surface_certified"] is False
    assert qa["prospective_5x5_max_abs_error_meV_per_GaN"] < 0.02
