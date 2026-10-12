"""Real live snapshots: residual diagnosis is not material convergence."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
CASE = REPO / "benchmarks/hfo2_channels/20261008/clamped_flip_residual_E139_20261012"
SNAPSHOTS = [
    ("28709790", "28737465", 37,
     "a4593fcf2c0f1613cced5e884afcf58518b1285c38ffd3aa9336f043e5e4ce31",
     0.11760457679798154, 0.059113965810629034,
     2.8811660508892304e-7, 0.005878231381024376),
    ("28709791", "28764704", 15,
     "c26a3c8f783c8c758fdf5716bb3aea6dff2885f5f19c5fcde3629ca43f55282f",
     0.14367414108221876, 0.010533619598725376,
     3.606716254545368e-12, 0.0053451093461657495),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("snapshot", SNAPSHOTS, ids=["preserving37", "reversing15"])
def test_native_observation_is_complete_not_converged(snapshot):
    from examples.hfo2_fixed_input_factory import CONTRACT

    parent, job, step, pin, *_ = snapshot
    directory = CASE / "completed_HF" / parent / f"step_{step:04d}"
    assert sha(directory / "observation.json") == pin
    report = json.loads((directory / "observation.json").read_text())
    assert report["status"] == "complete_observation_not_final_result"
    assert report["source_job_id"] == job and report["snapshot_step"] == step
    assert report["new_DFT_calls"] == 0 and not report["ordinary_residual_pass"]
    assert report["pressure_GPa"] == 0 and report["climb"] is False
    assert report["fmax_target_eV_A"] == .10 and report["k_eV_A2"] == .2
    assert report["n_total_images"] == 9 and report["n_active_images"] == 7
    assert report["formula_units"] == 4
    assert report["mechanical_boundary"]["fixed_cell_rows"] == [0, 1]
    assert report["mechanical_boundary"]["allow_tilt"] is True
    assert len(report["raw_image_evaluations"]) == 9
    for row in report["raw_image_evaluations"]:
        assert {name: row["input_sha256"][name] for name in CONTRACT} == CONTRACT
        assert not row["new_SCF_for_export"]


@pytest.mark.parametrize("snapshot", SNAPSHOTS, ids=["preserving37", "reversing15"])
def test_retained_raw_scfs_and_ordered_geometry(snapshot):
    from ase.io import read
    from examples.hfo2_fixed_input_factory import read_fixed_hfo2_stru, same_ordered_geometry
    from scripts.audit_hfo2_static_replica import audited_results

    parent, _, step, pin, *_ = snapshot
    directory = CASE / "completed_HF" / parent / f"step_{step:04d}"
    assert sha(directory / "observation.json") == pin
    report = json.loads((directory / "observation.json").read_text())
    images = read(directory / "evaluated_chain.traj", index=":")
    assert sha(directory / "evaluated_chain.traj") == report["evaluated_chain_sha256"]
    assert len(images) == 9
    for image, row in zip(images, report["raw_image_evaluations"]):
        raw = directory / "raw" / f"image_{row['image_index']:04d}"
        assert sha(raw / "OUT.ABACUS/running_scf.log") == row["raw_log_sha256"]
        assert sha(raw / "source_audit.json") == row["audit_sha256"]
        for name in ("INPUT", "KPT", "STRU"):
            assert sha(raw / name) == row["input_sha256"][name]
        assert same_ordered_geometry(image, read_fixed_hfo2_stru(raw / "STRU"))
        parsed = audited_results(raw)
        for key, expected in (("energy", row["energy_eV_cell"]),
                              ("forces", row["forces_eV_A"]),
                              ("stress", row["stress_ASE_voigt_eV_A3"])):
            np.testing.assert_allclose(parsed[key], expected, atol=1e-12, rtol=0)


@pytest.mark.parametrize("snapshot", SNAPSHOTS, ids=["preserving37", "reversing15"])
def test_native_residual_replay_excludes_spring_dominance(snapshot):
    from scripts.analyze_hfo2_clamped_residual import analyze

    parent, job, step, _, fmax, cell_max, spring_projection, chain_spring_max = snapshot
    directory = CASE / "completed_HF" / parent / f"step_{step:04d}"
    result = analyze(directory)
    assert result["source_job_id"] == job and result["snapshot_step"] == step
    assert result["new_DFT_calls"] == 0 and not result["ordinary_residual_pass"]
    assert result["fmax_eV_A"] == pytest.approx(fmax, abs=1e-12)
    assert result["limiting_image_index"] == 4 and result["limiting_vector_block"] == "atoms"
    winner = next(row for row in result["images"] if row["image_index"] == 4)
    assert winner["residual_open_cell_max_vector_eV_A"] == pytest.approx(cell_max, abs=1e-12)
    assert winner["spring_projection_on_worst_residual_eV_A"] == pytest.approx(
        spring_projection, abs=1e-12)
    assert winner["component_sum_max_abs_error_eV_A"] < 1e-10
    assert (winner["perpendicular_projection_on_worst_residual_eV_A"]
            + winner["spring_projection_on_worst_residual_eV_A"]) == pytest.approx(fmax, abs=1e-12)
    assert result["spring_only_fmax_same_geometry_eV_A"] == pytest.approx(
        chain_spring_max, abs=1e-12)
