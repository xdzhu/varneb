"""Calculator-free global mode-subspace preflight in the universal ASE driver."""

from __future__ import annotations

import json
import hashlib
import sys

import numpy as np
import pytest
from ase import Atoms
from ase.io import read, write

from examples import run_vcneb_ase
from vcneb import (
    ModePlane, ReferenceCellCoordinates,
    endpoint_completed_mode_subspace_for_vcneb, interpolate_vcneb,
    make_vcneb_subspace_artifact,
)


def _case(tmp_path):
    reference = Atoms(
        "He", scaled_positions=[[0.21, 0.34, 0.43]],
        cell=[[4.0, 0.0, 0.0], [0.4, 4.8, 0.0], [0.1, 0.2, 5.1]], pbc=True,
    )
    chart = ReferenceCellCoordinates(reference)
    axes = np.zeros((9, 2))
    axes[0, 0], axes[1, 1] = 1.0, 1.0
    plane = ModePlane(
        reference=np.zeros(9), basis=np.eye(9), metric_weights=np.ones(9),
        axis_weights=axes, axis_labels=("x", "y"),
        amplitude_unit="toy length", reference_id="toy/ase-driver",
    )
    initial = chart.to_atoms(plane.frozen_coordinates([-0.1, 0.0]))
    final_chart = plane.frozen_coordinates([0.1, 0.0])
    final_chart[-4] = 0.02  # endpoint-completion cell direction
    final = chart.to_atoms(final_chart)
    scale = abs(np.linalg.det(initial.cell.array)) ** (1.0 / 3.0)
    subspace = endpoint_completed_mode_subspace_for_vcneb(
        plane, chart, initial, final, cell_scale_A=scale,
    )
    artifact = make_vcneb_subspace_artifact(
        subspace, initial, final, source_sha256={"toy_reference": "a" * 64},
    )
    paths = [tmp_path / name for name in ("initial.traj", "final.traj", "mode_subspace.json")]
    write(paths[0], initial)
    write(paths[1], final)
    paths[2].write_text(json.dumps(artifact), encoding="utf-8")
    return initial, final, paths


def _arguments(paths, workdir):
    return [
        "run_vcneb_ase.py", "--initial", str(paths[0]), "--final", str(paths[1]),
        "--subspace-artifact", str(paths[2]), "--workdir", str(workdir),
        "--n-images", "5", "--no-align-cells", "--cell-interpolation", "linear",
        "--mapping", "identity", "--validate-only",
    ]


def test_universal_ase_driver_accepts_audited_subspace_without_calculator(tmp_path, monkeypatch):
    initial, final, paths = _case(tmp_path)
    workdir = tmp_path / "validate"
    monkeypatch.setattr(sys, "argv", _arguments(paths, workdir))
    run_vcneb_ase.main()
    report = json.loads((workdir / "vcneb_preflight.json").read_text(encoding="utf-8"))
    assert report["mode_subspace"]["n_directions"] == 3
    assert report["mode_subspace"]["subspace_kind"] == "endpoint_completed"
    assert report["calculator"] is None
    assert report["calculator_validation"] == "not_instantiated_validate_only"
    assert (workdir / "initial-vcneb-unprojected.traj").is_file()
    images = read(workdir / "initial-vcneb.traj", index=":")
    assert len(images) == 5
    assert np.allclose(images[0].positions, initial.positions)
    assert np.allclose(images[-1].positions, final.positions)


def test_universal_ase_driver_rejects_off_subspace_resume(tmp_path, monkeypatch):
    initial, final, paths = _case(tmp_path)
    resumed = interpolate_vcneb(initial, final, 5, align_cells=False)
    resumed[2].positions[0, 2] += 0.03
    resume_path = tmp_path / "resume.traj"
    write(resume_path, resumed)
    monkeypatch.setattr(
        sys, "argv", _arguments(paths, tmp_path / "reject")
        + ["--resume-snapshot", str(resume_path)],
    )
    with pytest.raises(ValueError, match="outside the declared global mode subspace"):
        run_vcneb_ase.main()


def test_universal_ase_driver_projects_explicit_unconstrained_initial_chain(tmp_path, monkeypatch):
    initial, final, paths = _case(tmp_path)
    raw = interpolate_vcneb(initial, final, 5, align_cells=False)
    raw[2].positions[0, 2] += 0.03
    raw_path = tmp_path / "raw_chain.traj"
    write(raw_path, raw)
    workdir = tmp_path / "projected"
    monkeypatch.setattr(
        sys, "argv", _arguments(paths, workdir) + ["--initial-chain", str(raw_path)],
    )
    run_vcneb_ase.main()
    report = json.loads((workdir / "vcneb_preflight.json").read_text(encoding="utf-8"))
    assert report["mode_subspace"]["initial_projection_displacements_vcneb_A"][2] > 0.02
    projected = read(workdir / "initial-vcneb.traj", index=":")
    original = read(workdir / "initial-vcneb-unprojected.traj", index=":")
    assert np.linalg.norm(projected[2].positions - original[2].positions) > 0.02
    assert np.allclose(projected[0].positions, original[0].positions)
    assert np.allclose(projected[-1].positions, original[-1].positions)


def test_production_subspace_requires_a_pinned_artifact_hash(tmp_path, monkeypatch):
    _, _, paths = _case(tmp_path)
    args = _arguments(paths, tmp_path / "production")
    args.remove("--validate-only")
    args += ["--calculator", "ase.calculators.emt:EMT"]
    monkeypatch.setattr(sys, "argv", args)
    with pytest.raises(ValueError, match="production subspace runs require"):
        run_vcneb_ase.main()


def test_subspace_preflight_rejects_a_changed_artifact_hash(tmp_path, monkeypatch):
    _, _, paths = _case(tmp_path)
    monkeypatch.setattr(
        sys, "argv", _arguments(paths, tmp_path / "bad_pin")
        + ["--subspace-artifact-sha256", "0" * 64],
    )
    with pytest.raises(ValueError, match="does not match the pinned digest"):
        run_vcneb_ase.main()
    monkeypatch.setattr(
        sys, "argv", _arguments(paths, tmp_path / "good_pin")
        + ["--subspace-artifact-sha256", hashlib.sha256(paths[2].read_bytes()).hexdigest()],
    )
    run_vcneb_ase.main()
    report = json.loads((tmp_path / "good_pin" / "vcneb_preflight.json").read_text(encoding="utf-8"))
    assert report["mode_subspace"]["artifact_sha256_pinned"] is True
