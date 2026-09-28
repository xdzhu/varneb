"""Tests for the calculator-free configuration workflow."""

from __future__ import annotations

import json

import pytest

from ase import Atoms
from ase.io import write

from vcneb.config import RunConfig, prepare_run


def test_prepare_run_uses_safe_geometry_defaults(tmp_path) -> None:
    initial = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]], cell=[5, 5, 5], pbc=True)
    final = initial.copy()
    final.set_cell([5.2, 5.2, 5.2], scale_atoms=True)
    write(tmp_path / "initial.vasp", initial, format="vasp", direct=True, vasp5=True)
    write(tmp_path / "final.vasp", final, format="vasp", direct=True, vasp5=True)
    config_path = tmp_path / "varneb.json"
    config_path.write_text(json.dumps({
        "schema_version": 1,
        "backend": "abinit",
        "initial": "initial.vasp",
        "final": "final.vasp",
        "workdir": "run",
        "n_images": 5,
    }), encoding="utf-8")

    config, report_path = prepare_run(config_path)

    assert isinstance(config, RunConfig)
    assert config.mapping == "auto"
    assert config.cell_interpolation == "log_strain"
    assert (config.workdir / "initial-vcneb.traj").is_file()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["calculator_attached"] is False
    assert report["initial_path_geometry"]["valid"] is True
    _, repeated_report = prepare_run(config_path)
    assert repeated_report == report_path


def test_prepare_refuses_to_replace_a_changed_chain(tmp_path) -> None:
    initial = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]], cell=[5, 5, 5], pbc=True)
    final = initial.copy()
    final.positions[1, 2] = 1.0
    write(tmp_path / "initial.vasp", initial, format="vasp", direct=True, vasp5=True)
    write(tmp_path / "final.vasp", final, format="vasp", direct=True, vasp5=True)
    config_path = tmp_path / "varneb.json"
    config_path.write_text(json.dumps({
        "schema_version": 1, "backend": "ase", "initial": "initial.vasp",
        "final": "final.vasp", "workdir": "run",
    }), encoding="utf-8")
    prepare_run(config_path)
    trajectory = tmp_path / "run" / "initial-vcneb.traj"
    saved = trajectory.read_bytes()
    final.positions[1, 2] = 1.2
    write(tmp_path / "final.vasp", final, format="vasp", direct=True, vasp5=True)
    with pytest.raises(FileExistsError, match="prepared path"):
        prepare_run(config_path)
    assert trajectory.read_bytes() == saved


@pytest.mark.parametrize("key,value", [
    ("fmax_ev_per_angstrom", float("nan")),
    ("k", 0.0),
    ("pressure_gpa", float("inf")),
    ("image_workers", -1),
])
def test_run_config_rejects_nonphysical_or_unsafe_values(tmp_path, key, value) -> None:
    kwargs = {
        "backend": "ase", "initial": tmp_path / "is.vasp",
        "final": tmp_path / "fs.vasp", "workdir": tmp_path / "run",
        key: value,
    }
    with pytest.raises(ValueError):
        RunConfig(**kwargs)


def test_identity_mapping_rejects_reordered_endpoint(tmp_path) -> None:
    initial = Atoms("NaCl", positions=[[0, 0, 0], [1, 1, 1]], cell=[4, 4, 4], pbc=True)
    final = Atoms("ClNa", positions=[[0, 0, 0], [1, 1, 1]], cell=[4, 4, 4], pbc=True)
    write(tmp_path / "initial.vasp", initial, format="vasp", direct=True, vasp5=True)
    write(tmp_path / "final.vasp", final, format="vasp", direct=True, vasp5=True)
    config_path = tmp_path / "varneb.json"
    config_path.write_text(json.dumps({
        "schema_version": 1,
        "backend": "vasp",
        "initial": "initial.vasp",
        "final": "final.vasp",
        "workdir": "run",
        "mapping": "identity",
    }), encoding="utf-8")
    try:
        prepare_run(config_path)
    except ValueError as exc:
        assert "identity mapping" in str(exc)
    else:  # pragma: no cover - assertion guard
        raise AssertionError("identity mapping unexpectedly accepted reordered endpoints")


@pytest.mark.parametrize("field,value,message", [
    ("pressure_GPa", 45.7, "unknown config field"),
    ("n_images", 7.5, "n_images.*integer"),
    ("n_images", True, "n_images.*integer"),
    ("pressure_gpa", "45.7", "pressure_gpa.*number"),
    ("pressure_gpa", 10**400, "pressure_gpa.*numeric range"),
    ("mic", "false", "mic.*boolean"),
    ("align_translation", 0, "align_translation.*boolean"),
    ("climb", "false", "climb.*boolean"),
    ("steps", 10.5, "steps.*integer"),
])
def test_config_file_rejects_silent_coercions_and_typo(
    tmp_path, field, value, message
) -> None:
    payload = {
        "schema_version": 1,
        "backend": "ase",
        "initial": "initial.vasp",
        "final": "final.vasp",
        "workdir": "run",
    }
    payload[field] = value
    path = tmp_path / "varneb.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        RunConfig.from_file(path)


@pytest.mark.parametrize("payload,message", [
    (["not", "an", "object"], "JSON object"),
    ({"schema_version": True}, "schema_version must be 1"),
    ({"schema_version": 1, "backend": "ase", "initial": "is.vasp",
      "final": "fs.vasp", "workdir": "run",
      "calculator": {"commnad": "vasp_std"}}, "unknown calculator field"),
])
def test_config_file_rejects_invalid_shape_or_calculator_typo(
    tmp_path, payload, message
) -> None:
    path = tmp_path / "varneb.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        RunConfig.from_file(path)
