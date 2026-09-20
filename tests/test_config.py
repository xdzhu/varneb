"""Tests for the calculator-free configuration workflow."""

from __future__ import annotations

import json

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
