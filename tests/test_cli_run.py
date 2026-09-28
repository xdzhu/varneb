"""Installed-package style config-to-execution regression without DFT."""

from __future__ import annotations

import json
from pathlib import Path
import sys

from ase import Atoms
from ase.io import write

from vcneb.cli import main
from vcneb import material_runner


def _toy_config(tmp_path: Path) -> Path:
    initial = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.75]], cell=[6, 6, 6], pbc=True)
    final = Atoms("H2", positions=[[0, 0, 0], [0, 0, 1.25]], cell=[6, 6, 6], pbc=True)
    write(tmp_path / "initial.vasp", initial, format="vasp", direct=True, vasp5=True)
    write(tmp_path / "final.vasp", final, format="vasp", direct=True, vasp5=True)
    path = tmp_path / "varneb.json"
    path.write_text(json.dumps({
        "schema_version": 1,
        "backend": "ase",
        "initial": "initial.vasp",
        "final": "final.vasp",
        "workdir": "run",
        "n_images": 3,
        "cell_mode": "fixed",
        "mapping": "identity",
        "steps": 1,
        "calculator": {
            "kind": "ase_class",
            "symbol": "ase.calculators.emt:EMT",
            "parameters": {},
        },
    }), encoding="utf-8")
    return path


def test_run_requires_explicit_execution(monkeypatch, tmp_path, capsys):
    config = _toy_config(tmp_path)
    monkeypatch.setattr(sys, "argv", ["varneb", "run", str(config)])
    assert main() == 2
    assert "requires --execute" in capsys.readouterr().err
    assert not (tmp_path / "run").exists()


def test_run_executes_same_config_with_ase_calculator(monkeypatch, tmp_path, capsys):
    config = _toy_config(tmp_path)
    monkeypatch.setattr(sys, "argv", ["varneb", "prepare", str(config)])
    assert main() == 0
    prepared = tmp_path / "run" / "initial-vcneb.traj"
    prepared_bytes = prepared.read_bytes()

    monkeypatch.setattr(sys, "argv", ["varneb", "run", str(config), "--execute"])
    assert main() == 0
    capsys.readouterr()
    summary = json.loads((tmp_path / "run" / "vcneb_summary.json").read_text(encoding="utf-8"))
    assert summary["backend_label"] == "ase"
    assert summary["requires_stress"] is False
    assert summary["n_interior_images"] == 1
    assert summary["endpoint_evaluation_policy"] == "fixed_cached_once"
    assert summary["config_source_sha256"] is not None
    assert prepared.read_bytes() == prepared_bytes

    # An accidental second run must not overwrite the successful case.
    monkeypatch.setattr(sys, "argv", ["varneb", "run", str(config), "--execute"])
    assert main() == 2
    assert "previous execution" in capsys.readouterr().err


def test_run_rejects_implicit_dft_launcher(monkeypatch, tmp_path, capsys):
    config = _toy_config(tmp_path)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["backend"] = "vasp"
    config.write_text(json.dumps(payload), encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["varneb", "run", str(config), "--execute"])
    assert main() == 2
    assert "explicit calculator.command" in capsys.readouterr().err
    assert not (tmp_path / "run").exists()


def test_config_dispatch_keeps_factory_and_optimizer_independent(monkeypatch, tmp_path):
    config = _toy_config(tmp_path)
    payload = json.loads(config.read_text(encoding="utf-8"))
    payload["backend"] = "cp2k"
    payload["optimizer"] = "SplitFIRE"
    payload["image_workers"] = 2
    payload["calculator"] = {
        "kind": "factory",
        "symbol": "vcneb.backends:make_ase_cp2k_factory",
        "parameters": {"cutoff_ry": 400},
        "command": "srun -n 32 cp2k_shell.psmp",
        "factory_kwargs": {},
    }
    config.write_text(json.dumps(payload), encoding="utf-8")
    seen: list[str] = []
    monkeypatch.setattr(material_runner, "main", lambda argv: seen.extend(argv))
    monkeypatch.setattr(sys, "argv", ["varneb", "run", str(config), "--execute"])
    assert main() == 0
    assert seen[seen.index("--factory") + 1] == "vcneb.backends:make_ase_cp2k_factory"
    assert seen[seen.index("--optimizer") + 1] == "SplitFIRE"
    assert seen[seen.index("--image-workers") + 1] == "2"
    assert seen[seen.index("--backend-label") + 1] == "cp2k"
    assert json.loads(seen[seen.index("--parameters-json") + 1]) == {"cutoff_ry": 400}
