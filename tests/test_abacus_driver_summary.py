"""ABACUS driver results must distinguish a finished run from a converged band."""

from __future__ import annotations

import json
import runpy
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "examples" / "run_vcneb_abacus.py"


@pytest.mark.parametrize("force,status", [(0.08, "completed"), (0.12, "step_limit_reached")])
def test_abacus_neb_summary_truthfully_records_convergence(tmp_path, monkeypatch, force, status) -> None:
    namespace = runpy.run_path(str(DRIVER))
    initial = Atoms("Ba", cell=[4, 4, 4], pbc=True)
    final = Atoms("Ba", scaled_positions=[[0.1, 0.0, 0.0]], cell=[4.1, 4, 4], pbc=True)
    initial_path, final_path = tmp_path / "initial.vasp", tmp_path / "final.vasp"
    write(initial_path, initial, format="vasp")
    write(final_path, final, format="vasp")

    def fake_factory(*, parameters, command):
        def make(index, atoms, directory):
            calculator = SinglePointCalculator(
                atoms, energy=-12.5, forces=np.zeros((1, 3)), stress=np.zeros(6)
            )
            calculator.directory = str(directory)
            return calculator
        return make

    def fake_run(images, **kwargs):
        return SimpleNamespace(
            images=images, enthalpies=np.zeros(len(images)), climb=kwargs["climb"],
            plot_band=lambda *args: None, barrier=lambda: (0.0, 0.0),
            get_forces=lambda: np.zeros((1, 3)), gradient_norm=lambda *args: force,
            path_diagnostics=lambda: {}, saddle_diagnostics=lambda: {},
        ), None

    monkeypatch.setitem(namespace["main"].__globals__, "make_ase_abacus_factory", fake_factory)
    monkeypatch.setitem(namespace["main"].__globals__, "run_vcneb", fake_run)
    monkeypatch.setattr(sys, "argv", [
        str(DRIVER), "--initial", str(initial_path), "--final", str(final_path),
        "--workdir", str(tmp_path / "run"), "--n-images", "3", "--steps", "0",
    ])
    namespace["main"]()
    summary = json.loads((tmp_path / "run" / "vcneb_summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == status
    assert summary["converged"] == (status == "completed")
    assert summary["climbing_image_requested"] is False
    assert summary["climbing_image_active_final"] is False
