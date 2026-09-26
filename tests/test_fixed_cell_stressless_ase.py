"""Fixed-cell NEB needs energy and forces, but never a fabricated stress."""

from __future__ import annotations

import json
import sys

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import read, write

from examples import run_vcneb_ase
from vcneb import VCNEB, interpolate_vcneb, run_vcneb
from vcneb.calculator import CalculatorCapabilityError
from vcneb.executor import ThreadedCalculatorExecutor


class StresslessSlide(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x = float(atoms.positions[0, 0])
        energy = (x * (2.0 - x)) ** 2
        derivative = 4.0 * x * (2.0 - x) * (1.0 - x)
        forces = np.zeros((len(atoms), 3))
        forces[0, 0] = -derivative
        self.results = {"energy": energy, "forces": forces}

    def get_stress(self, *args, **kwargs):
        raise AssertionError("a fixed-cell path must not ask this calculator for stress")


def _images():
    initial = Atoms("H", positions=[[0.0, 1.0, 1.0]], cell=[6.0, 6.0, 20.0], pbc=True)
    final = initial.copy()
    final.positions[0, 0] = 2.0
    images = interpolate_vcneb(initial, final, 5, align_cells=False)
    for image in images:
        image.calc = StresslessSlide()
    return images


def test_fixed_cell_serial_runs_without_stress_and_reports_unavailable():
    images = _images()
    chain, _ = run_vcneb(
        images, cell_mask=np.zeros((3, 3)), climb=False, steps=0,
        logfile=None, trajectory=None, snapshot_dir=None,
    )
    diagnostics = chain.path_diagnostics()
    assert len(diagnostics["images"]) == len(images)
    assert all(row["stress_available"] is False and row["stress_eV_per_A3"] is None
               for row in diagnostics["images"])
    assert all(np.array_equal(image.cell.array, images[0].cell.array) for image in images)


def test_fixed_cell_threaded_cache_and_snapshots_do_not_claim_stress(tmp_path):
    images = _images()
    executor = ThreadedCalculatorExecutor(
        2, cache_dir=tmp_path / "cache", cache_namespace="stressless-fixed",
        require_stress=False,
    )
    chain, _ = run_vcneb(
        images, cell_mask=np.zeros((3, 3)), image_executor=executor,
        climb=False, steps=0, logfile=None,
        trajectory=tmp_path / "chain.traj", snapshot_dir=tmp_path / "snapshots",
    )
    assert all(not value.stress_available for value in chain._last_evaluations)
    assert all("stress" not in image.calc.results for image in read(tmp_path / "chain.traj", index=":"))
    again = _images()
    # The controller evaluates fixed endpoints once outside the worker cache.
    cached = executor.evaluate(again[1:-1], indices=range(1, len(again) - 1))
    assert all(not value.stress_available for value in cached)
    assert all(attempt == 0 for attempt in executor.last_attempts.values())
    with pytest.raises(ValueError, match="cache metadata"):
        ThreadedCalculatorExecutor(
            1, cache_dir=tmp_path / "cache", cache_namespace="stressless-fixed",
            require_stress=True,
        )


def test_variable_cell_still_rejects_stressless_calculators_before_dft():
    images = _images()
    with pytest.raises(CalculatorCapabilityError, match="stress"):
        run_vcneb(images, steps=0, logfile=None, trajectory=None, snapshot_dir=None)
    assert all(not image.calc.results for image in images)


def test_fixed_cell_rejects_different_image_cells():
    images = _images()
    images[2].cell[0, 0] += 0.01
    with pytest.raises(ValueError, match="same cell"):
        VCNEB(images, cell_mask=np.zeros((3, 3)))


def test_ase_cli_static_fixed_cell_accepts_stressless_calculator(tmp_path, monkeypatch):
    initial, final = _images()[0], _images()[-1]
    initial_path = tmp_path / "initial.vasp"
    final_path = tmp_path / "final.vasp"
    write(initial_path, initial, format="vasp")
    write(final_path, final, format="vasp")
    monkeypatch.setattr(run_vcneb_ase, "_load_symbol", lambda _: StresslessSlide)
    monkeypatch.setattr(sys, "argv", [
        "run_vcneb_ase.py", "--initial", str(initial_path), "--final", str(final_path),
        "--workdir", str(tmp_path / "run"), "--calculator", "toy:StresslessSlide",
        "--n-images", "3", "--cell-mode", "fixed", "--static-only",
        "--mapping", "identity", "--cell-interpolation", "linear", "--no-align-cells",
    ])
    run_vcneb_ase.main()
    summary = json.loads((tmp_path / "run" / "ase_static_summary.json").read_text())
    assert summary["cell_mode"] == "fixed"
    assert summary["requires_stress"] is False
    assert len(summary["static_endpoints"]) == 2
    assert all(row["stress_available"] is False and row["stress_eV_per_A3"] is None
               for row in summary["static_endpoints"])


def test_ase_cli_threaded_fixed_cell_completes_without_stress(tmp_path, monkeypatch):
    images = _images()
    initial_path = tmp_path / "initial.vasp"
    final_path = tmp_path / "final.vasp"
    write(initial_path, images[0], format="vasp")
    write(final_path, images[-1], format="vasp")
    monkeypatch.setattr(run_vcneb_ase, "_load_symbol", lambda _: StresslessSlide)
    monkeypatch.setattr(sys, "argv", [
        "run_vcneb_ase.py", "--initial", str(initial_path), "--final", str(final_path),
        "--workdir", str(tmp_path / "run"), "--calculator", "toy:StresslessSlide",
        "--n-images", "5", "--cell-mode", "fixed", "--image-workers", "2",
        "--steps", "0", "--no-climb", "--mapping", "identity",
        "--cell-interpolation", "linear", "--no-align-cells",
    ])
    run_vcneb_ase.main()
    summary = json.loads((tmp_path / "run" / "vcneb_summary.json").read_text())
    assert summary["cell_mode"] == "fixed"
    assert summary["requires_stress"] is False
    assert all(row["stress_available"] is False and row["stress_eV_per_A3"] is None
               for row in summary["path_diagnostics"]["images"])
    assert all("stress" not in image.calc.results
               for image in read(tmp_path / "run" / "vcneb.traj", index=":"))


def test_ase_cli_fixed_cell_rejects_pressure_and_mismatched_cells():
    images = _images()
    with pytest.raises(ValueError, match="pressure_gpa"):
        run_vcneb_ase._cell_mask_for_mode("fixed", images, 1.0)
    images[-1].cell[0, 0] += 0.01
    with pytest.raises(ValueError, match="identical cells"):
        run_vcneb_ase._cell_mask_for_mode("fixed", images, 0.0)
