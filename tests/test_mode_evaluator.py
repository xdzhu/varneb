"""Recoverable ASE energy/gradient evaluation without a DFT backend."""

from __future__ import annotations

import json
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from vcneb.mode_evaluator import CalculatorModeEvaluator
from vcneb.reference_cell import ReferenceCellCoordinates


def _chart() -> ReferenceCellCoordinates:
    reference = Atoms(
        "HeLi", scaled_positions=[[0.1, 0.2, 0.3], [0.6, 0.5, 0.4]],
        cell=np.diag([4.0, 5.0, 6.0]), pbc=True,
    )
    return ReferenceCellCoordinates(reference)


def _factory(counter: list[int]):
    def make(index: int, atoms: Atoms, directory):
        counter.append(index)
        positions, cell = atoms.positions, atoms.cell.array
        energy = 0.5 * (np.sum(positions**2) + 0.04 * np.sum(cell**2))
        forces = -positions
        stress = (positions.T @ positions + 0.04 * cell.T @ cell) / atoms.get_volume()
        return SinglePointCalculator(atoms, energy=energy, forces=forces, stress=stress)
    return make


def test_evaluator_persists_and_resumes_exact_completed_coordinates(tmp_path) -> None:
    chart = _chart()
    counter: list[int] = []
    evaluator = CalculatorModeEvaluator(
        chart, _factory(counter), workdir=tmp_path / "run",
        calculator_id="analytic:test:v1", pressure_eV_per_A3=0.002,
        minimum_distance_A=1.0,
    )
    zero = np.zeros(chart.coordinate_count)
    first_energy, first_gradient = evaluator(zero)
    repeat_energy, repeat_gradient = evaluator(zero)
    assert first_energy == repeat_energy
    assert np.allclose(first_gradient, repeat_gradient)
    assert counter == [0]
    changed = zero.copy()
    changed[0] = 0.01
    evaluator(changed)
    assert counter == [0, 1]
    assert evaluator.n_new_evaluations == 2
    assert len(list((tmp_path / "run").glob("eval-*/result.json"))) == 2

    resumed_counter: list[int] = []
    resumed = CalculatorModeEvaluator(
        chart, _factory(resumed_counter), workdir=tmp_path / "run",
        calculator_id="analytic:test:v1", pressure_eV_per_A3=0.002,
        minimum_distance_A=1.0, resume=True,
    )
    resumed_energy, resumed_gradient = resumed(zero)
    assert resumed_energy == first_energy
    assert np.allclose(resumed_gradient, first_gradient)
    assert resumed_counter == []
    with pytest.raises(ValueError, match="resume contract changed"):
        CalculatorModeEvaluator(
            chart, _factory([]), workdir=tmp_path / "run",
            calculator_id="different-functional", pressure_eV_per_A3=0.002,
            minimum_distance_A=1.0, resume=True,
        )


def test_unsafe_geometry_rejected_before_calculator_and_failed_eval_preserved(tmp_path) -> None:
    chart = _chart()
    zero = np.zeros(chart.coordinate_count)
    counter: list[int] = []
    too_strict = CalculatorModeEvaluator(
        chart, _factory(counter), workdir=tmp_path / "strict",
        calculator_id="analytic:test:v1", minimum_distance_A=10.0,
    )
    with pytest.raises(ValueError, match="calculator was not called"):
        too_strict(zero)
    assert counter == []
    assert not list((tmp_path / "strict").glob("eval-*"))

    def fail(index, atoms, directory):
        raise RuntimeError("injected calculator failure")

    failing = CalculatorModeEvaluator(
        chart, fail, workdir=tmp_path / "recover", calculator_id="analytic:test:v1",
    )
    with pytest.raises(RuntimeError, match="injected"):
        failing(zero)
    failed = list((tmp_path / "recover").glob("eval-*/failure.json"))
    assert len(failed) == 1
    assert json.loads(failed[0].read_text(encoding="utf-8"))["status"] == "failed"
    recovered_counter: list[int] = []
    resumed = CalculatorModeEvaluator(
        chart, _factory(recovered_counter), workdir=tmp_path / "recover",
        calculator_id="analytic:test:v1", resume=True,
    )
    resumed(zero)
    assert recovered_counter == [1]
    assert len(list((tmp_path / "recover").glob("eval-*/result.json"))) == 1


def test_result_validator_blocks_unreviewed_calculator_output(tmp_path) -> None:
    chart = _chart()
    zero = np.zeros(chart.coordinate_count)

    def reject(index, atoms, directory):
        raise ValueError("SCF convergence marker missing")

    validator = CalculatorModeEvaluator(
        chart, _factory([]), workdir=tmp_path / "validated",
        calculator_id="analytic:test:v1", result_validator=reject,
        validation_id="scf-log-audit:v1",
    )
    with pytest.raises(ValueError, match="SCF convergence"):
        validator(zero)
    assert len(list((tmp_path / "validated").glob("eval-*/failure.json"))) == 1
    assert not list((tmp_path / "validated").glob("eval-*/result.json"))

    def accept(index, atoms, directory):
        return {"scf_converged": True, "index": index}

    resumed = CalculatorModeEvaluator(
        chart, _factory([]), workdir=tmp_path / "validated",
        calculator_id="analytic:test:v1", result_validator=accept,
        validation_id="scf-log-audit:v1", resume=True,
    )
    resumed(zero)
    result_path = next((tmp_path / "validated").glob("eval-*/result.json"))
    assert json.loads(result_path.read_text(encoding="utf-8"))["validation"]["scf_converged"] is True
    with pytest.raises(ValueError, match="resume contract changed"):
        CalculatorModeEvaluator(
            chart, _factory([]), workdir=tmp_path / "validated",
            calculator_id="analytic:test:v1", result_validator=accept,
            validation_id="different-validator", resume=True,
        )
