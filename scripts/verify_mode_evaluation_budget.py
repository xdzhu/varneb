"""Bounded analytic budget/recovery check; no DFT or scheduler calls."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

import vcneb.mode_evaluator as implementation
from vcneb import CalculatorModeEvaluator, ReferenceCellCoordinates


def verify(workdir: Path) -> dict:
    if workdir.exists():
        raise FileExistsError("verification requires a fresh workdir")
    workdir.mkdir(parents=True)
    reference = Atoms("HeLi", scaled_positions=[[.1, .2, .3], [.6, .5, .4]],
                      cell=np.diag([4., 5., 6.]), pbc=True)
    chart = ReferenceCellCoordinates(reference)
    calls = []

    def factory(index, atoms, directory):
        calls.append(index)
        assert (directory.parent/f"attempt-{index:06d}.json").is_file()
        if index == 1:
            raise RuntimeError("intentional analytic factory failure")
        positions, cell = atoms.positions, atoms.cell.array
        energy = .5*(np.sum(positions**2)+.04*np.sum(cell**2))
        stress = (positions.T@positions+.04*cell.T@cell)/atoms.get_volume()
        return SinglePointCalculator(atoms, energy=energy, forces=-positions, stress=stress)

    options = dict(workdir=workdir/"capped", calculator_id="analytic:test:v1", max_evaluation_attempts=2)
    evaluate = CalculatorModeEvaluator(chart, factory, **options)
    zero = np.zeros(chart.coordinate_count)
    expected = evaluate(zero)
    for _ in range(3):
        cached = evaluate(zero)
        assert cached[0] == expected[0]
        np.testing.assert_array_equal(cached[1], expected[1])
    shifted = zero.copy()
    shifted[0] = .01
    try:
        evaluate(shifted)
    except RuntimeError as error:
        assert "intentional" in str(error)
    else:
        raise AssertionError("injected failure was not retained")
    resumed = CalculatorModeEvaluator(chart, factory, **options, resume=True)
    actual = resumed(zero)
    assert actual[0] == expected[0]
    np.testing.assert_array_equal(actual[1], expected[1])
    try:
        resumed(shifted)
    except RuntimeError as error:
        assert "attempt budget exhausted" in str(error)
    else:
        raise AssertionError("spent budget admitted a fresh attempt")
    changed_policy_refusals = 0
    for limit in (None, 3):
        try:
            CalculatorModeEvaluator(chart, factory, **dict(options, max_evaluation_attempts=limit), resume=True)
        except ValueError as error:
            assert "budget changed" in str(error)
            changed_policy_refusals += 1
        else:
            raise AssertionError("resume changed the registered budget")
    legacy = CalculatorModeEvaluator(chart, factory, workdir=workdir/"legacy", calculator_id="analytic:test:v1")
    assert legacy.contract == resumed.contract
    assert calls == [0, 1]
    assert evaluate.n_new_evaluations == 1 and resumed.n_new_evaluations == 0
    assert evaluate.n_evaluation_attempts == resumed.n_evaluation_attempts == 2
    source = Path(implementation.__file__).resolve()
    return dict(scope="analytic implementation replay; not material prediction",
                factory_invocations=calls, successful_new_current_process=evaluate.n_new_evaluations,
                successful_new_after_resume=resumed.n_new_evaluations,
                cumulative_reserved_attempts=resumed.n_evaluation_attempts,
                marker_files=len(list(resumed.workdir.glob("attempt-*.json"))),
                completed_result_files=len(list(resumed.workdir.glob("eval-*/result.json"))),
                failed_result_files=len(list(resumed.workdir.glob("eval-*/failure.json"))),
                exact_cache_recovery=True, changed_or_removed_budget_refusals=changed_policy_refusals,
                physical_contract_unchanged_from_legacy=True,
                physical_contract_sha256=resumed.contract["contract_sha256"],
                imported_implementation=str(source),
                imported_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                new_material_DFT_calls=0, scheduler_mutations=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    report = verify(args.workdir)
    with (args.workdir/"verification.json").open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps(report))


if __name__ == "__main__":
    main()
