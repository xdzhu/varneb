"""Synthetic arithmetic/coverage evidence, not HfO2 held-out predictions."""
from dataclasses import replace

import numpy as np
import pytest

from vcneb.prediction_evaluation import (
    AuditedBarrierReference, FrozenChannelForecast, absolute_error_interval,
    evaluate_frozen_panel, shared_reference_gain_interval,
)


CONTRACT = 'a'*64


def panel(*, near_baseline=False):
    references = [AuditedBarrierReference('switch', .2, .005, CONTRACT, 'b'*64),
                  AuditedBarrierReference('escape', .4, .005, CONTRACT, 'c'*64)]
    values = {'B5': [.2, .4], 'B0': [.201, .399] if near_baseline else [.23, .43],
              'B2_Cmma_gauge1': [.23, .43], 'B2_Cmma_gauge2': [.24, .44]}
    forecasts = [FrozenChannelForecast(model, reference.case_id, value, None, CONTRACT, str(i+1)*64)
                 for i, (model, predictions) in enumerate(values.items())
                 for reference, value in zip(references, predictions)]
    return forecasts, references


def score(forecasts, references):
    return evaluate_frozen_panel(forecasts, references, candidate_model='B5',
        baseline_models=['B0', 'B2_Cmma_gauge1', 'B2_Cmma_gauge2'])


@pytest.mark.parametrize('prediction,interval,expected', [
    (.2, (.19, .21), (0, .01)), (.25, (.19, .21), (.04, .06)),
    (.1, (.19, .21), (.09, .11)), (.2, (.2, .2), (0, 0))])
def test_absolute_error_interval(prediction, interval, expected):
    assert absolute_error_interval(prediction, interval) == pytest.approx(expected)


def test_shared_reference_is_not_two_independent_error_draws():
    # Same prediction, same unknown reference: gain is identically zero,
    # although each absolute error separately spans [0,.01].
    assert shared_reference_gain_interval(.2, [.2], [.19, .21]) == (0, 0)
    assert shared_reference_gain_interval(.2, [.22], [.19, .21]) == pytest.approx((0, .02), abs=1e-15)


def test_strong_baseline_envelope_needs_midpoint_breakpoints():
    # Candidate .2, baselines .1/.3, y in[.1,.3]: gain -0.1 at endpoints
    # and +0.1 at their midpoint, not at a baseline prediction.
    assert shared_reference_gain_interval(.2, [.1, .3], [.1, .3]) == pytest.approx((-.1, .1))


def test_48_random_intervals_independent_dense_grid_encloses_and_approaches_extrema():
    rng = np.random.default_rng(61463)
    for _ in range(48):
        lower, upper = sorted(rng.uniform(-.5, .5, 2))
        candidate = rng.uniform(-.6, .6)
        baselines = rng.uniform(-.6, .6, 4)
        lo, hi = shared_reference_gain_interval(candidate, baselines, (lower, upper))
        grid = np.linspace(lower, upper, 10001)
        gain = np.min(np.abs(baselines[:, None]-grid), axis=0)-np.abs(candidate-grid)
        gap = 2*(upper-lower)/10000+1e-14  # Lipschitz2, not a statistical error
        assert lo-1e-14 <= gain.min() <= lo+gap
        assert hi-gap <= gain.max() <= hi+1e-14


def test_resolved_gain_requires_every_registered_baseline_and_case():
    forecasts, references = panel()
    original = [(f.model_id, f.case_id, f.prediction_eV_fu) for f in forecasts]
    report = score(forecasts, references)
    assert report['full_panel_strict_advantage']
    assert all(r['complete_registered_panel_strict_gain'] for r in report['comparisons'])
    assert [(f.model_id, f.case_id, f.prediction_eV_fu) for f in forecasts] == original
    assert report['new_DFT_calls'] == 0 and report['energy_unit'] == 'eV/formula_unit'
    assert report['models'][0]['raw_maximum_absolute_error_eV_fu'] == 0


def test_one_strong_reference_defeats_a_weak_reference_advantage_claim():
    forecasts, references = panel(near_baseline=True)
    report = score(forecasts, references)
    assert not report['full_panel_strict_advantage']
    for row in report['comparisons']:
        assert row['pairwise_error_gain_interval_eV_fu']['B2_Cmma_gauge1'][0] > 0
        assert row['best_available_baseline_oracle_gain_interval_eV_fu'][0] < 0


def test_abstention_preserves_signed_value_and_reports_coverage_not_zero_error():
    forecasts, references = panel()
    forecasts[1] = replace(forecasts[1], prediction_eV_fu=-.01, abstention_reason='branch domain failed')
    report = score(forecasts, references)
    candidate = report['models'][0]
    assert candidate['coverage'] == .5
    assert candidate['cases'][1]['prediction_eV_fu'] == -.01
    assert candidate['cases'][1]['raw_absolute_error_eV_fu'] is None
    assert candidate['cases'][1]['absolute_error_interval_eV_fu'] is None
    assert not report['full_panel_strict_advantage']
    assert report['comparisons'][1]['best_available_baseline_oracle_gain_interval_eV_fu'] is None


def test_unavailable_strong_reference_is_not_evidence_it_performs_worse():
    forecasts, references = panel()
    forecasts[6] = replace(forecasts[6], prediction_eV_fu=None, abstention_reason='unmeasured curvature')
    report = score(forecasts, references)
    row = report['comparisons'][0]
    assert row['strict_gain_over_all_available_baselines']
    assert not row['complete_registered_panel_strict_gain']
    assert not report['full_panel_strict_advantage']


def test_all_candidate_abstentions_have_unknown_error_not_perfect_score():
    forecasts, references = panel()
    forecasts[:2] = [replace(f, prediction_eV_fu=None, abstention_reason='unresolved') for f in forecasts[:2]]
    report = score(forecasts, references)
    assert report['models'][0]['coverage'] == 0
    assert report['models'][0]['raw_maximum_absolute_error_eV_fu'] is None
    assert report['models'][0]['maximum_absolute_error_interval_eV_fu'] is None
    assert not report['full_panel_strict_advantage']


@pytest.mark.parametrize('fault', ['missing', 'duplicate', 'extra_case', 'extra_model',
    'contract', 'reference_contract', 'duplicate_reference', 'changed_model'])
def test_incomplete_or_mixed_panel_rejected(fault):
    forecasts, references = panel()
    if fault == 'missing': forecasts.pop()
    elif fault == 'duplicate': forecasts.append(forecasts[0])
    elif fault == 'extra_case': forecasts.append(replace(forecasts[0], case_id='posthoc'))
    elif fault == 'extra_model': forecasts.append(replace(forecasts[0], model_id='posthoc'))
    elif fault == 'contract': forecasts[0] = replace(forecasts[0], comparison_contract_sha256='d'*64)
    elif fault == 'reference_contract': references[0] = replace(references[0], comparison_contract_sha256='d'*64)
    elif fault == 'duplicate_reference': references.append(references[0])
    elif fault == 'changed_model': forecasts[1] = replace(forecasts[1], frozen_model_sha256='d'*64)
    with pytest.raises(ValueError): score(forecasts, references)


@pytest.mark.parametrize('value', [True, np.nan, np.inf, '0.2', 1+2j])
def test_nonphysical_numeric_types_rejected(value):
    with pytest.raises(ValueError): absolute_error_interval(value, [.1, .2])
    with pytest.raises(ValueError): AuditedBarrierReference('case', value, .01, CONTRACT, 'b'*64)


@pytest.mark.parametrize('interval', [[.2, .1], [np.nan, .2], [.1], [True, .2]])
def test_invalid_reference_intervals_rejected(interval):
    with pytest.raises(ValueError): shared_reference_gain_interval(.2, [.3], interval)


def test_no_silent_zero_or_negative_prediction_clipping():
    with pytest.raises(ValueError): FrozenChannelForecast('B5', 'case', None, None, CONTRACT, 'b'*64)
    with pytest.raises(ValueError): FrozenChannelForecast('B5', 'case', -.1, None, CONTRACT, 'b'*64)
    record = FrozenChannelForecast('B5', 'case', -.1, 'negative signed gap', CONTRACT, 'b'*64)
    assert record.prediction_eV_fu == -.1


def test_computational_roundoff_guard_is_not_a_DFT_error_bound():
    forecasts, references = panel()
    references = [replace(r, error_bound_eV_fu=0) for r in references]
    forecasts = [replace(f, prediction_eV_fu=references[i % 2].barrier_eV_fu+5e-13)
                 if f.model_id != 'B5' else f for i, f in enumerate(forecasts)]
    report = score(forecasts, references)
    assert not report['full_panel_strict_advantage']
    assert report['floating_tolerance_eV_fu'] == 1e-12


def test_finite_input_arithmetic_overflow_cannot_be_reported_as_valid_error():
    with pytest.raises(ValueError, match='overflow'):
        AuditedBarrierReference('case', 1e308, 1e308, CONTRACT, 'b'*64)
    with pytest.raises(ValueError, match='overflow'):
        absolute_error_interval(1e308, [-1e308, 1e308])
    with pytest.raises(ValueError, match='overflow'):
        shared_reference_gain_interval(1e308, [-1e308], [-1e308, 1e308])


def test_bounded_experiment_rejects_false_harmonic_or_weak_reference_advantage():
    from scripts.check_prediction_evaluation import check
    report = check()
    assert report['random_shared_reference_intervals_checked'] == 48
    assert report['non_grid_aligned_cusp_checked']
    assert report['cusp_grid_extremum_distance_eV_fu'] > 0
    assert not report['example_complete_panel']['full_panel_strict_advantage']
    assert not report['HfO2_training_or_holdout_read']
    assert not report['B2_to_B5_material_models_frozen']
