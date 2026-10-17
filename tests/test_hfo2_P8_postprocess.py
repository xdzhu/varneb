"""Synthetic score/admission checks, never actual P8 test labels."""
from pathlib import Path
from copy import deepcopy
import importlib.util
import json
import pytest

BASE=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/G3_P8_execution_E266_20261018'
spec=importlib.util.spec_from_file_location('E266_frozen_score_tests',BASE/'postprocess.py')
post=importlib.util.module_from_spec(spec);spec.loader.exec_module(post)
RAW='28935677|COMPLETED|0:0|3600|32|node13|hfacnormal01|1\n28935679|COMPLETED|0:0|1800|32|node40|hfacnormal01|1\n'


def synthetic():
    forecasts=json.loads((BASE/'forecasts.json').read_bytes())
    values={p['id']:p['prediction_raw_eV_cell'] for p in forecasts['new_models'][0]['predictions']}
    native=dict(complete_cohort_audited=True,forecast_sha256=post.FORECAST_SHA,
        rows=[dict(id=i,label=dict(energy_change_eV_cell=values[i])) for i in range(1,9)],
        native_SCFs_reparsed=16,allocation_core_hours=17.)
    return forecasts,native


def test_synthetic_score_retains_all23_arms_184_points_and_the_frozen_refusal():
    forecasts,native=synthetic();result=post.score(forecasts,native)
    assert len(result['arms'])==23 and sum(len(r['points']) for r in result['arms'])==184
    assert result['arms'][0]['raw']['max_absolute_meV_fu']==0
    assert result['arms'][2]['model_abstention_ids']==[8]
    assert result['primary_four_arm_common_ids']==list(range(1,8))
    assert all(r['denominator']==8 for r in result['arms'])
    assert result['old_2D_refusal_preserved'] and not result['stable_sheet_TS_resolved_advantage_or_JCTC_readiness_claimed']


def test_missing_native_label_never_becomes_zero_error_or_shrunk_denominator():
    forecasts,native=synthetic();native['rows'][0]['label']=None
    result=post.score(forecasts,native);row=result['arms'][0]
    assert row['denominator']==8 and row['valid_labels']==7 and row['raw']['n']==7
    assert row['missing_label_ids']==[1] and row['points'][0]['signed_raw_error_eV_cell'] is None
    assert 1 not in result['primary_four_arm_common_ids']


def test_signed_cell_to_formula_unit_errors_and_old_calibration_not_refitted():
    forecasts,native=synthetic()
    for row in native['rows']:row['label']['energy_change_eV_cell']+=.004
    result=post.score(forecasts,native)
    assert result['arms'][0]['raw']['max_absolute_meV_fu']==pytest.approx(1.)
    assert result['arms'][0]['points'][0]['signed_raw_error_eV_cell']==pytest.approx(-.004)
    assert result['no_refit_step_selection_or_rebias']


@pytest.mark.parametrize('change',[
    lambda x:x.update(complete_cohort_audited=False),
    lambda x:x.update(forecast_sha256='0'*64),
    lambda x:x.update(rows=x['rows'][:-1])])
def test_wrong_or_incomplete_native_scope_refuses(change):
    forecasts,native=synthetic();change(native)
    with pytest.raises(ValueError,match='complete native'):post.score(forecasts,native)


def test_both_exact_terminal_allocations_required_and_cost_retained():
    result=post.accounting(RAW)
    assert sum(r['allocation_core_hours'] for r in result)==48
    with pytest.raises(ValueError):post.accounting(RAW.replace('COMPLETED','RUNNING',1))
    with pytest.raises(ValueError):post.accounting(RAW.replace('|32|','|64|',1))
    with pytest.raises(ValueError):post.accounting(RAW.splitlines()[0])
