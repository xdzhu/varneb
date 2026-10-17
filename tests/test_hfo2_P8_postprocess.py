"""Synthetic refusal tests plus the externally pinned completed P8 replay."""
from pathlib import Path
from copy import deepcopy
import importlib.util
import json
import hashlib
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


def actual():
    paths={'native_complete_E267.json':'dd245ac8f7335442ef5355b484b5c01dcba03ae3201ec1ad02131a3ae16ac654',
           'score_E267.json':'c51dec44d4e92e05ad824bd7952a9e1cee16d662476684675c9ae0c7489c2607',
           'forecasts.json':post.FORECAST_SHA}
    result={}
    for name,digest in paths.items():
        raw=(BASE/name).read_bytes();assert hashlib.sha256(raw).hexdigest()==digest
        result[name]=json.loads(raw)
    return result


def test_actual_complete_native24_cost_inputs_MPI_and_both_gradients():
    values=actual();native=values['native_complete_E267.json']
    assert native['complete_cohort_audited'] and native['valid_candidate_labels']==8
    assert len(native['native_records'])==native['native_SCFs_reparsed']==24
    assert [r['native_SCFs'] for r in native['rows']]==[2,2,2,2,6,4,3,3]
    assert native['allocation_core_hours']==pytest.approx((1069+2075)*32/3600)
    config=json.loads((BASE/'config.json').read_bytes())
    for record in native['native_records']:
        assert record['validation']['genuine_MPI_ranks']==32
        assert all(record['validation']['input_sha256'][k]==v for k,v in config['physical_input_sha256'].items())
    for row in native['rows']:
        assert row['label']['released_gradient_norm_eV_A']<=.03
        assert abs(row['label']['zeta_gradient_eV_A'])<=1e-4
        assert not row['label']['stability_certified']


def test_actual_all23_score_replays_exactly_without_refit_or_hidden_selection():
    values=actual();stored=values['score_E267.json']
    recomputed=post.score(values['forecasts.json'],values['native_complete_E267.json'])
    assert recomputed=={k:v for k,v in stored.items() if k not in ('source_audit_sha256','forecast_sha256','source_LF_sha256')}
    assert sum(len(a['points']) for a in recomputed['arms'])==184
    assert recomputed['primary_four_arm_common_ids']==list(range(1,8))
    assert recomputed['all23_common_ids']==[]


def test_actual_max_error_and_RMSE_do_not_support_a_metric_independent_nonlinear_win():
    score=actual()['score_E267.json'];a=score['arms']
    assert [r['prediction_support'] for r in a[:4]]==[8,8,7,8]
    for harmonic,nonlinear in ((a[0],a[1]),(a[2],a[3])):
        assert harmonic['raw_on_primary_common_ids']['max_absolute_meV_fu'] < nonlinear['raw_on_primary_common_ids']['max_absolute_meV_fu']
        assert harmonic['raw_on_primary_common_ids']['RMSE_meV_fu'] > nonlinear['raw_on_primary_common_ids']['RMSE_meV_fu']
    assert not score['stable_sheet_TS_resolved_advantage_or_JCTC_readiness_claimed']


def test_actual_main_and_complete_appendix_numbers_match_frozen_scores():
    paper=(BASE.parents[3]/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    score=actual()['score_E267.json']
    def table_after(caption):
        text=paper.split(caption,1)[1];lines=text.splitlines();start=next(i for i,l in enumerate(lines) if l.startswith('|'))
        table=[]
        for line in lines[start:]:
            if not line.startswith('|'):break
            table.append([c.strip() for c in line.strip('|').split('|')])
        return table[2:]
    main=table_after('Table 5. Revised preserving errors')
    assert len(main)==4
    for printed,arm in zip(main,score['arms'][:4]):
        assert printed[2]==f"{arm['prediction_support']}/8"
        assert list(map(float,printed[-2:]))==pytest.approx([
            arm['raw_on_primary_common_ids']['max_absolute_meV_fu'],
            arm['raw_on_primary_common_ids']['RMSE_meV_fu']],abs=5.01e-7)
    appendix=table_after('Table A19. All23frozen arms')
    assert len(appendix)==23
    for printed,arm in zip(appendix,score['arms']):
        assert int(printed[1])==arm['raw']['n']
        expected=[arm[m][k] for m in ('raw','origin_aligned') for k in ('max_absolute_meV_fu','RMSE_meV_fu')]
        for cell,value in zip(printed[2:],expected):
            assert cell=='--' if value is None else float(cell)==pytest.approx(value,abs=5.01e-7)
