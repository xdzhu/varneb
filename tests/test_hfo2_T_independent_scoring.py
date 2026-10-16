"""Fixed scoring and native-audit refusal logic; toy labels are not DFT evidence."""
import importlib.util
import json
from pathlib import Path
import math
import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_T_independent_E258_20261017'


@pytest.fixture(scope='module')
def setup():
    s=importlib.util.spec_from_file_location('original_T_scoring',CASE/'postprocess.py')
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    return m,m.pinned(CASE/'config.json',m.CONFIG_SHA)


def test_every_original_standard_reference_and_step_is_scored(setup):
    m,c=setup;arms=m.models(c)
    assert len(arms)==19 and len({a['name'] for a in arms})==19
    assert [a['name'] for a in arms[:5]]==list(c['points'][0]['models'])
    assert [a['name'] for a in arms[5:]]==[f"{r['reference']}@{r['step_A']:.2f}" for r in c['reference_predictions']]
    for arm in arms:
        assert [p['id'] for p in arm['points']]==list(range(1,9))
        assert [p['q_A'] for p in arm['points']]==[p['q_A'] for p in c['points']]


def test_missing_labels_are_not_zero_error_or_success(setup):
    m,c=setup;r=m.score_labels(c,{})
    assert r['T_denominator']==r['preserving_denominator']==8
    assert r['total_original_independent_inventory']==16 and r['preserving_valid_labels']==0
    assert r['all_19_arm_common_ids']==[]
    assert all(a['calibrated'] is None and a['uncalibrated'] is None for a in r['arms'])
    assert all(a['missing_label_ids']==list(range(1,9)) and not a['scored_ids'] for a in r['arms'])
    assert not r['resolved_advantage_TS_or_JCTC_readiness_claimed']


def test_error_sign_units_and_RMSE_are_independent_of_material_labels(setup):
    m,_=setup;r=m.metrics([.001,-.003])
    assert r['mean_signed_meV_cell']==-1 and r['max_absolute_meV_cell']==3
    assert r['RMSE_meV_cell']==pytest.approx(math.sqrt(5))
    assert r['max_absolute_meV_per_fu']==.75 and r['RMSE_meV_per_fu']==pytest.approx(math.sqrt(5)/4)
    assert m.metrics([]) is None
    with pytest.raises(ValueError,match='finite'):m.metrics([float('nan')])


def test_calibrated_and_raw_errors_remain_separate_without_rebias(setup):
    m,c=setup;p=c['points'][0]['models']['original_H_0.01']
    toy_truth=p['seen_origin_calibrated_energy_change_eV_cell']+.002
    result=m.score_labels(c,{1:toy_truth})
    arm=next(a for a in result['arms'] if a['name']=='original_H_0.01')
    row=arm['points'][0]
    assert row['signed_calibrated_error_eV_cell']==pytest.approx(-.002)
    assert row['signed_raw_error_eV_cell']==pytest.approx(p['value']-toy_truth)
    assert arm['origin_bias_eV_cell']==p['origin_bias_eV_cell']
    assert result['observed_origin_energy_change_eV_cell']==c['standard_predictions']['observed_origin_energy_change_eV_cell']
    assert arm['calibrated']['max_absolute_meV_per_fu']==pytest.approx(.5)
    assert arm['denominator']==8 and len(arm['missing_label_ids'])==7


def test_all_arm_common_set_and_full_geometry_support_are_fixed_without_labels(setup):
    m,c=setup;result=m.score_labels(c,{i:0. for i in range(1,9)})
    assert result['all_19_arm_common_ids']==[1,2,3,4]
    assert result['full_support_common_label_ids']==list(range(1,9))
    assert len(result['full_prediction_support_arm_names'])==9
    assert all('T_all_optical@'+step in result['full_prediction_support_arm_names'] for step in ('0.01','0.02'))
    assert all('T_two_complete_optical_doublets@'+step in result['full_prediction_support_arm_names'] for step in ('0.01','0.02'))
    assert all(a['common_all_arms']['n']==4 for a in result['arms'])


def test_geometry_abstention_keeps_formal_value_but_no_scored_error(setup):
    m,c=setup;triplet=next(a for a in m.models(c) if a['name']=='T_geometric_triplet@0.01')
    unsupported=next(p for p in triplet['points'] if p['calibrated_eV_cell'] is None)
    result=m.score_labels(c,{unsupported['id']:0.})
    arm=next(a for a in result['arms'] if a['name']==triplet['name'])
    row=next(p for p in arm['points'] if p['id']==unsupported['id'])
    assert math.isfinite(row['formal_raw_eV_cell']) and row['signed_calibrated_error_eV_cell'] is None
    assert arm['calibrated'] is None and not arm['scored_ids']
    assert unsupported['id'] in arm['model_abstention_ids']


@pytest.mark.parametrize('bad',[{0:0.},{9:0.},{1:float('inf')}])
def test_wrong_target_or_nonfinite_label_is_refused(setup,bad):
    m,c=setup
    with pytest.raises(ValueError):m.score_labels(c,bad)


def test_accounting_requires_both_original_complete_allocations(setup):
    m,_=setup
    raw='28913247|COMPLETED|0:0|800|32|node3|hfacnormal01|1\n28913248|COMPLETED|0:0|900|32|node3|hfacnormal01|1\n'
    result=m.accounting(raw)
    assert sum(r['allocation_core_hours'] for r in result)==pytest.approx(1700*32/3600)
    for changed in (raw.replace('COMPLETED','RUNNING',1),raw.replace('|32|','|1|',1),raw.splitlines()[0],
                    raw.replace('hfacnormal01','wrong',1),raw.replace('0:0','1:0',1)):
        with pytest.raises(ValueError):m.accounting(changed)


def test_live_native_audit_refuses_before_worker_replay_or_any_DFT(setup,monkeypatch):
    m,_=setup
    raw='28913247|RUNNING|0:0|800|32|node3|hfacnormal01|1\n28913248|RUNNING|0:0|900|32|node3|hfacnormal01|1\n'
    monkeypatch.setattr(m.subprocess,'check_output',lambda *a,**kw:raw)
    with pytest.raises(ValueError,match='complete normally'):m.native_audit(CASE)


def test_changed_pinned_source_is_refused_without_refit(setup,tmp_path):
    m,_=setup;p=tmp_path/'modified.json';p.write_text(json.dumps(dict(changed=True)))
    with pytest.raises(ValueError,match='strict source'):m.pinned(p,m.CONFIG_SHA)
