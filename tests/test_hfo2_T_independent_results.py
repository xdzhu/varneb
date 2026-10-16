"""Actual native T labels and all-arm errors, with source-derived paper tables."""
from collections import Counter
import hashlib,importlib.util,json
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_T_independent_E258_20261017'
PAPER=ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'
NAMES=['Bilinear','Tensor quadratic','Refitted total quadratic','Original H/g, 0.01','Original H/g, 0.02',
 'T triplet, 0.01','T two optical doublets, 0.01','T all optical, 0.01',
 'Cmma 103, 0.01','Cmma 106, 0.01','Cmma 114, 0.01','Cmma 127, 0.01',
 'T triplet, 0.02','T two optical doublets, 0.02','T all optical, 0.02',
 'Cmma 103, 0.02','Cmma 106, 0.02','Cmma 114, 0.02','Cmma 127, 0.02']


@pytest.fixture(scope='module')
def evidence():
    pins={'native_independent_E259.json':'5dac550e6a802225465b0717ec9b2bb42bb55fc457e3adddc2f05fdd5d455ab5',
          'score_E259.json':'272bb11781be401e74ad333e4a211973410217ccce68a69fc9f2bd41f41a14cc',
          'config.json':'d78a39c60a6b7d2af52ab8a9cfabfe829635f4fcd94d6c682e0f63c98826506e'}
    for n,h in pins.items():assert hashlib.sha256((CASE/n).read_bytes()).hexdigest()==h
    a,s,c=[json.loads((CASE/n).read_bytes()) for n in pins]
    spec=importlib.util.spec_from_file_location('T_actual_scoring',CASE/'postprocess.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    assert hashlib.sha256((CASE/'postprocess.py').read_bytes().replace(b'\r\n',b'\n')).hexdigest()==a['audit_LF_code_sha256']
    return a,s,c,m


def test_both_actual_allocations_and_all_eighteen_native_SCFs_are_complete(evidence):
    a,_,_,_=evidence
    assert [(r['job_id'],r['state'],r['exit'],r['elapsed_seconds']) for r in a['allocations']]==[
        ('28913247','COMPLETED','0:0',837),('28913248','COMPLETED','0:0',1320)]
    assert a['valid_gradient_candidate_labels']==a['original_targets_processed']==8
    assert a['new_native_SCFs_reparsed']==len(a['native_records'])==len(a['native_sources'])==18
    assert len({p['result_path'] for p in a['native_sources']})==18
    fixed={k:v for k,v in a['native_records'][0]['validation']['input_sha256'].items() if k!='STRU'}
    assert len(fixed)==6
    for row in a['native_records']:
        assert row['validation']['genuine_MPI_ranks']==32 and row['validation']['raw_SCFs']==1
        assert {k:v for k,v in row['validation']['input_sha256'].items() if k!='STRU'}==fixed
        assert row['coordinates'][-1]==0
    assert a['allocation_core_hours']==pytest.approx((837+1320)*32/3600)
    assert a['completed_call_proxy_core_hours']==pytest.approx(sum(r['validation']['elapsed_seconds'] for r in a['native_records'])*32/3600)


def test_all_original_points_retain_full_native_labels_without_training_replacement(evidence):
    a,_,c,_=evidence
    for row,point in zip(a['rows'],c['points']):
        assert row['id']==point['id'] and row['q_A']==point['q_A']
        assert row['status']=='single_start_gradient_candidate_stability_unchecked'
        assert row['label'] is not None and 0<=row['label']['released_gradient_norm_eV_A']<=.03
        assert 1<=row['reserved_attempts']<=12 and row['new_native_SCFs']==row['reserved_attempts']
    assert sum(r['reserved_attempts'] for r in a['rows'])==18<=96
    assert a['preserving_independent_labels']==0 and a['total_original_independent_inventory']==16
    assert not a['physical_plus005_or_private_author_arrays_accessed']
    assert not a['full_stability_TS_or_whole_path_predictive_gain_certified']


def test_all_nineteen_signed_metrics_are_reconstructed_without_refit(evidence):
    a,s,c,m=evidence
    regenerated=m.score(c,a);regenerated['source_native_audit_sha256']=hashlib.sha256((CASE/'native_independent_E259.json').read_bytes()).hexdigest()
    assert regenerated==s and len(s['arms'])==19
    peak=c['standard_predictions']['raw_peak_reference_energy_eV_cell']
    labels={r['id']:r['label']['energy_eV_cell']-peak for r in a['rows']}
    for arm in s['arms']:
        for field,key in [('signed_calibrated_error_eV_cell','calibrated'),('signed_raw_error_eV_cell','uncalibrated')]:
            values=[p[field]*1000 for p in arm['points'] if p[field] is not None]
            assert arm[key]['max_absolute_meV_cell']==pytest.approx(max(abs(v) for v in values),abs=2e-12)
            assert arm[key]['RMSE_meV_per_fu']==pytest.approx(np.sqrt(np.mean(np.square(values)))/4,abs=2e-12)
        for p in arm['points']:
            assert p['observed_change_eV_cell']==labels[p['id']]
            expected=None if p['calibrated_eV_cell'] is None else p['calibrated_eV_cell']-labels[p['id']]
            assert p['signed_calibrated_error_eV_cell']==expected
    assert s['all_19_arm_common_ids']==[1,2,3,4]


def test_observed_ranking_depends_on_registered_alignment_not_model_selection(evidence):
    _,s,_,_=evidence;arms={r['name']:r for r in s['arms']}
    assert arms['tensor_quadratic']['calibrated']['max_absolute_meV_per_fu']<arms['original_H_0.01']['calibrated']['max_absolute_meV_per_fu']
    assert arms['tensor_quadratic']['uncalibrated']['max_absolute_meV_per_fu']>arms['original_H_0.01']['uncalibrated']['max_absolute_meV_per_fu']
    for step in ('0.01','0.02'):
        for metric in ('calibrated','uncalibrated'):
            assert arms['T_all_optical@'+step][metric]['max_absolute_meV_per_fu']==pytest.approx(
                arms['original_H_'+step][metric]['max_absolute_meV_per_fu'],abs=2e-12)
    assert not s['resolved_advantage_TS_or_JCTC_readiness_claimed']


def test_both_new_manuscript_tables_match_actual_fixed_scoring(evidence):
    _,s,_,_=evidence;lines=Counter(l.strip() for l in PAPER.read_text(encoding='utf-8').splitlines() if l.startswith('|'))
    for i,row in enumerate(s['arms']):
        a,u,co=row['calibrated'],row['uncalibrated'],row['common_all_arms']
        appendix=f"| {NAMES[i]} | {len(row['scored_ids'])} | {a['max_absolute_meV_per_fu']:.6f} | {a['RMSE_meV_per_fu']:.6f} | {co['max_absolute_meV_per_fu']:.6f} | {u['max_absolute_meV_per_fu']:.6f} |"
        assert lines[appendix]==1
        if i<5:
            name=NAMES[i]+(' Angstrom' if i in (3,4) else '')
            main=f"| {name} | {a['max_absolute_meV_per_fu']:.6f} | {a['RMSE_meV_per_fu']:.6f} | {u['max_absolute_meV_per_fu']:.6f} |"
            assert lines[main]==1
    assert sum(lines.values())==182
