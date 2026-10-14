"""Complete measured E159 data replay; no new DFT, stability or TS assertion."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np
import pytest

from scripts import hfo2_G3_conditional_canary as canary
from scripts import summarize_hfo2_G3_curvature as measured

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/hfo2_channels/20261008'
CASE=BASE/'G3_complete_curvature_E159_20261014'
PINS=dict(T='a182e4c386889880c77c2069e939cd65665724d98d3838e251a077204cad1c3b',
          preserving='5d137fd0e676ad4fb7e87076d76b53e618ab9de25bc86228e02d4bcf0285a778')
CONFIGS=dict(T='f77ebe3708f89b14dbbda053505f56fbf26672a2830bd0231680d6798d59a739',
             preserving='e8ce5429c24fb443f4a651eef14aac8cf0190a7399cb19a0abef838658be12b5')


def arithmetic_same(actual,expected):
    """Only derived real arithmetic has a BLAS-roundoff allowance."""
    assert type(actual) is type(expected)
    if isinstance(expected,dict):
        assert actual.keys()==expected.keys()
        for key in expected:arithmetic_same(actual[key],expected[key])
    elif isinstance(expected,list):
        assert len(actual)==len(expected)
        for a,e in zip(actual,expected):arithmetic_same(a,e)
    elif isinstance(expected,float):assert actual==pytest.approx(expected,abs=2e-12,rel=0)
    else:assert actual==expected


@pytest.fixture(scope='module')
def material():
    data=(BASE/'G3_curvature_E150_20261014/manifest.json').read_bytes()
    assert hashlib.sha256(data).hexdigest()==canary.PARENT_SHA
    reports={}
    with zipfile.ZipFile(CASE/'native_reports.zip') as archive:
        assert set(archive.namelist())=={channel+'_native.json' for channel in PINS}
        for channel,pin in PINS.items():
            raw=archive.read(channel+'_native.json')
            assert hashlib.sha256(raw).hexdigest()==pin
            reports[channel]=json.loads(raw)
    return json.loads(data),reports


def test_both_original_amplitudes_and_no_native_failures(material):
    manifest,reports=material
    for channel,report in reports.items():
        assert report['raw_SCFs_reparsed']==report['completed_records']==148
        assert report['conservative_attempt_reservations']==report['observed_native_SCF_directories']==148
        assert report['failures_present']==[]
        assert [(b['step_A'],b['complete_pairs']) for b in report['blocks']]==[(.01,37),(.02,37)]
        assert all(b['full_stencil_at_this_step'] for b in report['blocks'])
        assert report['center_source']==manifest['centers'][channel]['raw_source']
        for source in report['native_record_sources']:
            validation=source['validation']
            assert validation['genuine_MPI_ranks']==32 and validation['raw_SCFs']==1
            assert {name:validation['input_sha256'][name] for name in canary.study.CONTRACT}==canary.study.CONTRACT
        assert not report['physical_error_bound_or_TS_certificate']
        assert not report['physical_reserved_geometry_generated_or_read']


def test_complete_summary_replays_without_filling_or_reparsing(material):
    manifest,reports=material
    result=measured.summarize(manifest,reports)
    expected=json.loads((CASE/'complete_summary.json').read_bytes())
    arithmetic_same(result,{k:v for k,v in expected.items() if k!='input_sha256'})
    assert expected['input_sha256']==dict(manifest=canary.PARENT_SHA,**PINS)
    assert all(row['complete_amplitudes']==[.01,.02] and row['two_step_comparison_available'] for row in result['rows'])
    assert not result['new_raw_SCF_reparse'] and result['new_DFT_or_scheduler_operations']==0


@pytest.mark.parametrize('channel',['T','preserving'])
def test_real_complete_canary_policy_and_same_positive_complement_screen(material,channel):
    manifest,reports=material
    raw=(CASE/(channel+'_config.json')).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==CONFIGS[channel]
    config=json.loads(raw)
    _,_,_,matrices=canary.validate_config(config,manifest,channel,reports[channel],PINS[channel])
    axes=np.asarray(config['internal_axes_36_by_2'])
    complement=np.linalg.svd(axes.T,full_matrices=True)[2][2:].T
    assert all(np.count_nonzero(np.linalg.eigvalsh(h)<0)==2 for h in matrices)
    np.testing.assert_allclose([np.linalg.eigvalsh(complement.T@h@complement)[0] for h in matrices],
        config['raw_complement_minima_eV_A2'],atol=2e-12,rtol=0)
    assert min(config['raw_complement_minima_eV_A2'])>0
    assert config['settings']['maximum_attempts_per_arm']==12
    assert config['settings']['maximum_attempts_per_center']==24
    assert config['settings']['modeling_nodes_used']==1 and config['settings']['independent_test_points_used']==0
    assert not config['full_internal_index_or_TS_certified']
    assert not config['final_surface_or_forecast_domain_locked']


@pytest.mark.parametrize('channel',['T','preserving'])
def test_one_native_tangent_still_leaves_raw_transverse_instability(material,channel):
    manifest,reports=material
    row=next(r for r in measured.summarize(manifest,reports)['rows'] if r['channel']==channel)
    for block in row['blocks']:
        assert block['raw_negative_count']==2
        assert block['tangent_complement_raw_negative_count']==1
        assert block['minimum_retained_dimension_for_raw_positive_complement']==2
        assert not block['physical_error_bound_or_stationary_index_or_TS_certificate']


def test_second_amplitude_keeps_all_strong_references_and_optical_control():
    result=json.loads((CASE/'reference_complements_complete.json').read_bytes())
    assert len(result['rows'])==28 and result['native_report_sha256']==PINS
    expected=[('T_geometric_triplet',3),('T_two_complete_optical_doublets',4),('T_all_optical',33),
              *((f'Cmma_frame_{i}',4) for i in (103,106,114,127))]
    for channel in PINS:
        for step in (.01,.02):
            rows=[r for r in result['rows'] if r['channel']==channel and r['step_A']==step]
            assert [(r['reference'],r['retained_rank']) for r in rows]==expected
            for row in rows:
                joint=row['joint_release']
                assert joint['released_dimension']==36-row['retained_rank']
                if row['reference']=='T_all_optical':
                    assert joint['raw_negative_count']==0 and joint['raw_minimum_eigenvalue_eV_A2']>2.9
                    assert joint['raw_positive_block_Schur_retained_negative_count']==2
                else:assert joint['raw_negative_count']>=1
    assert not result['raw_or_transformed_author_arrays_exported']
    assert not result['source_frequencies_used_as_PBE_curvature']
    assert not result['physical_sign_floor_error_bound_or_stationary_response_certified']


def test_missing_native_pair_cannot_be_promoted_to_a_canary(material):
    manifest,reports=material;report=deepcopy(reports['T'])
    report['blocks'][1]['complete_pairs']=36
    with pytest.raises(ValueError,match='both original37-column'):
        canary.prepare(manifest,'T',report,PINS['T'])
