"""Actual four-point energies, forces, provenance and bounded manuscript claims."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_soft_probe_E256_20261016'


@pytest.fixture(scope='module')
def evidence():
    pins={
        'native_complete_audit.json':'41321300ae00dc46a526ea14069fb541f6d02811fcac8dcccdf5fda72ace9928',
        'analysis_E257.json':'761bb4bd2784fca366cfbb32abc4800647b8545e8281b5ccd03cf93028ceb318',
        'config.json':'f7dc1053bff96bbd141d0a7e2e5c14d949f91f69832852d61bf7af9b932de7a7',
        'soft_direction.json':'1e5dcd01ec5c7731646c9e1f6b0f230fc8bcf45bbcad5a785125534a0d2ee065',
    }
    for name,digest in pins.items():
        assert hashlib.sha256((CASE/name).read_bytes()).hexdigest()==digest
    audit=json.loads((CASE/'native_complete_audit.json').read_bytes())
    analysis=json.loads((CASE/'analysis_E257.json').read_bytes())
    config=json.loads((CASE/'config.json').read_bytes())
    soft=json.loads((CASE/'soft_direction.json').read_bytes())
    spec=importlib.util.spec_from_file_location('four_statics_postprocess',CASE/'postprocess.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    text=(CASE/'postprocess.py').read_bytes().replace(b'\r\n',b'\n')
    assert hashlib.sha256(text).hexdigest()==audit['audit_text_sha256']
    return audit,analysis,config,soft,module


def test_all_four_native_records_keep_original_EFS_MPI_and_source_spelling(evidence):
    audit,_,config,_,_=evidence
    assert audit['new_probe_SCFs_reparsed']==len(audit['records'])==len(audit['native_sources'])==4
    assert audit['old_center_SCFs_reparsed']==1 and audit['job_id']=='28912687'
    assert audit['accounting']['normal_completion'] and audit['accounting']['elapsed_seconds']==478
    fixed={k:v for k,v in audit['records'][0]['validation']['input_sha256'].items() if k!='STRU'}
    assert len(fixed)==6
    for p,r,s in zip(config['points'],audit['records'],audit['native_sources']):
        np.testing.assert_array_equal(p['coordinates_full_A'],r['coordinates'])
        assert r['validation']==s['validation']
        assert '/bundle/../probes/evaluations/' in s['result_path']
        assert r['validation']['genuine_MPI_ranks']==32 and r['validation']['raw_SCFs']==1
        assert {k:v for k,v in r['validation']['input_sha256'].items() if k!='STRU'}==fixed


def test_signed_lowerings_and_formula_unit_normalization_are_not_refitted(evidence):
    audit,analysis,config,soft,module=evidence
    assert module.analyze(audit,config,soft)==analysis
    expected=[-.05053285167377908,-.05026035250921268,-.0686314615450101,-.06797588139306754]
    for row,value,p,r in zip(analysis['rows'],expected,config['points'],audit['records']):
        assert row['delta_meV_cell']==pytest.approx(value,abs=2e-12,rel=0)
        assert row['delta_meV_per_fu']==row['delta_meV_cell']/4
        assert row['delta_eV_cell']==r['enthalpy_eV']-analysis['base_energy_eV_cell']
        assert row['original_two_quadratic_delta_eV_cell']==p['quadratic_delta_eV_cell']
    assert not analysis['predictions_refitted']
    assert not analysis['amplitudes_or_signs_selected_after_labels']


@pytest.mark.parametrize('index,h',[(0,.05),(1,.10)])
def test_energy_and_gradient_secants_remain_distinct_observables(evidence,index,h):
    _,a,_,_,_=evidence
    plus=next(r for r in a['rows'] if r['t_A']==h);minus=next(r for r in a['rows'] if r['t_A']==-h)
    pair=a['pairs'][index]
    ke=(plus['delta_eV_cell']+minus['delta_eV_cell'])/(h*h)
    kg=(plus['gradient_along_direction_eV_A']-minus['gradient_along_direction_eV_A'])/(2*h)
    assert pair['energy_central_curvature_eV_A2']==ke
    assert pair['gradient_central_curvature_eV_A2']==kg
    assert ke<0
    assert (kg<0) if index==0 else (kg>0)  # not coerced into a harmonic identity


def test_retrospective_adjacent_work_does_not_hide_center_interval_discrepancy(evidence):
    _,a,_,_,_=evidence
    points=[(r['t_A'],r['delta_eV_cell'],r['gradient_along_direction_eV_A']) for r in a['rows']]
    points.append((0.,0.,a['base_soft_gradient_eV_A']));points.sort()
    residual=[]
    for left,right in zip(points,points[1:]):
        residual.append(((right[1]-left[1])-(right[0]-left[0])*(left[2]+right[2])/2)*1000)
    np.testing.assert_allclose(residual,[.00011524165910889689,.029244940950440206,
        -.029491290548395308,-.0004607514555066109],atol=2e-12,rtol=0)
    assert not a['physical_error_bound_stable_branch_TS_or_independent_accuracy_certified']


def test_literal_path_alias_is_kept_only_after_exact_resolved_identity(tmp_path,evidence):
    module=evidence[-1];namespace=tmp_path/'probes';namespace.mkdir()
    expected=namespace/'result.json';expected.write_text('{}',encoding='utf-8')
    (tmp_path/'bundle').mkdir()
    alias=tmp_path/'bundle'/'..'/'probes'/'result.json'
    assert module.checked_literal_path(expected,namespace,{'result_path':str(alias)})==alias
    other=tmp_path/'other.json';other.write_text('{}',encoding='utf-8')
    with pytest.raises(ValueError,match='exact expected file'):
        module.checked_literal_path(expected,namespace,{'result_path':str(other)})


@pytest.mark.parametrize('fault',['missing_record','changed_coordinate','changed_amplitude'])
def test_incomplete_or_changed_target_cannot_be_scored(evidence,fault):
    audit,_,config,soft,module=evidence;a=deepcopy(audit);c=deepcopy(config)
    if fault=='missing_record':a['records']=a['records'][:3]
    elif fault=='changed_coordinate':a['records'][0]['coordinates'][0]+=.001
    else:c['points'][0]['t_A']=.2
    with pytest.raises(ValueError):module.analyze(a,c,soft)


def test_table_A17_is_sourced_and_original_test_refusal_remains(evidence):
    _,a,_,_,_=evidence
    text=(ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    block=text.split('Table A17. Four frozen directional training statics',1)[1].split('\n### ',1)[0]
    rows=[line.strip().strip('|').split('|') for line in block.splitlines()
          if line.startswith('|') and not line.startswith('|---')][1:]
    assert len(rows)==4
    for row,source in zip(rows,a['rows']):
        assert float(row[0])==source['t_A']
        assert float(row[1])==pytest.approx(source['delta_meV_cell'],abs=5.01e-7,rel=0)
        assert float(row[2])==pytest.approx(source['gradient_along_direction_eV_A'],abs=5.01e-10,rel=0)
    assert '0.029245/-0.029491 meV/cell' in text
    section=text.split('### 3.4.',1)[1].split('## 4.',1)[0]
    assert 'its eight independent points remain uncomputed' in ' '.join(section.split()) and 'original refusal' in section
    assert a['independent_points_computed']==0 and a['original_independent_admission_refusal_preserved']
