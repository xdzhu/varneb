"""Actual 36 portable training frames and independent native work contractions."""
from pathlib import Path
import numpy as np
import pytest
from ase.io import read
from scripts.analyze_hfo2_clamped_path_work import analyze,partition,path_arrays
from scripts.analyze_hfo2_clamped_reference import terminal_images,original_chart_images,pinned_json
from scripts.analyze_hfo2_clamped_external_work import SEEDS,SEEDS_SHA
from scripts.analyze_hfo2_cmma_path_reference import MASS
from vcneb.path_work import sampled_path_work
from scripts.audit_hfo2_static_replica import sha256
import json

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/hfo2_channels/20261008'
CASES=[
    ('clamped_M_terminal_E068_20261011','56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c'),
    ('clamped_T_terminal_E079_20261011','c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e'),
    ('clamped_reversing_terminal_E083_20261011','632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41'),
    ('clamped_plus1_T_terminal_E087_20261011','32a6e49a618e6e62632ae71a5e0ed093095254f6858221be2ade3d83cc194a28'),
]


def test_unknown_or_holdout_receipt_rejects_before_any_source_access():
    with pytest.raises(ValueError,match='before label access'):
        analyze(None,None,None,None,None,terminal_audit_sha256='0'*64,
                clamped_seed_manifest=None,clamped_seed_manifest_sha256=None)


@pytest.mark.parametrize('case,digest',CASES)
def test_actual_physical_native_work_equals_dual_reference_atom_and_cell_work(case,digest):
    images,report,_,_=terminal_images(BASE/case,digest)
    seeds=pinned_json(SEEDS,SEEDS_SHA)
    t=read(BASE/'reference_variants/T.vasp')
    original=original_chart_images(images,seeds['orientation'],t)
    rotation=np.asarray(seeds['orientation']['cartesian_rotation_rows'])
    x,g,energy,direct,_,chart=path_arrays(images,original,t,rotation,report)
    random=np.linalg.qr(np.random.default_rng(92).normal(size=(36,4)))[0]
    parts,diagnostic=partition(random,MASS,chart.cell_dofs)
    result=sampled_path_work(x,g,energy,parts)
    np.testing.assert_allclose(result['segment_full_work'],[r['total_eV_per_control'] for r in direct],atol=1e-12,rtol=0)
    np.testing.assert_allclose(result['segment_partition_work']['open_cell'],[r['cell_eV_per_control'] for r in direct],atol=1e-12,rtol=0)
    atom=np.sum([result['segment_partition_work'][n] for n in ('retained_atomic','omitted_atomic','rigid_translation')],axis=0)
    np.testing.assert_allclose(atom,[r['atomic_eV_per_control'] for r in direct],atol=1e-12,rtol=0)
    assert diagnostic['retained_rank']==4 and diagnostic['omitted_atomic_rank']==29
    assert all(a.calc is None for a in original)
    assert all(a.calc.__class__.__name__=='SinglePointCalculator' for a in images)
    assert result['physical_error_bound_measured'] is False


def test_translation_only_registered_span_is_rejected_not_deleted():
    translation=np.zeros((36,1));translation[::3,0]=np.sqrt(MASS)/np.linalg.norm(np.sqrt(MASS))
    with pytest.raises(ValueError,match='rank lost'):partition(translation,MASS,3)


@pytest.mark.parametrize('name,job,peak,energy_minus_work',[
    ('zero_M','28722320',3,.08359097282891609),
    ('zero_T','28692775',4,.03384274886041727),
    ('zero_reversing','28692776',4,.14110129045160047),
    ('plus1_T','28709788',4,.03767302780244475),
])
def test_actual_export_keeps_signed_parts_all_strong_frames_and_observed_coarse_residual(name,job,peak,energy_minus_work):
    result=json.loads((BASE/'clamped_path_work_E092_20261011'/(name+'.json')).read_text())
    assert result['source_job_id']==job and result['sampled_peak_image']==peak
    assert result['script_sha256']==sha256(ROOT/'scripts/analyze_hfo2_clamped_path_work.py')
    assert result['core_source_sha256']==sha256(ROOT/'vcneb/path_work.py')
    assert result['all_E023_Cmma_frames_retained'] is True
    assert result['portable_native_frames_reparsed']==9
    assert result['complete_six_physical_bytes_rechecked_locally'] is False
    assert result['new_DFT_calls']==result['scheduler_mutations']==0
    assert result['G3_selection_or_holdout_access'] is False
    expected=['T_geometric_triplet','T_two_lowest_complete_optical_doublets',
              'Cmma_frame_103','Cmma_frame_106','Cmma_frame_114','Cmma_frame_127']
    assert list(result['representations'])==expected
    first=None
    for key,span in result['representations'].items():
        assert span['retained_rank']==(3 if key==expected[0] else 4)
        assert span['source_modes_or_forces_ASR_corrected'] is False
        assert span['physical_error_bound_measured'] is False
        assert span['independent_native_work_max_defect_eV_cell']<1e-12
        assert span['independent_cell_work_max_defect_eV_cell']<1e-12
        assert span['cumulative_energy_minus_work'][peak]==pytest.approx(energy_minus_work,abs=1e-12)
        parts=np.sum(list(span['segment_partition_work'].values()),axis=0)
        np.testing.assert_allclose(parts,span['segment_full_work'],atol=1e-12,rtol=0)
        if first is not None:
            np.testing.assert_array_equal(span['segment_full_work'],first['segment_full_work'])
            np.testing.assert_array_equal(span['segment_partition_work']['open_cell'],first['segment_partition_work']['open_cell'])
        first=span
