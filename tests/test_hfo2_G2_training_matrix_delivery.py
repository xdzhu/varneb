"""Replay six real zero-DFT prepared G2 seeds; never access remote assets."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from ase.io import read

from examples.hfo2_fixed_input_factory import CONTRACT,same_ordered_geometry
from scripts.audit_hfo2_G1_gate import EXPECTED
from scripts.prepare_hfo2_clamped_chains import CHANNEL_FINAL,endpoint_evidence
from vcneb import validate_path_geometry,validate_periodic_path_lift


ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/hfo2_channels/20261008'
CASE=BASE/'clamped_G2_remaining_E056_20261010'
MATRIX=[('strain_0000','PO_to_T'),('strain_0000','PO_flip_T_pattern_reversing')]
MATRIX += [('strain_p0100',c) for c in CHANNEL_FINAL]


def load(path):return json.loads(path.read_text(encoding='utf-8'))
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize('condition,channel',MATRIX)
def test_remaining_starting_chain_keeps_registered_contract_geometry_and_cache(condition,channel):
    seed=CASE/'seeds'/condition/channel
    m=load(seed/'manifest.json')
    assert m['status']=='G2_geometry_prepared_DFT_pending' and m['channel']==channel
    assert m['strain']==(0. if condition=='strain_0000' else .01)
    assert m['physical_contract_sha256']==CONTRACT and not m['physical_inputs_changed']
    assert m['new_DFT_calls']==0 and not m['holdout_generated'] and not m['climb']
    assert (m['n_total_images'],m['n_internal_images'],m['n_fixed_endpoints'])==(9,7,2)
    assert m['cached_interior_evaluations']==0 and m['fmax_eV_A']==.10
    assert m['source_observation_sha256']==dict((c,h) for c,h,_ in EXPECTED)[channel]
    assert not m['lift_preparation']['source_results_reused'] and not m['lift_preparation']['atom_permutation_applied']
    assert all(digest(seed/p)==h for p,h in m['files_sha256'].items())
    assert m['preparation_script_sha256']==digest(CASE/'executed_source/scripts/prepare_hfo2_clamped_chains.py')
    images=read(seed/'seed.traj',index=':')
    assert len(images)==9 and all(a.calc is None for a in images)
    validate_periodic_path_lift(images)
    validate_path_geometry(images,cell_scale=m['cell_scale_A'],minimum_distance=1.6,maximum_deformation=.25)
    for a in images:
        assert a.get_chemical_symbols()==['Hf']*4+['O']*8
        np.testing.assert_allclose(a.cell.array[:2],np.asarray(m['reference_cell_A'])[:2],atol=1e-10,rtol=0.)
    for a,phase in zip((images[0],images[-1]),('PO_plus',CHANNEL_FINAL[channel])):
        assert same_ordered_geometry(a,endpoint_evidence(BASE,condition,phase)[1])
    p=load(seed/'factory_parameters.json')
    assert p['seed_cache_records'][1:8]==[None]*7 and p['seed_cache_records'][0] is not None
    g=load(seed/'vcneb_preflight_HF.json')
    assert g['n_images']==9 and g['n_interior_images']==7 and g['requires_stress']
    assert g['mechanical_boundary']['kind']=='clamped_plane' and g['mechanical_boundary']['cell_dofs']==3
    assert not g['climbing_image_requested'] and g['fmax_target_eV_per_A']==.10
    assert g['calculator_validation']=='not_instantiated_validate_only'


def test_receipt_requires_actual_six_channel_cache_checks_and_no_submission():
    receipt=load(CASE/'preparation_receipt.json')
    assert receipt['prepared_chains']==6 and receipt['existing_independent_G2_chains']==2
    assert receipt['total_registered_G2_matrix_chains']==8 and receipt['new_independent_DFT_chains_started']==0
    assert receipt['new_DFT_calls']==0 and not receipt['physical_parameters_changed']
    assert not receipt['old_source_or_jobs_mutated'] and not receipt['held_out_condition_generated_or_labels_read']
    assert receipt['external_executable_launch_prohibited_during_preparation']
    assert receipt['execution_source_code_files_byte_checked']>100
    assert {(r['condition'],r['channel']) for r in receipt['prepared_channels']}==set(MATRIX)
    for r in receipt['prepared_channels']:
        seed=CASE/'seeds'/r['condition']/r['channel']
        assert digest(seed/'manifest.json')==r['manifest_sha256']
        assert digest(seed/'vcneb_preflight_HF.json')==r['geometry_preflight_sha256']
        assert r['readiness']=='seed_and_cache_preflight_passed_NOT_submitted'
        assert r['all_six_original_physical_bytes_and_native_raw_EFS_checked'] and r['DFT_calls']==0
        assert len(r['cached_endpoints'])==2
        for point in r['cached_endpoints']:
            path=CASE/'cache_preflight'/r['condition']/r['channel']/f"image_{point['image_index']:04d}/seed_cache_audit.json"
            assert digest(path)==point['cache_audit_sha256'] and point['new_DFT_calls']==0
            assert load(path)['raw_source']==point['actual_raw_source']
    assert not any(p.suffix in ('.upf','.orb') for p in CASE.rglob('*') if p.is_file())


def test_preparation_and_first_pilots_complete_only_the_finite_training_start_matrix():
    prior=BASE/'clamped_G2_E053_20261010/strain_0000'
    all_seeds={}
    for condition,channel in MATRIX:
        all_seeds[condition,channel]=load(CASE/'seeds'/condition/channel/'manifest.json')
    for channel in ('PO_to_M','PO_flip_T_pattern_preserving'):
        all_seeds['strain_0000',channel]=load(prior/channel/'manifest.json')
    assert set(all_seeds)=={(s,c) for s in ('strain_0000','strain_p0100') for c in CHANNEL_FINAL}
    cells={}
    for (condition,channel),m in all_seeds.items():
        plane=np.asarray(m['reference_cell_A'])[:2]
        if condition in cells:np.testing.assert_allclose(plane,cells[condition],atol=1e-12,rtol=0.)
        else:cells[condition]=plane
    np.testing.assert_allclose(cells['strain_p0100'],1.01*cells['strain_0000'],atol=1e-12,rtol=0.)
    assert not (CASE/'seeds/strain_p0050').exists()
