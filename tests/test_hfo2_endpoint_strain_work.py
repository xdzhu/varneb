"""Actual training endpoints, no holdout path labels or extra calculations."""
from pathlib import Path
import hashlib
import json
import pytest
from scripts.analyze_hfo2_endpoint_strain_work import analyse


CASE=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008'


@pytest.fixture(scope='module')
def report():return analyse(CASE)


def test_actual_pairs_are_same_boundary_original_contract_and_not_forecasts(report):
    assert report['status']=='seen_training_endpoint_response_NOT_barrier_forecast'
    assert report['training_strains']==[0.,.01] and len(report['endpoints'])==5
    assert report['formula_units']==4 and report['pressure_eV_A3']==0.
    assert not report['new_DFT_calls'] and not report['holdout_generated_or_read']
    assert not report['physical_inputs_changed'] and not report['forecast_batch_frozen']
    assert not report['full_relaxed_branch_derivative_or_Hessian_stability_certified']
    assert not report['gradient_consistency_or_barrier_error_bound_established']
    for endpoints in report['endpoints'].values():
        for e in endpoints:
            assert e['max_atomic_force_eV_A']<.03 and e['open_traction_kbar']<2
            assert not e['full_pseudo_or_basis_bytes_rechecked_offline']
            assert e['source_audit']['raw_full_physical_bytes_checked_on_HF']
    assert all(p['unmeasured_interior_work_or_quadrature_error_bound'] is None for p in report['pairs'])


def test_B1_training_nulls_keep_initial_shift_and_distinguish_final_wells(report):
    shifts={p['phase_representation']:p['energy_shift_meV_fu'] for p in report['pairs']}
    assert shifts['PO_plus']==pytest.approx(6.8065023770031985,abs=1e-8)
    expected={'PO_to_T':1.7859082190625486,'PO_to_M':-22.206478976841027,
              'PO_flip_T_pattern_preserving':.03525880219967803,
              'PO_flip_T_pattern_reversing':.03525880129018333}
    for row in report['B1_training_response_implications']:
        assert row['B1_fixed_bottleneck_shift_meV_fu']==pytest.approx(-shifts['PO_plus'])
        assert row['B1_final_following_bottleneck_shift_meV_fu']==pytest.approx(expected[row['channel']],abs=1e-8)
        assert 'absolute barrier requires' in row['interpretation']


def test_delivered_HF_replay_checks_original_ten_terminal_assets_without_new_DFT(report):
    delivery=CASE/'endpoint_strain_work_E059_20261010'
    receipt=json.loads((delivery/'replay_receipt.json').read_text())
    actual=json.loads((delivery/'endpoint_work_HF.json').read_text())
    text_checks=json.loads((delivery/'source_text_checks.json').read_text())
    root=CASE.parents[2]
    assert receipt['status']=='original_HF_environment_and_ten_raw_endpoints_replay_passed'
    assert receipt['source_files_byte_checked']==8
    exact={'vcneb/strain_work.py','scripts/analyze_hfo2_endpoint_strain_work.py',
           'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md','paper/VARNEB_JCTC/METHODS_DRAFT.md'}
    assert text_checks['actual_source_files']==8 and text_checks['raw_source_still_matches_original_replay']
    assert not text_checks['source_files_modified'] and not text_checks['new_DFT_calls']
    for name,h in receipt['source_sha256'].items():
        raw=(root/name).read_bytes()
        c=text_checks['checks'][name]
        assert c['raw_executed_sha256']==h
        if name in exact:
            assert hashlib.sha256(raw).hexdigest()==h
        else:
            # Compare BOTH legacy representations to the explicitly exported
            # canonical text. The exact executed raw hash remains separate;
            # Git archives and Windows checkouts need not share their newlines.
            raw=raw.replace(b'\r\n',b'\n')
            assert hashlib.sha256(raw).hexdigest()==c['LF_canonical_text_sha256']
    assert receipt['maximum_numeric_replay_difference']<1e-9
    assert receipt['numeric_tolerance_is_NOT_DFT_error_bound']
    assert not receipt['new_DFT_calls'] and receipt['external_executable_launch_prohibited']
    assert not receipt['environment_installed_or_upgraded'] and not receipt['HF_pytest_run']
    assert not receipt['physical_parameters_changed'] and not receipt['holdout_generated_or_read']
    assert not receipt['production_source_or_job_changed']
    assert len(receipt['actual_terminal_checks'])==10
    assert all(r['original_six_physical_bytes_and_STRU_log_checked'] and r['native_32MPI_raw_EFS_checked']
               for r in receipt['actual_terminal_checks'])
    assert actual['status']==report['status'] and not actual['forecast_batch_frozen']
    for local,remote in zip(report['B1_training_response_implications'],actual['B1_training_response_implications']):
        assert local['channel']==remote['channel']
        assert local['B1_final_following_bottleneck_shift_meV_fu']==pytest.approx(remote['B1_final_following_bottleneck_shift_meV_fu'],abs=1e-9)
