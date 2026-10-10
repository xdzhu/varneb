"""Actual training endpoints, no holdout path labels or extra calculations."""
from pathlib import Path
from copy import deepcopy
import hashlib
from io import BytesIO
import json
import math
from numbers import Real
import zipfile
import pytest
from scripts.analyze_hfo2_endpoint_strain_work import analyse


CASE=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008'


@pytest.fixture(scope='module')
def report():return analyse(CASE)


def assert_historical_sources(bundle, receipt, text_checks, source_proof):
    """Exact executed bytes, not mutable paper/library working-copy bytes."""
    assert hashlib.sha256(bundle).hexdigest() == source_proof['source_bundle_sha256']
    assert source_proof['exact_raw_source_members'] == receipt['source_sha256']
    with zipfile.ZipFile(BytesIO(bundle)) as archive:
        assert len(archive.namelist()) == 8 and set(archive.namelist()) == set(receipt['source_sha256'])
        for name,h in receipt['source_sha256'].items():
            raw = archive.read(name)
            assert hashlib.sha256(raw).hexdigest() == h
            c = text_checks['checks'][name]
            assert c['raw_executed_sha256'] == h
            assert hashlib.sha256(raw.replace(b'\r\n',b'\n')).hexdigest() == c['LF_canonical_text_sha256']


def numeric_fields(value, prefix=()):
    if isinstance(value,dict):
        return {p:n for k,v in value.items() for p,n in numeric_fields(v,prefix+(k,)).items()}
    if isinstance(value,(list,tuple)):
        return {p:n for i,v in enumerate(value) for p,n in numeric_fields(v,prefix+(i,)).items()}
    if isinstance(value,Real) and not isinstance(value,bool):
        assert math.isfinite(value)
        return {prefix:float(value)}
    return {}


def assert_numeric_replay(local, actual):
    left,right = numeric_fields(local),numeric_fields(actual)
    assert left.keys() == right.keys()
    assert all(abs(left[p]-right[p]) < 1e-9 for p in left)
    return len(left)


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
    source_delivery=CASE/'replay_delivery_E061_20261011'
    source_proof=json.loads((source_delivery/'source_export_receipt.json').read_text())
    assert receipt['status']=='original_HF_environment_and_ten_raw_endpoints_replay_passed'
    assert receipt['source_files_byte_checked']==8
    assert text_checks['actual_source_files']==8 and text_checks['raw_source_still_matches_original_replay']
    assert not text_checks['source_files_modified'] and not text_checks['new_DFT_calls']
    assert_historical_sources((source_delivery/'source_bundle.zip').read_bytes(),receipt,text_checks,source_proof)
    assert source_proof['historical_E059_receipt_sha256'] == hashlib.sha256((delivery/'replay_receipt.json').read_bytes()).hexdigest()
    assert not source_proof['new_numerical_HF_replay'] and not source_proof['new_DFT_calls']
    assert not source_proof['DFT_inputs_or_holdout_read'] and not source_proof['physical_parameters_or_running_sources_or_jobs_changed']
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
    assert_numeric_replay(report,actual)
    for local,remote in zip(report['B1_training_response_implications'],actual['B1_training_response_implications']):
        assert local['channel']==remote['channel']
        assert local['B1_final_following_bottleneck_shift_meV_fu']==pytest.approx(remote['B1_final_following_bottleneck_shift_meV_fu'],abs=1e-9)


def test_editable_manuscript_changes_do_not_invalidate_historical_source_evidence(monkeypatch):
    delivery = CASE/'endpoint_strain_work_E059_20261010'
    sources = CASE/'replay_delivery_E061_20261011'
    receipt = json.loads((delivery/'replay_receipt.json').read_text())
    checks = json.loads((delivery/'source_text_checks.json').read_text())
    proof = json.loads((sources/'source_export_receipt.json').read_text())
    read_bytes = Path.read_bytes
    root = CASE.parents[2]
    drafts = {root/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md',root/'paper/VARNEB_JCTC/METHODS_DRAFT.md'}
    def revised(self):
        body = read_bytes(self)
        return body+b'\nNoncomputational drafting changes.\n' if self in drafts else body
    monkeypatch.setattr(Path,'read_bytes',revised)
    # The old working-copy contract would fail; no actual paper file is edited.
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() != receipt['source_sha256'][p.relative_to(root).as_posix()] for p in drafts)
    assert_historical_sources((sources/'source_bundle.zip').read_bytes(),receipt,checks,proof)


@pytest.mark.parametrize('field',['partial','atomic','released_cell','chord','well','B1'])
def test_full_material_replay_rejects_changed_work_or_well_terms(report,field):
    changed = deepcopy(report)
    if field == 'partial':changed['endpoints']['PO_plus'][0]['imposed_plane_partial_eV_fu_per_strain'] += .001
    elif field == 'atomic':changed['pairs'][0]['configuration_chord_endpoint_terms'][0]['atomic_eV_fu_per_strain'] += .001
    elif field == 'released_cell':changed['pairs'][0]['configuration_chord_endpoint_terms'][0]['released_cell_eV_fu_per_strain'] += .001
    elif field == 'chord':changed['pairs'][0]['configuration_chord_defect_meV_fu'] += .001
    elif field == 'well':changed['pairs'][0]['energy_shift_meV_fu'] += .001
    else:changed['B1_training_response_implications'][0]['B1_final_following_bottleneck_shift_meV_fu'] += .001
    with pytest.raises(AssertionError):assert_numeric_replay(changed,report)


def test_forged_bundle_digest_does_not_hide_changed_historical_source():
    delivery = CASE/'endpoint_strain_work_E059_20261010'
    sources = CASE/'replay_delivery_E061_20261011'
    receipt = json.loads((delivery/'replay_receipt.json').read_text())
    checks = json.loads((delivery/'source_text_checks.json').read_text())
    proof = json.loads((sources/'source_export_receipt.json').read_text())
    forged = BytesIO()
    with zipfile.ZipFile(sources/'source_bundle.zip') as old, zipfile.ZipFile(forged,'w') as out:
        for name in old.namelist():
            body = old.read(name)
            if name == 'vcneb/strain_work.py':body += b'\n# changed source\n'
            out.writestr(name,body)
    raw = forged.getvalue();proof['source_bundle_sha256'] = hashlib.sha256(raw).hexdigest()
    with pytest.raises(AssertionError):assert_historical_sources(raw,receipt,checks,proof)
