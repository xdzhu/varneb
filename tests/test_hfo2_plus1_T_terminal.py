"""Reparse nine real portable +1% terminal frames, not all42 earlier SCFs."""
import importlib.util
import json
from pathlib import Path, PurePosixPath

import numpy as np
import pytest
from ase.io import read

from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.analyze_hfo2_clamped_residual import analyze
from scripts.audit_hfo2_static_replica import audited_results, sha256

REPO = Path(__file__).resolve().parents[1]
BASE = REPO/'benchmarks/hfo2_channels/20261008'
CASE = BASE/'clamped_plus1_T_terminal_E087_20261011'
OBS = CASE/'observations/step_0005'
AUDIT_SHA = '32a6e49a618e6e62632ae71a5e0ed093095254f6858221be2ade3d83cc194a28'


@pytest.fixture(scope='module')
def evidence():
    assert sha256(CASE/'audit_receipt.json') == AUDIT_SHA
    return (json.loads((CASE/'audit_receipt.json').read_text()),
            json.loads((OBS/'observation.json').read_text()),
            read(OBS/'evaluated_chain.traj', index=':'))


def test_actual_terminal_receipt_has_registered_native_scope(evidence):
    audit, obs, images = evidence
    assert audit['job_id'] == obs['source_job_id'] == '28709788'
    assert audit['ordinary_converged'] and obs['ordinary_residual_pass']
    assert audit['terminal_step'] == obs['snapshot_step'] == 5
    assert audit['fresh_interior_SCFs'] == len(audit['calls']) == 42
    assert audit['initial_fresh_interior_SCFs'] == 7
    assert len(audit['runtime_code_sha256']) == 338
    assert audit['scheduler_row'] == '28709788|COMPLETED|0:0|32|hfacnormal01|iai806'
    assert audit['new_DFT_calls'] == audit['scheduler_mutations'] == 0
    assert not audit['holdout_generated_or_read']
    assert not audit['continuous_TS_or_barrier_error_bound_certified']
    boundary = obs['mechanical_boundary']
    assert boundary['reference_file_sha256'] == 'f1b88a625b851149201149faa91f5d4d5565b1e77022a0419974272d6c9cfa74'
    assert boundary['fixed_cell_rows'] == [0, 1] and boundary['allow_tilt']
    assert obs['formula_units'] == 4 and obs['pressure_GPa'] == 0 and not obs['climb']
    assert len(images) == 9 and sha256(OBS/'evaluated_chain.traj') == obs['evaluated_chain_sha256']
    for call in audit['calls']:
        assert all(call['input_sha256'][name] == h for name, h in CONTRACT.items())


@pytest.mark.parametrize('index', range(9))
def test_portable_raw_frame_matches_native_record_and_geometry(evidence, index):
    audit, obs, images = evidence
    raw = OBS/f'raw/image_{index:04d}'
    record, image = obs['raw_image_evaluations'][index], images[index]
    assert record['image_index'] == index and record['ordered_periodic_geometry_matched']
    assert all(record['input_sha256'][name] == h for name, h in CONTRACT.items())
    for name in ('INPUT', 'KPT', 'STRU'):
        assert sha256(raw/name) == record['input_sha256'][name]
    assert sha256(raw/'source_audit.json') == record['audit_sha256']
    assert sha256(raw/'OUT.ABACUS/running_scf.log') == record['raw_log_sha256']
    assert sha256(OBS/f'POSCAR_{index:02d}') == record['snapshot_POSCAR_sha256']
    assert same_ordered_geometry(image, read_fixed_hfo2_stru(raw/'STRU'))
    parsed = audited_results(raw)  # one native DSIZE32 and electronic convergence
    for key in ('energy', 'forces', 'stress'):
        np.testing.assert_allclose(parsed[key], image.calc.results[key], atol=1e-12, rtol=0)


def test_complete_ordinary_force_replay_matches_native_components(evidence):
    result = analyze(OBS)
    native = json.loads((OBS/'residual.json').read_text())
    assert result['source_job_id'] == '28709788' and result['ordinary_residual_pass']
    assert result['fmax_eV_A'] == pytest.approx(.0949479822061884, abs=1e-12)
    assert result['limiting_image_index'] == 7 and result['limiting_vector_block'] == 'atoms'
    assert result['new_DFT_calls'] == 0
    for a, b in zip(result['images'], native['images']):
        assert a.keys() == b.keys()
        for key in a:
            if isinstance(a[key], (int, float, list)):
                np.testing.assert_allclose(a[key], b[key], atol=1e-10, rtol=0)
            else:
                assert a[key] == b[key]


def test_signed_peak_identity_and_seen_training_contrast(evidence):
    audit, obs, images = evidence
    energies = np.array([image.get_potential_energy() for image in images])
    relative = (energies-energies[0])*1000/4
    np.testing.assert_allclose(relative, obs['relative_enthalpy_meV_fu'], atol=1e-9, rtol=0)
    assert np.argmax(relative) == 4
    forward, reverse = obs['sampled_forward_barrier_meV_fu'], obs['sampled_reverse_barrier_meV_fu']
    assert forward == pytest.approx(96.76365039467782, abs=1e-9)
    assert reverse == pytest.approx(82.61776706149249, abs=1e-9)
    assert forward-reverse == pytest.approx(relative[-1], abs=1e-9)
    old = json.loads((BASE/'clamped_T_terminal_E079_20261011/observations/step_0010/observation.json').read_text())
    assert old['formula_units'] == obs['formula_units'] and old['pressure_GPa'] == 0
    assert old['source_job_id'] == '28692775' and old['ordinary_residual_pass']
    for condition in (old, obs):
        boundary = condition['mechanical_boundary']
        assert boundary['fixed_cell_rows'] == [0, 1] and boundary['allow_tilt']
        assert boundary['cell_dofs'] == 3
        for frame in condition['raw_image_evaluations']:
            assert all(frame['input_sha256'][name] == h for name, h in CONTRACT.items())
    np.testing.assert_allclose(np.array(obs['mechanical_boundary']['reference_cell_A'])[:2],
                               1.01*np.array(old['mechanical_boundary']['reference_cell_A'])[:2],
                               atol=1e-10, rtol=0)
    assert old['mechanical_boundary']['reference_file_sha256'] != obs['mechanical_boundary']['reference_file_sha256']
    peak_change = forward-old['sampled_forward_barrier_meV_fu']
    well_change = relative[-1]-old['relative_enthalpy_meV_fu'][-1]
    assert peak_change == pytest.approx(12.347026513452874, abs=1e-9)
    assert well_change == pytest.approx(1.7859082190625486, abs=1e-9)
    # Numeric serialization checks are not a DFT/continuous-sampling error bound.


def test_reviewed_parent_and_accepted_once_only_M_transaction(evidence):
    path = CASE/'executed_source/release_existing_M_once.py'
    spec = importlib.util.spec_from_file_location('E086_parent_check_only', path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    # The actual helper executed on Linux. Replay its POSIX path declaration
    # on Windows without changing the native receipt or executed source bytes.
    helper.WAIT = PurePosixPath(helper.WAIT.as_posix())
    helper.verified_parent(evidence[0])  # no native dispatch function is invoked
    journal = CASE/'admission_journal.jsonl'
    receipt = json.loads((CASE/'admission_receipt.json').read_text())
    rows = [json.loads(line) for line in journal.read_text().splitlines()]
    assert [r['event'] for r in rows] == ['admission_verified', 'release_intent', 'release_accepted']
    assert rows[0]['parent_audit_sha256'] == AUDIT_SHA
    assert rows[1]['job_id'] == rows[2]['job_id'] == receipt['target_job'] == '28709789'
    assert receipt['journal_sha256'] == sha256(journal)
    assert receipt['executed_helper_sha256'] == sha256(path)
    assert receipt['release_writes'] == 1
    assert receipt['new_submissions'] == receipt['dependency_writes'] == receipt['new_independent_chains'] == 0
    assert not receipt['startup_claimed'] and not receipt['holdout_generated_or_read']
    startup = json.loads((CASE/'startup_observation.json').read_text())
    assert startup['actual_first_interior_header_DSIZE'] == startup['allocated_CPUs'] == 32
    assert startup['active_study_allocations'] == 2
    assert not startup['first_electronic_or_step_zero_convergence_claimed']
