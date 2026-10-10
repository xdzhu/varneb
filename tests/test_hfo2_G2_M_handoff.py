import ast
import importlib.util
import json
from pathlib import Path
import pytest
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
from scripts.audit_hfo2_static_replica import audited_results, sha256
from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from vcneb import VCNEB, clamped_plane_vcneb_boundary

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('M_handoff', ROOT/'benchmarks/hfo2_channels/20261008/clamped_M_continue_E062_20261011/queue_continuation.py')
handoff = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(handoff)


def waiter(dependency='afterok:28661019(unfulfilled)', state='PENDING'):
    return ('JobId=28692775 JobName=hfo2-G2-E058-PO_to_T UserId=iai806(16284) '
            f'JobState={state} Partition=hfacnormal01 NumCPUs=32 Dependency={dependency}')


def test_exact_afterok_dependency_handoff():
    handoff.validate_waiter(waiter(), '28692775', {'28661019'})
    handoff.validate_waiter(waiter('afterok:28661019(unfulfilled):999(unfulfilled)'),
                            '28692775', {'28661019','999'})
    handoff.validate_waiter(waiter('afterok:28661019(unfulfilled),afterok:999(unfulfilled)'),
                            '28692775', {'28661019','999'})


@pytest.mark.parametrize('dependency', ['afterany:28661019', 'afterok:28661019?afterok:999',
    'afterok:28661019,afterany:999', '(null)', 'afterok:123'])
def test_changed_or_missing_dependency_stops_before_mutation(dependency):
    with pytest.raises(ValueError):
        handoff.validate_waiter(waiter(dependency), '28692775', {'28661019'})


def test_running_successor_cannot_be_rewired():
    with pytest.raises(ValueError):
        handoff.validate_waiter(waiter(state='RUNNING'), '28692775', {'28661019'})


def test_spare_slot_check_ignores_other_project_but_rejects_duplicate_study():
    body = ('28661019|hfo2-G2-cont|RUNNING|hfacnormal01|32\n'
            '28692775|hfo2-G2-E058|PENDING|hfacnormal01|32\n'
            '800|pc-other|RUNNING|hfacnormal01|64\n')
    assert handoff.validate_active_queue(body, handoff.ALLOWED) == ['28661019']
    with pytest.raises(ValueError):
        handoff.validate_active_queue(body+'999|hfo2-duplicate|RUNNING|hfacnormal01|32\n', handoff.ALLOWED)
    with pytest.raises(ValueError):
        handoff.validate_active_queue(body.replace('28692775|hfo2-G2-E058|PENDING',
                                                   '28692775|hfo2-G2-E058|RUNNING'), handoff.ALLOWED)


def test_single_use_journal_hold_verify_then_release_and_no_retry():
    source = Path(handoff.__file__).read_text()
    assert "open('x')" in source and "'--hold'" in source and 'os.fsync' in source
    assert source.index("record('dependency_update_verified'") < source.index("run(['scontrol','release'")
    assert 'reconcile_do_not_retry' in source
    assert not any(isinstance(node, ast.While) for node in ast.walk(ast.parse(source)))
    assert 'sleep(' not in source and 'scancel' not in source


def test_reconciliation_uses_existing_handle_without_any_submission():
    source = (Path(handoff.__file__).parent/'reconcile_held.py').read_text()
    assert "'28722320'" in source and "'accepted_held_handle'" in source
    assert 'new_submissions_during_reconciliation=0' in source and "open('x')" in source
    assert 'sbatch' not in source and 'submit_once(' not in source and 'scancel' not in source
    assert source.index("record('dependency_update_verified'") < source.index("rules.run(['scontrol','release'")


def test_read_acknowledgement_handles_stale_dependency_without_repeating_update():
    replies = iter([waiter(),waiter('afterok:28661019(unfulfilled),afterok:999(unfulfilled)')])
    queries, pauses = [], []
    def lookup(argv):
        queries.append(argv)
        return next(replies)
    actual = handoff.confirm_waiter('28692775', {'28661019','999'}, lookup=lookup, pause=pauses.append)
    assert handoff.dependency_handles(actual['Dependency']) == {'28661019','999'}
    assert len(queries) == 2 and pauses == [1]
    assert all(q[:3] == ['scontrol','show','job'] for q in queries)
    with pytest.raises(TimeoutError):
        handoff.confirm_waiter('28692775', {'28661019','999'}, lookup=lambda a:waiter(), pause=lambda n:None)
    with pytest.raises(ValueError):
        handoff.confirm_waiter('28692775', {'28661019','999'}, lookup=lambda a:waiter(state='RUNNING'), pause=lambda n:None)


def test_partial_handoff_finish_is_single_use_and_cannot_resubmit():
    source = (Path(handoff.__file__).parent/'finish_handoff.py').read_text()
    assert 'sbatch' not in source and 'submit_once(' not in source and 'scancel' not in source
    assert 'first_dependency_update_repeated=False' in source and "open('x')" in source
    assert source.count("rules.run(['scontrol','update'") == 1


def test_separate_verified_release_has_no_submission_or_dependency_mutation():
    source = (Path(handoff.__file__).parent/'release_verified.py').read_text()
    assert 'sbatch' not in source and "'update'" not in source and 'scancel' not in source
    assert 'execution_started_claimed=False' in source and "open('x')" in source
    assert source.count("rules.run(['scontrol','release'") == 1


@pytest.mark.parametrize('step', [0, 20])
def test_actual_HF_nine_image_raw_EFS_replays_same_boundary_locally(step):
    case = Path(handoff.__file__).parent
    observation = case/'observations'/f'step_{step:04d}'
    report = json.loads((observation/'observation.json').read_text())
    images = read(observation/'evaluated_chain.traj', index=':')
    assert sha256(observation/'evaluated_chain.traj') == report['evaluated_chain_sha256']
    assert len(images) == 9 and report['source_job_id'] == '28661020'
    for i, image in enumerate(images):
        raw = observation/'raw'/f'image_{i:04d}'
        pinned = report['raw_image_evaluations'][i]
        assert pinned['input_sha256'].keys() == {*CONTRACT,'STRU'}
        assert all(pinned['input_sha256'][n] == h for n,h in CONTRACT.items())
        # Only portable INPUT/KPT/STRU/log bytes are physically present here.
        # The six-file/140-call physical check was actually performed on HF.
        for n in ('INPUT','KPT','STRU'):
            assert sha256(raw/n) == pinned['input_sha256'][n]
        assert sha256(raw/'OUT.ABACUS/running_scf.log') == pinned['raw_log_sha256']
        assert same_ordered_geometry(image, read_fixed_hfo2_stru(raw/'STRU'))
        parsed = audited_results(raw)
        for key, expected in (('energy',pinned['energy_eV_cell']), ('forces',pinned['forces_eV_A']),
                              ('stress',pinned['stress_ASE_voigt_eV_A3'])):
            assert np.allclose(parsed[key],expected,atol=1e-12,rtol=0)
        image.calc = SinglePointCalculator(image,**parsed)
    reference = read(case/'seed/substrate.vasp').cell.array
    boundary = clamped_plane_vcneb_boundary(12,reference,allow_tilt=True)
    chain = VCNEB(images,cell_scale=report['cell_scale_A'],k=.2,climb=False,pressure=0,
                  **boundary.vcneb_kwargs(images))
    assert np.linalg.norm(chain.get_forces(),axis=1).max() == pytest.approx(report['replayed_fmax_eV_A'],abs=1e-10,rel=0)
    np.testing.assert_allclose((chain.enthalpies-chain.enthalpies[0])*1000/4,
                               report['relative_enthalpy_meV_fu'],atol=1e-9,rtol=0)
    assert not report['ordinary_residual_pass'] and report['new_DFT_calls'] == 0


def test_actual_HF_dispatch_and_audit_receipts_pin_one_unchanged_chain():
    case = Path(handoff.__file__).parent
    audit = json.loads((case/'audit_receipt.json').read_text())
    receipt = json.loads((case/'submission_receipt.json').read_text())
    original = [json.loads(r) for r in (case/'submission_journal.jsonl').read_text().splitlines()]
    assert [r['job_id'] for r in original if r['event']=='accepted_held_handle'] == ['28722320']
    assert receipt['job_id'] == '28722320' and receipt['total_E062_submissions'] == 1
    assert receipt['new_submissions_during_reconciliation'] == 0
    assert receipt['waiting_jobs_verified'] == ['28692775','28692776']
    assert receipt['audit_receipt_sha256'] == sha256(case/'audit_receipt.json')
    assert audit['fresh_parent_interior_SCFs'] == 140
    assert audit['current_frame_exact_caches_checked'] == 9
    assert not receipt['physical_inputs_changed'] and not receipt['holdout_generated_or_read']
    assert sha256(case/'executed_queue_continuation.py') == 'fd6ee53777bc8b13b4fad37e6b572ad01824de916295ad9aef797d8d4f30ae5b'
    assert receipt['executed_release_sha256'] == sha256(case/'release_verified.py')
    assert receipt['reviewed_handoff_rules_sha256'] == sha256(case/'queue_continuation.py')
    assert not receipt['execution_started_claimed']  # independent accounting is separate evidence
