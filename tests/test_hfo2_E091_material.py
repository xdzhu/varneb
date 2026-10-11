"""Reparse 18 portable native frames, not all 357 historical HF SCFs locally."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.analyze_hfo2_clamped_residual import analyze
from scripts.audit_hfo2_static_replica import audited_results, sha256

CASE = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/clamped_rolling_native_E091_20261011'
PARENTS = {
    '28709789': ('terminal',10,77,.5701098697954412,'28727368',
        '1a05735749d8893b2c6f992f0ab90f347da9cd65b674a8c1d920227dd8565178'),
    '28723655': ('terminal-r2',40,280,.14343717607417714,'28727382',
        '61edc4548d85984e78bbafc29d8b9d366157cfef48ec8175cabf3adfc7ce7e93'),
}


@pytest.mark.parametrize('parent',list(PARENTS))
def test_actual_terminal_latest_force_raw_frames_seed_and_once_only_journal(parent):
    subdir,step,count,force,job,digest=PARENTS[parent]
    root=CASE/parent
    terminal=root/subdir
    assert sha256(terminal/'audit_receipt.json')==digest
    receipt=json.loads((terminal/'audit_receipt.json').read_text())
    assert receipt['ordinary_converged'] is False
    assert receipt['fresh_interior_SCFs']==count and receipt['current_frame_exact_caches_checked']==9
    assert receipt['all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked'] is True
    observation=terminal/'observations'/f'step_{step:04d}'
    r=json.loads((observation/'observation.json').read_text())
    for image in r['raw_image_evaluations']:
        directory=observation/'raw'/f"image_{image['image_index']:04d}"
        for name in ('INPUT','KPT','STRU'):
            assert sha256(directory/name)==image['input_sha256'][name]
        assert sha256(directory/'OUT.ABACUS/running_scf.log')==image['raw_log_sha256']
        raw=audited_results(directory)
        for key,field in (('energy','energy_eV_cell'),('forces','forces_eV_A'),('stress','stress_ASE_voigt_eV_A3')):
            np.testing.assert_allclose(raw[key],image[field],atol=1e-12,rtol=0)
    residual=analyze(observation)
    assert residual['fmax_eV_A']==pytest.approx(force,abs=1e-12)
    assert not residual['ordinary_residual_pass'] and residual['new_DFT_calls']==0
    assert sha256(terminal/'seed/seed.traj')==r['evaluated_chain_sha256']
    seed=json.loads((terminal/'seed/manifest.json').read_text())
    assert seed['source_step']==step and seed['optimizer_state_restored'] is False
    for name,h in seed['files_sha256'].items():assert sha256(terminal/'seed'/name)==h
    dispatch=json.loads((root/'dispatch_receipt.json').read_text())
    assert dispatch['job_id']==job and dispatch['audit_receipt_sha256']==digest
    assert dispatch['new_independent_chains']==dispatch['dependency_writes']==0
    assert dispatch['total_submissions']==dispatch['release_writes']==1
    journal=root/'dispatch_journal.jsonl'
    assert sha256(journal)==dispatch['journal_sha256']
    records=[json.loads(line) for line in journal.read_text().splitlines()]
    assert [e['event'] for e in records]==['one_held_submission_intent','one_accepted_handle',
        'one_release_intent','one_existing_handle_released']
    assert all(e['job_id']==job for e in records[1:])


def test_actual_old_prepare_rejects_registered_continuation_before_writes(tmp_path):
    source=CASE/'original_executed_prepare_v1.py'
    assert sha256(source)=='2dc91260be19a6bc9ee44fcc39d319199bbf634dd4f1d90396427b735d72b7b8'
    spec=importlib.util.spec_from_file_location('original_E054_prepare_snapshot',source)
    legacy=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    observation=CASE/'28723655/terminal/observations/step_0040'
    with pytest.raises(ValueError,match='audited incomplete ordinary clamped G2 pilot'):
        legacy.prepare(observation,tmp_path/'seed')
    assert not (tmp_path/'seed').exists()
