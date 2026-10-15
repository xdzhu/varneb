"""Complete native proof and all-hole accounting; no new material calculations."""
from copy import deepcopy
from collections import Counter
import importlib.util
import json
from pathlib import Path
import re

import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_preserving_terminal_E179_20261015'
spec = importlib.util.spec_from_file_location('complete_preserving_holes', CASE / 'replay_holes.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


@pytest.fixture
def proof():
    path = CASE / 'preserving_complete_native_audit.json'
    assert a.sha(path) == a.PROOF_SHA
    return json.loads(path.read_bytes())


def test_actual_complete_proof_has_all_eight_holes_not_a_complete_converged_surface(proof):
    assert a.validate_proof(proof) == {13, 14, 19, 20, 21, 22, 23, 24}
    assert proof['gradient_candidates'] == 17 and proof['new_SCFs_reparsed'] == 172
    assert all(r['replay']['original_order_cache_only_replay_passed'] for r in proof['candidate_rows'])
    assert proof['completed_call_proxy_core_hours'] == pytest.approx(193.85907278201645, abs=1e-12, rel=0)
    assert proof['allocation_core_hours_at_observation'] == 195.12
    assert not proof['released_stability_branch_uniqueness_contour_TS_or_forecast_certified']
    holes = [r for r in proof['candidate_rows'] if r['candidate'] is None]
    assert {r['id']: r['completed_SCFs'] for r in holes} == {
        13: 12, 14: 12, 19: 12, 20: 12, 21: 12, 22: 12, 23: 10, 24: 12}
    assert sum(r['completed_SCFs'] for r in holes) == 94


@pytest.mark.parametrize('fault', ['prefix', 'job', 'running', 'CPUs', 'pending', 'missing_point',
    'candidate_count', 'SCF_count', 'stability', 'failure', 'runtime'])
def test_false_completion_or_incomplete_source_counts_rejected(proof, fault):
    p = deepcopy(proof)
    if fault == 'prefix': p['audit_mode'] = 'prefix'
    elif fault == 'job': p['job_id'] = 'different'
    elif fault == 'running': p['sacct'] = p['sacct'].replace('COMPLETED', 'RUNNING')
    elif fault == 'CPUs': p['sacct'] = p['sacct'].replace('|32|', '|1|')
    elif fault == 'pending': p['pending_node'] = dict(point_id=24)
    elif fault == 'missing_point': p['candidate_rows'].pop()
    elif fault == 'candidate_count': p['gradient_candidates'] += 1
    elif fault == 'SCF_count': p['new_SCFs_reparsed'] += 1
    elif fault == 'stability': p['released_stability_branch_uniqueness_contour_TS_or_forecast_certified'] = True
    elif fault == 'failure': p['unfinished_node_failure_paths'] = ['example_failure']
    else: p['runtime_manifest_sha256'] = '0' * 64
    with pytest.raises(ValueError): a.validate_proof(p)


def test_existing_output_refuses_before_proof_or_native_access(tmp_path, monkeypatch):
    path = tmp_path / 'old.json'
    path.write_text('preserve')
    monkeypatch.setattr(a.sys, 'argv', ['holes', '--proof', 'unused', '--output', str(path)])
    monkeypatch.setattr(a, 'observe', lambda *args: pytest.fail('must refuse before access'))
    with pytest.raises(FileExistsError): a.main()
    assert path.read_text() == 'preserve'


def test_changed_proof_bytes_refuse_before_loading_original_runtime(tmp_path, monkeypatch):
    path = tmp_path / 'fake.json'
    path.write_text('{}')
    for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        monkeypatch.setenv(key, '1')
    monkeypatch.setattr(a, 'module', lambda *args: pytest.fail('must refuse before imports'))
    with pytest.raises(ValueError, match='exact complete native proof'): a.observe(path)


def test_actual_all_eight_guard_replays_preserve_holes_and_original_decisions():
    path = CASE / 'completed_holes.json'
    assert a.sha(path) == 'd7c89c5125e398f08aa856dd222da759770ecb04ea7c1defcaa94dacd0d094ab'
    result = json.loads(path.read_bytes())
    assert [h['point_id'] for h in result['holes']] == [13, 14, 19, 20, 21, 22, 23, 24]
    assert result['raw_SCFs_reparsed_by_original_complete_audit'] == 172
    assert result['new_raw_SCF_reparse_this_diagnostic'] == result['new_DFT_or_scheduler_writes'] == 0
    assert not result['optimizer_guard_force_target_or_physics_changed']
    assert not result['conditional_stability_boundary_KKT_unreachable_minimum_TS_or_prediction_certified']
    for h in result['holes']:
        assert h['replay'] == dict(original_order_cache_only_replay_passed=True, new_DFT=0)
        assert h['reported_final_candidate'] is None and h['last_trial_not_promoted']
        assert not h['optimizer_or_guard_decisions_changed']
        assert dict(Counter(e['reason'] for e in h['guard_events'])) == h['guard_counts']
        assert set(h['guard_counts']) == {'accepted', 'original_internal_radius'}
        rejected = [e for e in h['guard_events'] if not e['safe']]
        assert all(e['minimum_distance_A'] is None for e in rejected)
        assert not any(t['below_gradient_target'] for t in h['completed_trial_rows'])
        assert min(t['released_gradient_norm_eV_A'] for t in h['completed_trial_rows']) > .03
    assert sum(h['guard_counts']['original_internal_radius'] for h in result['holes']) == 436


def test_paper_retains_all_previous_tables_figures_citations_and_wrapper():
    baseline = json.loads((CASE / 'preservation.json').read_bytes())
    paper = ROOT / 'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'
    text = paper.read_text(encoding='utf-8')
    # Approved main/appendix movement changes order, not content or multiplicity.
    current = Counter(line for line in text.splitlines() if line.startswith('|'))
    assert not Counter(baseline['table_lines'])-current, 'previous table row removed or changed'
    assert len(baseline['table_lines']) == 119
    assert re.findall(r'!\[[^\]]*\]\(([^)]+)\)', text) == baseline['figure_targets']
    assert re.findall(r'\]\((https?://[^)]+)\)', text) == baseline['external_links']
    assert a.sha(ROOT / 'paper/VARNEB_JCTC/varneb_JCTC_draft.tex') == baseline['wrapper_sha256']
    for name, digest in baseline['figure_sha256'].items():
        assert a.sha(paper.parent / name) == digest
    assert '17, including the' in text and 'They consumed 94 of the 172 SCFs' in text
    assert 'not nonexistence of a solution' in text
