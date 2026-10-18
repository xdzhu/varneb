"""Six independent final-candidate native replays, not a TS certificate."""
import importlib.util
import json
from pathlib import Path
import shutil

import pytest

BASE = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008'
spec = importlib.util.spec_from_file_location('complete_T_controls', BASE/'matched_pattern_seeded_E271_20261018/analyze_complete.py')
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


def test_all_six_candidates_replayed_and_all_costs_preserved():
    r = module.analyze()
    assert len(r['native_candidate_replays']) == 6
    assert sum(x['local_native_SCFs_reparsed'] for x in r['native_candidate_replays']) == 6
    assert r['costs']['total_new_SCF_calls'] == 55
    assert r['costs']['allocated_core_hours'] == pytest.approx(52.3733333333)
    assert r['costs']['all_earlier_curvature_and_path_preparation_costs_additional']
    assert [x['dE_J_minus_F_meV_fu'] for x in r['rows']] == pytest.approx([-.17737459666,-.08609264614,-.24002449709])
    assert [x['candidate_source_batch'] for x in r['native_candidate_replays']] == ['E271','E271','E270','E271','E271','E271']
    assert sum(x['new_repair_attempts'] for x in r['native_candidate_replays']) == 14
    for x in r['native_candidate_replays']:
        assert x['released_gradient_norm_eV_A'] <= .03
        assert x['physical_retained_defect_A'] < 1e-15
        assert not x['energy_partition_or_stationary_TS_or_minimum_certified']
        if x['arm'] == 'A': assert x['cell_coordinate_norm_A'] == 0
    assert r['old_labels_and_scores_unchanged'] and r['new_DFT_calls'] == 0
    assert not r['full_stationarity_stability_or_corrected_barrier_certified']


def test_corrupt_archive_is_rejected_before_extraction(tmp_path):
    source = BASE/module.ARCHIVES['E271'][0]; case = tmp_path/'case'; case.mkdir()
    shutil.copyfile(source/'native_terminal_receipt.json',case/'native_terminal_receipt.json')
    (case/'native_terminal.tar.gz').write_bytes(b'not-the-committed-archive')
    destination = tmp_path/'out'
    with pytest.raises(ValueError,match='archive changed'):
        module.restore(case,module.ARCHIVES['E271'][1],destination)
    assert not destination.exists()
