"""The measured atomic treatment is not a joint-stability or prediction claim."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import re

import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'benchmarks/hfo2_channels/20261008/manuscript_atomic_treatment_E195_20261015'
PAPER = ROOT / 'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'
spec = importlib.util.spec_from_file_location('manuscript_actual_atomic', CASE / 'check_tables.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


def test_four_new_rows_preserve_both_centers_and_steps():
    data, native = a.source()
    assert len(a.require_rows(PAPER.read_text(encoding='utf-8'), data, native)) == 4
    assert [(r['channel'], r['probe_step_A']) for r in data['rows']] == [
        ('T', .01), ('T', .02), ('preserving', .01), ('preserving', .02)]


@pytest.mark.parametrize('index', range(4))
def test_changed_missing_or_duplicate_complete_treatment_row_refused(index):
    data, native = a.source()
    text = PAPER.read_text(encoding='utf-8')
    row = a.expected_rows(data, native)[index]
    for broken in (text.replace(row, ''), text.replace(row, row+'\n'+row),
                   text.replace(row, row.replace('|', '| wrong ', 1))):
        with pytest.raises(ValueError, match='source-derived'):
            a.require_rows(broken, data, native)


@pytest.mark.parametrize('fault', ['prefix', 'failed', 'missing', 'unequal_cohort'])
def test_incomplete_or_unequal_data_cannot_support_complete_table(fault):
    data, native = (deepcopy(x) for x in a.source())
    if fault == 'prefix': native['audit_mode'] = 'prefix'
    elif fault == 'failed': native['sacct'] = native['sacct'].replace('COMPLETED', 'FAILED')
    elif fault == 'missing': native['gradient_candidates'] = 9
    else: data['rows'][0]['non_origin_coverage_statistics']['model_minus_observed']['atomic']['shape_eV_cell']['ids'] = [1]
    with pytest.raises(ValueError):
        a.expected_rows(data, native)


def test_prose_is_local_signed_and_not_unique_partition_or_joint_stability():
    text = ' '.join(PAPER.read_text(encoding='utf-8').split())
    for phrase in ('Across the five registered coordinates', '44.06--48.35 meV/cell',
                   '0.80--1.14 meV/cell lower', 'not the whole switching path',
                   'not a unique strain-energy term', 'neither joint stationarity nor transverse stability',
                   'Both original quadratic models were fixed before the atomic labels',
                   'independent material predictions are still missing',
                   '39.749558 core-hours', '40.684444 core-hours', 'not a speedup measurement'):
        assert phrase in text
    ledger = (PAPER.parent / 'TERMINOLOGY_LEDGER.md').read_text(encoding='utf-8')
    assert '| Atomic-only treatment |' in ledger and '| Additional joint-release difference |' in ledger


def test_current_source_checker_does_not_claim_full_goal():
    result = a.check()
    assert result['new_quantitative_rows'] == 4
    assert result['new_DFT_or_scheduler_writes'] == 0
    assert result['structural_reorganization_applied']
    assert not result['stable_surface_TS_independent_accuracy_or_JCTC_readiness_established']


def test_approved_four_part_results_and_six_appendix_diagnostics():
    text = PAPER.read_text(encoding='utf-8')
    main = text.split('## 3.', 1)[1].split('## 4.', 1)[0]
    assert re.findall(r'^### (3\.\d+)\.', main, re.M) == ['3.1', '3.2', '3.3', '3.4']
    assert 'Across the five registered coordinates' in ' '.join(main.split())
    assert 'An independent prediction advantage has not yet been established' in main
    assert 'Eight prospective local test coordinates per center' in main
    assert 'They remain uncomputed' in main
    assert 'failed\nthe necessary post-release sign condition' in main
    assert 'original refusal' in main
    # Six historical diagnostics remain; the actual four-point result adds E.7.
    assert re.findall(r'^### (E\.\d+)\.', text, re.M) == [f'E.{i}' for i in range(1, 8)]
    assert not re.search(r'Section 3\.[5-9]\b|Table [4-9]\b', text)
    for number in range(11, 18):
        assert f'Table A{number}' in text
