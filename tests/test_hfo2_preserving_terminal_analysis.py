"""Original full signed native batch: local successes and bounded failures."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/preserving_pattern_branches_E272_20261018'
spec=importlib.util.spec_from_file_location('E272_whole_native_test',CASE/'analyze_terminal.py')
reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)


@pytest.fixture(scope='module')
def report():
    return reader.analyze()


def test_all_54_native_attempts_and_costs_retained(report):
    assert report['local_native_SCFs_reparsed']==report['costs']['new_SCF_calls']==54
    assert sum(r['total_new_attempts'] for r in report['rows'])==54
    assert report['costs']['allocated_core_hours']==pytest.approx(46.026666666666664)
    assert report['costs']['native_proxy_core_hours']==pytest.approx(45.64812823663155)
    assert report['costs']['all_earlier_costs_additional']
    assert report['new_DFT_calls']==0 and report['old_labels_and_scores_unchanged']
    assert not report['stationary_TS_or_full_stability_or_corrected_barrier_certified']
    assert not report['full_orbital_pseudopotential_bytes_rechecked_locally']


def test_zero_signs_have_opposite_transverse_structural_proxies_not_new_barriers(report):
    rows=report['rows'][:2]
    assert [r['sign'] for r in rows]==['minus','plus'] and all(r['strain_percent']==0 for r in rows)
    energies=[]
    for row,sign in zip(rows,(1,-1)):
        assert row['total_new_attempts']==15
        for arm in row['arms']:
            assert arm['status']=='gradient_candidate_stability_unchecked' and arm['released_gradient_eV_A']<.03
            geometry=arm['actual_geometry']
            assert .12<sign*geometry['structural_polar_proxy_change_A'][1]<.123
            assert abs(geometry['structural_polar_proxy_change_A'][2])<1e-12
            assert geometry['structural_proxy_is_not_electronic_polarization']
        energies.append(row['arms'][1]['dE_minus_F_meV_fu'])
        secants=row['directional_secants_at_A']
        assert secants['at']=='A' and not secants['full_stability_certified']
        assert secants['gradient_eV_A2']>0 and secants['energy_eV_A2']>0
    assert abs(energies[0]-energies[1])<.001


def test_plus1_budget_failures_have_no_imputed_candidate_or_J_or_stability_check(report):
    rows=report['rows'][2:]
    assert len(rows)==2 and all(r['strain_percent']==1 for r in rows)
    for row in rows:
        a,j=row['arms']
        assert row['total_new_attempts']==12 and a['new_attempts']==12 and j['new_attempts']==0
        assert a['status']=='attempt_budget_exhausted' and j['status']=='not_attempted_atomic_candidate_unavailable'
        for arm in (a,j):
            assert arm['energy_eV_cell'] is arm['dE_minus_F_meV_fu'] is arm['actual_geometry'] is None
        assert row['directional_secants_at_A'] is None
        assert .4999<row['last_attempt_coordinate_norm_A']<.5
        assert row['last_attempt_atomic_released_gradient_eV_A']>.37
        assert row['last_attempt_is_not_automatically_a_candidate']
