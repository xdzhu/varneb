"""Actual four-corner energies and conditional interventions, not synthetic gains."""
import importlib.util
from pathlib import Path

import pytest

CASE=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/midpoint_factorial_E274_20261019'
spec=importlib.util.spec_from_file_location('E274_terminal_results',CASE/'analyze_terminal.py')
reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)


@pytest.fixture(scope='module')
def data():
    return reader.analyze()


def test_preserving_contrast_persists_at_both_fixed_cells(data):
    p=data['rows'][0]
    assert p['channel']=='preserving'
    assert p['diagonal_change_meV_fu']==pytest.approx(10.820843716373929,abs=1e-9)
    assert p['internal_at_A_cell_meV_fu']==pytest.approx(11.398469141568057,abs=1e-9)
    assert p['internal_at_B_cell_meV_fu']==pytest.approx(10.36220439709723,abs=1e-9)
    assert p['cell_at_A_internal_meV_fu']==pytest.approx(.4586393192766991,abs=1e-9)
    assert p['cell_at_B_internal_meV_fu']==pytest.approx(-.5776254251941282,abs=1e-9)
    assert p['finite_interaction_meV_fu']==pytest.approx(-1.0362647444708273,abs=1e-9)


def test_T_remains_small_under_both_orders_of_the_frozen_intervention(data):
    t=data['rows'][1]
    assert t['channel']=='T'
    assert t['diagonal_change_meV_fu']==pytest.approx(-.13193540735301212,abs=1e-9)
    for key in ('internal_at_A_cell_meV_fu','internal_at_B_cell_meV_fu',
                'cell_at_A_internal_meV_fu','cell_at_B_internal_meV_fu','finite_interaction_meV_fu'):
        assert abs(t[key])<.1
    for row in data['rows']:
        assert row['internal_at_A_cell_meV_fu']+row['cell_at_B_internal_meV_fu']==pytest.approx(row['diagonal_change_meV_fu'],abs=1e-8)
        assert row['cell_at_A_internal_meV_fu']+row['internal_at_B_cell_meV_fu']==pytest.approx(row['diagonal_change_meV_fu'],abs=1e-8)
        assert not row['unique_elastic_or_phonon_energy_partition']


def test_four_actual_new_calls_reuse_four_diagonal_energies_without_rescoring(data):
    assert data['source_job']=='28950839' and data['new_SCF_calls']==4
    assert data['all_four_new_native_EFS_reparsed_locally']
    assert data['native_diagonal_energies_replayed_locally']
    assert len(data['corners'])==8 and sum(r['new_SCF'] for r in data['corners'])==4
    assert data['native_proxy_core_hours']==pytest.approx(3.466703908211655)
    assert data['allocated_core_hours']==pytest.approx(3.493333333333333)
    assert data['new_DFT_calls_for_analysis']==0 and data['old_labels_and_forecasts_unchanged']
    assert not data['new_relaxed_path_or_TS_certified']
    assert not data['full_pseudo_orbital_bytes_rechecked_locally']
