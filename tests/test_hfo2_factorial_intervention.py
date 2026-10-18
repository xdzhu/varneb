"""Frozen four-state interventions: geometry and contrast, not DFT results."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/midpoint_reverse_E269_20261018'
spec = importlib.util.spec_from_file_location('E269_factorial', CASE/'prepare_factorial.py')
study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)


def test_affine_exchange_keeps_fractions_not_Cartesian_atoms():
    a = Atoms('HfO', cell=np.eye(3)*6, scaled_positions=[[.2, .3, .4], [.6, .5, .2]], pbc=True)
    b = a.copy(); cell = b.cell.array.copy(); cell[2] = [.02, -.01, 6.1]
    b.set_cell(cell, scale_atoms=True)
    result = study.crossover(a, b)
    np.testing.assert_allclose(result.cell.array, b.cell.array)
    np.testing.assert_allclose(result.get_scaled_positions(wrap=False), a.get_scaled_positions(wrap=False))
    assert not np.allclose(result.positions, a.positions)
    assert result.calc is None
    bad = b.copy(); cell = bad.cell.array.copy(); cell[0, 0] += .01; bad.set_cell(cell)
    with pytest.raises(ValueError, match='substrate'): study.crossover(a, bad)


def test_finite_interaction_is_not_double_counted_or_called_mode_energy():
    # Synthetic arithmetic only: eV/cell -> meV/fu at four formula units.
    energies = {'AA': -100., 'BA': -99.96, 'AB': -99.98, 'BB': -99.92}
    result = study.contrasts(energies)
    assert result['diagonal_change_meV_fu'] == pytest.approx(20.)
    assert result['internal_at_A_cell_meV_fu'] == pytest.approx(10.)
    assert result['cell_at_A_internal_meV_fu'] == pytest.approx(5.)
    assert result['finite_interaction_meV_fu'] == pytest.approx(5.)
    assert result['diagonal_change_meV_fu'] == pytest.approx(result['internal_at_A_cell_meV_fu']+
        result['cell_at_A_internal_meV_fu']+result['finite_interaction_meV_fu'])
    shifted = study.contrasts({k: v-10000 for k, v in energies.items()})
    assert shifted['finite_interaction_meV_fu'] == pytest.approx(result['finite_interaction_meV_fu'], abs=1e-9)
    assert not result['unique_elastic_or_phonon_energy_partition']


@pytest.mark.parametrize('value', [None, float('nan'), True])
def test_unmeasured_corner_is_not_scored_as_zero(value):
    with pytest.raises(ValueError):
        study.contrasts({'AA': -100., 'AB': value, 'BA': -99.96, 'BB': -99.92})


def test_actual_four_decks_are_inert_and_keep_fixed_inputs(tmp_path):
    output = tmp_path/'prepared'
    result = study.prepare(ROOT, CASE/'initial_complete_chains', output, initial=True)
    assert result['additional_SCFs_required'] == len(result['decks']) == 4
    assert result['new_DFT_calls'] == result['jobs_submitted'] == 0
    assert not result['execution_admitted'] and not result['diagonals_recalculated']
    assert result['upper_stage'] == 'initial_evaluation'
    for row in result['channels']:
        assert row['corner_energy_eV_cell']['AB'] is row['corner_energy_eV_cell']['BA'] is None
        assert not row['measured_crossover_energies_available']
    for row in result['decks']:
        assert all(row['input_sha256'][n] == CONTRACT[n] for n in ('INPUT', 'KPT'))
        atoms = read_fixed_hfo2_stru(output/row['path']/'STRU')
        assert atoms.get_chemical_symbols() == ['Hf']*4+['O']*8
        assert atoms.calc is None and not row['actual_HF_ASE_preflight_done']
    original = (output/'preparation.json').read_bytes()
    with pytest.raises(FileExistsError): study.prepare(ROOT, CASE/'initial_complete_chains', output, initial=True)
    assert original == (output/'preparation.json').read_bytes()
    with pytest.raises(ValueError, match='initial snapshots'):
        study.prepare(ROOT, CASE/'initial_complete_chains', tmp_path/'refused', initial=False)
