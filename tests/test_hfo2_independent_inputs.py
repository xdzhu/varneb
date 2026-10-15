"""Technical checks only; no independent material label or calculator run."""
from pathlib import Path
import copy
import importlib.util
import json

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_independent_inputs_E186_20261015'
SOURCE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_local_forecast_E185_20261015'
spec = importlib.util.spec_from_file_location('E186_inputs', CASE / 'prepare.py')
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


def test_exact_eight_targets_each_and_reserved_budget_are_inherited():
    c = p.prepare(SOURCE)
    assert sum(len(r['points']) for r in c['channels']) == 16
    assert c['future_attempt_limit_total'] == 192
    assert c['independent_points_computed_each'] == 0
    assert c['production_test_admitted'] is False
    assert c['six_asset_sha256'] == p.study.CONTRACT
    assert c['maximum_internal_radius_A'] == .5 and c['future_gradient_target_eV_A'] == .03


def test_literal_full_start_is_preserved_not_reconstructed():
    locked = p.forecast(SOURCE)
    for row in locked['channels']:
        data = p.faces.load_inputs(row['channel'])
        for point in row['targets']:
            atoms = p.checked_start(data, point)
            original = data['chart'].to_atoms(point['starting_full_coordinates_A'])
            np.testing.assert_array_equal(atoms.positions, original.positions)
            np.testing.assert_array_equal(atoms.cell.array, original.cell.array)
            assert atoms.get_chemical_symbols() == ['Hf']*4 + ['O']*8


@pytest.mark.parametrize('field,value',[
    ('q_A',[.049,0.]), ('training_parent_id',1),
    ('starting_full_coordinates_A',[float('nan')]*40),
    ('starting_internal_coordinates_A',[.51]*36)])
def test_inconsistent_start_is_rejected_without_repair(field,value):
    point = copy.deepcopy(p.forecast(SOURCE)['channels'][0]['targets'][0])
    point[field] = value
    with pytest.raises(ValueError, match='frozen start violates'):
        p.checked_start(p.faces.load_inputs('T'), point)


def test_native_input_preflight_cannot_be_used_from_local_workspace(tmp_path):
    with pytest.raises(ValueError, match='immutable HF runtime'):
        p.native_preflight(p.prepare(SOURCE), SOURCE, tmp_path/'initial_decks')


def test_config_reproduces_portable_frozen_file():
    assert json.loads(json.dumps(p.prepare(SOURCE))) == json.loads((CASE/'config.json').read_bytes())


def test_entrypoint_has_no_submit_or_DFT_run_mode(monkeypatch):
    monkeypatch.setattr('sys.argv',['prepare.py','run','--bundle','x','--config','x'])
    with pytest.raises(SystemExit) as exc:
        p.main()
    assert exc.value.code == 2
