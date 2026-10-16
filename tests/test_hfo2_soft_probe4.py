"""Four finite training statics must not become an expanded or independent run."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_soft_probe_E256_20261016'


@pytest.fixture(scope='module')
def worker():
    spec = importlib.util.spec_from_file_location('four_soft_training_probes',CASE/'worker.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_four_literal_training_points_have_locked_predictions_and_budget(worker):
    config = worker.prepare(CASE)
    assert config == json.loads((CASE/'config.json').read_bytes())
    data,anchor,direction,soft = worker.validate(config,CASE)
    assert len(config['points']) == 4
    assert config['policy']['t_A'] == [.05,-.05,.10,-.10]
    assert config['policy']['maximum_new_attempts'] == 4
    assert config['policy']['amended_training_committed'] == {'T':34,'preserving':39}
    assert config['policy']['independent_computed_each'] == 0
    assert config['original_independent_admission_refusal_preserved']
    assert not config['physical_plus005_generated_or_read']
    for p,t in zip(config['points'],config['policy']['t_A']):
        x = np.asarray(p['coordinates_full_A'])
        np.testing.assert_array_equal(x,anchor+t*direction)
        atoms = worker.checked_geometry(data,x)
        assert atoms.get_chemical_symbols() == ['Hf']*4+['O']*8
        assert np.linalg.norm(x) < .5
        np.testing.assert_allclose(data['axes'].T@(data['basis'].T@x),0.,atol=2e-12,rtol=0)
        expected = [soft['soft_gradient_eV_A']*t+.5*k*t*t for k in soft['first_direction_Rayleigh_eV_A2']]
        np.testing.assert_allclose(p['quadratic_delta_eV_cell'],expected,atol=2e-12,rtol=0)


@pytest.mark.parametrize('fault',['attempts','amplitude','independent','coordinate'])
def test_changed_budget_targets_or_test_status_are_refused(worker,fault):
    config = worker.prepare(CASE)
    if fault == 'attempts': config['policy']['maximum_new_attempts'] = 5
    elif fault == 'amplitude': config['policy']['t_A'][0] = .2
    elif fault == 'independent': config['policy']['independent_computed_each'] = 1
    else: config['points'][0]['coordinates_full_A'][0] += .001
    with pytest.raises(ValueError,match='contract differs'):
        worker.validate(config,CASE)
    assert worker.POLICY['maximum_new_attempts'] == 4  # no mutable policy alias


def test_literal_geometry_guard_does_not_repair_control_or_radius(worker):
    data,anchor,_,_ = worker.inputs(CASE)
    bad = anchor.copy(); bad[-1] = .001
    with pytest.raises(ValueError): worker.checked_geometry(data,bad)
    with pytest.raises(ValueError): worker.checked_geometry(data,anchor*3)
    bad = anchor + data['basis']@data['axes'][:,0]*.05
    with pytest.raises(ValueError): worker.checked_geometry(data,bad)


def test_run_refuses_local_shell_before_any_DFT(worker,monkeypatch,tmp_path):
    for name in ('RUN_DFT','SLURM_JOB_ID','SLURM_JOB_PARTITION'):
        monkeypatch.delenv(name,raising=False)
    with pytest.raises(ValueError,match='allocation required'):
        worker.run(worker.prepare(CASE),CASE,tmp_path/'must_not_exist')
    assert not (tmp_path/'must_not_exist').exists()


def test_preflight_refuses_existing_or_unapproved_output_namespace(worker,tmp_path):
    with pytest.raises(ValueError,match='namespace required'):
        worker.preflight(worker.prepare(CASE),CASE,tmp_path)
