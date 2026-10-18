"""Physical covectors, exact fixed-cell embedding and finite attempt accounting."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from ase.io import read

from scripts.plot_hfo2_mechanism import original_chart_images, rotated_t_triplet, path_coordinates
from vcneb.mode_surface import _orthogonal_directions
from vcneb import ModePlane

REPOSITORY = Path(__file__).resolve().parents[1]
CASE = REPOSITORY/'benchmarks/hfo2_channels/20261008/matched_pattern_release_E270_20261018'
spec = importlib.util.spec_from_file_location('E270_patterns', CASE/'study.py')
study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)


@pytest.fixture(scope='module')
def prepared(tmp_path_factory):
    folder = tmp_path_factory.mktemp('E270')/'prepared'
    study.prepare(REPOSITORY, folder)
    return folder


def test_three_native_sources_and_fixed_contract_without_new_DFT(prepared):
    config = json.loads((prepared/'config.json').read_bytes())
    assert config['policy'] == study.POLICY
    assert config['runtime_root'] == '/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1'
    assert config['new_DFT_calls'] == 0 and not config['preserving_admitted']
    assert [r['source_job'] for r in config['frozen']] == ['28692775', '28938773', '28709788']
    assert [r['substrate_strain'] for r in config['frozen']] == [0., .005, .01]
    for i in range(3):
        seed = prepared/'seeds'/str(i); chart = study.chart_for(seed)
        energy, gradient = study.cached_frozen(seed, chart)
        assert energy == config['frozen'][i]['frozen_energy_eV_cell']
        assert gradient.shape == (39,) and np.isfinite(gradient).all()


def test_physical_covectors_equal_common_pattern_and_polar_changes_at_all_strains(prepared):
    t = read(REPOSITORY/study.BASE/'reference_variants/T.vasp')
    parent, patterns, _ = rotated_t_triplet(t)
    orientation = json.loads((REPOSITORY/study.BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json').read_bytes())['orientation']
    rng = np.random.default_rng(270)
    for i in range(3):
        seed = prepared/'seeds'/str(i); chart = study.chart_for(seed)
        c = np.asarray(json.loads((seed/'seed.json').read_bytes())['covectors_39_by_4'])
        for x in rng.normal(scale=.002, size=(4, 39)):
            images = original_chart_images([chart.reference_atoms, chart.to_atoms(x)], orientation, t)
            q, polar, _, _ = path_coordinates(images, parent, patterns)
            measured = np.r_[(q[1]-q[0])[[0, 2, 1]], polar[1, 2]-polar[0, 2]]
            np.testing.assert_allclose(c.T@x, measured, atol=1e-12, rtol=0)


def test_release_ranks_exact_cell_and_four_quantities_are_preserved(prepared):
    seed = prepared/'seeds/0'; chart = study.chart_for(seed)
    data = json.loads((seed/'seed.json').read_bytes())
    q = np.asarray(data['orthonormal_39_by_4'])
    translations = np.tile(np.eye(3), (12, 1))
    for atomic, rank in ((True, 29), (False, 32)):
        axes = q[:36, :2] if atomic else q[:, :2]
        gauge = translations if atomic else np.vstack([translations, np.zeros((3, 3))])
        frozen = np.column_stack([q[:len(axes), 2:], gauge])
        plane = ModePlane(np.zeros(len(axes)), axes, np.ones(len(axes)), np.eye(2), ('1', '2'), 'A', 'test')
        assert _orthogonal_directions(plane, frozen).shape == (len(axes), rank)
    rng = np.random.default_rng(729)
    gradient = rng.normal(scale=.025, size=39)
    def quadratic(x):
        return float(gradient@x+.5*x@x), gradient+x
    cache = study.ExactCache(quadratic, 0., gradient)
    recorded = []
    rows = study.treatments(chart, data, cache, record=recorded.append)
    assert recorded == rows
    assert [r['treatment'] for r in rows] == ['A', 'J']
    c = np.asarray(data['covectors_39_by_4'])
    for row in rows:
        assert row['candidate'] is not None and row['new_attempts'] <= 12
        point = row['candidate']; x = np.asarray(point['coordinates'])
        np.testing.assert_allclose(c.T@x, 0., atol=1e-10, rtol=0)
        assert point['released_gradient_norm_eV_A'] <= .03
        assert not point['stationary_TS_or_stability_certified']
        if row['treatment'] == 'A':
            assert not np.any(x[36:])
            assert np.array_equal(chart.to_atoms(x).cell.array, chart.reference_atoms.cell.array)
    assert rows[1]['candidate']['energy_eV_cell'] <= rows[0]['candidate']['energy_eV_cell'] <= 0
    assert cache.fresh == sum(r['new_attempts'] for r in rows)


def test_exact_cache_and_cap_do_not_substitute_nearby_geometry():
    calls = []
    def evaluate(x):
        calls.append(x.copy()); return float(x@x), 2*x
    cache = study.ExactCache(evaluate, 0., np.zeros(39))
    cache(np.zeros(39)); assert cache.fresh == 0 and len(calls) == 0
    for i in range(12):
        x = np.zeros(39); x[0] = (i+1)*1e-15
        cache(x); cache(x)
    assert len(calls) == cache.fresh == 12
    x = np.zeros(39); x[1] = 1e-15
    with pytest.raises(study.AttemptCap): cache(x)
    assert len(calls) == 12


def test_no_atomic_candidate_means_no_joint_attempt(prepared):
    seed = prepared/'seeds/0'; chart = study.chart_for(seed)
    data = json.loads((seed/'seed.json').read_bytes())
    q = np.asarray(data['orthonormal_39_by_4'])
    plane = ModePlane(np.zeros(36), q[:36, :2], np.ones(36), np.eye(2), ('1','2'), 'A', 'linear-test')
    frozen = np.column_stack([q[:36, 2:], np.tile(np.eye(3), (12, 1))])
    direction = np.r_[_orthogonal_directions(plane, frozen)[:, 0], np.zeros(3)]
    g = .2*direction
    cache = study.ExactCache(lambda x: (float(g@x), g), 0., g)
    rows = study.treatments(chart, data, cache)
    assert rows[0]['candidate'] is None and rows[0]['new_attempts'] <= 12
    assert rows[1]['status'] == 'not_attempted_atomic_candidate_unavailable'
    assert rows[1]['new_attempts'] == 0 and cache.fresh <= 12


def test_prepare_never_calls_calculator_and_refuses_existing_target(tmp_path, monkeypatch):
    from examples.hfo2_fixed_input_factory import FixedHfo2Calculator
    def forbidden(*args, **kwargs):
        raise AssertionError('prepare cannot launch DFT')
    monkeypatch.setattr(FixedHfo2Calculator, 'calculate', forbidden)
    destination = tmp_path/'prepared'
    study.prepare(REPOSITORY, destination)
    before = (destination/'config.json').read_bytes()
    with pytest.raises(FileExistsError): study.prepare(REPOSITORY, destination)
    assert (destination/'config.json').read_bytes() == before


def test_duplicate_submission_intent_blocks_sbatch(tmp_path, monkeypatch):
    monkeypatch.setattr(study, 'ROOT', tmp_path)
    monkeypatch.setattr(study, 'validate', lambda pin: {})
    monkeypatch.setattr(study, 'queue', lambda: ['28946180|hfo2-E269-midpoint-0|RUNNING|hfacnormal01|32'])
    study.save(tmp_path/'preflight.json', dict(bundle_pin='pin', rows=[0,1,2], new_DFT_calls=0))
    study.save(tmp_path/'intent.json', dict(command=['already-attempted']))
    def forbidden(*args, **kwargs):
        raise AssertionError('duplicate intent cannot submit')
    monkeypatch.setattr(study.subprocess, 'check_output', forbidden)
    with pytest.raises(ValueError, match='existing submission intent'):
        study.dispatch('submit-held', 'pin')
