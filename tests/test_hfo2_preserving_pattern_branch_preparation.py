"""Measured omitted direction, matched physical controls and no new DFT."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/preserving_pattern_branches_E272_20261018'
spec = importlib.util.spec_from_file_location('P_branch_preparation', CASE/'prepare.py')
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


@pytest.fixture(scope='module')
def prepared(tmp_path_factory):
    destination = tmp_path_factory.mktemp('E272')/'prepared'
    result = module.prepare(destination)
    return destination, result


def test_only_original_extreme_strains_prepared_and_midpoint_not_admitted(prepared):
    folder, result = prepared
    assert [r['strain'] for r in result['rows']] == [0., .01]
    assert not result['execution_admitted'] and result['new_DFT_calls'] == result['new_submissions'] == 0
    assert not result['policy']['midpoint_admitted']
    assert not result['policy']['new_full_Hessian']
    assert result['policy']['maximum_prospective_attempts'] == 104
    assert result['policy']['maximum_prospective_optimization_attempts'] == 96
    assert result['policy']['maximum_directional_check_attempts'] == 8
    assert not (folder/'seeds/1').exists()


def test_both_measured_probe_amplitudes_show_same_omitted_atomic_direction(prepared):
    _, result = prepared
    rows = result['diagnostics']
    assert [r['probe_step_A'] for r in rows] == [.01, .02]
    assert [r['atomic_min_eigenvalue_eV_A2'] for r in rows] == pytest.approx([-.4474490104, -.455585944])
    assert all(r['atomic_negative_count'] == 1 for r in rows)
    assert all(abs(r['gradient_along_signed_direction_eV_A']) < 1e-6 for r in rows)
    assert all(.78 < r['oxygen_squared_norm_fraction'] < .79 for r in rows)
    assert result['signed_direction_probe_overlap'] > .999999
    assert all(not r['physical_stationarity_or_stability_certified'] for r in rows)


def test_signed_starts_keep_all_four_pattern_values_and_exact_cells(prepared):
    folder, _ = prepared
    for index in (0, 2):
        seed = folder/'seeds'/str(index)
        data = json.loads((seed/'seed.json').read_bytes())
        branch = json.loads((seed/'branch_preparation.json').read_bytes())
        c = np.asarray(data['covectors_39_by_4']); starts = np.asarray(branch['signed_starts'])
        assert starts.shape == (2, 39)
        assert not np.any(starts[:, 36:])
        np.testing.assert_allclose(starts[0], -starts[1], atol=1e-15, rtol=0)
        np.testing.assert_allclose(np.linalg.norm(starts, axis=1), .1, atol=1e-12, rtol=0)
        np.testing.assert_allclose(starts@c, 0., atol=1e-12, rtol=0)
        np.testing.assert_allclose(starts[:, :36].reshape(2,12,3).sum(axis=1), 0., atol=1e-12, rtol=0)


def test_numerical_seed_positive_without_removing_physical_negative_curvature(prepared):
    folder, result = prepared
    assert not result['policy']['eigenvalue_clipping']
    for index in (0, 2):
        seed = folder/'seeds'/str(index)
        data = json.loads((seed/'seed.json').read_bytes())
        branch = json.loads((seed/'branch_preparation.json').read_bytes())
        h = np.asarray(branch['numerical_hessian_39'])
        np.linalg.cholesky(h)
        for n in (36, 39):
            u = module.release(data, n)
            np.linalg.cholesky(u.T@h[:n, :n]@u)
    assert result['diagnostics'][0]['atomic_min_eigenvalue_eV_A2'] < 0


def test_covectors_equal_measured_original_T_pattern_and_polar_changes(prepared):
    from ase.io import read
    from scripts.plot_hfo2_mechanism import original_chart_images, rotated_t_triplet, path_coordinates
    study = module.module(module.BASE/'matched_pattern_release_E270_20261018/study.py', 'E272_chart_test')
    t = read(module.BASE/'reference_variants/T.vasp'); parent, patterns, _ = rotated_t_triplet(t)
    orientation = json.loads((module.BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json').read_bytes())['orientation']
    folder, _ = prepared; rng = np.random.default_rng(272)
    for index in (0, 2):
        seed = folder/'seeds'/str(index); chart = study.chart_for(seed)
        c = np.asarray(json.loads((seed/'seed.json').read_bytes())['covectors_39_by_4'])
        for x in rng.normal(scale=.002, size=(3,39)):
            images = original_chart_images([chart.reference_atoms, chart.to_atoms(x)], orientation, t)
            q, polar, _, _ = path_coordinates(images, parent, patterns)
            np.testing.assert_allclose(c.T@x, np.r_[(q[1]-q[0])[[0,2,1]], polar[1,2]-polar[0,2]], atol=1e-12, rtol=0)


def test_preparation_cannot_launch_DFT_or_overwrite(tmp_path, monkeypatch):
    from examples.hfo2_fixed_input_factory import FixedHfo2Calculator
    monkeypatch.setattr(FixedHfo2Calculator, 'calculate', lambda *a, **k: pytest.fail('no DFT admitted'))
    destination = tmp_path/'prepared'
    module.prepare(destination)
    with pytest.raises(FileExistsError): module.prepare(destination)
