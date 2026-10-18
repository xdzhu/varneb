"""Known curvature is a possible numerical initializer, not a new prediction."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

CASE = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/matched_pattern_release_E270_20261018'
spec = importlib.util.spec_from_file_location('E270_block_seed', CASE/'numerical_seed.py')
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


@pytest.mark.parametrize('probe', [0, 1])
@pytest.mark.parametrize('index', [0, 1, 2])
def test_both_probe_steps_are_positive_on_all_target_release_spaces(index, probe):
    seed, report = module.seed_for(CASE, index, probe_index=probe)
    np.testing.assert_allclose(seed, seed.T, atol=1e-12, rtol=0)
    assert report['measured_probe_step_A'] == [.01, .02][probe]
    assert report['target_strain'] == [0., .005, .01][index]
    assert [r['rank'] for r in report['restrictions']] == [29, 32]
    assert all(r['minimum_seed_eigenvalue_eV_A2'] > 0 for r in report['restrictions'])
    assert not report['eigenvalue_clipping'] and report['new_DFT_calls'] == 0
    assert report['source_DFT_mixed_block_still_present']
    assert not report['execution_admitted'] and not report['physical_strain_Hessian_or_stability_prediction']


def test_transport_matches_actual_fractional_and_open_cell_increments():
    spec = importlib.util.spec_from_file_location('E270_chart_test', CASE/'study.py')
    study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)
    old = study.chart_for(CASE/'prepared/seeds/0')
    rng = np.random.default_rng(271)
    for i in range(3):
        target = study.chart_for(CASE/'prepared/seeds'/str(i))
        transform = module.chart_transport(old, target)
        cs, ct = [np.asarray(json.loads((CASE/'prepared/seeds'/str(j)/'seed.json').read_bytes())['covectors_39_by_4'])
                  for j in (0, i)]
        np.testing.assert_allclose(transform.T@cs, ct, atol=1e-12, rtol=0)
        for x in rng.normal(scale=.002, size=(3,39)):
            source_increment = transform@x
            a, b = old.to_atoms(source_increment), target.to_atoms(x)
            np.testing.assert_allclose(a.get_scaled_positions(wrap=False)-old.reference_atoms.get_scaled_positions(wrap=False),
                b.get_scaled_positions(wrap=False)-target.reference_atoms.get_scaled_positions(wrap=False), atol=1e-12, rtol=0)
            np.testing.assert_allclose(a.cell.array-old.reference_cell, b.cell.array-target.reference_cell, atol=1e-12, rtol=0)
