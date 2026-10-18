import numpy as np
import pytest
from scripts.analyze_hfo2_oxygen_motion import analyze


@pytest.fixture(scope='module')
def report(): return analyze()


def test_zero_species_mean_does_not_mean_straight_oxygen_motion(report):
    old = report['bands']['original']
    assert abs(old['center_O_mean_relative_y_A']) < 1e-4
    assert old['center_O_max_abs_relative_y_A'] > .32
    assert old['center_O_relative_y_nonuniform_norm_A'] > .7
    assert not report['literature_CN_ST_family_assigned']
    assert not report['electronic_polarization_computed']


def test_lower_candidate_is_mixed_motion_and_retains_ordered_layer_bands(report):
    delta = np.array(report['center_lower_minus_original_translation_free_A'])
    np.testing.assert_allclose(delta.mean(axis=0), 0., atol=1e-16)
    assert np.max(np.abs(delta[:4, 0])) > .067
    assert np.max(delta[4:, 2]) > .20 and np.min(delta[4:, 2]) < -.20
    assert report['new_DFT_calls'] == 0
    for band in report['bands'].values():
        assert band['no_Hf_layer_band_crossing_at_samples']
        assert band['minimum_layer_band_margin_fractional_z'] > .14
    # This is relative to the initial endpoint, not the absolute T-parent proxy.
    assert report['bands']['lower']['center_O_mean_relative_y_A'] == pytest.approx(.1210794665, abs=1e-9)
