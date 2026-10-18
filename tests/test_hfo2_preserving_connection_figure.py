"""Native first-update structure/force evidence, not a final barrier claim."""
from pathlib import Path

import numpy as np
import pytest
from scripts.plot_hfo2_preserving_connection import build_data,draw


@pytest.fixture(scope='module')
def data():
    return build_data(1)


def test_original_seed_and_one_actual_complete_update_replay(data):
    old,seed,new=data['bands']
    assert data['step']==1 and data['selected_new_SCF_calls']==7 and data['new_DFT_calls']==0
    assert old['sampled_polygon_peak_meV_fu']==pytest.approx(96.82818300279905,abs=1e-9)
    assert seed['sampled_polygon_peak_meV_fu']==pytest.approx(85.49988851655144,abs=1e-9)
    assert new['sampled_polygon_peak_meV_fu']==pytest.approx(85.26929141089568,abs=1e-9)
    assert old['ordinary_residual_pass']
    assert not seed['ordinary_residual_pass'] and not new['ordinary_residual_pass']
    assert new['fmax_eV_A']==pytest.approx(.1237397955929541,abs=1e-12)
    assert data['old_labels_and_forecast_scores_unchanged']
    assert not data['stationary_TS_or_sampling_bound_certified']


def test_same_endpoint_caches_and_actual_transverse_structural_proxy(data):
    old,seed,new=data['bands']
    for band in (seed,new):
        assert len(band['rows'])==9
        for i in (0,8):
            assert band['rows'][i]['energy_meV_fu']==old['rows'][i]['energy_meV_fu']
            assert band['rows'][i]['NEB_residual_eV_A'] is None
    for key in ('Qx_A','Qy_A','Qz_A','polar_proxy_z_A'):
        assert seed['rows'][4][key]==pytest.approx(old['rows'][4][key],abs=1e-12)
    assert new['rows'][4]['polar_proxy_y_A']>.12
    assert data['structural_proxies_are_not_electronic_polarization']
    assert data['one_chain_not_statistical_replicates']


def test_limiting_neighbor_force_is_not_a_spring_only_spacing_error(data):
    force=data['bands'][-1]['residual_decomposition']
    assert force['limiting_image_index']==3 and force['limiting_vector_block']=='atoms'
    winner=force['images'][2]
    p=winner['perpendicular_projection_on_worst_residual_eV_A']
    s=winner['spring_projection_on_worst_residual_eV_A']
    assert p>.12 and -.002<s<0
    assert p+s==pytest.approx(winner['residual_max_vector_eV_A'],abs=1e-12)
    assert all(r['component_sum_max_abs_error_eV_A']<1e-12 for r in force['images'])
    assert force['perpendicular_only_fmax_same_geometry_eV_A']>.12
    assert force['spring_only_fmax_same_geometry_eV_A']<.008


def test_editable_aligned_four_panel_exports_respect_author_style(data,tmp_path):
    qa=draw(data,tmp_path)
    assert qa['panel_letter_weight']=='normal' and qa['panel_letter_size_pt']==13
    assert qa['top_and_right_spines'] and not qa['subpanel_titles']
    assert qa['legend_white_background_alpha']==.85
    assert not qa['interpolation_or_error_bars_or_synthetic_points']
    b=np.asarray(qa['axes_bounds'])
    np.testing.assert_allclose(b[:2,1],b[0,1],atol=1e-12,rtol=0)
    np.testing.assert_allclose(b[[0,2],0],b[0,0],atol=1e-12,rtol=0)
    svg=(tmp_path/'connection_step_0001.svg').read_text(encoding='utf-8')
    assert '<text' in svg and 'Free-chain step 1' in svg
    assert (tmp_path/'connection_step_0001.pdf').stat().st_size>1000
