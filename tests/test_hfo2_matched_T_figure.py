"""Completed native T controls only; geometry is not an energy partition."""
import numpy as np
import pytest

from scripts.plot_hfo2_matched_T_release import build_data, draw


@pytest.fixture(scope='module')
def data(): return build_data()


def test_native_energies_geometry_and_actual_physical_controls(data):
    rows=data['rows']
    assert [r['strain_percent'] for r in rows] == [0., .5, 1.]
    assert data['native_candidates_independently_replayed'] == 6 and data['new_DFT_calls'] == 0
    assert [r['dE_J_minus_F_meV_fu'] for r in rows] == pytest.approx([-.17737459666,-.08609264614,-.24002449709])
    assert [r['atomic_displacement_F_J_norm_A'] for r in rows] == pytest.approx([.07033733872,.02186298047,.09724638223])
    assert [r['joint_minus_frozen_open_y_stretch_percent'] for r in rows] == pytest.approx([.1409941389,.0993404427,.1013058837])
    assert all(r['atomic_displacement_A_J_norm_A'] > .019 for r in rows)
    assert max(r['measured_four_pattern_defect_A'] for r in rows) < 1e-12
    assert data['older_two_axis_controls_not_nested_with_this_four_pattern_control']
    assert data['geometric_norms_are_not_RMS_or_energy_fractions']
    assert data['local_energy_lowering_is_not_barrier_error_bound_or_stability_certificate']


def test_editable_aligned_exports_follow_author_style(data,tmp_path):
    from pypdf import PdfReader
    qa=draw(data,tmp_path)
    positions=np.asarray(qa['axes_bounds'])
    np.testing.assert_allclose(positions[:,1],positions[0,1],atol=1e-12,rtol=0)
    np.testing.assert_allclose(positions[:,3],positions[0,3],atol=1e-12,rtol=0)
    assert qa['panel_letter_weight']=='normal' and qa['panel_letter_size_pt']==13
    assert qa['top_and_right_spines'] and qa['framed_white_legend_alpha']==.85
    assert not qa['smoothing_or_error_bars_or_synthetic_data']
    svg=(tmp_path/'matched_T_release.svg').read_text(encoding='utf-8')
    assert svg.count('<text') > 10
    assert len(PdfReader(tmp_path/'matched_T_release.pdf').pages[0].extract_text()) > 100
    assert (tmp_path/'matched_T_release.png').stat().st_size > 10000
