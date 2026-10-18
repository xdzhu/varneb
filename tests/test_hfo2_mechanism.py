"""Mechanism figures replay archived paths; no DFT or scheduler calls."""
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from ase.io import read

from scripts.plot_hfo2_mechanism import (
    BASE, NETWORK, build_data, load, path_coordinates, response_rows,
    rotated_t_triplet, competition_figure, geometry_figure, treatment_figure,
)

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def data():
    return build_data(ROOT)


def test_common_initial_offset_and_opposite_switching_response(data):
    rows = {r['channel']: r for r in data['responses']}
    for r in rows.values():
        assert r['delta_initial_meV_fu'] == pytest.approx(6.806502377)
        assert r['delta_absolute_sampled_peak_meV_fu']-r['delta_initial_meV_fu'] == pytest.approx(r['delta_barrier_meV_fu'])
    assert rows['preserving']['delta_barrier_meV_fu'] == pytest.approx(28.596248155)
    assert rows['reversing']['delta_barrier_meV_fu'] == pytest.approx(-44.024330642)
    assert rows['preserving']['endpoint_only_delta_barrier_meV_fu'] == pytest.approx(rows['reversing']['endpoint_only_delta_barrier_meV_fu'])
    assert data['new_DFT_calls'] == 0 and not data['full_native_audit_repeated']


@pytest.mark.parametrize('fault', ['missing', 'duplicate', 'initial'])
def test_invalid_response_comparison_refused(fault):
    rows = deepcopy(load(ROOT/NETWORK)['rows'])
    if fault == 'missing': rows.pop()
    elif fault == 'duplicate': rows[-1] = deepcopy(rows[0])
    else: rows[0]['matched_initial_energy_eV_cell'] += .001
    with pytest.raises(ValueError): response_rows(rows)


def test_all_eight_bands_keep_ordered_structural_distinctions(data):
    assert len(data['geometry']) == 72 and len(data['bands']) == 8
    basis = np.array(data['basis_cartesian_A_unit_vectors']).reshape(3, 36)
    np.testing.assert_allclose(basis@basis.T, np.eye(3), atol=1e-12)
    for band in data['bands']:
        rows = band['rows']
        assert len(rows) == 9 and [r['image'] for r in rows] == list(range(9))
        assert np.all(np.diff([r['s'] for r in rows]) > 0)
        np.testing.assert_allclose([r['F_xx'] for r in rows], 1+band['strain'], atol=1e-10)
        np.testing.assert_allclose([r['F_zz'] for r in rows], 1+band['strain'], atol=1e-10)
        if band['channel'] == 'preserving':
            assert all(r['Q_x_A'] < 0 for r in rows)
        if band['channel'] == 'reversing':
            assert rows[0]['Q_x_A'] < 0 < rows[-1]['Q_x_A']
            assert rows[band['peak_image']]['residual_atomic_A'] > 1.3
        if band['channel'] in ('preserving', 'reversing'):
            assert rows[0]['polar_dz_A'] > .27 and rows[-1]['polar_dz_A'] < -.27
        if band['channel'] == 'T':
            assert abs(rows[-1]['Q_x_A']) < 1e-6
            assert rows[-1]['Q_z_A'] > .8


def test_projection_keeps_global_translation_and_integer_gauge_invariance():
    t = read(ROOT/BASE/'reference_variants/T.vasp')
    parent, basis, _ = rotated_t_triplet(t)
    images = [t.copy(), t.copy()]
    q, polar, _, residual = path_coordinates(images, parent, basis)
    for atoms in images:
        scaled = atoms.get_scaled_positions(wrap=False)
        scaled += [0.01, -.02, .03]
        scaled[0] += [1, -2, 0]
        atoms.set_scaled_positions(scaled)
    other_q, other_polar, _, other_residual = path_coordinates(images, parent, basis)
    np.testing.assert_allclose(q, other_q, atol=1e-12)
    np.testing.assert_allclose(polar, other_polar, atol=1e-12)
    np.testing.assert_allclose(residual, other_residual, atol=1e-12)


def test_treatment_points_not_counted_twice_as_curvature_replicates(data):
    assert len(data['treatments']) == 10
    for channel in ('preserving', 'T'):
        rows = [r for r in data['treatments'] if r['channel'] == channel]
        assert len(rows) == 5 and len({(r['Q1_A'], r['Q2_A']) for r in rows}) == 5
        for r in rows:
            assert r['atomic_minus_frozen_meV_fu']+r['joint_minus_atomic_meV_fu'] == pytest.approx(r['joint_minus_frozen_meV_fu'])
            assert r['atomic_minus_frozen_meV_fu'] < 0 and r['joint_minus_atomic_meV_fu'] < 0
        lowerings = [-r['atomic_minus_frozen_meV_fu'] for r in rows]
        assert (min(lowerings), max(lowerings)) == pytest.approx(
            (11.015956, 12.087173) if channel == 'preserving' else (.067850, .166184), abs=1e-6)


@pytest.mark.parametrize('builder', [competition_figure, geometry_figure, treatment_figure])
def test_style_alignment_and_all_data_inside_axes(data, builder):
    fig = builder(data)
    try:
        renderer = fig.canvas.get_renderer()
        canvas = fig.bbox
        for i, ax in enumerate(fig.axes):
            assert not ax.get_title() and all(s.get_visible() for s in ax.spines.values())
            panel = next(t for t in ax.texts if t.get_text() == f'({chr(97+i)})')
            assert panel.get_fontweight() == 'normal' and panel.get_fontsize() >= 13
            for text in [panel, ax.xaxis.label, ax.yaxis.label, *ax.get_xticklabels(), *ax.get_yticklabels()]:
                bbox = text.get_window_extent(renderer)
                assert bbox.x0 >= canvas.x0-1 and bbox.x1 <= canvas.x1+1
                assert bbox.y0 >= canvas.y0-1 and bbox.y1 <= canvas.y1+1
            ymin, ymax = ax.get_ylim()
            for line in ax.lines:
                if line.get_transform() == ax.transData:
                    y = np.asarray(line.get_ydata(), dtype=float)
                    assert np.min(y) >= ymin and np.max(y) <= ymax
            for patch in ax.patches:
                assert ymin <= patch.get_y() <= ymax
                assert ymin <= patch.get_y()+patch.get_height() <= ymax
        for legend in fig.legends:
            assert 0 < legend.get_frame().get_alpha() < 1
        rows = [fig.axes] if builder != geometry_figure else np.array(fig.axes).reshape(3, 3)
        for row in rows:
            assert np.ptp([ax.get_position().y0 for ax in row]) < 1e-12
            assert np.ptp([ax.get_position().height for ax in row]) < 1e-12
    finally:
        import matplotlib.pyplot as plt
        plt.close(fig)


def test_main_text_connects_three_evidence_levels_and_preserves_old_figures():
    text = (ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    main = text.split('## 3.', 1)[1].split('## Appendix A.', 1)[0]
    for name in ('hfo2_competition.pdf', 'hfo2_path_geometry.pdf', 'hfo2_accommodation.pdf'):
        assert name in main
    assert (ROOT/'paper/VARNEB_JCTC/figures/hfo2_mechanism_20261018/hfo2_relaxation.pdf').is_file()
    assert 'these are not energy percentages' in main
    for n in range(1, 5): assert f'Figure {n}.' in main
    for n in range(1, 9): assert f'Figure A{n}.' in text
    assert 'correlation of a stretch with a peak' in main
    assert 'not an isolated causal test' in main
    assert 'no metric-independent' in text or 'not show a\nmetric-independent' in text
    assert 'Appendix F. Reference coverage' in text
    for n in range(20, 23): assert f'Table A{n}.' in text
