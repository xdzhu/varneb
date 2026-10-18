"""Existing local energies acquire a common geometric interpretation, not new labels."""
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.analyze_hfo2_local_accommodation import analyze, displacement_parts, export, figure, POINTS

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def data():
    return analyze(ROOT)


def test_all_ten_physical_points_once_not_two_probe_step_replicates(data):
    assert data['physical_points'] == 10 and data['geometries'] == 30
    assert len(data['states']) == 30 and len(data['changes']) == 20
    assert data['new_DFT_calls'] == 0 and not data['new_native_audit']
    assert not data['independent_prediction_test']
    assert not data['energy_partition_or_TS_or_stable_surface_certified']
    assert [r['source_job'] for r in data['centers']] == ['28727382', '28692775']
    for channel in ('preserving', 'T'):
        states = [r for r in data['states'] if r['channel'] == channel]
        assert [(r['point_id'], r['treatment']) for r in states] == [
            (i, t) for i in range(5) for t in ('F', 'A', 'J')]
        assert [[r['Q1_A'], r['Q2_A']] for r in states if r['treatment'] == 'F'] == POINTS
        for i in range(5):
            f, a, _ = [r for r in states if r['point_id'] == i]
            assert a['open_y_stretch_percent'] == f['open_y_stretch_percent']


def test_original_energy_contrasts_and_observed_geometry_not_refitted(data):
    frozen_origins = {r['channel']: r for r in data['states']
                      if r['point_id'] == 0 and r['treatment'] == 'F'}
    assert frozen_origins['preserving']['Qx_A'] == pytest.approx(-.7992885087386647)
    assert frozen_origins['T']['Qz_A'] == pytest.approx(.13662528423958056)
    assert frozen_origins['preserving']['open_y_stretch_percent'] == pytest.approx(-.057769867354418025)
    origins = {r['channel']: r for r in data['changes']
               if r['point_id'] == 0 and r['contrast'] == 'A_minus_F'}
    p, t = origins['preserving'], origins['T']
    assert p['energy_change_meV_fu'] == pytest.approx(-12.087172655810718)
    assert t['energy_change_meV_fu'] == pytest.approx(-.06784954257454956)
    assert p['delta_Qx_A'] == pytest.approx(.0627939257813048)
    assert p['total_norm_A'] == pytest.approx(.3097750967954315)
    assert p['complement_norm_A'] == pytest.approx(.3033439191994945)
    assert p['oxygen_squared_norm_fraction'] == pytest.approx(.8114238621059985)
    assert t['triplet_norm_A'] == pytest.approx(.024477471787310284)
    for row in data['changes']:
        assert row['total_norm_A']**2 == pytest.approx(row['Hf_norm_A']**2+row['oxygen_norm_A']**2)
        assert row['total_norm_A']**2 == pytest.approx(row['triplet_norm_A']**2+row['complement_norm_A']**2)
        assert row['triplet_norm_A']**2 == pytest.approx(sum(row[f'delta_Q{d}_A']**2 for d in 'xyz'))
        assert row['geometric_fractions_are_not_energy_fractions']
    preserving = [r for r in data['changes'] if r['channel'] == 'preserving' and r['contrast'] == 'A_minus_F']
    assert all(.958 < r['complement_squared_norm_fraction'] < .971 for r in preserving)
    assert all(.808 < r['oxygen_squared_norm_fraction'] < .813 for r in preserving)


def test_geometry_translation_and_complete_subspace_rotation(data):
    rng = np.random.default_rng(421)
    first, second = rng.normal(size=(2, 12, 3))
    patterns = np.asarray(data['pattern_basis'])
    a = displacement_parts(first, second, patterns)
    b = displacement_parts(first+np.array([1, 2, 3]), second+np.array([-2, .4, 1]), patterns)
    mixing = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    c = displacement_parts(first, second, np.einsum('ab,bij->aij', mixing, patterns))
    for key in ('total_norm_A', 'Hf_norm_A', 'oxygen_norm_A', 'triplet_norm_A', 'complement_norm_A',
                'oxygen_squared_norm_fraction', 'complement_squared_norm_fraction'):
        assert a[key] == pytest.approx(b[key])
        assert a[key] == pytest.approx(c[key])
    zero = displacement_parts(first, first, patterns)
    assert zero['total_norm_A'] == 0
    assert zero['oxygen_squared_norm_fraction'] is None and zero['complement_squared_norm_fraction'] is None


@pytest.mark.parametrize('fault', ['shape', 'nonfinite', 'nonorthogonal', 'translation'])
def test_bad_displacement_or_reference_refused(data, fault):
    first, second, modes = np.zeros((12, 3)), np.ones((12, 3)), np.array(data['pattern_basis'])
    if fault == 'shape': first = first[:11]
    elif fault == 'nonfinite': second[0, 0] = np.nan
    elif fault == 'nonorthogonal': modes[0] *= 2
    else: modes[0] += 1
    with pytest.raises(ValueError):
        displacement_parts(first, second, modes)


def test_analysis_cannot_launch_backend_or_scheduler(monkeypatch):
    import subprocess
    from examples.hfo2_fixed_input_factory import FixedHfo2Calculator
    def forbidden(*args, **kwargs):
        raise AssertionError('a geometric replay must not run DFT or query/mutate the scheduler')
    monkeypatch.setattr(subprocess, 'run', forbidden)
    monkeypatch.setattr(subprocess, 'check_output', forbidden)
    monkeypatch.setattr(FixedHfo2Calculator, 'calculate', forbidden)
    assert analyze(ROOT)['physical_points'] == 10


def test_changed_archived_source_pin_is_refused(monkeypatch):
    from scripts import analyze_hfo2_local_accommodation as module
    pins = dict(module.PINS)
    pins[next(iter(pins))] = '0'*64
    monkeypatch.setattr(module, 'PINS', pins)
    with pytest.raises(ValueError, match='archived treatment source changed'):
        module.analyze(ROOT)


def test_export_preserves_numeric_source_and_refuses_existing_target(data, tmp_path):
    output = tmp_path/'output'
    export(data, output, draw=False)
    assert json.loads((output/'analysis.json').read_bytes()) == data
    for name, count in (('states', 30), ('changes', 20)):
        with (output/f'{name}.csv').open() as handle:
            assert len(list(csv.DictReader(handle))) == count
    before = (output/'analysis.json').read_bytes()
    with pytest.raises(FileExistsError):
        export(data, output, draw=False)
    assert (output/'analysis.json').read_bytes() == before


def test_four_panels_keep_points_alignment_labels_and_explicitly_unequal_scales(data):
    import matplotlib.pyplot as plt
    fig = figure(data)
    try:
        assert len(fig.axes) == 4 and len(fig.legends) == 2
        for ax, label in zip(fig.axes, ('(a)', '(b)', '(c)', '(d)')):
            assert len(ax.patches) == 10 and not ax.get_title()
            assert all(s.get_visible() for s in ax.spines.values())
            text = next(t for t in ax.texts if t.get_text() == label)
            assert text.get_fontweight() == 'normal' and text.get_fontsize() == 13
            assert text.get_position()[0] == 0
        assert fig.axes[0].get_position().x0 == fig.axes[2].get_position().x0
        assert fig.axes[1].get_position().x0 == fig.axes[3].get_position().x0
        assert fig.axes[0].get_position().y0 == fig.axes[1].get_position().y0
        assert fig.axes[2].get_position().y0 == fig.axes[3].get_position().y0
        assert fig.axes[0].get_ylim() != fig.axes[1].get_ylim()
    finally:
        plt.close(fig)
