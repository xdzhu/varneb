"""Execute new comparison definitions on archived real reference bands only."""
import copy
import importlib.util
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019'


def reader():
    spec = importlib.util.spec_from_file_location('E275_analyze_test', CASE/'analyze_terminal.py')
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


@pytest.fixture(scope='module')
def reference(): return reader().reference_data()


def test_real_reference_bands_keep_old_labels_and_resolve_only_zero_candidate(reference):
    assert reference['new_DFT_calls'] == 0 and reference['old_labels_and_scores_unchanged']
    assert len(reference['bands']) == 8
    zero = {b['kind']: b for b in reference['bands'] if b['strain'] == 0.}
    assert zero['registered_preserving']['sampled_polygon_peak_meV_fu'] == pytest.approx(96.8281830028)
    assert zero['released_seed']['sampled_polygon_peak_meV_fu'] == pytest.approx(83.1662809401)
    assert zero['T_escape']['sampled_polygon_peak_meV_fu'] == pytest.approx(84.4166238812)
    for b in reference['bands']:
        assert b['ordinary_residual_pass']
        assert b['rows'][0]['NEB_residual_eV_A'] is b['rows'][-1]['NEB_residual_eV_A'] is None
    # Actual geometric labels, with the original tolerance sweep; not a new phase.
    old_symmetry = zero['registered_preserving']['peak_structure_audit']['symmetry_sweep']
    new_symmetry = zero['released_seed']['peak_structure_audit']['symmetry_sweep']
    assert [r['symprec_A'] for r in old_symmetry] == [.001, .01, .05]
    assert [r['number'] for r in old_symmetry] == [60]*3
    assert [r['number'] for r in new_symmetry] == [41]*3
    assert zero['released_seed']['peak_structure_not_a_stable_phase_or_electronic_polarization']
    assert not any(b['kind'] == 'released_seed' and b['strain'] != 0. for b in reference['bands'])
    upper = next(b for b in reference['bands'] if b['kind'] == 'upper_seed_midpoint')
    assert upper['sampled_polygon_peak_meV_fu'] == pytest.approx(109.3517206327)


def test_competition_requires_actual_common_native_well(reference):
    m = reader(); bands = [b for b in reference['bands'] if b['strain'] == 0.]
    original = next(b for b in bands if b['kind'] == 'registered_preserving')
    new = next(b for b in bands if b['kind'] == 'released_seed')
    m.require_same_well(original, new, final=True)
    wrong = copy.deepcopy(new); wrong['endpoint_native_records'][0]['energy_eV_cell'] += .001
    with pytest.raises(ValueError, match='common native well'): m.require_same_well(original, wrong)


def test_missing_original_targets_are_not_filled_from_a_reference(tmp_path, monkeypatch):
    output = tmp_path/'analysis'
    # An empty namespace exercises missing-file behavior without moving or
    # deleting the real terminal archives or inventing material observations.
    m = reader(); monkeypatch.setattr(m, 'CASE', tmp_path/'missing_case')
    with pytest.raises(FileNotFoundError, match='await original terminal archive'):
        m.analyze(output)
    assert not output.exists()


def test_partial_summary_reports_only_actual_completed_strains(reference):
    rows = reader().summarize(reference['bands'], require_all=False)
    assert len(rows) == 1 and rows[0]['strain'] == 0.
    assert rows[0]['new_T_minus_preserving_meV_fu'] == pytest.approx(1.2503429411)
    with pytest.raises(KeyError): reader().summarize(reference['bands'])
    duplicate = reference['bands']+[copy.deepcopy(reference['bands'][0])]
    with pytest.raises(ValueError, match='duplicate'): reader().summarize(duplicate, require_all=False)
