"""Analytic work/curvature checks for the eight-point pilot."""

import numpy as np
from ase import Atoms

from scripts.prepare_hfo2_channel_work_probes import probe_directions, write_structure
from scripts.score_hfo2_channel_work_probes import score_pairs


def test_directions_remove_translation_and_separate_cell():
    center = Atoms("Hf4O8", cell=np.eye(3) * 5.2, pbc=True)
    previous, following = center.copy(), center.copy()
    following.positions[4, 2] += 0.05
    following.positions += 0.01
    directions = probe_directions(center, previous, following)
    np.testing.assert_allclose(directions @ directions.T, np.eye(2), atol=1e-12)
    np.testing.assert_allclose(directions[0, :36].reshape(12, 3).sum(axis=0), 0, atol=1e-12)
    assert np.count_nonzero(directions[1, :36]) == 0


def test_exact_cubic_diagnostic_shows_finite_step_bias():
    # E(t)=E0+g*t+.5*k*t²+a*t³. Energy gradient has a*h² bias.
    manifest = {"center_energy_eV_cell": -100, "predicted_center_directional_gradients_eV_A": [0.2, -0.1], "points": []}
    results = {}
    for direction in (0, 1):
        g, k, a = manifest["predicted_center_directional_gradients_eV_A"][direction], -0.5, 2.0
        for h in (0.01, 0.02):
            for sign in (-1, 1):
                i, t = len(manifest["points"]), sign * h
                manifest["points"].append({"index": i, "direction": direction, "step_A": h, "sign": sign})
                results[i] = {"energy_eV_cell": -100 + g * t + 0.5 * k * t * t + a * t ** 3,
                              "directional_gradient_eV_A": g + k * t + 3 * a * t ** 2}
    pairs = score_pairs(manifest, results)
    assert len(pairs) == 4
    for pair in pairs:
        np.testing.assert_allclose(pair["gradient_difference_eV_A"], 2 * pair["step_A"] ** 2, atol=1e-11)
        np.testing.assert_allclose(pair["energy_curvature_eV_A2"], -0.5, atol=1e-8)
        np.testing.assert_allclose(pair["gradient_curvature_eV_A2"], -0.5, atol=1e-11)


def test_missing_pair_is_not_filled_by_interpolation():
    manifest = {"points": [{"direction": 0, "step_A": 0.01, "sign": -1, "index": 0}]}
    assert score_pairs(manifest, {0: {"energy_eV_cell": 0}}) == []


def test_writer_uses_full_precision_fixed_species_and_orbital_names(tmp_path):
    atoms = Atoms("Hf4O8", cell=np.eye(3) * 5.123456789012345, pbc=True)
    atoms.set_scaled_positions(np.linspace(0.05, 0.9, 36).reshape(12, 3))
    write_structure(tmp_path / "STRU", atoms)
    body = (tmp_path / "STRU").read_text()
    assert "100Ry_4s2p2d1f.orb" in body and "100Ry_2s2p1d.orb" in body
    assert "5.1234567890123452" in body or "5.123456789012345" in body
    assert "Direct" in body
