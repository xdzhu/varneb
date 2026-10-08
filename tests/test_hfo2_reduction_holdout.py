import json
from pathlib import Path

import numpy as np
import pytest

from scripts.prepare_hfo2_reduction_holdout import holdout_vectors
from scripts.score_hfo2_reduction_holdout import score_pairs


def test_response_keeps_Qx_amplitude_without_renormalization():
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"
    model = json.loads((root / "T_atomic_reduction.json").read_text())
    vectors = holdout_vectors(model)
    q = np.array(model["retained_basis_cartesian"])
    assert np.linalg.norm(vectors["responded"]) > 1
    np.testing.assert_allclose(q.T @ vectors["responded"], [1, 0, 0], atol=1e-12)
    model["families"]["d0.01"]["full_atomic_response"][0][0] = float("nan")
    with pytest.raises(ValueError):
        holdout_vectors(model)


def test_preregistered_gate_and_missing_pair():
    manifest = {"prediction_fixed_before_DFT": {"frozen": 4., "responded": 2.},
                "relative_curvature_acceptance": .1, "T_energy_eV_cell": -10., "points": []}
    results = {}
    for kind, k in manifest["prediction_fixed_before_DFT"].items():
        for h in (.05, .10):
            for sign in (-1, 1):
                i = len(manifest["points"])
                manifest["points"].append({"index": i, "kind": kind, "amplitude_A": h, "sign": sign})
                x = h * sign
                # Deliberately nonzero reference gradient cancels in curvatures.
                results[i] = {"energy_eV_cell": -10 + .3*x + .5*k*x*x,
                              "directional_gradient_eV_A": .3 + k*x}
    pairs = score_pairs(manifest, results)
    assert len(pairs) == 4 and all(p["passes_preregistered_10percent_gate"] for p in pairs)
    results[7]["directional_gradient_eV_A"] += 1
    assert not score_pairs(manifest, results)[-1]["passes_preregistered_10percent_gate"]
    del results[7]
    assert len(score_pairs(manifest, results)) == 3
