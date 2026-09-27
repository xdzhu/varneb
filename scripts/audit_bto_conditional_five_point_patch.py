"""Recompute the BTO five-point conditional-patch geometry screen offline.

This compact snapshot is downstream of separate raw ABACUS audits. It does not
certify a continuous conditional PES, a Hessian minimum, or a VCNEB barrier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


EXPECTED_Q = {
    "q060_q000": (0.6, 0.0),
    "q090_q000": (0.9, 0.0),
    "q060_q030": (0.6, 0.3),
    "q090_q030": (0.9, 0.3),
    "q075_q015_holdout": (0.75, 0.15),
}
MIRRORED_LABELS = {"q090_q000", "q090_q030"}


def _vector(point: dict, key: str) -> list[float]:
    values = point.get(key)
    if (not isinstance(values, list) or len(values) != 21
            or any(not isinstance(x, (int, float)) or not math.isfinite(x) for x in values)):
        raise ValueError(f"{point['label']}: malformed {key}")
    return [float(x) for x in values]


def _weighted_square(a: list[float], b: list[float], weights: list[float]) -> float:
    return math.fsum(weight * (x - y) ** 2 for weight, x, y in zip(weights, a, b))


def audit_patch(source: Path) -> dict:
    """Check declared mirror evidence and report measured-patch interpolation.

    The Q axes are assumed metric-orthonormal as independently established by
    the cited preflight. Subtracting the known Q step removes its prescribed
    part from each neighboring structural difference.
    """

    raw = source.read_bytes()
    data = json.loads(raw)
    contract = data.get("calculator_contract", {})
    if (contract.get("ecutwfc_Ry") != 100
            or contract.get("electronic_kpoints") != [4, 4, 4]
            or contract.get("gamma_force_constant_supercell") != [1, 1, 1]
            or "10 au DZP" not in contract.get("orbitals", "")):
        raise ValueError("BTO fixed calculator/phonon contract changed")
    weights = data.get("metric_weights", [])
    if (len(weights) != 21 or any(not isinstance(w, (int, float))
                                  or not math.isfinite(w) or w <= 0 for w in weights)):
        raise ValueError("expected 21 positive finite atomic-plus-strain metric weights")
    weights = [float(w) for w in weights]
    mirror_indices = data.get("mirror_y_negated_coordinate_indices_zero_based")
    if mirror_indices != [1, 4, 7, 10, 13, 18, 20]:
        raise ValueError("the declared BTO y-mirror action changed")
    points = data.get("points", [])
    by_label = {point["label"]: point for point in points}
    if len(points) != len(EXPECTED_Q) or set(by_label) != set(EXPECTED_Q):
        raise ValueError("expected exactly four corners and the preregistered holdout")

    selected: dict[str, list[float]] = {}
    aligned: dict[str, list[float]] = {}
    energies: dict[str, float] = {}
    qy: dict[str, float] = {}
    mirror_checks: dict[str, dict] = {}
    for label, q in EXPECTED_Q.items():
        point = by_label[label]
        actual_q = point.get("q", [])
        if len(actual_q) != 2 or any(abs(a - b) > 1e-12 for a, b in zip(actual_q, q)):
            raise ValueError(f"{label}: wrong fixed-Q location")
        digest = point.get("source_sha256", "")
        if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            raise ValueError(f"{label}: missing source-file digest")
        energy = float(point["selected_energy_eV_per_BTO"])
        selected_qy = float(point["selected_qy_sqrt_amu_A"])
        if not math.isfinite(energy) or not math.isfinite(selected_qy):
            raise ValueError(f"{label}: nonfinite energy or third-mode amplitude")
        coords = _vector(point, "selected_coordinates_u_A_eta_voigt")
        selected[label] = coords
        energies[label] = energy
        if label in MIRRORED_LABELS:
            reflected = [(-x if i in mirror_indices else x)
                         for i, x in enumerate(coords)]
            alternative = _vector(point, "aligned_positive_qy_coordinates_u_A_eta_voigt")
            max_coordinate_difference = max(abs(x - y)
                                            for x, y in zip(reflected, alternative))
            alternative_energy = float(point["aligned_positive_qy_energy_eV_per_BTO"])
            alternative_qy = float(point["aligned_positive_qy_sqrt_amu_A"])
            energy_difference = abs(energy - alternative_energy)
            if (selected_qy >= 0 or alternative_qy <= 0
                    or max_coordinate_difference > 1e-8
                    or energy_difference > 1e-8
                    or abs(selected_qy + alternative_qy) > 1e-8):
                raise ValueError(f"{label}: alternative is not verified mirror-equivalent")
            aligned[label] = alternative
            qy[label] = alternative_qy
            mirror_checks[label] = {
                "max_coordinate_difference_mixed_units": max_coordinate_difference,
                "energy_difference_eV_per_BTO": energy_difference,
            }
        else:
            if selected_qy <= 0 or "aligned_positive_qy_coordinates_u_A_eta_voigt" in point:
                raise ValueError(f"{label}: expected the directly selected positive branch")
            aligned[label] = coords
            qy[label] = selected_qy

    corner_labels = ("q060_q000", "q090_q000", "q060_q030", "q090_q030")
    holdout_label = "q075_q015_holdout"
    prediction_energy = math.fsum(energies[label] for label in corner_labels) / 4
    prediction_qy = math.fsum(qy[label] for label in corner_labels) / 4
    prediction_coordinates = [math.fsum(aligned[label][i] for label in corner_labels) / 4
                              for i in range(21)]
    holdout = selected[holdout_label]
    atomic_square = _weighted_square(holdout[:15], prediction_coordinates[:15], weights[:15])
    strain_square = _weighted_square(holdout[15:], prediction_coordinates[15:], weights[15:])

    edges = (("q060_q000", "q090_q000"), ("q060_q030", "q090_q030"),
             ("q060_q000", "q060_q030"), ("q090_q000", "q090_q030"))

    def neighbor_jumps(coordinates: dict[str, list[float]]) -> list[float]:
        jumps = []
        for left, right in edges:
            q_left, q_right = EXPECTED_Q[left], EXPECTED_Q[right]
            prescribed_square = math.fsum((a - b) ** 2 for a, b in zip(q_left, q_right))
            off_plane_square = _weighted_square(
                coordinates[left], coordinates[right], weights,
            ) - prescribed_square
            if off_plane_square < -1e-8:
                raise ValueError("coordinates and declared fixed-Q metric are inconsistent")
            jumps.append(math.sqrt(max(0.0, off_plane_square)))
        return jumps

    raw_prediction = [math.fsum(selected[label][i] for label in corner_labels) / 4
                      for i in range(21)]
    return {
        "status": "five_point_metrics_recomputed_not_conditional_PES_certificate",
        "source_snapshot_sha256": hashlib.sha256(raw).hexdigest(),
        "corner_order": list(corner_labels),
        "holdout_q": list(EXPECTED_Q[holdout_label]),
        "predicted_holdout_energy_eV_per_BTO": prediction_energy,
        "measured_holdout_energy_eV_per_BTO": energies[holdout_label],
        "measured_minus_predicted_energy_meV_per_BTO": 1000 * (
            energies[holdout_label] - prediction_energy),
        "predicted_positive_branch_Qy_sqrt_amu_A": prediction_qy,
        "measured_positive_branch_Qy_sqrt_amu_A": qy[holdout_label],
        "atomic_coordinate_error_sqrt_amu_A": math.sqrt(atomic_square),
        "strain_coordinate_error_sqrt_amu_A": math.sqrt(strain_square),
        "full_coordinate_error_sqrt_amu_A": math.sqrt(atomic_square + strain_square),
        "raw_selected_holdout_coordinate_error_sqrt_amu_A": math.sqrt(
            _weighted_square(holdout, raw_prediction, weights)),
        "aligned_neighbor_off_plane_jumps_sqrt_amu_A": neighbor_jumps(aligned),
        "raw_selected_neighbor_off_plane_jumps_sqrt_amu_A": neighbor_jumps(selected),
        "mirror_checks": mirror_checks,
        "limitations": [
            "Snapshot digests locate original reports; this script does not re-audit remote raw DFT.",
            "No material-specific interpolation thresholds or rigorous Hessian uncertainty are supplied.",
            "One cell and one holdout cannot certify a continuous or global conditional PES.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="compact five-point source snapshot")
    args = parser.parse_args()
    print(json.dumps(audit_patch(args.source), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
