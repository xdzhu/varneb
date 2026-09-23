"""Check BaTiO3 T/C endpoint identity before an expensive VC-NEB run.

The ordered five-atom conventional perovskite cells used by the BTO examples
must contain a genuinely tetragonal, polar T endpoint and a nearly cubic,
non-polar C endpoint.  This geometry-only check is calculator-independent; it
does not replace a force/stress gate or a phonon stability calculation.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import endpoint_structure_record  # noqa: E402


def _wrapped_delta(value: np.ndarray) -> np.ndarray:
    return (value + 0.5) % 1.0 - 0.5


def bto_phase_metrics(atoms) -> dict:
    """Return cell anisotropy and Ti displacement from the equatorial O plane."""

    if Counter(atoms.get_chemical_symbols()) != Counter({"Ba": 1, "Ti": 1, "O": 3}):
        raise ValueError("BTO phase check requires exactly one Ba, one Ti, and three O atoms")
    a, b, c = (float(value) for value in atoms.cell.lengths())
    angles = [float(value) for value in atoms.cell.angles()]
    if max(abs(value - 90.0) for value in angles) > 1.0:
        raise ValueError("BTO T/C case requires an approximately orthogonal conventional cell")
    fractional = atoms.get_scaled_positions(wrap=True)
    symbols = atoms.get_chemical_symbols()
    titanium = fractional[symbols.index("Ti")]
    oxygen = fractional[[index for index, symbol in enumerate(symbols) if symbol == "O"]]
    offsets = _wrapped_delta(oxygen - titanium)
    # The two equatorial oxygen sites have the largest in-plane offset from Ti.
    equatorial = np.argsort(np.linalg.norm(offsets[:, :2], axis=1))[-2:]
    polar_offset = abs(float(np.mean(offsets[equatorial, 2]))) * c
    return {
        "n_atoms": len(atoms),
        "c_over_mean_ab": c / ((a + b) / 2.0),
        "ab_relative_mismatch": abs(a - b) / ((a + b) / 2.0),
        "ti_equatorial_o_z_offset_A": polar_offset,
        "cell_lengths_A": [a, b, c],
        "endpoint": endpoint_structure_record(atoms),
    }


def audit_bto_phase_pair(
    tetragonal,
    cubic,
    *,
    minimum_tetragonality: float = 0.005,
    minimum_t_polar_offset_A: float = 0.01,
    maximum_cubic_anisotropy: float = 0.002,
    maximum_c_polar_offset_A: float = 0.01,
) -> dict:
    """Validate declared BTO phase identities without evaluating a calculator."""

    if min(
        minimum_tetragonality,
        minimum_t_polar_offset_A,
        maximum_cubic_anisotropy,
        maximum_c_polar_offset_A,
    ) <= 0.0:
        raise ValueError("phase identity thresholds must be positive")
    t = bto_phase_metrics(tetragonal)
    c = bto_phase_metrics(cubic)
    issues = []
    if t["ab_relative_mismatch"] > maximum_cubic_anisotropy:
        issues.append("T endpoint a and b are not equivalent within the declared tolerance")
    if t["c_over_mean_ab"] < 1.0 + minimum_tetragonality:
        issues.append("T endpoint has no declared tetragonal c/a distortion")
    if t["ti_equatorial_o_z_offset_A"] < minimum_t_polar_offset_A:
        issues.append("T endpoint has no declared Ti-O polar displacement")
    if c["ab_relative_mismatch"] > maximum_cubic_anisotropy or abs(c["c_over_mean_ab"] - 1.0) > maximum_cubic_anisotropy:
        issues.append("C endpoint is not cubic within the declared tolerance")
    if c["ti_equatorial_o_z_offset_A"] > maximum_c_polar_offset_A:
        issues.append("C endpoint retains a Ti-O polar displacement")
    return {
        "status": "passed" if not issues else "failed",
        "tetragonal": t,
        "cubic": c,
        "thresholds": {
            "minimum_tetragonality": minimum_tetragonality,
            "minimum_t_polar_offset_A": minimum_t_polar_offset_A,
            "maximum_cubic_anisotropy": maximum_cubic_anisotropy,
            "maximum_c_polar_offset_A": maximum_c_polar_offset_A,
        },
        "issues": issues,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tetragonal", required=True, help="ASE-readable T endpoint")
    parser.add_argument("--cubic", required=True, help="ASE-readable C endpoint")
    parser.add_argument("--output", help="optional JSON audit path")
    parser.add_argument("--minimum-tetragonality", type=float, default=0.005)
    parser.add_argument("--minimum-t-polar-offset-A", type=float, default=0.01)
    parser.add_argument("--maximum-cubic-anisotropy", type=float, default=0.002)
    parser.add_argument("--maximum-c-polar-offset-A", type=float, default=0.01)
    args = parser.parse_args()
    report = audit_bto_phase_pair(
        read(args.tetragonal),
        read(args.cubic),
        minimum_tetragonality=args.minimum_tetragonality,
        minimum_t_polar_offset_A=args.minimum_t_polar_offset_A,
        maximum_cubic_anisotropy=args.maximum_cubic_anisotropy,
        maximum_c_polar_offset_A=args.maximum_c_polar_offset_A,
    )
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(payload, encoding="utf-8")
    print(payload, end="")
    if report["issues"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
