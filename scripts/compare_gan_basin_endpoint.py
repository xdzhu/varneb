"""Compare a pressure-relaxed GaN basin candidate with a mapped VCNEB endpoint.

This is a structural identity check, not an energy comparison or full
thermodynamic phase certificate. Both cells must have the same atom order.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import spglib
from ase.io import read

from scripts.audit_gan_ts_basin_pilot import ga_n_coordination
from scripts.prepare_gan_ts_newton_probe import sha256


def relative_displacement_metrics(candidate, reference) -> dict:
    if (len(candidate) != len(reference)
            or candidate.get_chemical_symbols() != reference.get_chemical_symbols()
            or len(candidate) != 4):
        raise ValueError("candidate/reference must be mapped Ga2N2 four-atom cells")
    difference = candidate.get_scaled_positions() - reference.get_scaled_positions()
    common_translation = difference[0].copy()
    difference -= common_translation
    difference -= np.rint(difference)
    cartesian = difference @ reference.cell.array
    norms = np.linalg.norm(cartesian, axis=1)
    return {
        "common_fractional_translation_mod_1": np.mod(common_translation, 1.0).tolist(),
        "relative_displacement_rms_A": float(np.sqrt(np.mean(norms**2))),
        "relative_displacement_max_A": float(np.max(norms)),
        "relative_volume_difference": float(candidate.get_volume() / reference.get_volume() - 1.0),
    }


def spacegroup_symbol(atoms, symprec_A: float) -> str:
    dataset = spglib.get_symmetry_dataset(
        (atoms.cell.array, atoms.get_scaled_positions(), atoms.get_atomic_numbers()),
        symprec=symprec_A,
    )
    if dataset is None:
        raise ValueError("spglib could not classify the mapped GaN cell")
    return str(dataset.international)


def compare(candidate_path: Path, source_trajectory: Path, output: Path,
            endpoint_index: int = -1) -> dict:
    if output.exists():
        raise FileExistsError(output)
    if endpoint_index not in (0, -1):
        raise ValueError("only original VCNEB B4/B1 endpoints are allowed")
    candidate = read(candidate_path, format="vasp")
    reference = read(source_trajectory, index=endpoint_index)
    metrics = relative_displacement_metrics(candidate, reference)
    symprec = 0.01
    result = {
        "status": "GaN_basin_candidate_endpoint_structural_comparison_only",
        "endpoint": "B4" if endpoint_index == 0 else "B1",
        "same_order_species": candidate.get_chemical_symbols(),
        "candidate_coordination_GaN_2p4A": ga_n_coordination(candidate),
        "reference_coordination_GaN_2p4A": ga_n_coordination(reference),
        "candidate_spacegroup_at_0p01A": spacegroup_symbol(candidate, symprec),
        "reference_spacegroup_at_0p01A": spacegroup_symbol(reference, symprec),
        "candidate_volume_A3": float(candidate.get_volume()),
        "reference_volume_A3": float(reference.get_volume()),
        **metrics,
        "source_sha256": {"candidate": sha256(candidate_path),
                          "trajectory": sha256(source_trajectory),
                          "analyzer": sha256(Path(__file__))},
        "limitations": [
            "A common translation is removed; atom permutation is not guessed.",
            "A 0.01 Å symmetry classification and 2.4 Å coordination screen do not prove stress convergence.",
            "The original 600 eV VCNEB endpoint and this 1000 eV diagnostic have incomparable absolute enthalpies.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"endpoint": result["endpoint"],
                      "rms_A": metrics["relative_displacement_rms_A"],
                      "candidate_spacegroup": result["candidate_spacegroup_at_0p01A"]}))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--source-trajectory", type=Path, required=True)
    parser.add_argument("--endpoint-index", type=int, choices=[0, -1], default=-1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    compare(args.candidate, args.source_trajectory, args.output,
            args.endpoint_index)


if __name__ == "__main__":
    main()
