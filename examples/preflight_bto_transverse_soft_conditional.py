"""Preflight one BTO variable-cell conditional point on the physical soft/soft plane.

This creates no calculator and submits no job. The cubic Gamma mode source is
the archived 1x1x1 phonon calculation, not a 4x4x4 phonon supercell. The latter
number is the electronic k mesh of the unchanged 100 Ry / 10 au DZP protocol.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.bto_q1q2_reference import (
    bto_transverse_soft_plane, load_bto_q1q2_reference, sha256,
    validate_bto_gamma_source,
)
from vcneb.mode_surface import ModePlane
from vcneb.phonons import GammaModes


def remaining_soft_y_direction(modes: GammaModes, plane: ModePlane) -> np.ndarray:
    """Orient the third unstable cubic mode with the Ti_y - Ba_y anchor."""

    if modes.n_atoms != 5 or plane.reference.size != 3 * modes.n_atoms + 6:
        raise ValueError("expected five-atom BTO with six open strain coordinates")
    anchor = np.zeros(3 * modes.n_atoms)
    anchor[4], anchor[1] = 1.0, -1.0  # Ti_y - Ba_y in Ba,Ti,O,O,O order
    subspace = modes.eigenvectors[:, :3]
    mass_root = np.repeat(np.sqrt(modes.masses_amu), 3)
    coefficients = subspace.T @ (mass_root * anchor)
    fixed = plane.axis_weights[:3, :]
    coefficients -= fixed @ (fixed.T @ coefficients)
    norm = float(np.linalg.norm(coefficients))
    if norm < 1e-8:
        raise ValueError("Ti_y - Ba_y has no independent projection in the soft triplet")
    direction = np.zeros(plane.reference.size)
    direction[:3 * modes.n_atoms] = modes.cartesian_eigenvectors()[:, :3] @ (coefficients / norm)
    metric_norm = float(np.sqrt(direction @ (plane.metric_weights * direction)))
    if (abs(metric_norm - 1.0) > 1e-8
            or np.max(np.abs(plane.axis_vectors.T @ (plane.metric_weights * direction))) > 1e-8):
        raise ValueError("third soft direction is not metric-unit and orthogonal to fixed Q")
    return direction


def _minimum_distance(atoms) -> float:
    distances = atoms.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    return float(np.min(distances))


def build_preflight(args: argparse.Namespace) -> dict:
    q = np.array([args.q_parallel, args.q_transverse], dtype=float)
    if (not np.all(np.isfinite(q)) or not np.isfinite(args.seed_amplitude)
            or args.seed_amplitude <= 0.0 or not np.isfinite(args.minimum_distance_A)
            or args.minimum_distance_A <= 0.0):
        raise ValueError("Q, seed amplitude and geometry bound must be finite and valid")
    gamma_hashes = validate_bto_gamma_source(
        args.gamma_provenance, args.force_sets, args.eigenpairs_provenance,
    )
    loaded = load_bto_q1q2_reference(
        args.report, args.reference, args.force_constants, args.phonopy_eigenpairs,
    )
    chart, modes = loaded.chart, loaded.modes
    cell = chart.reference_cell
    total_mass = float(np.sum(modes.masses_amu))
    strain_weights = total_mass * float(np.linalg.det(cell) ** (2.0 / 3.0)) * np.array(
        [1.0, 1.0, 1.0, 2.0, 2.0, 2.0]
    )
    plane = bto_transverse_soft_plane(
        loaded, strain_metric_weights_amu_A2=strain_weights,
    )
    if plane.reference.size != chart.coordinate_count:
        raise ValueError("mode plane and reference-cell chart disagree")
    gauge = chart.rigid_translation_directions()
    q_gauge_overlap = plane.axis_vectors.T @ (plane.metric_weights[:, None] * gauge)
    if np.max(np.abs(q_gauge_overlap)) > 1e-8:
        raise ValueError("Q axes overlap the rigid-translation gauge")
    third = remaining_soft_y_direction(modes, plane)
    if np.max(np.abs(gauge.T @ (plane.metric_weights * third))) > 1e-8:
        raise ValueError("third soft direction overlaps rigid translation")
    frozen = plane.frozen_coordinates(q)
    starts = [frozen, frozen + args.seed_amplitude * third,
              frozen - args.seed_amplitude * third]
    seed_rows = []
    for label, start in zip(("frozen", "+Q_y", "-Q_y"), starts):
        atoms = chart.to_atoms(start)
        distance = _minimum_distance(atoms)
        if not np.allclose(plane.project(start), q, rtol=0.0, atol=1e-10):
            raise ValueError("one branch seed changes fixed Q coordinates")
        if distance < args.minimum_distance_A or atoms.get_volume() <= 0.0:
            raise ValueError(f"unsafe conditional starting geometry: {label}")
        seed_rows.append({
            "label": label,
            "coordinates_u_A_eta_voigt": start.tolist(),
            "minimum_atomic_distance_A": distance,
            "volume_A3": float(atoms.get_volume()),
        })
    grid = json.loads(args.grid_result_manifest.read_text(encoding="utf-8"))
    parameters = grid.get("calculator_parameters", {})
    if (grid.get("status") != "converged" or grid.get("mpi_ranks") != 32
            or parameters.get("calculation") != "scf"
            or parameters.get("basis_type") != "lcao"
            or parameters.get("dft_functional") != "pbe"
            or parameters.get("ecutwfc") != 100.0
            or parameters.get("kpts") != [4, 4, 4]
            or parameters.get("scf_thr") != 1e-8
            or "Orb-DZP-10au" not in parameters.get("basis_dir", "")
            or not all("10au_100Ry" in name for name in parameters.get("basis", {}).values())
            or len(parameters.get("basis", {})) != 3
            or not grid.get("abacus_binary_sha256") or not grid.get("asset_sha256")):
        raise ValueError("archived grid does not establish the unchanged BTO 100Ry/10au contract")
    cubic = [point for point in grid.get("points", [])
             if abs(point.get("q1", np.inf)) < 1e-12 and abs(point.get("q2", np.inf)) < 1e-12]
    if len(cubic) != 1 or not np.isfinite(cubic[0].get("energy_eV", np.nan)):
        raise ValueError("archived grid lacks one finite cubic energy reference")
    identity = {"calculator_parameters": parameters, "asset_sha256": grid["asset_sha256"],
                "abacus_binary_sha256": grid["abacus_binary_sha256"]}
    calculator_hash = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    sources = {
        **loaded.source_hashes, **gamma_hashes,
        "grid_result_manifest": sha256(args.grid_result_manifest),
    }
    return {
        "kind": "bto_transverse_soft_conditional_preflight_no_dft",
        "axis_kind": "cubic_Gamma_unstable_triplet_Ti_z_Ba_z_and_Ti_x_Ba_x_anchored",
        "q_parallel_q_transverse_sqrt_amu_A": q.tolist(),
        "reference_id": plane.reference_id,
        "source_sha256": sources,
        "calculator_id": f"abacus-bto-pbe100-dzp10au-4x4x4-sha256:{calculator_hash}",
        "phonon_supercell": [1, 1, 1],
        "electronic_kpoints": [4, 4, 4],
        "strain_metric_weights_amu_A2": strain_weights.tolist(),
        "strain_metric_rule": "M_total * V_C^(2/3) * (1,1,1,2,2,2); optimizer coordinates only",
        "n_total_coordinates": chart.coordinate_count,
        "n_fixed_q_axes": 2,
        "n_fixed_rigid_translations": 3,
        "n_relaxed_orthogonal_coordinates": chart.coordinate_count - 5,
        "remaining_soft_y_metric_unit_direction": third.tolist(),
        "seed_amplitude_sqrt_amu_A": args.seed_amplitude,
        "minimum_allowed_atomic_distance_A": args.minimum_distance_A,
        "branch_starts": seed_rows,
        "cubic_reference_energy_eV": float(cubic[0]["energy_eV"]),
        "no_dft_launched": True,
        "not_a_conditional_PES_or_T_to_C_barrier": True,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("report", "reference", "force-constants", "phonopy-eigenpairs",
                 "gamma-provenance", "force-sets", "eigenpairs-provenance",
                 "grid-result-manifest", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--q-parallel", type=float, required=True)
    parser.add_argument("--q-transverse", type=float, required=True)
    parser.add_argument("--seed-amplitude", type=float, default=0.4)
    parser.add_argument("--minimum-distance-A", type=float, default=1.6)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite existing preflight: {args.output}")
    report = build_preflight(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "q": report["q_parallel_q_transverse_sqrt_amu_A"],
                      "minimum_seed_distance_A": min(row["minimum_atomic_distance_A"] for row in report["branch_starts"]),
                      "calculator_id": report["calculator_id"], "no_dft_launched": True}, indent=2))


if __name__ == "__main__":
    main()
