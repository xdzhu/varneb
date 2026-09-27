"""Test whether two fixed 600-eV image-15 joint modes represent the GaN path.

This is a geometry-only gate, not a transition-state certificate. Both the
archived VCNEB chain and the local image-15 Hessian use the original 600-eV
electronic contract. No DFT job is launched.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.io.trajectory import Trajectory

from scripts.analyze_gan_ts_mode_path_overlap import joint_delta
from vcneb.joint_curvature import JointCurvatureCoordinates


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plane_coverage(displacements: np.ndarray, modes: np.ndarray) -> dict:
    """Return absolute and relative residuals for orthonormal joint modes."""

    y = np.asarray(displacements, dtype=float)
    e = np.asarray(modes, dtype=float)
    if (y.ndim != 2 or e.ndim != 2 or y.shape[1] != e.shape[0]
            or not 1 <= e.shape[1] <= y.shape[1]
            or not np.isfinite(y).all() or not np.isfinite(e).all()
            or not np.allclose(e.T @ e, np.eye(e.shape[1]), atol=1e-8, rtol=0)):
        raise ValueError("invalid displacements or nonorthonormal mode basis")
    coordinates = y @ e
    residual = y - coordinates @ e.T
    residual_norm = np.linalg.norm(residual, axis=1)
    total = float(np.sum(y * y))
    return {
        "coordinates_A": coordinates.tolist(),
        "off_plane_norm_A": residual_norm.tolist(),
        "maximum_off_plane_norm_A": float(np.max(residual_norm)),
        "endpoint_off_plane_norm_A": [float(residual_norm[0]), float(residual_norm[-1])],
        "captured_squared_path_fraction": (
            float(1.0 - np.sum(residual * residual) / total) if total > 0 else 1.0
        ),
    }


def analyze(trajectory: Path, chain_audit: Path, center_poscar: Path,
            hessian_manifest: Path, hessian_audit: Path, hessian_npz: Path) -> dict:
    chain = json.loads(chain_audit.read_text(encoding="utf-8"))
    manifest = json.loads(hessian_manifest.read_text(encoding="utf-8"))
    audit = json.loads(hessian_audit.read_text(encoding="utf-8"))
    vasp = chain.get("backends", {}).get("vasp", {})
    if (chain.get("status") != "five_hashed_final_trajectories_match_plotted_enthalpy_curves"
            or chain.get("n_images_total") != 29
            or vasp.get("compact_final_chain_sha256") != sha256(trajectory)
            or manifest.get("center_POSCAR_sha256") != sha256(center_poscar)
            or audit.get("source_sha256", {}).get("manifest") != sha256(hessian_manifest)
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("coordinate_count") != 18
            or audit.get("n_static_displacements") != 36
            or audit.get("status") != "GaN_image15_local_joint_curvature_at_nonstationary_NEB_candidate"
            or audit.get("center_gradient_translation_free_eV_per_A", 0) <= 0.05):
        raise ValueError("archived GaN path or Hessian provenance is not the audited case")
    expected_inputs = {
        "INCAR": "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772",
        "KPOINTS": "b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee",
        "POTCAR": "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168",
    }
    if (len(audit.get("cases", [])) != 36
            or any({key: case.get("input_sha256", {}).get(key) for key in expected_inputs}
                   != expected_inputs for case in audit["cases"])):
        raise ValueError("image-15 Hessian does not use the production 600-eV inputs")
    with Trajectory(trajectory) as source:
        if len(source) != 29:
            raise ValueError("expected the compact final 29-image GaN chain")
        frames = list(source)
    if any(frame.get_chemical_symbols() != ["Ga", "Ga", "N", "N"] for frame in frames):
        raise ValueError("GaN atom count, order, or species changed")
    center = read(center_poscar, format="vasp")
    scale = float(manifest["cell_scale_A"])
    chart = JointCurvatureCoordinates(center, cell_scale_A=scale)
    free = chart.translation_free_basis()
    with np.load(hessian_npz) as data:
        eig = np.asarray(data["eigenvalues"], dtype=float)
        modes = np.asarray(data["eigenvectors"], dtype=float)
    if (eig.shape != (15,) or modes.shape != (18, 15)
            or not np.allclose(eig, audit["translation_free_eigenvalues_eV_per_A2"], atol=1e-8)
            or not np.allclose(modes.T @ modes, np.eye(15), atol=1e-8)
            or not eig[0] < 0 < eig[1]):
        raise ValueError("local Hessian eigensystem no longer matches its audit")
    shifts = []
    rotations = []
    margins = []
    q0 = center.get_scaled_positions(wrap=False)
    for frame in frames:
        delta, rotation = joint_delta(center, frame, scale)
        shifts.append(free @ (free.T @ delta))
        rotations.append(rotation)
        dq = frame.get_scaled_positions(wrap=False) - q0
        dq -= np.rint(dq)
        margins.append(float(np.min(0.5 - np.abs(dq))))
    y = np.asarray(shifts)
    if np.linalg.norm(y[15]) > 1e-5:
        raise ValueError("Hessian center does not coincide with final image 15")
    unstable = modes[:, 0].copy()
    if np.dot(unstable, y[-1] - y[0]) < 0:
        unstable *= -1
    after_unstable = y - np.outer(y @ unstable, unstable)
    _, residual_singular_values, vt = np.linalg.svd(after_unstable, full_matrices=False)
    secondary = vt[0].copy()
    if residual_singular_values[0] <= 1e-10:
        raise ValueError("no independent path direction remains after unstable mode")
    secondary -= unstable * np.dot(secondary, unstable)
    secondary /= np.linalg.norm(secondary)
    if np.dot(secondary, y[-1] - y[0]) < 0:
        secondary *= -1
    _, singular_values, empirical_vt = np.linalg.svd(y, full_matrices=False)
    anchored = plane_coverage(y, np.column_stack((unstable, secondary)))
    local_eigenpair = plane_coverage(y, modes[:, :2])
    empirical = plane_coverage(y, empirical_vt[:2].T)
    peak_center_offset = float(np.linalg.norm(y[15]))
    max_rotation = float(np.max(rotations))
    min_margin = float(np.min(margins))
    # A global fixed plane is promoted only if both end-to-end geometric
    # coverage and the periodic/lattice chart pass predeclared margins.
    full_path_fixed_plane_viable = bool(
        anchored["maximum_off_plane_norm_A"] <= 0.05
        and anchored["captured_squared_path_fraction"] >= 0.99
        and min_margin >= 0.02 and max_rotation <= 1e-5
    )
    return {
        "kind": "GaN_45p7_600eV_image15_two_joint_mode_global_coverage_gate_geometry_only",
        "source_protocols": {
            "path": "VASP PBE 600 eV, Ga_d/N, 8x8x6, 29 images, 45.7 GPa",
            "local_modes": "VASP PBE 600 eV, Ga_d/N, 8x8x6, image-15 Hessian, 45.7 GPa",
        },
        "n_images": 29,
        "n_atoms": 4,
        "hessian_lowest_two_eV_per_A2": eig[:2].tolist(),
        "image15_to_center_joint_distance_A": peak_center_offset,
        "image15_translation_free_gradient_eV_per_A": audit["center_gradient_translation_free_eV_per_A"],
        "max_cell_rotation_residual_frobenius": max_rotation,
        "min_periodic_half_cell_margin_fractional": min_margin,
        "max_successive_joint_step_A": float(np.max(np.linalg.norm(np.diff(y, axis=0), axis=1))),
        "image15_negative_mode_plus_path_empirical_transverse": anchored,
        "image15_lowest_two_eigenmodes": local_eigenpair,
        "unconstrained_empirical_rank2": empirical,
        "secondary_axis_overlap_with_next_image15_mode": float(abs(np.dot(secondary, modes[:, 1]))),
        "path_joint_singular_values_A": singular_values.tolist(),
        "full_path_fixed_plane_viable_under_predeclared_gate": full_path_fixed_plane_viable,
        "predeclared_gate": {
            "max_off_plane_A": 0.05,
            "min_captured_squared_fraction": 0.99,
            "min_periodic_half_cell_margin_fractional": 0.02,
            "max_cell_rotation_residual_frobenius": 1e-5,
        },
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "chain_audit": sha256(chain_audit),
            "center_POSCAR": sha256(center_poscar),
            "hessian_manifest": sha256(hessian_manifest),
            "hessian_audit": sha256(hessian_audit),
            "hessian_npz": sha256(hessian_npz),
            "analysis_script": sha256(Path(__file__)),
        },
        "limitations": [
            "A geometric projection is not a computed two-dimensional enthalpy surface.",
            "Image 15 has a nonzero full gradient; its Hessian is not a TS certificate.",
            "The secondary anchored axis is an empirical combination, not a phonon eigenmode.",
            "The nearest-image coordinate chart must not be extrapolated through an ambiguous periodic mapping.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "chain-audit", "center-poscar", "hessian-manifest",
                 "hessian-audit", "hessian-npz", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = analyze(args.trajectory, args.chain_audit, args.center_poscar,
                     args.hessian_manifest, args.hessian_audit, args.hessian_npz)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "fixed_plane_viable": report["full_path_fixed_plane_viable_under_predeclared_gate"],
        "anchored_max_off_plane_A": report["image15_negative_mode_plus_path_empirical_transverse"]["maximum_off_plane_norm_A"],
        "anchored_captured_fraction": report["image15_negative_mode_plus_path_empirical_transverse"]["captured_squared_path_fraction"],
    }))


if __name__ == "__main__":
    main()
