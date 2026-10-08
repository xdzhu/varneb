"""Register published Cmma directions in audited HfO2 paths, without DFT.

Reference orientation is selected from endpoint geometry, never mode coverage,
energy or an unseen strained path. All geometry-equivalent choices are retained.
Projected spans are descriptive and do not establish a predictive advantage.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from ase.io import read
import numpy as np
from scipy.optimize import linear_sum_assignment

from scripts.audit_hfo2_cmma_reference import AUTHOR_COMMIT, audit as audit_author
from scripts.audit_hfo2_static_replica import sha256
from scripts.analyze_hfo2_chain_observations import analyze_observation
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from vcneb.continuous_projection import continuous_reference_coordinates, project_reference_basis
from vcneb.mode_subspaces import real_mode_subspace, mass_translation_overlaps
from vcneb.phonons import identify_acoustic_modes
from vcneb.qe_modes import read_qe_gamma_modes


T_AUTHOR_BLOB = "b8cb27541fa7d1fe166b37fa9c432d31cf664047"
T_AUTHOR_SHA256 = "afcb23ec69f99d4ef9fa546795222d4c28eac4249659eb35c528fe5f82131a53"
MASS = np.array([178.49] * 4 + [15.9994] * 8)
COST_TIE_A = 1e-7


def proper_signed_frames():
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            rotation = np.eye(3, dtype=int)[list(permutation)] * np.array(signs)[:, None]
            if round(np.linalg.det(rotation)) == 1:
                yield rotation


def endpoint_assignment(source, target, rotation, origin, metric_cell, *, second_best=False):
    """Nearest-site *registration*, not same-phase/equivalent-cell matching.

    The orthogonal original-T metric measures all candidates identically.
    Source/target cells may differ physically; that is not silently fitted away.
    Only the endpoint chooses a species assignment. No path image is remapped.
    """
    metric = np.asarray(metric_cell, float)
    rotation, origin = np.asarray(rotation, float), np.asarray(origin, float)
    if (metric.shape != (3, 3) or not np.isfinite(metric).all()
            or np.linalg.det(metric) <= 0
            or not np.allclose(metric @ metric.T, np.diag(np.diag(metric @ metric.T)), atol=1e-10, rtol=0)
            or rotation.shape != (3, 3) or origin.shape != (3,)
            or not np.isfinite(np.r_[rotation.ravel(), origin]).all()
            or not np.allclose(rotation, np.rint(rotation), atol=1e-12, rtol=0)
            or not np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-12, rtol=0)
            or not np.isclose(np.linalg.det(rotation), 1, atol=1e-12, rtol=0)):
        raise ValueError("proper signed frame, finite origin and orthogonal positive metric required")
    if (not 2 <= len(source) <= 128 or len(source) != len(target)
            or not all(source.pbc) or not all(target.pbc)
            or sorted(source.numbers) != sorted(target.numbers)):
        raise ValueError("compatible bounded periodic compositions required")
    for atoms in (source, target):
        if (not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell.array).all()
                or np.linalg.det(atoms.cell.array) <= 0):
            raise ValueError("finite right-handed source and endpoint required")
    q = source.get_scaled_positions(wrap=False) @ rotation.T + origin
    delta = q[:, None, :] - target.get_scaled_positions(wrap=False)[None, :, :]
    integers = np.rint(delta)
    distances = np.linalg.norm((delta - integers) @ metric, axis=2)
    distances[source.numbers[:, None] != target.numbers[None, :]] = np.inf
    row, permutation = linear_sum_assignment(distances)
    matched = distances[row, permutation]
    result = {"cost_sum_A": float(matched.sum()), "maximum_site_offset_A": float(matched.max()),
              "source_to_target": permutation.tolist(), "site_offsets_A": matched.tolist(),
              "integer_shifts_source_minus_target": integers[row, permutation].astype(int).tolist()}
    if second_best:
        alternatives = []
        for i, j in enumerate(permutation):
            forbidden = distances.copy()
            forbidden[i, j] = np.inf
            try:
                r, p = linear_sum_assignment(forbidden)
            except ValueError:
                continue
            alternatives.append(float(forbidden[r, p].sum()))
        if not alternatives:
            raise ValueError("assignment ambiguity cannot be measured")
        result["second_best_assignment_cost_gap_A"] = min(alternatives) - result["cost_sum_A"]
    return result


def register_frames(author_t, cmma, tetragonal, po):
    records = []
    for rotation in proper_signed_frames():
        for origin in itertools.product((-.25, .25), repeat=3):
            t_fit = endpoint_assignment(author_t, tetragonal, rotation, origin, tetragonal.cell.array)
            po_fit = endpoint_assignment(cmma, po, rotation, origin, tetragonal.cell.array)
            records.append({"id": len(records), "rotation_and_integer_basis": rotation.tolist(),
                            "origin_fractional": list(origin), "T_cost_sum_A": t_fit["cost_sum_A"],
                            "Cmma_to_PO_cost_sum_A": po_fit["cost_sum_A"]})
    if len(records) != 192:
        raise ValueError("bounded signed-frame/origin enumeration incomplete")
    t_best = min(r["T_cost_sum_A"] for r in records)
    compatible = [r for r in records if r["T_cost_sum_A"] <= t_best + COST_TIE_A]
    po_best = min(r["Cmma_to_PO_cost_sum_A"] for r in compatible)
    selected = [r.copy() for r in compatible if r["Cmma_to_PO_cost_sum_A"] <= po_best + COST_TIE_A]
    for record in selected:
        assignment = endpoint_assignment(cmma, po, record["rotation_and_integer_basis"],
                                         record["origin_fractional"], tetragonal.cell.array, second_best=True)
        if assignment["second_best_assignment_cost_gap_A"] <= COST_TIE_A:
            raise ValueError("common PO anchor has ambiguous atom correspondence")
        record["registered_PO_assignment"] = assignment
    return records, selected, {"T_minimum_cost_sum_A": t_best,
        "T_compatible_frame_count": len(compatible), "PO_minimum_cost_sum_A": po_best,
        "selected_frame_count": len(selected), "cost_tie_tolerance_A": COST_TIE_A,
        "selection_data": "known G1 T and common PO+ geometry only, before G2/holdout labels",
        "selection_rule": "T-compatible proper frames first; closest Cmma-to-PO registration second; retain all ties",
        "not_used_for_selection": ["mode coverage", "energy", "forces", "G2", "holdout"]}


def prepare_frame(cm_modes, cmma, tetragonal, author_mapping, frame):
    rotation = np.asarray(frame["rotation_and_integer_basis"], int)
    origin = np.asarray(frame["origin_fractional"])
    assignment = np.asarray(frame["registered_PO_assignment"]["source_to_target"], int)
    qe_to_cm = np.asarray(author_mapping["source_to_target_indices"], int)
    qe_rotation = np.asarray(author_mapping["rotation_Cartesian_proper"])
    values = cm_modes.mass_weighted_eigenvectors[:, :4].reshape(12, 3, 4)
    ordered = np.empty_like(values)
    ordered[assignment[qe_to_cm]] = np.einsum("ab,nbk->nak", rotation @ qe_rotation, values)
    cell = rotation @ cmma.cell.array @ rotation.T
    # A declared common fractional-displacement chart, not a physical strain
    # applied to any source calculation or a new dynamical-matrix eigenbasis.
    pullback = np.linalg.solve(cell, tetragonal.cell.array)
    weighted_common = np.einsum("nck,cd->ndk", ordered, pullback).reshape(36, 4)
    span = real_mode_subspace(weighted_common, expected_rank=4)
    reference = cmma[np.argsort(assignment)]
    q = cmma.get_scaled_positions(wrap=False) @ rotation.T + origin
    reference.set_cell(tetragonal.cell.array)
    reference.set_scaled_positions(q[np.argsort(assignment)])
    reference.calc = None
    diagnostics = {
        "Cartesian_direction_to_common_T_chart": pullback.tolist(),
        "source_QE_to_production_indices": assignment[qe_to_cm].tolist(),
        "real_rank": span.expected_rank, "real_singular_values": span.real_singular_values.tolist(),
        "relative_rank_tolerance": span.relative_rank_tolerance,
        "maximum_relative_complex_projection_error": float(span.relative_column_projection_errors.max()),
        "orthogonality_max_defect": span.orthogonality_max_defect,
        "rigid_translation_overlap_subspace_Frobenius_norm": float(np.linalg.norm(mass_translation_overlaps(span.basis_columns, MASS))),
        "ASR_or_translation_removal_from_basis": False,
        "basis_coordinates": "real SVD axes, not named phonon eigenvectors; only total subspace weights compared",
    }
    return reference, span.basis_columns, diagnostics


def lowest_complete_doublet_indices(frequencies, optical):
    frequencies = np.asarray(frequencies, float)
    optical = np.asarray(optical)
    if (frequencies.ndim != 1 or not np.isfinite(frequencies).all()
            or optical.ndim != 1 or not np.issubdtype(optical.dtype, np.integer)
            or len(optical) < 4 or len(set(optical.tolist())) != len(optical)
            or np.any(optical < 0) or np.any(optical >= len(frequencies))):
        raise ValueError("complete finite optical frequency indices required")
    lowest = optical[np.argsort(frequencies[optical])[:4]]
    groups = [optical[np.isclose(frequencies[optical], frequencies[lowest[i]], atol=1e-6, rtol=0)]
              for i in (0, 2)]
    if (any(len(group) != 2 for group in groups)
            or len(set(groups[0].tolist()) | set(groups[1].tolist())) != 4):
        raise ValueError("registered four-mode baseline requires two complete distinct doublets")
    return lowest


def analyze(source_root, dataset, variants, gamma_path):
    author = audit_author(source_root)
    path = source_root / "Tetragonal.vasp"
    if path.stat().st_size != 919:
        raise ValueError("published T frame anchor byte count changed")
    data = path.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    if blob != T_AUTHOR_BLOB or hashlib.sha256(data).hexdigest() != T_AUTHOR_SHA256:
        raise ValueError("published T frame anchor changed")
    t, po, cmma, author_t = [read(p, format="vasp") for p in
        (variants / "T.vasp", variants / "PO.vasp", source_root / "Cmma.vasp", path)]
    if any(a.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8 for a in (t, po, cmma, author_t)):
        raise ValueError("ordered Hf4O8 references required")
    scores, selected, registration = register_frames(author_t, cmma, t, po)
    modes = read_qe_gamma_modes(source_root / "matdyn.modes", MASS)
    if not np.all(modes.frequencies_THz[:4] < -.1):
        raise ValueError("reviewed four negative reference candidates absent")
    parent, patterns, _ = rotated_t_triplet(t)
    weighted_patterns = patterns.reshape(3, 36).T * np.repeat(np.sqrt(MASS), 3)[:, None]
    t_span = real_mode_subspace(weighted_patterns, expected_rank=3)
    frames = [prepare_frame(modes, cmma, t, author["QE_to_published_POSCAR"], frame) for frame in selected]
    folders = sorted((p.parent for p in dataset.glob("*/observation.json")), key=lambda p: p.name.casefold())
    if not folders or len(folders) > 16:
        raise ValueError("one to sixteen complete frozen observations required")
    observations = []
    with np.load(gamma_path, allow_pickle=False) as gamma:
        acoustic, _ = identify_acoustic_modes(gamma["eigenvectors_mass_weighted"], MASS)
        optical = np.setdiff1d(np.arange(36), acoustic)
        lowest = lowest_complete_doublet_indices(gamma["frequencies_thz"], optical)
        low_basis = gamma["eigenvectors_mass_weighted"][:, lowest]
        for folder in folders:
            audit = analyze_observation(folder, t, parent, patterns, gamma)
            images = read(folder / "evaluated_chain.traj", index=":")
            chart = continuous_reference_coordinates(images, t)
            triplet = project_reference_basis(chart["displacements_A"], t_span.basis_columns, metric_weights=MASS)
            low = project_reference_basis(chart["displacements_A"], low_basis, metric_weights=MASS)
            frame_results = []
            for frame_record, (reference, basis, diagnostics) in zip(selected, frames):
                common = project_reference_basis(chart["displacements_A"], basis, metric_weights=MASS)
                native_chart = continuous_reference_coordinates(images, reference)
                difference = native_chart["displacements_A"] - chart["displacements_A"]
                if not np.allclose(difference, difference[0], atol=1e-10, rtol=0):
                    raise ValueError("Cmma affine origin is not constant across the continuous path")
                native = project_reference_basis(native_chart["displacements_A"], basis, metric_weights=MASS)
                frame_results.append({"frame_id": frame_record["id"], "basis_audit": diagnostics,
                    "fixed_initial_integer_gauge": native_chart["fixed_integer_lattice_shifts_by_atom"].tolist(),
                    "common_T_origin_fraction": common["captured_squared_norm_fraction"].tolist(),
                    "common_T_origin_residual_sqrt_amu_A": common["residual_norm"].tolist(),
                    "Cmma_affine_origin_residual_sqrt_amu_A": native["residual_norm"].tolist(),
                    "Cmma_affine_origin_fraction_not_same_denominator": native["captured_squared_norm_fraction"].tolist()})
            observations.append({"source": folder.name, "source_observation_SHA256": sha256(folder / "observation.json"),
                "source_trajectory_SHA256": sha256(folder / "evaluated_chain.traj"),
                "audited_images": len(images), "replayed_fmax_eV_A": audit["replayed_fmax_eV_A"],
                "highest_image_index": audit["highest_image_index"],
                "common_T_origin_zero_displacement": triplet["zero_displacement"].tolist(),
                "T_triplet_mass_metric_fraction_rank3": triplet["captured_squared_norm_fraction"].tolist(),
                "T_triplet_residual_sqrt_amu_A": triplet["residual_norm"].tolist(),
                "T_lowest_two_optical_doublets_fraction_rank4": low["captured_squared_norm_fraction"].tolist(),
                "T_lowest_two_optical_doublets_indices_zero_based": lowest.tolist(),
                "Cmma_reference_frames": frame_results})
    return {"schema_version": 1, "status": "descriptive_strong_reference_registration_not_prediction",
        "analysis_script_SHA256": sha256(Path(__file__)),
        "mode_subspaces_module_SHA256": sha256(Path(__file__).resolve().parents[1] / "vcneb/mode_subspaces.py"),
        "author_commit": AUTHOR_COMMIT, "author_source_audit": author,
        "additional_T_anchor": {"URL": f"https://raw.githubusercontent.com/yuboqiuab/unstableflatband/{AUTHOR_COMMIT}/figure1/Tetragonal.vasp",
                                "git_blob": blob, "bytes": len(data), "SHA256": T_AUTHOR_SHA256},
        "production_T_SHA256": sha256(variants / "T.vasp"), "production_PO_SHA256": sha256(variants / "PO.vasp"),
        "production_T_Gamma_SHA256": sha256(gamma_path),
        "frame_scores_geometry_only": scores, "registration": registration,
        "selected_frame_registrations": selected, "observations": observations,
        "new_DFT_calls": 0, "physical_parameters_changed": False,
        "metric": "original production T fractional-displacement chart; mass metric Hf178.49/O15.9994amu; sqrt(amu)A residuals",
        "limitations": ["G1-informed descriptive registration, not unseen-condition prediction or blind reference choice",
            "nearest-site registration under a declared chart is not proof of unique physical atom correspondence",
            "four Cmma negative directions and T doublets are not complete optical spaces; T-triplet dimension differs",
            "common-origin increment fractions and affine-origin fractions have different denominators and are not conflated",
            "SVD basis represents a declared geometric span; source translation admixture retained, no ASR/new eigenpairs",
            "literature LDA vectors transformed to a common chart are not PBE force constants or curvatures",
            "coverage is not energy partition, a conditional surface, a whole-variable TS or predictive accuracy",
            "raw and transformed author eigenvector arrays are not redistributed; source hashes and scalar evidence only"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source-root", "dataset", "variants", "gamma", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing existing Cmma path analysis")
    result = analyze(args.source_root, args.dataset, args.variants, args.gamma)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"registration": result["registration"],
                      "audited_observations": len(result["observations"]), "new_DFT_calls": 0}))


if __name__ == "__main__":
    main()
