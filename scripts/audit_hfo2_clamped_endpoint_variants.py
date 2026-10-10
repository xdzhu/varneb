"""Read-only ordered-gauge review of the five registered zero-strain endpoints.

Parent-site translation/permutation equivalence is not ordered cache identity
or proof that two entire paths are equivalent. No atom or input is rewritten.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import read

from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.analyze_hfo2_network_update import structure_audit
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.prepare_hfo2_clamped_endpoints import AXES, geometry_guard, orient_long_axis_x
from scripts.relax_clamped_ase_endpoint import load_seed
from vcneb.reference_variants import apply_parent_operation


PHASES = ("T", "PO_plus", "M", "PO_minus_T_preserving", "PO_minus_T_reversing")


def declared_translation_pair(left, right, translation, permutation, *, tolerance_A=1e-8):
    """Test one declared same-frame operation, not a best-fit remapping.

    permutation[i] is the target index for source i. The periodic lift is an
    analysis record only and must not be substituted into a production chain.
    """
    geometry_guard([left, right])
    if (len(left) != 12 or left.get_chemical_symbols() != ["Hf"]*4+["O"]*8
            or right.get_chemical_symbols() != left.get_chemical_symbols()):
        raise ValueError("requires the same ordered Hf4O8 species")
    t = np.asarray(translation, dtype=float)
    p = np.asarray(permutation)
    if (t.shape != (3,) or not np.isfinite(t).all() or p.shape != (12,)
            or not np.issubdtype(p.dtype, np.integer) or sorted(p.tolist()) != list(range(12))
            or not np.array_equal(left.numbers, right.numbers[p])
            or not np.isfinite(tolerance_A) or tolerance_A <= 0):
        raise ValueError("finite translation, species-preserving bijection and positive tolerance required")
    f, g = left.get_scaled_positions(wrap=False), right.get_scaled_positions(wrap=False)
    ordered_vectors, ordered_distances = find_mic((g-f)@left.cell.array, left.cell, pbc=True)
    difference = g[p]-f-t
    vectors, distances = find_mic(difference@left.cell.array, left.cell, pbc=True)
    # Recover the lattice lift from the true MIC, not component-wise rounding
    # that can give the wrong image in a skewed cell.
    lift_float = difference-vectors@np.linalg.inv(left.cell.array)
    lift = np.rint(lift_float).astype(int)
    if not np.allclose(lift_float, lift, atol=1e-9, rtol=0):
        raise ValueError("minimum-image displacement lacks a consistent integer lift")
    cell_difference = float(np.abs(left.cell.array-right.cell.array).max())
    return {
        "declared_fractional_translation": t.tolist(), "source_to_target": p.tolist(),
        "target_periodic_lift_relative_to_translated_source": lift.tolist(),
        "cell_max_absolute_difference_A": cell_difference,
        "ordered_periodic_max_distance_A": float(ordered_distances.max()),
        "declared_operation_max_distance_A": float(distances.max()),
        "geometry_tolerance_A": tolerance_A,
        "same_ordered_periodic_geometry": same_ordered_geometry(left, right),
        "same_periodic_geometry_under_declared_translation_and_permutation": bool(
            cell_difference < tolerance_A and distances.max() < tolerance_A),
        "production_atoms_wrapped_translated_permuted_or_modified": False,
        "whole_path_equivalence_topology_or_electronic_polarization_certified": False,
    }


def geometric_patterns(atoms, parent, basis):
    delta = atoms.get_scaled_positions(wrap=False)-parent.get_scaled_positions(wrap=False)
    delta -= np.rint(delta)
    margin = float(.5-np.abs(delta).max())
    if margin <= .05:
        return {"nearest_reference_chart_accepted": False, "chart_margin_fractional": margin,
                "Q_A": None, "reason": "outside the registered local quarter-site chart"}
    vector = delta@parent.cell.array
    vector -= vector.mean(axis=0)
    q = np.einsum("aij,ij->a", basis, vector)
    residual = vector-np.einsum("a,aij->ij", q, basis)
    return {"nearest_reference_chart_accepted": True, "chart_margin_fractional": margin,
            "Q_A": q.tolist(), "orthogonal_residual_norm_A": float(np.linalg.norm(residual)),
            "oxygen_minus_hafnium_mean_displacement_A": (vector[4:].mean(0)-vector[:4].mean(0)).tolist(),
            "electronic_polarization_or_energy_decomposition": False}


def audit(case_root):
    matrix = case_root/"clamped_endpoint_matrix_20261009/strain_0000"
    t_source = case_root/"reference_variants/T.vasp"
    variant_path = case_root/"reference_variants/variant_manifest.json"
    variant = json.loads(variant_path.read_text())
    seeds_path = case_root/"clamped_endpoint_seeds/clamped_seed_manifest.json"
    seeds = json.loads(seeds_path.read_text())
    registered = {row["phase_label"]:row["manifest_sha256"] for row in seeds["seeds"]
                  if row["strain"] == 0.}
    if set(registered) != set(PHASES) or seeds["physical_contract_sha256"] != CONTRACT:
        raise ValueError("original finite common-substrate seed matrix required")
    if sha256(t_source) != variant["generated_structure_sha256"]["T.vasp"]:
        raise ValueError("registered T parent changed")
    t = read(t_source, format="vasp")
    parent_old, basis_old, _ = rotated_t_triplet(t)
    parent = orient_long_axis_x(parent_old)
    basis = basis_old[..., list(AXES)]
    # Composition of the two previously registered inversion operations gives
    # this parent translation. It is fixed before inspecting endpoint fits.
    tp = np.asarray(variant["variants"]["T_preserving_inversion"]["operation"]["translation_fractional"])
    tr = np.asarray(variant["variants"]["T_reversing_inversion"]["operation"]["translation_fractional"])
    translation = np.mod((tr-tp)[list(AXES)], 1.)
    operation = apply_parent_operation(parent, parent, np.eye(3), translation)
    records, initial, terminal, results = {}, {}, {}, {}
    for phase in PHASES:
        seed = case_root/"clamped_endpoint_seeds/strain_0000"/phase/"endpoint_seed.json"
        if sha256(seed) != registered[phase]:
            raise ValueError("registered endpoint seed manifest changed")
        atoms, boundary, manifest = load_seed(seed)
        endpoint = (case_root/"clamped_PO_continuation_E046_20261009/completed_HF/endpoint"
                    if phase == "PO_plus" else matrix/phase/"completed_HF/endpoint")
        summary_path = endpoint/"endpoint_relax_summary.json"
        summary = json.loads(summary_path.read_text())
        if (manifest["strain"] != 0. or manifest["pressure_gpa"] != 0.
                or not manifest["allow_tilt"] or manifest["physical_contract_sha256"] != CONTRACT
                or summary["strain"] != 0.
                or summary["external_pressure_gpa"] != 0. or not summary["allow_tilt"]
                or summary["physical_contract_sha256"] != CONTRACT
                or summary["fmax_target_eV_per_A"] != .03 or summary["open_stress_target_kbar"] != 2.
                or not summary["converged"] or summary["max_atomic_force_eV_per_A"] >= .03
                or summary["open_traction_norm_kbar"] >= 2.):
            raise ValueError("same-boundary physically screened zero-strain endpoint required")
        calls = sorted((endpoint/"calculator/image_0000").glob("scf_*"))
        static = calls[-1]
        call = json.loads((static/"call_audit.json").read_text())
        if (any(sha256(static/name) != CONTRACT[name] for name in ("INPUT","KPT"))
                or any(call["input_sha256"][name] != digest for name,digest in CONTRACT.items())
                or sha256(static/"STRU") != call["input_sha256"]["STRU"]
                or sha256(static/"OUT.ABACUS/running_scf.log") != call["raw_log_sha256"]):
            raise ValueError("terminal raw input/log provenance changed")
        final = read_fixed_hfo2_stru(static/"STRU")
        boundary.validate_images([final])
        if not same_ordered_geometry(final, read(endpoint/"CONTCAR", format="vasp")):
            raise ValueError("terminal ordered geometry changed")
        raw = audited_results(static)
        if any(not np.allclose(raw[k],call["results"][k],atol=1e-12,rtol=0) for k in raw):
            raise ValueError("terminal E/F/stress changed")
        initial[phase], terminal[phase], results[phase] = atoms, final, raw
        records[phase] = {"seed_sha256": sha256(seed), "summary_sha256": sha256(summary_path),
                          "terminal_STRU_sha256": sha256(static/"STRU"), "raw_log_sha256": call["raw_log_sha256"],
                          "energy_eV_Hf4O8": raw["energy"], "physical_endpoint_screen_passed": True,
                          "structure_audit": structure_audit(final),
                          "registered_seed_patterns": geometric_patterns(atoms, parent, basis),
                          "terminal_patterns": geometric_patterns(final, parent, basis)}
    preserving, reversing = "PO_minus_T_preserving", "PO_minus_T_reversing"
    pairs = {name: declared_translation_pair(values[preserving], values[reversing],
                                             translation, operation.source_to_target)
             for name,values in (("registered_seed",initial),("relaxed_endpoint",terminal))}
    a, b = results[preserving], results[reversing]
    p = np.asarray(operation.source_to_target)
    covariance = {"terminal_energy_difference_eV_cell": b["energy"]-a["energy"],
                  "terminal_force_covariance_max_absolute_error_eV_A": float(np.abs(b["forces"][p]-a["forces"]).max()),
                  "terminal_stress_covariance_max_absolute_error_eV_A3": float(np.abs(b["stress"]-a["stress"]).max()),
                  "evaluations_compared": "two independent actual zero-strain endpoint SCFs",
                  "not_a_general_error_bound_or_plus_one_percent_cache_certificate": True}
    return {"status": "five_zero_strain_endpoint_representations_reviewed_not_G2_barriers",
            "new_DFT_calls": 0, "analysis_script_sha256": sha256(Path(__file__)),
            "variant_manifest_sha256": sha256(variant_path),
            "finite_seed_matrix_sha256": sha256(seeds_path),
            "terminal_INPUT_KPT_bytes_checked_here": True,
            "licensed_pseudo_orbital_bytes_checked_here": False,
            "software": {name: importlib.metadata.version(name) for name in ("ase","pymatgen","spglib")},
            "fixed_reference": "original T quarter-site scaffold, proper axes(2,0,1), original atom order",
            "pattern_basis_order": ["old_x_geometric_T_pattern","old_y_geometric_T_pattern","old_z_geometric_T_pattern"],
            "pattern_metric": "translation-free Cartesian Angstrom in oriented original T lattice; no phonon/irrep claim",
            "substrate_strain": 0., "pressure_GPa": 0., "holdout_generated_or_read": False,
            "endpoints": records, "minus_pair": pairs, "raw_covariance": covariance,
            "production_source_runtime_or_geometry_modified": False,
            "limitations": ["parent translation equivalence is not identity at ordered path indices",
                            "retain both registered path lifts/correspondences; do not merge paths or count different winding sectors",
                            "geometric mean displacement is not electronic polarization",
                            "phase labels/force screens are not full Hessian stability or TS certification",
                            "no new physical variant, plus-one endpoint cache eligibility, G2 barrier or holdout prediction certified"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", type=Path, default=Path("benchmarks/hfo2_channels/20261008"))
    args = parser.parse_args()
    print(json.dumps(audit(args.case_root), sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
