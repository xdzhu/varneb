"""Prepare one registered G2 channel from audited G1 geometry and clamped wells.

No DFT or submission. Original ordered endpoint lifts are fixed against the
registered *unrelaxed* seeds, never inferred by fitting a new relaxed endpoint.
The free-cell polygon is a starting mechanism, not a clamped MEP or TS proof.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.audit_hfo2_G1_gate import EXPECTED
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.prepare_hfo2_clamped_endpoints import geometry_guard, orient_long_axis_x, _json
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from scripts.relax_clamped_ase_endpoint import load_seed
from vcneb.periodic_path import minimum_image_path_lift, validate_periodic_path_lift


CHANNEL_FINAL = {"PO_to_T": "T", "PO_to_M": "M",
                 "PO_flip_T_pattern_preserving": "PO_minus_T_preserving",
                 "PO_flip_T_pattern_reversing": "PO_minus_T_reversing"}
HF_BASE = "/public/home/iai806/abacus/agent-runs/"
HF_PHYSICS = HF_BASE + ("20260912-vcneb-hf/validation/hfo2_t_to_po/"
                       "vcneb_n7_fire_distributed_resume27678004_v4/00")
ZERO_RUNTIME = {
    "PO_plus": "20261009-varneb-clamped-PO-E046-r1/endpoint",
    "T": "20261009-varneb-clamped-T-E047-r1/endpoint",
    "M": "20261009-varneb-clamped-pair-E048-r1/M/endpoint",
    "PO_minus_T_preserving": "20261009-varneb-clamped-pair-E048-r1/PO_minus_T_preserving/endpoint",
    "PO_minus_T_reversing": "20261009-varneb-clamped-reversing-minus-E049-r1/endpoint",
}


def endpoint_evidence(case, condition, phase):
    seed_path = case / "clamped_endpoint_seeds" / condition / phase / "endpoint_seed.json"
    seed, boundary, manifest = load_seed(seed_path)
    if manifest["physical_contract_sha256"] != CONTRACT or manifest["pressure_gpa"] != 0:
        raise ValueError("registered P=0 physical contract required")
    folder = case / "clamped_endpoint_matrix_20261009" / condition / phase
    if condition == "strain_0000" and phase == "PO_plus":
        folder = case / "clamped_PO_continuation_E046_20261009"
        audit_path = folder / "completed_HF/endpoint/continuation_audit.json"
    elif condition == "strain_0000" and phase == "T":
        audit_path = folder / "completed_HF/T_endpoint_audit.json"
    else:
        audit_path = folder / "raw_audit_HF.json"
    root = folder / "completed_HF/endpoint"
    audit = json.loads(audit_path.read_text())
    summary_path = root / "endpoint_relax_summary.json"
    summary = json.loads(summary_path.read_text())
    if (not audit["endpoint_converged"] or not audit["raw_full_physical_bytes_checked_here"]
            or audit["summary_sha256"] != sha256(summary_path)
            or not summary["converged"] or summary["status"] != "completed"
            or summary["physical_contract_sha256"] != CONTRACT
            or summary["strain"] != manifest["strain"] or summary["external_pressure_gpa"] != 0
            or not summary["allow_tilt"] or summary["fmax_target_eV_per_A"] != .03
            or summary["open_stress_target_kbar"] != 2. or summary["maxstep"] != .02
            or summary["max_atomic_force_eV_per_A"] >= .03 or summary["open_traction_norm_kbar"] >= 2.):
        raise ValueError("raw-audited same-boundary converged endpoint required")
    calls = sorted((root / "calculator/image_0000").glob("scf_*"))
    terminal = calls[-1]
    point = audit["rows"][-1]
    call = json.loads((terminal / "call_audit.json").read_text())
    hashes = call["input_sha256"]
    if (any(hashes[name] != h for name, h in CONTRACT.items())
            or any(sha256(terminal/name) != hashes[name] for name in ("INPUT", "KPT", "STRU"))
            or hashes != point["input_sha256"]
            or sha256(terminal/"OUT.ABACUS/running_scf.log") != point["raw_log_sha256"]
            or call["raw_log_sha256"] != point["raw_log_sha256"]):
        raise ValueError("exported endpoint raw hashes changed")
    atoms = read_fixed_hfo2_stru(terminal / "STRU")
    raw = audited_results(terminal)
    if (not same_ordered_geometry(atoms, read(root/"CONTCAR", format="vasp"))
            or any(not np.allclose(raw[k], call["results"][k], atol=1e-12, rtol=0) for k in raw)
            or not np.isclose(raw["energy"], summary["potential_energy_eV"], atol=1e-10, rtol=0)):
        raise ValueError("terminal geometry or E/F/stress changed")
    boundary.validate_images([atoms])
    # The fixed calculator writes wrapped Direct STRU. Explicitly reconstruct
    # and record its short-step lift against the original registered seed.
    # For E046 include the canary before the geometry-only continuation reset.
    parent_calls = []
    if condition == "strain_0000" and phase == "PO_plus":
        parent = case / "G1_review_clamped_canary_20261009/completed_HF/endpoint/calculator/image_0000"
        parent_calls = sorted(parent.glob("scf_*"))
    history = [seed] + [read_fixed_hfo2_stru(p/"STRU") for p in parent_calls + calls]
    lifted_history, endpoint_lift = minimum_image_path_lift(history)
    validate_periodic_path_lift(lifted_history)
    atoms = lifted_history[-1]
    if condition == "strain_0000":
        remote = HF_BASE + ZERO_RUNTIME[phase]
    else:
        run = {"PO_plus": "20261010-varneb-clamped-p0100-E050-r1/PO_plus/endpoint",
               "T": "20261010-varneb-clamped-p0100-E050-r1/T/endpoint",
               "M": "20261010-varneb-clamped-p0100-M-E051-r1/endpoint"}.get(phase)
        if run is None:
            run = f"20261010-varneb-clamped-p0100-minus-E052-r1/{phase}/endpoint"
        remote = HF_BASE + run
    cache = {"directory": remote + "/calculator/image_0000/" + terminal.name,
             "input_sha256": hashes, "raw_log_sha256": point["raw_log_sha256"]}
    provenance = {"phase_label_seed": phase, "strain": manifest["strain"],
                  "seed_manifest_sha256": sha256(seed_path), "summary_sha256": sha256(summary_path),
                  "raw_audit_sha256": sha256(audit_path), "cache": cache,
                  "wrapped_raw_endpoint_short_step_lift": endpoint_lift,
                  "raw_full_physical_bytes_checked_on_HF": True,
                  "full_pseudo_or_basis_bytes_rechecked_offline": False}
    return seed, atoms, boundary, manifest, provenance


def deform_polygon(images, seeds, endpoints, boundary, n_images=9):
    """Affine substrate change + smooth endpoint correction, no remapping."""
    if n_images != 9:
        raise ValueError("registered G2 recipe is seven internal images plus two endpoints")
    source = [orient_long_axis_x(a) for a in images]
    validate_periodic_path_lift(source)
    for a in source:
        cell = a.cell.array.copy()
        cell[:2] = boundary.reference_cell[:2]
        a.set_cell(cell, scale_atoms=True)
    lifts, corrections, cell_corrections = [], [], []
    for original, seed, final in zip((source[0], source[-1]), seeds, endpoints):
        if not same_ordered_geometry(original, seed):
            raise ValueError("G1 ordered endpoint does not match registered unrelaxed seed")
        lift = np.rint(original.get_scaled_positions(wrap=False) - seed.get_scaled_positions(wrap=False)).astype(int)
        target = final.get_scaled_positions(wrap=False) + lift
        lifts.append(lift.tolist())
        corrections.append(target - original.get_scaled_positions(wrap=False))
        cell_corrections.append(final.cell.array - original.cell.array)
    result, segments = [], []
    for t in np.linspace(0., 1., n_images):
        x = t * (len(source) - 1)
        j = min(int(np.floor(x)), len(source) - 2)
        w = x - j
        q = (1-w)*source[j].get_scaled_positions(wrap=False) + w*source[j+1].get_scaled_positions(wrap=False)
        cell = (1-w)*source[j].cell.array + w*source[j+1].cell.array
        q += (1-t)*corrections[0] + t*corrections[1]
        cell += (1-t)*cell_corrections[0] + t*cell_corrections[1]
        cell[:2] = boundary.reference_cell[:2]
        a = source[j].copy()
        a.calc = None
        a.set_cell(cell, scale_atoms=False)
        a.set_scaled_positions(q)
        result.append(a)
        segments.append({"normalized_source_index": float(t), "source_segment": j, "segment_fraction": float(w)})
    boundary.validate_images(result)
    geometry_guard(result)
    validate_periodic_path_lift(result)
    if not all(same_ordered_geometry(a,b) for a,b in zip((result[0],result[-1]), endpoints)):
        raise ValueError("prepared endpoints changed ordered raw geometry")
    return result, {"fixed_endpoint_integer_lifts": lifts, "resampling": segments,
                    "source_total_images": len(source), "atom_permutation_applied": False,
                    "source_results_reused": False,
                    "policy": "source_index_polygon_with_affine_plane_and_linear_endpoint_corrections"}


def prepare(case, condition, channel, output):
    if output.exists():
        raise FileExistsError("refusing existing G2 seed namespace")
    if condition not in ("strain_0000", "strain_p0100") or channel not in CHANNEL_FINAL:
        raise ValueError("only registered training strains and four G2 channels allowed; no holdout")
    spec_path = case / "reversing_peak_sampling_20261009/network_specification.json"
    spec = json.loads(spec_path.read_text())
    entry = next(e for e in spec["channels"] if e["name"] == channel)
    source_path = (spec_path.parent / entry["observation"]).resolve()
    images, observation, digest = read_evaluated_observation(source_path)
    expected_digest = {name:digest for name,digest,_ in EXPECTED}[channel]
    if digest != expected_digest or observation["replayed_fmax_eV_A"] >= .10:
        raise ValueError("frozen G1 ordinary-converged source changed")
    if entry["reverse"]:
        images = list(reversed(images))
    left = endpoint_evidence(case, condition, "PO_plus")
    right = endpoint_evidence(case, condition, CHANNEL_FINAL[channel])
    if left[3]["cell_scale_A"] != right[3]["cell_scale_A"]:
        raise ValueError("one common cell metric required")
    left[2].validate_images([right[1]])
    result, lift = deform_polygon(images, (left[0],right[0]), (left[1],right[1]), left[2])
    output.mkdir(parents=True, exist_ok=False)
    write_clean_poscar(output/"initial.vasp", result[0])
    write_clean_poscar(output/"final.vasp", result[-1])
    write_clean_poscar(output/"substrate.vasp", left[0])
    write(output/"seed.traj", result)
    caches = [left[4]["cache"]] + [None]*7 + [right[4]["cache"]]
    _json(output/"factory_parameters.json", {"source_directory": HF_PHYSICS, "seed_cache_records": caches})
    manifest = {"format_version": 1, "status": "G2_geometry_prepared_DFT_pending", "channel": channel,
                "mechanical_family": "common_T_substrate_clamped_plane_P0_E0",
                "strain": left[3]["strain"], "holdout_generated": False, "pressure_GPa": 0.,
                "allow_tilt": True, "cell_scale_A": left[3]["cell_scale_A"],
                "reference_cell_A": left[3]["reference_cell_A"],
                "source_observation_sha256": digest, "source_reversed": entry["reverse"],
                "endpoint_evidence": [left[4],right[4]], "lift_preparation": lift,
                "n_total_images": 9, "n_internal_images": 7, "n_fixed_endpoints": 2,
                "fmax_eV_A": .10, "k_eV_A2": .2, "climb": False,
                "physical_inputs_changed": False, "physical_contract_sha256": CONTRACT,
                "cached_clamped_endpoint_evaluations": 2, "cached_interior_evaluations": 0,
                "new_DFT_calls": 0, "preparation_script_sha256": sha256(Path(__file__)),
                "files_sha256": {n: sha256(output/n) for n in ("initial.vasp","final.vasp","substrate.vasp","seed.traj","factory_parameters.json")},
                "limitations": ["clamped starting polygon is not the optimized G2 path",
                                "no TS, Hessian, electronic polarization or distinct winding certification",
                                "full physical bytes must be checked again on HF before using the two endpoint caches"]}
    _json(output/"manifest.json", manifest)
    return manifest


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case-root", type=Path, default=Path("benchmarks/hfo2_channels/20261008"))
    p.add_argument("--condition", choices=("strain_0000","strain_p0100"), required=True)
    p.add_argument("--channel", choices=tuple(CHANNEL_FINAL), required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    m = prepare(a.case_root, a.condition, a.channel, a.output)
    print(json.dumps({k:m[k] for k in ("status","channel","strain","n_total_images","new_DFT_calls")}))
