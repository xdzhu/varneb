"""Resolve native training gradients in all registered fixed-reference spans.

No new DFT, curvature, model selection, held-out labels or source-mode export.
Atomic covectors are transported dually to the original reference chart;
open-cell gradients and substrate reactions remain separate diagnostics.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.analyze_hfo2_clamped_external_work import ADMITTED_TRAINING_AUDITS
from scripts.analyze_hfo2_clamped_reference import (
    analyze as reference_analysis, original_chart_images, pinned_json, terminal_images,
)
from scripts.analyze_hfo2_cmma_path_reference import MASS, prepare_frame
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.audit_hfo2_static_replica import sha256
from vcneb.continuous_projection import project_reference_atomic_gradient
from vcneb.mode_subspaces import real_mode_subspace
from vcneb.qe_modes import read_qe_gamma_modes

ROOT = Path(__file__).resolve().parents[1]


def analyze(case, source_root, registration, variants, gamma, *, terminal_audit_sha256,
            clamped_seed_manifest, clamped_seed_manifest_sha256):
    if terminal_audit_sha256 not in ADMITTED_TRAINING_AUDITS:
        raise ValueError("admitted completed training audit required before label access")
    # Reuse the complete existing source/registration/boundary audit, not just
    # a caller-supplied coverage JSON. No new frame or mode is chosen here.
    coverage = reference_analysis(case, source_root, registration, variants, gamma,
        terminal_audit_sha256=terminal_audit_sha256,
        clamped_seed_manifest=clamped_seed_manifest,
        clamped_seed_manifest_sha256=clamped_seed_manifest_sha256)
    images, _, residual, receipt = terminal_images(case, terminal_audit_sha256)
    registered = pinned_json(registration, coverage["original_registration_sha256"])
    t = read(variants/"T.vasp", format="vasp")
    seeds = pinned_json(clamped_seed_manifest, clamped_seed_manifest_sha256)
    original = original_chart_images(images, seeds["orientation"], t)
    rotation = np.asarray(seeds["orientation"]["cartesian_rotation_rows"])
    _, patterns, _ = rotated_t_triplet(t)
    triplet = real_mode_subspace(patterns.reshape(3, 36).T*np.repeat(np.sqrt(MASS), 3)[:, None], expected_rank=3)
    bases = [("T_geometric_triplet", triplet.basis_columns)]
    with np.load(gamma, allow_pickle=False) as source:
        bases.append(("T_two_lowest_complete_optical_doublets", source["eigenvectors_mass_weighted"][:,
                      coverage["T_lowest_two_complete_doublets_rank4_indices"]]))
    modes = read_qe_gamma_modes(source_root/"matdyn.modes", MASS)
    cmma = read(source_root/"Cmma.vasp", format="vasp")
    for frame in registered["selected_frame_registrations"]:
        _, basis, _ = prepare_frame(modes, cmma, t,
            registered["author_source_audit"]["QE_to_published_POSCAR"], frame)
        bases.append((f"Cmma_frame_{frame['id']}", basis))
    rows = []
    for i, (native, image) in enumerate(zip(images, original)):
        # Only the covector array is transported; no incorrectly framed
        # calculator results are attached to the geometry-only copies.
        force = native.get_forces(apply_constraint=False) @ rotation
        spans = {name: project_reference_atomic_gradient(force, image.cell.array,
                    t.cell.array, basis, masses_amu=MASS) for name, basis in bases}
        rows.append(dict(image_index=i, relative_energy_meV_fu=coverage["relative_energy_meV_fu"][i],
                         raw_atomic_force_max_eV_A=float(np.linalg.norm(force, axis=1).max()),
                         representations=spans))
    peak = coverage["highest_image_index"]
    return dict(status="training_atomic_gradient_coverage_not_model_accuracy_or_TS",
        script_sha256=sha256(Path(__file__)),
        projection_source_sha256=sha256(ROOT/"vcneb/continuous_projection.py"),
        terminal_audit_sha256=terminal_audit_sha256,
        registration_sha256=coverage["original_registration_sha256"],
        author_source_sha256=coverage["author_source_sha256"],
        clamped_seed_manifest_sha256=clamped_seed_manifest_sha256,
        source_job_id=receipt["job_id"], terminal_step=receipt["terminal_step"],
        ordinary_fmax_eV_A=residual["fmax_eV_A"], strain_training_condition=coverage["strain_training_condition"],
        fixed_frame_transport=coverage["fixed_G2_orientation"],
        sampled_peak_image=peak, sampled_peak=rows[peak], images=rows,
        all_E023_Cmma_frames_retained=True, portable_native_frames_reparsed=9,
        complete_six_physical_bytes_rechecked_locally=False,
        new_DFT_calls=0, scheduler_mutations=0, G3_selection_or_holdout_access=False,
        physical_inputs_production_modes_geometry_or_optimizer_changed=False,
        limitations=[
            "atomic fixed-cell energy covectors only; open-cell gradients are not atomic modes",
            "translation-free gauge restriction is recorded, not an ASR repair of source modes/forces",
            "mass/reference-chart gradient fractions are not energy partitions or ordinary eV/A thresholds",
            "an omitted movable atomic direction is not a mechanical support reaction",
            "native force samples and ordinary convergence do not certify local curvature, TS index or a conditional branch",
            "ranks3/4, reference choices and all four Cmma frames are unchanged; no strong-model selection/accuracy claim",
        ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "source-root", "registration", "variants", "gamma", "clamped-seed-manifest", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--terminal-audit-sha256", required=True)
    parser.add_argument("--clamped-seed-manifest-sha256", required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("fresh gradient analysis output required")
    result = analyze(args.case, args.source_root, args.registration, args.variants, args.gamma,
        terminal_audit_sha256=args.terminal_audit_sha256,
        clamped_seed_manifest=args.clamped_seed_manifest,
        clamped_seed_manifest_sha256=args.clamped_seed_manifest_sha256)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write("\n")
    print(json.dumps(dict(job=result["source_job_id"], peak=result["sampled_peak"], new_DFT_calls=0)))


if __name__ == "__main__":
    main()
