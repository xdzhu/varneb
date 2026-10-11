"""Native-cache stationarity diagnostics in the original clamped active space."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.analyze_hfo2_clamped_reference import terminal_images
from scripts.audit_hfo2_static_replica import sha256
from vcneb import VCNEB, clamped_plane_vcneb_boundary

FIELDS = (
    "allowed_true_generalized_force_norm_eV_per_A",
    "allowed_true_generalized_force_max_vector_eV_per_A",
    "constraint_reaction_force_norm_eV_per_A",
    "constraint_reaction_force_max_vector_eV_per_A",
)


def analyze(case, audit_sha256):
    # This reparses the nine portable INPUT/KPT/STRU/logs and checks their
    # ordered E/F/stress, not just an already rendered numeric figure report.
    images, observation, residual, receipt = terminal_images(case, audit_sha256)
    reference = np.asarray(observation["mechanical_boundary"]["reference_cell_A"])
    boundary = clamped_plane_vcneb_boundary(12, reference, allow_tilt=True)
    chain = VCNEB(images, cell_scale=observation["cell_scale_A"], k=.2, climb=False,
                  pressure=0, **boundary.vcneb_kwargs(images))
    before = chain.get_forces().copy()
    saddle, path = chain.saddle_diagnostics(), chain.path_diagnostics()
    after = chain.get_forces()
    if not np.array_equal(before, after):
        raise ValueError("diagnostic replay changed the ordinary force")
    records = []
    for row in path["images"]:
        i = row["image_index"]
        item = dict(image_index=i, energy_relative_PO_meV_fu=observation["relative_enthalpy_meV_fu"][i],
            **{k: row[k] for k in FIELDS},
            true_tangential_force_eV_per_A=row["true_tangential_force_eV_per_A"],
            true_perpendicular_force_eV_per_A=row["true_perpendicular_force_eV_per_A"],
            ordinary_residual_max_vector_eV_per_A=row["neb_residual_generalized_force_eV_per_A"])
        # Legacy path tangents/perpendiculars already use the allowed force.
        # Legacy saddle norms did not, so preserve both meanings explicitly.
        if i == saddle["image_index"]:
            if (any(abs(row[k]-saddle[k]) > 1e-10 for k in FIELDS)
                    or abs(row["true_tangential_force_eV_per_A"]-
                           saddle["allowed_true_tangential_force_eV_per_A"]) > 1e-10
                    or abs(row["true_perpendicular_force_eV_per_A"]-
                           saddle["allowed_true_perpendicular_force_eV_per_A"]) > 1e-10):
                raise ValueError("saddle and path active-space diagnostics disagree")
        records.append(item)
    if abs(np.linalg.norm(before, axis=1).max()-residual["fmax_eV_A"]) > 1e-12:
        raise ValueError("ordinary residual differs from the native terminal audit")
    return dict(status="ordinary_clamped_candidate_gradient_audit_not_TS_certificate",
        source_job_id=receipt["job_id"], snapshot_step=receipt["terminal_step"],
        terminal_audit_sha256=audit_sha256,
        source_observation_sha256=sha256(case/"observations"/f"step_{receipt['terminal_step']:04d}"/"observation.json"),
        script_sha256=sha256(Path(__file__)), core_source_sha256=sha256(Path(__file__).parents[1]/"vcneb/core.py"),
        native_portable_frames_reparsed=9, full_six_physical_bytes_rechecked_locally=False,
        original_clamped_reference_cell_A=reference.tolist(), cell_scale_A=observation["cell_scale_A"],
        ordinary_fmax_eV_A=residual["fmax_eV_A"], sampled_peak=saddle, images=records,
        diagnostics_leave_ordinary_forces_byte_identical=True, new_DFT_calls=0,
        geometry_parameters_or_optimizer_changed=False, G3_selection_or_holdout_access=False,
        limitations=["fixed cell reactions are not open-space stationarity failures",
            "mode/subspace projection excludes prescribed coordinates only, not arbitrary unstable directions",
            "allowed generalized force norms use the registered cell metric, not raw stress units",
            "nonzero sampled tangent and discrete path curvature do not establish Hessian index or a continuous saddle",
            "no Hessian, stationary response derivative, G3 selection or forecast certified"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--audit-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("fresh analysis output required")
    result = analyze(args.case, args.audit_sha256)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"job": result["source_job_id"], "sampled_peak": result["sampled_peak"]["image_index"],
        "allowed_max_vector_eV_A": result["sampled_peak"]["allowed_true_generalized_force_max_vector_eV_per_A"],
        "new_DFT_calls": 0}))


if __name__ == "__main__":
    main()
