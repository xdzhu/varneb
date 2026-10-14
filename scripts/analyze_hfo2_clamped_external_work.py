"""Raw-stress external work on an audited training band, not a barrier slope.

Reuses the existing biaxial chart and nine portable native SCFs. The two
reported partial derivatives hold different open-row coordinates fixed;
their difference must not be interpreted as a relaxed barrier prediction.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.analyze_hfo2_clamped_reference import pinned_json, terminal_images
from scripts.audit_hfo2_static_replica import sha256
from vcneb import BiaxialClampedCurvatureCoordinates, clamped_plane_vcneb_boundary
from vcneb.epitaxial_boundary import cell_work_derivative

ROOT = Path(__file__).resolve().parents[1]
SEEDS = ROOT / "benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/clamped_seed_manifest.json"
SEEDS_SHA = "46c19e628b71e0a44a595be5f0ec9a2518e0719db0d1b286b294f5160f4af9a4"
ADMITTED_TRAINING_AUDITS = {
    # Actual completed zero-strain M, T, both flips and +0.01 T/M/preserving receipts;
    # never arbitrary caller digests or the reserved +0.005 prediction condition.
    "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c",
    "c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e",
    "632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41",
    "32a6e49a618e6e62632ae71a5e0ed093095254f6858221be2ade3d83cc194a28",
    "46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202",
    "b86cdb72eab9e3d8ff2125e16e6fbb2d0954b34d54c0f05099a581c066108389",
    "0c8d91532bc9051decdfd600eae12ce409fd84a86c677158077948fa352c9d8a",
}


def image_work(image, boundary, *, strain, cell_scale_A):
    """Frozen fractional atoms; row2 fixed versus its in-plane affine motion."""
    boundary.validate_images([image])
    chart = BiaxialClampedCurvatureCoordinates(
        image, boundary=boundary, cell_scale_A=cell_scale_A,
        anchor_strain=strain, strain_scale_A=cell_scale_A,
    )
    cell = image.cell.array.copy()
    stress = image.get_stress(voigt=False, apply_constraint=False)
    forces = image.get_forces(apply_constraint=False)
    gradient = chart.enthalpy_gradient(np.zeros(chart.size), image, forces, stress, 0.)
    fixed_open_row_jacobian = np.zeros((3, 3))
    fixed_open_row_jacobian[:2] = cell[:2] / (1 + strain)
    exact = cell_work_derivative(stress, cell, fixed_open_row_jacobian)
    if not np.isclose(gradient[-1] * chart.strain_scale_A, exact, atol=1e-11, rtol=0):
        raise ValueError("controlled chart and raw-stress Jacobian disagree")
    projector = np.eye(3) - np.outer(boundary.normal, boundary.normal)
    affine_jacobian = cell @ projector / (1 + strain)
    affine = cell_work_derivative(stress, cell, affine_jacobian)
    planar = image.get_volume() * np.sum(stress * projector) / (1 + strain)
    open_row_work = cell_work_derivative(stress, cell, affine_jacobian - fixed_open_row_jacobian)
    if (not np.isclose(planar, affine, atol=1e-11, rtol=0)
            or not np.isclose(affine - exact, open_row_work, atol=1e-11, rtol=0)):
        raise ValueError("affine reparameterization work identity failed")
    return dict(
        fixed_open_row_partial_eV_cell_per_strain=exact,
        affine_open_row_partial_eV_cell_per_strain=affine,
        affine_minus_fixed_open_row_eV_cell_per_strain=open_row_work,
        fixed_open_row_cell_jacobian_A_per_strain=fixed_open_row_jacobian.tolist(),
        affine_cell_jacobian_A_per_strain=affine_jacobian.tolist(),
        raw_stress_ASE_eV_A3=stress.tolist(),
        volume_A3=image.get_volume(),
    )


def analyze(case, audit_sha256):
    # Register actual completed training receipts before opening their labels.
    # Boundary rejection after reading a held-out receipt would be too late.
    if audit_sha256 not in ADMITTED_TRAINING_AUDITS:
        raise ValueError("an admitted completed training audit is required before label access")
    images, observation, residual, receipt = terminal_images(case, audit_sha256)
    seeds = pinned_json(SEEDS, SEEDS_SHA)
    if (seeds["strain_conditions"] != [0., .01]
            or seeds["holdout_strain_reserved_not_generated"] != .005
            or observation["pressure_GPa"] != 0.):
        raise ValueError("only the registered P=0 training study is admitted")
    reference = np.asarray(observation["mechanical_boundary"]["reference_cell_A"])
    unstrained = np.asarray(seeds["unstrained_T_reference_cell_A"])
    matches = [e for e in seeds["strain_conditions"]
               if np.allclose(reference[:2], (1 + e) * unstrained[:2], atol=1e-9, rtol=0)]
    if len(matches) != 1 or observation["mechanical_boundary"]["allow_tilt"] is not True:
        raise ValueError("a registered training substrate with open tilts is required")
    strain = matches[0]
    boundary = clamped_plane_vcneb_boundary(len(images[0]), reference, allow_tilt=True)
    boundary.validate_images(images)
    rows = [dict(image_index=i, **image_work(image, boundary, strain=strain,
                  cell_scale_A=observation["cell_scale_A"])) for i, image in enumerate(images)]
    nfu = observation["formula_units"]
    for row in rows:
        for prefix in ("fixed_open_row", "affine_open_row"):
            value = row[f"{prefix}_partial_eV_cell_per_strain"]
            row[f"{prefix}_partial_meV_fu_per_unit_strain"] = value * 1000 / nfu
            row[f"{prefix}_relative_PO_partial_meV_fu_per_unit_strain"] = (
                value - rows[0][f"{prefix}_partial_eV_cell_per_strain"]) * 1000 / nfu
    peak = int(np.argmax(observation["relative_enthalpy_meV_fu"]))
    return dict(
        status="native_training_band_frozen_external_work_not_stationary_barrier_response",
        source_job_id=receipt["job_id"], terminal_step=receipt["terminal_step"],
        terminal_audit_sha256=audit_sha256, training_seed_manifest_sha256=SEEDS_SHA,
        script_sha256=sha256(Path(__file__)),
        core_source_sha256={name: sha256(ROOT/name) for name in
                            ("vcneb/biaxial_curvature.py", "vcneb/epitaxial_boundary.py")},
        strain_training_condition=strain, pressure_GPa=0., formula_units=nfu,
        native_portable_frames_reparsed=len(images),
        complete_six_physical_bytes_rechecked_locally=False,
        ordinary_fmax_eV_A=residual["fmax_eV_A"], sampled_peak_image=peak,
        controlled_chart="first two row vectors scale with1+epsilon; open row2 fixed; fractional atoms fixed",
        comparison_chart="homogeneous in-plane affine motion also carries row2 tilt; fractional atoms fixed",
        sampled_peak=rows[peak], images=rows,
        new_DFT_calls=0, scheduler_mutations=0,
        geometry_parameters_optimizer_or_running_source_changed=False,
        G3_selection_or_holdout_access=False,
        limitations=[
            "raw substrate reaction stress is retained, not an internal optimization residual",
            "two partial derivatives differ by actual open-row work; different held internal coordinates",
            "their equality at a stationary released branch is not established by an ordinary NEB pass",
            "relative-PO partial work is not a derivative of a relaxed MEP maximum or a prospective prediction",
            "no finite-strain energy differences, Hessian index, uncertainty bound or H1/H2 pass measured",
        ],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--audit-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("fresh external-work output required")
    result = analyze(args.case, args.audit_sha256)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(job=result["source_job_id"], peak=result["sampled_peak_image"],
        peak_work=result["sampled_peak"], new_DFT_calls=0)))


if __name__ == "__main__":
    main()
