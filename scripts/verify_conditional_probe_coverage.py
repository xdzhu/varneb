"""Verify finite probe coverage on synthetic grids, with zero calculator calls.

This is an implementation example, not an HfO2 surface, prospectively selected
material bottleneck, continuous interpolation bound or improvement benchmark.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np

from vcneb import (
    ConditionalHoldout, ConditionalModeSurface, ConditionalPointEvidence,
    ConditionalProbePlan, ConditionalProbeRecord, ConditionalProbeSpec, ModePlane,
    audit_conditional_probe_coverage, register_conditional_probe_plan,
)


def write_new(path, data):
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n")


def verify(workdir):
    workdir = Path(workdir)
    if workdir.exists():
        raise FileExistsError("fresh synthetic verification namespace required")
    plane = ModePlane(reference=np.zeros(3), basis=np.diag([.5, 1., 1./3.]),
        metric_weights=np.array([4., 1., 9.]), axis_weights=np.eye(3)[:, :2],
        axis_labels=("Q1", "Q2"), amplitude_unit="sqrt(amu)*angstrom", reference_id="synthetic-reference")
    axes = np.linspace(0., 1., 9)
    x, y = np.meshgrid(axes, axes)
    surface = ConditionalModeSurface(q1=axes, q2=axes, energies=2*x + 3*y + .5*x*y,
        axis_labels=plane.axis_labels, amplitude_unit=plane.amplitude_unit,
        reference_id=plane.reference_id, energy_unit="eV/cell", potential_kind="E",
        boundary_condition="synthetic fixed-Q local branch; not material data",
        orthogonal_gradient_norms=np.zeros((9, 9)), n_evaluations=np.ones((9, 9), dtype=int),
        selected_starts=np.zeros((9, 9), dtype=int), coordinates=np.stack([x/2., y, .03*x + .02*y], axis=-1))
    evidence = ConditionalPointEvidence(calculator_id="synthetic-only", energy_reference_id="common-zero",
        raw_audit_sha256="a"*64, branch_audit_sha256="b"*64, curvature_audit_sha256="c"*64,
        curvature_resolution_sha256="d"*64, minimum_orthogonal_curvature=.3, curvature_uncertainty=.05)
    grid_evidence = [[evidence for _ in axes] for _ in axes]
    probes = [ConditionalProbeSpec(f"p{i}", i, i, ((axes[i]+axes[i+1])/2.,)*2) for i in range(8)]
    plan = register_conditional_probe_plan(plane, surface, grid_evidence, probes,
        training_evaluation_ids=[f"training-{i}" for i in range(81)],
        protocol_id="synthetic-example", start_policy_id="declared-common-starts",
        selection_rule_id="diagonal-cell-centers", probe_labels_unread=True,
        energy_zero_id="common-zero", gradient_tolerance=.01, maximum_neighbor_orthogonal_jump=.1,
        maximum_holdout_energy_error=.02, maximum_holdout_coordinate_error=.05)
    workdir.mkdir(parents=True, exist_ok=False)
    # Persist the training-only registration before constructing any test labels.
    write_new(workdir / "probe_manifest.json", plan.to_dict())
    plan = ConditionalProbePlan.from_dict(json.loads((workdir / "probe_manifest.json").read_text(encoding="utf-8")))
    observations = {}
    for probe in plan.probes:
        a, b = probe.q
        observations[probe.probe_id] = ConditionalProbeRecord(
            f"new-{probe.probe_id}", ConditionalHoldout(probe.cell_i, probe.cell_j, probe.q,
                2*a + 3*b + .5*a*b, np.array([a/2., b, .03*a + .02*b]), evidence, 0.))
    perfect = audit_conditional_probe_coverage(plane, surface, grid_evidence, plan, observations)
    if perfect.n_passed != 8 or perfect.n_cells != 64 or perfect.whole_grid_screen.ready_for_exploratory_interpolation:
        raise ValueError("eight finite checks incorrectly promoted to whole-grid readiness")
    write_new(workdir / "perfect_8_of_64.json", perfect.to_dict())
    point = observations["p0"].point
    coordinates = point.coordinates.copy(); coordinates[2] += .04
    observations["p0"] = replace(observations["p0"], point=replace(point, energy=point.energy-.07, coordinates=coordinates))
    observations["p1"] = ConditionalProbeRecord("failed-calculation", failure="synthetic SCF failure")
    del observations["p2"]
    failed = audit_conditional_probe_coverage(plane, surface, grid_evidence, plan, observations)
    if (failed.n_passed, failed.n_failed_checks, failed.n_failed_calculations, failed.n_pending) != (5, 1, 1, 1):
        raise ValueError("failed/pending tests were silently removed")
    if not np.isclose(failed.maximum_observed_absolute_energy_error, .07) or not np.isclose(failed.maximum_observed_metric_coordinate_error, .12):
        raise ValueError("signed energy or mass-metric coordinate error mismatch")
    write_new(workdir / "failure_and_pending.json", failed.to_dict())
    rejected = []
    for label, altered_surface, altered_observations in (
        ("refit_after_registration", replace(surface, energies=surface.energies+.001), observations),
        ("training_evaluation_reused", surface, dict(observations, p3=replace(observations["p3"], evaluation_id="training-0"))),
        ("unregistered_observation", surface, dict(observations, extra=observations["p3"]))):
        try:
            audit_conditional_probe_coverage(plane, altered_surface, grid_evidence, plan, altered_observations)
        except ValueError as error:
            rejected.append(dict(case=label, reason=str(error)))
        else:
            raise ValueError(f"negative gate unexpectedly accepted {label}")
    result = dict(scope="synthetic-implementation-only", grid_nodes=81, grid_cells=64,
        registered_probes=8, passing_point_checks=8, cells_with_passing_points=8,
        whole_grid_ready=False, continuous_region_bound=False, independent_blinding_proven=False,
        calculator_calls=0, negative_gates=rejected,
        source_sha256={str(path.name): hashlib.sha256(path.read_bytes()).hexdigest() for path in (
            Path(__file__), Path(__file__).resolve().parents[1] / "vcneb/conditional_coverage.py",
            Path(__file__).resolve().parents[1] / "vcneb/conditional_evidence.py")})
    write_new(workdir / "verification.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.workdir)
    print(json.dumps({"grid_cells": result["grid_cells"], "passing_points": result["passing_point_checks"],
                      "whole_grid_ready": False, "calculator_calls": 0}))


if __name__ == "__main__":
    main()
