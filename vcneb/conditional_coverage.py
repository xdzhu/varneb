"""Finite, preregistered conditional-surface checks; no whole-area certificate.

This supplements, and never relaxes, ``screen_conditional_interpolation``.
Registration uses training data only. Publish the manifest before measuring
probes: a hash or caller attestation alone cannot prove prospective blinding.
No calculator, interpolation fit or scheduler is invoked by this module.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
from typing import Mapping, Sequence

import numpy as np

from .conditional_evidence import (
    ConditionalHoldout, ConditionalInterpolationScreen, ConditionalPointEvidence,
    screen_conditional_interpolation,
)
from .mode_surface import ConditionalModeSurface, ModePlane


def _digest(data: dict) -> str:
    raw = json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _snapshot(plane, surface, evidence) -> str:
    def arrays(record):
        return {key: value.tolist() if isinstance(value, np.ndarray) else value
                for key, value in asdict(record).items()}

    return _digest({"plane": arrays(plane), "surface": arrays(surface),
                    "point_evidence": [[asdict(item) for item in row] for row in evidence],
                    "predictor": "four-corner-bilinear-v1"})


def _validate_locations(plan, surface):
    if len(plan.training_evaluation_ids) != surface.energies.size:
        raise ValueError("one training evaluation ID per measured grid node is required")
    for probe in plan.probes:
        i, j = probe.cell_i, probe.cell_j
        if not (i < len(surface.q1) - 1 and j < len(surface.q2) - 1):
            raise ValueError("registered probe cell lies outside the measured grid")
        if not (surface.q1[i] < probe.q[0] < surface.q1[i + 1]
                and surface.q2[j] < probe.q[1] < surface.q2[j + 1]):
            raise ValueError("registered probe must lie strictly inside its cell")


@dataclass(frozen=True)
class ConditionalProbeSpec:
    """One fixed, strictly interior location, selected without its label."""

    probe_id: str
    cell_i: int
    cell_j: int
    q: tuple[float, float]

    def __post_init__(self):
        if not isinstance(self.probe_id, str) or not self.probe_id.strip():
            raise ValueError("a nonempty probe_id is required")
        for value in (self.cell_i, self.cell_j):
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 0:
                raise ValueError("cell indices must be nonnegative integers")
        q = np.asarray(self.q, dtype=float)
        if q.shape != (2,) or not np.all(np.isfinite(q)):
            raise ValueError("a probe requires two finite Q coordinates")
        object.__setattr__(self, "q", tuple(float(value) for value in q))
        object.__setattr__(self, "cell_i", int(self.cell_i))
        object.__setattr__(self, "cell_j", int(self.cell_j))


@dataclass(frozen=True)
class ConditionalProbePlan:
    """Training snapshot, test locations and thresholds, with a JSON manifest.

    The supplied evaluation IDs and policy IDs are declarations; source audits
    must establish actual independence, physical inputs and branch starts.
    """

    training_snapshot_sha256: str
    training_evaluation_ids: tuple[str, ...]
    probes: tuple[ConditionalProbeSpec, ...]
    protocol_id: str
    start_policy_id: str
    selection_rule_id: str
    energy_zero_id: str
    gradient_tolerance: float
    maximum_neighbor_orthogonal_jump: float
    maximum_holdout_energy_error: float
    maximum_holdout_coordinate_error: float

    def __post_init__(self):
        digest = self.training_snapshot_sha256
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("training snapshot requires a SHA256 digest")
        ids = tuple(self.training_evaluation_ids)
        probes = tuple(self.probes)
        if (not ids or any(not isinstance(value, str) or not value.strip() for value in ids)
                or len(set(ids)) != len(ids)):
            raise ValueError("training evaluation IDs must be nonempty and unique")
        if not probes or any(not isinstance(item, ConditionalProbeSpec) for item in probes):
            raise ValueError("a nonempty sequence of ConditionalProbeSpec is required")
        if len({item.probe_id for item in probes}) != len(probes) or len({item.q for item in probes}) != len(probes):
            raise ValueError("probe IDs and Q positions must be unique")
        for name in ("protocol_id", "start_policy_id", "selection_rule_id", "energy_zero_id"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ValueError(f"nonempty {name} is required")
        for name in ("gradient_tolerance", "maximum_neighbor_orthogonal_jump",
                     "maximum_holdout_energy_error", "maximum_holdout_coordinate_error"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not np.isfinite(value) or value <= 0:
                raise ValueError("probe thresholds must be finite and positive")
            object.__setattr__(self, name, float(value))
        object.__setattr__(self, "training_evaluation_ids", ids)
        object.__setattr__(self, "probes", probes)

    def settings(self) -> dict:
        return {key: getattr(self, key) for key in (
            "energy_zero_id", "gradient_tolerance", "maximum_neighbor_orthogonal_jump",
            "maximum_holdout_energy_error", "maximum_holdout_coordinate_error")}

    def to_dict(self) -> dict:
        data = asdict(self)
        data["training_evaluation_ids"] = list(self.training_evaluation_ids)
        data["probes"] = [dict(asdict(item), q=list(item.q)) for item in self.probes]
        data.update(format_version=1, caller_attests_probe_labels_unread=True,
                    independent_blinding_proven=False, scope="finite-pointwise-bilinear-checks")
        return dict(data, registration_sha256=_digest(data))

    @classmethod
    def from_dict(cls, data: dict) -> "ConditionalProbePlan":
        payload = dict(data)
        expected = payload.pop("registration_sha256", None)
        if _digest(payload) != expected:
            raise ValueError("registered probe manifest checksum mismatch")
        for name, value in (("format_version", 1), ("caller_attests_probe_labels_unread", True),
                            ("independent_blinding_proven", False), ("scope", "finite-pointwise-bilinear-checks")):
            if payload.pop(name, None) != value:
                raise ValueError("unsupported probe manifest semantics")
        payload["probes"] = tuple(ConditionalProbeSpec(**item) for item in payload["probes"])
        plan = cls(**payload)
        if plan.to_dict() != data:
            raise ValueError("noncanonical registered probe manifest")
        return plan


def register_conditional_probe_plan(
    plane: ModePlane, surface: ConditionalModeSurface,
    point_evidence: Sequence[Sequence[ConditionalPointEvidence]],
    probes: Sequence[ConditionalProbeSpec], *, training_evaluation_ids: Sequence[str],
    protocol_id: str, start_policy_id: str, selection_rule_id: str,
    probe_labels_unread: bool = False, **settings,
) -> ConditionalProbePlan:
    """Bind a plan to training data only; persist/publish before probe labels.

    ``settings`` are the five explicit settings of the unchanged whole-grid
    screen. Energy error uses surface.energy_unit; coordinate error and jump
    use the plane's metric/amplitude unit, not raw Cartesian Angstrom by default.
    This function neither selects modes/probes nor reads test observations.
    """
    if probe_labels_unread is not True:
        raise ValueError("explicit caller attestation of unread probe labels required")
    # Validate the same coordinate/evidence contract without any test labels.
    screen_conditional_interpolation(plane, surface, point_evidence, (), **settings)
    if len(training_evaluation_ids) != surface.energies.size:
        raise ValueError("one training evaluation ID per measured grid node is required")
    plan = ConditionalProbePlan(
        training_snapshot_sha256=_snapshot(plane, surface, point_evidence),
        training_evaluation_ids=tuple(training_evaluation_ids), probes=tuple(probes),
        protocol_id=protocol_id, start_policy_id=start_policy_id,
        selection_rule_id=selection_rule_id, **settings,
    )
    _validate_locations(plan, surface)
    return plan


@dataclass(frozen=True)
class ConditionalProbeRecord:
    """One new evaluation, or an explicit failed calculation (not zero energy)."""

    evaluation_id: str
    point: ConditionalHoldout | None = None
    failure: str | None = None

    def __post_init__(self):
        if not isinstance(self.evaluation_id, str) or not self.evaluation_id.strip():
            raise ValueError("a nonempty probe evaluation_id is required")
        if self.point is not None:
            if not isinstance(self.point, ConditionalHoldout) or self.failure is not None:
                raise ValueError("a measured holdout cannot also be a failed calculation")
        elif not isinstance(self.failure, str) or not self.failure.strip():
            raise ValueError("missing measurements require an explicit calculation failure")


@dataclass(frozen=True)
class ConditionalProbeResult:
    probe_id: str
    cell_i: int
    cell_j: int
    q: tuple[float, float]
    status: str
    measured_minus_predicted_energy: float | None
    metric_coordinate_error: float | None
    blockers: tuple[str, ...]


@dataclass(frozen=True)
class ConditionalProbeCoverage:
    """Observed test-point coverage, never a continuous-cell error bound."""

    registration_sha256: str
    energy_unit: str
    coordinate_error_unit: str
    n_cells: int
    n_registered: int
    n_measured: int
    n_passed: int
    n_failed_checks: int
    n_failed_calculations: int
    n_pending: int
    cells_with_measurements: tuple[tuple[int, int], ...]
    cells_with_passed_point_checks: tuple[tuple[int, int], ...]
    maximum_observed_absolute_energy_error: float | None
    maximum_observed_metric_coordinate_error: float | None
    point_results: tuple[ConditionalProbeResult, ...]
    whole_grid_screen: ConditionalInterpolationScreen

    def to_dict(self) -> dict:
        return dict(asdict(self), scope="finite-pointwise-bilinear-checks",
                    continuous_region_error_bound_established=False,
                    prospective_blinding_independently_proven=False,
                    new_calculator_calls=0)


def audit_conditional_probe_coverage(
    plane: ModePlane, surface: ConditionalModeSurface,
    point_evidence: Sequence[Sequence[ConditionalPointEvidence]],
    plan: ConditionalProbePlan, observations: Mapping[str, ConditionalProbeRecord],
) -> ConditionalProbeCoverage:
    """Report every registered point, including failures and missing results.

    Each measured point uses the existing screen on its four-corner patch.
    The full-grid screen is also returned, unchanged. A passing point cannot
    certify its cell, unsampled space, global lower envelope, TS or barrier.
    Evaluation IDs catch declared reuse; hashes/IDs do not authenticate raw
    calculation provenance or establish that labels were prospectively hidden.
    """
    if not isinstance(plan, ConditionalProbePlan):
        raise TypeError("a registered ConditionalProbePlan is required")
    if _snapshot(plane, surface, point_evidence) != plan.training_snapshot_sha256:
        raise ValueError("training prediction changed after probe registration")
    _validate_locations(plan, surface)
    ids = {item.probe_id for item in plan.probes}
    if set(observations) - ids:
        raise ValueError("unregistered probe observations must not be silently included")
    evaluation_ids = set(plan.training_evaluation_ids)
    for record in observations.values():
        if not isinstance(record, ConditionalProbeRecord):
            raise TypeError("observations require ConditionalProbeRecord values")
        if record.evaluation_id in evaluation_ids:
            raise ValueError("probe evaluation ID repeats a training or test evaluation")
        evaluation_ids.add(record.evaluation_id)

    rows, measured, passed_cells = [], [], set()
    for probe in plan.probes:
        record = observations.get(probe.probe_id)
        common = dict(probe_id=probe.probe_id, cell_i=probe.cell_i, cell_j=probe.cell_j, q=probe.q)
        if record is None or record.point is None:
            rows.append(ConditionalProbeResult(
                **common, status="pending" if record is None else "calculation_failed",
                measured_minus_predicted_energy=None, metric_coordinate_error=None,
                blockers=("no observation",) if record is None else (record.failure,),
            ))
            continue
        point = record.point
        if (point.cell_i != probe.cell_i or point.cell_j != probe.cell_j
                or not np.array_equal(np.asarray(point.q), probe.q)):
            raise ValueError("observation does not match its registered cell/Q location")
        i, j = probe.cell_i, probe.cell_j
        patch_fields = {
            name: getattr(surface, name)[j:j + 2, i:i + 2]
            for name in ("energies", "orthogonal_gradient_norms", "n_evaluations", "selected_starts", "coordinates")
        }
        patch = replace(surface, q1=surface.q1[i:i + 2], q2=surface.q2[j:j + 2], **patch_fields)
        evidence = [row[i:i + 2] for row in point_evidence[j:j + 2]]
        check = screen_conditional_interpolation(
            plane, patch, evidence, [replace(point, cell_i=0, cell_j=0)], **plan.settings(),
        )
        tx = (probe.q[0] - patch.q1[0]) / np.diff(patch.q1)[0]
        ty = (probe.q[1] - patch.q2[0]) / np.diff(patch.q2)[0]
        weights = np.array([[(1 - tx)*(1 - ty), tx*(1 - ty)], [(1 - tx)*ty, tx*ty]])
        error = float(point.energy - np.sum(weights * patch.energies))
        status = "point_checks_passed" if check.ready_for_exploratory_interpolation else "point_checks_failed"
        rows.append(ConditionalProbeResult(
            **common, status=status, measured_minus_predicted_energy=error,
            metric_coordinate_error=check.maximum_holdout_coordinate_error, blockers=check.blockers,
        ))
        measured.append(point)
        if check.ready_for_exploratory_interpolation:
            passed_cells.add((i, j))
    whole = screen_conditional_interpolation(plane, surface, point_evidence, measured, **plan.settings())
    return ConditionalProbeCoverage(
        registration_sha256=plan.to_dict()["registration_sha256"],
        energy_unit=surface.energy_unit, coordinate_error_unit=plane.amplitude_unit,
        n_cells=whole.n_cells, n_registered=len(plan.probes), n_measured=len(measured),
        n_passed=sum(row.status == "point_checks_passed" for row in rows),
        n_failed_checks=sum(row.status == "point_checks_failed" for row in rows),
        n_failed_calculations=sum(row.status == "calculation_failed" for row in rows),
        n_pending=sum(row.status == "pending" for row in rows),
        cells_with_measurements=tuple(sorted({(p.cell_i, p.cell_j) for p in measured})),
        cells_with_passed_point_checks=tuple(sorted(passed_cells)),
        maximum_observed_absolute_energy_error=max((abs(row.measured_minus_predicted_energy)
                                                   for row in rows if row.measured_minus_predicted_energy is not None), default=None),
        maximum_observed_metric_coordinate_error=whole.maximum_holdout_coordinate_error,
        point_results=tuple(rows), whole_grid_screen=whole,
    )


__all__ = ["ConditionalProbeSpec", "ConditionalProbePlan", "ConditionalProbeRecord",
           "ConditionalProbeResult", "ConditionalProbeCoverage",
           "register_conditional_probe_plan", "audit_conditional_probe_coverage"]
