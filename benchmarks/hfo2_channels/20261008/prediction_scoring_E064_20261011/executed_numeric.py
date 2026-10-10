"""Deterministic, shared-reference scoring of an already frozen forecast panel.

Conventional error bookkeeping, not a new predictive theory. This module
never fits/selects models, certifies a material reference, proves blinding,
or authorizes a held-out calculation. Values and bounds are eV/formula unit.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Real
from typing import Sequence

import numpy as np


def _real(value, name, *, nonnegative=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not np.isfinite(value):
        raise ValueError(f"finite real {name} required")
    result = float(value)
    if nonnegative and result < 0:
        raise ValueError(f"nonnegative {name} required")
    return result


def _name(value, name):
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"nonempty, trimmed {name} required")
    return value


def _sha(value, name):
    if not isinstance(value, str) or len(value) != 64 or any(c not in '0123456789abcdef' for c in value):
        raise ValueError(f"lowercase SHA256 {name} required")
    return value


def _interval(value):
    if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != 2:
        raise ValueError("two finite ordered reference bounds required")
    lower, upper = (_real(v, "reference bound") for v in value)
    if lower > upper:
        raise ValueError("ordered reference bounds required")
    return lower, upper


def absolute_error_interval(prediction_eV_fu, reference_interval_eV_fu):
    """Exact range of |prediction-y| for one deterministic reference interval."""
    prediction = _real(prediction_eV_fu, "prediction")
    lower, upper = _interval(reference_interval_eV_fu)
    result = (max(lower-prediction, prediction-upper, 0.),
              max(abs(prediction-lower), abs(prediction-upper)))
    if not np.isfinite(result).all():
        raise ValueError("absolute-error arithmetic overflow")
    return result


def shared_reference_gain_interval(candidate_eV_fu, baseline_predictions_eV_fu,
                                   reference_interval_eV_fu):
    """Range of min_b |b-y| - |candidate-y| on the SAME interval for y.

    Positive means smaller candidate error. For multiple fixed baselines the
    nearest-baseline envelope is a conservative post-test oracle, not a
    newly fitted deployable selector. Its breakpoints are each prediction,
    the candidate and adjacent-baseline midpoints. Endpoint/breakpoint
    evaluation gives the extrema (up to floating-point arithmetic), without
    pretending two comparisons see independent draws of the DFT reference.
    """
    candidate = _real(candidate_eV_fu, "candidate prediction")
    if not isinstance(baseline_predictions_eV_fu, (list, tuple, np.ndarray)) or not len(baseline_predictions_eV_fu):
        raise ValueError("at least one fixed baseline prediction required")
    baselines = sorted(set(_real(v, "baseline prediction") for v in baseline_predictions_eV_fu))
    lower, upper = _interval(reference_interval_eV_fu)
    nodes = [lower, upper]
    switches = [a/2+b/2 for a, b in zip(baselines[:-1], baselines[1:])]
    nodes += [v for v in [candidate, *baselines, *switches] if lower <= v <= upper]
    gain = [min(abs(b-y) for b in baselines)-abs(candidate-y) for y in nodes]
    if not np.isfinite(gain).all():
        raise ValueError("shared-reference arithmetic overflow")
    return min(gain), max(gain)


@dataclass(frozen=True)
class FrozenChannelForecast:
    """One frozen model/case record, including explicit abstention.

    A signed raw prediction can be preserved with an abstention reason; it
    is never clipped into a valid barrier or counted as a numerical success.
    Hashes identify caller-audited artifacts, not independently prove them.
    """
    model_id: str
    case_id: str
    prediction_eV_fu: float | None
    abstention_reason: str | None
    comparison_contract_sha256: str
    frozen_model_sha256: str

    def __post_init__(self):
        for name in ("model_id", "case_id"):
            _name(getattr(self, name), name)
        for name in ("comparison_contract_sha256", "frozen_model_sha256"):
            _sha(getattr(self, name), name)
        if self.abstention_reason is not None:
            _name(self.abstention_reason, "abstention reason")
        if self.prediction_eV_fu is None:
            if self.abstention_reason is None:
                raise ValueError("missing prediction needs an explicit abstention reason")
        else:
            value = _real(self.prediction_eV_fu, "barrier prediction")
            if value < 0 and self.abstention_reason is None:
                raise ValueError("negative signed prediction must be an explicit abstention")
            object.__setattr__(self, "prediction_eV_fu", value)


@dataclass(frozen=True)
class AuditedBarrierReference:
    """A caller-audited directional barrier and a measured deterministic bound.

    Its numerical bound is not a Gaussian standard deviation. Convergence,
    sampling, phase, physical input and source audits must precede this object.
    Ordinary NEB fmax cannot be supplied as a barrier-error bound by this API.
    """
    case_id: str
    barrier_eV_fu: float
    error_bound_eV_fu: float
    comparison_contract_sha256: str
    source_audit_sha256: str

    def __post_init__(self):
        _name(self.case_id, "case_id")
        for name in ("comparison_contract_sha256", "source_audit_sha256"):
            _sha(getattr(self, name), name)
        for name in ("barrier_eV_fu", "error_bound_eV_fu"):
            object.__setattr__(self, name, _real(getattr(self, name), name, nonnegative=True))
        if not np.isfinite(self.interval_eV_fu).all():
            raise ValueError("reference-bound arithmetic overflow")

    @property
    def interval_eV_fu(self):
        return (self.barrier_eV_fu-self.error_bound_eV_fu,
                self.barrier_eV_fu+self.error_bound_eV_fu)


def evaluate_frozen_panel(forecasts: Sequence[FrozenChannelForecast],
                          references: Sequence[AuditedBarrierReference], *,
                          candidate_model: str, baseline_models: Sequence[str],
                          floating_tolerance_eV_fu: float = 1e-12):
    """Score every registered case/model; abstention is coverage, not zero error.

    The full rectangular panel is required, including unavailable records.
    Missing a baseline is not evidence that it performs worse. Strict full-
    panel advantage requires every model/case available and an error gain
    above the floating arithmetic tolerance throughout each reference bound.
    This tolerance is NOT a DFT uncertainty or changed physical threshold.
    """
    _name(candidate_model, "candidate model")
    if (isinstance(baseline_models, str) or not isinstance(baseline_models, Sequence)
            or not baseline_models or len(set(baseline_models)) != len(baseline_models)):
        raise ValueError("explicit unique registered baseline models required")
    baselines = tuple(_name(m, "baseline model") for m in baseline_models)
    if candidate_model in baselines:
        raise ValueError("candidate cannot be its own baseline")
    tolerance = _real(floating_tolerance_eV_fu, "floating tolerance", nonnegative=True)
    if not references or any(not isinstance(r, AuditedBarrierReference) for r in references):
        raise ValueError("nonempty audited reference panel required")
    by_case = {r.case_id:r for r in references}
    if len(by_case) != len(references):
        raise ValueError("duplicate reference case")
    models = (candidate_model, *baselines)
    if any(not isinstance(f, FrozenChannelForecast) for f in forecasts):
        raise TypeError("frozen forecast records required")
    by_key = {(f.model_id, f.case_id):f for f in forecasts}
    if (len(by_key) != len(forecasts)
            or set(by_key) != {(model, case) for model in models for case in by_case}):
        raise ValueError("complete rectangular registered panel required; no omitted failures or extra cases")
    contract = next(iter(by_case.values())).comparison_contract_sha256
    if (any(r.comparison_contract_sha256 != contract for r in references)
            or any(f.comparison_contract_sha256 != contract for f in forecasts)):
        raise ValueError("same comparison contract required")
    if any(len({by_key[m,c].frozen_model_sha256 for c in by_case}) != 1 for m in models):
        raise ValueError("one unchanged frozen model artifact per model required")
    summaries = []
    for model in models:
        rows = []
        for case, reference in by_case.items():
            forecast = by_key[model,case]
            row = dict(case_id=case, prediction_eV_fu=forecast.prediction_eV_fu,
                       abstention_reason=forecast.abstention_reason,
                       available=forecast.abstention_reason is None,
                       raw_absolute_error_eV_fu=None, absolute_error_interval_eV_fu=None,
                       reference_source_audit_sha256=reference.source_audit_sha256)
            if row["available"]:
                row["raw_absolute_error_eV_fu"] = abs(forecast.prediction_eV_fu-reference.barrier_eV_fu)
                row["absolute_error_interval_eV_fu"] = absolute_error_interval(
                    forecast.prediction_eV_fu, reference.interval_eV_fu)
            rows.append(row)
        valid = [r for r in rows if r["available"]]
        summaries.append(dict(model_id=model, frozen_model_sha256=by_key[model,next(iter(by_case))].frozen_model_sha256,
            n_registered_cases=len(rows), n_available=len(valid), coverage=len(valid)/len(rows),
            raw_maximum_absolute_error_eV_fu=max((r["raw_absolute_error_eV_fu"] for r in valid), default=None),
            maximum_absolute_error_interval_eV_fu=(
                [max(r["absolute_error_interval_eV_fu"][i] for r in valid) for i in (0,1)] if valid else None),
            cases=rows))
    comparisons = []
    for case, reference in by_case.items():
        candidate = by_key[candidate_model,case]
        available = [m for m in baselines if by_key[m,case].abstention_reason is None]
        all_available = candidate.abstention_reason is None and len(available) == len(baselines)
        pairwise = {m:(shared_reference_gain_interval(candidate.prediction_eV_fu,
            [by_key[m,case].prediction_eV_fu], reference.interval_eV_fu)
            if candidate.abstention_reason is None else None) for m in available}
        gain = (shared_reference_gain_interval(candidate.prediction_eV_fu,
                    [by_key[m,case].prediction_eV_fu for m in available], reference.interval_eV_fu)
                if available and candidate.abstention_reason is None else None)
        comparisons.append(dict(case_id=case, available_baseline_models=available,
            pairwise_error_gain_interval_eV_fu=pairwise,
            best_available_baseline_oracle_gain_interval_eV_fu=gain,
            all_registered_models_available=all_available,
            strict_gain_over_all_available_baselines=bool(gain is not None and gain[0]>tolerance),
            complete_registered_panel_strict_gain=bool(all_available and gain[0]>tolerance)))
    return dict(format_version=1, energy_unit="eV/formula_unit", comparison_contract_sha256=contract,
        candidate_model=candidate_model, baseline_models=list(baselines), floating_tolerance_eV_fu=tolerance,
        models=summaries, comparisons=comparisons,
        full_panel_strict_advantage=all(r["complete_registered_panel_strict_gain"] for r in comparisons),
        new_DFT_calls=0, limitations=[
            "caller must independently audit physical/sampling errors and frozen-model provenance",
            "hash declarations do not prove chronological blinding or an unseen-label audit",
            "strongest-available baseline envelope is a conservative scoring oracle, not fitted selection",
            "per-model error maxima with different coverage are not a fair performance comparison",
            "deterministic reference intervals are not statistical confidence or model-error guarantees",
            "scores alone do not establish equal-information/cost or cross-material generalization",
            "this scoring does not freeze B2-B5, select material branches, or establish JCTC readiness"])
