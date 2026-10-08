"""Simple prospective controls for a registered, same-boundary path study.

These conventional baselines are not a new prediction method. A training
interval propagated through interpolation is NOT a bound on an unseen path.
The caller must audit source reports and register selection before test labels.
"""

from __future__ import annotations

import copy
import json
import math
from numbers import Integral


def _finite(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"finite numeric {name} required")
    return float(value)


def _hash(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _network(report):
    required = {"formula_units", "energy_unit", "physical_contract", "mechanical_family",
                "mechanical_parameters", "pressure_eV_A3", "required_channels", "channels",
                "common_initial_raw_energy_eV_cell", "ready_for_bounded_discrete_comparison",
                "blockers", "missing_required_channels", "fmax_target_eV_A"}
    if not isinstance(report, dict) or not required <= report.keys():
        raise ValueError("complete audited channel-network report required")
    if (report["ready_for_bounded_discrete_comparison"] is not True
            or report["blockers"] or report["missing_required_channels"]):
        raise ValueError("training needs residual, sampling and measured-error audits; incomplete network rejected")
    fu = report["formula_units"]
    if isinstance(fu, bool) or not isinstance(fu, Integral) or fu < 1 or report["energy_unit"] != "eV/fu":
        raise ValueError("positive integer formula normalization and eV/fu required")
    contract = report["physical_contract"]
    if not isinstance(contract, dict) or not contract or any(not isinstance(k, str) or not _hash(v) for k, v in contract.items()):
        raise ValueError("explicit physical-input SHA256 contract required")
    if not isinstance(report["mechanical_family"], str) or not report["mechanical_family"].strip():
        raise ValueError("explicit mechanical boundary family required")
    _finite(report["common_initial_raw_energy_eV_cell"], "initial raw energy in eV/cell")
    _finite(report["pressure_eV_A3"], "pressure in eV/A3")
    target = _finite(report["fmax_target_eV_A"], "force target in eV/A")
    if target <= 0:
        raise ValueError("positive ordinary force target required")
    declared, records = report["required_channels"], report["channels"]
    if (not isinstance(declared, dict) or set(declared.values()) != {"switching", "decay"}
            or any(not isinstance(n, str) or not n for n in declared)
            or not isinstance(records, list) or len(records) != len(declared)):
        raise ValueError("complete declared switching/decay coverage required")
    result = {}
    row_fields = {"name", "role", "source_id", "NEB_residual_passed", "source_NEB_fmax_eV_A",
                  "sampling_audit_sha256", "barrier_error_audit_sha256", "barrier_interval_eV_fu",
                  "barrier_from_common_initial_eV_fu", "reaction_energy_from_common_initial_eV_fu",
                  "relative_profile_from_common_initial_eV_fu"}
    for row in records:
        if (not isinstance(row, dict) or not row_fields <= row.keys()
                or row.get("name") not in declared or row["name"] in result):
            raise ValueError("unique declared channel records required")
        name = row["name"]
        if (row.get("role") != declared[name] or row.get("NEB_residual_passed") is not True
                or not isinstance(row.get("source_id"), str) or not row["source_id"].strip()
                or not _hash(row.get("sampling_audit_sha256"))
                or not _hash(row.get("barrier_error_audit_sha256"))):
            raise ValueError("role, provenance, residual, sampling and energy-error evidence required")
        residual = _finite(row["source_NEB_fmax_eV_A"], "source ordinary residual")
        barrier = _finite(row["barrier_from_common_initial_eV_fu"], "barrier in eV/fu")
        reaction = _finite(row["reaction_energy_from_common_initial_eV_fu"], "reaction energy in eV/fu")
        interval = row.get("barrier_interval_eV_fu")
        if not isinstance(interval, list) or len(interval) != 2:
            raise ValueError("measured training-barrier interval required")
        lower, upper = [_finite(v, "training bound in eV/fu") for v in interval]
        profile = row.get("relative_profile_from_common_initial_eV_fu")
        if not isinstance(profile, list) or len(profile) < 3:
            raise ValueError("complete relative energy profile required")
        energies = [_finite(v, "profile energy in eV/fu") for v in profile]
        if (residual < 0 or residual > target or barrier < 0 or not 0 <= lower <= barrier <= upper
                or abs(energies[0]) > 1e-10 or abs(energies[-1] - reaction) > 1e-10
                or abs(max(energies) - barrier) > 1e-10):
            raise ValueError("numeric residual/barrier/profile contradicts training audit")
        result[name] = row
    return result


def _representative(records, role):
    group = [r for r in records.values() if r["role"] == role]
    minimum_upper = min(r["barrier_interval_eV_fu"][1] for r in group)
    possible = sorted(r["name"] for r in group if r["barrier_interval_eV_fu"][0] <= minimum_upper)
    return {"selected": possible[0], "possible_lowest_within_training_bounds": possible,
            "unique_lowest_resolved": len(possible) == 1,
            "selection": "lexical first among possible minima; no target labels used"}


def _point(value, *, endpoint_reaction=None):
    result = {"raw_signed_prediction_eV_fu": float(value),
              "prediction_error_bound_eV_fu": None,
              "status": "point_prediction_not_error_certified"}
    if value < 0 or (endpoint_reaction is not None and value < max(0., endpoint_reaction) - 1e-10):
        result["status"] = "abstain_unphysical_endpoint_ordering"
        result["barrier_prediction_eV_fu"] = None
    else:
        result["barrier_prediction_eV_fu"] = float(value)
    return result


def simple_prediction_controls(lower, upper, *, parameter, target, endpoint_features=None):
    """Select training-only edges; compute B0 and two optional B1 controls.

    The conditions must differ only in one declared parameter. Select the
    representative at the lower condition even if a different channel is
    lowest at the upper condition. Endpoint features, if supplied, are visible
    test inputs, not complete test-path labels or zero-cost DFT information.
    No fitting, calculator, path optimization or unseen-error estimate occurs.
    """
    a, b = _network(lower), _network(upper)
    for key in ("formula_units", "energy_unit", "physical_contract", "mechanical_family",
                "pressure_eV_A3", "required_channels", "fmax_target_eV_A"):
        if lower[key] != upper[key]:
            raise ValueError(f"matched training contract required: {key}")
    p0, p1 = lower["mechanical_parameters"], upper["mechanical_parameters"]
    if (not isinstance(parameter, str) or not isinstance(p0, dict) or not isinstance(p1, dict)
            or parameter not in p0 or set(p0) != set(p1)
            or {k: v for k, v in p0.items() if k != parameter} != {k: v for k, v in p1.items() if k != parameter}):
        raise ValueError("one parameter in the same mechanical boundary definition may change")
    json.dumps([p0, p1], allow_nan=False)
    x0, x1, x = [_finite(v, "registered condition") for v in (p0[parameter], p1[parameter], target)]
    if not x0 < x < x1:
        raise ValueError("strictly interior unseen condition, not extrapolation, required")
    weight = (x - x0) / (x1 - x0)
    expected_target = copy.deepcopy(p0)
    expected_target[parameter] = x
    fu = int(lower["formula_units"])
    features = copy.deepcopy(endpoint_features)
    if features is not None:
        if lower["pressure_eV_A3"] != 0.:
            raise ValueError("raw-E-only endpoint controls require P=0; finite pressure needs audited volumes/enthalpy")
        expected = {"physical_contract", "mechanical_family", "formula_units", "pressure_eV_A3",
                    "anchor_parameters", "target_parameters", "initial", "final_by_channel"}
        if not isinstance(features, dict) or set(features) != expected:
            raise ValueError("only registered endpoint features, not target path labels, are accepted")
        for key in ("physical_contract", "mechanical_family", "formula_units", "pressure_eV_A3"):
            if features[key] != lower[key]:
                raise ValueError(f"endpoint feature contract differs: {key}")
        _finite(features["pressure_eV_A3"], "endpoint pressure in eV/A3")
        if features["anchor_parameters"] != p0 or features["target_parameters"] != expected_target:
            raise ValueError("endpoint features use another condition")
        finals = features["final_by_channel"]
        if not isinstance(finals, dict) or not set(finals) <= set(a):
            raise ValueError("registered final endpoint channels required")
        for row in [features["initial"], *finals.values()]:
            if not isinstance(row, dict) or set(row) != {"anchor_energy_eV_cell", "target_energy_eV_cell",
                                                       "anchor_audit_sha256", "target_audit_sha256",
                                                       "target_endpoint_DFT_calls"}:
                raise ValueError("endpoint energies, audit hashes and actual feature cost required")
            for key in ("anchor_energy_eV_cell", "target_energy_eV_cell"):
                _finite(row[key], key)
            calls = row["target_endpoint_DFT_calls"]
            if (not _hash(row["anchor_audit_sha256"]) or not _hash(row["target_audit_sha256"])
                    or isinstance(calls, bool) or not isinstance(calls, Integral) or calls < 1):
                raise ValueError("audited endpoint sources and positive recorded DFT cost required")
        if abs(features["initial"]["anchor_energy_eV_cell"] - lower["common_initial_raw_energy_eV_cell"]) > 1e-8:
            raise ValueError("anchor initial feature is not the common training well")
        for name, row in finals.items():
            expected_e = lower["common_initial_raw_energy_eV_cell"] + fu * a[name]["reaction_energy_from_common_initial_eV_fu"]
            if abs(row["anchor_energy_eV_cell"] - expected_e) > 1e-8:
                raise ValueError("anchor final feature does not match the training channel")
    outputs = []
    for role in ("switching", "decay"):
        selection, upper_selection = _representative(a, role), _representative(b, role)
        name = selection["selected"]
        r0, r1 = a[name], b[name]
        anchor = r0["barrier_from_common_initial_eV_fu"]
        endpoint_reaction = None
        if features is not None and name in features["final_by_channel"]:
            endpoint_reaction = (features["final_by_channel"][name]["target_energy_eV_cell"]
                                 - features["initial"]["target_energy_eV_cell"]) / fu
        b0 = _point((1-weight)*anchor + weight*r1["barrier_from_common_initial_eV_fu"], endpoint_reaction=endpoint_reaction)
        b0["propagated_training_reference_interval_eV_fu"] = [
            (1-weight)*r0["barrier_interval_eV_fu"][i] + weight*r1["barrier_interval_eV_fu"][i] for i in (0, 1)]
        missing = {"status": "unavailable_endpoint_features_missing", "barrier_prediction_eV_fu": None}
        b1is, b1fs = dict(missing), dict(missing)
        if features is not None:
            initial = features["initial"]
            delta_initial = (initial["target_energy_eV_cell"] - initial["anchor_energy_eV_cell"]) / fu
            b1is = _point(anchor - delta_initial, endpoint_reaction=endpoint_reaction)
            if name in features["final_by_channel"]:
                final = features["final_by_channel"][name]
                delta_final = (final["target_energy_eV_cell"] - final["anchor_energy_eV_cell"]) / fu
                b1fs = _point(anchor + delta_final - delta_initial, endpoint_reaction=endpoint_reaction)
        outputs.append({"role": role, "selected_channel": name, "lower_training_selection": selection,
                        "upper_training_selection": upper_selection,
                        "same_representative_at_both_conditions": name == upper_selection["selected"],
                        "B0_direct_barrier_interpolation": b0,
                        "B1_fixed_absolute_bottleneck": b1is,
                        "B1_bottleneck_follows_final": b1fs})
    calls = (sum(row["target_endpoint_DFT_calls"] for row in [features["initial"], *features["final_by_channel"].values()])
             if features is not None else 0)
    return {"format_version": 1, "energy_unit": "eV/fu", "formula_units": fu,
            "physical_contract": copy.deepcopy(lower["physical_contract"]),
            "pressure_eV_A3": float(lower["pressure_eV_A3"]),
            "mechanical_family": lower["mechanical_family"], "target_parameters": expected_target,
            "lower_parameters": copy.deepcopy(p0), "upper_parameters": copy.deepcopy(p1),
            "parameter": parameter, "interpolation_weight": weight, "predictions": outputs,
            "endpoint_features": features, "visible_target_endpoint_DFT_calls": int(calls),
            "new_DFT_calls_by_this_analysis": 0, "target_path_labels_read_by_this_API": False,
            "limitations": ["conventional controls, not a new predictor or novelty claim",
                            "propagated training intervals do not bound unseen-condition model error",
                            "endpoint-only features have visible DFT cost and must be shared by all controls",
                            "raw-E endpoint feature schema supports P=0 only; B0 uses training thermodynamic barriers",
                            "two selected edges do not predict the full network minimum",
                            "source audit hashes and caller label-access attestations are not independently certified here",
                            "no target accuracy, selectivity window, TS, rate or device-lifetime verdict"]}
