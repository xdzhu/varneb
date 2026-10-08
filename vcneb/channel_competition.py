"""Calculator-free, same-initial-state comparisons of declared path channels.

Discrete image maxima are not certified saddles or global minimum barriers.
Measured barrier-error bounds and explicit channel coverage are mandatory
before comparing selectivity with interval arithmetic. No error is inferred
from a NEB force threshold, and no kinetic lifetime is inferred from energy.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from typing import Mapping, Sequence

import numpy as np
from ase import Atoms

from .periodic_path import validate_periodic_path_lift


@dataclass(frozen=True)
class ChannelPath:
    name: str
    role: str  # "switching" or "decay"
    images: Sequence[Atoms]
    energies_eV_cell: Sequence[float]
    physical_contract: Mapping[str, str]
    mechanical_family: str
    mechanical_parameters: Mapping[str, object]
    pressure_eV_A3: float
    neb_fmax_eV_A: float
    source_id: str
    reverse: bool = False
    barrier_error_eV_cell: float | None = None
    sampling_audit_sha256: str | None = None
    barrier_error_audit_sha256: str | None = None


def _number(value, label, *, minimum=None):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{label} must be a finite number")
    result = float(value)
    if not np.isfinite(result) or (minimum is not None and result < minimum):
        raise ValueError(f"{label} must be finite and >= {minimum}")
    return result


def _hash(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _same_initial(left, right, tolerance):
    if (left.get_chemical_symbols() != right.get_chemical_symbols()
            or not np.array_equal(left.pbc, right.pbc)
            or not np.allclose(left.cell.array, right.cell.array, rtol=0, atol=tolerance)):
        return False
    delta = left.get_scaled_positions(wrap=False) - right.get_scaled_positions(wrap=False)
    delta[:, left.pbc] -= np.rint(delta[:, left.pbc])
    return bool(np.max(np.linalg.norm(delta @ left.cell.array, axis=1)) <= tolerance)


def _canonical_parameters(value):
    if not isinstance(value, Mapping) or not value:
        raise ValueError("explicit nonempty mechanical parameters required")
    try:
        return json.loads(json.dumps(dict(value), sort_keys=True, allow_nan=False))
    except (ValueError, TypeError) as exc:
        raise ValueError("mechanical parameters must be finite JSON data") from exc


def _difference(left, right):
    return [left[0] - right[1], left[1] - right[0]]


def summarize_competing_paths(
    paths: Sequence[ChannelPath], *, formula_units: int,
    required_channels: Mapping[str, str], fmax_target_eV_A: float = .10,
    geometry_tolerance_A: float = 1e-9, initial_energy_tolerance_eV: float = 1e-8,
) -> dict:
    """Report every path with one physical initial state and one ensemble.

    Input energies are raw E per simulation cell; E+PV is evaluated using
    each image's own volume. ``reverse`` only reverses the view, not atomic
    mappings, periodic lifts, source optimization metrics or raw data.
    Sampling checks/error bounds are caller-supplied audited evidence, not
    generated here. Required channel names prevent a missing leakage path
    from being silently omitted when a minimum is reported.
    """
    if (isinstance(formula_units, bool) or not isinstance(formula_units, (int, np.integer))
            or formula_units < 1):
        raise ValueError("formula_units must be a positive integer")
    formula_units = int(formula_units)  # keep derived flags/intervals JSON-native
    if (not paths or not isinstance(required_channels, Mapping) or not required_channels
            or set(required_channels.values()) != {"switching", "decay"}
            or any(not isinstance(n, str) or not n.strip() for n in required_channels)):
        raise ValueError("declare nonempty switching and decay channel coverage")
    target = _number(fmax_target_eV_A, "fmax target", minimum=0)
    geometry_tolerance = _number(geometry_tolerance_A, "geometry tolerance", minimum=0)
    energy_tolerance = _number(initial_energy_tolerance_eV, "initial energy tolerance", minimum=0)
    if target == 0 or geometry_tolerance == 0:
        raise ValueError("force and geometry tolerances must be positive")
    records, names, first = [], set(), None
    for path in paths:
        if (not isinstance(path, ChannelPath) or path.name in names
                or path.name not in required_channels or path.role != required_channels[path.name]
                or not isinstance(path.reverse, bool) or not isinstance(path.source_id, str)
                or not path.source_id.strip() or not isinstance(path.mechanical_family, str)
                or not path.mechanical_family.strip()):
            raise ValueError("unique declared channel, role, direction and provenance required")
        names.add(path.name)
        contract = dict(path.physical_contract)
        if not contract or any(not isinstance(k, str) or not k or not _hash(v) for k, v in contract.items()):
            raise ValueError("physical contract needs explicit input SHA256 records")
        parameters = _canonical_parameters(path.mechanical_parameters)
        pressure = _number(path.pressure_eV_A3, "pressure")
        fmax = _number(path.neb_fmax_eV_A, "source NEB residual", minimum=0)
        images = list(path.images)
        energy = np.array(path.energies_eV_cell, dtype=float, copy=True)
        if len(images) < 3 or energy.shape != (len(images),) or not np.isfinite(energy).all():
            raise ValueError("a complete finite image/energy chain required")
        for image in images:
            if (not isinstance(image, Atoms) or not len(image) or not image.pbc.all()
                    or not np.isfinite(image.positions).all() or not np.isfinite(image.cell.array).all()
                    or not np.isfinite(np.linalg.det(image.cell.array))
                    or np.linalg.det(image.cell.array) <= 0
                    or image.get_chemical_symbols() != images[0].get_chemical_symbols()):
                raise ValueError("finite ordered periodic positive-volume images required")
        validate_periodic_path_lift(images)
        enthalpy = energy + pressure * np.array([a.get_volume() for a in images])
        if not np.isfinite(enthalpy).all():
            raise ValueError("finite E+PV profile required")
        viewed = images[::-1] if path.reverse else images
        viewed_enthalpy = enthalpy[::-1] if path.reverse else enthalpy
        initial_energy = float(energy[-1] if path.reverse else energy[0])
        identity = (viewed[0], initial_energy, contract, path.mechanical_family, parameters, pressure)
        if first is None:
            first = identity
        elif (contract != first[2] or path.mechanical_family != first[3]
              or parameters != first[4] or pressure != first[5]):
            raise ValueError("mixed physical inputs or mechanical conditions cannot be compared")
        elif not _same_initial(viewed[0], first[0], geometry_tolerance):
            raise ValueError("channels do not share the same ordered periodic initial structure")
        elif abs(initial_energy - first[1]) > energy_tolerance:
            raise ValueError("common initial raw energy differs; do not shift energy zeros to hide it")
        peak = float(enthalpy.max())
        barrier = float(peak - viewed_enthalpy[0])
        interval = None
        if path.barrier_error_eV_cell is not None:
            error = _number(path.barrier_error_eV_cell, "measured barrier error", minimum=0)
            interval = [max(0., barrier - error) / formula_units, (barrier + error) / formula_units]
        sampling = path.sampling_audit_sha256 is not None
        if sampling and not _hash(path.sampling_audit_sha256):
            raise ValueError("sampling audit needs its SHA256, not an implicit pass flag")
        if path.barrier_error_audit_sha256 is not None and not _hash(path.barrier_error_audit_sha256):
            raise ValueError("barrier error audit needs its SHA256")
        records.append({
            "name": path.name, "role": path.role, "source_id": path.source_id,
            "source_direction_used": "reverse" if path.reverse else "forward",
            "n_total_images": len(images), "source_NEB_fmax_eV_A": fmax,
            "NEB_residual_passed": fmax <= target,
            "source_forward_barrier_eV_fu": float(peak - enthalpy[0]) / formula_units,
            "source_reverse_barrier_eV_fu": float(peak - enthalpy[-1]) / formula_units,
            "barrier_from_common_initial_eV_fu": barrier / formula_units,
            "reaction_energy_from_common_initial_eV_fu": float(viewed_enthalpy[-1] - viewed_enthalpy[0]) / formula_units,
            "highest_image_in_view": int(np.argmax(viewed_enthalpy)),
            "relative_profile_from_common_initial_eV_fu": ((viewed_enthalpy - viewed_enthalpy[0]) / formula_units).tolist(),
            "barrier_interval_eV_fu": interval,
            "sampling_audit_sha256": path.sampling_audit_sha256,
            "barrier_error_audit_sha256": path.barrier_error_audit_sha256,
        })
    missing = sorted(set(required_channels) - names)
    minima, intervals = {}, {}
    for role in ("switching", "decay"):
        group = [r for r in records if r["role"] == role]
        complete = all(name in names for name, required_role in required_channels.items() if required_role == role)
        minima[role] = min(r["barrier_from_common_initial_eV_fu"] for r in group) if complete else None
        intervals[role] = ([min(r["barrier_interval_eV_fu"][i] for r in group) for i in (0, 1)]
                           if complete and all(r["barrier_interval_eV_fu"] is not None for r in group) else None)
    blockers = (["missing required channels: " + ", ".join(missing)] if missing else [])
    blockers += [f"{r['name']}: {reason}" for r in records for reason, failed in (
        ("ordinary residual above target", not r["NEB_residual_passed"]),
        ("sampling audit missing", r["sampling_audit_sha256"] is None),
        ("measured barrier error missing", r["barrier_interval_eV_fu"] is None),
        ("barrier error provenance missing", r["barrier_error_audit_sha256"] is None)) if failed]
    return {
        "format_version": 1, "formula_units": int(formula_units), "energy_unit": "eV/fu",
        "required_channels": dict(required_channels), "missing_required_channels": missing,
        "physical_contract": first[2], "mechanical_family": first[3],
        "mechanical_parameters": first[4], "pressure_eV_A3": first[5],
        "common_initial_raw_energy_eV_cell": first[1], "fmax_target_eV_A": target,
        "channels": records, "provisional_minimum_barriers_eV_fu": minima,
        "minimum_barrier_intervals_eV_fu": intervals,
        "provisional_selectivity_eV_fu": minima["decay"] - minima["switching"] if not missing else None,
        "selectivity_interval_eV_fu": _difference(intervals["decay"], intervals["switching"])
        if all(v is not None for v in intervals.values()) else None,
        "ready_for_bounded_discrete_comparison": not blockers, "blockers": blockers,
        "limitations": ["comparison covers only declared candidate paths, not all possible channels",
                        "discrete image maxima are not stationary-saddle or sampling certificates",
                        "error bounds are supplied evidence, not inferred from the force tolerance",
                        "reversal retains the original optimization metric and source residual",
                        "no DFT, rates, coercive field, retention lifetime or polarization computed"],
    }


def compare_channel_selectivity(base: dict, changed: dict, *, decay_noninferiority_margin_eV_fu: float) -> dict:
    """Conservative interval changes within one mechanical boundary family.

    The noninferiority margin must be specified *before* examining responses.
    It is not guessed from overlapping error bars. Independent-error or
    Gaussian assumptions are not used; interval arithmetic can be conservative
    when barriers share correlated errors. Only declared candidate paths count.
    """
    margin = _number(decay_noninferiority_margin_eV_fu, "predeclared noninferiority margin", minimum=0)
    for key in ("energy_unit", "formula_units", "physical_contract", "mechanical_family", "pressure_eV_A3", "required_channels"):
        if base[key] != changed[key]:
            raise ValueError(f"incompatible comparison: {key}")
    if set(base["mechanical_parameters"]) != set(changed["mechanical_parameters"]):
        raise ValueError("mechanical parameter definitions differ")
    blockers = [f"{label}: {reason}" for label, record in (("base", base), ("changed", changed))
                for reason in record["blockers"]]
    if blockers:
        return {"status": "unresolved", "blockers": blockers, "claims_evaluated": False}
    changes = {role: _difference(changed["minimum_barrier_intervals_eV_fu"][role],
                                base["minimum_barrier_intervals_eV_fu"][role]) for role in ("switching", "decay")}
    selectivity = _difference(changed["selectivity_interval_eV_fu"], base["selectivity_interval_eV_fu"])
    lower_switch = changes["switching"][1] < 0
    noninferior_decay = changes["decay"][0] >= -margin
    return {
        "status": "bounded_declared_path_comparison", "claims_evaluated": True,
        "barrier_change_intervals_eV_fu": changes, "selectivity_change_interval_eV_fu": selectivity,
        "relative_selectivity_improves": selectivity[0] > 0,
        "switching_barrier_lowers": lower_switch,
        "decay_barrier_lowers": changes["decay"][1] < 0,
        "decay_barrier_increases": changes["decay"][0] > 0,
        "decay_noninferiority_margin_eV_fu": margin, "decay_noninferior_within_margin": noninferior_decay,
        "easier_switching_without_decay_loss_beyond_margin": lower_switch and noninferior_decay,
        "limitations": "bounded comparison of declared discrete paths; not global kinetics or proof that an unexamined decay channel is absent",
    }


__all__ = ["ChannelPath", "summarize_competing_paths", "compare_channel_selectivity"]
