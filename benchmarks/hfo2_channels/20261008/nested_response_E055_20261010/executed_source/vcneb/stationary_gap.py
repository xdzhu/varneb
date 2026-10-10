"""Pair restricted stationary energies under one prescribed external parameter.

This is local harmonic energy-gap response, not a certified activation barrier,
branch selector or independent material prediction. Contracts and domains are
caller-audited declarations; hashes cannot certify their physical correctness.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .stationary_branch import StationaryQuadraticBranch, StationaryQuadraticPoint


def _real(value, name):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not np.isfinite(value)):
        raise ValueError(f"finite real {name} required")
    return float(value)


def _sha(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


@dataclass(frozen=True)
class StationaryResponseContract:
    """One audited Hamiltonian, boundary family, energy zero and control.

    The boundary digest must bind the common substrate reference, open freedom,
    loads/field and other fixed boundary settings, not a strain-dependent cell.
    Coordinates are in Angstrom and the potential is in eV per simulation cell.
    Distinct local centres may use distinct coordinate scales.
    """

    physical_contract: Mapping[str, str]
    mechanical_family: str
    boundary_definition_sha256: str
    parameter: str
    parameter_unit: str
    anchor_parameter: float
    energy_zero_id: str
    formula_units: int
    pressure_eV_A3: float
    potential: str

    def __post_init__(self):
        hashes = self.physical_contract
        if (not isinstance(hashes, Mapping) or not hashes
                or any(not isinstance(k, str) or not k.strip() or not _sha(v) for k, v in hashes.items())
                or not _sha(self.boundary_definition_sha256)):
            raise ValueError("physical-input and fixed-boundary SHA256 declarations required")
        for key in ('mechanical_family', 'parameter', 'parameter_unit', 'energy_zero_id'):
            if not isinstance(getattr(self, key), str) or not getattr(self, key).strip():
                raise ValueError(f"explicit {key} required")
        if (isinstance(self.formula_units, (bool, np.bool_))
                or not isinstance(self.formula_units, Integral) or self.formula_units < 1):
            raise ValueError("positive integer formula normalization required")
        pressure = _real(self.pressure_eV_A3, 'pressure')
        if self.potential not in ('E', 'H=E+PV') or (self.potential == 'E' and pressure != 0.):
            raise ValueError("E supports P=0 only; finite pressure requires same-reference H=E+PV")
        object.__setattr__(self, 'physical_contract', MappingProxyType(dict(sorted(hashes.items()))))
        object.__setattr__(self, 'anchor_parameter', _real(self.anchor_parameter, 'anchor parameter'))
        object.__setattr__(self, 'pressure_eV_A3', pressure)
        object.__setattr__(self, 'formula_units', int(self.formula_units))


@dataclass(frozen=True)
class ControlledStationaryModel:
    """Explicit local scaling and caller-declared, not automatically proved, domain.

    t = scale_A_per_parameter_unit * (x-anchor). Bounds and displacement radius
    are necessary checks, not proof of probe-convex-hull coverage, anharmonic
    accuracy or transverse stability. No unseen point is used to select them.
    """

    branch: StationaryQuadraticBranch
    contract: StationaryResponseContract
    scale_A_per_parameter_unit: float
    parameter_interval: tuple[float, float]
    maximum_full_displacement_A: float
    source_audit_sha256: str

    def __post_init__(self):
        if not isinstance(self.branch, StationaryQuadraticBranch) or not isinstance(self.contract, StationaryResponseContract):
            raise TypeError("stationary quadratic and explicit response contract required")
        scale = _real(self.scale_A_per_parameter_unit, 'control scale')
        radius = _real(self.maximum_full_displacement_A, 'declared displacement radius')
        if scale <= 0 or radius <= 0 or not _sha(self.source_audit_sha256):
            raise ValueError("positive control scale/radius and source audit SHA256 required")
        if not isinstance(self.parameter_interval, (tuple, list)) or len(self.parameter_interval) != 2:
            raise ValueError("two declared parameter bounds required")
        lo, hi = [_real(v, 'domain bound') for v in self.parameter_interval]
        if not lo < hi or not lo <= self.contract.anchor_parameter <= hi:
            raise ValueError("ordered parameter bounds must include the anchor")
        object.__setattr__(self, 'scale_A_per_parameter_unit', scale)
        object.__setattr__(self, 'maximum_full_displacement_A', radius)
        object.__setattr__(self, 'parameter_interval', (lo, hi))

    def evaluate(self, parameter):
        x = _real(parameter, 'prescribed parameter')
        if not self.parameter_interval[0] <= x <= self.parameter_interval[1]:
            raise ValueError("outside declared parameter domain; no extrapolation repair")
        point = self.branch.evaluate(self.scale_A_per_parameter_unit * (x-self.contract.anchor_parameter))
        if np.linalg.norm(point.full_displacement) > self.maximum_full_displacement_A:
            raise ValueError("stationary displacement leaves declared radius; no automatic recentering")
        return point


@dataclass(frozen=True)
class StationaryGapResponse:
    """Signed gap and its response; no full-DFT TS/barrier certificate is implied."""

    parameter: float
    reference_gap_eV_cell: float
    stationary_anchor_correction_eV_cell: float
    stationary_anchor_gap_eV_fu: float
    initial_response_eV_cell: float
    bottleneck_response_eV_cell: float
    gap_response_eV_fu: float
    signed_stationary_gap_eV_fu: float
    gap_derivative_eV_fu_per_parameter_unit: float
    gap_curvature_eV_fu_per_parameter_unit2: float
    initial_point: StationaryQuadraticPoint
    bottleneck_point: StationaryQuadraticPoint


def restricted_stationary_gap_response(
    initial: ControlledStationaryModel, bottleneck: ControlledStationaryModel,
    parameter: float, *, reference_gap_eV_cell: float,
    reference_gap_identity_tolerance_eV_cell: float,
) -> StationaryGapResponse:
    """Pair an internal minimum and index-one saddle in matched declarations.

    Supply the independently recorded centre-to-centre gap without subtracting
    large raw totals to calculate response. Cross-check those reference totals
    with the caller's explicit identity tolerance; this tolerance is not an
    unseen prediction error or a new DFT setting. Keep reference relaxation
    corrections. Compute each energy increment about its own stationary anchor
    and apply its own external-coordinate chain rule. Negative gaps stay signed.

    Omitted admissible gradients and clamped reactions remain in both points.
    Resolved *represented* indices do not certify an unsampled full DFT saddle,
    continuity, the lowest channel or the highest point of a complete path.
    """
    if not isinstance(initial, ControlledStationaryModel) or not isinstance(bottleneck, ControlledStationaryModel):
        raise TypeError("two controlled stationary models required")
    if initial.contract != bottleneck.contract:
        raise ValueError("same physical, mechanical, control and energy-reference contract required")
    if initial.branch.expected_index != 0 or bottleneck.branch.expected_index != 1:
        raise ValueError("initial model must be index0 and bottleneck model index1")
    gap = _real(reference_gap_eV_cell, 'recorded reference gap')
    tolerance = _real(reference_gap_identity_tolerance_eV_cell, 'reference gap identity tolerance')
    if tolerance < 0:
        raise ValueError("reference gap identity tolerance must be nonnegative")
    raw_gap = bottleneck.branch.conditional.reference_energy - initial.branch.conditional.reference_energy
    if abs(raw_gap-gap) > tolerance:
        raise ValueError("recorded gap contradicts the paired reference energies")
    x = _real(parameter, 'prescribed parameter')
    contract = initial.contract
    p0, s0 = initial.evaluate(contract.anchor_parameter), bottleneck.evaluate(contract.anchor_parameter)
    p, s = initial.evaluate(x), bottleneck.evaluate(x)
    ti, ts = p.control_value, s.control_value
    di = p0.control_gradient*ti + .5*initial.branch.control_curvature*ti*ti
    ds = s0.control_gradient*ts + .5*bottleneck.branch.control_curvature*ts*ts
    correction = s0.energy_change - p0.energy_change
    fu = contract.formula_units
    derivative = (bottleneck.scale_A_per_parameter_unit*s.control_gradient
                  - initial.scale_A_per_parameter_unit*p.control_gradient) / fu
    curvature = (bottleneck.scale_A_per_parameter_unit*bottleneck.scale_A_per_parameter_unit*bottleneck.branch.control_curvature
                 - initial.scale_A_per_parameter_unit*initial.scale_A_per_parameter_unit*initial.branch.control_curvature) / fu
    values = (gap, correction, (gap+correction)/fu, di, ds, (ds-di)/fu,
              (gap+correction+ds-di)/fu, derivative, curvature)
    if not all(np.isfinite(v) for v in values):
        raise ValueError("nonfinite paired response; no repair or clipping")
    return StationaryGapResponse(x, *map(float, values), p, s)
