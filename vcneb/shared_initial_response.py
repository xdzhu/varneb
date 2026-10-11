"""Pair a local bottleneck with shared, externally audited initial energies.

No initial Hessian or calculator is required. This is a feature-conditioned
energy-gap estimate, not a blind/zero-DFT forecast or a saddle certificate.
The current raw-energy feature contract deliberately supports P=0 only.
"""

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .stationary_branch import StationaryQuadraticPoint
from .stationary_gap import ControlledStationaryModel, StationaryResponseContract, _real, _sha


@dataclass(frozen=True)
class SharedInitialEnergyResponse:
    """One visible initial-well feature, shared unchanged across all controls.

    Record the small energy increment independently; the raw totals are only
    an identity check with an explicit tolerance, not the response calculation.
    Audit geometry, phase, force/open traction and the energy reference outside
    this class. Hashes cannot certify those checks or unread path labels.
    Preparation calls include the measured endpoint's cost, even when reused.
    """

    contract: StationaryResponseContract
    parameter: float
    anchor_energy_eV_cell: float
    target_energy_eV_cell: float
    recorded_energy_change_eV_cell: float
    energy_identity_tolerance_eV_cell: float
    anchor_audit_sha256: str
    target_audit_sha256: str
    target_endpoint_DFT_calls: int

    def __post_init__(self):
        if not isinstance(self.contract, StationaryResponseContract):
            raise TypeError("explicit stationary response contract required")
        if self.contract.potential != 'E' or self.contract.pressure_eV_A3 != 0.:
            raise ValueError("raw initial-energy features support P=0 only; enthalpy/volume audit required")
        for name in ('parameter', 'anchor_energy_eV_cell', 'target_energy_eV_cell',
                     'recorded_energy_change_eV_cell', 'energy_identity_tolerance_eV_cell'):
            object.__setattr__(self, name, _real(getattr(self, name), name))
        if self.energy_identity_tolerance_eV_cell < 0.:
            raise ValueError("energy identity tolerance must be nonnegative")
        if not _sha(self.anchor_audit_sha256) or not _sha(self.target_audit_sha256):
            raise ValueError("anchor and target endpoint audit SHA256 declarations required")
        calls = self.target_endpoint_DFT_calls
        if isinstance(calls, (bool, np.bool_)) or not isinstance(calls, Integral) or calls < 1:
            raise ValueError("positive integer target endpoint preparation calls required")
        object.__setattr__(self, 'target_endpoint_DFT_calls', int(calls))
        difference = self.target_energy_eV_cell - self.anchor_energy_eV_cell
        if (not np.isfinite(difference)
                or abs(difference-self.recorded_energy_change_eV_cell) > self.energy_identity_tolerance_eV_cell):
            raise ValueError("recorded initial response contradicts the raw endpoint energies")
        if self.parameter == self.contract.anchor_parameter and self.recorded_energy_change_eV_cell != 0.:
            raise ValueError("same-condition initial response must be zero")


@dataclass(frozen=True)
class SharedInitialGapResponse:
    """Signed feature-conditioned gap; endpoint derivatives remain unknown.

    Two measured energies do not determine the initial branch derivative or
    curvature. They are deliberately not supplied as zero or a secant slope.
    The feature's preparation cost is shared, not incurred once per model.
    """

    parameter: float
    reference_gap_eV_cell: float
    stationary_anchor_correction_eV_cell: float
    stationary_anchor_gap_eV_fu: float
    initial_response_eV_cell: float
    bottleneck_response_eV_cell: float
    gap_response_eV_fu: float
    signed_stationary_gap_eV_fu: float
    shared_initial_feature: SharedInitialEnergyResponse
    bottleneck_anchor_point: StationaryQuadraticPoint
    bottleneck_point: StationaryQuadraticPoint


def stationary_gap_with_shared_initial_energy(
    initial: SharedInitialEnergyResponse, bottleneck: ControlledStationaryModel,
    *, reference_gap_eV_cell: float, reference_gap_identity_tolerance_eV_cell: float,
) -> SharedInitialGapResponse:
    """Use one visible initial response with any preregistered bottleneck arm.

    G is the raw bottleneck-centre energy minus the *measured* initial anchor
    energy. Retain the bottleneck's stationary-anchor correction; do not also
    apply an initial-model correction to an already measured initial energy.
    Evaluate the target in the bottleneck's own control scale/domain and keep
    its omitted-gradient and clamped-reaction diagnostics. Signed gaps are not
    clipped into barriers, calibrated to NEB or certified as connected TSs.

    This calculation invokes no DFT, but the shared feature has preparation
    cost. All B0--B5 must receive the same allowed feature and that cost; full
    dataset accounting, prospective freezing and independent labels are still
    external responsibilities. No feature availability/choice is frozen here.
    """
    if not isinstance(initial, SharedInitialEnergyResponse) or not isinstance(bottleneck, ControlledStationaryModel):
        raise TypeError("shared initial-energy feature and controlled bottleneck required")
    if initial.contract != bottleneck.contract:
        raise ValueError("same physical, mechanical, control and energy-reference contract required")
    if bottleneck.branch.expected_index != 1:
        raise ValueError("bottleneck model must be index1")
    gap = _real(reference_gap_eV_cell, 'recorded reference gap')
    tolerance = _real(reference_gap_identity_tolerance_eV_cell, 'reference gap identity tolerance')
    if tolerance < 0.:
        raise ValueError("reference gap identity tolerance must be nonnegative")
    raw_gap = bottleneck.branch.conditional.reference_energy - initial.anchor_energy_eV_cell
    if not np.isfinite(raw_gap) or abs(raw_gap-gap) > tolerance:
        raise ValueError("recorded gap contradicts the bottleneck and initial anchor energies")
    anchor = bottleneck.evaluate(initial.contract.anchor_parameter)
    target = bottleneck.evaluate(initial.parameter)
    t = target.control_value
    ds = anchor.control_gradient*t + .5*bottleneck.branch.control_curvature*t*t
    di = initial.recorded_energy_change_eV_cell
    correction = anchor.energy_change
    fu = initial.contract.formula_units
    values = (gap, correction, (gap+correction)/fu, di, ds,
              (ds-di)/fu, (gap+correction+ds-di)/fu)
    if not all(np.isfinite(value) for value in values):
        raise ValueError("nonfinite shared-initial response; no repair or clipping")
    return SharedInitialGapResponse(initial.parameter, *map(float, values), initial, anchor, target)
