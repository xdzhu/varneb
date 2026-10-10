"""Nested stationary-response controls from measured Hessian columns.

This is standard restricted harmonic algebra and training-only instability
promotion, not a full anharmonic branch search or a material barrier predictor.
Unmeasured Hessian blocks never authorize release or certify stability.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np

from .stationary_branch import StationaryQuadraticBranch, stationary_quadratic_branch
from .stationary_gap import (
    ControlledStationaryModel, StationaryResponseContract, StationaryGapResponse,
    restricted_stationary_gap_response,
)


def _array(value, name, *, ndim):
    if np.iscomplexobj(value):
        raise ValueError(f"{name} requires a declared real chart; complex phases cannot be discarded")
    a = np.array(value, dtype=float, copy=True)
    if a.ndim != ndim or not np.isfinite(a).all():
        raise ValueError(f"finite {ndim}-dimensional {name} required")
    a.setflags(write=False)
    return a


def _real(value, name):
    if isinstance(value, (bool,np.bool_)) or not isinstance(value, Real) or not np.isfinite(value):
        raise ValueError(f"finite real {name} required")
    return float(value)


def _basis(value, n, name, *, allow_empty=False):
    b = _array(value,name,ndim=2)
    if (b.shape[0]!=n or (not allow_empty and not b.shape[1])
            or not np.allclose(b.T@b,np.eye(b.shape[1]),atol=1e-10,rtol=0)):
        raise ValueError(f"orthonormal {name} in the declared chart required")
    return b


def _inside(a,b):
    return bool(np.allclose(b@(b.T@a),a,atol=1e-10,rtol=0))


def _complement(outer, inner):
    """Orthonormal complement inside outer; never project a forbidden basis."""
    _, singular, vt = np.linalg.svd(inner.T@outer,full_matrices=True)
    rank = int(np.count_nonzero(singular>1e-10))
    return outer@vt[rank:].T


@dataclass(frozen=True)
class MeasuredStationaryCentre:
    """Same-centre full gradient and directional Hessian data, in Å/eV.

    M is an orthonormal *measured* basis. D[:,j] is d(physical gradient)/dM_j,
    including rows outside M; it is not a NEB/projected force derivative.
    Q and atomic/internal spaces are explicit training-audited declarations.
    No full Hessian, mode orthogonalization, source audit or cost is invented.
    """

    measured_basis: np.ndarray
    hessian_columns: np.ndarray
    gradient: np.ndarray
    reference_potential_eV_cell: float
    retained_internal_basis: np.ndarray
    admissible_internal_basis: np.ndarray
    atomic_internal_basis: np.ndarray
    control_direction: np.ndarray
    expected_index: int
    stability_floor_eV_A2: float
    internal_curvature_floor_eV_A2: float
    symmetry_tolerance_eV_A2: float
    source_audit_sha256: str
    measurement_ids: tuple[str,...]

    def __post_init__(self):
        g = _array(self.gradient,"physical gradient",ndim=1)
        n = g.size
        m = _basis(self.measured_basis,n,"measured basis")
        q = _basis(self.retained_internal_basis,n,"retained internal basis")
        b = _basis(self.admissible_internal_basis,n,"admissible internal basis")
        a = _basis(self.atomic_internal_basis,n,"atomic internal basis",allow_empty=True)
        c = _array(self.control_direction,"control",ndim=1)
        d = _array(self.hessian_columns,"Hessian columns",ndim=2)
        if (c.shape!=(n,) or d.shape!=m.shape or abs(np.linalg.norm(c)-1)>1e-10
                or not np.allclose(b.T@c,0,atol=1e-10,rtol=0)
                or not _inside(q,b) or not _inside(a,b)
                or not _inside(m,np.column_stack((b,c)))):
            raise ValueError("retained/atomic/measured spaces must obey the actual mechanical/control chart")
        if (isinstance(self.expected_index,(bool,np.bool_)) or not isinstance(self.expected_index,Integral)
                or self.expected_index not in (0,1)):
            raise ValueError("expected internal stationary index must be integer0 or1")
        floor = _real(self.stability_floor_eV_A2,"release stability floor")
        internal = _real(self.internal_curvature_floor_eV_A2,"internal curvature floor")
        tolerance = _real(self.symmetry_tolerance_eV_A2,"measured symmetry tolerance")
        if floor<0 or internal<0 or tolerance<=0:
            raise ValueError("nonnegative measured curvature floors and positive symmetry tolerance required")
        projected=m.T@d
        if np.max(np.abs(projected-projected.T))>tolerance:
            raise ValueError("measured Hessian columns contradict reciprocity within the declared tolerance")
        sha=self.source_audit_sha256
        ids=self.measurement_ids
        if (not isinstance(sha,str) or len(sha)!=64 or any(s not in '0123456789abcdef' for s in sha)
                or not isinstance(ids,(list,tuple)) or not ids
                or any(not isinstance(s,str) or not s.strip() for s in ids) or len(set(ids))!=len(ids)):
            raise ValueError("source audit SHA256 and unique actual measurement IDs required")
        for name,value in (("gradient",g),("measured_basis",m),("hessian_columns",d),
                           ("retained_internal_basis",q),("admissible_internal_basis",b),
                           ("atomic_internal_basis",a),("control_direction",c)):
            object.__setattr__(self,name,value)
        object.__setattr__(self,"reference_potential_eV_cell",_real(self.reference_potential_eV_cell,"reference potential"))
        object.__setattr__(self,"stability_floor_eV_A2",floor)
        object.__setattr__(self,"internal_curvature_floor_eV_A2",internal)
        object.__setattr__(self,"symmetry_tolerance_eV_A2",tolerance)
        object.__setattr__(self,"measurement_ids",tuple(ids))

    def action_extension(self):
        """Symmetric action on M only; the M-orthogonal block is UNKNOWN.

        A zero computational completion outside M is never used as a physical
        Hessian: every retained/released/control column must lie inside M.
        Full physical-gradient rows remain available for omission diagnostics.
        """
        m,d=self.measured_basis,self.hessian_columns
        projected=m.T@d
        symmetric=.5*(projected+projected.T)
        corrected=d+m@(symmetric-projected)
        return corrected@m.T+m@corrected.T-m@symmetric@m.T


@dataclass(frozen=True)
class NestedControlResult:
    control: str
    status: str
    reason: str | None
    branch: StationaryQuadraticBranch | None
    retained_dimension: int
    released_dimension: int
    promoted_dimension: int
    unmeasured_internal_dimension: int
    measurement_ids: tuple[str,...]
    source_audit_sha256: str


def nested_stationary_controls(centre: MeasuredStationaryCentre) -> tuple[NestedControlResult,...]:
    """B2 frozen / B3 atomic / B4 joint / B5 training-instability promotion.

    All controls see the SAME measured columns and gradient. B5 promotes the
    complete R eigensubspace at/below the registered floor using training H,
    then checks the combined internal index without pseudoinversion. This
    implements B5's local increment-dimension arm, NOT discovery/selection of
    multiple anharmonic branches. A failure is coverage/assumption failure,
    never a zero barrier or evidence of poorer material prediction accuracy.
    """
    if not isinstance(centre,MeasuredStationaryCentre):
        raise TypeError("measured stationary centre required")
    x=centre
    h=x.action_extension()
    q,b,a,c=x.retained_internal_basis,x.admissible_internal_basis,x.atomic_internal_basis,x.control_direction
    joint=_complement(b,q)
    empty=np.empty((len(c),0))
    unmeasured=int(np.linalg.matrix_rank(b-x.measured_basis@(x.measured_basis.T@b),tol=1e-10))
    outputs=[]
    for name in ('B2_frozen','B3_atomic_release','B4_joint_release','B5_training_instability_promotion'):
        retained=q
        released=empty
        promoted=0
        reason=None
        status='available_restricted_quadratic_only'
        if name=='B3_atomic_release':
            if not _inside(q,a):
                status='unavailable_incompatible_atomic_ablation'
                reason='retained coordinates include cell motion; cannot call this an atomic-only control'
            else:
                released=_complement(a,q)
        elif name in ('B4_joint_release','B5_training_instability_promotion'):
            released=joint
        needed=np.column_stack((retained,released,c))
        if reason is None and not _inside(needed,x.measured_basis):
            status='unavailable_unmeasured_curvature'
            reason='required Hessian action outside the measured space; no stable or zero-block assumption'
        if reason is None and name=='B5_training_instability_promotion' and released.shape[1]:
            values,vectors=np.linalg.eigh(released.T@h@released)
            unstable=values<=x.stability_floor_eV_A2
            promoted=int(np.count_nonzero(unstable))
            retained=np.column_stack((q,released@vectors[:,unstable]))
            released=released@vectors[:,~unstable]
        branch=None
        if reason is None:
            try:
                branch=stationary_quadratic_branch(h,x.gradient,x.reference_potential_eV_cell,
                    retained,released,c,admissible_internal_basis=b,expected_index=x.expected_index,
                    stability_floor=x.stability_floor_eV_A2,
                    internal_curvature_floor=x.internal_curvature_floor_eV_A2,
                    symmetry_tolerance=x.symmetry_tolerance_eV_A2)
            except ValueError as error:
                status='unavailable_stationary_or_release_assumption'
                reason=str(error)
        outputs.append(NestedControlResult(name,status,reason,branch,retained.shape[1],released.shape[1],
                       promoted,unmeasured,x.measurement_ids,x.source_audit_sha256))
    return tuple(outputs)


@dataclass(frozen=True)
class NestedGapResult:
    control: str
    status: str
    reason: str | None
    initial: NestedControlResult
    bottleneck: NestedControlResult
    response: StationaryGapResponse | None
    available_measurement_ids: tuple[str,...]


def nested_stationary_gap_responses(
    initial: MeasuredStationaryCentre, bottleneck: MeasuredStationaryCentre,
    parameter: float, *, contract: StationaryResponseContract,
    initial_domain: dict, bottleneck_domain: dict,
    reference_gap_eV_cell: float, reference_gap_identity_tolerance_eV_cell: float,
) -> tuple[NestedGapResult,...]:
    """Pair all four controls with explicit distinct domains/control scales.

    Domain dictionaries contain only scale_A_per_parameter_unit,
    parameter_interval and maximum_full_displacement_A. The supplied common
    contract must be externally verified for BOTH measurement datasets.
    Signed responses, affine corrections and omission diagnostics are preserved;
    unavailable controls have no numeric prediction. No network selection,
    calibration, prediction-error bound, cost timing or label freezing occurs.
    """
    if not isinstance(initial,MeasuredStationaryCentre) or not isinstance(bottleneck,MeasuredStationaryCentre):
        raise TypeError("two measured stationary centres required")
    if not isinstance(contract,StationaryResponseContract):
        raise TypeError("shared physical/energy/boundary response contract required")
    if initial.expected_index!=0 or bottleneck.expected_index!=1:
        raise ValueError("initial index0 and bottleneck index1 required")
    fields={'scale_A_per_parameter_unit','parameter_interval','maximum_full_displacement_A'}
    if any(not isinstance(d,dict) or set(d)!=fields for d in (initial_domain,bottleneck_domain)):
        raise ValueError("only explicit control scale, parameter bounds and displacement radius accepted")
    for d in (initial_domain,bottleneck_domain):
        if (_real(d['scale_A_per_parameter_unit'],'control scale')<=0
                or _real(d['maximum_full_displacement_A'],'displacement radius')<=0):
            raise ValueError("positive control scale and displacement radius required")
        interval=d['parameter_interval']
        if not isinstance(interval,(list,tuple)) or len(interval)!=2:
            raise ValueError("two explicit parameter bounds required")
        lo,hi=(_real(v,'domain bound') for v in interval)
        if not lo<hi or not lo<=contract.anchor_parameter<=hi:
            raise ValueError("ordered domain bounds must include the stationary anchor")
    # Validate global scalar declarations even if every control later abstains.
    x=_real(parameter,'parameter')
    gap=_real(reference_gap_eV_cell,'reference gap')
    tolerance=_real(reference_gap_identity_tolerance_eV_cell,'reference gap identity tolerance')
    if tolerance<0 or abs(bottleneck.reference_potential_eV_cell-initial.reference_potential_eV_cell-gap)>tolerance:
        raise ValueError("recorded raw reference gap contradicts the same-centre energies")
    ids=tuple(sorted(set(initial.measurement_ids)|set(bottleneck.measurement_ids)))
    results=[]
    for pi,ps in zip(nested_stationary_controls(initial),nested_stationary_controls(bottleneck)):
        response=None
        reason=None
        status='available_restricted_quadratic_gap_only'
        if pi.branch is None or ps.branch is None:
            reason='; '.join(f'{role}: {r.status}: {r.reason}' for role,r in (('initial',pi),('bottleneck',ps)) if r.branch is None)
            status='unavailable_control_assumption_or_measurement'
        else:
            mi=ControlledStationaryModel(branch=pi.branch,contract=contract,
                source_audit_sha256=initial.source_audit_sha256,**initial_domain)
            ms=ControlledStationaryModel(branch=ps.branch,contract=contract,
                source_audit_sha256=bottleneck.source_audit_sha256,**bottleneck_domain)
            try:
                response=restricted_stationary_gap_response(mi,ms,x,reference_gap_eV_cell=gap,
                    reference_gap_identity_tolerance_eV_cell=tolerance)
            except ValueError as error:
                status='unavailable_local_domain'
                reason=str(error)
        results.append(NestedGapResult(pi.control,status,reason,pi,ps,response,ids))
    return tuple(results)
