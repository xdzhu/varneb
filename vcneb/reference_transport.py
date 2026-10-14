"""Declared atomic-reference spans in a local length-valued response chart.

This is a tangent/units/gauge conversion, not a new mode calculation. The
source atoms, vectors and physical gradients are never modified. A basis
change supplies neither curvature nor an accuracy/TS/prediction certificate.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ase import Atoms

from .reference_cell import ReferenceCellCoordinates
from .clamped_mode_coordinates import ClampedModeCoordinates
from .biaxial_mode_coordinates import BiaxialClampedModeCoordinates


def _owned(values):
    result=np.array(values,copy=True)
    result.flags.writeable=False
    return result


@dataclass(frozen=True)
class TransportedAtomicReference:
    """Same internal span, with its source-amplitude map and excluded gauge.

    Input amplitudes use sqrt(amu)*Angstrom. ``source_to_local`` maps them
    to the returned orthonormal length coordinates in Angstrom. Its columns
    are mixed geometric directions, not named phonon eigenvectors. Full raw
    work equals internal work plus the explicit translation work; neither
    part justifies repairing source forces or assuming ASR.
    """
    local_basis: np.ndarray
    raw_mapped_columns: np.ndarray
    excluded_translation_columns: np.ndarray
    source_to_local: np.ndarray
    internal_singular_values: np.ndarray

    @property
    def rank(self):
        return self.local_basis.shape[1]


def transport_atomic_reference_span(
    source_atoms: Atoms, mass_weighted_columns: np.ndarray, masses_amu: np.ndarray,
    target_chart, *, source_to_target_cell_basis: np.ndarray,
    relative_rank_tolerance: float=1e-10,
) -> TransportedAtomicReference:
    """Transport an already registered real atomic span to a chart origin.

    For ASE row cells, P is the *declared* integer cell-basis change, not a
    fitted Cartesian rotation: delta_s_target=delta_s_source@inv(P).
    Thus the Cartesian tangent is (V/sqrt(m))@inv(H_source)@inv(P)@H_target.
    A simultaneous Cartesian rotation is already encoded in H_target.
    Atom indices and periodic lifts are not inferred or rematched. Source
    and target origins may differ; the caller must audit their registration.

    Only atomic columns are transported. All target cell/control components
    stay zero. Restriction to the chart's three rigid translations is an
    explicit analysis-gauge choice, not a change to the raw source vectors.
    Both raw and excluded columns are returned. Any rank loss is refused;
    complete strong references cannot be silently shortened. Complex source
    directions require a previously declared real-span conversion.
    The rank tolerance is algebraic, not a curvature or DFT error threshold.
    """
    if not isinstance(target_chart,(ReferenceCellCoordinates,ClampedModeCoordinates,
                                    BiaxialClampedModeCoordinates)):
        raise TypeError('a supported local response chart is required')
    if not isinstance(source_atoms,Atoms):
        raise TypeError('ordered source Atoms required')
    target=target_chart.reference_atoms
    if (source_atoms.get_chemical_symbols()!=target.get_chemical_symbols()
            or not np.array_equal(source_atoms.pbc,target.pbc)
            or source_atoms.constraints or target.constraints
            or not np.isfinite(source_atoms.positions).all()):
        raise ValueError('same ordered atoms, periodicity and unconstrained source required; no rematching')
    values=np.asarray(mass_weighted_columns)
    if np.iscomplexobj(values):
        raise ValueError('a declared real span is required; do not discard imaginary components')
    values=np.array(values,dtype=float,copy=True)
    masses=np.asarray(masses_amu,dtype=float)
    n=len(source_atoms)
    if (values.ndim!=2 or values.shape[0]!=3*n or not 1<=values.shape[1]<=3*n
            or not np.isfinite(values).all() or masses.shape!=(n,)
            or not np.isfinite(masses).all() or np.any(masses<=0)
            or isinstance(relative_rank_tolerance,(bool,np.bool_))
            or not np.isfinite(relative_rank_tolerance) or not 0<relative_rank_tolerance<1):
        raise ValueError('finite real atomic columns, ordered positive masses and algebraic tolerance required')
    p=np.array(source_to_target_cell_basis,dtype=float,copy=True)
    if (p.shape!=(3,3) or not np.isfinite(p).all() or not np.array_equal(p,np.rint(p))
            or not np.isclose(np.linalg.det(p),1.,atol=1e-12,rtol=0)):
        raise ValueError('explicit right-handed unimodular integer cell basis required')
    source_cell=source_atoms.cell.array
    target_cell=target.cell.array
    if any(not np.isfinite(h).all() or not np.isfinite(np.linalg.det(h))
           or np.linalg.det(h)<=0 or np.linalg.cond(h)>1e12
           for h in (source_cell,target_cell,p)):
        raise ValueError('finite positive well-conditioned declared cells/basis required')
    jacobian=np.linalg.solve(source_cell,np.linalg.solve(p,target_cell))
    raw=np.zeros((target_chart.coordinate_count,values.shape[1]))
    raw[:3*n]=np.einsum('nik,ij->njk',values.reshape(n,3,-1)/np.sqrt(masses)[:,None,None],
                        jacobian).reshape(3*n,-1)
    translations=np.array(target_chart.rigid_translation_directions(),dtype=float,copy=True)
    # The chart returns Cartesian translation directions, not unit columns.
    translations/=np.linalg.norm(translations,axis=0)
    if (translations.shape!=(target_chart.coordinate_count,3)
            or not np.allclose(translations.T@translations,np.eye(3),atol=1e-12,rtol=0)):
        raise ValueError('three independent Cartesian translation directions required')
    excluded=translations@(translations.T@raw)
    internal=raw-excluded
    if not np.isfinite(raw).all():
        raise ValueError('transported reference is not numerically representable')
    u,s,_=np.linalg.svd(internal[:3*n],full_matrices=False)
    if (not np.isfinite(raw).all() or not s.size or s[0]<=0
            or np.count_nonzero(s>relative_rank_tolerance*np.linalg.norm(raw,ord=2))!=values.shape[1]):
        raise ValueError('transport/restriction loses declared rank; no reference direction may be dropped')
    basis=np.zeros_like(raw)
    basis[:3*n]=u[:,:values.shape[1]]
    coefficients=basis.T@internal
    return TransportedAtomicReference(*map(_owned,(basis,raw,excluded,coefficients,s)))
