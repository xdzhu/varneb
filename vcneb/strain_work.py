"""Energy/enthalpy directional work with ASE row cells and tensile stress.

This is a local gradient contraction, not a relaxed-branch derivative unless
the released variables are stationary. No calculator or optimizer is called.
"""
from __future__ import annotations

from numbers import Real
import numpy as np


def configuration_work(cell, stress, cell_direction, *, forces=None,
                       fractional_direction=None, pressure_eV_A3=0.0):
    """Return directional d(E+PV)/dt for H(t), s(t), r=sH.

    H and dH/dt are row matrices in Angstrom; ASE stress is a symmetric
    tensile-positive 3x3 tensor in eV/Angstrom^3. Forces are eV/Angstrom;
    ds/dt is dimensionless per unit t. The cell part is
    V (sigma+P I) : [(dH/dt)^T H^-T], at fixed fractional coordinates.
    The additional internal-coordinate part is -sum F_i dot [(ds_i/dt)H].
    Components are eV per unit t, without mass weighting or a cell-scale factor.
    Nonzero P is a fixed external pressure, not a strain-dependent pressure.
    """
    h, sigma, dh = [np.asarray(a, dtype=float) for a in (cell, stress, cell_direction)]
    if any(a.shape != (3,3) or not np.isfinite(a).all() for a in (h,sigma,dh)):
        raise ValueError('finite 3x3 cell, stress and cell direction required')
    if np.linalg.det(h) <= 1e-12 or np.linalg.cond(h) > 1e12:
        raise ValueError('positive-volume nonsingular cell required')
    if not np.allclose(sigma,sigma.T,atol=1e-12,rtol=1e-10):
        raise ValueError('symmetric ASE tensile stress required')
    if (isinstance(pressure_eV_A3,(bool,np.bool_)) or not isinstance(pressure_eV_A3,Real)
            or not np.isfinite(pressure_eV_A3)):
        raise ValueError('finite explicit pressure required')
    volume = float(np.linalg.det(h))
    velocity = dh.T @ np.linalg.inv(h.T)
    cell_work = float(volume*np.sum((sigma+pressure_eV_A3*np.eye(3))*velocity))
    atomic_work = 0.0
    if (forces is None) != (fractional_direction is None):
        raise ValueError('forces and fractional direction must be supplied together')
    if forces is not None:
        f, ds = np.asarray(forces,dtype=float), np.asarray(fractional_direction,dtype=float)
        if (f.ndim != 2 or f.shape[1] != 3 or f.shape[0] == 0 or ds.shape != f.shape
                or not np.isfinite(f).all() or not np.isfinite(ds).all()):
            raise ValueError('matching finite nonempty Nx3 forces and fractional direction required')
        atomic_work = float(-np.sum(f*(ds@h)))
    return dict(cell_eV_per_control=cell_work, atomic_eV_per_control=atomic_work,
                total_eV_per_control=cell_work+atomic_work)
