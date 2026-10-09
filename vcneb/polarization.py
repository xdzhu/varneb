"""Explicit polarization-branch bookkeeping, not a spontaneous-P estimator.

ABACUS logs can report a spin-paired modulo 2eR/V rather than eR/V.
Keep that reported modulus; never silently identify it with the physical
quantum or choose a branch merely by taking the smallest absolute value.
"""

from __future__ import annotations

import re
import io
import numpy as np


def quantum_lattice(cell_A):
    """Row lattice e*a_i/V in C/m^2 for row cell vectors in Angstrom."""
    cell = np.asarray(cell_A, dtype=float)
    if cell.shape != (3, 3) or not np.isfinite(cell).all():
        raise ValueError("finite 3x3 cell required")
    volume = float(np.linalg.det(cell))
    if volume <= 0 or np.linalg.cond(cell) > 1e8:
        raise ValueError("positive, nondegenerate right-handed cell required")
    return 16.02176634 * cell / volume


def modular_difference(value, reference, reported_modulus):
    """Signed nearest scalar difference modulo an explicitly supplied period.

    This does not select a physical path branch. At a half-period tie two
    choices are equally close; callers must not use this to unwrap a path.
    """
    values = np.asarray([value, reference, reported_modulus], dtype=float)
    if not np.isfinite(values).all() or reported_modulus <= 0:
        raise ValueError("finite values and positive explicit modulus required")
    delta = float(value - reference)
    return delta - reported_modulus * np.rint(delta / reported_modulus)


def sampled_reduced_branch(values_C_m2, quanta_C_m2, reported_periods_C_m2,
                           uncertainties_C_m2, *, initial_integer=0):
    """A conditional nearest-neighbour lift in the *actual-cell* reduced units.

    Keep the reported (possibly spin-paired) period, not an assumed eR/V.
    The first integer is an explicit gauge, not a spontaneous-P reference.
    Abstain on half-period ties or uncertainty-overlapping alternatives.
    Even a unique sampled lift cannot exclude winding between sampled images
    or certify full-BZ insulation. This routine does not close those gates.
    """
    arrays = [np.asarray(x, dtype=float) for x in
              (values_C_m2, quanta_C_m2, reported_periods_C_m2, uncertainties_C_m2)]
    if (any(x.ndim != 1 for x in arrays) or len(arrays[0]) < 2
            or any(x.shape != arrays[0].shape or not np.isfinite(x).all() for x in arrays)
            or np.any(arrays[1] <= 0) or np.any(arrays[2] <= 0) or np.any(arrays[3] < 0)
            or isinstance(initial_integer, bool) or not isinstance(initial_integer, int)):
        raise ValueError("aligned finite samples, positive quanta/periods and integer gauge required")
    values, quanta, periods, errors = arrays
    ratios = periods / quanta
    period = float(ratios[0])
    if not np.allclose(ratios, period, rtol=0, atol=2e-6):
        raise ValueError("reported reduced period changes along the path")
    reduced = values / quanta
    uncertainty = errors / quanta
    lifted = [float(reduced[0] + initial_integer * period)]
    integers = [initial_integer]
    links = []
    for i in range(1, len(values)):
        integer = int(np.rint((lifted[-1] - reduced[i]) / period))
        candidate = float(reduced[i] + integer * period)
        delta = candidate - lifted[-1]
        # This is a measured quadrature-sensitivity margin, not a rigorous
        # total DFT error bar or a bound on unsampled polarization winding.
        margin = period / 2 - abs(delta) - uncertainty[i] - uncertainty[i-1]
        unique = bool(margin > 1e-8)
        links.append({"left": i-1, "right": i, "delta_reduced": delta,
                      "half_period_margin_reduced": float(margin), "unique": unique})
        if not unique:
            return {"status": "ambiguous_sampled_link", "ambiguous_link": [i-1, i],
                    "reported_reduced_period": period, "links": links,
                    "lifted_reduced": None, "branch_integers": None,
                    "continuous_path_certified": False, "spontaneous_P_selected": False}
        lifted.append(candidate)
        integers.append(integer)
    return {"status": "conditional_sampled_lift", "reported_reduced_period": period,
            "initial_branch_integer": initial_integer, "raw_reduced": reduced.tolist(),
            "lifted_reduced": lifted, "branch_integers": integers, "links": links,
            "delta_reduced": lifted[-1] - lifted[0],
            "continuous_path_certified": False, "spontaneous_P_selected": False}


def parse_abacus_berry(body):
    """Read exactly one direction's native C/m^2 output without branch edits."""
    directions = re.findall(r"calculated polarization direction is in R([123]) direction", body)
    number = r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?"
    pattern = (rf"P\s*=\s*({number})\s*\(mod\s*({number})\)\s*"
               rf"\(\s*({number})\s*,\s*({number})\s*,\s*({number})\s*\)\s*C/m\^2")
    entries = re.findall(pattern, body)
    if len(directions) != 1 or len(entries) != 1:
        raise ValueError("expected exactly one complete Berry direction/C/m^2 result")
    values = np.asarray(entries[0], dtype=float)
    if not np.isfinite(values).all() or values[1] <= 0:
        raise ValueError("invalid Berry value or reported modulus")
    if not np.isclose(np.linalg.norm(values[2:]), abs(values[0]), atol=3e-7, rtol=0):
        raise ValueError("scalar Berry value and Cartesian projection disagree")
    return {"direction": int(directions[0]), "value_C_m2": float(values[0]),
            "reported_modulus_C_m2": float(values[1]),
            "cartesian_projection_C_m2": values[2:].tolist(),
            "branch_selected": False}


def sampled_band_gap(body, occupied_bands):
    """Indirect sampled gap from ABACUS istate.info, not a full BZ certificate.

    The occupation column is k-weighted in ABACUS; do not interpret its
    numeric value as a per-state electron occupation. The caller must fix
    the occupied-band count from the audited electron/pseudopotential count.
    """
    if not isinstance(occupied_bands, int) or occupied_bands < 1:
        raise ValueError("positive integer occupied-band count required")
    blocks = []
    current = None
    for line in body.splitlines():
        if "Kpoint =" in line and "BAND" in line:
            current = []
            blocks.append(current)
        elif current is not None and line.strip():
            parts = line.split()
            if len(parts) < 3:
                raise ValueError("malformed eigenvalue row")
            index, energy, occupation = int(parts[0]), float(parts[1]), float(parts[2])
            if index != len(current) + 1 or not np.isfinite([energy, occupation]).all():
                raise ValueError("nonfinite or unordered band table")
            current.append(energy)
    if not blocks or any(len(b) <= occupied_bands for b in blocks):
        raise ValueError("missing occupied/unoccupied states")
    if any(np.min(np.diff(b)) < -1e-7 for b in blocks):
        raise ValueError("eigenvalues are not ordered")
    vbm = max(b[occupied_bands - 1] for b in blocks)
    cbm = min(b[occupied_bands] for b in blocks)
    return {"n_kpoints": len(blocks), "occupied_bands": occupied_bands,
            "VBM_eV": vbm, "CBM_eV": cbm, "sampled_indirect_gap_eV": cbm - vbm}


def sampled_nscf_band_gap(body, occupied_bands):
    """Native BANDS_1.dat: k index, accumulated k distance, band energies(eV).

    ABACUS f7cb1d3 LCAO writes istate.info for SCF/MD/relax only. NSCF
    eigenvalues require explicit out_band=1; see module_io/nscf_band.cpp.
    """
    if not isinstance(occupied_bands, int) or occupied_bands < 1:
        raise ValueError("positive integer occupied-band count required")
    if not body.strip():
        raise ValueError("empty NSCF band table")
    try:
        table = np.loadtxt(io.StringIO(body), ndmin=2)
    except (ValueError, TypeError) as error:
        raise ValueError("malformed NSCF band table") from error
    if (table.shape[0] < 1 or table.shape[1] < occupied_bands + 3
            or not np.isfinite(table).all()
            or not np.array_equal(table[:, 0], np.arange(1, len(table) + 1))
            or np.min(np.diff(table[:, 2:], axis=1)) < -1e-7):
        raise ValueError("missing, nonfinite, incomplete or unordered NSCF bands")
    vbm = float(table[:, 1 + occupied_bands].max())
    cbm = float(table[:, 2 + occupied_bands].min())
    return {"n_kpoints": len(table), "occupied_bands": occupied_bands,
            "VBM_eV": vbm, "CBM_eV": cbm, "sampled_indirect_gap_eV": cbm - vbm,
            "source_format": "ABACUS_BANDS_1.dat"}
