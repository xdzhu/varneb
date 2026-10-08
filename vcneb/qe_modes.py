"""Read the first, explicitly Gamma, block of QE ``matdyn.modes``.

QE flvec contains normalized Cartesian displacements, not orthonormal
dynamical-matrix eigenvectors. Masses and ordered atoms must come from the
corresponding reference; this reader does not infer a structure, launch QE,
modify a mode, or certify reference stability. Later q blocks are not parsed.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import re

import numpy as np

from .phonons import _THZ_TO_CM1


_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
_Q = re.compile(rf"\s*q\s*=\s*({_NUMBER})\s+({_NUMBER})\s+({_NUMBER})\s*")
_FREQUENCY = re.compile(
    rf"\s*freq\s*\(\s*(\d+)\s*\)\s*=\s*({_NUMBER})\s*\[THz\]"
    rf"\s*=\s*({_NUMBER})\s*\[cm-1\]\s*"
)


def _float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


@dataclass(frozen=True)
class QEGammaModes:
    """Read-only complex columns in the supplied ordered-atom convention.

    ``cartesian_displacements`` are the printed flvec values. Multiplying
    each atom by sqrt(mass), then normalizing each column, produces
    ``mass_weighted_eigenvectors``. The measured Gram defect is retained:
    no QR, ASR, symmetrization or real-phase projection is applied here.
    Frequencies are signed reference frequencies, not a TS stability index.
    """

    frequencies_THz: np.ndarray
    frequencies_cm1: np.ndarray
    cartesian_displacements: np.ndarray
    mass_weighted_eigenvectors: np.ndarray
    masses_amu: np.ndarray
    printed_normalization_max_defect: float
    mass_weighted_gram_max_defect: float
    source_sha256: str


def read_qe_gamma_modes(
    path: str | Path,
    masses_amu,
    *,
    printed_norm_tolerance: float = 5e-5,
    gram_tolerance: float = 5e-4,
) -> QEGammaModes:
    """Read all 3N modes of the first Gamma block, with explicit masses.

    Bounds: 128 atoms, 12MiB input. Reject incomplete/duplicate rows, a first
    nonzero q, inconsistent THz/cm-1, non-unit flvec norms or inconsistent
    mass-weighted orthogonality. Six-decimal printed vectors require a finite
    tolerance; the measured defects, not exact orthogonality, are returned.
    This is the flvec format, not fleig. No nonzero-q unit is guessed.
    """

    masses = np.array(masses_amu, dtype=float, copy=True)
    if (masses.ndim != 1 or not 1 <= len(masses) <= 128
            or not np.isfinite(masses).all() or np.any(masses <= 0)):
        raise ValueError("one positive finite mass per ordered atom required, at most 128 atoms")
    if any(not np.isfinite(x) or not 0 < x < 1
           for x in (printed_norm_tolerance, gram_tolerance)):
        raise ValueError("finite norm/Gram tolerances must be within (0, 1)")
    source = Path(path)
    if source.stat().st_size > 12 * 1024**2:
        raise ValueError("QE mode input exceeds 12MiB bound")
    data = source.read_bytes()
    if len(data) > 12 * 1024**2:
        raise ValueError("QE mode input exceeds 12MiB bound")
    lines = data.decode("ascii").splitlines()
    started = False
    frequencies, wavenumbers, columns = [], [], []
    index = 0
    while index < len(lines):
        line = lines[index]
        q = _Q.fullmatch(line)
        if line.lstrip().startswith("q") and "=" in line and q is None:
            raise ValueError("malformed q header")
        if q is not None:
            if started:
                break
            if not np.allclose([_float(x) for x in q.groups()], 0, rtol=0, atol=1e-12):
                raise ValueError("first mode block is not Gamma")
            started = True
        elif started and line.strip() and set(line.strip()) != {"*"}:
            if (line.strip() == "diagonalizing the dynamical matrix ..."
                    and len(columns) == 3 * len(masses)):
                # QE repeats this separator before each subsequent q header.
                index += 1
                continue
            match = _FREQUENCY.fullmatch(line)
            if match is None or int(match[1]) != len(columns) + 1:
                raise ValueError("sequential frequency records required")
            thz, cm1 = _float(match[2]), _float(match[3])
            if (not np.isfinite([thz, cm1]).all()
                    or not np.isclose(cm1, thz * _THZ_TO_CM1, rtol=1e-7, atol=3e-5)):
                raise ValueError("THz and cm-1 frequency labels disagree")
            rows = []
            for _ in masses:
                index += 1
                if index >= len(lines):
                    raise ValueError("incomplete atom rows in Gamma mode")
                row = lines[index].strip()
                if not (row.startswith("(") and row.endswith(")")):
                    raise ValueError("one six-component row per atom required")
                tokens = row[1:-1].split()
                if len(tokens) != 6 or any(re.fullmatch(_NUMBER, token) is None for token in tokens):
                    raise ValueError("finite six-component displacement row required")
                values = np.array([_float(token) for token in tokens])
                if not np.isfinite(values).all():
                    raise ValueError("non-finite displacement row")
                rows.append(values[0::2] + 1j * values[1::2])
            frequencies.append(thz)
            wavenumbers.append(cm1)
            columns.append(np.asarray(rows).reshape(-1))
            if len(columns) > 3 * len(masses):
                raise ValueError("too many Gamma modes for supplied atom order")
        index += 1
    if not started or len(columns) != 3 * len(masses):
        raise ValueError("complete 3N Gamma mode block required")
    cartesian = np.stack(columns, axis=1)
    defect = float(np.max(np.abs(np.sum(np.abs(cartesian)**2, axis=0) - 1)))
    if defect > printed_norm_tolerance:
        raise ValueError("flvec displacement columns are not normalized")
    root_mass = np.sqrt(masses)
    # A common factor cancels in the column normalization; scale it out to
    # avoid overflow for otherwise finite user-supplied masses.
    weighted = cartesian * np.repeat(root_mass / root_mass.max(), 3)[:, None]
    norms = np.linalg.norm(weighted, axis=0)
    if not np.isfinite(norms).all() or np.any(norms <= 0):
        raise ValueError("mass-weighted displacement norms are not numerically resolvable")
    weighted /= norms
    gram_defect = float(np.max(np.abs(weighted.conj().T @ weighted - np.eye(len(columns)))))
    if gram_defect > gram_tolerance:
        raise ValueError("mass-weighted Gram defect exceeds tolerance; check masses/flvec convention")
    arrays = [np.asarray(frequencies), np.asarray(wavenumbers), cartesian, weighted, masses]
    for array in arrays:
        array.flags.writeable = False
    return QEGammaModes(*arrays, defect, gram_defect, hashlib.sha256(data).hexdigest())
