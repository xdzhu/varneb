"""Load one audited BTO cubic-reference Gamma-mode coordinate contract.

Both frozen and conditionally relaxed examples consume the same reference,
atom mapping, mode ordering and hashes. No calculator is constructed here.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import (
    ModePlane, PhonopyGammaEigenpairs, ReferenceCellCoordinates,
    anchored_gamma_axis_weights,
    gamma_modes_from_phonopy_eigenpairs, load_gamma_force_constants,
    load_phonopy_gamma_eigenpairs,
)
from vcneb.phonons import GammaModes


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_report_hash(report: dict, source: Path) -> str:
    digest = sha256(source)
    # Reports created on Windows retain backslash-separated path keys when
    # copied to a Linux cluster; pathlib on Linux does not parse those keys.
    matches = [value for key, value in report.get("input_sha256", {}).items()
               if str(key).replace("\\", "/").rsplit("/", 1)[-1] == source.name]
    if len(matches) != 1 or matches[0] != digest:
        raise ValueError(f"source hash is absent or mismatched in projection report: {source}")
    return digest


def validate_bto_gamma_source(
    provenance_path: Path,
    force_sets_path: Path,
    eigenpairs_provenance_path: Path,
) -> dict[str, str]:
    """Gate the archived BTO Gamma basis against its actual displacement data.

    The 1x1x1 phonon supercell is independent of the 4x4x4 electronic k mesh.
    This case-specific check must not be applied to non-Gamma phonon work.
    """

    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    eigenpairs_provenance = json.loads(eigenpairs_provenance_path.read_text(encoding="utf-8"))
    calculator = provenance.get("calculator", {})
    force_constants = provenance.get("force_constants", {})
    if (provenance.get("endpoint") != "cubic BaTiO3"
            or provenance.get("supercell") != [1, 1, 1]
            or provenance.get("n_signed_displacements") != 6
            or not np.isclose(provenance.get("finite_difference_amplitude_A", -1), 0.01)
            or calculator.get("code") != "ABACUS"
            or calculator.get("functional") != "PBE"
            or calculator.get("ecutwfc_Ry") != 100
            or "10 au" not in calculator.get("basis", "")
            or calculator.get("kpoints") != [4, 4, 4]
            or force_constants.get("raw_force_sets") != force_sets_path.name
            or force_constants.get("phonopy_gamma_eigenpairs_provenance") != eigenpairs_provenance_path.name):
        raise ValueError("archived BTO Gamma source violates the 1x1x1/100Ry/10au contract")
    with force_sets_path.open(encoding="utf-8") as source:
        header = [line.strip() for line in source if line.strip()][:2]
    if header != ["5", "6"]:
        raise ValueError("BTO Gamma FORCE_SETS must contain five atoms and six displacements")
    if (eigenpairs_provenance.get("force_sets_sha256") != sha256(force_sets_path)
            or eigenpairs_provenance.get("method")
            != "Phonopy run_qpoints([[0, 0, 0]], with_eigenvectors=True)"
            or eigenpairs_provenance.get("n_modes") != 15):
        raise ValueError("BTO Gamma eigenpair provenance is not tied to the archived FORCE_SETS")
    return {
        "gamma_provenance": sha256(provenance_path),
        "force_sets": sha256(force_sets_path),
        "eigenpairs_provenance": sha256(eigenpairs_provenance_path),
    }


@dataclass(frozen=True)
class BtoQ1Q2Reference:
    report: dict
    source_hashes: dict[str, str]
    chart: ReferenceCellCoordinates
    modes: GammaModes
    plane: ModePlane


def load_bto_q1q2_reference(
    report_path: Path,
    reference_path: Path,
    force_constants_path: Path,
    phonopy_eigenpairs_path: Path,
    *,
    strain_metric_weights_amu_A2: np.ndarray | None = None,
) -> BtoQ1Q2Reference:
    """Rebuild exactly the signed axes used by the archived path report."""

    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("kind") != "posthoc_q1q2_projection_not_PES":
        raise ValueError("report must be an audited posthoc Q1/Q2 projection")
    hashes = {
        "report": sha256(report_path),
        "reference": _require_report_hash(report, reference_path),
        "force_constants": _require_report_hash(report, force_constants_path),
        "phonopy_eigenpairs": _require_report_hash(report, phonopy_eigenpairs_path),
    }
    if report.get("reference_id") != f"sha256:{hashes['phonopy_eigenpairs']}":
        raise ValueError("reference_id differs from the reviewed eigenpair hash")
    reference = read(str(reference_path))
    translation = np.asarray(report.get("reference_fractional_translation"), dtype=float)
    permutation = report.get("reference_atom_permutation_in_path_order")
    if translation.shape != (3,) or not np.all(np.isfinite(translation)):
        raise ValueError("reference translation is missing or invalid")
    if permutation is None or sorted(permutation) != list(range(len(reference))):
        raise ValueError("reference permutation is invalid")
    reference.set_scaled_positions(reference.get_scaled_positions(wrap=False) + translation)
    reference = reference[permutation]
    chart = ReferenceCellCoordinates(reference)
    _, masses = load_gamma_force_constants(force_constants_path)
    if len(masses) != len(reference):
        raise ValueError("force-constant masses do not match the reference structure")
    dofs = [3 * atom + axis for atom in permutation for axis in range(3)]
    raw = load_phonopy_gamma_eigenpairs(phonopy_eigenpairs_path)
    reordered = PhonopyGammaEigenpairs(
        frequencies_thz=raw.frequencies_thz,
        eigenvectors_mass_weighted=raw.eigenvectors_mass_weighted[dofs, :],
    )
    modes = gamma_modes_from_phonopy_eigenpairs(reordered, masses[permutation])
    axis_weights = np.asarray(report["axis_mode_weights"], dtype=float)
    common = {
        "axis_labels": tuple(report["axis_labels"]),
        "reference_id": report["reference_id"],
    }
    if strain_metric_weights_amu_A2 is None:
        plane = ModePlane.from_gamma_modes(modes, axis_weights, **common)
    else:
        plane = ModePlane.from_gamma_modes_with_strain(
            modes, axis_weights,
            strain_metric_weights_amu_A2=strain_metric_weights_amu_A2,
            **common,
        )
    if plane.amplitude_unit != report.get("amplitude_unit"):
        raise ValueError("mode amplitude units differ from the archived report")
    return BtoQ1Q2Reference(report, hashes, chart, modes, plane)


def bto_transverse_soft_plane(
    reference: BtoQ1Q2Reference,
    *,
    strain_metric_weights_amu_A2: np.ndarray | None = None,
) -> ModePlane:
    """Orient the complete cubic unstable triplet with Ti--Ba z/x anchors.

    The resulting fixed-C-cell atomic plane is a *different* coordinate
    contract from the archived soft/stable Q1/Q2 plane. Existing off-axis
    DFT energies cannot be relabeled as samples of this new plane.
    """

    modes = reference.modes
    if reference.chart.reference_atoms.get_chemical_symbols() != ["Ba", "Ti", "O", "O", "O"]:
        raise ValueError("unexpected cubic BTO atom order")
    frequencies = modes.frequencies_cm1[:3]
    if np.any(frequencies >= 0.0) or np.ptp(frequencies) > 1e-4:
        raise ValueError("first three Gamma modes are not a degenerate unstable triplet")
    anchors = np.zeros((3 * modes.n_atoms, 2))
    anchors[5, 0], anchors[2, 0] = 1.0, -1.0  # Ti_z - Ba_z
    anchors[3, 1], anchors[0, 1] = 1.0, -1.0  # Ti_x - Ba_x
    weights = anchored_gamma_axis_weights(modes, np.arange(3), anchors)
    common = {
        "axis_labels": ("Q_parallel: cubic z polar soft mode", "Q_transverse: cubic x polar soft mode"),
        "reference_id": reference.report["reference_id"],
    }
    if strain_metric_weights_amu_A2 is None:
        return ModePlane.from_gamma_modes(modes, weights, **common)
    return ModePlane.from_gamma_modes_with_strain(
        modes, weights,
        strain_metric_weights_amu_A2=strain_metric_weights_amu_A2,
        **common,
    )


__all__ = [
    "BtoQ1Q2Reference", "bto_transverse_soft_plane", "load_bto_q1q2_reference",
    "sha256", "validate_bto_gamma_source",
]
