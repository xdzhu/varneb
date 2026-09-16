"""Quantum ESPRESSO calculator helpers for variable-cell NEB.

The VCNEB core sees only ASE energy, force, and stress calls.  This module
keeps the Quantum ESPRESSO-specific static ``pw.x`` setup at the adapter
boundary so a path image cannot accidentally invoke QE's own ionic/cell
optimizer.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
from pathlib import Path
import re
from typing import Callable, Mapping

from ase import Atoms
from ase.io import write


CalculatorFactory = Callable[[int, Atoms, Path], object]

_REQUIRED_CONTROL = {
    "calculation": "scf",
    "tstress": True,
    "tprnfor": True,
}
_UPF_ELEMENT_RE = re.compile(r"\belement\s*=\s*['\"]?([A-Za-z]{1,3})", re.IGNORECASE)


def validate_qe_pseudopotentials(
    pseudo_dir: str | Path,
    pseudopotentials: Mapping[str, str],
) -> list[dict]:
    """Validate and fingerprint QE UPF files before a VCNEB preflight.

    This is deliberately a provenance gate rather than a claim that the
    selected potentials are converged.  It rejects missing files, a filename
    that escapes ``pseudo_dir``, inconsistent element metadata, and files
    without a visible PBE marker.  The complete SHA256 is retained in the run
    report so a reviewed selection cannot silently change during a restart.
    """

    root = Path(pseudo_dir).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"QE pseudopotential directory does not exist: {root}")
    if not pseudopotentials:
        raise ValueError("at least one QE pseudopotential mapping is required")
    reports: list[dict] = []
    for species, filename in sorted(pseudopotentials.items()):
        label = str(species).strip()
        relative = Path(str(filename))
        if not label or relative.is_absolute() or relative.name != str(filename):
            raise ValueError(f"QE pseudopotential for {species!r} must be a filename inside {root}")
        if relative.suffix.lower() != ".upf":
            raise ValueError(f"QE pseudopotential for {label} is not a .UPF file: {filename!r}")
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(f"QE pseudopotential for {label} is missing: {path}")
        header = path.read_bytes()[:65536].decode("utf-8", errors="ignore")
        elements = {value.lower() for value in _UPF_ELEMENT_RE.findall(header)}
        if label.lower() not in elements:
            detail = ", ".join(sorted(elements)) if elements else "none"
            raise ValueError(f"QE UPF {path.name} element metadata {detail!r} does not match {label}")
        if "pbe" not in header.lower():
            raise ValueError(f"QE UPF {path.name} has no visible PBE functional marker")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        reports.append(
            {
                "species": label,
                "filename": path.name,
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": digest,
                "element_metadata": sorted(elements),
                "functional_marker": "PBE",
            }
        )
    return reports


def static_qe_input_data(input_data: Mapping | None = None) -> dict:
    """Return QE input data validated for a static VCNEB image evaluation.

    QE may perform ionic or variable-cell updates internally for calculations
    such as ``relax`` and ``vc-relax``.  Those updates conflict with the
    manager-owned VCNEB coordinates, so only ``scf`` is accepted here.  Stress
    and forces are mandatory because VCNEB transforms both into generalized
    forces.
    """

    data = deepcopy(dict(input_data or {}))
    control = data.get("control", {})
    if not isinstance(control, Mapping):
        raise ValueError("QE input_data['control'] must be a mapping")
    control = dict(control)

    for key, required in _REQUIRED_CONTROL.items():
        supplied = control.get(key, data.get(key, required))
        if key == "calculation":
            if str(supplied).strip().lower() != "scf":
                raise ValueError(
                    "VCNEB QE images must use static calculation='scf'; "
                    "do not use relax or vc-relax"
                )
        elif supplied is not True:
            raise ValueError(f"VCNEB QE images require {key}=True")
        control[key] = required
        data.pop(key, None)
    data["control"] = control
    return data


def attach_qe_calculators(
    images: list[Atoms],
    *,
    workdir: str | Path,
    factory: CalculatorFactory,
) -> None:
    """Attach one QE calculator in one immutable directory per path image."""

    root = Path(workdir)
    root.mkdir(parents=True, exist_ok=True)
    for image_index, image in enumerate(images):
        image_dir = root / f"{image_index:02d}"
        image_dir.mkdir(parents=True, exist_ok=True)
        write(image_dir / "POSCAR.start", image, format="vasp", direct=True, vasp5=True)
        image.calc = factory(image_index, image, image_dir)


def make_ase_espresso_factory(
    *,
    parameters: Mapping,
    command: str | None = None,
    pseudo_dir: str | Path | None = None,
    profile: object | None = None,
) -> CalculatorFactory:
    """Build an ASE ``Espresso`` factory with VCNEB-safe image parameters.

    Pass a prebuilt ASE ``EspressoProfile`` through ``profile`` or let this
    helper create one from ``command`` and ``pseudo_dir``.  One of these
    explicit launch paths is required to keep the calculator invocation in the
    run manifest rather than relying on a process-global environment variable.
    """

    if command is not None and profile is not None:
        raise ValueError("Pass either command or profile, not both")
    if profile is None and (command is None or pseudo_dir is None):
        raise ValueError("QE requires either profile or both command and pseudo_dir")

    supplied = dict(parameters)
    supplied["input_data"] = static_qe_input_data(supplied.get("input_data"))

    def factory(image_index: int, image: Atoms, image_dir: Path):
        try:
            from ase.calculators.espresso import Espresso, EspressoProfile
        except Exception as exc:  # pragma: no cover - optional ASE installation details
            raise ImportError(
                "This ASE installation does not provide ase.calculators.espresso. "
                "Install an ASE build with Quantum ESPRESSO support."
            ) from exc

        active_profile = profile
        if active_profile is None:
            active_profile = EspressoProfile(command, str(pseudo_dir))
        return Espresso(profile=active_profile, directory=str(image_dir), **deepcopy(supplied))

    return factory


__all__ = [
    "CalculatorFactory",
    "attach_qe_calculators",
    "make_ase_espresso_factory",
    "static_qe_input_data",
    "validate_qe_pseudopotentials",
]
