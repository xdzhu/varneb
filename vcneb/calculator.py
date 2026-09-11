"""Calculator capability checks and diagnostics for VC-NEB.

The VC-NEB core only needs the ASE calculator protocol.  This module provides
an explicit preflight check so a missing stress property is reported before an
optimization starts instead of being silently interpreted as zero cell force.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

from ase import Atoms


def _declared_properties(calculator: object) -> tuple[str, ...]:
    values = getattr(calculator, "implemented_properties", ())
    try:
        return tuple(sorted({str(value).lower() for value in values}))
    except TypeError:
        return ()


def _command_text(calculator: object) -> str | None:
    command = getattr(calculator, "command", None)
    if command is None:
        profile = getattr(calculator, "profile", None)
        command = getattr(profile, "command", None)
    if command is None:
        return None
    if isinstance(command, (list, tuple)):
        return " ".join(str(part) for part in command)
    return str(command)


@dataclass(frozen=True)
class CalculatorCapabilities:
    """Machine-readable capability report for one ASE calculator."""

    calculator_type: str
    module: str
    declared_properties: tuple[str, ...]
    has_energy: bool
    has_forces: bool
    has_stress: bool
    variable_cell: bool
    has_directory: bool
    directory: str | None
    command: str | None
    issues: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return not self.issues

    def to_dict(self) -> dict:
        data = asdict(self)
        data["ok"] = self.ok
        return data


class CalculatorCapabilityError(RuntimeError):
    """Raised when one or more calculators cannot support VC-NEB."""

    def __init__(self, message: str, reports: Sequence[CalculatorCapabilities] = ()):
        super().__init__(message)
        self.reports = tuple(reports)


def inspect_calculator(
    calculator: object,
    *,
    require_stress: bool = True,
    require_variable_cell: bool = True,
    require_directory: bool = False,
) -> CalculatorCapabilities:
    """Inspect an ASE-like calculator without executing an expensive calculation.

    A calculator with no ``implemented_properties`` declaration is accepted
    when its protocol methods are callable.  If the declaration is present,
    ``stress`` must be listed explicitly for variable-cell use.
    """

    declared = _declared_properties(calculator)
    has_energy = callable(getattr(calculator, "get_potential_energy", None))
    has_forces = callable(getattr(calculator, "get_forces", None))
    has_stress_method = callable(getattr(calculator, "get_stress", None))
    stress_declared = not declared or "stress" in declared
    has_stress = has_stress_method and stress_declared

    raw_directory = getattr(calculator, "directory", None)
    directory = None if raw_directory is None else str(Path(raw_directory))
    has_directory = directory is not None
    issues: list[str] = []
    if not has_energy:
        issues.append("calculator does not provide get_potential_energy()")
    if not has_forces:
        issues.append("calculator does not provide get_forces()")
    if require_stress and not has_stress:
        if not has_stress_method:
            issues.append("calculator does not provide get_stress(); VCNEB requires stress")
        elif not stress_declared:
            issues.append("calculator does not declare stress in implemented_properties")
    if require_variable_cell and not has_stress:
        issues.append("calculator cannot be verified for variable-cell force evaluation")
    if require_directory and not has_directory:
        issues.append("calculator has no per-image directory")

    return CalculatorCapabilities(
        calculator_type=type(calculator).__name__,
        module=type(calculator).__module__,
        declared_properties=declared,
        has_energy=has_energy,
        has_forces=has_forces,
        has_stress=has_stress,
        variable_cell=has_stress,
        has_directory=has_directory,
        directory=directory,
        command=_command_text(calculator),
        issues=tuple(issues),
    )


def validate_calculator(
    calculator: object,
    *,
    require_stress: bool = True,
    require_variable_cell: bool = True,
    require_directory: bool = False,
) -> CalculatorCapabilities:
    """Validate one calculator and return its capability report."""

    report = inspect_calculator(
        calculator,
        require_stress=require_stress,
        require_variable_cell=require_variable_cell,
        require_directory=require_directory,
    )
    if not report.ok:
        details = "; ".join(report.issues)
        raise CalculatorCapabilityError(
            f"{report.calculator_type} cannot be used for VC-NEB: {details}",
            [report],
        )
    return report


def validate_image_calculators(
    images: Iterable[Atoms],
    *,
    require_stress: bool = True,
    require_variable_cell: bool = True,
    require_directory: bool = False,
    require_unique_directories: bool = False,
) -> list[CalculatorCapabilities]:
    """Validate all image calculators and include image-local diagnostics."""

    reports: list[CalculatorCapabilities] = []
    failures: list[str] = []
    directories: dict[str, int] = {}
    for image_index, image in enumerate(images):
        calculator = getattr(image, "calc", None)
        if calculator is None:
            failures.append(f"image {image_index}: no calculator is attached")
            continue
        report = inspect_calculator(
            calculator,
            require_stress=require_stress,
            require_variable_cell=require_variable_cell,
            require_directory=require_directory,
        )
        reports.append(report)
        if report.issues:
            context = f"image {image_index}"
            if report.directory:
                context += f" ({report.directory})"
            failures.append(f"{context}: {'; '.join(report.issues)}")
        if require_unique_directories and report.directory is not None:
            previous = directories.get(report.directory)
            if previous is not None:
                failures.append(
                    f"image {image_index} and image {previous} share calculator directory {report.directory}"
                )
            else:
                directories[report.directory] = image_index

    if failures:
        detail = "\n".join(f"- {failure}" for failure in failures)
        raise CalculatorCapabilityError(f"VC-NEB calculator preflight failed:\n{detail}", reports)
    return reports


def calculator_context(calculator: object) -> str:
    """Return concise directory/command context for a runtime error."""

    report = inspect_calculator(
        calculator,
        require_stress=False,
        require_variable_cell=False,
        require_directory=False,
    )
    details = [f"calculator={report.calculator_type}"]
    if report.directory is not None:
        details.append(f"directory={report.directory}")
    if report.command is not None:
        details.append(f"command={report.command}")
    return ", ".join(details)


__all__ = [
    "CalculatorCapabilities",
    "CalculatorCapabilityError",
    "calculator_context",
    "inspect_calculator",
    "validate_calculator",
    "validate_image_calculators",
]
