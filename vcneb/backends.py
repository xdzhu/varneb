"""Common backend registry and ASE calculator factories.

The VCNEB core deliberately knows nothing about an external executable.  This
module is the small public boundary between the core and calculator adapters:
each factory receives an image index, an :class:`ase.Atoms` object and a
private image directory, and returns an ASE-compatible calculator exposing
energy, forces and stress.  Existing ``vcneb.vasp``, ``vcneb.abacus`` and
``vcneb.qe`` modules remain available for backwards compatibility; this module
only gives them a uniform discovery and user-facing entry point.

ABINIT, LAMMPS and CP2K are optional ASE integrations.  They are never imported at
module import time, so installing VARNEB does not require either executable.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from pathlib import Path
import shlex
import shutil
import tempfile
from typing import Callable, Mapping

from ase import Atoms
from ase.calculators.calculator import all_changes
from ase.io import write

CalculatorFactory = Callable[[int, Atoms, Path], object]


@dataclass(frozen=True)
class BackendSpec:
    """Static capability metadata used by ``varneb backends`` and docs."""

    name: str
    python_adapter: str
    executable: str
    variable_cell: bool
    status: str
    notes: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


_BACKENDS: tuple[BackendSpec, ...] = (
    BackendSpec(
        "abacus", "vcneb.abacus", "abacus", True, "validated",
        "ASE ABACUS adapter; use static force/stress evaluations per image.",
    ),
    BackendSpec(
        "vasp", "vcneb.vasp", "vasp_std", True, "validated",
        "VASP contract, exact POSCAR lattice round-trip and endpoint caching.",
    ),
    BackendSpec(
        "qe", "vcneb.qe", "pw.x", True, "adapter",
        "ASE Espresso static scf adapter; UPF provenance must be pinned.",
    ),
    BackendSpec(
        "lammps", "vcneb.backends", "lammps", True, "adapter",
        "ASE LAMMPS static force/stress adapter; validate potential and units.",
    ),
    BackendSpec(
        "cp2k", "vcneb.backends", "cp2k_shell", True, "adapter",
        "ASE CP2K static force/stress adapter; validate basis/potential files.",
    ),
    BackendSpec(
        "abinit", "vcneb.backends", "abinit", True, "adapter",
        "ASE ABINIT static force/stress adapter; validate pseudopotential paths and cutoffs.",
    ),
)


def backend_specs() -> tuple[BackendSpec, ...]:
    """Return the immutable list of supported backend specifications."""

    return _BACKENDS


def backend_capability_matrix() -> list[dict[str, object]]:
    """Return JSON-friendly backend metadata without importing executables."""

    return [spec.to_dict() for spec in _BACKENDS]


def get_backend_spec(name: str) -> BackendSpec:
    """Look up a backend by case-insensitive public name."""

    key = str(name).strip().lower()
    for spec in _BACKENDS:
        if spec.name == key:
            return spec
    choices = ", ".join(spec.name for spec in _BACKENDS)
    raise ValueError(f"unknown VARNEB backend {name!r}; choose one of: {choices}")


def attach_image_calculators(
    images: list[Atoms],
    *,
    workdir: str | Path,
    factory: CalculatorFactory,
    structure_filename: str = "structure.start.vasp",
) -> None:
    """Attach calculators in isolated, deterministic image directories.

    Endpoints may be cached by a caller and omitted from execution, but they
    still get a directory and input snapshot when this helper is used.  The
    manager can therefore apply the same provenance and path layout to every
    backend instead of maintaining backend-specific loops.
    """

    root = Path(workdir)
    root.mkdir(parents=True, exist_ok=True)
    for image_index, image in enumerate(images):
        image_dir = root / f"image_{image_index:04d}"
        image_dir.mkdir(parents=True, exist_ok=True)
        write(image_dir / structure_filename, image, format="vasp", direct=True, vasp5=True)
        image.calc = factory(image_index, image, image_dir)


def make_ase_lammps_factory(
    *,
    parameters: Mapping,
    command: str | None = None,
    workdir: str | Path | None = None,
) -> CalculatorFactory:
    """Create an ASE LAMMPS factory with one persistent directory per image.

    ``parameters`` must include a physically meaningful ``pair_style`` and
    ``pair_coeff`` (and, for nontrivial species, ``specorder``/``masses``).
    The adapter does not guess a potential.  ``command`` is normally a module
    provided executable on HF, e.g. ``lmp_mpi``; it is recorded in ASE's
    calculator parameters and can be launched through ``srun --exclusive``.
    """

    supplied = dict(parameters)
    if "pair_style" not in supplied or "pair_coeff" not in supplied:
        raise ValueError("LAMMPS requires explicit pair_style and pair_coeff")
    root = None if workdir is None else Path(workdir)

    def factory(image_index: int, image: Atoms, image_dir: Path):
        try:
            from ase.calculators.lammpsrun import LAMMPS
        except Exception as exc:  # pragma: no cover - optional ASE component
            raise ImportError("ASE LAMMPS support is unavailable") from exc
        params = dict(supplied)
        params.setdefault("keep_tmp_files", True)
        params.setdefault("tmp_dir", str(image_dir if root is None else root / f"image_{image_index:04d}"))
        if command is not None:
            params["command"] = command
        # LAMMPS uses the label as a filename prefix inside ``tmp_dir``;
        # embedding an absolute path would create slash-containing temporary
        # prefixes and fail before the executable is launched.
        return LAMMPS(label="lammps", **params)

    return factory


def make_ase_cp2k_factory(
    *,
    parameters: Mapping,
    command: str | None = None,
) -> CalculatorFactory:
    """Create an ASE CP2K factory with image-isolated output labels.

    CP2K's ASE calculator uses ``label`` as the output prefix rather than a
    separate ``directory`` keyword.  Keeping that prefix inside the image
    directory prevents concurrent VCNEB workers from sharing restart/output
    files.  CP2K input parameters are intentionally passed through unchanged.
    """

    supplied = dict(parameters)

    def output_label(image_dir: Path) -> tuple[str, str]:
        requested = image_dir / "cp2k"
        # CP2K's parser rejects PROJECT paths longer than 80 characters.  A
        # long shared HF path is common, so use a deterministic short local
        # alias and expose the original image directory in the report.  The
        # alias is per image and can be symlinked back for easy inspection.
        if len(str(requested)) <= 72:
            return str(requested), str(image_dir)
        digest = hashlib.sha256(str(image_dir.resolve()).encode()).hexdigest()[:12]
        short_dir = Path(tempfile.gettempdir()) / f"varneb-cp2k-{digest}"
        short_dir.mkdir(parents=True, exist_ok=True)
        for suffix in (".inp", ".out", ".pos"):
            link = image_dir / f"cp2k{suffix}"
            target = short_dir / f"cp2k{suffix}"
            if not link.exists():
                try:
                    link.symlink_to(target)
                except OSError:
                    pass
        return str(short_dir / "cp2k"), str(image_dir)

    def factory(image_index: int, image: Atoms, image_dir: Path):
        try:
            from ase.calculators.cp2k import CP2K
        except Exception as exc:  # pragma: no cover - optional ASE component
            raise ImportError("ASE CP2K support is unavailable") from exc
        image_dir.mkdir(parents=True, exist_ok=True)
        params = dict(supplied)
        params.setdefault("stress_tensor", True)
        label, original_directory = output_label(image_dir)
        class PersistentCP2K(CP2K):
            def _persist_outputs(self) -> None:
                source_prefix = Path(self.label)
                target_prefix = image_dir / "cp2k"
                for suffix in (".inp", ".out", ".pos"):
                    source = Path(f"{source_prefix}{suffix}")
                    target = Path(f"{target_prefix}{suffix}")
                    if not source.is_file():
                        continue
                    try:
                        if target.is_symlink() or target.exists():
                            target.unlink()
                        shutil.copy2(source, target)
                    except OSError:
                        continue

            def calculate(self, atoms=None, properties=None, system_changes=all_changes):
                try:
                    return super().calculate(atoms, properties, system_changes)
                finally:
                    self._persist_outputs()

        calculator = PersistentCP2K(
            label=label,
            command=command,
            **params,
        )
        calculator.varneb_directory = original_directory
        return calculator

    return factory


def make_ase_abinit_factory(
    *,
    parameters: Mapping,
    command: str | None = None,
    pp_paths: str | Path | list[str | Path] | None = None,
) -> CalculatorFactory:
    """Create an ASE ABINIT factory with isolated per-image directories.

    ``command`` is normally the module-provided ``abinit`` executable (or an
    ``srun`` wrapper).  ABINIT pseudopotentials are supplied through ASE's
    ``AbinitProfile`` and are deliberately not guessed by VARNEB.
    """

    supplied = dict(parameters)
    # ASE's ABINIT reader also consumes the eigenvalue sidecar generated by
    # ``prteig``; enabling it avoids a successful SCF being reported as a
    # calculator failure during result parsing.
    supplied.setdefault("prteig", 1)
    if pp_paths is None:
        profile_paths = None
    elif isinstance(pp_paths, (str, Path)):
        profile_paths = [str(pp_paths)]
    else:
        profile_paths = [str(path) for path in pp_paths]

    def factory(image_index: int, image: Atoms, image_dir: Path):
        try:
            from ase.calculators.abinit import Abinit, AbinitProfile
        except Exception as exc:  # pragma: no cover - optional ASE component
            raise ImportError("ASE ABINIT support is unavailable") from exc
        if command is None:
            raise ValueError("ABINIT requires an explicit command or srun wrapper")
        image_dir.mkdir(parents=True, exist_ok=True)
        # ABINIT 8.x consumes a three-line ``.files`` stream on stdin rather
        # than treating the input filename as a positional argument.  ASE's
        # generic profile supplies the filename as an argument, so insert a
        # tiny image-local launcher that bridges the two conventions and keeps
        # the ``.abo``/``_EIG`` outputs where ASE's reader expects them.
        tokens = shlex.split(command)
        if not tokens:
            raise ValueError("ABINIT command must not be empty")
        executable = tokens[-1]
        prefix = tokens[:-1]
        runner = image_dir / "varneb_abinit_runner.sh"
        runner.write_text(
            "#!/bin/sh\n"
            "set -eu\n"
            'input="${1:-abinit.in}"\n'
            'pseudo=$(grep -m1 "^[[:space:]]*pseudos" "$input" | sed \'s/.*pseudos[[:space:]]*"//; s/".*//\')\n'
            "test -n \"$pseudo\"\n"
            'sed \'/^[[:space:]]*pseudos[[:space:]]/d\' "$input" > abinit.varneb.in\n'
            f"printf '%s\\n%s\\n%s\\n%s\\n%s\\n%s\\n' abinit.varneb.in abinit.abo abinit.tmp abinito abinito \"$pseudo\" | exec {shlex.quote(executable)}\n",
            encoding="utf-8",
        )
        runner.chmod(0o755)
        launch = " ".join(shlex.quote(token) for token in [*prefix, str(runner)])
        profile = AbinitProfile(launch, pp_paths=profile_paths)
        calculator = Abinit(profile=profile, directory=str(image_dir), **supplied)
        calculator.varneb_directory = str(image_dir)
        return calculator

    return factory


__all__ = [
    "BackendSpec",
    "CalculatorFactory",
    "attach_image_calculators",
    "backend_capability_matrix",
    "backend_specs",
    "get_backend_spec",
    "make_ase_cp2k_factory",
    "make_ase_abinit_factory",
    "make_ase_lammps_factory",
]
