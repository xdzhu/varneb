"""ABACUS calculator adapter helpers for VC-NEB.

ASE's ABACUS calculator is version/distribution dependent, so this module uses
late imports and accepts either an installed ASE ABACUS calculator or a user
provided calculator factory.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Mapping, Optional

from ase import Atoms
from ase.io import write


CalculatorFactory = Callable[[int, Atoms, Path], object]


REQUIRED_VCNEB_PARAMETERS = {
    "cal_force": 1,
    "cal_stress": 1,
    "out_stru": 1,
}


def _install_parallel_safe_sort_writer() -> None:
    """Avoid ASE-ABACUS's process-global ``ase_sort.dat`` race for identity order.

    Recent ASE ABACUS adapters write ``ase_sort.dat`` in the Python process's
    current directory even though each calculator has an isolated directory.
    Concurrent image workers can therefore truncate one another's sort file.
    Our production fixtures keep atoms grouped by species, so the sort is the
    identity permutation and the file is unnecessary: the ABACUS force order
    already matches the ASE order.  Preserve the upstream writer for genuinely
    non-identity permutations (serial execution remains supported).
    """

    try:
        import ase.io.abacus as abacus_io
    except Exception:  # pragma: no cover - optional external adapter
        return
    if getattr(abacus_io, "_vcneb_sort_writer_patched", False):
        return
    original = abacus_io.write_input_stru_sort

    def safe_write_input_stru_sort(atoms_sort=None):
        if atoms_sort is not None:
            values = [int(value) for value in atoms_sort]
            if values == list(range(len(values))):
                return None
        return original(atoms_sort)

    abacus_io.write_input_stru_sort = safe_write_input_stru_sort
    abacus_io._vcneb_sort_writer_patched = True


_install_parallel_safe_sort_writer()


def attach_abacus_calculators(
    images: list[Atoms],
    *,
    workdir: str | Path,
    factory: CalculatorFactory,
) -> None:
    root = Path(workdir)
    root.mkdir(parents=True, exist_ok=True)
    for image_index, image in enumerate(images):
        image_dir = root / f"{image_index:02d}"
        image_dir.mkdir(parents=True, exist_ok=True)
        write(image_dir / "STRU.start.vasp", image, format="vasp", direct=True, vasp5=True)
        image.calc = factory(image_index, image, image_dir)


def make_ase_abacus_factory(
    *,
    parameters: Mapping,
    command: Optional[str] = None,
    profile: object = None,
    **kwargs,
) -> CalculatorFactory:
    """Create a factory for ASE builds that provide ``ase.calculators.abacus``."""

    if command is not None and profile is not None:
        raise ValueError("Pass either command or profile, not both")

    def factory(image_index: int, image: Atoms, image_dir: Path):
        try:
            from ase.calculators.abacus import Abacus
        except Exception as exc:  # pragma: no cover - depends on external ASE build
            raise ImportError(
                "This ASE installation does not provide ase.calculators.abacus. "
                "Install an ASE/ABACUS adapter or pass a custom factory."
            ) from exc

        params = dict(parameters)
        params.update(kwargs)
        params.update(REQUIRED_VCNEB_PARAMETERS)
        if profile is not None:
            return Abacus(directory=str(image_dir), profile=profile, **params)
        if command is not None:
            try:
                from ase.calculators.abacus import AbacusProfile
            except Exception as exc:  # pragma: no cover - depends on external ASE build
                raise ImportError(
                    "This ASE ABACUS adapter does not expose AbacusProfile; "
                    "pass an explicit profile object instead of command."
                ) from exc
            return Abacus(directory=str(image_dir), profile=AbacusProfile(command), **params)
        return Abacus(directory=str(image_dir), **params)

    return factory
