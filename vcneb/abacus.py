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


def _read_vcneb_results(directory: Path, *, output_suffix: str, calculation: str) -> dict:
    """Read only the ABACUS properties required by variable-cell NEB.

    ASE's ABACUS reader eagerly parses eigenvalues while building its complete
    result dictionary.  Some ABACUS builds write a harmless but irregular
    k-point eigenvalue block, which makes that optional parser fail before the
    energy, forces, and stress are returned.  VCNEB does not use eigenvalues,
    so keep this adapter on the smaller, explicit calculator contract.
    """

    output = Path(directory) / ("OUT." + output_suffix) / f"running_{calculation}.log"
    try:
        from ase.io.abacus import _get_abacus_chunks
    except ImportError:  # pragma: no cover - version-specific ASE fallback
        from ase.io.abacus import read_abacus_results

        with output.open(encoding="utf-8") as handle:
            return read_abacus_results(handle, index=-1)[0]

    with output.open(encoding="utf-8") as handle:
        chunk = _get_abacus_chunks(handle, index=-1, non_convergence_ok=False)[0]
    values = {
        "energy": chunk.energy,
        "free_energy": chunk.free_energy,
        "forces": chunk.forces_sort,
        "stress": chunk.stress,
        "magmom": chunk.magmom,
        "dipole": chunk.dipole,
    }
    return {key: value for key, value in values.items() if value is not None}


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
            calculator = Abacus(directory=str(image_dir), profile=profile, **params)
        elif command is not None:
            try:
                from ase.calculators.abacus import AbacusProfile
            except Exception as exc:  # pragma: no cover - depends on external ASE build
                raise ImportError(
                    "This ASE ABACUS adapter does not expose AbacusProfile; "
                    "pass an explicit profile object instead of command."
                ) from exc
            calculator = Abacus(
                directory=str(image_dir), profile=AbacusProfile(command), **params
            )
        else:
            calculator = Abacus(directory=str(image_dir), **params)

        template = getattr(calculator, "template", None)
        if template is not None:
            template.read_results = lambda directory: _read_vcneb_results(
                directory,
                output_suffix=template.out_suffix,
                calculation=template.cal_name,
            )
        return calculator

    return factory
