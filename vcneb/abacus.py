"""ABACUS calculator adapter helpers for VC-NEB.

ASE's ABACUS calculator is version/distribution dependent, so this module uses
late imports and accepts either an installed ASE ABACUS calculator or a user
provided calculator factory.
"""

from __future__ import annotations

from pathlib import Path
import math
import re
from typing import Callable, Mapping, Optional

from ase import Atoms
import numpy as np
from ase.io import write


CalculatorFactory = Callable[[int, Atoms, Path], object]

_ABACUS_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


def _prefer_machine_precision_energy(output: Path, parsed_energy: float) -> float:
    """Use ABACUS's full-precision final marker when it matches ASE's energy.

    The human-readable ``final etot is`` line commonly rounds to seven
    decimals. Small frozen-mode curvature differences can be of that order.
    A disagreement larger than print rounding is an audit failure, not a
    reason to silently replace a different energy definition.
    """

    parsed = float(parsed_energy)
    if not math.isfinite(parsed):
        raise ValueError(f"ABACUS parsed energy is non-finite: {output}")
    marker = None
    pattern = re.compile(rf"!FINAL_ETOT_IS\s+({_ABACUS_FLOAT})\s*eV", re.IGNORECASE)
    with output.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = pattern.search(line)
            if match is not None:
                marker = float(match.group(1).replace("D", "E").replace("d", "e"))
    if marker is None:
        return parsed
    if not math.isfinite(marker) or abs(marker - parsed) > 1e-4:
        raise ValueError(f"ABACUS full-precision energy disagrees with parsed result: {output}")
    return marker


def _minimal_abacus_results(output: Path) -> dict:
    """Parse the static force/stress contract from a damaged ABACUS log."""

    text = output.read_text(encoding="utf-8", errors="replace")
    # ABACUS writes both the human-readable ``final etot is`` line and the
    # machine-oriented ``!FINAL_ETOT_IS`` marker.  A long SCF summary can
    # occasionally wrap/truncate the former even though the latter and the
    # force/stress contract are complete.
    energy_patterns = (
        rf"final\s+etot\s+is\s+({_ABACUS_FLOAT})\s*eV",
        rf"!FINAL_ETOT_IS\s+({_ABACUS_FLOAT})\s*eV",
    )
    energy_matches = [
        match.group(1)
        for pattern in energy_patterns
        for match in re.finditer(pattern, text, flags=re.IGNORECASE)
    ]
    if not energy_matches:
        raise ValueError(f"ABACUS log has no final total energy: {output}")
    energy = float(energy_matches[-1].replace("D", "E").replace("d", "e"))

    force_body = text.rsplit("TOTAL-FORCE", 1)[-1]
    force_body = force_body.split("TOTAL-STRESS", 1)[0]
    forces: list[list[float]] = []
    for line in force_body.splitlines():
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            values = [float(value.replace("D", "E").replace("d", "e")) for value in fields[-3:]]
        except ValueError:
            continue
        if len(fields) == 4 or re.match(r"^[A-Za-z][A-Za-z0-9_]*$", fields[0]):
            forces.append(values)
    if not forces:
        raise ValueError(f"ABACUS log has no parseable final force block: {output}")

    stress_body = text.rsplit("TOTAL-STRESS", 1)[-1]
    stress_rows: list[list[float]] = []
    for line in stress_body.splitlines():
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            values = [float(value.replace("D", "E").replace("d", "e")) for value in fields[-3:]]
        except ValueError:
            continue
        stress_rows.append(values)
        if len(stress_rows) == 3:
            break
    if len(stress_rows) != 3:
        raise ValueError(f"ABACUS log has no parseable final stress block: {output}")

    from ase.stress import full_3x3_to_voigt_6_stress

    # ABACUS prints compressive-positive kbar; ASE expects tensile-positive
    # eV/A^3. 1 eV/A^3 = 1602.176634 kbar. Multiplication by a GPa factor
    # here would silently corrupt the fallback stress by orders of magnitude.
    stress = -np.asarray(stress_rows, dtype=float) / 1602.176634
    return {
        "energy": energy,
        "forces": np.asarray(forces, dtype=float),
        "stress": full_3x3_to_voigt_6_stress(stress),
    }


REQUIRED_VCNEB_PARAMETERS = {
    "cal_force": 1,
    "cal_stress": 1,
    "out_stru": 1,
}


def _reject_broken_mpi_launch(directory: Path) -> None:
    """Reject ABACUS output produced after Intel MPI fell back to singletons.

    A direct Slurm/PMIx launch can print this warning once per requested rank,
    then let every independent process write apparently usable energy and
    force blocks to the same image directory.  Parsing those blocks as one
    valid distributed calculation would silently contaminate the NEB chain.
    """

    for name in ("abacus.out", "abacus.err"):
        path = directory / name
        if not path.is_file():
            continue
        with path.open(encoding="utf-8", errors="replace") as handle:
            if any("PMI server not found" in line for line in handle):
                raise RuntimeError(
                    f"ABACUS MPI launch failed in {directory}: {path} reports "
                    "'PMI server not found'; discard this image calculation "
                    "and use a validated MPI launcher"
                )


def _read_vcneb_results(directory: Path, *, output_suffix: str, calculation: str) -> dict:
    """Read only the ABACUS properties required by variable-cell NEB.

    ASE's ABACUS reader eagerly parses eigenvalues while building its complete
    result dictionary.  Some ABACUS builds write a harmless but irregular
    k-point eigenvalue block, which makes that optional parser fail before the
    energy, forces, and stress are returned.  VCNEB does not use eigenvalues,
    so keep this adapter on the smaller, explicit calculator contract.
    """

    directory = Path(directory)
    _reject_broken_mpi_launch(directory)
    output = directory / ("OUT." + output_suffix) / f"running_{calculation}.log"
    try:
        from ase.io.abacus import _get_abacus_chunks
    except ImportError:  # pragma: no cover - version-specific ASE fallback
        from ase.io.abacus import read_abacus_results

        with output.open(encoding="utf-8") as handle:
            return read_abacus_results(handle, index=-1)[0]

    try:
        with output.open(encoding="utf-8") as handle:
            chunk = _get_abacus_chunks(handle, index=-1, non_convergence_ok=False)[0]
    except (AttributeError, IndexError, KeyError, TypeError, ValueError) as exc:
        # ABACUS can leave a usable final force/stress block while its verbose
        # header or optional eigenvalue block is malformed.  Keep the original
        # parser error in the chained exception if the contract fallback also
        # fails, but do not discard valid VCNEB data merely because eigenvalue
        # diagnostics are unreadable.
        try:
            return _minimal_abacus_results(output)
        except Exception as fallback_exc:
            raise ValueError(
                f"ABACUS result parsing failed for {output}: {exc}; "
                f"minimal contract parser also failed: {fallback_exc}"
            ) from exc
    values = {
        "energy": (_prefer_machine_precision_energy(output, chunk.energy)
                   if chunk.energy is not None else None),
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
