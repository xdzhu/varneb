"""Optional non-mutating lattice preflight using a separately licensed harness.

No VASP implementation or binary is bundled. A consistent classification is
not an electronic SCF result, nor a guarantee that VASP will accept every input.
The diagnostic wrapper in scripts/vasp_lattice_probe.f90 uses a user's installed
Intel-classic lattlib.o; it must first be checked against that VASP executable.
"""
from dataclasses import dataclass, asdict
from hashlib import sha256
from pathlib import Path
import subprocess

import numpy as np

from .step_control import CandidateStepRejected
from .vasp_contract import (
    VaspInputContractError, poscar_lattice_roundtrip, validate_vasp_image_geometry,
)


@dataclass(frozen=True)
class LatticeClassification:
    real_type: int
    reciprocal_type: int
    expected_reciprocal_type: int
    cell_sha256: str

    @property
    def consistent(self):
        return self.reciprocal_type == self.expected_reciprocal_type

    def to_dict(self):
        return dict(asdict(self), consistent=self.consistent)


class NativeVaspLatticeProbe:
    """Call a fixed executable and tolerance; never alter atoms or run DFT."""

    def __init__(self, executable, *, symprec=1e-4, timeout=10):
        self.executable = Path(executable).resolve(strict=True)
        if not self.executable.is_file():
            raise ValueError("lattice probe executable must be a file")
        self.executable_sha256 = sha256(self.executable.read_bytes()).hexdigest()
        self.symprec = float(symprec)
        self.timeout = float(timeout)
        if not np.isfinite(self.symprec) or self.symprec <= 0:
            raise ValueError("symprec must be finite and positive")
        if not np.isfinite(self.timeout) or self.timeout <= 0:
            raise ValueError("timeout must be finite and positive")
        self._cache = {}

    def classify(self, cell):
        if sha256(self.executable.read_bytes()).hexdigest() != self.executable_sha256:
            raise RuntimeError("lattice probe executable changed after its contract was recorded")
        cell = np.asarray(cell, dtype=float)
        if cell.shape != (3, 3) or not np.isfinite(cell).all() or np.linalg.det(cell) <= 0:
            raise ValueError("finite, positive-volume 3x3 cell required")
        encoded = (format(self.symprec, ".17g") + "\n" + "\n".join(
            " ".join(format(float(value), ".17g") for value in row) for row in cell
        ) + "\n")
        digest = sha256(encoded.encode("ascii")).hexdigest()
        if digest in self._cache:
            return self._cache[digest]
        result = subprocess.run([str(self.executable)], input=encoded, text=True,
                                capture_output=True, timeout=self.timeout, check=True)
        fields = result.stdout.split()
        if len(fields) != 3 or not all(field.isdigit() for field in fields):
            raise RuntimeError("invalid lattice probe output: expected exactly three type indices")
        indices = tuple(int(field) for field in fields)
        if any(value < 1 or value > 14 for value in indices):
            raise RuntimeError("lattice probe returned an unknown Bravais type")
        report = LatticeClassification(*indices, cell_sha256=digest)
        self._cache[digest] = report
        return report

    def descriptor(self):
        return {"executable": str(self.executable), "sha256": self.executable_sha256,
                "symprec": self.symprec, "timeout_seconds": self.timeout,
                "electronic_result": False, "geometry_modified": False}


class NativeVaspCandidateValidator:
    """Check all candidate geometries and writer-rounded lattices without DFT."""

    def __init__(self, probe, *, minimum_distance=None, maximum_deformation=None):
        if not callable(getattr(probe, "classify", None)):
            raise TypeError("probe must provide classify(cell)")
        for value in (minimum_distance, maximum_deformation):
            if value is not None and (not np.isfinite(value) or value <= 0):
                raise ValueError("geometry limits must be finite and positive")
        self.probe = probe
        self.minimum_distance = minimum_distance
        self.maximum_deformation = maximum_deformation

    def __call__(self, images):
        from .core import validate_path_geometry
        try:
            for image in images:
                validate_vasp_image_geometry(image, minimum_distance=self.minimum_distance)
            validate_path_geometry(images, maximum_deformation=self.maximum_deformation)
        except (VaspInputContractError, ValueError) as error:
            raise CandidateStepRejected(str(error), details=[{"category": "candidate_geometry"}]) from error
        failures = []
        for index, image in enumerate(images):
            # Match the exact-round-trip policy used by ExplicitPotcarVasp.
            # The preview is separate and never replaces manager atoms.
            serialized_cell = poscar_lattice_roundtrip(image.cell.array)
            for representation, cell in (("manager", image.cell.array),
                                         ("POSCAR_preview", serialized_cell)):
                report = self.probe.classify(cell)
                if not report.consistent:
                    failures.append({"image_index": index, "representation": representation,
                                     "category": "vasp_bravais_lattice_inconsistency", **report.to_dict()})
        if failures:
            raise CandidateStepRejected("Candidate lattice classification is inconsistent",
                                        details=failures)
