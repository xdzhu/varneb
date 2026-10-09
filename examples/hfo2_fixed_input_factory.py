"""Case-specific ASE transport preserving historical HfO2 inputs bytewise.

This is a bounded research adapter, not a replacement generic ABACUS backend.
Each call has its own raw SCF directory; endpoint BFGS and VARNEB use the same
energy/forces/stress interface without modifying the physical contract.
"""

from pathlib import Path
import json
import shlex
import shutil
import subprocess
import time

from ase.io import read
from ase import Atoms
from ase.units import Bohr
from ase.calculators.calculator import Calculator, all_changes
import numpy as np

from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_gamma_and_seeds import geometry_roundtrip


CONTRACT = {
    "INPUT": "dc6684ffa709bf3d953c647272589b05d844c1293477bd144226d495207cb9df",
    "KPT": "92b917107d9df11da28a465cc4900f3574956a057ddf9f7ad125e75c956ec508",
    "Hf.upf": "7c516a68310c779a4f572c9ea6b3a73f27ca3a61908b9ade0da3db6db727c8a5",
    "O.upf": "8c6886f443f3bfe7b8d18f2f6b28dd8749746d996499a542fff9666da34d7cf6",
    "Hf_gga_10au_100Ry_4s2p2d1f.orb": "0c72d33ee28f930f3426a0d81255a233c3f9caf1326b5679d861ea6bf83c0519",
    "O_gga_10au_100Ry_2s2p1d.orb": "af216fe56366f36583f7d0d3f6e8820381327ae6f80ce0efc86fde1a1b6bc1d3",
}


def read_fixed_hfo2_stru(path):
    """Read only this case's exact Direct/Hf4O8 writer contract.

    Stock ASE does not necessarily register ABACUS I/O. This deliberately
    narrow reader is not a general STRU parser and never changes coordinates.
    """
    lines = [s.strip() for s in Path(path).read_text(encoding="utf-8").splitlines() if s.strip()]
    expected = {0: "ATOMIC_SPECIES", 1: "Hf 178.49 Hf.upf", 2: "O 15.999 O.upf",
                3: "NUMERICAL_ORBITAL", 4: "Hf_gga_10au_100Ry_4s2p2d1f.orb",
                5: "O_gga_10au_100Ry_2s2p1d.orb", 6: "LATTICE_CONSTANT",
                7: "1.8897261258369282", 8: "LATTICE_VECTORS",
                12: "ATOMIC_POSITIONS", 13: "Direct", 14: "Hf", 15: "0.0",
                16: "4", 21: "O", 22: "0.0", 23: "8"}
    if len(lines) != 32 or any(lines[i] != s for i, s in expected.items()):
        raise ValueError("only the fixed Hf4O8 Direct STRU writer contract is supported")
    lattice = [s.split() for s in lines[9:12]]
    coordinates = [lines[i].split() for i in (*range(17, 21), *range(24, 32))]
    if any(len(row) != 3 for row in lattice) or any(
            len(row) != 6 or row[3:] != ["1", "1", "1"] for row in coordinates):
        raise ValueError("fixed STRU must release all atomic coordinates")
    cell = np.array(lattice, dtype=float) * float(lines[7]) * Bohr
    scaled = np.array([row[:3] for row in coordinates], dtype=float)
    if not np.isfinite(cell).all() or not np.isfinite(scaled).all() or np.linalg.det(cell) <= 0:
        raise ValueError("finite positive-volume fixed STRU required")
    return Atoms(["Hf"]*4+["O"]*8, cell=cell, scaled_positions=scaled, pbc=True)


class FixedHfo2Calculator(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def __init__(self, *, source, command, directory):
        super().__init__(directory=str(directory))
        self.source = Path(source)
        self.command = command
        self.argv = shlex.split(command)
        if self.argv[:3] != ["mpirun", "-np", "32"]:
            raise ValueError("case adapter requires the validated mpirun -np32 command")
        self._check_contract()
        # Every existing failed/completed call remains recoverable on resume.
        calls = [int(p.name.removeprefix("scf_")) for p in Path(self.directory).glob("scf_*")
                 if p.is_dir() and p.name.removeprefix("scf_").isdigit()]
        self.next_call = max(calls, default=-1) + 1

    def _check_contract(self):
        if any(sha256(self.source / name) != expected for name, expected in CONTRACT.items()):
            raise ValueError("HfO2 physical INPUT/KPT/pseudo/orbital contract changed")

    def calculate(self, atoms=None, properties=("energy", "forces", "stress"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self._check_contract()
        if self.atoms.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8:
            raise ValueError("case calculator requires ordered Hf4O8")
        directory = Path(self.directory) / f"scf_{self.next_call:06d}"
        directory.mkdir(parents=True, exist_ok=False)
        self.next_call += 1
        for name in CONTRACT:
            shutil.copyfile(self.source / name, directory / name)
        write_structure(directory / "STRU", self.atoms)
        geometry_roundtrip(directory / "STRU", self.atoms)
        hashes = {name: sha256(directory / name) for name in (*CONTRACT, "STRU")}
        (directory / "input_sha256.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")
        started = time.perf_counter()
        with (directory / "abacus.out").open("w") as stdout, (directory / "abacus.err").open("w") as stderr:
            subprocess.run(self.argv, cwd=directory, stdout=stdout, stderr=stderr, check=True)
        raw = audited_results(directory)
        self.results = raw
        audit = {"elapsed_seconds": time.perf_counter() - started, "input_sha256": hashes,
                 "raw_log_sha256": sha256(directory / "OUT.ABACUS/running_scf.log"),
                 "results": {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in raw.items()}}
        (directory / "call_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")


def make_factory(*, parameters, command):
    if set(parameters) != {"source_directory"}:
        raise ValueError("only source_directory accepted; physical parameter overrides prohibited")
    def factory(image_index, atoms, directory):
        return FixedHfo2Calculator(source=parameters["source_directory"], command=command, directory=directory)
    return factory


def make_endpoint_cached_factory(*, parameters, command):
    """Continue one audited endpoint geometry, not a free-cell/nearest cache.

    Source bytes and parent raw log/STRU are pinned. Movement invalidates ASE's
    cache normally and the same byte-preserving adapter performs the next SCF.
    This restarts BFGS's Hessian; it does not restore optimizer history.
    """
    if set(parameters) != {"source_directory", "seed_static_directory", "seed_input_sha256",
                            "seed_raw_log_sha256"}:
        raise ValueError("endpoint cache requires exact parent hashes; physical overrides prohibited")
    source = Path(parameters["seed_static_directory"])
    expected = parameters["seed_input_sha256"]
    if set(expected) != {*CONTRACT, "STRU"} or any(expected[n] != h for n, h in CONTRACT.items()):
        raise ValueError("endpoint seed physical contract changed")

    def factory(image_index, atoms, directory):
        if type(image_index) is not int or image_index != 0:
            raise ValueError("exactly one endpoint cache, index0, is allowed")
        calc = FixedHfo2Calculator(source=parameters["source_directory"], command=command, directory=directory)
        if ({n: sha256(source/n) for n in expected} != expected
                or sha256(source/"OUT.ABACUS/running_scf.log") != parameters["seed_raw_log_sha256"]):
            raise ValueError("endpoint cached input/raw log changed")
        if not same_ordered_geometry(atoms, read_fixed_hfo2_stru(source/"STRU")):
            raise ValueError("endpoint cache ordered geometry differs")
        raw = audited_results(source)
        calc.atoms = atoms.copy()
        calc.results = {k: v.copy() if isinstance(v, np.ndarray) else v for k, v in raw.items()}
        target = Path(directory)
        if target.exists() and any(p.name != "structure.start.vasp" for p in target.iterdir()):
            raise ValueError("fresh endpoint calculator directory required")
        target.mkdir(parents=True, exist_ok=True)
        record = {"policy": "identical_ordered_clamped_seed_only_no_optimizer_history",
                  "raw_source": str(source), "input_sha256": expected,
                  "raw_log_sha256": parameters["seed_raw_log_sha256"], "new_DFT_calls": 0}
        (target/"seed_cache_audit.json").write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8", newline="\n")
        return calc
    return factory


def same_ordered_geometry(left, right):
    """Periodic equivalence without atom remapping, rotation or force transforms."""
    if (left.get_chemical_symbols() != right.get_chemical_symbols()
            or not np.array_equal(left.pbc, right.pbc)
            or not np.allclose(left.cell.array, right.cell.array, atol=1e-10, rtol=0)):
        return False
    delta = left.get_scaled_positions(wrap=False) - right.get_scaled_positions(wrap=False)
    delta -= np.rint(delta)
    return bool(np.max(np.abs(delta @ left.cell.array)) < 1e-10)


def make_seed_cached_factory(*, parameters, command):
    """Reuse individually audited seed SCFs, especially unchanged endpoints.

    A null directory denotes a new, unevaluated interior geometry. Cached
    results only apply to an identical ordered geometry. ASE invalidates
    them after movement and the byte-preserving calculator then performs a new
    SCF. No nearest-neighbor/permutation cache lookup is performed.
    """
    if set(parameters) != {"source_directory", "seed_static_directories"}:
        raise ValueError("seed factory accepts source_directory and exact ordered seed_static_directories only")
    directories = parameters["seed_static_directories"]
    if not isinstance(directories, list) or len(directories) < 3:
        raise ValueError("at least three ordered seed statics required")

    def factory(image_index, atoms, directory):
        calc = FixedHfo2Calculator(source=parameters["source_directory"], command=command, directory=directory)
        if directories[image_index] is None:
            return calc
        source = Path(directories[image_index])
        if any(sha256(source / n) != h for n, h in CONTRACT.items()):
            raise ValueError("seed cached SCF electronic contract changed")
        if not same_ordered_geometry(atoms, read(source / "STRU", format="abacus")):
            raise ValueError("seed cache ordered geometry differs; refuse stale energy/forces/stress")
        raw = audited_results(source)
        calc.atoms = atoms.copy()
        calc.results = {k: v.copy() if isinstance(v, np.ndarray) else v for k, v in raw.items()}
        target = Path(directory)
        target.mkdir(parents=True, exist_ok=True)
        record = {"policy": "identical_ordered_seed_geometry_only", "raw_source": str(source),
                  "input_sha256": {n: sha256(source / n) for n in (*CONTRACT, "STRU")},
                  "raw_log_sha256": sha256(source / "OUT.ABACUS/running_scf.log")}
        (target / "seed_cache_audit.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        return calc
    return factory
