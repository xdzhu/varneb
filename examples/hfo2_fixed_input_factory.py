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
