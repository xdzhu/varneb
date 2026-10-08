import json
from pathlib import Path

from ase.io import read
import numpy as np
import pytest

import examples.hfo2_fixed_input_factory as transport
from scripts.audit_hfo2_static_replica import sha256


def test_transport_keeps_each_call_and_ase_cache(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    (source / "INPUT").write_text("fixed electronic input\n")
    (source / "KPT").write_text("fixed k points\n")
    monkeypatch.setattr(transport, "CONTRACT", {n: sha256(source / n) for n in ("INPUT", "KPT")})
    monkeypatch.setattr(transport, "geometry_roundtrip", lambda path, atoms: 2.)
    launches = []
    monkeypatch.setattr(transport.subprocess, "run", lambda argv, **kwargs: launches.append((argv, kwargs["cwd"])))
    result = {"energy": -1., "forces": np.zeros((12, 3)), "stress": np.zeros(6)}
    def fake_audit(directory):
        (directory / "OUT.ABACUS").mkdir()
        (directory / "OUT.ABACUS/running_scf.log").write_text("synthetic fixture only")
        return result
    monkeypatch.setattr(transport, "audited_results", fake_audit)
    atoms = read(Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/reference_variants/M_seed.vasp")
    atoms.calc = transport.FixedHfo2Calculator(source=source, command="mpirun -np 32 abacus", directory=tmp_path / "calc")
    assert atoms.get_potential_energy() == -1.
    atoms.get_forces()
    atoms.get_stress()
    assert len(launches) == 1
    atoms.positions[0, 0] += .001
    atoms.get_forces()
    assert len(launches) == 2
    assert launches[0][1] != launches[1][1]
    assert sha256(launches[0][1] / "INPUT") == sha256(source / "INPUT")
    restarted = transport.FixedHfo2Calculator(source=source, command="mpirun -np 32 abacus", directory=tmp_path / "calc")
    assert restarted.next_call == 2
    (source / "INPUT").write_text("changed input")
    with pytest.raises(ValueError, match="contract changed"):
        restarted.calculate(atoms)
    assert len(launches) == 2


def test_factory_rejects_hidden_physics_or_mpi_override():
    with pytest.raises(ValueError, match="overrides prohibited"):
        transport.make_factory(parameters={"source_directory": "x", "ecutwfc": 120}, command="mpirun -np 32 abacus")
    with pytest.raises(ValueError, match="mpirun"):
        transport.FixedHfo2Calculator(source="x", command="srun -n 32 abacus", directory="x")
