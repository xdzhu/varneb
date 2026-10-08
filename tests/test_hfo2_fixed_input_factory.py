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


def test_seed_cache_never_silently_relabels_geometry(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    (source / "INPUT").write_text("unchanged fixture input")
    monkeypatch.setattr(transport, "CONTRACT", {"INPUT": sha256(source / "INPUT")})
    atoms = read(Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/reference_variants/T.vasp")
    (source / "STRU").write_text("geometry handled by fixture read")
    (source / "OUT.ABACUS").mkdir()
    (source / "OUT.ABACUS/running_scf.log").write_text("synthetic only")
    monkeypatch.setattr(transport, "read", lambda *a, **k: atoms.copy())
    result = {"energy": -2., "forces": np.zeros((12, 3)), "stress": np.zeros(6)}
    monkeypatch.setattr(transport, "audited_results", lambda path: result)
    monkeypatch.setattr(transport.subprocess, "run", lambda *a, **k: pytest.fail("cached identical seed must not launch DFT"))
    factory = transport.make_seed_cached_factory(parameters={"source_directory": str(source),
                "seed_static_directories": [str(source)] * 3}, command="mpirun -np 32 abacus")
    cached = atoms.copy()
    cached.calc = factory(0, cached, tmp_path / "cached")
    assert cached.get_potential_energy() == -2
    np.testing.assert_array_equal(cached.get_forces(), result["forces"])
    assert (tmp_path / "cached/seed_cache_audit.json").exists()
    uncached_factory = transport.make_seed_cached_factory(parameters={"source_directory": str(source),
                    "seed_static_directories": [str(source), None, str(source)]}, command="mpirun -np 32 abacus")
    assert uncached_factory(1, atoms, tmp_path / "uncached").results == {}
    changed = atoms.copy()
    changed.positions[0, 0] += .001
    with pytest.raises(ValueError, match="geometry differs"):
        factory(1, changed, tmp_path / "changed")


def test_periodic_geometry_cache_is_ordered_not_permutation_invariant():
    atoms = read(Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/reference_variants/T.vasp")
    translated = atoms.copy()
    translated.positions[0] += atoms.cell[0]
    assert transport.same_ordered_geometry(atoms, translated)
    swapped = atoms.copy()
    swapped.positions[[0, 1]] = swapped.positions[[1, 0]]
    assert not transport.same_ordered_geometry(atoms, swapped)
