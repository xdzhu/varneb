from types import SimpleNamespace
import subprocess

import numpy as np
import pytest

from vcneb.vasp_lattice import NativeVaspLatticeProbe


def make_probe(tmp_path):
    executable = tmp_path / "diagnostic"
    executable.write_bytes(b"test fixture: not an executable")
    return NativeVaspLatticeProbe(executable)


def test_reports_conflict_without_mutation_and_caches(tmp_path, monkeypatch):
    probe = make_probe(tmp_path)
    calls = []
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: calls.append(kw) or SimpleNamespace(stdout="4 11 4\n"))
    cell = np.eye(3)
    original = cell.copy()
    report = probe.classify(cell)
    assert not report.consistent
    assert report.to_dict()["expected_reciprocal_type"] == 4
    assert probe.classify(cell) is report
    assert len(calls) == 1
    assert calls[0]["check"] is True
    assert np.array_equal(cell, original)
    assert probe.descriptor()["electronic_result"] is False


@pytest.mark.parametrize("output", ["4 4", "4 4 4 garbage", "-1 4 4", "15 4 4", "4.0 4 4"])
def test_invalid_output_is_fatal(tmp_path, monkeypatch, output):
    probe = make_probe(tmp_path)
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: SimpleNamespace(stdout=output))
    with pytest.raises(RuntimeError):
        probe.classify(np.eye(3))


@pytest.mark.parametrize("cell", [np.zeros((3, 3)), np.diag([-1, 1, 1]), np.full((3, 3), np.nan), np.eye(2)])
def test_invalid_cells_do_not_start_process(tmp_path, monkeypatch, cell):
    probe = make_probe(tmp_path)
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: pytest.fail("must not launch"))
    with pytest.raises(ValueError):
        probe.classify(cell)


def test_timeout_does_not_become_success_or_cached_result(tmp_path, monkeypatch):
    probe = make_probe(tmp_path)
    def fail(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], kwargs["timeout"])
    monkeypatch.setattr(subprocess, "run", fail)
    with pytest.raises(subprocess.TimeoutExpired):
        probe.classify(np.eye(3))
    assert not probe._cache


def test_changed_executable_invalidates_even_cached_classification(tmp_path, monkeypatch):
    probe = make_probe(tmp_path)
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: SimpleNamespace(stdout="4 4 4"))
    assert probe.classify(np.eye(3)).consistent
    probe.executable.write_bytes(b"different diagnostic")
    with pytest.raises(RuntimeError, match="changed"):
        probe.classify(np.eye(3))
