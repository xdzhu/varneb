"""Synthetic parser checks, not a physical hafnia benchmark."""

from pathlib import Path

import pytest

from scripts.audit_hfo2_static_replica import INPUT_FILES, audit_replica


def fixture_image(path: Path, *, delta: float = 0.0, ranks: int = 32):
    path.mkdir()
    for name in INPUT_FILES:
        (path / name).write_text(name, encoding="utf-8")
    output = path / "OUT.ABACUS"
    output.mkdir()
    body = f"DSIZE = {ranks}\ncharge density convergence is achieved\n"
    body += "TOTAL-FORCE (eV/Angstrom)\n"
    body += "".join(f"Hf{i} 0.01 0.02 -0.01\n" for i in range(12))
    body += "TOTAL-STRESS (KBAR)\n0.1 0 0\n0 0.2 0\n0 0 0.3\n"
    body += f"!FINAL_ETOT_IS {-100 + delta:.10f} eV\n"
    (output / "running_scf.log").write_text(body, encoding="utf-8")


def test_identical_contract_and_full_results(tmp_path):
    source, replica = tmp_path / "source", tmp_path / "replica"
    fixture_image(source)
    fixture_image(replica, delta=1e-6)
    report = audit_replica(source, replica)
    assert report["status"] == "passed"
    assert report["energy_delta_meV_fu"] == pytest.approx(0.00025)
    assert len(report["results"]["forces"]) == 12
    assert report["results"]["stress"][0] == pytest.approx(-0.1 / 1602.176634)


@pytest.mark.parametrize("failure", ["inputs", "ranks", "singleton", "scf"])
def test_invalid_evidence_is_rejected(tmp_path, failure):
    source, replica = tmp_path / "source", tmp_path / "replica"
    fixture_image(source)
    fixture_image(replica, ranks=1 if failure == "ranks" else 32)
    if failure == "inputs":
        (replica / "INPUT").write_text("ecutwfc 120", encoding="utf-8")
    elif failure == "singleton":
        (replica / "abacus.err").write_text("PMI server not found", encoding="utf-8")
    elif failure == "scf":
        log = replica / "OUT.ABACUS/running_scf.log"
        log.write_text(log.read_text().replace("charge density convergence is achieved", "not converged"))
    with pytest.raises(ValueError):
        audit_replica(source, replica)


def test_valid_but_changed_results_require_review(tmp_path):
    source, replica = tmp_path / "source", tmp_path / "replica"
    fixture_image(source)
    fixture_image(replica, delta=0.001)
    assert audit_replica(source, replica)["status"] == "review_required"
