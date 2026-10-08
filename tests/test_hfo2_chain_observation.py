import json
from pathlib import Path

import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

import scripts.export_hfo2_chain_observation as observer
from scripts.prepare_hfo2_switching_chains import ordered_seed


@pytest.fixture
def evaluated_workdir(tmp_path, monkeypatch):
    variants = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/reference_variants"
    images = ordered_seed(read(variants / "PO.vasp"), read(variants / "PO_minus_T_preserving_inversion.vasp"))
    work = tmp_path / "live"
    snapshot = work / "snapshots/step_0002"
    snapshot.mkdir(parents=True)
    for i, image in enumerate(images):
        image.calc = SinglePointCalculator(image, energy=-10 + .1 * np.sin(i * np.pi / 8),
                                          forces=np.zeros((12, 3)), stress=np.zeros(6))
        write(snapshot / f"POSCAR_{i:02d}", image, format="vasp", direct=True)
        source = work / f"image_{i:04d}/scf_000002"
        (source / "OUT.ABACUS").mkdir(parents=True)
        (source / "OUT.ABACUS/running_scf.log").write_text("synthetic complete SCF")
        (source / "INPUT").write_text("synthetic fixed physical input")
        write(source / "STRU", image, format="vasp", direct=True)
    # Only numeric cached data are mocked here. Actual contract/SCF parsing,
    # geometry matching and raw-log freshness have separate material tests.
    def exact(image, directory):
        i = int(directory.name.split("_")[-1])
        return directory / "scf_000002", dict(images[i].calc.results)

    monkeypatch.setattr(observer, "exact_cached_source", exact)
    monkeypatch.setattr(observer, "CONTRACT", {"INPUT": "synthetic"})
    chain, forces, _, _ = observer.replay(images)
    fmax = np.linalg.norm(forces, axis=1).max()
    (work / "vcneb.opt.log").write_text(f"CheckedFIRE: 2 18:00:00 {chain.enthalpies.max():.6f} {fmax:.6f}\n")
    return work


def test_complete_numeric_export_does_not_mutate_source(evaluated_workdir, tmp_path):
    work = evaluated_workdir
    before = {p.relative_to(work): observer.sha256(p) for p in work.rglob("*") if p.is_file()}
    out = tmp_path / "observation"
    report = observer.export(work, 2, out, "synthetic-no-DFT")
    assert report["status"] == "complete_observation_not_final_result"
    assert report["new_DFT_calls"] == 0
    assert report["n_active_images"] == 7
    assert len(report["extended_reaction_coordinate_A"]) == 9
    assert all(p["input_sha256"]["STRU"] for p in report["raw_image_evaluations"])
    cached = read(out / "evaluated_chain.traj", index=":")
    assert all(isinstance(a.calc, SinglePointCalculator) for a in cached)
    assert {p: observer.sha256(work / p) for p in before} == before
    assert json.loads((out / "observation.json").read_text()) == report
    with pytest.raises(FileExistsError):
        observer.export(work, 2, out, "synthetic")


@pytest.mark.parametrize("fault", ["missing", "wrong_energy", "wrong_force", "duplicate_row", "broken_lift", "wrong_order"])
def test_incomplete_or_inconsistent_observation_fails_closed(evaluated_workdir, tmp_path, fault):
    work = evaluated_workdir
    if fault == "missing":
        (work / "snapshots/step_0002/POSCAR_04").unlink()
    elif fault in ("wrong_energy", "wrong_force", "duplicate_row"):
        log = work / "vcneb.opt.log"
        row = log.read_text()
        if fault == "duplicate_row":
            log.write_text(row + row)
        else:
            fields = row.split()
            fields[3 if fault == "wrong_energy" else 4] = "1.0"
            log.write_text(" ".join(fields) + "\n")
    else:
        path = work / "snapshots/step_0002/POSCAR_04"
        a = read(path)
        if fault == "broken_lift":
            a.positions[0] += a.cell[0]
        else:
            a = a[[4, 1, 2, 3, 0, 5, 6, 7, 8, 9, 10, 11]]
        write(path, a, format="vasp", direct=True)
    out = tmp_path / "bad"
    with pytest.raises(ValueError):
        observer.export(work, 2, out, "synthetic")
    assert not out.exists()


@pytest.mark.parametrize("step", [-1, True, 1.5])
def test_invalid_snapshot_step_rejected(tmp_path, step):
    with pytest.raises(ValueError, match="nonnegative"):
        observer.export(tmp_path, step, tmp_path / "out", "synthetic")
