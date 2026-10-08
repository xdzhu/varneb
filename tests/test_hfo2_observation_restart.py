import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np
import pytest
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator

import scripts.prepare_hfo2_observation_restart as restart


@pytest.fixture
def stopped_observation(tmp_path, monkeypatch):
    case = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/chain_observations/PO_M_lifted_step10"
    observation = tmp_path / "observation"
    shutil.copytree(case, observation)
    work = tmp_path / "terminal_work"
    snapshot = work / "snapshots/step_0010"
    snapshot.mkdir(parents=True)
    report = json.loads((observation / "observation.json").read_text())
    report["source_workdir"] = str(work)
    records = {}
    for i, point in enumerate(report["raw_image_evaluations"]):
        poscar = observation / f"POSCAR_{i:02d}"
        shutil.copyfile(poscar, snapshot / poscar.name)
        source = work / f"image_{i:04d}/scf_000010"
        (source / "OUT.ABACUS").mkdir(parents=True)
        (source / "INPUT").write_text("synthetic fixed contract")
        (source / "OUT.ABACUS/running_scf.log").write_text("synthetic complete output")
        write(source / "STRU", read(poscar), format="vasp", direct=True)
        point["raw_source"] = str(source)
        point["input_sha256"] = {n: restart.sha256(source / n) for n in ("INPUT", "STRU")}
        point["raw_log_sha256"] = restart.sha256(source / "OUT.ABACUS/running_scf.log")
        records[source] = {"energy": point["energy_eV_cell"],
                           "forces": np.array(point["forces_eV_A"]),
                           "stress": np.array(point["stress_ASE_voigt_eV_A3"])}
    (observation / "observation.json").write_text(json.dumps(report))
    monkeypatch.setattr(restart, "CONTRACT", {"INPUT": report["raw_image_evaluations"][0]["input_sha256"]["INPUT"]})
    monkeypatch.setattr(restart, "audited_results", lambda p: records[p])
    real_read = restart.read
    monkeypatch.setattr(restart, "read", lambda p, **kw: real_read(p, **{**kw, "format": "vasp"})
                        if kw.get("format") == "abacus" else real_read(p, **kw))
    return observation, work, report


def terminal(job_id):
    return {"job_id": job_id, "state": "COMPLETED", "exit_code": "0:0", "end": "2026-10-08T18:42:50"}


def test_already_converged_residual_is_not_resubmitted(stopped_observation, tmp_path, monkeypatch):
    observation, _, report = stopped_observation
    report["replayed_fmax_eV_A"] = .05
    (observation / "observation.json").write_text(json.dumps(report))
    monkeypatch.setattr(restart, "replay", lambda images: (None, np.array([[.05, 0, 0]]), None, {}))
    with pytest.raises(ValueError, match="already converged"):
        restart.prepare(observation, tmp_path / "unnecessary", job_reader=terminal)
    assert not (tmp_path / "unnecessary").exists()


def test_exact_restart_preserves_geometry_cache_and_source(stopped_observation, tmp_path):
    observation, work, original = stopped_observation
    before = {p: restart.sha256(p) for p in work.rglob("*") if p.is_file()}
    output = tmp_path / "restart"
    manifest = restart.prepare(observation, output, job_reader=terminal)
    assert manifest["new_DFT_calls"] == 0 and manifest["all_seed_SCFs_reused"]
    assert not manifest["physical_inputs_changed"] and not manifest["old_optimizer_state_reused"]
    assert manifest["next_optimizer_step_cap"] == 80 and manifest["wall_cap_hours"] == 24
    assert manifest["restart_fmax_eV_A"] == pytest.approx(.2933253052250083, abs=1e-10)
    images = read(output / "seed.traj", index=":")
    assert all(a.calc is None for a in images)
    for a, old in zip(images, read(observation / "evaluated_chain.traj", index=":")):
        assert np.allclose(a.get_scaled_positions(wrap=False), old.get_scaled_positions(wrap=False), atol=1e-14)
    parameters = json.loads((output / "factory_parameters.json").read_text())
    assert parameters["seed_static_directories"] == [p["raw_source"] for p in original["raw_image_evaluations"]]
    assert all(restart.sha256(p) == h for p, h in before.items())
    with pytest.raises(FileExistsError):
        restart.prepare(observation, output, job_reader=terminal)


@pytest.mark.parametrize("fault", ["running", "newer_snapshot", "input", "raw_log", "raw_structure", "old_snapshot_changed", "raw_energy"])
def test_unsafe_or_stale_restart_fails_before_output(stopped_observation, tmp_path, monkeypatch, fault):
    observation, work, report = stopped_observation
    point = report["raw_image_evaluations"][3]
    source = Path(point["raw_source"])
    reader = terminal
    if fault == "running":
        reader = lambda job_id: {"job_id": job_id, "state": "RUNNING"}
    elif fault == "newer_snapshot":
        shutil.copytree(work / "snapshots/step_0010", work / "snapshots/step_0011")
    elif fault in ("input", "raw_log", "old_snapshot_changed"):
        p = {"input": source / "INPUT", "raw_log": source / "OUT.ABACUS/running_scf.log",
             "old_snapshot_changed": work / "snapshots/step_0010/POSCAR_03"}[fault]
        p.write_text(p.read_text() + "changed")
    elif fault == "raw_structure":
        a = read(source / "STRU", format="vasp")
        a.positions[0, 0] += .01
        write(source / "STRU", a, format="vasp", direct=True)
    else:
        raw_reader = restart.audited_results
        monkeypatch.setattr(restart, "audited_results", lambda p: {**raw_reader(p), "energy": raw_reader(p)["energy"] + .001})
    with pytest.raises(ValueError):
        restart.prepare(observation, tmp_path / "bad", job_reader=reader)
    assert not (tmp_path / "bad").exists()


@pytest.mark.parametrize("stdout,valid", [
    ("123|RUNNING|0:0|Unknown\n123.batch|COMPLETED|0:0|now\n", False),
    ("123|COMPLETED|0:0|2026-10-08T18:42:50\n123.batch|COMPLETED|0:0|now\n", True),
    ("123|TIMEOUT|0:0|2026-10-08T18:42:50\n", True),
    ("124|COMPLETED|0:0|now\n", False),
])
def test_sacct_main_allocation_is_authoritative_not_a_completed_step(monkeypatch, stdout, valid):
    monkeypatch.setattr(restart.subprocess, "run", lambda *a, **kw: SimpleNamespace(stdout=stdout))
    if valid:
        assert restart.terminal_job_record("123")["state"] in restart.TERMINAL
    else:
        with pytest.raises(ValueError):
            restart.terminal_job_record("123")


@pytest.mark.parametrize("case,count,expected", [
    ("gap", 10, .28499524288774647),
    ("PO_M", 9, .2933253052250083),
])
def test_published_terminal_restart_replays_without_remote_DFT(case, count, expected):
    """Check delivered numeric evidence; live raw-SCF verification is separate."""
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/ordinary_continuation"
    seed = root / f"{case}_restart"
    manifest = json.loads((seed / "manifest.json").read_text())
    preflight = json.loads((root / f"{case}_cache_preflight.json").read_text())
    assert manifest["source_job"]["state"] == "COMPLETED"
    assert manifest["source_snapshot_step"] == 10
    assert manifest["n_total_images"] == count
    assert manifest["all_seed_SCFs_reused"] and manifest["new_DFT_calls"] == 0
    assert not manifest["climb"] and not manifest["physical_inputs_changed"]
    assert not manifest["old_optimizer_state_reused"]
    assert manifest["next_optimizer_step_cap"] == 80 and manifest["wall_cap_hours"] == 24
    assert preflight["status"] == "passed_actual_production_factory_without_DFT"
    assert preflight["cached_evaluations"] == count and preflight["new_SCF_calls"] == 0
    assert preflight["seed_manifest_sha256"] == restart.sha256(seed / "manifest.json")
    assert preflight["seed_traj_sha256"] == restart.sha256(seed / "seed.traj")
    assert all(restart.sha256(seed / name) == value
               for name, value in manifest["seed_file_sha256"].items())
    images = read(seed / "seed.traj", index=":")
    assert len(images) == count and all(a.calc is None for a in images)
    for i, (image, point) in enumerate(zip(images, manifest["raw_image_evaluations"])):
        assert point["image_index"] == i and point["ordered_periodic_geometry_matched"]
        assert image.get_chemical_symbols() == ["Hf"] * 4 + ["O"] * 8
        assert all(point["input_sha256"][name] == value for name, value in restart.CONTRACT.items())
        image.calc = SinglePointCalculator(image, energy=point["energy_eV_cell"],
                                          forces=np.array(point["forces_eV_A"]),
                                          stress=np.array(point["stress_ASE_voigt_eV_A3"]))
    _, forces, _, _ = restart.replay(images)
    assert np.linalg.norm(forces, axis=1).max() == pytest.approx(expected, abs=1e-10)
    assert preflight["fmax_eV_A"] == pytest.approx(expected, abs=1e-10)
    assert manifest["restart_fmax_eV_A"] == pytest.approx(expected, abs=1e-10)


@pytest.mark.parametrize("case,job_id,expected", [
    ("PO_minus_T_preserving", "28288045", .980656798544163),
    ("PO_minus_T_reversing", "28288063", .7002788564042406),
])
def test_terminal_switching_restart_replays_complete_fixed_contract_cache(case, job_id, expected):
    from scripts.analyze_hfo2_channel_network import read_evaluated_observation

    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"
    seed = root / "switching_continuation" / f"{case}_restart"
    manifest = json.loads((seed / "manifest.json").read_text())
    preflight = json.loads((seed.parent / f"{case}_cache_preflight.json").read_text())
    observation, report, _ = read_evaluated_observation(root / "chain_observations" / f"{case}_step10")
    assert manifest["source_job"]["job_id"] == job_id == report["source_job_id"]
    assert manifest["source_job"]["state"] == "COMPLETED"
    assert manifest["source_snapshot_step"] == report["snapshot_step"] == 10
    assert manifest["n_total_images"] == 9 and manifest["n_active_images"] == 7
    assert manifest["new_DFT_calls"] == preflight["new_SCF_calls"] == 0
    assert manifest["all_seed_SCFs_reused"] and preflight["cached_evaluations"] == 9
    assert not manifest["climb"] and not manifest["physical_inputs_changed"]
    assert not manifest["old_optimizer_state_reused"] and manifest["fmax_eV_A"] == .10
    assert manifest["next_optimizer_step_cap"] == 80 and manifest["wall_cap_hours"] == 24
    assert preflight["status"] == "passed_actual_production_factory_without_DFT"
    assert preflight["seed_manifest_sha256"] == restart.sha256(seed / "manifest.json")
    assert preflight["seed_traj_sha256"] == restart.sha256(seed / "seed.traj")
    assert preflight["verification_script_sha256"] == restart.sha256(seed.parent / "verify_production_cache.py")
    assert all(restart.sha256(seed / name) == digest for name, digest in manifest["seed_file_sha256"].items())
    images = read(seed / "seed.traj", index=":")
    assert len(images) == len(observation) == 9 and all(image.calc is None for image in images)
    assert all(restart.same_ordered_geometry(image, original) for image, original in zip(images, observation))
    assert manifest["restart_fmax_eV_A"] == preflight["fmax_eV_A"] == pytest.approx(expected, abs=1e-10)
