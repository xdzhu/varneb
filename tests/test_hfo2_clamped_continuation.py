"""Actual five-SCF replay plus synthetic transport/refusal checks; no DFT."""
import json
from pathlib import Path
import shutil

import numpy as np
import pytest
from ase.io import read

import examples.hfo2_fixed_input_factory as transport
from scripts.audit_hfo2_clamped_canary import audit, replay_directories
from scripts.audit_hfo2_clamped_continuation import audit as audit_continuation, validate_summary
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_clamped_continuation import prepare
from scripts.relax_clamped_ase_endpoint import load_seed

REPOSITORY = Path(__file__).resolve().parents[1]
CASE = REPOSITORY/"benchmarks/hfo2_channels/20261008"
ROOT = CASE/"G1_review_clamped_canary_20261009/completed_HF/endpoint"
SEED = CASE/"clamped_endpoint_seeds/strain_0000/PO_plus/endpoint_seed.json"
STATIC = ROOT/"calculator/image_0000/scf_000004"
CONTINUATION = CASE/"clamped_PO_continuation_E046_20261009"
T_CASE = CASE/"clamped_endpoint_matrix_20261009/strain_0000/T"


def test_actual_canary_replayed_by_stock_ASE_without_general_ABACUS_IO():
    actual = audit(ROOT, SEED, full_physical_bytes=False)
    original = json.loads((ROOT/"canary_audit.json").read_text())
    assert actual["new_SCF_calls"] == 5 and not actual["endpoint_converged"]
    assert actual["raw_full_physical_bytes_checked_here"] is False
    assert original["raw_full_physical_bytes_checked_here"] is True
    assert actual["summary_sha256"] == original["summary_sha256"]
    assert actual["SCF_seconds"] == original["SCF_seconds"]
    for a, b in zip(actual["rows"], original["rows"]):
        assert a["raw_log_sha256"] == b["raw_log_sha256"] and a["input_sha256"] == b["input_sha256"]
        for key in ("energy_eV_cell", "max_atomic_force_eV_A", "open_traction_norm_kbar"):
            assert a[key] == pytest.approx(b[key], abs=1e-10, rel=0)
        assert [s["symbol"] for s in a["structure_audit"]["symmetry_sweep"]] == ["Pca2_1"]*3


def test_exact_writer_geometry_agrees_with_actual_terminal_POSCAR():
    atoms = transport.read_fixed_hfo2_stru(STATIC/"STRU")
    assert transport.same_ordered_geometry(atoms, read(ROOT/"CONTCAR", format="vasp"))


@pytest.mark.parametrize("before,after", [("Direct", "Cartesian"), ("Hf 178.49", "Hf 178.50"),
    ("1 1 1", "1 0 1"), ("5.2390601337432896", "nan"), ("5.2390601337432896", "-5.2390601337432896")])
def test_narrow_reader_refuses_unsupported_or_invalid_STRU(tmp_path, before, after):
    path = tmp_path/"STRU"
    path.write_text((STATIC/"STRU").read_text().replace(before, after, 1))
    with pytest.raises(ValueError):
        transport.read_fixed_hfo2_stru(path)


def test_actual_partial_seed_preparation_keeps_same_boundary_and_no_DFT(tmp_path):
    out = tmp_path/"seed"
    p = prepare(REPOSITORY, ROOT, out, full_physical_bytes=False)
    assert p["new_DFT_calls"] == 0 and p["planned_max_new_SCF_calls"] == 20
    assert not p["parent_review_full_physical_bytes_checked_here"]
    atoms, boundary, manifest = load_seed(out/"endpoint_seed.json")
    assert transport.same_ordered_geometry(atoms, transport.read_fixed_hfo2_stru(STATIC/"STRU"))
    boundary.validate_images([atoms])
    assert manifest["continuation"]["parent_BFGS_steps"] == 4
    assert manifest["continuation"]["restore_BFGS_Hessian"] is False
    parameters = json.loads((out/"factory_parameters.json").read_text())
    assert parameters["seed_input_sha256"]["STRU"] == sha256(STATIC/"STRU")
    with pytest.raises(FileExistsError):
        prepare(REPOSITORY, ROOT, out, full_physical_bytes=False)


def cache_fixture(tmp_path, monkeypatch):
    # Reduced two-file contract tests cache logic only. Real six-file bytes are
    # checked on HF, not inferred from this deliberately limited local fixture.
    source = tmp_path/"seed"
    source.mkdir()
    for name in ("INPUT", "KPT", "STRU"):
        shutil.copyfile(STATIC/name, source/name)
    (source/"OUT.ABACUS").mkdir()
    shutil.copyfile(STATIC/"OUT.ABACUS/running_scf.log", source/"OUT.ABACUS/running_scf.log")
    monkeypatch.setattr(transport, "CONTRACT", {n: sha256(source/n) for n in ("INPUT", "KPT")})
    return {"source_directory": str(source), "seed_static_directory": str(source),
            "seed_input_sha256": {n: sha256(source/n) for n in (*transport.CONTRACT, "STRU")},
            "seed_raw_log_sha256": sha256(source/"OUT.ABACUS/running_scf.log")}


def test_one_endpoint_cache_reused_until_geometry_changes(tmp_path, monkeypatch):
    p = cache_fixture(tmp_path, monkeypatch)
    factory = transport.make_endpoint_cached_factory(parameters=p, command="mpirun -np 32 abacus")
    atoms = transport.read_fixed_hfo2_stru(STATIC/"STRU")
    directory = tmp_path/"calculator"
    directory.mkdir()
    (directory/"structure.start.vasp").write_text("attachment snapshot fixture")
    def refused_launch(*args, **kwargs):
        raise RuntimeError("synthetic fixture reached next fresh SCF")
    monkeypatch.setattr(transport.subprocess, "run", refused_launch)
    monkeypatch.setattr(transport, "geometry_roundtrip", lambda path, a: None)
    atoms.calc = factory(0, atoms, directory)
    assert atoms.get_potential_energy() == pytest.approx(-9782.972556311173)
    assert atoms.get_forces().shape == (12, 3) and atoms.get_stress().shape == (6,)
    assert atoms.calc.next_call == 0
    record = json.loads((directory/"seed_cache_audit.json").read_text())
    assert record["new_DFT_calls"] == 0
    atoms.positions[0, 0] += .001
    with pytest.raises(RuntimeError, match="next fresh SCF"):
        atoms.get_potential_energy()
    assert atoms.calc.next_call == 1 and (directory/"scf_000000").exists()


@pytest.mark.parametrize("mutation", ["log", "geometry", "index", "existing"])
def test_stale_or_ambiguous_endpoint_cache_refused(tmp_path, monkeypatch, mutation):
    p = cache_fixture(tmp_path, monkeypatch)
    factory = transport.make_endpoint_cached_factory(parameters=p, command="mpirun -np 32 abacus")
    atoms = transport.read_fixed_hfo2_stru(STATIC/"STRU")
    index, target = 0, tmp_path/"target"
    if mutation == "log":
        (Path(p["seed_static_directory"])/"OUT.ABACUS/running_scf.log").write_text("changed raw log")
    elif mutation == "geometry":
        atoms.positions[0, 0] += .001
    elif mutation == "index":
        index = 1
    else:
        target.mkdir()
        (target/"old_result.json").write_text("preserve this old namespace")
    with pytest.raises(ValueError):
        factory(index, atoms, target)


def test_fixed_continuation_recipe_and_cache_call_budget(tmp_path):
    out = tmp_path/"seed"
    prepare(REPOSITORY, ROOT, out, full_physical_bytes=False)
    m = json.loads((out/"endpoint_seed.json").read_text())
    s = json.loads((ROOT/"endpoint_relax_summary.json").read_text())
    s.update(steps_requested=20, optimizer_steps=20)
    validate_summary(s, m, 20)
    with pytest.raises(ValueError):
        validate_summary(s, m, 21)  # Recomputed initial seed is not acceptable.
    s["open_stress_target_kbar"] = .1
    with pytest.raises(ValueError):
        validate_summary(s, m, 20)
    wrapper = (REPOSITORY/"cluster/hf_hfo2_clamped_PO_continue_E046_20261009.slurm").read_text()
    assert '--steps 20' in wrapper and 'make_endpoint_cached_factory' in wrapper
    assert '#SBATCH --time=02:00:00' in wrapper and 'sbatch' not in wrapper and 'srun' not in wrapper


def test_actual_two_SCF_continuation_replays_original_HF_audit():
    endpoint = CONTINUATION/"completed_HF/endpoint"
    actual = audit_continuation(endpoint, CONTINUATION/"seed_HF/endpoint_seed.json", ROOT,
                                full_physical_bytes=False)
    original = json.loads((endpoint/"continuation_audit.json").read_text())
    assert actual["endpoint_converged"] and actual["new_SCF_calls"] == 2
    assert actual["fresh_BFGS_steps"] == 2 and actual["reused_seed_SCFs"] == 1
    assert actual["raw_full_physical_bytes_checked_here"] is False
    assert original["raw_full_physical_bytes_checked_here"] is True
    for key in ("summary_sha256", "seed_manifest_sha256", "SCF_seconds", "SCF_core_hours"):
        assert actual[key] == original[key]
    for a, b in zip(actual["rows"], original["rows"]):
        assert a["input_sha256"] == b["input_sha256"] and a["raw_log_sha256"] == b["raw_log_sha256"]
        for key in ("energy_eV_cell", "max_atomic_force_eV_A", "open_traction_norm_kbar"):
            assert a[key] == pytest.approx(b[key], abs=1e-10, rel=0)
    assert actual["terminal_phase_symbols"] == ["Pca2_1"]*3
    assert actual["rows"][-1]["max_atomic_force_eV_A"] < .03
    assert actual["rows"][-1]["open_traction_norm_kbar"] < 2.


def test_actual_registered_T_endpoint_one_raw_SCF_not_optimizer_trajectory():
    endpoint = T_CASE/"completed_HF/endpoint"
    seed = CASE/"clamped_endpoint_seeds/strain_0000/T/endpoint_seed.json"
    _, boundary, _ = load_seed(seed)
    original = json.loads((T_CASE/"completed_HF/T_endpoint_audit.json").read_text())
    summary = json.loads((endpoint/"endpoint_relax_summary.json").read_text())
    calls = sorted((endpoint/"calculator/image_0000").glob("scf_*"))
    actual = replay_directories(calls, boundary, full_physical_bytes=False)
    assert len(actual) == original["new_SCF_calls"] == 1
    assert summary["optimizer_steps"] == 0 and summary["steps_requested"] == 4
    assert summary["converged"] and original["endpoint_converged"]
    assert sha256(seed) == summary["seed_manifest_sha256"] == original["seed_manifest_sha256"]
    assert sha256(endpoint/"endpoint_relax_summary.json") == original["summary_sha256"]
    assert actual[0]["input_sha256"] == original["rows"][0]["input_sha256"]
    assert actual[0]["raw_log_sha256"] == original["rows"][0]["raw_log_sha256"]
    assert [s["symbol"] for s in actual[0]["structure_audit"]["symmetry_sweep"]] == ["P4_2/nmc"]*3
    for key in ("energy_eV_cell", "max_atomic_force_eV_A", "open_traction_norm_kbar"):
        assert actual[0][key] == pytest.approx(original["rows"][0][key], abs=1e-10, rel=0)


def test_actual_common_substrate_well_gap_is_not_a_barrier():
    po = json.loads((CONTINUATION/"completed_HF/endpoint/endpoint_relax_summary.json").read_text())
    t = json.loads((T_CASE/"completed_HF/endpoint/endpoint_relax_summary.json").read_text())
    for s in (po, t):
        assert s["strain"] == s["external_pressure_gpa"] == 0.
        assert s["physical_contract_sha256"] == transport.CONTRACT
        assert s["allow_tilt"] and s["open_stress_target_kbar"] == 2.
        assert s["endpoint"]["composition"] == {"Hf": 4, "O": 8}
    np.testing.assert_allclose(np.asarray(po["cell_A"])[:2], np.asarray(t["cell_A"])[:2], atol=1e-12, rtol=0)
    gap = (t["potential_energy_eV"] - po["potential_energy_eV"])*1000/4
    assert gap == pytest.approx(12.35997511412279, abs=1e-9, rel=0)


@pytest.mark.parametrize("root,count", [(ROOT.parent, 49), (CONTINUATION/"completed_HF", 26),
                                       (T_CASE/"completed_HF", 17)])
def test_actual_observable_exports_are_complete_without_licensed_basis(root, count):
    files = [p for p in root.rglob("*") if p.is_file()]
    assert len(files) == count
    assert not any(p.suffix in (".upf", ".orb") or p.name in ("abacus", "SPIN1_CHG.cube") for p in files)
    if root == ROOT.parent:
        inventory = json.loads((ROOT.parent.parent/"observable_inventory.json").read_text())
        assert len(inventory["files"]) == len(files)
        for row in inventory["files"]:
            p = root/row["path"]
            assert p.stat().st_size == row["bytes"] and sha256(p) == row["sha256"]
