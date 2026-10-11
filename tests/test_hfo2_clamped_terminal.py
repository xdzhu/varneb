import json
import shlex
from pathlib import Path
import tarfile

import pytest

from scripts.audit_hfo2_clamped_terminal import completed_step, validate_terminal, verify_runtime, validate_job_source, validate_accounting_source
from scripts.audit_hfo2_clamped_terminal import fresh_interior_call_count
from scripts.export_hfo2_clamped_observation import PILOT_SCRIPT_SHA256, M_CONTINUE_SCRIPT_SHA256
from scripts.audit_hfo2_static_replica import sha256


def summary():
    return dict(status="converged", converged=True, termination="force_threshold",
                final_max_generalized_force_eV_per_A=.09, climbing_image_active_final=False)


def seed_preflight(resume=False):
    return dict(factory="examples.hfo2_fixed_input_factory:"+(
        "make_clamped_resume_cached_factory" if resume else "make_clamped_seed_cached_factory"),
        calculator_parameters=dict(source_directory="pinned/source",
            seed_cache_records=[{"pinned": True} if resume or i in (0, 8) else None for i in range(9)]))


@pytest.mark.parametrize("step", [0, 1, 6, 10])
@pytest.mark.parametrize("resume", [False, True])
def test_registered_initial_fresh_calls_are_not_lost(step, resume):
    initial, total = fresh_interior_call_count(step, seed_preflight(resume),
        M_CONTINUE_SCRIPT_SHA256 if resume else PILOT_SCRIPT_SHA256)
    assert initial == (0 if resume else 7)
    assert total == 7*step+initial


@pytest.mark.parametrize("case", ["unknown_recipe", "wrong_factory", "no_endpoint",
    "extra_cache", "partial_resume", "physical_override", "negative_step", "boolean_step"])
def test_fresh_count_refuses_inferred_or_changed_cache_policy(case):
    resume = case == "partial_resume"
    preflight = seed_preflight(resume)
    digest = M_CONTINUE_SCRIPT_SHA256 if resume else PILOT_SCRIPT_SHA256
    step = 6
    if case == "unknown_recipe": digest = "a"*64
    if case == "wrong_factory": preflight["factory"] = "other:calculator"
    if case == "no_endpoint": preflight["calculator_parameters"]["seed_cache_records"][0] = None
    if case == "extra_cache": preflight["calculator_parameters"]["seed_cache_records"][1] = {"pinned": True}
    if case == "partial_resume": preflight["calculator_parameters"]["seed_cache_records"][1] = None
    if case == "physical_override": preflight["calculator_parameters"]["ecutwfc"] = 120
    if case == "negative_step": step = -1
    if case == "boolean_step": step = True
    with pytest.raises(ValueError):
        fresh_interior_call_count(step, preflight, digest)


def test_numeric_convergence_not_scheduler_success():
    validate_terminal("123", "123|COMPLETED|0:0|32|hfacnormal01|iai806", summary(), False)


@pytest.mark.parametrize("field", [None, "JobId", "JobState", "JobName", "UserId",
    "Partition", "NumCPUs", "Command", "StdOut", "StdErr"])
def test_sacct_handle_cannot_be_attached_to_another_workdir(tmp_path, field):
    root = tmp_path/"run"
    work = root/"PO_to_M/band"
    script = root/"production.slurm"
    values = dict(JobId="123", JobState="COMPLETED", JobName="hfo2-G2-M",
        UserId="iai806(16284)", Partition="hfacnormal01", NumCPUs="32",
        Command=str(script), StdOut=str(root/"PO_to_M.slurm.out"), StdErr=str(root/"PO_to_M.slurm.err"))
    if field:
        values[field] = "different"
    body = " ".join(f"{k}={v}" for k,v in values.items())
    if field:
        with pytest.raises(ValueError, match="namespaces differ"):
            validate_job_source(body, "123", work, script)
    else:
        assert validate_job_source(body, "123", work, script)["JobId"] == "123"


@pytest.mark.parametrize("change", [None, "work", "script", "output", "job", "state", "multiple_export"])
def test_persistent_accounting_after_scontrol_expiry(tmp_path, change):
    work = tmp_path/"run/PO_to_M/band"
    script = tmp_path/"run/production.slurm"
    submit = "sbatch --parsable --hold --export="+shlex.quote("ALL,RUN_DFT=1,WORKDIR="+str(work))+" "+shlex.quote(str(script))
    values = ["123", "COMPLETED", "iai806", "32", "hfacnormal01", "hfo2-G2-M",
              str(tmp_path/"run/PO_to_M.slurm.out"), str(tmp_path/"run/PO_to_M.slurm.err"), submit]
    if change == "work": values[-1] = submit.replace(str(work), "/other/band")
    if change == "script": values[-1] = submit.replace(str(script), "/other/script")
    if change == "output": values[6] = "/other/output"
    if change == "job": values[0] = "124"
    if change == "state": values[1] = "RUNNING"
    if change == "multiple_export": values[-1] = submit+" --export=ALL,WORKDIR=/other"
    if change:
        with pytest.raises(ValueError):
            validate_accounting_source("|".join(values), "123", work, script)
    else:
        assert validate_accounting_source("|".join(values), "123", work, script)["JobId"] == "123"


@pytest.mark.parametrize("change", [dict(status="max_steps_reached"), dict(converged=False),
    dict(termination="step_limit"), dict(climbing_image_active_final=True),
    dict(final_max_generalized_force_eV_per_A=.100001),
    dict(final_max_generalized_force_eV_per_A=-.01),
    dict(final_max_generalized_force_eV_per_A=True),
    dict(final_max_generalized_force_eV_per_A=float("nan")),
    dict(final_max_generalized_force_eV_per_A=float("inf"))])
def test_refuses_nonordinary_incomplete_or_nonfinite(change):
    report = {**summary(), **change}
    with pytest.raises(ValueError, match="ordinary"):
        validate_terminal("123", "123|COMPLETED|0:0|32|hfacnormal01|iai806", report, False)


@pytest.mark.parametrize("row", ["123|RUNNING|0:0|32|hfacnormal01|iai806",
    "123|TIMEOUT|0:0|32|hfacnormal01|iai806", "124|COMPLETED|0:0|32|hfacnormal01|iai806",
    "123|COMPLETED|1:0|32|hfacnormal01|iai806", "123|COMPLETED|0:0|1|hfacnormal01|iai806",
    "123|COMPLETED|0:0|32|other|iai806", "123|COMPLETED|0:0|32|hfacnormal01|other"])
def test_registered_actual_terminal_identity(row):
    with pytest.raises(ValueError, match="actual completed"):
        validate_terminal("123", row, summary(), False)


def test_failure_report_not_ignored():
    with pytest.raises(ValueError, match="ordinary"):
        validate_terminal("123", "123|COMPLETED|0:0|32|hfacnormal01|iai806", summary(), True)


@pytest.mark.parametrize("steps", [[], [0, 2], [0, 1, 1], [1, 2]])
def test_completed_history_refuses_partial_duplicate_or_skipped_steps(tmp_path, steps):
    path = tmp_path/"optimizer.log"
    path.write_text("\n".join(f"FIRE: {n} 01:00:00 -10 .09" for n in steps))
    with pytest.raises(ValueError, match="consecutive"):
        completed_step(path)


def test_complete_history(tmp_path):
    path = tmp_path/"optimizer.log"
    path.write_text("Step Time Energy fmax\nFIRE: 0 01:00:00 -10 .11\nFIRE: 1 01:01:00 -11 .09\n")
    assert completed_step(path) == 1


def test_runtime_bytes_are_checked_not_just_manifest(tmp_path):
    source = tmp_path/"source"
    (source/"vcneb").mkdir(parents=True)
    script = source/"vcneb/material_runner.py"
    script.write_bytes(b"# immutable\n")
    archive = tmp_path/"runtime.tar"
    with tarfile.open(archive, "w") as bundle:
        bundle.add(script, arcname="vcneb/material_runner.py")
    code = verify_runtime(source, archive, sha256(archive))
    assert code == {"vcneb/material_runner.py": sha256(script)}
    script.write_bytes(b"# edited\n")
    with pytest.raises(ValueError, match="source changed"):
        verify_runtime(source, archive, sha256(archive))
    with pytest.raises(ValueError, match="archive changed"):
        verify_runtime(source, archive, "a"*64)


def test_empty_runtime_archive_is_not_source_validation(tmp_path):
    archive = tmp_path/"empty.tar"
    with tarfile.open(archive, "w"):
        pass
    with pytest.raises(ValueError, match="complete registered runtime"):
        verify_runtime(tmp_path, archive, sha256(archive))


def test_real_unfinished_material_is_not_promoted_to_terminal():
    case = Path(__file__).resolve().parents[1]/"benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011"
    report = json.loads((case/"validation.json").read_text())
    assert report["parent_ordinary_converged"] is False
    wrong = {**summary(), "final_max_generalized_force_eV_per_A": report["terminal_fmax_eV_A"]}
    with pytest.raises(ValueError, match="ordinary"):
        validate_terminal(report["parent_job"], report["parent_job"]+"|COMPLETED|0:0|32|hfacnormal01|iai806", wrong, False)


def test_actual_M_terminal_portable_native_replay():
    import numpy as np
    from ase.io import read
    from ase.calculators.singlepoint import SinglePointCalculator
    from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry, read_fixed_hfo2_stru
    from scripts.audit_hfo2_static_replica import audited_results
    from scripts.analyze_hfo2_clamped_residual import analyze
    from vcneb import VCNEB, clamped_plane_vcneb_boundary
    case = Path(__file__).resolve().parents[1]/"benchmarks/hfo2_channels/20261008/clamped_M_terminal_E068_20261011"
    receipt = json.loads((case/"audit_receipt.json").read_text())
    summary_path = case/"vcneb_summary.json"
    report = receipt["terminal_observation"]
    directory = case/"observations/step_0013"
    assert sha256(case/"audit_receipt.json") == "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c"
    assert sha256(summary_path) == receipt["terminal_summary_sha256"]
    validate_terminal("28722320", receipt["scheduler_row"], json.loads(summary_path.read_text()), False)
    # HF paths are POSIX even on a Windows replay host; source identity was
    # independently checked by native accounting on HF, not reconstructed here.
    assert "--job-name=hfo2-G2-M-E062" in receipt["scheduler_job_source_raw"]
    assert receipt["fresh_interior_SCFs"] == 91 and len(receipt["calls"]) == 91
    assert len(receipt["runtime_code_sha256"]) == 334
    assert receipt["new_DFT_calls"] == receipt["scheduler_mutations"] == 0
    assert receipt["continuous_TS_or_barrier_error_bound_certified"] is False
    images = read(directory/"evaluated_chain.traj", index=":")
    assert sha256(directory/"evaluated_chain.traj") == report["evaluated_chain_sha256"]
    for i,image in enumerate(images):
        raw = directory/"raw"/f"image_{i:04d}"
        pinned = report["raw_image_evaluations"][i]
        assert set(pinned["input_sha256"]) == {*CONTRACT,"STRU"}
        assert all(pinned["input_sha256"][n] == h for n,h in CONTRACT.items())
        # Original six complete bytes were checked on HF. Only redistributable
        # INPUT/KPT/STRU/native log and E/F/stress are replayed on this host.
        for name in ("INPUT","KPT","STRU"):
            assert sha256(raw/name) == pinned["input_sha256"][name]
        assert sha256(raw/"OUT.ABACUS/running_scf.log") == pinned["raw_log_sha256"]
        assert same_ordered_geometry(image, read_fixed_hfo2_stru(raw/"STRU"))
        parsed = audited_results(raw)
        for key,expected in (("energy",pinned["energy_eV_cell"]),("forces",pinned["forces_eV_A"]),
                             ("stress",pinned["stress_ASE_voigt_eV_A3"])):
            np.testing.assert_allclose(parsed[key],expected,atol=1e-12,rtol=0)
        image.calc = SinglePointCalculator(image,**parsed)
    boundary = clamped_plane_vcneb_boundary(12,np.asarray(report["mechanical_boundary"]["reference_cell_A"]),allow_tilt=True)
    chain = VCNEB(images,cell_scale=report["cell_scale_A"],k=.2,climb=False,pressure=0,**boundary.vcneb_kwargs(images))
    assert np.linalg.norm(chain.get_forces(),axis=1).max() == pytest.approx(.0992140881476476,abs=1e-10,rel=0)
    profile = (chain.enthalpies-chain.enthalpies[0])*1000/4
    np.testing.assert_allclose(profile,report["relative_enthalpy_meV_fu"],atol=1e-9,rtol=0)
    assert profile.max() == pytest.approx(98.45812195362669,abs=1e-9,rel=0)
    assert profile.max()-profile[-1] == pytest.approx(190.99861634595106,abs=1e-9,rel=0)
    actual = analyze(directory)
    historical = json.loads((directory/"residual.json").read_text())
    assert actual["ordinary_residual_pass"] is True
    for key in ("fmax_eV_A","perpendicular_only_fmax_same_geometry_eV_A","spring_only_fmax_same_geometry_eV_A"):
        assert actual[key] == pytest.approx(historical[key],abs=1e-10,rel=0)
