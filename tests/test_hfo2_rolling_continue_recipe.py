"""Exact embedded E090 startup checks, synthetic files only, zero physical DFT."""
import hashlib
import json
import pathlib
from pathlib import Path
import sys

import pytest

from scripts.export_hfo2_clamped_observation import REGISTERED_FACTORIES, ROLLING_E090_SCRIPT_SHA256
from scripts.audit_hfo2_clamped_terminal import fresh_interior_call_count, validate_accounting_source

REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY/"cluster/hf_hfo2_clamped_continue_E090_20261011.slurm"
RUNS = "/public/home/iai806/abacus/agent-runs/"
PARENTS = {
    "28723655": (40, 0, 334, RUNS+"20261010-varneb-clamped-G2-E054-r1/source-fixed",
                 RUNS+"20261011-varneb-clamped-flip-E069-r1/PO_flip_T_pattern_preserving/band"),
    "28709789": (10, 7, 338, RUNS+"20261010-varneb-nested-controls-E055-r1/source",
                 RUNS+"20261010-varneb-clamped-G2-E060-r1/PO_to_M/band"),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def fixture(tmp_path, parent):
    step, initial, count, source, work = PARENTS[parent]
    seed, audit, runtime = tmp_path/"seed/manifest.json", tmp_path/"audit/audit_receipt.json", tmp_path/"runtime"
    files = {}
    for name in ("initial.vasp", "final.vasp", "seed.traj", "substrate.vasp", "factory_parameters.json"):
        path = seed.parent/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("synthetic recipe fixture, not a physical structure", encoding="utf-8")
        files[name] = digest(path)
    observation = audit.parent/"observations"/f"step_{step:04d}"/"observation.json"
    write_json(observation, {"synthetic_recipe_fixture_not_DFT": True})
    manifest = dict(status="audited_G2_geometry_continuation_prepared", source_job_id=parent,
        source_step=step, current_frame_exact_caches=9, new_independent_chains=0,
        new_DFT_calls=0, physical_inputs_changed=False, climb=False, holdout_generated_or_read=False,
        optimizer_state_restored=False, pressure_GPa=0, ordinary_fmax_eV_A=.10,
        n_total_images=9, n_fixed_endpoints=2, n_active_images=7, cell_scale_A=5.12968067458423,
        source_observation_sha256=digest(observation), files_sha256=files)
    write_json(seed, manifest)
    code = {}
    for name in ["vcneb/material_runner.py", *[f"scripts/fixture_{i}.py" for i in range(count-1)]]:
        path = runtime/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# Synthetic recipe-only source hash\n", encoding="utf-8")
        code[name] = digest(path)
    receipt = dict(status="audited_step_cap_not_converged_same_chain_seed", job_id=parent, workdir=work,
        terminal_step=step, ordinary_converged=False, scheduler_row=parent+"|COMPLETED|0:0|32|hfacnormal01|iai806",
        fresh_interior_SCFs=7*step+initial, initial_fresh_interior_SCFs=initial,
        current_frame_exact_caches_checked=9, seed_manifest_sha256=digest(seed),
        all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True,
        physical_inputs_or_running_source_changed=False, holdout_generated_or_read=False,
        FIRE_state_restored=False, new_DFT_calls=0, scheduler_mutations=0, new_independent_chains=0,
        terminal_observation=json.loads(observation.read_text()), runtime_code_sha256=code)
    write_json(audit, receipt)
    return seed, audit, runtime, source


def run_validator(seed, audit, runtime, source, parent, monkeypatch):
    # Only pathlib's native-source lookup is routed to a temporary fixture;
    # execute the unchanged embedded program, not a rewritten validator.
    # Capture the platform's concrete class: patching the generic Path factory
    # changes the identity check in Path.__new__ on both Windows and POSIX.
    real_path = type(runtime)
    body = SCRIPT.read_text(encoding="utf-8").split(" 'import hashlib,json,pathlib,sys\n", 1)[1]
    body = "import hashlib,json,pathlib,sys\n"+body.split("' \"${seed_root}/manifest.json\"", 1)[0]
    with monkeypatch.context() as context:
        context.setattr(sys, "argv", ["recipe", str(seed), digest(seed), str(audit), digest(audit), parent, source])
        context.setattr(pathlib, "Path", lambda p: runtime if str(p)==source else real_path(p))
        exec(compile(body, str(SCRIPT), "exec"), {})


@pytest.mark.parametrize("parent", list(PARENTS))
def test_two_known_recipe_branches_no_structure_or_DFT_claim(tmp_path, parent, monkeypatch, capsys):
    run_validator(*fixture(tmp_path, parent), parent, monkeypatch)
    assert float(capsys.readouterr().out) == 5.12968067458423


@pytest.mark.parametrize("mutation", ["converged", "live", "failed", "historical_step", "cache_count",
    "new_chain", "DFT", "physical", "holdout", "FIRE", "manifest_type", "seed_file", "runtime",
    "observation", "receipt_observation", "SCF_count", "source", "unknown_parent"])
def test_unauthorized_or_changed_recipe_rejects(tmp_path, mutation, monkeypatch):
    parent = "28723655"
    seed, audit, runtime, source = fixture(tmp_path, parent)
    manifest, receipt = [json.loads(p.read_text()) for p in (seed, audit)]
    if mutation == "converged": receipt["ordinary_converged"] = True
    elif mutation in ("live", "failed"): receipt["scheduler_row"] = receipt["scheduler_row"].replace("COMPLETED", mutation.upper())
    elif mutation == "historical_step": manifest["source_step"] = 35
    elif mutation == "cache_count": receipt["current_frame_exact_caches_checked"] = 8
    elif mutation == "new_chain": manifest["new_independent_chains"] = 1
    elif mutation == "DFT": receipt["new_DFT_calls"] = 1
    elif mutation == "physical": manifest["physical_inputs_changed"] = True
    elif mutation == "holdout": receipt["holdout_generated_or_read"] = True
    elif mutation == "FIRE": receipt["FIRE_state_restored"] = True
    elif mutation == "manifest_type": manifest["new_DFT_calls"] = False  # not integer0
    elif mutation == "seed_file": (seed.parent/"seed.traj").write_text("changed")
    elif mutation == "runtime": (runtime/"vcneb/material_runner.py").write_text("changed")
    elif mutation == "observation": (audit.parent/"observations/step_0040/observation.json").write_text("changed")
    elif mutation == "receipt_observation": receipt["terminal_observation"] = {"changed": True}
    elif mutation == "SCF_count": receipt["fresh_interior_SCFs"] = 77
    elif mutation == "source": source += "/different"
    else: parent = "28709790"  # existing held starter is not a completed parent
    write_json(seed, manifest)
    receipt["seed_manifest_sha256"] = digest(seed)
    write_json(audit, receipt)
    with pytest.raises(ValueError):
        run_validator(seed, audit, runtime, source, parent, monkeypatch)


def test_recipe_registration_and_finite_original_driver_contract():
    assert digest(SCRIPT) == ROLLING_E090_SCRIPT_SHA256
    assert REGISTERED_FACTORIES[ROLLING_E090_SCRIPT_SHA256].endswith(":make_clamped_resume_cached_factory")
    text = SCRIPT.read_text(encoding="utf-8")
    for token in ("--steps 40", "--time=12:00:00", "--cpus-per-task=32", "mpirun -np 32",
        "--no-climb", "--fmax .10", "--maxstep .02", "--k .2", "--minimum-distance 1.6",
        "--maximum-deformation .25", "--maximum-cell-step .02", "--candidate-step-retries 4",
        "--mapping identity", "--mic", "--no-align-cells", "--cell-interpolation linear",
        "--require-continuous-periodic-lift", "--clamped-allow-tilt true"):
        assert token in text
    assert "sbatch" not in text and "scontrol" not in text and "srun" not in text
    assert text.index("actual completed unconverged") < text.index("module purge")


@pytest.mark.parametrize("parent", list(PARENTS))
def test_future_audit_recipe_uses_zero_initial_calls_and_channel_named_logs(parent):
    # Synthetic future accounting, not a submitted or accepted native job.
    parameters = {"source_directory": "original-six-byte-source",
                  "seed_cache_records": [{"directory": f"original-{i}"} for i in range(9)]}
    preflight = {"factory": REGISTERED_FACTORIES[ROLLING_E090_SCRIPT_SHA256],
                 "calculator_parameters": parameters}
    assert fresh_interior_call_count(40, preflight, ROLLING_E090_SCRIPT_SHA256) == (0, 280)
    channel = "PO_flip_T_pattern_preserving" if parent == "28723655" else "PO_to_M"
    root = Path("/synthetic-new-namespace")
    work = root/channel/"band"
    fields = ["99900001", "COMPLETED", "iai806", "32", "hfacnormal01",
              "hfo2-G2-continue-E090", str(root/(channel+".slurm.out")),
              str(root/(channel+".slurm.err")),
              f"sbatch --export=RUN_DFT=1,WORKDIR={work} {SCRIPT}"]
    # Quote the Windows fixture paths; native production paths contain no spaces.
    import shlex
    fields[-1] = "sbatch " + shlex.quote(f"--export=RUN_DFT=1,WORKDIR={work}") + " " + shlex.quote(str(SCRIPT))
    assert validate_accounting_source("|".join(fields), "99900001", work, SCRIPT)["NumCPUs"] == "32"
