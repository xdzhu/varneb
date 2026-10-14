"""Execute the frozen E170 same-input arm, after the eight E166 controls.

One fresh, separately capped training arm; no dispatcher or default change.
Original failed trials are exact cache entries, never accepted end points.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
import numpy as np

HERE = Path(__file__).resolve().parent
RUNS = Path("/public/home/iai806/abacus/agent-runs")
OLD = RUNS / "20261015-varneb-G3-faces-E161-r1"
CONTROLS = RUNS / "20261015-varneb-G3-frozen-R-E166-r1"
RUNTIME_SHA = "d07bb011c9e313afe23ca4e78207356ff98adfc43ca231bdd047e6d5d0f1be20"
MATCHED_SHA = "f937a2e31ce30ec9b15852fa6a8e06cb73f4abbd4aa9b7f1f26bd29d561b4fd3"
AUDITOR_SHA = "ab041335a93b8eaa98fa5cf8a812db5126931105cf972f586e4d15f8d06742da"
CONTROL_CONFIG_SHA = "5ffb1fe0d707aec5d0a6546a44e0ad0002cd0f0280f1ba17ca2d955ac2d4c930"
CONTROL_RUNTIME_SHA = "d5a02ba70d9db2c1cecdf0eecc82f7b367bd4fe185fdaa57af4819a19adcc5e8"
POLICY = dict(
    job_name="hfo2-G3-E171-feasible",
    partition="hfacnormal01",
    nodes=1,
    MPI_ranks=32,
    allocation_hours=1,
    maximum_study_allocations=2,
    maximum_new_attempts=12,
    training_reservation_already_charged_in="E170",
    additional_training_slots=0,
    E166_eight_controls_must_complete_first=True,
    baseline_recomputed=False,
    climb=False,
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def inputs(bundle):
    return (
        bundle / "native_prefix.json",
        bundle / "original_point13_cache.zip",
        bundle / "feasible_relaxation.py",
    )


def prepare(bundle):
    m = module("matched_E170_frozen_execution", bundle / "matched.py")
    if sha(bundle / "matched_config.json") != MATCHED_SHA:
        raise ValueError("exact E170 configuration bytes required")
    c = json.loads((bundle / "matched_config.json").read_bytes())
    m.validate(c, *inputs(bundle))
    return dict(
        schema_version=1,
        purpose="bounded_matched_material_optimizer_training_arm",
        policy=POLICY,
        matched_config_sha256=MATCHED_SHA,
        original_runtime_manifest_sha256=RUNTIME_SHA,
        bundle_text_sha256={
            n: m.faces.study.text_sha(bundle / n)
            for n in ("matched.py", "feasible_relaxation.py", "worker.py", "run.slurm")
        },
        auditor_sha256=AUDITOR_SHA,
        original_cache_zip_sha256=m.CACHE_SHA,
        native_prefix_sha256=m.NATIVE_SHA,
        old_existing_SCFs=12,
        old_curvature_SCFs_both_centers=296,
        preparation_and_old_cache_costs_additional=True,
        original_control_and_domain_and_start_unchanged=True,
        physical_reserved_geometry_generated_or_read=False,
    )


def validate(config, bundle):
    if json.dumps(config, sort_keys=True, allow_nan=False) != json.dumps(
        prepare(bundle), sort_keys=True, allow_nan=False
    ):
        raise ValueError("execution source, budget or allocation contract differs")
    m = module("matched_E170_frozen_execution", bundle / "matched.py")
    c = json.loads((bundle / "matched_config.json").read_bytes())
    (point, records), data = m.validate(c, *inputs(bundle))
    return m, c, point, records, data


def admission_ready(parents, control, queue, own_job=None):
    if set(parents) != {"28862319", "28862441"}:
        raise ValueError("both original parent handles required")
    if not control or control["state_exit"] != ["COMPLETED", "0:0"]:
        return False
    if any(
        list(v) not in (["COMPLETED", "0:0"], ["RUNNING", "0:0"], ["PENDING", "0:0"])
        for v in parents.values()
    ) or not any(list(v) == ["COMPLETED", "0:0"] for v in parents.values()):
        return False
    if any(r["CPUs"] != 32 or r["nodes"] != 1 for r in queue):
        raise ValueError("committed study resource identity differs")
    if own_job is None:
        return len(queue) < 2
    return (
        len(queue) <= 2
        and sum(
            r["job_id"] == own_job and r["name"] == POLICY["job_name"] for r in queue
        )
        == 1
    )


def scheduler_snapshot(own_job=None):
    receipt = CONTROLS / "accepted.json"
    if not receipt.is_file():
        return dict(ready=False, reason="E166_not_submitted", scheduler_writes=0)
    control_id = json.loads(receipt.read_bytes())["job_id"]
    if not isinstance(control_id, str) or not control_id.isdecimal():
        raise ValueError("exact E166 accepted handle required")
    ids = ["28862319", "28862441", control_id]
    body = subprocess.run(
        [
            "sacct",
            "-n",
            "-P",
            "-X",
            "-j",
            ",".join(ids),
            "--format=JobIDRaw,State,ExitCode",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    states = {}
    for line in body.splitlines():
        f = line.strip().split("|")
        if f[0] in ids:
            if f[0] in states or len(f) < 3:
                raise ValueError("ambiguous accounting handle")
            states[f[0]] = f[1:3]
    if set(states) != set(ids):
        raise ValueError(
            "all original/control accounting handles required; missing is not failure"
        )
    queue_text = subprocess.run(
        [
            "squeue",
            "-h",
            "-u",
            "iai806",
            "-p",
            "hfacnormal01",
            "-o",
            "%i|%.100j|%T|%C|%D",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    queue = []
    for line in queue_text.splitlines():
        f = [s.strip() for s in line.split("|")]
        if len(f) != 5:
            raise ValueError("unparsed queue observation")
        if f[1].startswith("hfo2-G3-"):
            queue.append(
                dict(
                    job_id=f[0], name=f[1], state=f[2], CPUs=int(f[3]), nodes=int(f[4])
                )
            )
    parents = {j: states[j] for j in ids[:2]}
    control = dict(job_id=control_id, state_exit=states[control_id])
    return dict(
        ready=admission_ready(parents, control, queue, own_job),
        parent_states=parents,
        control=control,
        committed_G3_queue=queue,
        own_job_id=own_job,
        scheduler_writes=0,
    )


def verify_runtime(m):
    if m.faces.study.ROOT.resolve() != (OLD / "source").resolve():
        raise ValueError("original immutable E161 Python runtime required")
    manifest = OLD / "source/runtime_files.json"
    if sha(manifest) != RUNTIME_SHA:
        raise ValueError("original runtime manifest differs")
    files = json.loads(manifest.read_bytes())
    if len(files) != 380 or any(sha(OLD / "source" / n) != h for n, h in files.items()):
        raise ValueError("original runtime bytes differ")


def audit_controls(m):
    """Reparse all eight actual E166 statics, not a completion-file assertion."""
    source = CONTROLS / "source"
    if sha(source / "runtime_files.json") != CONTROL_RUNTIME_SHA:
        raise ValueError("E166 runtime manifest differs")
    files = json.loads((source / "runtime_files.json").read_bytes())
    if len(files) != 383 or any(sha(source / n) != h for n, h in files.items()):
        raise ValueError("E166 immutable runtime differs")
    cp = (
        source
        / "benchmarks/hfo2_channels/20261008/G3_frozen_R_controls_E166_20261015/config.json"
    )
    if sha(cp) != CONTROL_CONFIG_SHA:
        raise ValueError("original eight-control configuration differs")
    config = json.loads(cp.read_bytes())
    directory = CONTROLS / "frozen_R"
    summary = json.loads((directory / "controls.json").read_bytes())
    accepted = json.loads((CONTROLS / "accepted.json").read_bytes())
    if (
        summary["status"]
        != "eight_frozen_R_statics_completed_not_PES_or_independent_forecast"
        or summary["job_id"] != accepted["job_id"]
        or summary["new_SCFs"] != 8
        or summary["total_attempts"] != 8
        or len(summary["rows"]) != 8
    ):
        raise ValueError("complete original E166 result and exact handle required")
    before = list(sys.path)
    try:
        native_worker = module(
            "immutable_E166_static_auditor",
            source
            / "benchmarks/hfo2_channels/20261008/G3_frozen_R_controls_E166_20261015/worker.py",
        )
    finally:
        sys.path[:] = before
    if native_worker.faces is not m.faces:
        raise ValueError("original physical worker identity differs")
    sources = []
    for channel in ("preserving", "T"):
        data = m.faces.load_inputs(channel)
        for p in config["points"][channel]:
            d = directory / channel / f"point_{p['id']:03d}"
            audited = dict(channel=channel, **native_worker.native_point(data, p, d))
            expected = [
                r
                for r in summary["rows"]
                if r["channel"] == channel and r["point_id"] == p["id"]
            ]
            if expected != [audited]:
                raise ValueError("actual native static differs from E166 summary")
            sources.append(audited)
    return dict(
        controls_summary_sha256=sha(directory / "controls.json"),
        job_id=accepted["job_id"],
        original_native_SCFs_reparsed=8,
        sources=sources,
        new_DFT=0,
    )


class CountedDelegate:
    """Separate new cache hits, completed calls and denied-budget requests."""

    def __init__(self, evaluator):
        self.evaluator = evaluator
        self.contract = evaluator.contract
        self.fixed_coordinate_covectors = evaluator.fixed_coordinate_covectors
        self.requests = []

    def __call__(self, x):
        before = self.evaluator.n_new_evaluations
        row = dict(
            coordinate_sha256=hashlib.sha256(np.asarray(x).tobytes()).hexdigest()
        )
        self.requests.append(row)
        try:
            value = self.evaluator(x)
        except Exception as exc:
            row.update(
                completed=False, exception_type=type(exc).__name__, message=str(exc)
            )
            raise
        row.update(completed=True, new_SCF=self.evaluator.n_new_evaluations - before)
        return value


class ExactStartingOracle:
    """Fail before DFT if runtime arithmetic no longer reaches the old start."""

    def __init__(self, cache, original_first):
        self.cache = cache
        self.original_first = np.asarray(original_first)
        self.fixed_coordinate_covectors = cache.fixed_coordinate_covectors
        self.called = False

    def __call__(self, x):
        if not self.called and not np.array_equal(x, self.original_first):
            raise ValueError(
                "first evaluation must exactly reuse the original native starting SCF"
            )
        self.called = True
        return self.cache(x)


def native_preflight(config, bundle, *, require_controls):
    m, c, point, records, data = validate(config, bundle)
    verify_runtime(m)
    m.faces.native_preflight(m.faces.prepare("preserving"), "preserving")
    auditor = bundle / "audit_faces.py"
    if sha(auditor) != AUDITOR_SHA:
        raise ValueError("existing E163 auditor bytes required")
    audit = module("original_face_auditor_for_execution", auditor)
    native, _ = audit.native_records(data, Path(m.POINT_PATH).parent / "node_013")
    if native != records:
        raise ValueError("original cached/native records differ")
    baseline = audit.replay_node(
        data, point, records, np.array(point["starting_internal_coordinates_A"])
    )
    probe = m.probe(c, *inputs(bundle))
    controls = audit_controls(m) if require_controls else None
    verify_runtime(m)
    return (m, c, point, records, data), dict(
        original_runtime_files=380,
        original_native_SCFs=12,
        old_baseline_replay=baseline,
        old_cache_hits_before_first_novel=probe["old_cache_hits_before_first_novel"],
        controls=controls,
        controls_required=require_controls,
        new_DFT=0,
    )


def compute(m, c, point, records, data, workdir, *, factory, solver):
    """Injectable calculator/solver for tests; production preflight is mandatory in run."""
    workdir.mkdir()
    study = m.faces.study
    delegate = m.faces.canary.CalculatorModeEvaluator(
        data["chart"],
        factory,
        workdir=workdir / "new_evaluations",
        calculator_id="original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2",
        minimum_distance_A=1.6,
        result_validator=study.result_validator,
        validation_id="E150-native-six-input-DSIZE32-EFS-one-SCF-v1",
        max_evaluation_attempts=12,
    )
    counted = CountedDelegate(delegate)
    oracle = m.SharedExactCache(records, counted)
    evaluate = m.faces.canary.InternalEvaluator(
        ExactStartingOracle(oracle, records[0]["coordinates"]),
        data["basis"],
        data["center_energy"],
        data["center_gradient"],
    )
    plane = m.faces.canary.ModePlane(
        np.zeros(36),
        data["axes"],
        np.ones(36),
        np.eye(2),
        ("coupled internal1", "coupled internal2"),
        "local length",
        "G3-first-conditional-sheet",
    )
    screens = dict(checks=0, rejected=0)

    def guard(x):
        screens["checks"] += 1
        ok = m.faces.canary.geometry_guard(data["chart"], data["basis"], x)
        screens["rejected"] += int(not ok)
        return ok

    result, reason = None, None
    try:
        result = solver.relax_orthogonal_in_ball(
            plane,
            m.POLICY["q_A"],
            evaluate,
            start=np.array(point["starting_internal_coordinates_A"]),
            orthogonal_radius_bound=c["released_radius_A"],
            trial_validator=guard,
            initial_hessian=data["hessian"],
            initial_trust_radius=0.05,
            gradient_tolerance=0.03,
            max_iterations=200,
        )
    except RuntimeError as exc:
        reason = str(exc)
        if (
            list((workdir / "new_evaluations").glob("eval-*/failure.json"))
            or reason
            not in (
                "mode evaluation attempt budget exhausted; calculator was not called",
                "no full released-gradient convergence within declared ball",
            )
            or (
                reason.startswith("mode evaluation attempt budget exhausted")
                and delegate.n_evaluation_attempts != 12
            )
        ):
            raise
    fresh = [
        json.loads(p.read_bytes())
        for p in sorted((workdir / "new_evaluations").glob("eval-*/result.json"))
    ]
    if (
        len(fresh) != delegate.n_new_evaluations
        or not 0 <= delegate.n_evaluation_attempts <= 12
    ):
        raise ValueError("new attempt/result accounting differs")
    candidate = None
    if result is not None:
        full = data["basis"] @ result.coordinates
        matches = [
            r for r in [*records, *fresh] if np.array_equal(r["coordinates"], full)
        ]
        if len(matches) != 1 or not guard(result.coordinates):
            raise ValueError(
                "candidate must match one exact complete physical result and original guard"
            )
        row = matches[0]
        g = data["basis"].T @ np.asarray(row["gradient"])
        residual = float(np.linalg.norm(g - data["axes"] @ (data["axes"].T @ g)))
        if (
            row["enthalpy_eV"] != result.energy
            or residual > 0.03
            or not np.isclose(
                residual, result.orthogonal_gradient_norm, atol=2e-12, rtol=0
            )
            or not np.allclose(
                data["axes"].T @ result.coordinates, m.POLICY["q_A"], atol=2e-12, rtol=0
            )
            or full[-1] != 0
        ):
            raise ValueError("candidate E/Q/full released-gradient/control differs")
        candidate = dict(
            coordinates_internal_A=result.coordinates.tolist(),
            coordinates_full_A=full.tolist(),
            energy_eV_cell=result.energy,
            released_gradient_norm_eV_A=residual,
            native_validation=row["validation"],
            released_stability_unchecked=True,
        )
    summary = dict(
        status="gradient_candidate_release_stability_unchecked"
        if candidate
        else "optimizer_not_converged",
        candidate=candidate,
        reason=reason,
        policy=POLICY,
        q_A=m.POLICY["q_A"],
        original_start=point["starting_internal_coordinates_A"],
        new_completed_SCFs=delegate.n_new_evaluations,
        new_reserved_attempts=delegate.n_evaluation_attempts,
        old_cache_hits=oracle.old_cache_hits,
        physical_requests=oracle.requests,
        internal_callback_requests=evaluate.callback_requests,
        center_cache_hits=evaluate.center_cache_hits,
        new_delegate_requests=counted.requests,
        new_backend_cache_hits=sum(
            r.get("completed") is True and r["new_SCF"] == 0 for r in counted.requests
        ),
        geometry_screens=screens,
        sum_new_completed_call_seconds=sum(
            r["validation"]["elapsed_seconds"] for r in fresh
        ),
        old_baseline_existing_SCFs=12,
        old_baseline_recomputed_SCFs=0,
        prior_curvature_both_centers_SCFs=296,
        net_speedup_ratio=None,
        old_capped_failure_is_censored_not_speedup_denominator=True,
        preparation_costs_additional=True,
        physical_reserved_geometry_generated_or_read=False,
        stable_conditional_minimum_or_TS_or_independent_prediction_certified=False,
        job_id=os.environ.get("SLURM_JOB_ID"),
    )
    study.exclusive_json(workdir / "arm.json", summary)
    return summary


def run(config, bundle, workdir):
    expected = dict(
        RUN_DFT="1",
        SLURM_JOB_PARTITION="hfacnormal01",
        SLURM_JOB_NUM_NODES="1",
        SLURM_NTASKS="1",
        SLURM_CPUS_PER_TASK="32",
        SLURM_JOB_NAME=POLICY["job_name"],
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
    )
    if any(os.environ.get(k) != v for k, v in expected.items()) or not os.environ.get(
        "SLURM_JOB_ID"
    ):
        raise ValueError(
            "authorized one-node32MPI single-thread HF allocation required"
        )
    if (
        workdir.resolve() != bundle.resolve().parent / "arm"
        or workdir.exists()
        or not workdir.resolve().is_relative_to(RUNS)
        or not bundle.resolve().is_relative_to(RUNS)
    ):
        raise ValueError("fresh dedicated arm sibling of immutable bundle required")
    admission = scheduler_snapshot(os.environ["SLURM_JOB_ID"])
    if not admission["ready"]:
        raise RuntimeError("E166/control/resource gate not ready; no DFT")
    (m, c, point, records, data), proof = native_preflight(
        config, bundle, require_controls=True
    )
    solver = module("vcneb.feasible_relaxation", bundle / "feasible_relaxation.py")
    m.faces.study.exclusive_json(
        bundle.parent / "run_native_preflight.json", dict(**proof, admission=admission)
    )
    result = compute(
        m,
        c,
        point,
        records,
        data,
        workdir,
        factory=m.faces.study.make_factory(
            parameters={"source_directory": str(m.faces.study.ASSETS)},
            command=m.faces.study.COMMAND,
        ),
        solver=solver,
    )
    verify_runtime(m)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("prepare", "preflight", "run"))
    p.add_argument("--bundle", type=Path, default=HERE)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--config-sha256")
    p.add_argument("--output", type=Path)
    p.add_argument("--workdir", type=Path)
    args = p.parse_args()
    if args.mode == "prepare":
        if args.config.exists():
            raise FileExistsError("fresh execution configuration required")
        result = prepare(args.bundle)
    else:
        if not args.config_sha256 or sha(args.config) != args.config_sha256:
            raise ValueError("exact frozen execution configuration required")
        config = json.loads(args.config.read_bytes())
        if args.mode == "preflight":
            if args.output is None or args.output.exists():
                raise FileExistsError("fresh preflight report required")
            _, result = native_preflight(config, args.bundle, require_controls=False)
        else:
            if args.workdir is None:
                raise ValueError("fresh arm output required")
            result = run(config, args.bundle, args.workdir)
            print(
                json.dumps(
                    {
                        k: result[k]
                        for k in (
                            "status",
                            "new_completed_SCFs",
                            "new_reserved_attempts",
                        )
                    }
                )
            )
            return
    target = args.config if args.mode == "prepare" else args.output
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")
    print(json.dumps(dict(status=args.mode + "_passed_no_DFT", new_DFT=0)))


if __name__ == "__main__":
    main()
