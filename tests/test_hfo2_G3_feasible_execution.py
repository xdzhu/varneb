"""Execution gates and injected analytic/cache-only tests; not new DFT evidence."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "benchmarks/hfo2_channels/20261008"
CASE = BASE / "G3_feasible_execution_E171_20261015"
OLDCASE = BASE / "G3_feasible_material_E170_20261015"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


w = load("feasible_execution_tests", CASE / "worker.py")
m = load("frozen_matched_for_execution_tests", OLDCASE / "matched.py")
NATIVE = (
    BASE / "G3_bound_diagnostic_E168_20261015/preserving_native_prefix_and_holes.json"
)


@pytest.fixture(scope="module")
def material():
    c = json.loads((OLDCASE / "config.json").read_bytes())
    (point, records), data = m.validate(
        c,
        NATIVE,
        OLDCASE / "original_point13_cache.zip",
        ROOT / "vcneb/feasible_relaxation.py",
    )
    return c, point, records, data


def test_e166_controls_and_at_most_two_committed_allocations_are_mandatory():
    p = {"28862319": ["COMPLETED", "0:0"], "28862441": ["RUNNING", "0:0"]}
    control = dict(job_id="test", state_exit=["COMPLETED", "0:0"])
    q = [dict(job_id="28862441", name="hfo2-G3-E161-preserving", CPUs=32, nodes=1)]
    assert not w.admission_ready(p, None, q)
    assert not w.admission_ready(p, dict(state_exit=["RUNNING", "0:0"]), q)
    assert w.admission_ready(p, control, q)
    own = dict(job_id="new", name=w.POLICY["job_name"], CPUs=32, nodes=1)
    assert not w.admission_ready(p, control, q + [own])
    assert w.admission_ready(p, control, q + [own], "new")
    assert not w.admission_ready(
        p, control, q + [own, dict(own, job_id="third")], "new"
    )
    assert not w.admission_ready(
        dict(p, **{"28862441": ["TIMEOUT", "0:0"]}), control, []
    )
    assert not w.admission_ready(
        dict(p, **{"28862319": ["FAILED", "1:0"]}), control, []
    )
    with pytest.raises(ValueError, match="both original"):
        w.admission_ready({}, control, [])
    with pytest.raises(ValueError, match="resource"):
        w.admission_ready(p, control, [dict(own, CPUs=64)])


def test_counter_does_not_call_budget_denial_a_new_backend_cache_hit():
    class Delegate:
        n_new_evaluations = 0
        contract = {}
        fixed_coordinate_covectors = np.empty((1, 0))

        def __call__(self, x):
            if x[0] == 2:
                raise RuntimeError("cap")
            self.n_new_evaluations += int(x[0] == 1)
            return 0.0, np.zeros(1)

    counted = w.CountedDelegate(Delegate())
    counted(np.array([1.0]))
    counted(np.array([0.0]))
    with pytest.raises(RuntimeError, match="cap"):
        counted(np.array([2.0]))
    assert counted.requests[0]["new_SCF"] == 1
    assert counted.requests[1]["new_SCF"] == 0
    assert counted.requests[2]["completed"] is False
    assert "new_SCF" not in counted.requests[2]


def test_runtime_roundoff_cannot_trigger_an_extra_starting_SCF():
    calls = []

    class Cache:
        fixed_coordinate_covectors = np.eye(2)

        def __call__(self, x):
            calls.append(x.copy())
            return 0.0, x.copy()

    oracle = w.ExactStartingOracle(Cache(), np.ones(2))
    with pytest.raises(ValueError, match="exactly reuse"):
        oracle(np.array([1.0, 1.0 + 1e-15]))
    assert calls == []
    oracle(np.ones(2))
    assert len(calls) == 1


def test_execution_rejects_unallocated_run_before_any_native_check(
    tmp_path, monkeypatch
):
    monkeypatch.delenv("RUN_DFT", raising=False)
    monkeypatch.setattr(
        w, "native_preflight", lambda *a, **k: pytest.fail("must refuse first")
    )
    with pytest.raises(ValueError, match="allocation"):
        w.run({}, tmp_path, tmp_path / "arm")


def test_execution_policy_does_not_reserve_a_second_training_slot():
    assert w.POLICY["training_reservation_already_charged_in"] == "E170"
    assert w.POLICY["additional_training_slots"] == 0
    assert w.POLICY["maximum_new_attempts"] == 12 and not w.POLICY["climb"]
    assert w.POLICY["E166_eight_controls_must_complete_first"]
    assert w.POLICY["MPI_ranks"] == 32 and w.POLICY["allocation_hours"] == 1


def test_native_probe_binds_original_runtime_and_novel_query_not_portable_cache_guesses(
    material,
):
    c, _, _, _ = material
    old = json.loads((OLDCASE / "native_matched_probe.json").read_bytes())
    assert old["original_runtime_Python_files_rechecked"] == 380
    assert old["old_cache_hits_before_first_novel"] == 7
    assert old["first_original_start_static_exact"]
    assert (
        old["first_novel_internal_norm_A"]
        <= c["policy"]["original_total_internal_radius_A"]
    )
    assert old["first_novel_full_A"][-1] == 0


def test_original_nonconverged_cache_is_not_promoted_by_the_execution_driver(
    material, tmp_path
):
    c, point, records, data = material

    class FailSolver:
        def relax_orthogonal_in_ball(self, plane, q, evaluate, **kw):
            assert (
                kw["gradient_tolerance"] == 0.03 and kw["initial_trust_radius"] == 0.05
            )
            assert kw["orthogonal_radius_bound"] == c["released_radius_A"]
            np.testing.assert_array_equal(
                kw["start"], point["starting_internal_coordinates_A"]
            )
            raise RuntimeError(
                "no full released-gradient convergence within declared ball"
            )

    out = w.compute(
        m,
        c,
        point,
        records,
        data,
        tmp_path / "arm",
        factory=lambda *a: pytest.fail("no DFT in injected failure test"),
        solver=FailSolver(),
    )
    assert out["candidate"] is None and out["status"] == "optimizer_not_converged"
    assert out["new_reserved_attempts"] == out["new_completed_SCFs"] == 0
    assert out["net_speedup_ratio"] is None
    assert json.loads((tmp_path / "arm/arm.json").read_bytes()) == out


def test_unknown_runtime_failure_is_never_relabelled_as_optimizer_nonconvergence(
    material, tmp_path
):
    c, point, records, data = material
    solver = SimpleNamespace(
        relax_orthogonal_in_ball=lambda *a, **k: (_ for _ in ()).throw(
            RuntimeError("native backend error")
        )
    )
    with pytest.raises(RuntimeError, match="native backend error"):
        w.compute(
            m,
            c,
            point,
            records,
            data,
            tmp_path / "arm",
            factory=lambda *a: None,
            solver=solver,
        )
    assert not (tmp_path / "arm/arm.json").exists()


def test_synthetic_zero_gradient_cache_candidate_is_verified_without_a_factory_call(
    material, tmp_path
):
    """Synthetic modification is solely a unit test, not a physical HfO2 result."""
    c, point, records, data = material
    altered = deepcopy(records)
    altered[0]["gradient"] = np.zeros(40).tolist()
    # BLAS roundoff can differ from Linux. Build an explicitly synthetic exact
    # local cache vector; NEVER approximate-match the immutable physical cache.
    from vcneb.mode_surface import ModePlane, _orthogonal_directions

    plane = ModePlane(
        np.zeros(36),
        data["axes"],
        np.ones(36),
        np.eye(2),
        ("coupled internal1", "coupled internal2"),
        "local length",
        "G3-first-conditional-sheet",
    )
    base = plane.frozen_coordinates(m.POLICY["q_A"])
    u = _orthogonal_directions(plane, None)
    x = data["basis"] @ (
        base + u @ (u.T @ (np.asarray(point["starting_internal_coordinates_A"]) - base))
    )
    altered[0]["coordinates"] = x.tolist()
    altered[0]["coordinate_sha256"] = hashlib.sha256(x.tobytes()).hexdigest()
    solver = load("vcneb.feasible_relaxation", ROOT / "vcneb/feasible_relaxation.py")
    out = w.compute(
        m,
        c,
        point,
        altered,
        data,
        tmp_path / "arm",
        factory=lambda *a: pytest.fail("cache-only analytic test"),
        solver=solver,
    )
    assert out["candidate"]["released_gradient_norm_eV_A"] == 0
    assert out["old_cache_hits"] == 2
    assert out["new_completed_SCFs"] == out["new_reserved_attempts"] == 0
    assert out["candidate"]["released_stability_unchecked"]
    assert not out[
        "stable_conditional_minimum_or_TS_or_independent_prediction_certified"
    ]


def test_synthetic_small_residual_is_rejected_even_if_solver_claims_convergence(
    material, tmp_path
):
    c, point, records, data = material
    # Exact first native vector, but the actual released gradient exceeds .03.
    solver = SimpleNamespace(
        relax_orthogonal_in_ball=lambda *a, **k: SimpleNamespace(
            coordinates=np.asarray(point["starting_internal_coordinates_A"]),
            energy=records[0]["enthalpy_eV"],
            orthogonal_gradient_norm=0.0,
        )
    )
    with pytest.raises(ValueError, match="candidate"):
        w.compute(
            m,
            c,
            point,
            records,
            data,
            tmp_path / "arm",
            factory=lambda *a: None,
            solver=solver,
        )


def test_existing_execution_configuration_is_refused_before_preparation(
    tmp_path, monkeypatch
):
    path = tmp_path / "old.json"
    path.write_text("preserve")
    monkeypatch.setattr("sys.argv", ["worker", "prepare", "--config", str(path)])
    monkeypatch.setattr(w, "prepare", lambda *a: pytest.fail("must refuse first"))
    with pytest.raises(FileExistsError):
        w.main()
    assert path.read_text() == "preserve"


def test_execution_manifest_reproduces_and_refuses_allocation_drift(tmp_path):
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    pairs = {
        "matched.py": OLDCASE / "matched.py",
        "matched_config.json": OLDCASE / "config.json",
        "original_point13_cache.zip": OLDCASE / "original_point13_cache.zip",
        "native_prefix.json": NATIVE,
        "feasible_relaxation.py": ROOT / "vcneb/feasible_relaxation.py",
        "worker.py": CASE / "worker.py",
        "run.slurm": CASE / "run.slurm",
    }
    for name, path in pairs.items():
        shutil.copyfile(path, bundle / name)
    config = w.prepare(bundle)
    w.validate(config, bundle)
    bad = deepcopy(config)
    bad["policy"]["maximum_new_attempts"] = 24
    with pytest.raises(ValueError, match="contract differs"):
        w.validate(bad, bundle)
    assert config["policy"]["additional_training_slots"] == 0


def test_control_gate_uses_original_point_id_schema_and_all_eight_rows(
    tmp_path, monkeypatch
):
    source = tmp_path / "source"
    source.mkdir()
    directory = tmp_path / "frozen_R"
    directory.mkdir()
    files = {f"file_{i}.py": "fake-source-digest" for i in range(383)}
    (source / "runtime_files.json").write_text(json.dumps(files))
    cp = (
        source
        / "benchmarks/hfo2_channels/20261008/G3_frozen_R_controls_E166_20261015/config.json"
    )
    cp.parent.mkdir(parents=True)
    points = {c: [dict(id=i) for i in range(1, 5)] for c in ("preserving", "T")}
    cp.write_text(json.dumps(dict(points=points)))
    (tmp_path / "accepted.json").write_text(json.dumps(dict(job_id="unit-job")))
    rows = [
        dict(channel=c, point_id=i, result_sha256="unit")
        for c in points
        for i in range(1, 5)
    ]
    summary = dict(
        status="eight_frozen_R_statics_completed_not_PES_or_independent_forecast",
        job_id="unit-job",
        new_SCFs=8,
        total_attempts=8,
        rows=rows,
    )
    path = directory / "controls.json"
    path.write_text(json.dumps(summary))
    monkeypatch.setattr(w, "CONTROLS", tmp_path)
    monkeypatch.setattr(
        w,
        "sha",
        lambda p: w.CONTROL_RUNTIME_SHA
        if p.name == "runtime_files.json"
        else w.CONTROL_CONFIG_SHA
        if p == cp
        else "fake-source-digest",
    )
    fakefaces = SimpleNamespace(load_inputs=lambda c: dict(channel=c))
    nativeworker = SimpleNamespace(
        faces=fakefaces,
        native_point=lambda data, p, d: dict(point_id=p["id"], result_sha256="unit"),
    )
    monkeypatch.setattr(w, "module", lambda *a: nativeworker)
    r = w.audit_controls(SimpleNamespace(faces=fakefaces))
    assert r["original_native_SCFs_reparsed"] == 8 and r["sources"] == rows
    summary["rows"][2]["result_sha256"] = "changed"
    path.write_text(json.dumps(summary))
    with pytest.raises(ValueError, match="native static differs"):
        w.audit_controls(SimpleNamespace(faces=fakefaces))


def test_real_evaluator_caps_twelve_synthetic_new_calls_and_separates_budget_denial(
    material, tmp_path, monkeypatch
):
    """Actual ASE evaluator/reservations, injected analytic results only."""
    c, point, records, data = material
    from vcneb.mode_surface import ModePlane, _orthogonal_directions

    plane = ModePlane(
        np.zeros(36),
        data["axes"],
        np.ones(36),
        np.eye(2),
        ("coupled internal1", "coupled internal2"),
        "local length",
        "G3-first-conditional-sheet",
    )
    base = plane.frozen_coordinates(m.POLICY["q_A"])
    u = _orthogonal_directions(plane, None)
    start = base + u @ (
        u.T @ (np.asarray(point["starting_internal_coordinates_A"]) - base)
    )
    altered = deepcopy(records)
    full = data["basis"] @ start
    altered[0]["coordinates"] = full.tolist()
    altered[0]["coordinate_sha256"] = hashlib.sha256(full.tobytes()).hexdigest()
    monkeypatch.setattr(
        m.faces.study, "result_validator", lambda *a: dict(elapsed_seconds=0.01)
    )
    calls = []

    def factory(index, atoms, directory):
        calls.append(index)
        return SinglePointCalculator(
            atoms, energy=0.0, forces=np.zeros((12, 3)), stress=np.zeros(6)
        )

    class Solver:
        def relax_orthogonal_in_ball(self, plane, q, evaluate, **kw):
            evaluate(start)
            for i in range(1, 14):
                evaluate(start + u[:, 0] * (i * 0.0001))

    out = w.compute(
        m, c, point, altered, data, tmp_path / "arm", factory=factory, solver=Solver()
    )
    assert calls == list(range(12))
    assert (
        out["candidate"] is None
        and out["new_reserved_attempts"] == out["new_completed_SCFs"] == 12
    )
    assert out["old_cache_hits"] == 1 and out["new_backend_cache_hits"] == 0
    assert (
        len(out["new_delegate_requests"]) == 13
        and not out["new_delegate_requests"][-1]["completed"]
    )
    assert len(list((tmp_path / "arm/new_evaluations").glob("attempt-*.json"))) == 12


def test_wrapper_keeps_original_mpi_and_environment_without_an_implicit_submit():
    text = (CASE / "run.slurm").read_text()
    assert "#SBATCH -p hfacnormal01" in text and "#SBATCH --cpus-per-task=32" in text
    assert 'PYTHONPATH="${old}"' in text and "MKL_NUM_THREADS=1" in text
    assert "sbatch " not in text and "srun " not in text
    assert "/scripts/env.sh" in text and '--workdir "${bundle}/../arm"' in text
    attributes = (CASE / ".gitattributes").read_text()
    assert "*.json -text -eol" in attributes and "*.xml -text -eol" in attributes


def test_actual_native_preflight_and_static_admission_do_not_claim_material_arm_results():
    assert (
        w.sha(CASE / "execution.json")
        == "b80dc0156ca0a7d6bb2dbc26912dd286cc6c6bcb5b2370132beb1fd4e38e5911"
    )
    assert (
        w.sha(CASE / "bundle_files.json")
        == "d4f7a2fbf1950001b725eed0ff1d286e98825c23d62d9a070ac5bf6d4aaa00d2"
    )
    assert (
        w.sha(CASE / "native_execution_preflight.json")
        == "1ffb76f6ea708f68979cdc2df986d213480baf8c1d5a616abfd266272ad60b54"
    )
    proof = json.loads((CASE / "native_execution_preflight.json").read_bytes())
    assert (
        proof["original_runtime_files"] == 380 and proof["original_native_SCFs"] == 12
    )
    assert proof["old_cache_hits_before_first_novel"] == 7
    assert (
        not proof["controls_required"]
        and proof["controls"] is None
        and proof["new_DFT"] == 0
    )
    native = json.loads((CASE / "T_complete_native_audit.json").read_bytes())
    assert (
        native["audit_mode"] == "complete" and native["completed_nodes_audited"] == 25
    )
    assert native["gradient_candidates"] == 25 and native["new_SCFs_reparsed"] == 118
    assert native["sacct"].split("|")[1:3] == ["COMPLETED", "0:0"]
    assert not native[
        "released_stability_branch_uniqueness_contour_TS_or_forecast_certified"
    ]
    receipt = json.loads((CASE / "E166_submission_receipts.json").read_bytes())
    assert receipt["accepted"]["job_id"] == "28865656"
    assert receipt["submission_intent"]["T_complete_native_audit_sha256"] == w.sha(
        CASE / "T_complete_native_audit.json"
    )
    snapshot = json.loads((CASE / "milestone_observation.json").read_bytes())
    assert (
        not snapshot["optimizer_arm_submitted"]
        and not snapshot["optimizer_arm_started"]
    )
    assert snapshot["admission"]["control"] == dict(
        job_id="28865656", state_exit=["PENDING", "0:0"]
    )
    assert "AssocGrpCpuLimit" in snapshot["queue"]
    assert (
        snapshot["preserving_processed"] == 19
        and snapshot["preserving_gradient_candidates"] == 17
    )
