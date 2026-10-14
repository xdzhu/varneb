"""Prepare/probe a matched material arm; no calculator or submission entry yet."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import zipfile

sys.dont_write_bytecode = True
import numpy as np
from scripts import hfo2_G3_conditional_faces as faces

CASE = Path(__file__).resolve().parent
NATIVE_SHA = "7c9581eb0c102f6110415ca590b7f045b4875e3d3498a7af5019d0cb8a14f57b"
CACHE_SHA = "02cf4f9b4bd169a775762411b798554c345be1496309b64b7ef6218c5ab99610"
POINT_PATH = "/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-faces-E161-r1/preserving/point_013.json"
POINT_SHA = "43150958a4e8f8a2e8d16e56f97f8ac6046a467a4359651336f8f87778fffcec"
POLICY = dict(
    reservation_id="E170-one-preserving-point13-feasible-optimizer-arm",
    channel="preserving",
    point_id=13,
    q_A=[-0.2, -0.1],
    new_training_slots=1,
    preserving_training_committed=30,
    preserving_training_remaining=51,
    T_training_committed=29,
    T_training_remaining=52,
    independent_slots_untouched_per_face=8,
    total_training_allowance_per_face=81,
    maximum_new_attempts=12,
    maximum_iterations=200,
    original_total_internal_radius_A=0.5,
    gradient_target_eV_A=0.03,
    initial_trust_radius_A=0.05,
    minimum_distance_A=1.6,
    pressure_GPa=0.0,
    strain=0.0,
    ordinary_NEB_fmax_eV_A=0.10,
    allocation_nodes=1,
    MPI_ranks=32,
    allocation_hours=1,
    maximum_committed_study_allocations=2,
    original_E166_controls_have_priority=True,
    climb=False,
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_sources(native_path, cache_path):
    if sha(native_path) != NATIVE_SHA or sha(cache_path) != CACHE_SHA:
        raise ValueError("original native/cache archive bytes differ")
    native = json.loads(native_path.read_bytes())
    if native["point_snapshot_sha256"][POINT_PATH] != POINT_SHA:
        raise ValueError("original failed point snapshot differs")
    pins = [r for r in native["native_record_sources"] if r["point_id"] == 13]
    if len(pins) != 12:
        raise ValueError("exact original twelve native records required")
    with zipfile.ZipFile(cache_path) as z:
        if sorted(z.namelist()) != [
            "point_013.json",
            *[f"result_{i:03d}.json" for i in range(12)],
        ]:
            raise ValueError("archive entries differ")
        raw = z.read("point_013.json")
        if hashlib.sha256(raw).hexdigest() != POINT_SHA:
            raise ValueError("original point bytes differ")
        point = json.loads(raw)
        records = []
        for i, pin in enumerate(pins):
            raw = z.read(f"result_{i:03d}.json")
            if hashlib.sha256(raw).hexdigest() != pin["result_sha256"]:
                raise ValueError("original result bytes differ")
            r = json.loads(raw)
            x = np.asarray(r["coordinates"])
            if (
                r["status"] != "complete"
                or r["validation"] != pin["validation"]
                or r["validation"]["genuine_MPI_ranks"] != 32
                or x.shape != (40,)
                or x[-1] != 0
                or not np.isfinite(x).all()
                or hashlib.sha256(x.tobytes()).hexdigest() != r["coordinate_sha256"]
                or any(
                    r["validation"]["input_sha256"][n] != h
                    for n, h in faces.study.CONTRACT.items()
                )
            ):
                raise ValueError("cached record physical/coordinate identity differs")
            records.append(r)
    if (
        point["id"] != 13
        or point["q_A"] != POLICY["q_A"]
        or point["point"] is not None
        or point["status"] != "optimizer_not_converged"
        or point["reserved_attempts"] != 12
        or point["reason"]
        != "mode evaluation attempt budget exhausted; calculator was not called"
    ):
        raise ValueError("the exact original capped nonconverged point is required")
    return point, records


class NewCoordinateRequired(RuntimeError):
    def __init__(self, x):
        super().__init__("new coordinate required; zero-DFT probe has no delegate")
        self.coordinates = np.array(x, copy=True)


class SharedExactCache:
    """Reuse original physical results exactly; delegate budget remains separate."""

    fixed_coordinate_covectors = np.eye(40)[:, -1:]

    def __init__(self, records, delegate=None):
        self.records = {r["coordinate_sha256"]: r for r in records}
        self.delegate = delegate
        if len(self.records) != len(records):
            raise ValueError("unique original cache coordinates required")
        contracts = {r["contract_sha256"] for r in records}
        if len(contracts) != 1:
            raise ValueError("one original physical contract required")
        if delegate is not None and (
            delegate.contract["contract_sha256"] not in contracts
            or not np.array_equal(
                delegate.fixed_coordinate_covectors, self.fixed_coordinate_covectors
            )
        ):
            raise ValueError(
                "new delegate must preserve the exact original physical contract"
            )
        self.requests = []
        self.old_cache_hits = 0

    def __call__(self, x):
        if np.iscomplexobj(x):
            raise ValueError("real physical coordinates required")
        x = np.asarray(x, dtype=float)
        if x.shape != (40,) or not np.isfinite(x).all() or x[-1] != 0:
            raise ValueError(
                "finite full coordinates and exact prescribed control required"
            )
        digest = hashlib.sha256(x.tobytes()).hexdigest()
        row = self.records.get(digest)
        self.requests.append(
            dict(coordinate_sha256=digest, original_cache_hit=row is not None)
        )
        if row is not None:
            if not np.array_equal(x, row["coordinates"]):
                raise ValueError("exact coordinate cache identity differs")
            self.old_cache_hits += 1
            return row["enthalpy_eV"], np.array(row["gradient"], copy=True)
        if self.delegate is None:
            raise NewCoordinateRequired(x)
        return self.delegate(x)


def prepare(native_path, cache_path, solver_path):
    point, records = load_sources(native_path, cache_path)
    data = faces.load_inputs("preserving")
    if not faces.canary.geometry_guard(
        data["chart"], data["basis"], np.array(point["starting_internal_coordinates_A"])
    ):
        raise ValueError("original matched start is not safe")
    radius = math.sqrt(
        POLICY["original_total_internal_radius_A"] ** 2
        - sum(v * v for v in POLICY["q_A"])
    )
    return dict(
        schema_version=1,
        purpose="matched_same_input_optimizer_training_arm_not_independent_forecast",
        policy=json.loads(json.dumps(POLICY)),
        native_prefix_sha256=NATIVE_SHA,
        original_cache_zip_sha256=CACHE_SHA,
        original_point_sha256=POINT_SHA,
        original_twelve_coordinate_sha256=[r["coordinate_sha256"] for r in records],
        original_input_artifact_sha256=list(faces.PINS["preserving"]),
        original_six_input_sha256=dict(faces.study.CONTRACT),
        original_binary_sha256=faces.study.BINARY_SHA,
        original_continuation_parent_id=point["continuation_parent_id"],
        original_requested_start_internal_A=point["starting_internal_coordinates_A"],
        original_first_evaluated_full_A=records[0]["coordinates"],
        released_radius_A=radius,
        full_radius_equivalence_requires_original_guard=True,
        initial_hessian="original_first_step_full36D_symmetrized_internal_matrix_without_clipping",
        original_code_text_sha256=faces.prepare("preserving")["source_text_sha256"],
        matched_worker_text_sha256=faces.study.text_sha(Path(__file__)),
        new_solver_text_sha256=faces.study.text_sha(solver_path),
        old_baseline_new_SCFs=0,
        old_baseline_existing_SCFs=12,
        old_baseline_final_candidate=None,
        cache_preparation_costs_not_free_for_net_speedup=True,
        physical_reserved_geometry_generated_or_read=False,
    )


def validate(config, native_path, cache_path, solver_path):
    if json.dumps(config, sort_keys=True, allow_nan=False) != json.dumps(
        prepare(native_path, cache_path, solver_path), sort_keys=True, allow_nan=False
    ):
        raise ValueError("matched source/start/domain/physical/budget contract differs")
    return load_sources(native_path, cache_path), faces.load_inputs("preserving")


def probe(config, native_path, cache_path, solver_path):
    (point, records), data = validate(config, native_path, cache_path, solver_path)
    spec = importlib.util.spec_from_file_location(
        "vcneb.feasible_relaxation", solver_path
    )
    solver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(solver)
    oracle = SharedExactCache(records)
    evaluate = faces.canary.InternalEvaluator(
        oracle, data["basis"], data["center_energy"], data["center_gradient"]
    )
    plane = faces.canary.ModePlane(
        np.zeros(36),
        data["axes"],
        np.ones(36),
        np.eye(2),
        ("coupled internal1", "coupled internal2"),
        "local length",
        "G3-first-conditional-sheet",
    )
    try:
        solver.relax_orthogonal_in_ball(
            plane,
            POLICY["q_A"],
            evaluate,
            start=np.array(point["starting_internal_coordinates_A"]),
            orthogonal_radius_bound=config["released_radius_A"],
            trial_validator=lambda x: faces.canary.geometry_guard(
                data["chart"], data["basis"], x
            ),
            initial_hessian=data["hessian"],
            initial_trust_radius=0.05,
            gradient_tolerance=0.03,
            max_iterations=200,
        )
    except NewCoordinateRequired as exc:
        if (
            not oracle.requests[0]["original_cache_hit"]
            or oracle.requests[0]["coordinate_sha256"]
            != records[0]["coordinate_sha256"]
        ):
            raise ValueError(
                "native probe must reuse the exact original starting SCF first"
            )
        full = exc.coordinates
        internal = data["basis"].T @ full
        if (
            full[-1] != 0
            or not np.allclose(
                data["axes"].T @ internal, POLICY["q_A"], atol=2e-12, rtol=0
            )
            or not faces.canary.geometry_guard(data["chart"], data["basis"], internal)
        ):
            raise ValueError("novel proposal leaves original Q/control/domain/geometry")
        return dict(
            status="matched_old_start_reused_first_novel_proposal_no_DFT",
            old_cache_hits_before_first_novel=oracle.old_cache_hits,
            requests=oracle.requests,
            first_novel_full_A=full.tolist(),
            first_novel_internal_norm_A=float(np.linalg.norm(internal)),
            first_original_start_static_exact=True,
            new_DFT_calls=0,
            physical_reserved_geometry_generated_or_read=False,
        )
    raise ValueError(
        "zero-DFT probe did not reach a new proposal; no material result inferred"
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("prepare", "probe"))
    p.add_argument("--native", type=Path, required=True)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--solver", type=Path, required=True)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    if args.mode == "prepare":
        if args.config.exists():
            raise FileExistsError("fresh unique material-arm reservation required")
        result = prepare(args.native, args.cache, args.solver)
        args.config.parent.mkdir(parents=True, exist_ok=True)
        faces.study.exclusive_json(args.config, result)
    else:
        if args.output is None or args.output.exists():
            raise FileExistsError("fresh zero-DFT probe output required")
        result = probe(
            json.loads(args.config.read_bytes()), args.native, args.cache, args.solver
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        faces.study.exclusive_json(args.output, result)
    print(
        json.dumps(dict(status="prepared_or_probed_same_input_arm_no_DFT", new_DFT=0))
    )


if __name__ == "__main__":
    main()
