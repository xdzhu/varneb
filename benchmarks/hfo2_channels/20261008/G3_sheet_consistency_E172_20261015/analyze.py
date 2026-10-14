"""Post-hoc work/geometry diagnostics of an audited candidate sheet, not PES error.

No new DFT, interpolation, branch selection, stability, independent test or TS.
Uses the existing sampled_path_work rather than introducing an energy partition.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scripts import hfo2_G3_conditional_faces as faces
from vcneb.mode_surface import _orthogonal_directions
from vcneb.path_work import sampled_path_work

ROOT = Path(__file__).resolve().parents[4]
T_AUDIT = (
    ROOT
    / "benchmarks/hfo2_channels/20261008/G3_feasible_execution_E171_20261015/T_complete_native_audit.json"
)
T_SHA = "b8e7836722e3c7cadfaed9771eeb90e1d53097a03fc5d19272109eca0a47c166"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_T():
    if sha(T_AUDIT) != T_SHA:
        raise ValueError("exact original completed native T audit required")
    report = json.loads(T_AUDIT.read_bytes())
    if (
        report["audit_mode"] != "complete"
        or report["channel"] != "T"
        or report["job_id"] != "28862319"
        or report["sacct"].split("|")[1:3] != ["COMPLETED", "0:0"]
        or report["completed_nodes_audited"] != 25
        or report["new_SCFs_reparsed"] != 118
        or report["gradient_candidates"] != 25
        or len(report["candidate_rows"]) != 24
    ):
        raise ValueError(
            "complete original T batch required, not a partial live snapshot"
        )
    cfg = faces.canary.pinned(
        faces.BASE / "G3_conditional_faces_E161_20261015/T_config.json",
        "44f3b34e77a307b5d965893309a11b9ebfd1ad8961b3035904544843c6ec94d6",
    )
    data = faces.validate(cfg, "T")
    origin = data["origin_point"]
    rows = [
        dict(
            id=0,
            q_A=[0.0, 0.0],
            coordinates_internal_A=origin["coordinates_internal_A"],
            coordinates_full_A=origin["coordinates_full_A"],
            energy_eV_cell=origin["energy_eV_cell"],
            full_physical_gradient_eV_A=data["origin_record"]["gradient"],
            released_gradient_norm_eV_A=origin["orthogonal_gradient_norm_eV_A"],
        )
    ]
    for row, spec in zip(report["candidate_rows"], faces.ordered_points()[1:]):
        if (
            row["id"] != spec["id"]
            or row["q_A"] != spec["q_A"]
            or row["candidate"] is None
        ):
            raise ValueError("all original T candidate IDs and Q required")
        if row["replay"] != dict(
            original_order_cache_only_replay_passed=True, new_DFT=0
        ):
            raise ValueError("original trace replay required")
        rows.append(dict(id=row["id"], q_A=row["q_A"], **row["candidate"]))
    return data, rows, report


def validate_nodes(data, rows):
    if len(rows) != 25:
        raise ValueError("exact complete 5x5 T candidate grid required")
    b, a = data["basis"], data["axes"]
    values = []
    for r, spec in zip(rows, faces.ordered_points()):
        if r["id"] != spec["id"] or r["q_A"] != spec["q_A"]:
            raise ValueError("original point order and coordinate labels required")
        x = np.asarray(r["coordinates_internal_A"])
        full = np.asarray(r["coordinates_full_A"])
        gf = np.asarray(r["full_physical_gradient_eV_A"])
        q = np.asarray(r["q_A"])
        if (
            x.shape != (36,)
            or full.shape != (40,)
            or gf.shape != (40,)
            or full[-1] != 0
            or not all(np.isfinite(v).all() for v in (x, full, gf, q))
            or not np.isfinite(r["energy_eV_cell"])
            or not np.allclose(b @ x, full, atol=2e-12, rtol=0)
            or not np.allclose(a.T @ x, q, atol=2e-12, rtol=0)
            or not faces.canary.geometry_guard(data["chart"], b, x)
        ):
            raise ValueError(
                "candidate embedding, fixed Q/control or original guard differs"
            )
        g = b.T @ gf
        gr = g - a @ (a.T @ g)
        residual = float(np.linalg.norm(gr))
        if residual > 0.03 or not np.isclose(
            residual, r["released_gradient_norm_eV_A"], atol=2e-12, rtol=0
        ):
            raise ValueError("full released-gradient target differs")
        values.append(
            dict(
                id=r["id"],
                q=q,
                x=x,
                full=full,
                gf=gf,
                g=g,
                gr=gr,
                energy=r["energy_eV_cell"],
                radius=float(np.linalg.norm(x)),
                residual=residual,
            )
        )
    if len(values) != 25:
        raise ValueError("exact complete 5x5 T candidate grid required")
    return values


def segment(left, right, axes, release):
    x = np.array([left["x"], right["x"]])
    g = np.array([left["g"], right["g"]])
    energy = np.array([left["energy"], right["energy"]])
    w = sampled_path_work(x, g, energy, dict(retained_Q=axes, released_R=release))
    full = sampled_path_work(
        np.array([left["full"], right["full"]]),
        np.array([left["gf"], right["gf"]]),
        energy,
        dict(
            atomic=np.eye(40)[:, :36],
            open_cell=np.eye(40)[:, 36:39],
            prescribed_control=np.eye(40)[:, 39:],
        ),
    )
    dq = right["q"] - left["q"]
    dr = (right["x"] - left["x"]) - axes @ dq
    mg = 0.5 * (left["g"] + right["g"])
    wr = float((mg - axes @ (axes.T @ mg)) @ dr)
    wq = float((axes.T @ mg) @ dq)
    if (
        not np.isclose(
            wq, w["segment_partition_work"]["retained_Q"][0], atol=2e-12, rtol=0
        )
        or not np.isclose(
            wr, w["segment_partition_work"]["released_R"][0], atol=2e-12, rtol=0
        )
        or not np.isclose(wq + wr, full["segment_full_work"][0], atol=2e-12, rtol=0)
    ):
        raise ValueError("independent internal/full native work identity differs")
    bound = 0.5 * (left["residual"] + right["residual"]) * float(np.linalg.norm(dr))
    if abs(wr) > bound + 2e-12:
        raise ValueError("endpoint Cauchy inequality differs")
    return dict(
        ids=[left["id"], right["id"]],
        q_from_A=left["q"].tolist(),
        q_to_A=right["q"].tolist(),
        delta_energy_meV_cell=float(1000 * (energy[1] - energy[0])),
        retained_endpoint_work_meV_cell=1000 * wq,
        released_endpoint_work_meV_cell=1000 * wr,
        retained_only_energy_minus_work_meV_cell=float(
            1000 * (energy[1] - energy[0] - wq)
        ),
        full_energy_minus_trapezoid_work_meV_cell=float(
            1000 * w["segment_energy_minus_work"][0]
        ),
        endpoint_released_work_Cauchy_bound_meV_cell=1000 * bound,
        released_step_norm_A=float(np.linalg.norm(dr)),
        atomic_full_endpoint_work_meV_cell=1000
        * full["segment_partition_work"]["atomic"][0],
        open_cell_full_endpoint_work_meV_cell=1000
        * full["segment_partition_work"]["open_cell"][0],
        prescribed_control_endpoint_work_meV_cell=1000
        * full["segment_partition_work"]["prescribed_control"][0],
        endpoint_Cauchy_bound_is_not_an_energy_error_bound=True,
    )


def diagnose(data, rows):
    nodes = validate_nodes(data, rows)
    a = data["axes"]
    plane = faces.canary.ModePlane(
        np.zeros(36), a, np.ones(36), np.eye(2), ("Q1", "Q2"), "A", "existing-T-sheet"
    )
    release = _orthogonal_directions(plane, None)
    edges = []
    for left in nodes:
        for right in nodes:
            dq = right["q"] - left["q"]
            if (
                np.count_nonzero(dq) == 1
                and np.count_nonzero(np.isclose(dq, 0.1, atol=1e-15, rtol=0)) == 1
            ):
                edges.append(segment(left, right, a, release))
    if len(edges) != 40:
        raise ValueError("exactly forty nearest axial edges required")
    lookup = {tuple(np.rint(n["q"] / 0.1).astype(int)): n for n in nodes}
    loops = []
    for i in range(-2, 2):
        for j in range(-2, 2):
            ordered = [
                lookup[k]
                for k in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1), (i, j))
            ]
            w = sampled_path_work(
                np.array([n["x"] for n in ordered]),
                np.array([n["g"] for n in ordered]),
                np.array([n["energy"] for n in ordered]),
                dict(retained_Q=a, released_R=release),
            )
            loops.append(
                dict(
                    ids=[n["id"] for n in ordered],
                    lower_left_Q_A=(np.array([i, j]) * 0.1).tolist(),
                    retained_endpoint_circulation_meV_cell=1000
                    * w["cumulative_partition_work"]["retained_Q"][-1],
                    released_endpoint_circulation_meV_cell=1000
                    * w["cumulative_partition_work"]["released_R"][-1],
                    full_endpoint_circulation_meV_cell=1000
                    * w["cumulative_full_work"][-1],
                    raw_closed_loop_energy_change_meV_cell=1000
                    * w["cumulative_energy_change"][-1],
                )
            )
    e = np.array([n["energy"] for n in nodes])
    norm = np.array([n["radius"] for n in nodes])

    def extremum(key):
        return max(edges, key=lambda r: abs(r[key]))

    summary = dict(
        candidate_nodes=25,
        nearest_axial_edges=40,
        elementary_closed_loops=16,
        source_new_SCFs=118,
        new_DFT_calls=0,
        source_curvature_SCFs_both_centers_additional=296,
        relative_to_released_origin_energy_range_meV_cell=[
            float(1000 * (e.min() - e[0])),
            float(1000 * (e.max() - e[0])),
        ],
        total_radius_range_A=[float(norm.min()), float(norm.max())],
        minimum_original_radius_margin_A=float(0.5 - norm.max()),
        released_gradient_range_eV_A=[
            min(n["residual"] for n in nodes),
            max(n["residual"] for n in nodes),
        ],
        largest_retained_only_energy_work_difference=extremum(
            "retained_only_energy_minus_work_meV_cell"
        ),
        largest_full_endpoint_quadrature_difference=extremum(
            "full_energy_minus_trapezoid_work_meV_cell"
        ),
        largest_released_endpoint_work=extremum("released_endpoint_work_meV_cell"),
        largest_released_displacement=extremum("released_step_norm_A"),
        largest_retained_loop_circulation=max(
            loops, key=lambda r: abs(r["retained_endpoint_circulation_meV_cell"])
        ),
        largest_full_loop_circulation=max(
            loops, key=lambda r: abs(r["full_endpoint_circulation_meV_cell"])
        ),
    )
    return dict(
        summary=summary,
        nodes=[
            dict(
                id=n["id"],
                q_A=n["q"].tolist(),
                energy_relative_origin_meV_cell=1000 * (n["energy"] - e[0]),
                total_internal_radius_A=n["radius"],
                released_gradient_norm_eV_A=n["residual"],
            )
            for n in nodes
        ],
        edges=edges,
        loops=loops,
    )


def analyze():
    data, rows, source = load_T()
    result = diagnose(data, rows)
    result.update(
        status="posthoc_complete_candidate_sheet_work_diagnostic_not_PES_error_or_prediction",
        source_complete_T_native_sha256=T_SHA,
        source_job_id=source["job_id"],
        source_audit_mode=source["audit_mode"],
        source_six_physical_inputs=faces.study.CONTRACT,
        source_pressure_GPa=0.0,
        source_physical_strain=0.0,
        formula_units=4,
        energy_work_unit="meV per Hf4O8 cell; divide by4 for meV per HfO2",
        derivative_chart="original fixed-strain36D translation-free internal length chart; unchanged joint metric",
        source_text_sha256={
            n: faces.study.text_sha(ROOT / n)
            for n in (
                "vcneb/path_work.py",
                "vcneb/mode_surface.py",
                "scripts/hfo2_G3_conditional_faces.py",
                "benchmarks/hfo2_channels/20261008/G3_sheet_consistency_E172_20261015/analyze.py",
            )
        },
        source_force_ASR_or_geometry_or_physics_changed=False,
        source_native_SCFs_newly_reparsed=False,
        physical_error_bound_or_stability_or_unique_branch_or_TS_or_independent_prediction_certified=False,
        physical_reserved_geometry_generated_or_read=False,
        additional_training_points=0,
        scheduler_writes=0,
        limitations=[
            "post-hoc training diagnostic, not an independent forecast or novelty theorem",
            "finite residuals and finite-spacing endpoint trapezoids are retained, not set to zero",
            "Cauchy inequality bounds only endpoint released trapezoid work, not energy or barrier error",
            "closed-loop endpoint circulation is a quadrature/branch-consistency diagnostic, not nonconservative DFT forces",
            "neighbor displacement size alone does not prove smoothness, a branch jump, stability or uniqueness",
            "signed atomic/cell work is not a state-function mode energy decomposition",
            "existing raw native audit is reused; no additional raw SCF validation is claimed",
        ],
    )
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError("fresh analysis output required")
    result = analyze()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")
    print(json.dumps(result["summary"]))


if __name__ == "__main__":
    main()
