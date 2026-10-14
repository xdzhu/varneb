"""Matched frozen/released training ablation; no fitted coefficients or new DFT."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from scripts import hfo2_G3_conditional_faces as faces
from vcneb.quadratic_reduction import condition_quadratic_energy

CONFIG = (
    ROOT
    / "benchmarks/hfo2_channels/20261008/G3_frozen_R_controls_E166_20261015/config.json"
)
CONFIG_SHA = "5ffb1fe0d707aec5d0a6546a44e0ad0002cd0f0280f1ba17ca2d955ac2d4c930"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def decompose_shape(
    model_frozen,
    model_released,
    dft_frozen,
    dft_released,
    *,
    model_released_origin,
    dft_released_origin,
):
    """All inputs are energy changes relative to the same raw original center."""
    values = np.array(
        [
            model_frozen,
            model_released,
            dft_frozen,
            dft_released,
            model_released_origin,
            dft_released_origin,
        ],
        dtype=float,
    )
    if not np.isfinite(values).all():
        raise ValueError("finite same-reference energies required")
    mf, mr, df, dr, m0, d0 = values
    er = (mr - m0) - (dr - d0)
    ef = mf - df
    seen = (dr - df) - d0
    predicted = (mr - mf) - m0
    missing = predicted - seen
    defect = er - ef - missing
    if abs(defect) > 2e-11:
        raise ValueError("exact signed shape-error identity differs")
    return dict(
        frozen_shape_error_meV_cell=float(1000 * ef),
        released_shape_error_meV_cell=float(1000 * er),
        DFT_released_minus_frozen_meV_cell=float(1000 * (dr - df)),
        DFT_release_change_relative_origin_meV_cell=float(1000 * seen),
        model_release_change_relative_origin_meV_cell=float(1000 * predicted),
        missed_release_shape_change_meV_cell=float(1000 * missing),
        error_decomposition_identity_defect_meV_cell=float(1000 * defect),
        signed_DFT_release_not_clipped=True,
    )


def analyze(summary, native):
    controls = native["controls"]
    if (
        summary["status"]
        != "eight_frozen_R_statics_completed_not_PES_or_independent_forecast"
        or summary["job_id"] != "28865656"
        or summary["new_SCFs"] != 8
        or summary["total_attempts"] != 8
        or summary["source_original_centers_and_released_origins_reused"] is not True
        or summary[
            "independent_material_prediction_or_released_stability_or_TS_certified"
        ]
        is not False
        or summary["physical_reserved_geometry_generated_or_read"] is not False
        or len(summary["rows"]) != 8
        or controls["job_id"] != "28865656"
        or controls["original_native_SCFs_reparsed"] != 8
        or controls["sources"] != summary["rows"]
        or native["original_runtime_files"] != 380
        or native["original_native_SCFs"] != 12
        or native["old_cache_hits_before_first_novel"] != 7
        or not native["controls_required"]
        or native["new_DFT"] != 0
    ):
        raise ValueError("complete exact eight-control native admission proof required")
    config = faces.canary.pinned(CONFIG, CONFIG_SHA)
    observed = summary["rows"]
    rows = []
    for channel in ("preserving", "T"):
        data = faces.load_inputs(channel)
        a = data["axes"]
        b = data["basis"]
        g = b.T @ data["center_gradient"]
        r = np.linalg.svd(a.T, full_matrices=True)[2][2:].T
        center = data["center_energy"]
        origin = data["origin_point"]["energy_eV_cell"]
        selected = [x for x in observed if x["channel"] == channel]
        if [x["point_id"] for x in selected] != [1, 2, 3, 4]:
            raise ValueError(
                "all four frozen axial controls required in original order"
            )
        for point, expected in zip(selected, config["points"][channel]):
            if (
                point["point_id"] != expected["id"]
                or point["q_A"] != expected["q_A"]
                or point["new_SCFs"] != 1
                or point["attempts"] != 1
                or point["native_validation"]["genuine_MPI_ranks"] != 32
                or any(
                    point["native_validation"]["input_sha256"][n] != h
                    for n, h in faces.study.CONTRACT.items()
                )
                or not np.array_equal(
                    point["coordinates_full_A"], expected["coordinates_full_A"]
                )
                or point["coordinates_full_A"][-1] != 0
                or point["matched_released_energy_eV_cell"]
                != expected["matched_released_energy_eV_cell"]
                or not np.isclose(
                    point["frozen_minus_released_meV_cell"],
                    1000
                    * (
                        point["frozen_energy_eV_cell"]
                        - point["matched_released_energy_eV_cell"]
                    ),
                    atol=2e-8,
                    rtol=0,
                )
            ):
                raise ValueError(
                    "fixed-Q/input/native source or signed release identity differs"
                )
        for block in data["curvature"]["blocks"]:
            raw = np.asarray(block["raw_projected_eV_A2"])[:36, :36]
            h = 0.5 * (raw + raw.T)
            model = condition_quadratic_energy(h, g, center, a, r, stability_floor=0.0)
            points = []
            for p in selected:
                q = np.asarray(p["q_A"])
                x = a @ q
                mf = float(g @ x + 0.5 * x @ h @ x)
                mr = model.energy_change(q)
                df = p["frozen_energy_eV_cell"] - center
                dr = p["matched_released_energy_eV_cell"] - center
                points.append(
                    dict(
                        point_id=p["point_id"],
                        q_A=p["q_A"],
                        DFT_frozen_energy_change_eV_cell=df,
                        DFT_released_energy_change_eV_cell=dr,
                        model_frozen_energy_change_eV_cell=mf,
                        model_released_energy_change_eV_cell=mr,
                        **decompose_shape(
                            mf,
                            mr,
                            df,
                            dr,
                            model_released_origin=model.reference_energy_change,
                            dft_released_origin=origin - center,
                        ),
                    )
                )
            rows.append(
                dict(
                    channel=channel,
                    curvature_probe_step_A=block["step_A"],
                    points=points,
                    raw_origin_DFT_release_meV_cell=1000 * (origin - center),
                    raw_origin_model_release_meV_cell=1000
                    * model.reference_energy_change,
                    max_abs_frozen_shape_error_meV_cell=max(
                        abs(p["frozen_shape_error_meV_cell"]) for p in points
                    ),
                    max_abs_released_shape_error_meV_cell=max(
                        abs(p["released_shape_error_meV_cell"]) for p in points
                    ),
                    max_abs_DFT_release_shape_change_meV_cell=max(
                        abs(p["DFT_release_change_relative_origin_meV_cell"])
                        for p in points
                    ),
                    max_abs_missed_release_shape_change_meV_cell=max(
                        abs(p["missed_release_shape_change_meV_cell"]) for p in points
                    ),
                )
            )
    return dict(
        status="matched_seen_frozen_released_shape_ablation_not_independent_prediction",
        rows=rows,
        all_four_axial_points_both_centers_and_both_probe_steps_retained=True,
        formula_units=4,
        energy_unit="meV per Hf4O8cell; divide by4for meV per HfO2",
        original_config_sha256=CONFIG_SHA,
        original_eight_frozen_statics=8,
        completed_static_call_proxy_core_hours=sum(
            p["native_validation"]["elapsed_seconds"] for p in observed
        )
        * 32
        / 3600,
        completed_static_attempts=sum(p["attempts"] for p in observed),
        new_DFT_calls=0,
        original_curvature_both_centers_SCFs_additional=296,
        all_models_same_original_affine_gradient_Hessian_no_fit=True,
        R_zero_does_not_mean_fully_fixed_cell=True,
        comparison_is_not_full_B2_B5_or_T_Cmma_control_matrix=True,
        physical_error_stability_TS_or_forecast_advantage_certified=False,
        physical_reserved_geometry_generated_or_read=False,
        source_text_sha256={
            n: faces.study.text_sha(ROOT / n)
            for n in (
                "vcneb/quadratic_reduction.py",
                "scripts/hfo2_G3_conditional_faces.py",
                "benchmarks/hfo2_channels/20261008/G3_frozen_release_E173_20261015/analyze.py",
            )
        },
        limitations=[
            "training-only four axial Q controls, not all25surface points",
            "signed decomposition identifies frozen versus released shape discrepancies, not a unique microscopic cause",
            "formal stable release at original raw matrices does not certify post-relaxation stability",
            "both finite-difference probe steps are retained, not interpreted as physical uncertainty bounds",
            "censored optimizer baseline is not a denominator for speedup",
        ],
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--native", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError("fresh ablation output required")
    summary = json.loads(args.summary.read_bytes())
    native = json.loads(args.native.read_bytes())
    if sha(args.summary) != native["controls"]["controls_summary_sha256"]:
        raise ValueError("exact original control summary bytes required")
    result = analyze(summary, native)
    result.update(
        controls_summary_sha256=sha(args.summary),
        native_admission_proof_sha256=sha(args.native),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")
    print(
        json.dumps(
            [{k: v for k, v in r.items() if k != "points"} for r in result["rows"]]
        )
    )


if __name__ == "__main__":
    main()
