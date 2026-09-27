"""Plot audited BTO fixed-Q candidates without implying a continuous PES.

Only four independently calculated coordinates are shown.  The open degrees
of freedom are relaxed, but branch continuity and a 2D interpolant are not
certified.  The frozen Q_y=0 branch is a *relaxed constrained branch*, not the
frozen cubic-cell contour plotted in the separate exploratory figure.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "paper/VARNEB_CPC/figures/bto_conditional_four_point_stage_2026-09-27"


def _read(relative: str) -> tuple[dict, str]:
    path = ROOT / relative
    data = path.read_bytes()
    return json.loads(data), hashlib.sha256(data).hexdigest()


def _replayed_point(qz: float, qx: float, result_path: str, replay_path: str) -> dict:
    result, result_hash = _read(result_path)
    replay, replay_hash = _read(replay_path)
    if not (result["status"] == "orthogonal_gradient_and_stress_converged_curvature_unchecked"
            and result["stress_target_passed"]
            and np.allclose(result["q1_q2_sqrt_amu_A"], [qz, qx], atol=1e-12)
            and replay["status"] == "all_three_branch_outcomes_reproduced_from_raw_audited_cache"
            and replay["source_sha256"]["summary"] == result_hash
            and replay["selected_start"] == result["selected_start"]
            and len(replay["branch_outcomes"]) == 3):
        raise ValueError("BTO fixed-Q result/replay contract failed")
    branches = replay["branch_outcomes"]
    chosen = branches[result["selected_start"]]
    if abs(chosen["energy_minus_c_eV_per_BTO"] - result["energy_minus_c_eV_per_BTO"]) > 1e-8:
        raise ValueError("selected result and audited-cache replay disagree")
    return {
        "Q_z_sqrt_amu_A": qz, "Q_x_sqrt_amu_A": qx,
        "selected_E_minus_C_eV_per_BTO": float(chosen["energy_minus_c_eV_per_BTO"]),
        "relaxed_Qy_zero_E_minus_C_eV_per_BTO": float(branches[0]["energy_minus_c_eV_per_BTO"]),
        "selected_Q_y_sqrt_amu_A": float(chosen["third_soft_y_amplitude_sqrt_amu_A"]),
        "result_sha256": result_hash, "replay_sha256": replay_hash,
    }


def _off_axis_point() -> dict:
    evidence, digest = _read(
        "benchmarks/numerical_integrity/bto_q060_q030_three_branch_evidence_2026-09-27.json"
    )
    if (evidence["fixed_Q_z_Q_x_sqrt_amu_A"] != [0.6, 0.3]
            or evidence["audit"]["raw_DFT_points_individually_validated"] != 76
            or len(evidence["branches"]) != 3):
        raise ValueError("off-axis branch evidence is incomplete")
    branches = evidence["branches"]
    chosen = min(branches, key=lambda branch: branch["energy_minus_c_eV_per_BTO"])
    return {
        "Q_z_sqrt_amu_A": 0.6, "Q_x_sqrt_amu_A": 0.3,
        "selected_E_minus_C_eV_per_BTO": float(chosen["energy_minus_c_eV_per_BTO"]),
        "relaxed_Qy_zero_E_minus_C_eV_per_BTO": float(branches[0]["energy_minus_c_eV_per_BTO"]),
        "selected_Q_y_sqrt_amu_A": float(chosen["final_Q_y_sqrt_amu_A"]),
        "result_sha256": digest, "replay_sha256": evidence["sha256"]["cache_only_replay"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-csv", type=Path,
                        help="plot from the committed source table in a fresh checkout")
    parser.add_argument("--output-prefix", type=Path, default=OUTPUT)
    args = parser.parse_args()
    output = args.output_prefix
    files = [output.with_suffix(extension) for extension in (".png", ".pdf", ".svg")]
    csv_path = Path(str(output) + "_source_data.csv")
    qa_path = Path(str(output) + "_qa.json")
    if any(path.exists() for path in [*files, csv_path, qa_path]):
        raise FileExistsError("refusing to overwrite a stage figure or its source data")
    points = [
        _replayed_point(
            0.6, 0.0,
            "outputs/batio3_t_to_c_pbe100_dzp10au/bto_transverse_soft_conditional_q060_result_job27783331.json",
            "outputs/batio3_t_to_c_pbe100_dzp10au/bto_transverse_soft_q060_branch_replay_job27783331.json",
        ),
        _off_axis_point(),
        _replayed_point(
            0.9, 0.0,
            "benchmarks/numerical_integrity/bto_q090_hf_2026-09-27/conditional_q090_q000_result.json",
            "benchmarks/numerical_integrity/bto_q090_hf_2026-09-27/replay-q090_q000-final.json",
        ),
        _replayed_point(
            0.9, 0.3,
            "benchmarks/numerical_integrity/bto_q090_hf_2026-09-27/conditional_q090_q030_result.json",
            "benchmarks/numerical_integrity/bto_q090_hf_2026-09-27/replay-q090_q030-final.json",
        ),
    ] if args.source_csv is None else []
    if args.source_csv is not None:
        with args.source_csv.open(newline="", encoding="utf-8") as stream:
            points = list(csv.DictReader(stream))
        if ({(float(row["Q_z_sqrt_amu_A"]), float(row["Q_x_sqrt_amu_A"]))
             for row in points} != {(0.6, 0.0), (0.6, 0.3), (0.9, 0.0), (0.9, 0.3)}
                or len(points) != 4):
            raise ValueError("source table is not the four audited coordinates")
        numeric = ("Q_z_sqrt_amu_A", "Q_x_sqrt_amu_A",
                   "selected_E_minus_C_eV_per_BTO",
                   "relaxed_Qy_zero_E_minus_C_eV_per_BTO",
                   "selected_Q_y_sqrt_amu_A")
        for row in points:
            for key in numeric:
                row[key] = float(row[key])
            if len(row["result_sha256"]) != 64 or len(row["replay_sha256"]) != 64:
                raise ValueError("source table lacks raw-evidence provenance hashes")
    for point in points:
        lowering = 1000 * (point["relaxed_Qy_zero_E_minus_C_eV_per_BTO"]
                           - point["selected_E_minus_C_eV_per_BTO"])
        if lowering <= 0 or abs(point["selected_Q_y_sqrt_amu_A"]) < 0.5:
            raise ValueError("selected branch does not reproduce the lower off-plane solution")
        point["selected_branch_lowering_meV_per_BTO"] = lowering

    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 11, "axes.labelsize": 12, "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5, "legend.fontsize": 10,
        "axes.linewidth": 0.9, "svg.fonttype": "none", "pdf.fonttype": 42,
    })
    colors = {0.0: "#245c9b", 0.3: "#c56e31"}
    fig, (ax, bars) = plt.subplots(1, 2, figsize=(10.0, 4.6))
    fig.subplots_adjust(left=0.095, right=0.975, bottom=0.18, top=0.85, wspace=0.44)
    for point in points:
        qz, qx = point["Q_z_sqrt_amu_A"], point["Q_x_sqrt_amu_A"]
        ax.scatter(qz, qx, s=255, marker="s", color=colors[qx],
                   edgecolors="#22313d", linewidths=0.8, zorder=3)
        energy = 1000 * point["selected_E_minus_C_eV_per_BTO"]
        ax.annotate(f"{energy:.1f}", (qz, qx), xytext=(0, 15 if qx == 0.0 else -24),
                    textcoords="offset points", ha="center", va="center", fontsize=10.5)
    ax.set(xlim=(0.51, 0.99), ylim=(-0.08, 0.38),
           xticks=[0.6, 0.75, 0.9], yticks=[0.0, 0.15, 0.3],
           xlabel=r"$Q_z$ ($\sqrt{\mathrm{amu}}$ Å)",
           ylabel=r"$Q_x$ ($\sqrt{\mathrm{amu}}$ Å)")
    ax.text(0.75, 0.15, "unsampled", ha="center", va="center", fontsize=10,
            color="#666666", style="italic")
    ax.text(0.73, -0.065, r"labels: $(E-E_C)$ in meV/BTO", ha="center", fontsize=9.5)

    ordered = sorted(points, key=lambda point: (point["Q_z_sqrt_amu_A"],
                                                point["Q_x_sqrt_amu_A"]))
    ys = np.arange(len(ordered))[::-1]
    widths = [point["selected_branch_lowering_meV_per_BTO"] for point in ordered]
    bars.barh(ys, widths, height=0.56,
              color=[colors[point["Q_x_sqrt_amu_A"]] for point in ordered],
              edgecolor="#22313d", linewidth=0.7)
    for y, width in zip(ys, widths):
        bars.text(width + 1.0, y, f"{width:.1f}", va="center", fontsize=10.5)
    bars.set(yticks=ys,
             yticklabels=[f"({p['Q_z_sqrt_amu_A']:.1f}, {p['Q_x_sqrt_amu_A']:.1f})"
                          for p in ordered],
             xlim=(0, 70), xlabel=r"Lowering from relaxed $Q_y=0$ (meV/BTO)",
             ylabel=r"Fixed $(Q_z,Q_x)$ ($\sqrt{\mathrm{amu}}$ Å)")
    for panel in (ax, bars):
        panel.tick_params(direction="in", top=True, right=True, length=4)
        for spine in panel.spines.values():
            spine.set_visible(True)
    fig.text(0.095, 0.88, "(a)", fontsize=15, fontweight="normal")
    fig.text(0.595, 0.88, "(b)", fontsize=15, fontweight="normal")
    output.parent.mkdir(parents=True, exist_ok=True)
    for path in files:
        fig.savefig(path, dpi=350, facecolor="white")
    plt.close(fig)

    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=points[0].keys())
        writer.writeheader()
        writer.writerows(ordered)
    qa = {
        "kind": "BTO_four_measured_fixed_Q_branch_candidates_not_continuous_PES_or_barrier",
        "n_measured_Q_points": 4,
        "interpolation_used": False,
        "phonon_supercell": [1, 1, 1],
        "electronic_kpoints": [4, 4, 4],
        "ABACUS_ecutwfc_Ry": 100,
        "ABACUS_orbitals": "10 au DZP",
        "source_csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "claim_limit": "Stationary branch candidates at four fixed-Q points; not a certified 2D conditional PES, a T-to-C barrier, or global branch continuity.",
    }
    qa_path.write_text(json.dumps(qa, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"figure": str(files[0]), "lowering_meV_per_BTO": widths}, indent=2))


if __name__ == "__main__":
    main()
