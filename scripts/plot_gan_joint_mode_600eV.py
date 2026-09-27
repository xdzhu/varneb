"""Publication panel from an audited, single-contract 600-eV GaN report.

Panel (a) compares local frozen-block and joint curvatures at one finite-
difference step. Panel (b) shows the geometric endpoint optical projection of
the atomic portion of that joint direction. This is not a TS certificate.
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
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "benchmarks/numerical_integrity/gan_600eV_joint_gamma_bridge_20260928.json"
DEFAULT_PREFIX = ROOT / "paper/VARNEB_CPC/figures/gan_joint_mode_600eV"
COLORS = {
    "blue": "#285D80", "orange": "#C47D45", "gray": "#9DAAB2", "ink": "#263238",
}
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 10.2, "axes.labelsize": 11.2,
    "xtick.labelsize": 10.0, "ytick.labelsize": 10.0,
    "axes.linewidth": 0.85,
    "axes.spines.top": True, "axes.spines.right": True,
})


def load_report(path: Path) -> tuple[dict, str]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    data = json.loads(path.read_text(encoding="utf-8"))
    if (data.get("status")
            != "GaN_same_600eV_joint_endpoint_Gamma_direction_audited_not_TS_certificate"
            or data.get("joint_hessian_step_A") != 0.02
            or data["joint_eigenvalues_first_two_eV_per_A2"][0] >= 0
            or data["joint_eigenvalues_first_two_eV_per_A2"][1] <= 0
            or data["numerical_integrity"]["energy_gradient_max_abs_difference_eV_per_A"] <= 0
            or any("1000" in str(value) for value in data["input_contract_sha256"].values())):
        raise ValueError("600-eV GaN figure report changed")
    blocks = data["joint_block_lowest_curvatures_eV_per_A2"]
    if not (blocks["lowest_frozen_atomic_eV_per_A2"] > 0
            and blocks["lowest_frozen_strain_eV_per_A2"] > 0
            and blocks["lowest_strain_relaxed_atomic_schur_eV_per_A2"] < 0):
        raise ValueError("atomic-strain sign pattern is absent")
    for phase in ("B4", "B1"):
        groups = [item for item in data["phases"][phase]["optical_groups"]
                  if item["fraction"] > 1e-4]
        if (len(groups) != 3 or sorted(len(item["mode_indices"]) for item in groups) != [1, 2, 2]
                or not np.isclose(sum(item["fraction"] for item in groups), 1.0, atol=1e-5)):
            raise ValueError(f"{phase} optical groups changed")
    return data, digest


def source_rows(data: dict, digest: str) -> list[dict]:
    blocks = data["joint_block_lowest_curvatures_eV_per_A2"]
    categories = (
        ("frozen atoms", blocks["lowest_frozen_atomic_eV_per_A2"]),
        ("frozen strain", blocks["lowest_frozen_strain_eV_per_A2"]),
        ("joint", blocks["lowest_joint_two_eV_per_A2"][0]),
        ("relaxed strain", blocks["lowest_strain_relaxed_atomic_schur_eV_per_A2"]),
    )
    rows = [{"panel": "a", "phase": "local candidate", "category": name,
             "value": value, "unit": "eV/A^2", "mode_indices": "",
             "mean_frequency_THz": "", "step_A": 0.02,
             "electronic_contract": "VASP 600 eV", "report_sha256": digest}
            for name, value in categories]
    for phase in ("B4", "B1"):
        groups = [item for item in data["phases"][phase]["optical_groups"]
                  if item["fraction"] > 1e-4]
        for name, item in zip(("lower pair", "single", "upper pair"), groups):
            rows.append({"panel": "b", "phase": phase, "category": name,
                         "value": 100 * item["fraction"], "unit": "% optical direction",
                         "mode_indices": ",".join(map(str, item["mode_indices"])),
                         "mean_frequency_THz": item["mean_frequency_THz"],
                         "step_A": 0.01, "electronic_contract": "VASP 600 eV",
                         "report_sha256": digest})
    return rows


def draw(rows: list[dict]) -> plt.Figure:
    fig = plt.figure(figsize=(7.20, 3.34), facecolor="white")
    ax_a = fig.add_axes((0.100, 0.235, 0.425, 0.615))
    ax_b = fig.add_axes((0.627, 0.235, 0.320, 0.615))
    fig.text(0.050, 0.91, "(a)", fontsize=13, fontweight="normal", color=COLORS["ink"])
    fig.text(0.578, 0.91, "(b)", fontsize=13, fontweight="normal", color=COLORS["ink"])
    keys = ("frozen atoms", "frozen strain", "joint", "relaxed strain")
    labels = ("Atomic\nblock", "Strain\nblock", "Joint\nHessian", "Strain-\nrelaxed")
    values = [next(float(row["value"]) for row in rows
                   if row["panel"] == "a" and row["category"] == key) for key in keys]
    x = np.arange(4)
    ax_a.axhline(0, color="#70777D", lw=0.95, ls="--", zorder=0)
    colors = [COLORS["blue"], COLORS["blue"], COLORS["orange"], COLORS["gray"]]
    ax_a.vlines(x, np.minimum(values, 0), np.maximum(values, 0), color=colors, lw=2.5)
    ax_a.scatter(x, values, s=65, c=colors, edgecolors="white", linewidths=0.7, zorder=4)
    ax_a.set_xlim(-0.4, 3.4)
    ax_a.set_ylim(-18.2, 6.2)
    ax_a.set_xticks(x, labels)
    ax_a.set_yticks([-15, -10, -5, 0, 5])
    ax_a.set_ylabel(r"Lowest curvature (eV Å$^{-2}$)")
    ax_a.tick_params(axis="both", direction="in", top=True, right=True, length=4, pad=5)

    categories = ("lower pair", "single", "upper pair")
    colors_b = (COLORS["blue"], COLORS["gray"], COLORS["orange"])
    for phase, ypos in (("B4", 1.0), ("B1", 0.0)):
        values_b = {row["category"]: float(row["value"]) for row in rows
                    if row["panel"] == "b" and row["phase"] == phase}
        left = 0.0
        for category, color in zip(categories, colors_b):
            width = values_b[category]
            ax_b.barh(ypos, width, left=left, height=0.38, color=color,
                      edgecolor="white", linewidth=0.8)
            if width > 10:
                ax_b.text(left + width / 2, ypos, f"{width:.1f}%", ha="center", va="center",
                          fontsize=10.0, color="white" if category == "lower pair" else COLORS["ink"])
            left += width
    ax_b.set_xlim(0, 100)
    ax_b.set_ylim(-0.46, 1.55)
    ax_b.set_xticks([0, 25, 50, 75, 100])
    ax_b.set_yticks([1, 0], ["B4", "B1"])
    ax_b.set_xlabel("Optical direction fraction (%)")
    ax_b.tick_params(axis="both", direction="in", top=True, right=True, length=4, pad=5)
    handles = [Patch(facecolor=color, edgecolor="white", label=label)
               for color, label in zip(colors_b, ("low pair", "single", "high pair"))]
    ax_b.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.99),
                ncol=3, fontsize=8.8, frameon=True, facecolor="white",
                edgecolor="#8B969C", framealpha=0.82, columnspacing=0.45,
                handlelength=0.8, handletextpad=0.28, borderpad=0.34)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_PREFIX)
    parser.add_argument("--overwrite", action="store_true",
                        help="replace only the explicitly selected generated figure bundle")
    args = parser.parse_args()
    data, digest = load_report(args.report)
    rows = source_rows(data, digest)
    prefix = args.output_prefix
    outputs = [prefix.with_suffix(suffix) for suffix in (".pdf", ".svg", ".png", ".tiff")]
    csv_path = prefix.with_name(prefix.name + "_source_data.csv")
    qa_path = prefix.with_name(prefix.name + "_qa.json")
    if not args.overwrite and any(path.exists() for path in (*outputs, csv_path, qa_path)):
        raise FileExistsError("refusing to overwrite an existing figure bundle")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    figure = draw(rows)
    for path in outputs:
        figure.savefig(path, dpi=600 if path.suffix in {".png", ".tiff"} else None,
                       bbox_inches="tight", facecolor="white")
        if path.suffix == ".svg":
            path.write_text("\n".join(line.rstrip() for line in
                                       path.read_text(encoding="utf-8").splitlines()) + "\n",
                            encoding="utf-8")
    plt.close(figure)
    qa = {
        "report_sha256": digest, "backend": "Python/matplotlib",
        "core_conclusion": "600-eV candidate has one coupling-driven joint negative direction",
        "archetype": "quantitative grid", "n_independent_material_runs": 1,
        "finite_difference_step_A": 0.02,
        "statistics": "Deterministic finite differences; no statistical error bars or hypothesis tests",
        "claim_limit": data["limitations"],
        "source_rows": len(rows),
        "output_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (*outputs, csv_path)},
    }
    qa_path.write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"figure": str(prefix), "source_rows": len(rows)}))


if __name__ == "__main__":
    main()
