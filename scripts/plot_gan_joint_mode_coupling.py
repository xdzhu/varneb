"""Plot the audited GaN atom--strain saddle coupling and endpoint-mode bridge.

Figure contract
---------------
Core conclusion: the local B4-to-B1 joint negative curvature comes from
atom--strain coupling although its frozen atomic and strain blocks are positive.
Archetype: quantitative grid, with the curvature-sign comparison as hero panel.
Backend: Python/matplotlib only; 183-mm vector PDF/SVG and 600-dpi raster.
Panel (a): two finite-difference steps of the 1000-eV local enthalpy Hessian.
Panel (b): the same unstable *direction* projected on 600-eV endpoint Gamma
optical subspaces; fractions are geometric, not energy/barrier contributions.
Review risk: the two electronic settings must not be mixed into an energy
comparison, and neither endpoint Gamma basis is a TS phonon calculation.
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


REPORT_SHA256 = "fc530fe491027c22eb4fe573007abf0ff2f6728b93c668b5d09ced1ee5baa20c"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "benchmarks" / "numerical_integrity" / "gan_ts_endpoint_gamma_bridge_2026-09-27.json"
DEFAULT_PREFIX = ROOT / "paper" / "VARNEB_CPC" / "figures" / "gan_joint_mode_coupling"

COLORS = {
    "step_small": "#285D80",
    "step_large": "#8FB5C8",
    "lower_pair": "#285D80",
    "single": "#B8C1C9",
    "upper_pair": "#C47D45",
    "ink": "#263238",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 9.5, "axes.labelsize": 10.5,
    "xtick.labelsize": 9.2, "ytick.labelsize": 9.2,
    "axes.linewidth": 0.85,
    "axes.spines.top": True, "axes.spines.right": True,
})


def load_evidence(path: Path) -> tuple[dict, str]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != REPORT_SHA256:
        raise ValueError("GaN coupling evidence differs from the audited pinned report")
    data = json.loads(path.read_text(encoding="utf-8"))
    if (data.get("status") != "audited_GaN_joint_unstable_direction_endpoint_Gamma_atomic_projection"
            or "1000 eV joint TS candidate and 600 eV endpoint Gamma bases" not in data.get("scope", "")):
        raise ValueError("GaN coupling report has a different scientific scope")
    blocks = data["joint_hessian_atomic_strain_block_diagnostics"]
    for step in ("0p01", "0p02"):
        record = blocks[step]
        if not (record["lowest_frozen_atomic_eV_per_A2"] > 0
                and record["lowest_frozen_strain_eV_per_A2"] > 0
                and record["lowest_joint_two_eV_per_A2"][0] < 0
                and record["lowest_strain_relaxed_atomic_schur_eV_per_A2"] < 0):
            raise ValueError("the reported atom--strain sign pattern is absent")
    if data["joint_unstable_direction_cross_step_absolute_overlap"] < 0.99:
        raise ValueError("unstable direction is not robust across steps")
    for phase in ("B4", "B1"):
        groups = data["phases"][phase]["d0.01"]["0p01"]["optical_groups"]
        selected = [item for item in groups if item["fraction"] > 1e-4]
        if (len(selected) != 3 or [len(item["mode_indices"]) for item in selected] != [2, 1, 2]
                or abs(sum(item["fraction"] for item in selected) - 1.0) > 1e-5):
            raise ValueError(f"{phase} optical grouping changed")
    return data, digest


def source_rows(data: dict, digest: str) -> list[dict[str, str | float]]:
    rows: list[dict[str, str | float]] = []
    fields = (
        ("frozen atoms", "lowest_frozen_atomic_eV_per_A2"),
        ("frozen strain", "lowest_frozen_strain_eV_per_A2"),
        ("joint", "lowest_joint_two_eV_per_A2"),
        ("relaxed strain", "lowest_strain_relaxed_atomic_schur_eV_per_A2"),
    )
    for step, step_A in (("0p01", 0.01), ("0p02", 0.02)):
        record = data["joint_hessian_atomic_strain_block_diagnostics"][step]
        for label, key in fields:
            value = record[key][0] if key == "lowest_joint_two_eV_per_A2" else record[key]
            rows.append({
                "panel": "a", "category": label, "phase": "local TS candidate",
                "step_A": step_A, "frequency_THz": "", "mode_indices": "",
                "value": value, "unit": "eV/A^2", "electronic_contract": "VASP 1000 eV",
                "report_sha256": digest,
            })
    for phase in ("B4", "B1"):
        groups = [item for item in data["phases"][phase]["d0.01"]["0p01"]["optical_groups"]
                  if item["fraction"] > 1e-4]
        for label, item in zip(("lower pair", "single", "upper pair"), groups):
            rows.append({
                "panel": "b", "category": label, "phase": phase,
                "step_A": 0.01, "frequency_THz": item["mean_frequency_THz"],
                "mode_indices": ",".join(map(str, item["mode_indices"])),
                "value": 100.0 * item["fraction"], "unit": "% atomic optical direction",
                "electronic_contract": "VASP 600 eV endpoint Gamma basis",
                "report_sha256": digest,
            })
    return rows


def draw_figure(rows: list[dict[str, str | float]]) -> plt.Figure:
    fig = plt.figure(figsize=(7.20, 3.25), facecolor="white")
    # Both plot boxes share the same top and bottom despite different y axes.
    ax_a = fig.add_axes((0.095, 0.225, 0.43, 0.625))
    ax_b = fig.add_axes((0.625, 0.225, 0.32, 0.625))
    fig.text(0.050, 0.905, "(a)", fontsize=12, fontweight="normal", color=COLORS["ink"])
    fig.text(0.580, 0.905, "(b)", fontsize=12, fontweight="normal", color=COLORS["ink"])

    labels = ["Atomic\nblock", "Strain\nblock", "Joint\nHessian", "Strain-\nrelaxed"]
    keys = ["frozen atoms", "frozen strain", "joint", "relaxed strain"]
    for step, dx, color, marker, legend in (
        (0.01, -0.085, COLORS["step_small"], "o", r"$\delta=0.01$ Å"),
        (0.02, +0.085, COLORS["step_large"], "s", r"$\delta=0.02$ Å"),
    ):
        observed = {row["category"]: float(row["value"]) for row in rows
                    if row["panel"] == "a" and row["step_A"] == step}
        x = np.arange(4, dtype=float) + dx
        y = np.asarray([observed[key] for key in keys])
        ax_a.vlines(x, np.minimum(0, y), np.maximum(0, y), color=color, lw=2.1, alpha=0.9)
        ax_a.scatter(x, y, s=43, marker=marker, color=color, edgecolor="white",
                     linewidth=0.65, zorder=4, label=legend)
    ax_a.axhline(0.0, color="#70777D", lw=0.9, ls="--", zorder=0)
    ax_a.set_xlim(-0.4, 3.4)
    ax_a.set_ylim(-18.2, 6.2)
    ax_a.set_xticks(range(4), labels)
    ax_a.set_yticks([-15, -10, -5, 0, 5])
    ax_a.set_ylabel(r"Lowest curvature (eV Å$^{-2}$)")
    ax_a.tick_params(axis="both", which="both", direction="in", top=True, right=True,
                     length=3.4, pad=4)
    ax_a.legend(loc="upper right", fontsize=8.7, ncol=1, frameon=True,
                facecolor="white", edgecolor="#8B969C", framealpha=0.85,
                borderpad=0.35, handletextpad=0.35)

    categories = ("lower pair", "single", "upper pair")
    colors = (COLORS["lower_pair"], COLORS["single"], COLORS["upper_pair"])
    for phase, ypos in (("B4", 1.0), ("B1", 0.0)):
        values = {row["category"]: float(row["value"]) for row in rows
                  if row["panel"] == "b" and row["phase"] == phase}
        left = 0.0
        for category, color in zip(categories, colors):
            width = values[category]
            ax_b.barh(ypos, width, left=left, height=0.39, color=color,
                      edgecolor="white", linewidth=0.8)
            if width > 10.0:
                ax_b.text(left + width / 2, ypos, f"{width:.1f}%", ha="center", va="center",
                          fontsize=9.0, color="white" if category == "lower pair" else COLORS["ink"])
            left += width
    ax_b.set_xlim(0, 100)
    ax_b.set_ylim(-0.45, 1.55)
    ax_b.set_xticks([0, 25, 50, 75, 100])
    ax_b.set_yticks([1.0, 0.0], ["B4", "B1"])
    ax_b.set_xlabel("Optical direction fraction (%)")
    ax_b.tick_params(axis="both", which="both", direction="in", top=True, right=True,
                     length=3.4, pad=4)
    handles = [Patch(facecolor=color, edgecolor="white", label=label)
               for color, label in zip(colors, ("low pair", "single", "high pair"))]
    ax_b.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.99),
                ncol=3, fontsize=8.0, frameon=True, facecolor="white",
                edgecolor="#8B969C", framealpha=0.85, columnspacing=0.48,
                handlelength=0.9, handletextpad=0.3, borderpad=0.35)
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_PREFIX)
    parser.add_argument("--overwrite", action="store_true",
                        help="replace only the explicitly selected generated figure bundle")
    args = parser.parse_args()
    data, digest = load_evidence(args.report)
    rows = source_rows(data, digest)
    prefix = args.output_prefix
    outputs = [prefix.with_suffix(suffix) for suffix in
               (".pdf", ".svg", ".png", ".tiff")]
    csv_path = prefix.with_name(prefix.name + "_source_data.csv")
    qa_path = prefix.with_name(prefix.name + "_qa.json")
    if not args.overwrite and any(path.exists() for path in (*outputs, csv_path, qa_path)):
        raise FileExistsError("refusing to overwrite an existing figure bundle")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    figure = draw_figure(rows)
    for path in outputs:
        figure.savefig(path, dpi=600 if path.suffix in {".png", ".tiff"} else None,
                       bbox_inches="tight", facecolor="white")
        if path.suffix == ".svg":
            # Matplotlib emits trailing spaces inside path data; normalize the
            # generated XML so the vector artifact passes repository hygiene.
            svg = path.read_text(encoding="utf-8")
            path.write_text("\n".join(line.rstrip() for line in svg.splitlines()) + "\n",
                            encoding="utf-8")
    plt.close(figure)
    qa = {
        "report_sha256": digest,
        "backend": "Python/matplotlib",
        "archetype": "quantitative grid",
        "core_conclusion": "Local negative GaN B4-to-B1 curvature arises from atomic-strain coupling",
        "numeric_steps_A": [0.01, 0.02],
        "joint_direction_absolute_overlap": data["joint_unstable_direction_cross_step_absolute_overlap"],
        "source_rows": len(rows),
        "n_independent_material_runs": 1,
        "statistics": "Deterministic finite differences; no statistical error bars or hypothesis tests",
        "claim_limit": "1000-eV local Hessian and 600-eV endpoint Gamma modes are compared only as directions, not energies or barrier fractions",
        "output_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (*outputs, csv_path)},
    }
    qa_path.write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"figure": str(prefix), "source_rows": len(rows),
                      "joint_direction_absolute_overlap": qa["joint_direction_absolute_overlap"]}, indent=2))


if __name__ == "__main__":
    main()
