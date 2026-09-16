"""Create editable source-data plots for the BTO Gamma-mode path analysis.

Core conclusion: the atomic component of the seven-image BTO T-to-C VCNEB
path is principally aligned with the cubic endpoint's unstable Gamma-mode
subspace, while energy and homogeneous volume remain separate path variables.
The inputs are deterministic single-path first-principles data, not a sampled
population; no error bars or statistical test are appropriate.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Editable vector text is intentional: the SVG is a manuscript working asset.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.linewidth": 0.8,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "legend.frameon": False,
    }
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUTS = ROOT / "outputs" / "batio3_t_to_c_pbe100_dzp10au"
DEFAULT_FIGURE = ROOT / "papar" / "VARNEB_CPC" / "bto_gamma_mode_path"

COLORS = {
    "soft": "#7C6CCF",
    "stable_1": "#33B5A5",
    "stable_2": "#767676",
    "energy": "#0F4D92",
    "volume": "#B64342",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--modes",
        type=Path,
        default=DEFAULT_OUTPUTS / "bto_tetragonal_to_cubic_n7_gamma_modes_phonopy_direct.json",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=DEFAULT_OUTPUTS / "bto_tetragonal_to_cubic_n7_static_audit_summary.json",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument(
        "--source-data",
        type=Path,
        default=ROOT / "papar" / "VARNEB_CPC" / "figures" / "bto_gamma_mode_path_source_data.csv",
    )
    return parser.parse_args()


def _label(group: dict) -> str:
    frequency = float(group["frequency_cm1"])
    if frequency < -1e-3:
        return rf"soft $\Gamma$ ({frequency:.0f} cm$^{{-1}}$)"
    return rf"stable $\Gamma$ ({frequency:.0f} cm$^{{-1}}$)"


def extract_source_data(mode_payload: dict, summary_payload: dict) -> tuple[list[dict], list[tuple[str, np.ndarray, str]]]:
    coordinates = np.asarray(mode_payload["normal_coordinates_sqrt_amu_A"], dtype=float)
    n_images = int(mode_payload["n_images"])
    if coordinates.shape[0] != n_images:
        raise ValueError("mode-coordinate image count does not match the report")
    diagnostics = summary_payload["path_diagnostics"]["images"]
    enthalpies = np.asarray(summary_payload["image_enthalpies_eV"], dtype=float)
    if len(diagnostics) != n_images or len(enthalpies) != n_images:
        raise ValueError("VCNEB summary image count does not match the mode report")
    curves: list[tuple[str, np.ndarray, str]] = []
    for group in mode_payload["degenerate_mode_subspaces"]:
        if group["translation_subspace"]:
            continue
        indices = np.asarray(group["mode_indices"], dtype=int)
        amplitude = np.linalg.norm(coordinates[:, indices], axis=1)
        if float(np.max(amplitude)) < 1e-10:
            continue
        frequency = float(group["frequency_cm1"])
        color = COLORS["soft"] if frequency < -1e-3 else (
            COLORS["stable_1"] if abs(frequency) < 200 else COLORS["stable_2"]
        )
        curves.append((_label(group), amplitude, color))
    curves.sort(key=lambda item: float(np.max(item[1])), reverse=True)
    rows = []
    initial_volume = float(diagnostics[0]["volume_A3"])
    for image_index in range(n_images):
        row = {
            "image_index": image_index,
            "relative_enthalpy_meV_per_formula_unit": 1000.0 * (enthalpies[image_index] - enthalpies[0]),
            "volume_A3": float(diagnostics[image_index]["volume_A3"]),
            "relative_volume_percent": 100.0 * (float(diagnostics[image_index]["volume_A3"]) / initial_volume - 1.0),
        }
        for label, values, _ in curves:
            row[f"{label}_norm_sqrt_amu_A"] = float(values[image_index])
        rows.append(row)
    return rows, curves


def write_source_data(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def add_panel_label(axis, label: str) -> None:
    axis.text(-0.17, 1.03, label, transform=axis.transAxes, fontweight="bold", fontsize=8)


def plot(rows: list[dict], curves: list[tuple[str, np.ndarray, str]], output: Path) -> list[Path]:
    image = np.asarray([row["image_index"] for row in rows])
    energy = np.asarray([row["relative_enthalpy_meV_per_formula_unit"] for row in rows])
    volume = np.asarray([row["relative_volume_percent"] for row in rows])
    figure = plt.figure(figsize=(7.2, 2.55))  # double-column width, 183 mm
    grid = figure.add_gridspec(2, 3, width_ratios=[1.15, 1.15, 0.9], wspace=0.55, hspace=0.6)
    hero = figure.add_subplot(grid[:, :2])
    energy_axis = figure.add_subplot(grid[0, 2])
    volume_axis = figure.add_subplot(grid[1, 2])

    for label, values, color in curves:
        hero.plot(image, values, "o-", color=color, lw=1.6, ms=4, label=label)
    hero.set_xlabel(r"VCNEB image index (T $\rightarrow$ C)")
    hero.set_ylabel(r"subspace coordinate norm ($\sqrt{\mathrm{amu}}$ \AA)")
    hero.set_xticks(image)
    hero.legend(loc="upper right", fontsize=6, handlelength=1.7)
    hero.text(
        0.02,
        0.04,
        r"cubic endpoint $\Gamma$ eigenvectors" "\nrigid translations removed",
        transform=hero.transAxes,
        fontsize=5.8,
        va="bottom",
        color="#4D4D4D",
    )
    add_panel_label(hero, "a")

    energy_axis.plot(image, energy, "o-", color=COLORS["energy"], lw=1.5, ms=3.5)
    energy_axis.set_ylabel("relative enthalpy\n(meV/f.u.)")
    energy_axis.set_xticks(image)
    energy_axis.set_xticklabels([])
    energy_axis.set_title("monotonic path", fontsize=7, pad=2)
    add_panel_label(energy_axis, "b")

    volume_axis.plot(image, volume, "o-", color=COLORS["volume"], lw=1.5, ms=3.5)
    volume_axis.set_xlabel("")
    volume_axis.set_ylabel(r"$\Delta V/V_0$ (%)")
    volume_axis.set_xticks(image)
    add_panel_label(volume_axis, "c")

    figure.text(
        0.995,
        0.005,
        "Single deterministic ABACUS/PBE path; no error bars.",
        ha="right",
        va="bottom",
        fontsize=5.5,
        color="#4D4D4D",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    saved = []
    for suffix, dpi in ((".svg", 600), (".pdf", 600), (".png", 600)):
        path = output.with_suffix(suffix)
        figure.savefig(path, dpi=dpi, bbox_inches="tight")
        if suffix == ".svg":
            # Matplotlib's wrapped SVG path data contains cosmetic trailing
            # spaces. Remove them so the editable vector asset passes the
            # repository whitespace gate without changing rendering.
            path.write_text(
                "\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines()) + "\n",
                encoding="utf-8",
            )
        saved.append(path)
    plt.close(figure)
    return saved


def main() -> None:
    args = parse_args()
    modes = json.loads(args.modes.read_text(encoding="utf-8"))
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    rows, curves = extract_source_data(modes, summary)
    write_source_data(rows, args.source_data)
    saved = plot(rows, curves, args.output)
    print("[DONE] wrote " + ", ".join(str(path) for path in saved))
    print(f"[DONE] source data={args.source_data}")


if __name__ == "__main__":
    main()
