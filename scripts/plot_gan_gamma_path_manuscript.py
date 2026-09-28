"""Plot the all-image GaN endpoint-Gamma decomposition for the CPC evidence set.

Figure contract: along one audited 600-eV B4-to-B1 VCNEB chain, three
endpoint optical subspaces describe the atomic displacement, while large cell
strain remains a distinct coordinate.  The separate reconstruction panel
shows what those three groups leave unexplained.  This is a quantitative-grid
figure, not an energy decomposition or a certified TS-mode analysis.
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

from scripts.analyze_gan_gamma_path_subspaces import audit


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODES = ROOT / "outputs/gan_b4_b1_gamma_1x1x1_20260926"
DEFAULT_PATH = ROOT / "benchmarks/numerical_integrity/gan_b4_b1_tetragonal_empirical_atomic_strain_20260926.csv"
DEFAULT_CHAIN = ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj"
DEFAULT_PREFIX = ROOT / "paper/VARNEB_CPC/figures/gan_gamma_path_600eV"
NAVY = "#284A68"
TEAL = "#258579"
OCHRE = "#BC7844"
GREY = "#566777"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.size": 10.0, "axes.labelsize": 10.7,
    "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "axes.linewidth": 0.9, "axes.spines.top": True,
    "axes.spines.right": True,
})


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assemble(directory: Path, path_csv: Path, chain: Path) -> tuple[dict, list[dict]]:
    report, records = audit(directory, path_csv)
    if report["status"] != "audited_GaN_endpoint_Gamma_atomic_subspaces_not_TS_modes":
        raise ValueError("endpoint Gamma subspace audit status changed")
    archived = json.loads((directory / "subspace_audit_v2.json").read_text(encoding="utf-8"))
    if archived != report:
        raise ValueError("endpoint Gamma subspace reconstruction differs from archived audit")
    original_chain = directory / "gan_vasp_tetragonal_final_29_images.traj"
    if sha256(chain) != sha256(original_chain):
        raise ValueError("Gamma projection trajectory is not the frozen 600-eV final chain")
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    dft_audit = json.loads((directory / "audit.json").read_text(encoding="utf-8"))
    expected_electronic = {
        "INCAR": "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772",
        "KPOINTS": "b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee",
        "POTCAR": "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168",
    }
    if (manifest.get("pressure_GPa") != 45.7
            or manifest.get("supercell_matrix") != np.eye(3, dtype=int).tolist()
            or dft_audit.get("source_manifest_sha256") != sha256(directory / "manifest.json")
            or any(sha256(directory / f"{phase}_CONTCAR")
                   != manifest["source"][phase]["sha256"]["CONTCAR"]
                   for phase in ("B4", "B1"))
            or any(manifest["source"][phase]["sha256"][key] != digest
                   for phase in ("B4", "B1")
                   for key, digest in expected_electronic.items())):
        raise ValueError("Gamma calculation differs from the VASP 600-eV production contract")
    phases = {phase: [row for row in records if row["reference_phase"] == phase]
              for phase in ("B4", "B1")}
    if any(len(phases[phase]) != 29 for phase in phases):
        raise ValueError("Gamma mode path must contain 29 images in each reference")
    joined = []
    for index, (b4, b1) in enumerate(zip(phases["B4"], phases["B1"])):
        if (b4["image_index"] != index or b1["image_index"] != index
                or abs(b4["reaction_coordinate_normalized"] - b1["reaction_coordinate_normalized"]) > 1e-12
                or abs(b4["relative_enthalpy_eV_per_GaN"] - b1["relative_enthalpy_eV_per_GaN"]) > 1e-12):
            raise ValueError("B4/B1 Gamma projections are not paired on the same 29-image chain")
        joined.append({
            "image_index": index,
            "path_s": b4["reaction_coordinate_normalized"],
            "relative_enthalpy_eV_per_GaN": b4["relative_enthalpy_eV_per_GaN"],
            **{f"B4_group{rank}_Q_sqrt_amu_A": b4[f"group{rank}_Q_norm_d0p01_sqrt_amu_A"]
               for rank in (1, 2, 3)},
            **{f"B1_group{rank}_Q_sqrt_amu_A": b1[f"group{rank}_Q_norm_d0p01_sqrt_amu_A"]
               for rank in (1, 2, 3)},
            "B4_three_group_residual_sqrt_amu_A": b4["three_group_atomic_residual_sqrt_amu_A"],
            "B1_three_group_residual_sqrt_amu_A": b1["three_group_atomic_residual_sqrt_amu_A"],
            **{f"B4_symmetric_strain_{axis}": b4[f"symmetric_strain_{axis}"]
               for axis in ("xx", "yy", "zz")},
        })
    if not np.all(np.diff([row["path_s"] for row in joined]) > 0):
        raise ValueError("path arc coordinate is not strictly increasing")
    return report, joined


def plot(report: dict, rows: list[dict], prefix: Path, sources: dict[str, Path]) -> dict:
    exports = [prefix.with_suffix(ext) for ext in (".svg", ".pdf", ".png", ".tiff")]
    data_path = prefix.with_name(prefix.name + "_source_data.csv")
    qa_path = prefix.with_name(prefix.name + "_qa.json")
    if any(path.exists() for path in (*exports, data_path, qa_path)):
        raise FileExistsError("refusing to overwrite a GaN manuscript figure or source data")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    s = np.asarray([row["path_s"] for row in rows], dtype=float)
    peak = int(np.argmax([row["relative_enthalpy_eV_per_GaN"] for row in rows]))
    if peak != 15:
        raise ValueError("600-eV final chain peak is no longer image 15")
    fig, axes = plt.subplots(2, 2, figsize=(7.22, 5.02), sharex=True)
    for panel, (phase, ax) in enumerate((("B4", axes[0, 0]), ("B1", axes[0, 1]))):
        for rank, color in enumerate((NAVY, TEAL, OCHRE), start=1):
            group = report["phases"][phase]["groups_ranked_by_path_amplitude"][rank - 1]
            indices = ",".join(str(item + 1) for item in group["mode_indices"])
            q = [row[f"{phase}_group{rank}_Q_sqrt_amu_A"] for row in rows]
            ax.plot(s, q, color=color, lw=1.75, label=rf"$\Gamma_{{{indices}}}$, {group['frequency_cm1_at_0p01']:.0f} cm$^{{-1}}$")
            ax.scatter(s[::4], np.asarray(q)[::4], s=10, color=color, zorder=3)
        ax.set_ylim(0, 4.7)
        ax.set_ylabel(fr"{phase} $|Q_\Gamma|$ ($\sqrt{{\mathrm{{amu}}}}\,\AA$)")
        ax.legend(loc="upper left" if panel == 0 else "upper right", fontsize=8.4,
                  frameon=True, facecolor="white", edgecolor="#A9B5BE", framealpha=0.83)
    strain = axes[1, 0]
    for axis, color, style in (("xx", NAVY, "-"), ("yy", TEAL, "--"), ("zz", OCHRE, "-.")):
        strain.plot(s, 100 * np.asarray([row[f"B4_symmetric_strain_{axis}"] for row in rows]),
                    color=color, lw=1.8, ls=style, label=rf"$\eta_{{{axis}}}$")
    strain.axhline(0, color="#A8B2BB", lw=0.7)
    strain.set_ylabel("Strain relative to B4 (%)")
    strain.legend(loc="lower left", bbox_to_anchor=(0.005, 1.025), ncol=3,
                  fontsize=8.5, columnspacing=0.9,
                  frameon=True, facecolor="white", edgecolor="#A9B5BE", framealpha=0.83)
    residual = axes[1, 1]
    for phase, color, style in (("B4", NAVY, "-"), ("B1", TEAL, "--")):
        value = 1000 * np.asarray([row[f"{phase}_three_group_residual_sqrt_amu_A"] for row in rows])
        residual.plot(s, value, color=color, ls=style, lw=1.75, label=phase)
    residual.set_ylim(bottom=0)
    residual.set_ylabel(r"Three-group residual ($10^{-3}\sqrt{\mathrm{amu}}\,\AA$)")
    residual.legend(loc="upper right", fontsize=8.8, frameon=True,
                    facecolor="white", edgecolor="#A9B5BE", framealpha=0.83)
    for index, ax in enumerate(axes.flat):
        ax.axvline(s[peak], color=GREY, lw=0.9, ls=":", zorder=0)
        ax.set_xlim(0, 1)
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.tick_params(direction="in", top=True, right=True, length=3.3)
        ax.text(-0.105, 1.035, f"({chr(97 + index)})", transform=ax.transAxes,
                fontsize=13, fontweight="normal", ha="left", va="bottom")
    for ax in axes[1]:
        ax.set_xlabel("Normalized VCNEB arc coordinate, $s$")
    fig.subplots_adjust(left=0.13, right=0.975, bottom=0.13, top=0.925,
                        wspace=0.31, hspace=0.32)
    for ext in (".svg", ".pdf"):
        fig.savefig(prefix.with_suffix(ext))
    svg_path = prefix.with_suffix(".svg")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines()) + "\n",
        encoding="utf-8",
    )
    fig.savefig(prefix.with_suffix(".png"), dpi=300)
    fig.savefig(prefix.with_suffix(".tiff"), dpi=600)
    plt.close(fig)
    with data_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    qa = {
        "status": "same_600eV_all_image_endpoint_Gamma_path_manuscript_figure",
        "core_conclusion": "Endpoint Gamma optical subspaces track atomic reconstruction while the GaN VCNEB path undergoes separate large cell strain.",
        "archetype": "quantitative_grid",
        "backend": "Python/Matplotlib",
        "panels": {"a": "B4-reference group amplitudes", "b": "B1-reference group amplitudes",
                   "c": "symmetric strain versus B4", "d": "three-group atomic reconstruction residual"},
        "n_images_total": 29,
        "pressure_GPa": 45.7,
        "highest_image_index_not_certified_TS": peak,
        "maximum_three_group_residual_sqrt_amu_A": {
            phase: report["phases"][phase]["max_three_group_atomic_reconstruction_residual_sqrt_amu_A"]
            for phase in ("B4", "B1")
        },
        "maximum_absolute_diagonal_strain": {
            axis: max(abs(row[f"B4_symmetric_strain_{axis}"]) for row in rows)
            for axis in ("xx", "yy", "zz")
        },
        "source_sha256": {name: sha256(path) for name, path in sources.items()},
        "source_data_sha256": sha256(data_path),
        "export_sha256": {path.name: sha256(path) for path in exports},
        "limitations": [
            "Gamma modes are fixed-cell endpoint bases, not joint atom-strain TS modes.",
            "The plotted amplitudes are geometry projections, not energy or barrier shares.",
            "The endpoint bases describe a 0-K path; they do not establish finite-temperature mechanisms.",
            "The dotted line marks the highest discrete image, not a certified stationary saddle.",
        ],
    }
    qa_path.write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    return qa


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modes", type=Path, default=DEFAULT_MODES)
    parser.add_argument("--path-csv", type=Path, default=DEFAULT_PATH)
    parser.add_argument("--chain", type=Path, default=DEFAULT_CHAIN)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_PREFIX)
    args = parser.parse_args()
    report, rows = assemble(args.modes, args.path_csv, args.chain)
    qa = plot(report, rows, args.output_prefix, {
        "endpoint_Gamma_audit": args.modes / "audit.json",
        "endpoint_Gamma_subspace_audit": args.modes / "subspace_audit_v2.json",
        "endpoint_Gamma_manifest": args.modes / "manifest.json",
        "600eV_final_chain": args.chain,
        "path_geometry_source": args.path_csv,
        "plotter": Path(__file__),
    })
    print(json.dumps({"status": qa["status"], "n_images": qa["n_images_total"],
                      "residual": qa["maximum_three_group_residual_sqrt_amu_A"]}))


if __name__ == "__main__":
    main()
