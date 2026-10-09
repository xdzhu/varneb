"""Plot six actual G1 reconstruction statics, never an interpolated MEP.

The original nine-image observation and all new raw SCF logs are replayed.
Physical template bytes were audited on HF; their recorded hashes must remain
the six fixed contract values. No new DFT or stationary-TS claim is made here.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import tempfile

import numpy as np
from ase.io import read
from ase.calculators.singlepoint import SinglePointCalculator

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.prepare_hfo2_sampling_bridge import SOURCE_OBSERVATION_SHA256, load_manifest
from scripts.prepare_hfo2_channel_work_probes import write_structure
from vcneb import VCNEB


def fixed_writer_geometry_matches(path, atoms):
    """Replay the case's actual fixed writer without requiring ASE-ABACUS I/O.

    HF already checks ASE-ABACUS round-trip equivalence before each SCF.
    Offline, only floating textual tails may differ between ASE versions;
    species/basis/order, grammar, units and every other token are retained.
    """
    with tempfile.TemporaryDirectory(prefix="hfo2-writer-replay-") as folder:
        expected = Path(folder) / "STRU"
        write_structure(expected, atoms)
        a, b = path.read_text().split(), expected.read_text().split()
    if len(a) != len(b):
        return False
    for left, right in zip(a, b):
        if left == right:
            continue
        try:
            x, y = float(left), float(right)
        except ValueError:
            return False
        if not np.isfinite(x) or not np.isfinite(y) or abs(x-y) > 1e-13:
            return False
    return True


def build_data(observation, completed):
    images, source, digest = read_evaluated_observation(observation)
    manifest = load_manifest(completed)
    summary = json.loads((completed / "summary.json").read_text())
    if (digest != SOURCE_OBSERVATION_SHA256 or manifest["source_observation_sha256"] != digest
            or source["snapshot_step"] != 69 or source["source_job_id"] != "28319570"
            or len(images) != 9 or manifest["physical_contract"] != CONTRACT
            or manifest["formula_units"] != 4 or manifest["pressure_GPa"] != 0
            or manifest["climb"] or manifest["ordinary_fmax_eV_A"] != .10
            or summary["status"] != "six_sampling_statics_audited_not_relaxed_MEP"
            or summary["new_SCF_calls"] != 6 or len(summary["points"]) != 6
            or summary["automatic_restart_or_G2_submission"]):
        raise ValueError("the complete frozen six-point G1 audit is required")
    band = VCNEB(images, pressure=0., climb=False, k=.2)
    # The fixed production metric, not an image-wise Cartesian MIC distance.
    arc = band.reaction_coordinate()
    arc /= arc[-1]
    energy0 = images[0].get_potential_energy()
    if abs(energy0-manifest["common_PO_energy_eV_cell"]) > 1e-10:
        raise ValueError("common PO reference changed")
    rows = [{"kind": "original_band_sample", "source_image": i, "segment": "",
             "fraction": None, "s_normalized": float(arc[i]),
             "energy_relative_PO_meV_fu": float((a.get_potential_energy()-energy0)*250),
             "raw_log_sha256": source["raw_image_evaluations"][i]["raw_log_sha256"]}
            for i, a in enumerate(images)]
    new = []
    added_atoms = {}
    for i, (point, record) in enumerate(zip(manifest["points"], summary["points"])):
        rawdir = completed / "calculations" / f"{i:02d}" / "scf_000000"
        actual = audited_results(rawdir)
        input_hashes = json.loads((rawdir / "input_sha256.json").read_text())
        if (point["index"] != i or record["index"] != i
                or record["segment"] != point["segment"] or record["fraction"] != point["fraction"]
                or record["manifest_sha256"] != sha256(completed / "manifest.json")
                or record["physical_input_sha256"] != CONTRACT
                or {n: input_hashes[n] for n in CONTRACT} != CONTRACT
                or sha256(rawdir / "STRU") != record["input_STRU_sha256"]
                or sha256(rawdir / "STRU") != input_hashes["STRU"]
                or sha256(rawdir / "OUT.ABACUS/running_scf.log") != record["raw_log_sha256"]
                or sha256(completed / point["geometry"]) != record["geometry_sha256"]
                or record["geometry_sha256"] != point["geometry_sha256"]
                or not fixed_writer_geometry_matches(rawdir / "STRU",
                                                     read(completed / point["geometry"], format="vasp"))
                or abs(actual["energy"]-record["results"]["energy"]) > 1e-10
                or not np.allclose(actual["forces"], record["results"]["forces"], atol=1e-12, rtol=0)
                or not np.allclose(actual["stress"], record["results"]["stress"], atol=1e-12, rtol=0)):
            raise ValueError("actual raw static provenance changed")
        energy = float((actual["energy"]-energy0)*250)
        if abs(energy-record["relative_energy_meV_fu"]) > 1e-10:
            raise ValueError("static energy conversion changed")
        left, right = point["segment"]
        fraction = point["fraction"]
        new.append({"kind": "new_linear_reconstruction_static", "source_image": None,
                    "segment": f"{left}->{right}", "fraction": fraction,
                    "s_normalized": float((1-fraction)*arc[left]+fraction*arc[right]),
                    "energy_relative_PO_meV_fu": energy, "raw_log_sha256": record["raw_log_sha256"]})
        atoms = read(completed / point["geometry"], format="vasp")
        atoms.calc = SinglePointCalculator(atoms, **actual)
        added_atoms.setdefault(left, []).append((fraction, atoms))
    old_max = max(r["energy_relative_PO_meV_fu"] for r in rows)
    new_max = max(r["energy_relative_PO_meV_fu"] for r in new)
    if (abs(old_max-summary["old_sampled_band_maximum_meV_fu"]) > 1e-10
            or abs(new_max-summary["highest_new_sample_meV_fu"]) > 1e-10
            or abs(new_max-old_max-summary["difference_meV_fu"]) > 1e-10):
        raise ValueError("summary sampled maxima changed")
    # Audit the added points' residuals without optimizing or making new SCFs.
    refined = []
    for i, atoms in enumerate(images):
        refined.append(atoms)
        refined.extend(a for _, a in sorted(added_atoms.get(i, []), key=lambda x: x[0]))
    reconstructed = VCNEB(refined, pressure=0., climb=False, k=.2)
    residual = float(np.linalg.norm(reconstructed.get_forces().reshape(-1, 3), axis=1).max())
    return {"rows": rows+new, "old_max_meV_fu": old_max, "new_max_meV_fu": new_max,
            "difference_meV_fu": new_max-old_max, "source_observation_sha256": digest,
            "new_manifest_sha256": sha256(completed / "manifest.json"),
            "summary_sha256": sha256(completed / "summary.json"),
            "reconstructed_band": {"total_images": len(refined), "moving_images": len(refined)-2,
                                   "optimizer_steps": 0, "new_SCF_calls": 0,
                                   "replayed_ordinary_fmax_eV_A": residual,
                                   "ordinary_target_eV_A": .10,
                                   "ordinary_residual_passed": residual <= .10,
                                   "status": "linear inserted cached band, not an optimized new MEP"},
            "new_static_SCF_count": 6, "new_DFT_calls_for_plot": 0,
            "stationary_TS_or_relaxed_MEP_certified": False,
            "interpretation": "actual samples of the stated linear reconstruction, not a relaxed MEP or blind prediction"}


def make_figure(data):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 10, "axes.labelsize": 11, "legend.fontsize": 9.3,
                         "axes.linewidth": .9, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(183/25.4, 99/25.4), sharey=True)
    fig.subplots_adjust(left=.12, right=.98, bottom=.17, top=.75, wspace=.20)
    old = [r for r in data["rows"] if r["kind"] == "original_band_sample"]
    new = [r for r in data["rows"] if r["kind"] != "original_band_sample"]
    axes[0].plot([r["s_normalized"] for r in old], [r["energy_relative_PO_meV_fu"] for r in old],
                 "o-", color="#75818A", markersize=4, linewidth=1.2, label="Original 9-image samples")
    axes[0].scatter([r["s_normalized"] for r in new], [r["energy_relative_PO_meV_fu"] for r in new],
                    marker="D", color="#C97738", s=24, zorder=3, label="6 new static checks")
    axes[0].legend(loc="lower left", bbox_to_anchor=(-.04, 1.09), frameon=True,
                   facecolor="white", edgecolor="#9CA3AA", framealpha=.85, borderpad=.4)
    for segment, endpoints, color, marker in (("2->3", (2, 3), "#256B91", "o"),
                                             ("5->6", (5, 6), "#C97738", "^")):
        points = [(0., old[endpoints[0]]["energy_relative_PO_meV_fu"]),
                  *[(r["fraction"], r["energy_relative_PO_meV_fu"]) for r in new if r["segment"] == segment],
                  (1., old[endpoints[1]]["energy_relative_PO_meV_fu"])]
        axes[1].plot(*np.asarray(points).T, marker=marker, color=color, linewidth=1.3,
                     markersize=4, label=f"Segment {segment.replace('->', chr(8594))}")
    axes[1].legend(loc="lower left", bbox_to_anchor=(-.04, 1.09), frameon=True,
                   facecolor="white", edgecolor="#9CA3AA", framealpha=.85, borderpad=.4)
    axes[0].set_xlabel("Normalized original path arc")
    axes[1].set_xlabel(r"Segment fraction $\lambda$")
    axes[0].set_ylabel(r"$E-E_{\rm PO^+}$ (meV/f.u.)")
    for i, axis in enumerate(axes):
        axis.axhline(data["old_max_meV_fu"], color="#A6AEB5", linestyle=":", linewidth=1., zorder=0)
        axis.set_xlim(-.025, 1.025)
        axis.set_ylim(-22, max(46, data["new_max_meV_fu"]+7))
        axis.set_xticks([0, .25, .5, .75, 1.])
        axis.tick_params(direction="in", top=True, right=True)
        for spine in axis.spines.values():
            spine.set_visible(True)
        axis.text(-.16, 1.03, f"({chr(97+i)})", transform=axis.transAxes,
                  fontsize=13, fontweight="normal", va="bottom")
    fig.canvas.draw()
    return fig


def export(data, output):
    if output.exists():
        raise FileExistsError("refusing an existing sampling figure")
    fig = make_figure(data)
    output.mkdir(parents=True)
    stem = output / "hfo2_G1_sampling_bridge"
    for ext in ("svg", "pdf", "png", "tiff"):
        kwargs = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        fig.savefig(stem.with_suffix("."+ext), dpi=600 if ext == "tiff" else 300,
                    facecolor="white", **kwargs)
    csv_path = output / "source_data.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data["rows"][0]))
        writer.writeheader()
        writer.writerows(data["rows"])
    qa = {k: v for k, v in data.items() if k != "rows"}
    qa.update({"plot_script_sha256": sha256(Path(__file__)), "source_data_sha256": sha256(csv_path),
               "panel_labels": ["(a)", "(b)"], "panel_label_weight": "normal",
               "panel_titles_empty": all(not a.get_title() for a in fig.axes),
               "all_spines_visible": all(s.get_visible() for a in fig.axes for s in a.spines.values()),
               "panel_positions_normalized": [a.get_position().bounds for a in fig.axes],
               "figure_dimensions_mm": (fig.get_size_inches()*25.4).tolist(),
               "line_interpretation": "straight sample-joining guides, not smooth MEP interpolation",
               "statistics": "9 reused band samples and 6 new SCFs; no independent replicates or error bars",
               "visual_review": "pending actual pixel review",
               "file_sha256": {p.name: sha256(p) for p in sorted(output.glob(stem.name+".*"))}})
    (output / "qa.json").write_text(json.dumps(qa, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return qa


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--observation", type=Path, required=True)
    p.add_argument("--completed", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing an existing sampling figure")
    qa = export(build_data(args.observation, args.completed), args.output)
    print(json.dumps({"new_static_SCF_count": 6, "old_max_meV_fu": qa["old_max_meV_fu"],
                      "new_max_meV_fu": qa["new_max_meV_fu"], "new_DFT_calls_for_plot": 0}))


if __name__ == "__main__":
    main()
