"""Rebuild the four archived BTO/GaN landscape figures without running DFT.

Each renderer reads committed audit data. The reconstructed source CSV must
match the manuscript version byte for byte; PNG identity is optional because
fonts and Matplotlib versions can change raster output across environments.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
FIGURES = (
    (
        "bto_restricted27", "scripts.plot_bto_qy0_restricted_sheet27",
        ("--report", "benchmarks/numerical_integrity/bto_qy0_27_node_sheet_20260929.json"),
        "bto_qy0_restricted_sheet_27_20260929", 27,
    ),
    (
        "bto_frozen289", "scripts.plot_bto_transverse_soft_landscape",
        ("--assembled", "outputs/batio3_t_to_c_pbe100_dzp10au/bto_transverse_soft_frozen25_assembled_job27781068.json",
         "--path-audit", "outputs/batio3_t_to_c_pbe100_dzp10au/bto_transverse_soft_plane_offline_audit_2026-09-26.json",
         "--path-report", "outputs/batio3_t_to_c_pbe100_dzp10au/bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json",
         "--analysis41", "outputs/batio3_t_to_c_pbe100_dzp10au/bto_transverse_soft_frozen41_center_holdout_score_job27781133.json",
         "--analysis59", "outputs/batio3_t_to_c_pbe100_dzp10au/bto_transverse_soft_frozen59_edge_holdout_score_job27781166.json",
         "--analysis81", "benchmarks/numerical_integrity/bto_transverse_soft_frozen81_20260928.json",
         "--analysis289", "benchmarks/numerical_integrity/bto_transverse_soft_frozen289_20260929.json",
         "--allow-exploratory-contours"),
        "bto_frozen_soft_mode_289_20260929", 289,
    ),
    (
        "gan_joint289", "scripts.plot_gan_600eV_local_joint_cut17",
        ("--pilot-audit", "benchmarks/numerical_integrity/gan_600eV_ts_2d_pilot_20260928.json",
         "--refinement5-audit", "benchmarks/numerical_integrity/gan_600eV_ts_2d_5x5_refinement_20260928.json",
         "--dense9-audit", "benchmarks/numerical_integrity/gan_600eV_ts_2d_9x9_refinement_20260928.json",
         "--dense17-audit", "benchmarks/numerical_integrity/gan_600eV_ts_2d_17x17_refinement_20260930.json"),
        "gan_600eV_local_joint_dft289_20260930_v2", 289,
    ),
    (
        "gan_transverse162", "scripts.plot_gan_600eV_atomic_transverse_q9",
        ("--report", "benchmarks/numerical_integrity/gan_600eV_atomic_tube_q9_20260929.json",
         "--refinement", "benchmarks/numerical_integrity/gan_600eV_atomic_tube_refinement_20260928.json",
         "--trajectory", "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj"),
        "gan_600eV_atomic_transverse_162_20260929_v2", 162,
    ),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rebuild(output_dir: Path, strict_png: bool = False) -> list[dict]:
    if output_dir.exists():
        raise FileExistsError(f"output directory already exists: {output_dir}")
    output_dir.mkdir(parents=True)
    records = []
    for name, module, options, archive_name, expected_nodes in FIGURES:
        prefix = output_dir / name
        command = [sys.executable, "-m", module, *options,
                   "--output-prefix", str(prefix)]
        subprocess.run(command, cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        archived = ROOT / "paper" / "VARNEB_CPC" / "figures" / archive_name
        csv_output = prefix.with_name(prefix.name + "_source_data.csv")
        csv_archived = archived.with_name(archive_name + "_source_data.csv")
        if sha256(csv_output) != sha256(csv_archived):
            raise ValueError(f"{name}: rebuilt source CSV differs from archived manuscript data")
        qa = json.loads(prefix.with_name(prefix.name + "_qa.json").read_text(encoding="utf-8"))
        node_count = (qa.get("n_measured_DFT_nodes") or qa.get("n_audited_DFT_points")
                      or qa.get("n_measured_DFT_points")
                      or qa.get("interpolation", {}).get("n_measured_DFT_samples"))
        if node_count != expected_nodes:
            raise ValueError(f"{name}: rebuilt QA node count is not {expected_nodes}")
        same_png = sha256(prefix.with_suffix(".png")) == sha256(archived.with_suffix(".png"))
        if strict_png and not same_png:
            raise ValueError(f"{name}: rebuilt PNG differs in this environment")
        records.append({"figure": name, "DFT_nodes": expected_nodes,
                        "source_csv_byte_identical": True,
                        "png_byte_identical": same_png})
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="fresh directory for regenerated PDF/SVG/PNG/CSV/QA files")
    parser.add_argument("--strict-png", action="store_true",
                        help="also require exact PNG bytes (environment-sensitive)")
    args = parser.parse_args()
    print(json.dumps(rebuild(args.output_dir, args.strict_png), indent=2))


if __name__ == "__main__":
    main()
