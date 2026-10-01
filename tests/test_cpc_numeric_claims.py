"""Keep headline CPC numbers tied to committed figure and launch sources."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

from ase.units import Ry
from scripts.audit_bto_conditional_five_point_patch import audit_patch


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "VARNEB_CPC"
FIGURES = PAPER / "figures"
MANUSCRIPT = re.sub(
    r"\s+", " ", (PAPER / "varneb_CPC.tex").read_text(encoding="utf-8")
)


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_five_gan_barriers_and_force_gate_match_manuscript() -> None:
    rows = _csv_rows(FIGURES / "gan_multibackend_validation_source_data.csv")
    order = ("ABACUS", "VASP", "QE", "ABINIT", "CP2K")
    barriers = []
    for backend in order:
        curve = sorted(
            (row for row in rows
             if row["backend"] == backend and row["status"] == "converged"),
            key=lambda row: int(row["image_index"]),
        )
        assert len(curve) == 29
        assert [int(row["image_index"]) for row in curve] == list(range(29))
        assert float(curve[0]["final_fmax_eV_per_A"]) <= 0.10
        peak = max(curve, key=lambda row: float(row["relative_enthalpy_eV_per_GaN"]))
        assert int(peak["image_index"]) == 15
        forward = float(peak["relative_enthalpy_eV_per_GaN"])
        assert abs(forward - float(curve[0]["forward_barrier_eV_per_GaN"])) < 1e-9
        barriers.append(f"{forward:.4f}")

    reported = ", ".join(barriers[:-1]) + ", and " + barriers[-1] + " eV/GaN"
    assert reported in MANUSCRIPT
    assert MANUSCRIPT.count(reported) == 2  # abstract and GaN results


def test_gan_plot_source_matches_archived_chain_evidence_and_cp2k_caveat() -> None:
    rows = _csv_rows(FIGURES / "gan_multibackend_validation_source_data.csv")
    evidence = json.loads(
        (PAPER / "evidence/gan_45p7_multibackend_vcneb_20260924.json")
        .read_text(encoding="utf-8")
    )
    assert evidence["contract"]["n_images_total"] == 29
    assert evidence["contract"]["pressure_gpa"] == 45.7
    for backend, source in evidence["backends"].items():
        curve = sorted(
            (row for row in rows if row["backend"].lower() == backend),
            key=lambda row: int(row["image_index"]),
        )
        assert len(curve) == 29
        assert source["status"] == "converged"
        assert source["highest_image_index"] == 15
        for row in curve:
            assert abs(
                float(row["forward_barrier_eV_per_GaN"])
                - source["barrier_eV_per_GaN"]
            ) < 1e-9
            assert abs(
                float(row["final_fmax_eV_per_A"])
                - source["final_max_generalized_force_eV_per_A"]
            ) < 1e-9
        archived_profile = source.get("relative_enthalpy_eV_per_cell")
        if archived_profile is not None:
            assert len(archived_profile) == 29
            for row, per_cell in zip(curve, archived_profile):
                assert abs(
                    float(row["relative_enthalpy_eV_per_GaN"]) - per_cell / 2
                ) < 1e-9

    cp2k_audit = json.loads(
        (ROOT / "benchmarks/numerical_integrity/"
         "gan_cp2k_final_chain_raw_energy_stress_20260928.json")
        .read_text(encoding="utf-8")
    )
    assert cp2k_audit["n_interior_images"] == 27
    assert cp2k_audit["n_raw_last_scf_converged"] == 27
    assert cp2k_audit["n_raw_last_program_end"] == 0
    assert cp2k_audit["max_absolute_cache_chain_difference"] == {
        "energy_ev": 0, "force_ev_per_angstrom": 0,
        "stress_ev_per_angstrom3": 0,
    }
    assert "original-text atomic forces and final run-end markers are unavailable" in MANUSCRIPT


def test_gan_calculator_table_matches_generated_peak_inputs() -> None:
    inputs = PAPER / "evidence/gan_45p7_final_peak_inputs_20260927"
    source_index = (
        PAPER / "evidence/gan_45p7_remote_source_index_20260927.md"
    ).read_text(encoding="utf-8")
    expected_hashes = {
        "qe.pwi": "b558fff1f3b16379a009f905514fa92c141062f24bb2859d1b323e57b190c767",
        "abinit.in": "b4d4388996655c4d13aae3c786dd3bd15f01c81ab13a23198a2d74e13b2069f3",
        "cp2k.inp": "935170d769aa368ddc8f1bc3089223953b275ea5b52bc70168b84285ab93cff1",
    }
    for filename, expected_hash in expected_hashes.items():
        assert hashlib.sha256((inputs / filename).read_bytes()).hexdigest() == expected_hash
        assert expected_hash in source_index

    qe = (inputs / "qe.pwi").read_text(encoding="utf-8")
    assert re.search(r"\becutwfc\s*=\s*100\b", qe)
    assert re.search(r"\becutrho\s*=\s*600\b", qe)
    assert "K_POINTS automatic\n4 4 3  0 0 0" in qe
    assert "PseudoDojo-NC-SR-PBE-v0.4-standard" in qe

    abinit = (inputs / "abinit.in").read_text(encoding="utf-8")
    assert "ecut 1400 eV" in abinit
    assert "ngkpt 4 4 3" in abinit
    assert "PseudoDojo-NC-SR-PBE-v0.4-standard" in abinit

    cp2k = (inputs / "cp2k.inp").read_text(encoding="utf-8")
    cutoff = re.search(r"\bCUTOFF \[eV\] ([\d.e+\-]+)", cp2k)
    assert cutoff is not None
    assert abs(float(cutoff.group(1)) / Ry - 800) < 1e-9
    assert "REL_CUTOFF 80" in cp2k
    assert "SCHEME MONKHORST-PACK 4 4 3" in cp2k
    assert "BASIS_SET DZVP-MOLOPT-SR-GTH" in cp2k
    assert "POTENTIAL GTH-PBE" in cp2k

    assert "QE 7.0 & PBE & PseudoDojo NC-SR UPF & 100/600 Ry" in MANUSCRIPT
    assert "ABINIT 8.6.1 & PBE & PseudoDojo NC-SR PSP8 & 1400 eV" in MANUSCRIPT
    assert "CP2K 2024.1 & PBE & GTH-PBE/DZVP-MOLOPT-SR & 800/80 Ry" in MANUSCRIPT
    assert "independently converged calculator contracts" not in MANUSCRIPT


def test_bto_prospective_errors_and_monotonic_rise_match_manuscript() -> None:
    points = _csv_rows(
        FIGURES / "bto_qy0_restricted_even_mode_sheet_points_source_data.csv"
    )
    blind = sorted(
        (row for row in points if row["class"] == "prospective"),
        key=lambda row: float(row["qz"]),
    )
    assert len(blind) == 2
    errors = [abs(float(row["observed_minus_predicted_meV_per_BTO"])) for row in blind]
    assert max(errors) < 2.0
    assert f"errors were {errors[0]:.3f} and {errors[1]:.3f} meV/BTO" in MANUSCRIPT

    paths = _csv_rows(FIGURES / "vcneb_path_comparison_source_data.csv")
    for n_images in (5, 7, 9):
        curve = sorted(
            (row for row in paths if row["material"] == "BaTiO3"
             and int(row["total_images"]) == n_images),
            key=lambda row: int(row["image_index"]),
        )
        assert len(curve) == n_images
        energies = [float(row["relative_enthalpy_eV_per_formula_unit"]) for row in curve]
        assert all(right >= left for left, right in zip(energies, energies[1:]))
        assert f"${energies[-1]:.7f}$ eV/BTO" in MANUSCRIPT


def test_bto_branch_alignment_and_holdout_claims_match_source() -> None:
    result = audit_patch(
        ROOT / "benchmarks/numerical_integrity/"
        "bto_conditional_five_point_patch_2026-09-28.json"
    )
    assert result["status"] == "five_point_metrics_recomputed_not_conditional_PES_certificate"
    assert f"${result['raw_selected_holdout_coordinate_error_sqrt_amu_A']:.3f}\\sqrt" in MANUSCRIPT
    assert f"${result['full_coordinate_error_sqrt_amu_A']:.3f}\\sqrt" in MANUSCRIPT
    assert f"${result['measured_minus_predicted_energy_meV_per_BTO']:.3f}$ meV/BTO" in MANUSCRIPT
    assert "not a continuous lower-envelope or local-minimum certificate" in MANUSCRIPT


def test_acceleration_and_hfo2_claims_match_audited_sources() -> None:
    audit = json.loads(
        (ROOT / "benchmarks/convergence/hf_slurm_abacus_launch_audit_20260928.json")
        .read_text(encoding="utf-8")
    )
    routes = audit["routes"]
    def count(name: str) -> int:
        return int(routes[name]["completed_launches_by_cutoff_second"])
    bto_before, bto_after = count("bto_global"), count("bto_block001")
    hfo2_before = count("hfo2_global")
    hfo2_after = count("hfo2_staged_coarse") + count("hfo2_staged_refine")
    bto_saved = 100 * (bto_before - bto_after) / bto_before
    hfo2_saved = 100 * (hfo2_before - hfo2_after) / hfo2_before
    assert f"{bto_saved:.1f}\\% and {hfo2_saved:.1f}\\%" in MANUSCRIPT
    assert (bto_before, bto_after, hfo2_before, hfo2_after) == (205, 44, 217, 134)

    barriers = _csv_rows(FIGURES / "vcneb_barrier_literature_comparison.csv")
    ci = next(
        row for row in barriers
        if row["material"] == "HfO2" and row["method"] == "CI n7"
    )
    per_cell = 4 * float(ci["barrier_eV_per_formula_unit"])
    assert f"{per_cell:.7f} eV per 12-atom cell" in MANUSCRIPT


def test_gan_original_endpoint_stress_claim_is_not_basin_return() -> None:
    from scripts.audit_gan_45p7_vasp_endpoints import audit

    original = audit()["endpoints"]
    assert [row["phase"] for row in original] == ["B4", "B1"]
    values = "/".join(f"{row['max_raw_stress_residual_kbar']:.3f}"
                      for row in original)
    assert values == "2.912/0.239"
    assert f"{values} kbar" in MANUSCRIPT
    assert "original B4 misses this diagnostic pressure gate" in MANUSCRIPT
    assert "separate signed basin-return" in MANUSCRIPT


def test_gan_strict_endpoint_diagnostic_is_not_barrier_replacement() -> None:
    from scripts.audit_gan_45p7_vasp_b4_refine import audit

    report = audit(PAPER / "evidence/gan_vasp_b4_endpoint_strict_20261001",
                   slurm_job=27812112)
    assert report["n_evaluated_ionic_frames"] == 5
    assert report["stress_pass_2kbar"]
    assert f"{report['final_max_raw_stress_residual_kbar']:.3f} kbar" in MANUSCRIPT
    assert "do not replace image~0 or recertify the reported VCNEB barrier" in MANUSCRIPT


def test_gan_normal_strain_step_check_is_claimed_within_audited_limits() -> None:
    report = json.loads((ROOT / "benchmarks/numerical_integrity/"
                         "gan_600eV_normal_strain_step_dependence_20261001.json")
                        .read_text(encoding="utf-8"))
    assert report["n_new_statics"] == 12
    magnitudes = [abs(value) for row in report["comparison_to_prior_0p02_A"]
                  for value in row["residual_eV_per_A"]]
    assert round(min(magnitudes), 4) == 0.0205
    assert round(max(magnitudes), 4) == 0.0238
    assert "$0.0205$--$0.0238$ eV/\\AA{}" in MANUSCRIPT
    assert "plausible contributor, not a uniquely established cause" in MANUSCRIPT


def test_gan_two_frozen_transverse_lowerings_match_path_force_slopes() -> None:
    numerical = ROOT / "benchmarks/numerical_integrity"
    grid = json.loads((numerical / "gan_600eV_atomic_tube_q9_20260929.json")
                      .read_text(encoding="utf-8"))
    components = json.loads((numerical / "gan_600eV_atomic_tube_components_20260928.json")
                            .read_text(encoding="utf-8"))
    trajectory = PAPER / "evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj"
    assert hashlib.sha256(trajectory.read_bytes()).hexdigest() == grid["source_sha256"]["trajectory"]
    assert grid["n_new_raw_audited_statics"] == 72
    q = grid["q_atom_A"]
    minus, center, plus = (q.index(value) for value in (-0.0125, 0.0, 0.0125))
    assert grid["central_image_indices"][10] == 15
    peak_h = grid["path_enthalpy_eV_per_cell"][10]
    forces = {row["image_index"]: row for row in components["anchor_components"]}
    for image, expected_drop, expected_gap in ((19, -0.236568957, 117.541829),
                                               (20, -0.205159309, 158.500804)):
        row = grid["central_image_indices"].index(image)
        values = grid["excess_enthalpy_meV_per_GaN"][row]
        assert abs(values[center]) < 1e-9
        assert abs(values[plus] - expected_drop) < 1e-6
        raw_cases = [item for item in grid["new_cases"]
                     if item["image_index"] == image and item["q_atom_A"] in (-0.0125, 0.0125)]
        assert len(raw_cases) == 2 and all(len(item["outcar_sha256"]) == 64
                                           for item in raw_cases)
        secant = (values[plus] - values[minus]) * 2 / 1000 / (q[plus] - q[minus])
        force_slope = forces[image]["path_force_transverse_slope_eV_per_A"]
        assert abs(secant - force_slope) < 0.0003
        gap = (peak_h - grid["path_enthalpy_eV_per_cell"][row]) * 500
        assert abs(gap - expected_gap) < 1e-5
        assert f"{secant:.5f}" in MANUSCRIPT
        assert f"{gap:.2f}" in MANUSCRIPT
    assert "No orthogonal atomic or cell degrees of freedom were released" in MANUSCRIPT


def test_quickstart_and_manual_describe_current_mode_surfaces() -> None:
    readme = re.sub(r"\s+", " ", (ROOT / "README.md").read_text(encoding="utf-8"))
    manual = re.sub(
        r"\s+", " ", (ROOT / "docs/USER_MANUAL.md").read_text(encoding="utf-8")
    )
    assert "27 audited conditional-DFT nodes on a 9×3 grid" in readme
    assert "17×17 grid of 289 *static DFT* points" in readme
    assert "144 off-path static evaluations (162 measured coordinates total)" in readme
    assert "9×3 sheet contains 27 audited conditional-DFT nodes" in manual
    assert "289 measured static ABACUS points on a complete 17×17 grid" in manual
    assert "18 path centers and 144 off-path VASP/600-eV statics" in manual
    assert "does not establish a lower relaxed MEP" in readme
    assert "neither certifies a lower fully relaxed MEP" in manual
