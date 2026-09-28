"""Keep headline CPC numbers tied to committed figure and launch sources."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


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
