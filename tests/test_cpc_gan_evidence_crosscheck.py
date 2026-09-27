"""Bind the published five-backend GaN panel to its evidence manifest."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "VARNEB_CPC"
MANIFEST = PAPER / "evidence" / "gan_45p7_multibackend_vcneb_20260924.json"
SOURCE = PAPER / "figures" / "gan_multibackend_validation_source_data.csv"


def test_gan_plot_matches_each_backend_evidence_record() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    with SOURCE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    contract = manifest["contract"]
    assert contract["n_images_total"] == 29
    assert contract["n_images_interior"] == 27
    assert contract["pressure_gpa"] == 45.7
    assert contract["fmax_target_eV_per_A"] == 0.10
    assert set(row["backend"] for row in rows) == {
        "ABACUS", "VASP", "QE", "ABINIT", "CP2K", "Qian et al."
    }

    for backend, record in manifest["backends"].items():
        curve = sorted(
            (row for row in rows if row["backend"] == backend.upper()),
            key=lambda row: int(row["image_index"]),
        )
        assert record["status"] == "converged"
        assert len(curve) == contract["n_images_total"]
        assert [int(row["image_index"]) for row in curve] == list(range(29))
        assert all(row["status"] == "converged" for row in curve)
        assert all(str(record["job_id"]) in row["provenance"] for row in curve)

        relative = np.array([float(row["relative_enthalpy_eV_per_GaN"]) for row in curve])
        assert np.isfinite(relative).all()
        np.testing.assert_allclose(relative[0], 0.0, atol=1e-11)
        assert int(np.argmax(relative)) == record["highest_image_index"] == 15
        np.testing.assert_allclose(
            relative.max(), record["barrier_eV_per_GaN"], atol=5e-10
        )
        np.testing.assert_allclose(
            relative[-1], record["reaction_enthalpy_eV_per_cell"] / 2, atol=5e-10
        )

        reverse = relative.max() - relative[-1]
        for row in curve:
            np.testing.assert_allclose(
                float(row["forward_barrier_eV_per_GaN"]), relative.max(), atol=5e-10
            )
            np.testing.assert_allclose(
                float(row["reverse_barrier_eV_per_GaN"]), reverse, atol=5e-10
            )
            np.testing.assert_allclose(
                float(row["final_fmax_eV_per_A"]),
                record["final_max_generalized_force_eV_per_A"],
                atol=5e-10,
            )

        if "relative_enthalpy_eV_per_cell" in record:
            np.testing.assert_allclose(
                relative,
                np.asarray(record["relative_enthalpy_eV_per_cell"]) / 2,
                atol=5e-10,
            )

    literature = [row for row in rows if row["backend"] == "Qian et al."]
    assert len(literature) == 29
    assert all(not row["reverse_barrier_eV_per_GaN"] for row in literature)
