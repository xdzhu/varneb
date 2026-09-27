"""Pin all five GaN figure curves to portable evaluated-chain snapshots."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "paper" / "VARNEB_CPC" / "evidence"
FIGURE = ROOT / "paper" / "VARNEB_CPC" / "figures" / "gan_multibackend_validation_source_data.csv"
AUDIT = EVIDENCE / "gan_45p7_final_chain_audit_20260927.json"


def test_five_compact_chains_reproduce_the_published_enthalpy_curves() -> None:
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    compact = json.loads((EVIDENCE / "gan_45p7_multibackend_vcneb_20260924.json").read_text())
    with FIGURE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert audit["status"] == "five_hashed_final_trajectories_match_plotted_enthalpy_curves"
    assert audit["pressure_gpa"] == compact["contract"]["pressure_gpa"] == 45.7
    assert audit["n_images_total"] == compact["contract"]["n_images_total"] == 29
    assert set(audit["backends"]) == {"abacus", "vasp", "qe", "abinit", "cp2k"}
    for name, record in audit["backends"].items():
        path = EVIDENCE / "gan_45p7_final_chains_20260927" / record["compact_final_chain_file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["compact_final_chain_sha256"]
        images = read(path, index=":")
        assert len(images) == 29
        assert all(len(image) == 4 and image.get_chemical_formula() == "Ga2N2"
                   and np.isfinite(image.get_forces()).all()
                   and np.isfinite(image.get_stress()).all() for image in images)
        enthalpy = np.asarray([image.get_potential_energy() + 45.7 * GPa * image.get_volume()
                               for image in images])
        relative = enthalpy - enthalpy[0]
        plotted = np.asarray([float(row["relative_enthalpy_eV_per_GaN"]) * 2
                              for row in rows if row["backend"] == name.upper()])
        np.testing.assert_allclose(relative, plotted, atol=3e-9, rtol=0)
        assert int(np.argmax(relative)) == record["peak_image_index"] == 15
        np.testing.assert_allclose(relative.max() / 2, record["barrier_eV_per_GaN"], atol=1e-9)
        np.testing.assert_allclose(relative[-1] / 2, record["reaction_enthalpy_eV_per_GaN"], atol=1e-9)
        np.testing.assert_allclose(relative.max() / 2,
                                   compact["backends"][name]["barrier_eV_per_GaN"], atol=1e-9)
        assert record["maximum_figure_curve_difference_eV_per_cell"] < 3e-9
        assert record["all_final_image_energies_forces_stresses_finite"] is True
