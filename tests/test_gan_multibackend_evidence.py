"""Keep the published CP2K GaN curve tied to its evaluated final chain."""

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


def test_cp2k_gan_final_chain_matches_multibackend_evidence() -> None:
    payload = json.loads((EVIDENCE / "gan_45p7_multibackend_vcneb_20260924.json").read_text())
    record = payload["backends"]["cp2k"]
    snapshot = EVIDENCE / "gan_cp2k_45p7_chain_step_0052.traj"
    images = read(snapshot, index=":")

    assert record["status"] == "converged"
    assert record["final_step"] == 52
    assert record["final_max_generalized_force_eV_per_A"] < payload["contract"]["fmax_target_eV_per_A"]
    assert record["geometry_valid"]
    assert hashlib.sha256(snapshot.read_bytes()).hexdigest() == record["chain_snapshot_sha256"]
    assert len(images) == payload["contract"]["n_images_total"] == 29
    assert all(len(image) == 4 and image.get_chemical_formula() == "Ga2N2" for image in images)
    assert all(np.isfinite(image.get_forces()).all() and np.isfinite(image.get_stress()).all()
               for image in images)

    enthalpies = np.asarray([
        image.get_potential_energy() + payload["contract"]["pressure_gpa"] * GPa * image.get_volume()
        for image in images
    ])
    relative = enthalpies - enthalpies[0]
    np.testing.assert_allclose(relative, record["relative_enthalpy_eV_per_cell"], atol=1e-10)
    assert int(np.argmax(relative)) == record["highest_image_index"] == 15
    np.testing.assert_allclose(relative.max(), record["barrier_eV_per_cell"], atol=1e-10)
    np.testing.assert_allclose(relative.max() / 2, record["barrier_eV_per_GaN"], atol=1e-10)


def test_multibackend_source_data_has_five_converged_paths() -> None:
    source = ROOT / "paper" / "VARNEB_CPC" / "figures" / "gan_multibackend_validation_source_data.csv"
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for backend in ("ABACUS", "VASP", "QE", "ABINIT", "CP2K"):
        path_rows = [row for row in rows if row["backend"] == backend]
        assert len(path_rows) == 29
        assert all(row["status"] == "converged" for row in path_rows)
        assert all(row["relative_enthalpy_eV_per_GaN"] for row in path_rows)
