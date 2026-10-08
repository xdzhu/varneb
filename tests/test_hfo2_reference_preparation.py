from pathlib import Path
import json
import shutil

import numpy as np
import pytest
from ase.io import read

from scripts.prepare_hfo2_gamma_and_seeds import gamma_family
from scripts.prepare_hfo2_reference_variants import distortion, fluorite_scaffold, translation_sectors
from scripts.audit_hfo2_static_replica import sha256
import scripts.audit_hfo2_gamma_and_seeds as scorer
from vcneb.reference_variants import apply_parent_operation
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.prepare_hfo2_sparse_gap import gap_geometries


DATA = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/reference_variants"


def test_two_sizes_are_both_central_and_do_not_expand_or_relabel():
    t = read(DATA / "T.vasp", format="vasp")
    for size in (.01, .02):
        family = gamma_family(t, size)
        assert len(family.supercells_with_displacements) == 8
        np.testing.assert_array_equal(family.supercell_matrix, np.eye(3))
        for i in range(0, 8, 2):
            plus, minus = family.dataset["first_atoms"][i:i+2]
            assert plus["number"] == minus["number"]
            np.testing.assert_allclose(plus["displacement"], -np.asarray(minus["displacement"]), atol=1e-12)


def test_parent_translation_sectors_reconstruct_and_are_orthogonal():
    t, po = [read(DATA / f"{label}.vasp") for label in ("T", "PO")]
    parent = fluorite_scaffold(t)
    values = distortion(po, parent)
    sectors = translation_sectors(parent, values)
    np.testing.assert_allclose(sum(sectors.values()), values, atol=1e-12)
    vectors = np.stack(list(sectors.values())).reshape(4, 36)
    np.testing.assert_allclose(vectors @ vectors.T, np.diag(np.sum(vectors**2, axis=1)), atol=1e-12)
    t_sectors = translation_sectors(parent, distortion(t, parent))
    assert np.linalg.norm(t_sectors["X_z"]) > .8
    assert max(np.linalg.norm(t_sectors[k]) for k in ("Gamma", "X_x", "X_y")) < 1e-10


def test_inverted_candidate_reverses_gamma_not_assumed_all_x_irreps():
    t, po = [read(DATA / f"{label}.vasp") for label in ("T", "PO")]
    parent = fluorite_scaffold(t)
    sectors = translation_sectors(parent, distortion(po, parent))
    flipped = apply_parent_operation(parent, po, -np.eye(3), [0, 0, 0]).atoms
    negative = translation_sectors(parent, distortion(flipped, parent))
    np.testing.assert_allclose(negative["Gamma"], -sectors["Gamma"], atol=2e-5)
    assert not np.allclose(negative["X_x"], -sectors["X_x"], atol=1e-4)


def test_unreviewed_steps_and_unsupported_t_embedding_fail():
    t = read(DATA / "T.vasp")
    with pytest.raises(ValueError, match="two declared"):
        gamma_family(t, .03)
    t.positions[0, 0] += .01
    with pytest.raises(ValueError, match="embedding"):
        fluorite_scaffold(t)


def test_full_scorer_on_synthetic_harmonic_forces(tmp_path, monkeypatch):
    # Analytic fixture only: no DFT result or HfO2 physical claim.
    shutil.copyfile(DATA / "T.vasp", tmp_path / "T_reference.vasp")
    t = read(tmp_path / "T_reference.vasp")
    records, results = [], []
    for distance in (.01, .02):
        family = gamma_family(t, distance)
        for local_index, item in enumerate(family.dataset["first_atoms"]):
            index = len(records)
            records.append({"index": index, "kind": "T_Gamma", "distance_A": distance,
                            "local_index": local_index, "atom_index": item["number"],
                            "displacement_A": item["displacement"]})
            delta = np.zeros((12, 3))
            delta[item["number"]] = item["displacement"]
            forces = -2 * (delta - delta.mean(axis=0))
            results.append({"index": index, "forces_eV_A": forces.tolist()})
    for name in ("minus_a", "minus_b", "M"):
        index = len(records)
        records.append({"index": index, "kind": "endpoint_candidate_static", "name": name})
        results.append({"index": index, "forces_eV_A": np.zeros((12, 3)).tolist()})
    manifest = {"n_static_points": 19, "points": records,
                "T_reference_sha256": sha256(tmp_path / "T_reference.vasp"), "Gamma_boundary": "synthetic"}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    monkeypatch.setattr(scorer, "point_audit", lambda root, index: results[index])
    report = scorer.audit(tmp_path, tmp_path / "analysis")
    assert len(report["points"]) == 19
    assert report["step_size_comparison"]["max_raw_fc_difference_eV_A2"] < 1e-12
    assert all(not family["negative_optical_indices"] for family in report["families"].values())
    assert all(family["minimum_translation_principal_overlap"] > .999 for family in report["families"].values())
    with pytest.raises(FileExistsError):
        scorer.audit(tmp_path, tmp_path / "analysis")


def test_rotated_t_triplet_resolves_axis_without_irrep_assumption():
    t = read(DATA / "T.vasp")
    parent, basis, _ = rotated_t_triplet(t)
    vector = distortion(t, parent)
    projections = np.einsum("aij,ij->a", basis, vector)
    np.testing.assert_allclose(projections[:2], 0, atol=1e-10)
    assert projections[2] > .8
    po = distortion(read(DATA / "PO.vasp"), parent)
    assert abs(np.sum(po * basis[0])) > .5


def test_sparse_gap_adds_only_three_fixed_order_geometries():
    t = read(DATA / "T.vasp")
    last = t.copy()
    last.positions[4, 2] += .04
    added = gap_geometries(t, last)
    assert len(added) == 3
    for fraction, atoms in zip((.25, .5, .75), added):
        assert atoms.get_chemical_symbols() == t.get_chemical_symbols()
        np.testing.assert_allclose(atoms.cell.array, t.cell.array)
        np.testing.assert_allclose(atoms.positions - t.positions,
                                   fraction * (last.positions - t.positions), atol=1e-12)
