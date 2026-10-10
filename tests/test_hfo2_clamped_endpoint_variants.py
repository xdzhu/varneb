"""Declared parent mapping must not erase actual ordered path identity."""
from pathlib import Path

from ase.io import read
import numpy as np
import pytest

from scripts.audit_hfo2_clamped_endpoint_variants import audit, declared_translation_pair, geometric_patterns
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.analyze_hfo2_network_update import structure_audit
from examples.hfo2_fixed_input_factory import read_fixed_hfo2_stru
from scripts.prepare_hfo2_clamped_endpoints import AXES, orient_long_axis_x

CASE = Path(__file__).resolve().parents[1]/"benchmarks/hfo2_channels/20261008"
P = [2,3,0,1,9,8,11,10,5,4,7,6]


def seed(phase):
    return read(CASE/"clamped_endpoint_seeds/strain_0000"/phase/"POSCAR.seed", format="vasp")


def test_actual_declared_translation_is_not_an_ordered_cache_match():
    result = declared_translation_pair(seed("PO_minus_T_preserving"), seed("PO_minus_T_reversing"), [.5,.5,0], P)
    assert result["same_periodic_geometry_under_declared_translation_and_permutation"]
    assert result["declared_operation_max_distance_A"] < 1e-12
    assert result["ordered_periodic_max_distance_A"] > 1.
    assert not result["same_ordered_periodic_geometry"]
    assert not result["production_atoms_wrapped_translated_permuted_or_modified"]
    assert not result["whole_path_equivalence_topology_or_electronic_polarization_certified"]


def test_wrong_parent_operation_not_silently_best_fit():
    a,b = seed("PO_minus_T_preserving"),seed("PO_minus_T_reversing")
    r = declared_translation_pair(a,b,[0,0,0],list(range(12)))
    assert not r["same_periodic_geometry_under_declared_translation_and_permutation"]


@pytest.mark.parametrize("translation,permutation,tolerance", [
    ([np.nan,0,0],P,1e-8), ([0,0],P,1e-8), ([.5,.5,0],[0]*12,1e-8),
    ([.5,.5,0],[float(x) for x in P],1e-8), ([.5,.5,0],[4,1,2,3,0,5,6,7,8,9,10,11],1e-8),
    ([.5,.5,0],P,0), ([.5,.5,0],P,np.nan)])
def test_invalid_declared_operation_is_refused(translation,permutation,tolerance):
    with pytest.raises(ValueError):
        declared_translation_pair(seed("PO_minus_T_preserving"),seed("PO_minus_T_reversing"),
                                  translation,permutation,tolerance_A=tolerance)


def test_geometric_pattern_chart_boundary_is_not_smoothed_or_rewrapped():
    t = read(CASE/"reference_variants/T.vasp",format="vasp")
    parent,basis,_ = rotated_t_triplet(t)
    parent = orient_long_axis_x(parent)
    outside = parent.copy()
    coordinates = outside.get_scaled_positions(wrap=False)
    coordinates[4,0] += .46
    outside.set_scaled_positions(coordinates)
    r = geometric_patterns(outside,parent,basis[...,list(AXES)])
    assert not r["nearest_reference_chart_accepted"] and r["Q_A"] is None
    assert r["chart_margin_fractional"] == pytest.approx(.04)


def test_actual_five_endpoint_review_replays_raw_covariance_without_new_DFT():
    r = audit(CASE)
    assert r["new_DFT_calls"] == 0 and not r["holdout_generated_or_read"]
    assert len(r["endpoints"]) == 5 and r["substrate_strain"] == r["pressure_GPa"] == 0.
    for pair in r["minus_pair"].values():
        assert pair["same_periodic_geometry_under_declared_translation_and_permutation"]
        assert not pair["same_ordered_periodic_geometry"]
        assert pair["source_to_target"] == P
    covariance = r["raw_covariance"]
    assert abs(covariance["terminal_energy_difference_eV_cell"]) < 1e-10
    assert covariance["terminal_force_covariance_max_absolute_error_eV_A"] < 1e-9
    assert covariance["terminal_stress_covariance_max_absolute_error_eV_A3"] < 1e-11
    e = r["endpoints"]
    assert e["PO_plus"]["terminal_patterns"]["oxygen_minus_hafnium_mean_displacement_A"][0] > 0
    for phase in ("PO_minus_T_preserving","PO_minus_T_reversing"):
        assert e[phase]["terminal_patterns"]["oxygen_minus_hafnium_mean_displacement_A"][0] < 0
        assert [x["symbol"] for x in e[phase]["structure_audit"]["symmetry_sweep"]] == ["Pca2_1"]*3
    assert e["PO_minus_T_preserving"]["terminal_patterns"]["Q_A"][0] < -.9
    assert e["PO_minus_T_reversing"]["terminal_patterns"]["Q_A"][0] > .9
    assert r["terminal_INPUT_KPT_bytes_checked_here"]
    assert not r["licensed_pseudo_orbital_bytes_checked_here"]


def test_registered_plus_one_T_is_metric_lowered_before_DFT_not_a_forced_tetragonal_phase():
    old = seed("T")
    strained = read(CASE/"clamped_endpoint_seeds/strain_p0100/T/POSCAR.seed",format="vasp")
    assert abs(np.linalg.norm(old.cell[1])-np.linalg.norm(old.cell[2])) < 1e-10
    assert abs(np.linalg.norm(strained.cell[1])-np.linalg.norm(strained.cell[2])) > .05
    assert [s["symbol"] for s in structure_audit(strained)["symmetry_sweep"]] == ["Ccce"]*3
    t = read(CASE/"reference_variants/T.vasp",format="vasp")
    parent,basis,_ = rotated_t_triplet(t)
    parent = orient_long_axis_x(parent);basis = basis[...,list(AXES)]
    final = read_fixed_hfo2_stru(CASE/"clamped_endpoint_matrix_20261009/strain_p0100/T/completed_HF/endpoint/calculator/image_0000/scf_000006/STRU")
    for atoms in (strained,final):
        r = geometric_patterns(atoms,parent,basis)
        assert r["nearest_reference_chart_accepted"]
        assert np.max(np.abs(r["Q_A"][:2])) < 1e-12
        assert r["orthogonal_residual_norm_A"] < 2e-12
    r = geometric_patterns(final,parent,basis)
    assert r["Q_A"][2] == pytest.approx(.8857526497026532,abs=1e-12,rel=0)
