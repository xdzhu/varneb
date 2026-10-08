from dataclasses import replace
import json

from ase import Atoms
from ase.units import GPa
from scipy.constants import elementary_charge
import numpy as np
import pytest

from vcneb import ChannelPath, summarize_competing_paths, compare_channel_selectivity


REQUIRED = {"switch_a": "switching", "switch_b": "switching", "decay_T": "decay", "decay_M": "decay"}


def channel(name, barrier, *, reverse=False, **kwargs):
    images = [Atoms("Ar", scaled_positions=[[.2 + .1 * i, .3, .4]], cell=np.eye(3) * 5,
                    pbc=True) for i in range(3)]
    energies = [0., barrier, -.1]
    if reverse:
        images, energies = images[::-1], energies[::-1]
    fields = dict(name=name, role=REQUIRED[name], images=images, energies_eV_cell=energies,
                  physical_contract={"INPUT": "a" * 64, "KPT": "b" * 64},
                  mechanical_family="same_substrate_tilt_open",
                  mechanical_parameters={"biaxial_strain": 0., "field": 0.}, pressure_eV_A3=0.,
                  neb_fmax_eV_A=.08, source_id="synthetic_audited_cache", reverse=reverse,
                  barrier_error_eV_cell=.002, sampling_audit_sha256="c" * 64,
                  barrier_error_audit_sha256="d" * 64)
    fields.update(kwargs)
    return ChannelPath(**fields)


def paths(switch=.20, decay=.50):
    return [channel("switch_a", switch), channel("switch_b", switch + .10),
            channel("decay_T", decay + .1, reverse=True), channel("decay_M", decay)]


def summary(values, **kwargs):
    return summarize_competing_paths(values, formula_units=1, required_channels=REQUIRED, **kwargs)


def test_common_initial_reverse_direction_units_and_barrier_identity():
    report = summarize_competing_paths(paths(), formula_units=4, required_channels=REQUIRED)
    assert report["ready_for_bounded_discrete_comparison"]
    assert report["provisional_minimum_barriers_eV_fu"] == pytest.approx({"switching": .05, "decay": .125})
    assert report["provisional_selectivity_eV_fu"] == pytest.approx(.075)
    reverse = next(r for r in report["channels"] if r["name"] == "decay_T")
    assert reverse["source_direction_used"] == "reverse"
    assert reverse["source_forward_barrier_eV_fu"] == pytest.approx(.175)
    assert reverse["barrier_from_common_initial_eV_fu"] == pytest.approx(.15)
    # Source T->PO runs opposite to the displayed PO->T view.
    assert reverse["source_forward_barrier_eV_fu"] - reverse["source_reverse_barrier_eV_fu"] == pytest.approx(.025)
    assert reverse["reaction_energy_from_common_initial_eV_fu"] == pytest.approx(-.025)


def test_integer_periodic_initial_equivalence_is_not_atom_remapping():
    values = paths()
    moved = [a.copy() for a in values[1].images]
    for a in moved:
        a.positions += a.cell[0]  # same integer gauge along the entire chain
    values[1] = replace(values[1], images=moved)
    assert summary(values)["ready_for_bounded_discrete_comparison"]


@pytest.mark.parametrize("fault", ["energy_zero", "geometry", "contract", "condition", "family", "pressure", "order", "lift"])
def test_mixed_or_bad_chain_cannot_be_hidden_by_rebaselining(fault):
    values = paths()
    candidate = values[1]
    if fault == "energy_zero":
        candidate = replace(candidate, energies_eV_cell=np.array(candidate.energies_eV_cell) + 1e-4)
    elif fault in ("geometry", "order", "lift"):
        images = [a.copy() for a in candidate.images]
        if fault == "geometry": images[0].positions[0, 0] += .001
        if fault == "order": images[1].numbers[0] = 2
        if fault == "lift": images[1].positions += images[1].cell[0]
        candidate = replace(candidate, images=images)
    else:
        change = {"contract": {"physical_contract": {"INPUT": "e" * 64}},
                  "condition": {"mechanical_parameters": {"biaxial_strain": .01, "field": 0.}},
                  "family": {"mechanical_family": "free_cell"}, "pressure": {"pressure_eV_A3": .001}}[fault]
        candidate = replace(candidate, **change)
    values[1] = candidate
    with pytest.raises(ValueError): summary(values)


def test_pressure_work_uses_each_volume_and_can_move_discrete_peak():
    values = paths()
    values = [replace(p, pressure_eV_A3=.01) for p in values]
    images = [a.copy() for a in values[0].images]
    for a, side in zip(images, (5., 6., 7.)):
        a.set_cell(np.eye(3) * side, scale_atoms=True)
    values[0] = replace(values[0], images=images, energies_eV_cell=[0., 1., .2])
    report = summary(values)
    first = report["channels"][0]
    assert first["highest_image_in_view"] == 2
    assert first["barrier_from_common_initial_eV_fu"] == pytest.approx(.2 + .01 * (343 - 125))
    assert first["source_reverse_barrier_eV_fu"] == 0.


def test_pressure_volume_normalization_matches_independent_SI_work():
    values = [replace(p, pressure_eV_A3=GPa) for p in paths()]
    images = [a.copy() for a in values[0].images]
    for a, side in zip(images, (5., 6., 7.)):
        a.set_cell(np.eye(3) * side, scale_atoms=True)
    values[0] = replace(values[0], images=images, energies_eV_cell=[0., 0., 0.])
    report = summarize_competing_paths(values, formula_units=4, required_channels=REQUIRED)
    # SI pressure*volume is joules; division by the library electron charge
    # converts to eV. ASE's pinned CODATA convention can differ in the 8th
    # digit; this dimensional check does not retune that production constant.
    independent_eV_fu = 1e9 * ((7e-10)**3 - (5e-10)**3) / elementary_charge / 4
    assert report["channels"][0]["barrier_from_common_initial_eV_fu"] == pytest.approx(independent_eV_fu, rel=5e-8)


def test_missing_leakage_channel_blocks_minimum_and_scientific_comparison():
    report = summary(paths()[:-1])
    assert report["missing_required_channels"] == ["decay_M"]
    assert report["provisional_minimum_barriers_eV_fu"]["decay"] is None
    assert report["provisional_selectivity_eV_fu"] is None
    comparison = compare_channel_selectivity(report, report, decay_noninferiority_margin_eV_fu=0.)
    assert comparison["status"] == "unresolved" and not comparison["claims_evaluated"]


@pytest.mark.parametrize("omission", ["sampling", "error", "error_provenance", "convergence"])
def test_force_threshold_does_not_supply_sampling_or_error_evidence(omission):
    values = paths()
    update = {"sampling": {"sampling_audit_sha256": None}, "error": {"barrier_error_eV_cell": None},
              "error_provenance": {"barrier_error_audit_sha256": None}, "convergence": {"neb_fmax_eV_A": .100001}}[omission]
    values[0] = replace(values[0], **update)
    report = summary(values)
    assert not report["ready_for_bounded_discrete_comparison"]
    assert not compare_channel_selectivity(report, summary(paths()), decay_noninferiority_margin_eV_fu=0.)["claims_evaluated"]


def test_relative_selectivity_improvement_can_make_decay_easier():
    # Both barriers fall, switching more strongly: selectivity improves but
    # the stronger noninferiority claim is false. No retention inference.
    base = summary(paths(switch=.3, decay=.5))
    changed = summary(paths(switch=.1, decay=.4))
    result = compare_channel_selectivity(base, changed, decay_noninferiority_margin_eV_fu=0.)
    assert result["relative_selectivity_improves"] and result["switching_barrier_lowers"]
    assert result["decay_barrier_lowers"]
    assert not result["easier_switching_without_decay_loss_beyond_margin"]


def test_numpy_formula_count_still_returns_JSON_serializable_flags():
    base = summarize_competing_paths(paths(.3, .5), formula_units=np.int64(4), required_channels=REQUIRED)
    changed = summarize_competing_paths(paths(.1, .6), formula_units=np.int64(4), required_channels=REQUIRED)
    result = compare_channel_selectivity(base, changed, decay_noninferiority_margin_eV_fu=0.)
    assert json.loads(json.dumps(result))["easier_switching_without_decay_loss_beyond_margin"] is True


def test_stronger_decoupling_and_uncertainty_overlap_are_distinct():
    base = summary(paths(switch=.3, decay=.5))
    result = compare_channel_selectivity(base, summary(paths(switch=.1, decay=.6)), decay_noninferiority_margin_eV_fu=0.)
    assert result["decay_barrier_increases"] and result["easier_switching_without_decay_loss_beyond_margin"]
    unchanged = compare_channel_selectivity(base, base, decay_noninferiority_margin_eV_fu=0.)
    assert not unchanged["decay_noninferior_within_margin"]  # overlap is not noninferiority proof
    declared = compare_channel_selectivity(base, base, decay_noninferiority_margin_eV_fu=.005)
    assert declared["decay_noninferior_within_margin"]
    assert not declared["relative_selectivity_improves"]


def test_minimum_interval_allows_the_lowest_channel_to_change():
    values = paths()
    values[0] = replace(values[0], energies_eV_cell=[0., .20, -.1], barrier_error_eV_cell=.10)
    values[1] = replace(values[1], energies_eV_cell=[0., .21, -.1], barrier_error_eV_cell=.001)
    report = summary(values)
    assert report["minimum_barrier_intervals_eV_fu"]["switching"] == pytest.approx([.10, .211])


@pytest.mark.parametrize("fault", ["family", "contract", "pressure", "coverage", "parameters"])
def test_free_and_clamped_results_are_not_a_smooth_strain_comparison(fault):
    base = summary(paths())
    changed = summary(paths())
    updates = {"family": ("mechanical_family", "free_cell"), "contract": ("physical_contract", {"INPUT": "f" * 64}),
               "pressure": ("pressure_eV_A3", .01), "coverage": ("required_channels", {"switch_a": "switching"}),
               "parameters": ("mechanical_parameters", {"different_parameter": .01})}
    key, value = updates[fault]
    changed[key] = value
    with pytest.raises(ValueError):
        compare_channel_selectivity(base, changed, decay_noninferiority_margin_eV_fu=0.)


@pytest.mark.parametrize("kwargs", [{"formula_units": True}, {"formula_units": 1.5},
                                    {"formula_units": 0}, {"geometry_tolerance_A": 0.},
                                    {"fmax_target_eV_A": float("nan")}])
def test_invalid_units_or_tolerances_are_rejected(kwargs):
    fields = dict(formula_units=1, required_channels=REQUIRED)
    fields.update(kwargs)
    with pytest.raises(ValueError): summarize_competing_paths(paths(), **fields)
