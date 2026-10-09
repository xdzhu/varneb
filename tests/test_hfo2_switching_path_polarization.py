"""Bounded property recipe and real frozen-band provenance; no local DFT."""
import copy
import json
from pathlib import Path

from ase import Atoms
import numpy as np
import pytest

from scripts import hfo2_switching_path_polarization as flow
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_channel_work_probes import write_structure
from vcneb.polarization import sampled_reduced_branch


REPOSITORY = Path(__file__).resolve().parents[1]


def test_actual_blueprint_reuses_fourteen_interiors_and_three_endpoints():
    plan = flow.blueprint(REPOSITORY)
    assert len(plan["points"]) == 14
    assert plan["maximum_new_SCF_calls"] == 14 and plan["maximum_new_NSCF_calls"] == 42
    assert [(p["channel"], p["image_index"]) for p in plan["points"]] == [
        (c, i) for c in ("preserving", "reversing") for i in range(1, 8)]
    assert all(c["ordinary_fmax_eV_A"] < .1 and c["max_R3_transverse_fraction"] < 2e-5
               for c in plan["channels"])
    assert not plan["SCF_physical_settings_changed"] and not plan["NSCF_energies_used_for_barriers"]
    assert plan["new_DFT_calls_for_preparation"] == 0 and not plan["automatic_restart_or_G2_submission"]


def test_reduced_lift_uses_actual_quantum_not_fixed_cartesian_period():
    quanta = np.array([1., .8, .6, .5])
    reduced = np.array([.8, -.9, -.6, -.4])
    result = sampled_reduced_branch(reduced*quanta, quanta, 2*quanta, np.zeros(4))
    np.testing.assert_allclose(result["lifted_reduced"], [.8, 1.1, 1.4, 1.6])
    assert result["branch_integers"] == [0, 1, 1, 1]
    assert result["delta_reduced"] == pytest.approx(.8)
    assert not result["continuous_path_certified"] and not result["spontaneous_P_selected"]
    # The selected start is explicit; changing it changes the gauge, not delta.
    other = sampled_reduced_branch(reduced*quanta, quanta, 2*quanta, np.zeros(4), initial_integer=-1)
    np.testing.assert_allclose(np.array(result["lifted_reduced"])-other["lifted_reduced"], 2.)
    assert other["delta_reduced"] == pytest.approx(result["delta_reduced"])


@pytest.mark.parametrize("values,errors", [([0., 1.], [0., 0.]), ([0., .99], [.01, .01])])
def test_branch_abstains_on_tie_or_measured_uncertainty_overlap(values, errors):
    result = sampled_reduced_branch(values, [1., 1.], [2., 2.], errors)
    assert result["status"] == "ambiguous_sampled_link"
    assert result["lifted_reduced"] is None and result["branch_integers"] is None
    assert result["ambiguous_link"] == [0, 1]


@pytest.mark.parametrize("quanta,periods,errors", [([0., 1.], [2., 2.], [0., 0.]),
    ([1., 1.], [2., 1.], [0., 0.]), ([1., 1.], [2., 2.], [-.1, 0.]),
    ([1., 1.], [2., 2.], [np.nan, 0.])])
def test_branch_invalid_units_and_periods_refused(quanta, periods, errors):
    with pytest.raises(ValueError):
        sampled_reduced_branch([0., .2], quanta, periods, errors)


@pytest.fixture
def staged_fixture(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    base = b"INPUT_PARAMETERS\ncalculation scf\nbasis_type lcao\necutwfc 100.0\ncal_force 1\ncal_stress 1\nout_stru 1\n"
    (source / "INPUT").write_bytes(base)
    (source / "KPT").write_bytes(b"K_POINTS\n0\nGamma\n2 2 2 0 0 0\n")
    for name in flow.INPUT_FILES:
        if name not in ("INPUT", "KPT", "STRU"):
            (source / name).write_text("fixture "+name)
    atoms = Atoms(["Hf"]*4+["O"]*8, positions=np.arange(36).reshape(12, 3)*.1,
                  cell=np.eye(3)*5, pbc=True)
    write_structure(source / "STRU", atoms)
    (source / "OUT.ABACUS").mkdir()
    (source / "OUT.ABACUS/running_scf.log").write_text("fixture only, not physical SCF")
    raw = {"energy": -10., "forces": np.zeros((12, 3)), "stress": np.zeros(6)}
    contract = {n: sha256(source/n) for n in flow.CONTRACT}
    monkeypatch.setattr(flow, "CONTRACT", contract)
    plan = {"purpose": "fixture", "physical_contract": contract, "tolerances": flow.TOLERANCES,
            "points": [], "maximum_new_SCF_calls": 14, "maximum_new_NSCF_calls": 42,
            "endpoint_summary_sha256": None}
    endpoint = tmp_path / "endpoint.json"
    endpoint.write_text("{}\n")
    plan["endpoint_summary_sha256"] = sha256(endpoint)
    for i in range(14):
        plan["points"].append({"index": i, "label": f"test_{i:02d}", "source_directory": str(source),
            "baseline_log_sha256": sha256(source/"OUT.ABACUS/running_scf.log"),
            "baseline_input_sha256": {n: sha256(source/n) for n in flow.INPUT_FILES},
            "baseline_results": {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in raw.items()},
            "cell_A": atoms.cell.array.tolist(), "positions_A": atoms.positions.tolist(),
            "quantum_lattice_C_m2": flow.quantum_lattice(atoms.cell.array).tolist()})
    monkeypatch.setattr(flow, "blueprint", lambda _: copy.deepcopy(plan))
    monkeypatch.setattr(flow.native, "read", lambda *a, **k: atoms.copy())
    monkeypatch.setattr(flow, "audited_results", lambda _: raw)

    def endpoint_replay(repository, root):
        (root/"endpoint_summary.json").write_bytes(endpoint.read_bytes())
        return {"status": "synthetic protocol fixture", "new_DFT_calls": 0}
    monkeypatch.setattr(flow, "replay_endpoints", endpoint_replay)
    root = tmp_path / "prepared"
    manifest = flow.prepare(tmp_path, root)
    return root, manifest, base


def test_preparation_only_changes_observable_flags_and_quadrature(staged_fixture):
    root, manifest, base = staged_fixture
    for point in manifest["points"]:
        scf = root/"points"/point["label"]/"scf"
        assert (scf/"INPUT").read_bytes() == base + flow.native.OUTPUT_ADDITION
        assert point["baseline_input_sha256"]["STRU"] == point["scf_input_sha256"]["STRU"]
        assert [n["nz"] for n in point["nscf_inputs"]] == [2, 4, 8]
    with pytest.raises(FileExistsError):
        flow.prepare(root.parent, root)


@pytest.mark.parametrize("mutation", ["tolerance", "budget", "geometry", "cutoff", "cache"])
def test_manifest_and_hash_consistent_parameter_tampering_refused(staged_fixture, mutation):
    root, m, base = staged_fixture
    if mutation == "tolerance":
        m["tolerances"]["SCF_energy_eV_cell"] = .1
    elif mutation == "budget":
        m["maximum_new_NSCF_calls"] = 420
    elif mutation == "cache":
        m["points"][0]["cached_scf_directory"] = "arbitrary"
    else:
        p = m["points"][0]
        path = root/"points"/p["label"]/"scf"/("STRU" if mutation == "geometry" else "INPUT")
        path.write_bytes(path.read_bytes()+b"\n" if mutation == "geometry" else
                         path.read_bytes().replace(b"ecutwfc 100.0", b"ecutwfc 80.0"))
        p["scf_input_sha256"][path.name] = sha256(path)
    flow.save(root/"manifest.json", m)
    with pytest.raises(ValueError):
        flow.load_manifest(root.parent, root)


def test_dft_requires_checked_allocation(staged_fixture, monkeypatch):
    root, _, _ = staged_fixture
    monkeypatch.delenv("RUN_DFT", raising=False)
    with pytest.raises(RuntimeError, match="allocation"):
        flow.run(root.parent, root, 0)
    with pytest.raises(ValueError, match="index"):
        flow.run(root.parent, root, 14)


def test_postprocessing_compatibility_is_one_exact_historical_driver_not_a_hash_bypass():
    expected = flow.blueprint(REPOSITORY)
    actual = copy.deepcopy(expected)
    actual["driver_sha256"] = flow.HISTORICAL_DRIVER_SHA256
    flow.validate_blueprint_metadata(actual, expected)
    actual["tolerances"]["SCF_energy_eV_cell"] = .1
    with pytest.raises(ValueError, match="tolerances"):
        flow.validate_blueprint_metadata(actual, expected)
    actual = copy.deepcopy(expected)
    actual["driver_sha256"] = "0"*64
    with pytest.raises(ValueError, match="driver"):
        flow.validate_blueprint_metadata(actual, expected)


def test_exact_known_native_CRLF_record_is_not_arbitrary_auditor_compatibility():
    expected = flow.blueprint(REPOSITORY)
    assert expected["native_adapter_sha256"] == flow.NATIVE_ADAPTER_LF_SHA256
    actual = copy.deepcopy(expected)
    actual["native_adapter_sha256"] = flow.NATIVE_ADAPTER_HISTORICAL_CRLF_SHA256
    flow.validate_blueprint_metadata(actual, expected)
    actual["native_adapter_sha256"] = "0"*64
    with pytest.raises(ValueError, match="native_adapter"):
        flow.validate_blueprint_metadata(actual, expected)
    expected["native_adapter_sha256"] = "1"*64
    actual["native_adapter_sha256"] = flow.NATIVE_ADAPTER_HISTORICAL_CRLF_SHA256
    with pytest.raises(ValueError, match="native_adapter"):
        flow.validate_blueprint_metadata(actual, expected)


def test_end_to_end_synthetic_summary_serializes_numpy_comparison_result(tmp_path, monkeypatch):
    points = [{"label": f"{name}_{i:02d}"} for name in ("preserving", "reversing") for i in range(1, 8)]
    channels = [{"name": name, "endpoint_labels": ["PO_plus", f"PO_minus_T_{name}"]}
                for name in ("preserving", "reversing")]
    monkeypatch.setattr(flow, "load_manifest", lambda *a: {"points": points, "channels": channels})
    (tmp_path/"manifest.json").write_text("{}\n")

    def records(value):
        return [{"nz": n, "value_modern_SI_C_m2": value, "physical_quantum_C_m2": 1.,
                 "reported_modulus_modern_SI_C_m2": 2.} for n in flow.native.GRIDS]
    endpoints = {l: records(.4 if l == "PO_plus" else -.4) for l in flow.native.LABELS}
    (tmp_path/"endpoint_summary.json").write_text(json.dumps({"native_Berry_audits": endpoints}))
    monkeypatch.setattr(flow.native, "audit_scf", lambda *a: {})
    monkeypatch.setattr(flow.native, "audit_nscf", lambda root, index, n, *a:
                        records(.3-.1*(index%7))[flow.native.GRIDS.index(n)])
    for index, point in enumerate(points):
        folder = tmp_path/"calculations"/point["label"]
        folder.mkdir(parents=True)
        (folder/"scf_audit.json").write_text("{}")
        for record in records(.3-.1*(index%7)):
            (folder/f"berry_22{record['nz']}.json").write_text(json.dumps(record))
    result = flow.summarize(REPOSITORY, tmp_path)
    assert json.loads((tmp_path/"summary.json").read_text()) == result
    assert result["all_longitudinal_gates_passed"] is True
    assert result["all_sampled_lifts_unique"] is True
    assert not result["full_G1_certified"]


def test_metadata_writer_serializes_native_numpy_scalars_without_accepting_nonfinite(tmp_path):
    path = tmp_path/"scalar.json"
    flow.save(path, {"passed": np.bool_(True), "count": np.int64(14), "value": np.float32(.5)})
    assert json.loads(path.read_text()) == {"passed": True, "count": 14, "value": .5}
    for bad in (np.nan, np.float32(np.inf), np.array([1.])):
        with pytest.raises((ValueError, TypeError)):
            flow.save(path, {"invalid": bad})
