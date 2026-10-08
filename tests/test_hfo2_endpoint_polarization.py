import json
from ase import Atoms
import numpy as np
import pytest

from scripts import hfo2_endpoint_polarization as flow
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_channel_work_probes import write_structure


BASE = b"INPUT_PARAMETERS\ncalculation scf\nbasis_type lcao\ndft_functional pbe\necutwfc 100.0\ncal_force 1\ncal_stress 1\nout_stru 1\nscf_thr 1e-08\nntype 2\n"


def test_nscf_explicit_delta_and_duplicate_rejection():
    old = flow.input_values(BASE.decode())
    new = flow.input_values(flow.berry_input(BASE).decode())
    changes = {k: v for k, v in new.items() if old.get(k) != v}
    assert changes == {"calculation": "nscf", "cal_force": "0", "cal_stress": "0", "out_stru": "0",
                       "init_chg": "file", "symmetry": "-1", "berry_phase": "1", "gdir": "3", "out_band": "1"}
    assert new["ecutwfc"] == old["ecutwfc"] == "100.0"
    with pytest.raises(ValueError, match="duplicate"):
        flow.input_values("INPUT_PARAMETERS\necutwfc 100\necutwfc 80\n")


def test_prepare_preserves_sources_and_output_only_scf(tmp_path, monkeypatch):
    source = tmp_path / "baseline"
    source.mkdir()
    (source / "INPUT").write_bytes(BASE)
    (source / "KPT").write_bytes(b"K_POINTS\n0\nGamma\n2 2 2 0 0 0\n")
    for symbol, valence in (("Hf", 12), ("O", 6)):
        (source / f"{symbol}.upf").write_text(f'<PP_HEADER z_valence=" {valence}"/>')
    for name in flow.INPUT_FILES:
        if name not in ("STRU", "INPUT", "KPT", "Hf.upf", "O.upf"):
            (source / name).write_text(name)
    atoms = Atoms(["Hf"]*4 + ["O"]*8, positions=np.arange(36).reshape(12, 3)*.1, cell=np.eye(3)*5, pbc=True)
    write_structure(source / "STRU", atoms)
    # Upstream local ASE has no ABACUS IO. Real unmocked read is a separate
    # preflight on hf's ABACUS-enabled ASE, not a claim from this fixture.
    monkeypatch.setattr(flow, "read", lambda path, format: atoms.copy())
    (source / "OUT.ABACUS").mkdir()
    (source / "OUT.ABACUS/running_scf.log").write_text("fixture")
    contract = {n: sha256(source / n) for n in flow.CONTRACT}
    monkeypatch.setattr(flow, "CONTRACT", contract)
    monkeypatch.setattr(flow, "audited_results", lambda _: {"energy": -10., "forces": np.zeros((12, 3)), "stress": np.zeros(6)})
    sources = {label: str(source) for label in flow.LABELS}
    root = tmp_path / "fresh"
    report = flow.prepare(sources, root)
    for point in report["points"]:
        scf = root / "points" / point["label"] / "scf"
        assert (scf / "INPUT").read_bytes() == BASE + flow.OUTPUT_ADDITION
        for name in flow.INPUT_FILES[1:]:
            assert sha256(scf / name) == sha256(source / name)
    assert (source / "INPUT").read_bytes() == BASE
    assert not report["SCF_physical_settings_changed"] and not report["NSCF_energies_used_for_barriers"]
    with pytest.raises(FileExistsError):
        flow.prepare(sources, root)
    (source / "INPUT").write_bytes(BASE.replace(b"100.0", b"80.0"))
    with pytest.raises(ValueError, match="contract changed"):
        flow.prepare(sources, tmp_path / "wrong")


def test_no_dft_outside_explicit_allocation(tmp_path, monkeypatch):
    monkeypatch.delenv("RUN_DFT", raising=False)
    with pytest.raises(RuntimeError, match="allocation"):
        flow.run_endpoint(tmp_path, 0)


def test_summary_preserves_branches_and_rejects_self_inverse(tmp_path):
    (tmp_path / "manifest.json").write_text(json.dumps({"tolerances": {
        "longitudinal_P_convergence_C_m2": .01, "inversion_modular_residual_C_m2": .01}, "limitations": "fixture"}))
    for label, p in zip(flow.LABELS, (.5, -.5, .7)):
        folder = tmp_path / "calculations" / label
        folder.mkdir(parents=True)
        (folder / "scf_audit.json").write_text("{}")
        for nz in flow.GRIDS:
            (folder / f"berry_22{nz}.json").write_text(json.dumps({"value_C_m2": p, "reported_modulus_C_m2": 1.2}))
    summary = flow.summarize(tmp_path)
    assert summary["status"] == "passed"
    assert summary["spontaneous_polarization_C_m2"] is None
    assert not summary["switching_path_branch_selected"]
    for label in flow.LABELS:
        for nz in flow.GRIDS:
            (tmp_path / "calculations" / label / f"berry_22{nz}.json").write_text(json.dumps({"value_C_m2": .6, "reported_modulus_C_m2": 1.2}))
    assert flow.summarize(tmp_path)["status"] == "review_required"


def test_actual_audit_contract_uses_different_scf_and_nscf_eigenvalue_files(tmp_path, monkeypatch):
    work = tmp_path / "work"
    work.mkdir()
    (work / "INPUT").write_text("fixture input")
    hashes = {"INPUT": sha256(work / "INPUT")}
    out = work / "OUT.ABACUS"
    out.mkdir()
    (out / "SPIN1_CHG.cube").write_text("immutable charge")
    (out / "running_scf.log").write_text("fixture SCF")
    # Eight kpoints, 48 occupied bands and at least one empty band.
    energies = list(np.arange(48)/48 - 1) + [3.]
    istate = "".join("BAND Energy(ev) Occupation Kpoint = %d\n" % (k+1) +
                     "".join(f"{i+1} {e} .25\n" for i, e in enumerate(energies)) for k in range(8))
    (out / "istate.info").write_text(istate)
    point = {"label": "PO_plus", "scf_input_sha256": hashes, "occupied_bands": 48,
             "baseline_results": {"energy": -10., "forces": np.zeros((12, 3)).tolist(), "stress": np.zeros(6).tolist()},
             "nscf_inputs": [{"nz": 2, "input_sha256": hashes}], "cell_A": (np.eye(3)*5).tolist(),
             "quantum_lattice_C_m2": flow.quantum_lattice(np.eye(3)*5).tolist()}
    manifest = {"points": [point], "tolerances": {"SCF_energy_eV_cell": 1e-5, "SCF_force_eV_A": 1e-4,
        "SCF_stress_kbar": .02, "sampled_gap_min_eV": .1}}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    monkeypatch.setattr(flow, "audited_results", lambda _: {"energy": -10., "forces": np.zeros((12, 3)), "stress": np.zeros(6)})
    audit = flow.audit_scf(tmp_path, 0, work)
    assert not (out / "BANDS_1.dat").exists()
    (out / "istate.info").unlink()  # NSCF does not write this file in actual f7cb1d3.
    (out / "BANDS_1.dat").write_text("".join(f"{k+1} 0 " + " ".join(map(str, energies)) + "\n" for k in range(8)))
    modulus = np.linalg.norm(np.array(point["quantum_lattice_C_m2"])[2])*2
    (out / "running_nscf.log").write_text(f"DSIZE = 32\nThe calculated polarization direction is in R3 direction\n"
        f"P = 0.5 (mod {modulus:.7f}) (0.0, 0.0, 0.5) C/m^2\n Total Time : 28\n")
    result = flow.audit_nscf(tmp_path, 0, 2, work, audit)
    assert result["bands"]["source_format"] == "ABACUS_BANDS_1.dat"
    (out / "SPIN1_CHG.cube").write_text("changed charge")
    with pytest.raises(ValueError, match="charge changed"):
        flow.audit_nscf(tmp_path, 0, 2, work, audit)
