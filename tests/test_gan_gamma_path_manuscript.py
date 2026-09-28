"""Guards for the all-image, same-protocol GaN Gamma figure input join."""

from __future__ import annotations

import json

import pytest

from scripts import plot_gan_gamma_path_manuscript as figure


def fixture_inputs(tmp_path, monkeypatch):
    modes = tmp_path / "modes"
    modes.mkdir()
    chain = tmp_path / "chain.traj"
    original = modes / "gan_vasp_tetragonal_final_29_images.traj"
    chain.write_bytes(b"original-600-eV-29-image-chain")
    original.write_bytes(chain.read_bytes())
    electronic = {
        "INCAR": "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772",
        "KPOINTS": "b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee",
        "POTCAR": "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168",
    }
    manifest = {
        "pressure_GPa": 45.7,
        "supercell_matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "source": {phase: {"sha256": dict(electronic)} for phase in ("B4", "B1")},
    }
    for phase in ("B4", "B1"):
        endpoint = modes / f"{phase}_CONTCAR"
        endpoint.write_text(f"{phase} endpoint", encoding="utf-8")
        manifest["source"][phase]["sha256"]["CONTCAR"] = figure.sha256(endpoint)
    (modes / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    report = {"status": "audited_GaN_endpoint_Gamma_atomic_subspaces_not_TS_modes"}
    (modes / "subspace_audit_v2.json").write_text(json.dumps(report), encoding="utf-8")
    (modes / "audit.json").write_text(json.dumps({
        "source_manifest_sha256": figure.sha256(modes / "manifest.json"),
    }), encoding="utf-8")
    rows = []
    for phase in ("B4", "B1"):
        for index in range(29):
            rows.append({
                "reference_phase": phase,
                "image_index": index,
                "reaction_coordinate_normalized": index / 28,
                "relative_enthalpy_eV_per_GaN": 0.338 if index == 15 else 0.0,
                **{f"group{rank}_Q_norm_d0p01_sqrt_amu_A": rank * 0.1
                   for rank in (1, 2, 3)},
                "three_group_atomic_residual_sqrt_amu_A": 0.001,
                **{f"symmetric_strain_{axis}": 0.01 for axis in ("xx", "yy", "zz")},
            })
    monkeypatch.setattr(figure, "audit", lambda _directory, _path_csv: (report, rows))
    return modes, chain, rows


def test_same_chain_and_contract_join_29_rows(tmp_path, monkeypatch):
    modes, chain, _ = fixture_inputs(tmp_path, monkeypatch)
    _, joined = figure.assemble(modes, tmp_path / "unused.csv", chain)
    assert len(joined) == 29
    assert joined[15]["relative_enthalpy_eV_per_GaN"] == 0.338
    assert joined[28]["B4_group3_Q_sqrt_amu_A"] == pytest.approx(0.3)


def test_rejects_different_final_chain(tmp_path, monkeypatch):
    modes, chain, _ = fixture_inputs(tmp_path, monkeypatch)
    chain.write_bytes(b"different-chain")
    with pytest.raises(ValueError, match="not the frozen 600-eV final chain"):
        figure.assemble(modes, tmp_path / "unused.csv", chain)


def test_rejects_gamma_path_phase_mismatch(tmp_path, monkeypatch):
    modes, chain, rows = fixture_inputs(tmp_path, monkeypatch)
    rows[29 + 8]["reaction_coordinate_normalized"] += 0.01
    with pytest.raises(ValueError, match="not paired"):
        figure.assemble(modes, tmp_path / "unused.csv", chain)


def test_rejects_changed_electronic_contract(tmp_path, monkeypatch):
    modes, chain, _ = fixture_inputs(tmp_path, monkeypatch)
    path = modes / "manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["source"]["B1"]["sha256"]["INCAR"] = "other"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="600-eV production contract"):
        figure.assemble(modes, tmp_path / "unused.csv", chain)
