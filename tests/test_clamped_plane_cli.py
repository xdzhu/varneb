"""Explicit clamped-plane config/production entry; no electronic calculation."""

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.io import read, write

from vcneb import material_runner
from vcneb.cli import main as cli_main
from vcneb.config import RunConfig, prepare_run
from vcneb.epitaxial_boundary import clamped_plane_vcneb_boundary


def _case(tmp_path, allow_tilt=True):
    initial = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = initial.cell.array.copy()
    cell[0] += [0.0, 0.10, 0.15]
    cell[1] += [0.05, 0.0, 0.08]
    initial.set_cell(cell, scale_atoms=True)
    initial.positions[0, 2] += 0.03
    final = initial.copy()
    normal = np.cross(cell[0], cell[1])
    normal /= np.linalg.norm(normal)
    delta = .04 * (np.array([.3, -.2, 1.]) if allow_tilt else normal)
    final.cell[2] += delta
    final.positions[0, 2] += .06
    for name, atoms in (("initial.traj", initial), ("final.traj", final), ("plane.traj", initial)):
        write(tmp_path / name, atoms)
    return initial, final


def _args(tmp_path, allow_tilt=True):
    return [
        "--initial", str(tmp_path / "initial.traj"),
        "--final", str(tmp_path / "final.traj"),
        "--workdir", str(tmp_path / "run"), "--n-images", "3",
        "--clamped-plane-reference", str(tmp_path / "plane.traj"),
        "--clamped-allow-tilt", str(allow_tilt).lower(), "--cell-scale", "3.6",
        "--no-align-cells", "--cell-interpolation", "linear", "--mapping", "identity",
    ]


def _config(tmp_path, allow_tilt=True):
    payload = {
        "schema_version": 1, "backend": "ase", "initial": "initial.traj",
        "final": "final.traj", "workdir": "run", "n_images": 3, "mapping": "identity",
        "clamped_plane_reference": "plane.traj", "clamped_allow_tilt": allow_tilt,
        "cell_scale_A": 3.6, "align_cells": False, "cell_interpolation": "linear",
        "align_translation": False, "mic": False, "steps": 3,
        "fmax_ev_per_angstrom": 1e-8,
        "calculator": {"kind": "ase_class", "symbol": "ase.calculators.emt:EMT",
                       "parameters": {}},
    }
    path = tmp_path / "varneb.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _forbidden_load(_):
    raise AssertionError("preflight/rejection must not load a calculator")


@pytest.mark.parametrize("allow_tilt", [True, False])
def test_preflight_retains_raw_geometry_and_mechanical_provenance(tmp_path, allow_tilt):
    initial, final = _case(tmp_path, allow_tilt)
    material_runner.main(_args(tmp_path, allow_tilt) + ["--validate-only"],
                         symbol_loader=_forbidden_load)
    report = json.loads((tmp_path / "run/vcneb_preflight.json").read_text())
    record = report["mechanical_boundary"]
    assert record["reference_file_sha256"] == hashlib.sha256((tmp_path / "plane.traj").read_bytes()).hexdigest()
    assert record["allow_tilt"] is allow_tilt
    assert record["cell_dofs"] == (3 if allow_tilt else 1)
    assert record["raw_chain_policy"] == "reject_incompatible_no_projection"
    assert record["cell_scale_A"] == 3.6
    assert report["requires_stress"] is True
    assert report["initial_path_geometry"]["cell_scale_A"] == 3.6
    assert report["mode_subspace"] is None
    assert report["climbing_image_requested"] is False
    assert not (tmp_path / "run/initial-vcneb-unprojected.traj").exists()
    images = read(tmp_path / "run/initial-vcneb.traj", index=":")
    np.testing.assert_array_equal(images[0].cell, initial.cell)
    # Interpolation transforms coordinates through fractional space: roundoff,
    # not a physical geometry correction or a relaxed convergence tolerance.
    np.testing.assert_allclose(images[-1].positions, final.positions, atol=1e-14, rtol=0)
    clamped_plane_vcneb_boundary(4, initial.cell.array, allow_tilt=allow_tilt).validate_images(images)


@pytest.mark.parametrize("seed_flag", ["--initial-chain", "--resume-snapshot"])
@pytest.mark.parametrize("bad_index", [0, 1, 2])
def test_bad_raw_plane_is_rejected_before_projection_or_calculator(tmp_path, seed_flag, bad_index):
    initial, final = _case(tmp_path)
    middle = initial.copy()
    middle.positions += .01
    images = [initial, middle, final]
    images[bad_index].cell[0, 0] += .001
    # Match endpoint hashes so a mechanical rejection is tested, not identity.
    if bad_index != 1:
        write(tmp_path / ("initial.traj" if bad_index == 0 else "final.traj"), images[bad_index])
    raw = tmp_path / "raw.traj"
    write(raw, images)
    before = raw.read_bytes()
    with pytest.raises(ValueError, match="substrate vectors"):
        material_runner.main(_args(tmp_path) + [seed_flag, str(raw), "--validate-only"],
                             symbol_loader=_forbidden_load)
    assert raw.read_bytes() == before
    assert not (tmp_path / "run/initial-vcneb.traj").exists()


def test_normal_only_rejects_interior_tilt(tmp_path):
    initial, final = _case(tmp_path, False)
    middle = initial.copy()
    middle.cell[2] += .01 * initial.cell[0]
    raw = tmp_path / "raw.traj"
    write(raw, [initial, middle, final])
    with pytest.raises(ValueError, match="tilt is not allowed"):
        material_runner.main(_args(tmp_path, False) + ["--resume-snapshot", str(raw), "--validate-only"],
                             symbol_loader=_forbidden_load)


def test_existing_preflight_cannot_silently_change_boundary(tmp_path):
    _case(tmp_path, False)
    material_runner.main(_args(tmp_path, False) + ["--validate-only"], symbol_loader=_forbidden_load)
    saved = (tmp_path / "run/vcneb_preflight.json").read_bytes()
    with pytest.raises(FileExistsError, match="different mechanical boundary"):
        material_runner.main(_args(tmp_path, True) + ["--validate-only"], symbol_loader=_forbidden_load)
    assert (tmp_path / "run/vcneb_preflight.json").read_bytes() == saved


@pytest.mark.parametrize("clamped_first", [True, False])
def test_existing_directory_cannot_switch_between_clamped_and_free(tmp_path, clamped_first):
    _case(tmp_path, False)
    clamped = _args(tmp_path, False) + ["--validate-only"]
    free = clamped.copy()
    for option in ("--clamped-plane-reference", "--clamped-allow-tilt"):
        index = free.index(option)
        del free[index:index + 2]
    first, second = (clamped, free) if clamped_first else (free, clamped)
    material_runner.main(first, symbol_loader=_forbidden_load)
    saved = (tmp_path / "run/vcneb_preflight.json").read_bytes()
    with pytest.raises(FileExistsError, match="different mechanical boundary"):
        material_runner.main(second, symbol_loader=_forbidden_load)
    assert (tmp_path / "run/vcneb_preflight.json").read_bytes() == saved


def test_extra_atomic_constraint_is_not_silently_ignored(tmp_path):
    initial, _ = _case(tmp_path)
    initial.set_constraint(FixAtoms(indices=[0]))
    write(tmp_path / "initial.traj", initial)
    with pytest.raises(ValueError, match="extra ASE constraints"):
        material_runner.main(_args(tmp_path) + ["--validate-only"], symbol_loader=_forbidden_load)


@pytest.mark.parametrize("missing,message", [
    ("--clamped-allow-tilt", "explicit boolean"),
    ("--cell-scale", "explicit cell_scale_A"),
    ("--no-align-cells", "align_cells false"),
    ("--clamped-plane-reference", "requires clamped_plane_reference"),
])
def test_partial_cli_boundary_is_rejected_before_file_read(tmp_path, missing, message):
    args = _args(tmp_path)
    index = args.index(missing)
    del args[index:index + (1 if missing == "--no-align-cells" else 2)]
    with pytest.raises(ValueError, match=message):
        material_runner.main(args + ["--validate-only"], symbol_loader=_forbidden_load)
    assert not (tmp_path / "run").exists()


@pytest.mark.parametrize("extra,message", [
    (["--cell-mode", "fixed"], "cell_mode full"),
    (["--subspace-artifact", "unread.json"], "cannot be combined"),
    (["--cell-scale", "nan"], "finite and positive"),
    (["--cell-interpolation", "log_strain"], "requires --cell-interpolation linear"),
])
def test_unsupported_compositions_fail_explicitly(tmp_path, extra, message):
    with pytest.raises(ValueError, match=message):
        material_runner.main(_args(tmp_path) + extra + ["--validate-only"], symbol_loader=_forbidden_load)


@pytest.mark.parametrize("allow_tilt", [True, False])
@pytest.mark.parametrize("workers", [0, 2])
def test_actual_few_step_execution_keeps_plane_and_endpoint_cache(tmp_path, allow_tilt, workers):
    initial, _ = _case(tmp_path, allow_tilt)
    calculators, backend_options = {}, {}

    class CountingEMT(EMT):
        def calculate(self, atoms=None, properties=("energy",), system_changes=None):
            self.calls += 1
            self.geometries.append(atoms.copy())
            super().calculate(atoms, properties, system_changes)

    def builder(**kwargs):
        backend_options.update(kwargs)
        def factory(index, image, directory):
            calc = CountingEMT(directory=str(directory))
            calc.calls, calc.geometries = 0, []
            calculators[index] = calc
            return calc
        return factory

    material_runner.main(_args(tmp_path, allow_tilt) + [
        "--factory", "toy:builder", "--parameters-json", '{"unchanged_contract":100}',
        "--factory-kwargs-json", '{"unchanged_basis":"10au_DZP"}',
        "--command", "unchanged-launcher", "--steps", "3", "--fmax", "1e-8",
        "--maxstep", ".01", "--maximum-cell-step", ".02", "--candidate-step-retries", "2",
        "--minimum-distance", "1.5", "--image-workers", str(workers), "--summary-format", "brief",
    ], symbol_loader=lambda _: builder)
    assert backend_options == {"parameters": {"unchanged_contract": 100},
                               "unchanged_basis": "10au_DZP", "command": "unchanged-launcher"}
    assert calculators[0].calls == calculators[2].calls == 1
    assert calculators[1].calls >= 3
    boundary = clamped_plane_vcneb_boundary(4, initial.cell.array, allow_tilt=allow_tilt)
    for calc in calculators.values():
        boundary.validate_images(calc.geometries)
    snapshots = read(tmp_path / "run/vcneb.traj", index=":")
    boundary.validate_images(snapshots)
    assert len(snapshots) == 12  # initial plus 3 steps, each with 3 images
    assert not np.array_equal(snapshots[1].positions, snapshots[-2].positions)
    summary = json.loads((tmp_path / "run/vcneb_summary.json").read_text())
    assert summary["requires_stress"] is True
    assert summary["endpoint_evaluation_policy"] == "fixed_cached_once"
    assert summary["cell_scale_A"] == 3.6
    assert not summary["climbing_image_active_final"]


@pytest.mark.parametrize("allow_tilt", [True, False])
def test_public_config_prepare_and_run_use_same_plane(tmp_path, monkeypatch, allow_tilt):
    initial, _ = _case(tmp_path, allow_tilt)
    path = _config(tmp_path, allow_tilt)
    config, report_path = prepare_run(path)
    assert config.clamped_plane_reference == (tmp_path / "plane.traj").resolve()
    report = json.loads(report_path.read_text())
    assert report["mechanical_boundary"]["allow_tilt"] is allow_tilt
    assert report["initial_path_geometry"]["cell_scale_A"] == 3.6
    saved = (tmp_path / "run/initial-vcneb.traj").read_bytes()
    assert prepare_run(path)[1] == report_path
    monkeypatch.setattr(sys, "argv", ["varneb", "run", str(path), "--execute"])
    assert cli_main() == 0
    assert (tmp_path / "run/initial-vcneb.traj").read_bytes() == saved
    boundary = clamped_plane_vcneb_boundary(4, initial.cell.array, allow_tilt=allow_tilt)
    boundary.validate_images(read(tmp_path / "run/vcneb.traj", index=":"))


@pytest.mark.parametrize("field,value,message", [
    ("clamped_allow_tilt", "false", "must be a boolean"),
    ("clamped_allow_tilt", None, "explicit boolean"),
    ("cell_scale_A", None, "explicit cell_scale_A"),
    ("cell_scale_A", 0., "finite and positive"),
    ("align_cells", True, "align_cells false"),
    ("align_cells", "false", "must be a boolean"),
    ("cell_mode", "fixed", "cell_mode full"),
    ("cell_interpolation", "log_strain", "cell_interpolation linear"),
])
def test_config_does_not_coerce_or_infer_boundary(tmp_path, field, value, message):
    path = _config(tmp_path)
    payload = json.loads(path.read_text())
    payload[field] = value
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        RunConfig.from_file(path)
    assert not (tmp_path / "run").exists()


def test_public_prepare_rejects_wrong_plane_before_creating_workdir(tmp_path):
    _, final = _case(tmp_path)
    final.cell[1, 1] += .001
    write(tmp_path / "final.traj", final)
    with pytest.raises(ValueError, match="substrate vectors"):
        prepare_run(_config(tmp_path))
    assert not (tmp_path / "run").exists()


def test_actual_hafnia_seeds_pass_but_released_endpoints_do_not(tmp_path):
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"
    seeds = root / "clamped_endpoint_seeds"
    manifest = json.loads((seeds / "clamped_seed_manifest.json").read_text())
    assert len(manifest["seeds"]) == 10
    # Reuse the actual manifest, without turning these unrelaxed seeds into DFT evidence.
    for index, entry in enumerate(manifest["seeds"]):
        endpoint_manifest = seeds / entry["manifest_file"]
        details = json.loads(endpoint_manifest.read_text())
        final = endpoint_manifest.parent / details["seed_file"]
        reference = endpoint_manifest.parent.parent / "T/POSCAR.seed"
        assert hashlib.sha256(final.read_bytes()).hexdigest() == entry["seed_sha256"]
        args = ["--initial", str(reference), "--final", str(final),
                "--clamped-plane-reference", str(reference), "--clamped-allow-tilt", "true",
                "--cell-scale", str(manifest["cell_scale_A"]),
                "--workdir", str(tmp_path / f"seed_{index}"), "--n-images", "3",
                "--no-align-cells", "--cell-interpolation", "linear", "--mapping", "identity",
                "--validate-only"]
        material_runner.main(args, symbol_loader=_forbidden_load)
        assert len(read(tmp_path / f"seed_{index}/initial-vcneb.traj", index=":")) == 3
    with pytest.raises(ValueError, match="substrate vectors"):
        material_runner.main([
            "--initial", str(root / "reference_variants/T.vasp"),
            "--final", str(root / "reference_variants/PO.vasp"),
            "--clamped-plane-reference", str(seeds / "strain_0000/T/POSCAR.seed"),
            "--clamped-allow-tilt", "true", "--cell-scale", str(manifest["cell_scale_A"]),
            "--workdir", str(tmp_path / "released_rejected"),
            "--no-align-cells", "--cell-interpolation", "linear", "--validate-only",
        ], symbol_loader=_forbidden_load)
    assert not (tmp_path / "released_rejected").exists()
