import json
from pathlib import Path
import shutil

import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read

from scripts.analyze_hfo2_clamped_residual import analyze, decompose
from scripts.export_hfo2_clamped_observation import CELL_SCALE
from scripts.audit_hfo2_static_replica import audited_results, sha256
from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from vcneb import VCNEB, clamped_plane_vcneb_boundary

ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT/"benchmarks/hfo2_channels/20261008/clamped_G2_E053_20261010/strain_0000/PO_to_M"


def synthetic_chain():
    images = read(SEED/"seed.traj", index=":")
    for i, image in enumerate(images):
        force = np.zeros((12, 3)); force[0, 2] = .03*i
        image.calc = SinglePointCalculator(image, energy=-10+np.sin(i*np.pi/8), forces=force,
                                          stress=np.array([.1, .1, .02, .003, .004, 0]))
    boundary = clamped_plane_vcneb_boundary(12, read(SEED/"substrate.vasp").cell.array, allow_tilt=True)
    return VCNEB(images, cell_scale=CELL_SCALE, k=.2, **boundary.vcneb_kwargs(images))


def check_sums(result):
    for row in result["images"]:
        p, s, r = (np.array(row[n+"_generalized_vector_eV_A"]) for n in ("perpendicular", "spring", "residual"))
        np.testing.assert_allclose(p+s, r, atol=1e-12, rtol=0)
        assert row["perpendicular_spring_dot_eV2_A2"] == pytest.approx(0, abs=1e-10)
        assert row["perpendicular_projection_on_worst_residual_eV_A"]+row["spring_projection_on_worst_residual_eV_A"] == pytest.approx(row["residual_max_vector_eV_A"], abs=1e-12)
        assert max(row["residual_atom_max_vector_eV_A"], row["residual_open_cell_max_vector_eV_A"]) == pytest.approx(row["residual_max_vector_eV_A"], abs=1e-12)


def test_exact_projected_vector_identity_and_blocks():
    chain = synthetic_chain()
    result = decompose(chain)
    check_sums(result)
    assert result["fmax_eV_A"] == pytest.approx(np.linalg.norm(chain.get_forces(), axis=1).max(), abs=1e-12)
    assert chain.dynamic_relaxation == 1 and not chain.climb


@pytest.mark.parametrize("option", ["climb", "weights"])
def test_nonordinary_or_weighted_snapshots_rejected(option):
    chain = synthetic_chain()
    if option == "climb":
        chain.climb = True
    else:
        chain.dynamic_relaxation = .5
    with pytest.raises(ValueError, match="unweighted ordinary"):
        decompose(chain)


@pytest.mark.parametrize("step", [0, 20])
def test_previous_native_M_frame_matches_portable_replay(step):
    directory = ROOT/f"benchmarks/hfo2_channels/20261008/clamped_M_continue_E062_20261011/observations/step_{step:04d}"
    result = analyze(directory)
    check_sums(result)
    assert result["new_DFT_calls"] == 0 and not result["ordinary_residual_pass"]


@pytest.mark.parametrize("step", [13, 14, 15])
def test_actual_flip_frame_independently_reproduces_HF_components(step):
    directory = ROOT/f"benchmarks/hfo2_channels/20261008/clamped_flip_residual_E063_20261011/observations/step_{step:04d}"
    result = analyze(directory)
    native = json.loads((directory/"residual.json").read_text())
    check_sums(result)
    assert result.keys() == native.keys()
    for left, right in zip(result["images"], native["images"]):
        assert left.keys() == right.keys()
        for key in left:
            if isinstance(left[key], (int, float, list)):
                np.testing.assert_allclose(left[key], right[key], atol=1e-10, rtol=0)
            else:
                assert left[key] == right[key]
    assert result["fmax_eV_A"] == pytest.approx(native["fmax_eV_A"], abs=1e-10)
    observation = json.loads((directory/"observation.json").read_text())
    images = read(directory/"evaluated_chain.traj", index=":")
    for i, image in enumerate(images):
        raw_directory = directory/f"raw/image_{i:04d}"
        pinned = observation["raw_image_evaluations"][i]
        assert all(pinned["input_sha256"][name] == value for name, value in CONTRACT.items())
        # Six physical bytes were checked on HF; only nonproprietary inputs
        # and native observables are present in the portable export.
        for name in ("INPUT", "KPT", "STRU"):
            assert sha256(raw_directory/name) == pinned["input_sha256"][name]
        assert sha256(raw_directory/"OUT.ABACUS/running_scf.log") == pinned["raw_log_sha256"]
        assert sha256(raw_directory/"source_audit.json") == pinned["audit_sha256"]
        assert same_ordered_geometry(image, read_fixed_hfo2_stru(raw_directory/"STRU"))
        parsed = audited_results(raw_directory)
        for key in ("energy", "forces", "stress"):
            np.testing.assert_allclose(parsed[key], image.calc.results[key], atol=1e-12, rtol=0)


@pytest.mark.parametrize("fault", ["trajectory", "native_force", "metric", "boundary", "weighted_contract"])
def test_changed_observation_refused(tmp_path, fault):
    source = ROOT/"benchmarks/hfo2_channels/20261008/clamped_M_continue_E062_20261011/observations/step_0020"
    directory = tmp_path/"observation"
    directory.mkdir()
    for name in ("observation.json", "evaluated_chain.traj"):
        shutil.copyfile(source/name, directory/name)
    path = directory/"observation.json"
    record = json.loads(path.read_text())
    if fault == "trajectory":
        record["evaluated_chain_sha256"] = "0"*64
    elif fault == "native_force":
        record["raw_image_evaluations"][4]["forces_eV_A"][0][0] += .1
    elif fault == "metric":
        record["cell_scale_A"] += .1
    elif fault == "boundary":
        record["mechanical_boundary"]["allow_tilt"] = False
    else:
        record["climb"] = True
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError):
        analyze(directory)
