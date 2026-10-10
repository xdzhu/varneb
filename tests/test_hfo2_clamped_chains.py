"""Registered material polygons and negative cache/geometry cases; no DFT."""
import json
from pathlib import Path
import shutil

import numpy as np
import pytest
from ase.io import read

import examples.hfo2_fixed_input_factory as transport
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_clamped_chains import CHANNEL_FINAL, deform_polygon, endpoint_evidence, prepare
from vcneb.periodic_path import validate_periodic_path_lift

REPOSITORY = Path(__file__).resolve().parents[1]
CASE = REPOSITORY / "benchmarks/hfo2_channels/20261008"


@pytest.mark.parametrize("condition,channel", [("strain_0000",c) for c in CHANNEL_FINAL]
                         + [("strain_p0100",c) for c in CHANNEL_FINAL])
def test_actual_audited_chain_preparation_has_two_exact_clamped_wells_and_no_free_EFS(tmp_path,condition,channel):
    output = tmp_path / "seed"
    m = prepare(CASE,condition,channel,output)
    assert m["new_DFT_calls"] == 0 and not m["holdout_generated"]
    assert m["n_total_images"] == 9 and m["n_internal_images"] == 7 and not m["climb"]
    assert m["cached_clamped_endpoint_evaluations"] == 2 and m["cached_interior_evaluations"] == 0
    assert not m["lift_preparation"]["atom_permutation_applied"]
    images = read(output/"seed.traj", index=":")
    assert len(images) == 9 and all(a.calc is None for a in images)
    validate_periodic_path_lift(images)
    for a in images:
        np.testing.assert_allclose(a.cell.array[:2],np.array(m["reference_cell_A"])[:2],atol=1e-12,rtol=0)
    for a,phase in zip((images[0],images[-1]),("PO_plus",CHANNEL_FINAL[channel])):
        assert transport.same_ordered_geometry(a,endpoint_evidence(CASE,condition,phase)[1])
    p = json.loads((output/"factory_parameters.json").read_text())
    assert p["seed_cache_records"][1:-1] == [None]*7
    assert "clamped" in p["seed_cache_records"][0]["directory"]
    with pytest.raises(FileExistsError):
        prepare(CASE,condition,channel,output)


def test_unregistered_holdout_or_channel_is_not_generated(tmp_path):
    for condition,channel in (("strain_p0050","PO_to_M"),("strain_0000","new_phase")):
        with pytest.raises(ValueError,match="registered"):
            prepare(CASE,condition,channel,tmp_path/"refused")
        assert not (tmp_path/"refused").exists()


def test_preparation_refuses_ordered_endpoint_remapping(tmp_path):
    seed,final,boundary,_,_ = endpoint_evidence(CASE,"strain_0000","PO_plus")
    # The transform expects old axes; use its proper inverse to build a fixture.
    from scripts.prepare_hfo2_clamped_endpoints import orient_long_axis_x
    original = orient_long_axis_x(orient_long_axis_x(seed))
    changed = seed.copy()
    changed.positions[[0,1]] = changed.positions[[1,0]]
    with pytest.raises(ValueError,match="ordered endpoint"):
        deform_polygon([original,original],(changed,seed),(final,final),boundary)


def cache_fixture(tmp_path, monkeypatch):
    seed,final,_,_,_ = endpoint_evidence(CASE,"strain_0000","M")
    root = CASE/"clamped_endpoint_matrix_20261009/strain_0000/M/completed_HF/endpoint/calculator/image_0000/scf_000017"
    source = tmp_path/"raw"
    source.mkdir()
    for n in ("INPUT","KPT","STRU","call_audit.json"):
        shutil.copyfile(root/n,source/n)
    (source/"OUT.ABACUS").mkdir()
    shutil.copyfile(root/"OUT.ABACUS/running_scf.log",source/"OUT.ABACUS/running_scf.log")
    monkeypatch.setattr(transport,"CONTRACT",{n:sha256(source/n) for n in ("INPUT","KPT")})
    record = {"directory":str(source),"input_sha256":{n:sha256(source/n) for n in (*transport.CONTRACT,"STRU")},
              "raw_log_sha256":sha256(source/"OUT.ABACUS/running_scf.log")}
    p = {"source_directory":str(source),"seed_cache_records":[record]+[None]*7+[record]}
    monkeypatch.setattr(transport.subprocess,"run",lambda *a,**k:pytest.fail("cache preflight must not launch DFT"))
    return p,final


def test_new_clamped_cache_uses_narrow_reader_and_only_two_wells(tmp_path,monkeypatch):
    p,atoms = cache_fixture(tmp_path,monkeypatch)
    monkeypatch.setattr(transport,"read",lambda *a,**k:pytest.fail("optional ABACUS I/O must not be needed"))
    factory = transport.make_clamped_seed_cached_factory(parameters=p,command="mpirun -np 32 abacus")
    for i in (0,8):
        a = atoms.copy()
        a.positions[0] += a.cell[0]  # pinned integer lift, not a different geometry
        a.calc = factory(i,a,tmp_path/f"calc{i}")
        assert a.get_potential_energy() == pytest.approx(-9783.343992919006,abs=1e-10,rel=0)
        assert a.get_forces().shape == (12,3) and a.get_stress().shape == (6,)
        assert a.calc.next_call == 0
    assert factory(1,atoms,tmp_path/"fresh").results == {}


@pytest.mark.parametrize("mutation",["log","STRU","geometry","interior","index","existing","override"])
def test_clamped_cache_rejects_mutations(tmp_path,monkeypatch,mutation):
    p,a = cache_fixture(tmp_path,monkeypatch)
    target = tmp_path/"calc"
    if mutation == "interior":
        p["seed_cache_records"][1] = p["seed_cache_records"][0]
    elif mutation == "override":
        p["ecutwfc"] = 120
    elif mutation in ("log","STRU"):
        path = Path(p["seed_cache_records"][0]["directory"])/("STRU" if mutation == "STRU" else "OUT.ABACUS/running_scf.log")
        path.write_text(path.read_text()+"changed")
    elif mutation == "geometry":
        a.positions[0,0] += .001
    elif mutation == "existing":
        target.mkdir()
        (target/"keep.json").write_text("preserve")
    with pytest.raises(ValueError):
        factory = transport.make_clamped_seed_cached_factory(parameters=p,command="mpirun -np 32 abacus")
        factory(9 if mutation == "index" else 0,a,target)


def test_bounded_clamped_template_does_not_use_free_cell_or_old_cache_contract():
    text = (REPOSITORY/"cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm").read_text()
    for required in ('--clamped-plane-reference','--clamped-allow-tilt true','--cell-scale',
                     '--n-images 9','--steps 10','--fmax .10','--no-climb','make_clamped_seed_cached_factory',
                     '#SBATCH --time=04:00:00','mpirun -np 32','--require-continuous-periodic-lift'):
        assert required in text
    assert 'sbatch' not in text and 'srun' not in text and 'make_seed_cached_factory' not in text


def test_clamped_resume_reuses_all_nine_exact_SCFS_then_invalidates_moved_interior(tmp_path,monkeypatch):
    p,atoms = cache_fixture(tmp_path,monkeypatch)
    p["seed_cache_records"] = [p["seed_cache_records"][0]]*9
    factory = transport.make_clamped_resume_cached_factory(parameters=p,command="mpirun -np 32 abacus")
    for i in range(9):
        a = atoms.copy()
        a.calc = factory(i,a,tmp_path/f"resume{i}")
        assert np.isfinite(a.get_potential_energy()) and a.get_forces().shape==(12,3)
        assert a.calc.next_call==0
        receipt = json.loads((tmp_path/f"resume{i}/seed_cache_audit.json").read_text())
        assert receipt["policy"]=="identical_ordered_clamped_resume_hash_pinned"
        assert receipt["new_DFT_calls"]==0
        if i==4:
            a.positions[0,0] += .001
            assert a.calc.check_state(a)==["positions"]


def test_clamped_resume_requires_every_geometry_pinned(tmp_path,monkeypatch):
    p,_ = cache_fixture(tmp_path,monkeypatch)
    with pytest.raises(ValueError,match="nine exact caches"):
        transport.make_clamped_resume_cached_factory(parameters=p,command="mpirun -np 32 abacus")
