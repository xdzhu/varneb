import json
from pathlib import Path

import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

import scripts.export_hfo2_clamped_observation as observer
from scripts.prepare_hfo2_channel_work_probes import write_structure
from vcneb import VCNEB, clamped_plane_vcneb_boundary


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT/"benchmarks/hfo2_channels/20261008/clamped_G2_E053_20261010/strain_0000/PO_to_M"


@pytest.mark.parametrize("channel",["PO_to_M","PO_flip_T_pattern_preserving"])
def test_archived_prepared_G2_artifacts_preserve_registered_raw_bytes(channel):
    seed = SEED.parent/channel
    manifest = json.loads((seed/"manifest.json").read_text())
    assert all(observer.sha256(seed/n)==h for n,h in manifest["files_sha256"].items())
    preflight = json.loads((seed/"vcneb_preflight_HF.json").read_text())
    assert observer.sha256(seed/"substrate.vasp")==preflight["mechanical_boundary"]["reference_file_sha256"]


@pytest.fixture
def evaluated(tmp_path, monkeypatch):
    work = tmp_path/"live"
    snapshot = work/"snapshots/step_0010"
    snapshot.mkdir(parents=True)
    images = read(SEED/"seed.traj", index=":")
    raw_by_source = {}
    caches = [None]*9
    contract = {"INPUT": None}
    for i, image in enumerate(images):
        source = work/f"image_{i:04d}/scf_000000"
        (source/"OUT.ABACUS").mkdir(parents=True)
        (source/"INPUT").write_text("fixed physical input")
        contract["INPUT"] = observer.sha256(source/"INPUT")
        (source/"OUT.ABACUS/running_scf.log").write_text("synthetic raw SCF")
        write_structure(source/"STRU", image)
        raw = {"energy": -10 + .1*np.sin(i*np.pi/8), "forces": np.zeros((12,3)),
               "stress": np.array([.1,.1,0,0,0,0])}
        raw["forces"][0,2] = .3  # intentionally nonconverged synthetic transverse force
        image.calc = SinglePointCalculator(image, **raw)
        raw_by_source[source] = raw
        write(snapshot/f"POSCAR_{i:02d}", image, format="vasp", direct=True)
        pinned = {"input_sha256": {n:observer.sha256(source/n) for n in ("INPUT","STRU")},
                  "raw_log_sha256": observer.sha256(source/"OUT.ABACUS/running_scf.log")}
        if i in (0,8):
            caches[i] = {"directory":str(source), **pinned}
            audit = {"policy":"identical_ordered_clamped_endpoint_hash_pinned",
                     "raw_source":str(source), **pinned}
            (source.parent/"seed_cache_audit.json").write_text(json.dumps(audit))
        else:
            (source/"call_audit.json").write_text(json.dumps({**pinned,
                "results":{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in raw.items()}}))
    monkeypatch.setattr(observer,"CONTRACT",contract)
    monkeypatch.setattr(observer,"audited_results",lambda p: raw_by_source[p])
    substrate = tmp_path/"substrate.vasp"
    substrate.write_bytes((SEED/"substrate.vasp").read_bytes())
    metadata = json.loads((SEED/"vcneb_preflight_HF.json").read_text())
    metadata["calculator_validation"] = "runtime_instantiated"
    metadata["factory"] = "examples.hfo2_fixed_input_factory:make_clamped_seed_cached_factory"
    metadata["calculator_parameters"] = {"source_directory":"unused-read-only", "seed_cache_records":caches}
    metadata["mechanical_boundary"]["reference_file"] = str(substrate)
    (work/"vcneb_preflight.json").write_text(json.dumps(metadata))
    boundary = clamped_plane_vcneb_boundary(12,read(substrate).cell.array,allow_tilt=True)
    chain = VCNEB(images,cell_scale=observer.CELL_SCALE,**boundary.vcneb_kwargs(images))
    force = np.linalg.norm(chain.get_forces(),axis=1).max()
    (work/"vcneb.opt.log").write_text(f"CheckedFIRE: 10 12:00:00 {chain.enthalpies.max():.6f} {force:.6f}\n")
    return work, images


def run(work, output):
    return observer.export(work,10,output,"synthetic-not-DFT",
                           production_script=ROOT/"cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm")


def test_clamped_replay_excludes_reaction_stress_and_preserves_source(evaluated,tmp_path):
    work, images = evaluated
    before = {p:observer.sha256(p) for p in work.rglob("*") if p.is_file()}
    report = run(work,tmp_path/"export")
    free = VCNEB(images,cell_scale=observer.CELL_SCALE)
    assert np.linalg.norm(free.get_forces(),axis=1).max() > report["replayed_fmax_eV_A"]+.5
    assert report["new_DFT_calls"] == 0
    assert report["mechanical_boundary"]["cell_dofs"] == 3
    assert np.isclose(report["sampled_forward_barrier_meV_fu"],25)
    assert all(observer.sha256(p)==h for p,h in before.items())
    assert all(isinstance(a.calc,SinglePointCalculator)
               for a in read(tmp_path/"export/evaluated_chain.traj",index=":"))
    with pytest.raises(FileExistsError):
        run(work,tmp_path/"export")


@pytest.mark.parametrize("fault",["metric","free_cell","reference","cache_pin","raw_log",
                                  "STRU","missing_snapshot","wrong_force","wrong_energy","duplicate_row",
                                  "incompatible_cell","broken_lift","unregistered_script","geometry_cache"])
def test_clamped_observation_rejects_incomplete_or_mixed_evidence(evaluated,tmp_path,fault):
    work,_ = evaluated
    path = work/"vcneb_preflight.json"
    metadata = json.loads(path.read_text())
    if fault == "metric":
        metadata["cell_scale_A"] += .01
    elif fault == "free_cell":
        metadata["mechanical_boundary"] = None
    elif fault == "reference":
        Path(metadata["mechanical_boundary"]["reference_file"]).write_text("changed")
    elif fault == "cache_pin":
        metadata["calculator_parameters"]["seed_cache_records"][0]["raw_log_sha256"] = "0"*64
    elif fault in ("raw_log","STRU"):
        name = "OUT.ABACUS/running_scf.log" if fault == "raw_log" else "STRU"
        (work/"image_0004/scf_000000"/name).write_text("changed")
    elif fault == "missing_snapshot":
        (work/"snapshots/step_0010/POSCAR_04").unlink()
    elif fault in ("wrong_force","wrong_energy","duplicate_row"):
        log = work/"vcneb.opt.log"
        row = log.read_text()
        if fault == "duplicate_row":
            log.write_text(row+row)
        else:
            fields = row.split()
            fields[4 if fault=="wrong_force" else 3] = "10.0"
            log.write_text(" ".join(fields)+"\n")
    elif fault in ("incompatible_cell","broken_lift","geometry_cache"):
        p = work/"snapshots/step_0010/POSCAR_04"
        a = read(p)
        if fault == "incompatible_cell":
            cell = a.cell.array.copy(); cell[0,0] += .01; a.set_cell(cell,scale_atoms=True)
        elif fault == "broken_lift":
            a.positions[0] += a.cell[0]
        else:
            a.positions[0,0] += .001
        write(p,a,format="vasp",direct=True)
    elif fault == "unregistered_script":
        script = tmp_path/"wrong.slurm"; script.write_text("altered contract")
        with pytest.raises(ValueError,match="unregistered"):
            observer.export(work,10,tmp_path/"bad","synthetic",production_script=script)
        return
    path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError):
        run(work,tmp_path/"bad")
    assert not (tmp_path/"bad").exists()


@pytest.mark.parametrize("step",[-1,True,1.5])
def test_clamped_export_rejects_invalid_step(tmp_path,step):
    with pytest.raises(ValueError,match="nonnegative"):
        observer.export(tmp_path,step,tmp_path/"bad","synthetic",production_script=Path("unused"))


def test_prepare_continuation_has_nine_pins_without_new_independent_chain(evaluated,tmp_path):
    from scripts.prepare_hfo2_clamped_resume import prepare
    work,_ = evaluated
    report = run(work,tmp_path/"observation")
    assert not report["ordinary_residual_pass"]
    manifest = prepare(tmp_path/"observation",tmp_path/"resume")
    assert manifest["current_frame_exact_caches"]==9 and manifest["new_independent_chains"]==0
    assert manifest["new_DFT_calls"]==0 and not manifest["optimizer_state_restored"]
    parameters = json.loads((tmp_path/"resume/factory_parameters.json").read_text())
    assert all(r is not None and "STRU" in r["input_sha256"] for r in parameters["seed_cache_records"])
    assert all(observer.sha256(tmp_path/"resume"/n)==h for n,h in manifest["files_sha256"].items())
    with pytest.raises(FileExistsError):
        prepare(tmp_path/"observation",tmp_path/"resume")


@pytest.mark.parametrize("fault",["trajectory","snapshot","audit","preflight","converged"])
def test_prepare_continuation_rejects_altered_or_converged_observation(evaluated,tmp_path,fault):
    from scripts.prepare_hfo2_clamped_resume import prepare
    work,_ = evaluated
    root = tmp_path/"observation"
    report = run(work,root)
    if fault=="trajectory":
        (root/"evaluated_chain.traj").write_bytes(b"changed")
    elif fault=="snapshot":
        (root/"POSCAR_04").write_text("changed")
    elif fault=="audit":
        Path(report["raw_image_evaluations"][4]["audit_path"]).write_text("changed")
    elif fault=="preflight":
        (root/"runtime_preflight.json").write_text("changed")
    else:
        report["ordinary_residual_pass"] = True
        (root/"observation.json").write_text(json.dumps(report))
    with pytest.raises(ValueError):
        prepare(root,tmp_path/"bad")
    assert not (tmp_path/"bad").exists()
