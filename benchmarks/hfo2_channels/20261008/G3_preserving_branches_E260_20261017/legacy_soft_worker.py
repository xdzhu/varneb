"""Exactly four new training statics along the frozen E255 soft direction."""
from pathlib import Path
from copy import deepcopy
import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
sys.dont_write_bytecode = True
LOCAL = Path(__file__).resolve().parents[4]
if (LOCAL / 'vcneb').is_dir():
    sys.path.insert(0, str(LOCAL))
import numpy as np
from scripts import hfo2_G3_conditional_faces as faces
from examples.hfo2_fixed_input_factory import write_structure, geometry_roundtrip

study = faces.study
EXEC = study.RUNS / '20261016-varneb-G3-soft4-E256-r1'
READER = study.RUNS / '20261015-varneb-G3-released-audit-E190-r1/audit.py'
READER_SHA = '51d0452968394ff03d26efface7705ee31c7679215cb2204f6c69e13a5391d4f'
CONFIG_SHA = '71a2e1b8a46b2ef903e4023d329109c5e74975aaf5bab54f4d4638d15c1add39'
AUDIT_SHA = 'a389154749a07109a3a91499343b86f06a2c185fd310611570637d71193e3381'
SOFT_SHA = '1e5dcd01ec5c7731646c9e1f6b0f230fc8bcf45bbcad5a785125534a0d2ee065'
POLICY = dict(t_A=[.05,-.05,.10,-.10], maximum_new_attempts=4,
    original_training_committed=dict(T=34,preserving=35),
    amended_training_committed=dict(T=34,preserving=39), training_cap_each=81,
    independent_reserved_each=8, independent_computed_each=0,
    maximum_internal_radius_A=.5, minimum_distance_A=1.6,
    partition='hfacnormal01', nodes=1, MPI_ranks=32, allocation_hours=1,
    maximum_committed_study_allocations=2, job_name='hfo2-G3-E256-soft4',
    automatic_retries=0, optimizer=None, climb=False)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def pinned(path, digest):
    if sha(path) != digest:
        raise ValueError('frozen source bytes differ')
    return json.loads(path.read_bytes())


def inputs(bundle):
    soft = pinned(bundle/'soft_direction.json', SOFT_SHA)
    if (LOCAL/'vcneb').is_dir():
        base = LOCAL/'benchmarks/hfo2_channels/20261008'
        cp = base/'G3_released_curvature_E181_20261015/config.json'
        ap = base/'G3_released_complete_E254_20261016/native_audit.json'
    else:
        cp = study.RUNS/'20261015-varneb-G3-release-curvature-E181-r1/bundle/config.json'
        ap = study.RUNS/'20261015-varneb-G3-released-audit-E190-r1/complete_20261016_2331.json'
    old = pinned(cp, CONFIG_SHA)
    audit = pinned(ap, AUDIT_SHA)
    if audit['mode'] != 'complete' or not audit['accounting']['normal_completion']:
        raise ValueError('actual complete native curvature required')
    data = faces.load_inputs('preserving')
    x0 = np.asarray(old['channels']['preserving']['literal_native_anchor_40_A'])
    direction = np.asarray(soft['soft_direction_full40'])
    if (direction.shape != (40,) or direction[-1] != 0.
            or not np.allclose(direction@direction,1.,atol=2e-12,rtol=0)
            or not np.array_equal(x0,data['origin_record']['coordinates'])
            or soft['original_independent_sign_admission_passed'] is not False):
        raise ValueError('original literal direction/anchor and refusal required')
    return data, x0, direction, soft


def checked_geometry(data, x):
    b,a = data['basis'],data['axes']
    if (x.shape != (40,) or not np.isfinite(x).all() or x[-1] != 0.
            or np.linalg.norm(b.T@x) > .5
            or not np.allclose(b@(b.T@x),x,atol=2e-12,rtol=0)
            or not np.allclose(a.T@(b.T@x),0.,atol=2e-12,rtol=0)):
        raise ValueError('literal probe violates original Q/control/span/radius; no repair')
    atoms = data['chart'].to_atoms(x)
    cell = atoms.cell.array
    dist = atoms.get_all_distances(mic=True); np.fill_diagonal(dist,np.inf)
    if (np.linalg.det(cell) <= 0 or np.linalg.cond(cell) > 1e12
            or dist.min() < 1.6 or not np.isfinite(atoms.positions).all()
            or not np.array_equal(cell[:2],data['chart'].reference_atoms.cell.array[:2])):
        raise ValueError('probe violates original geometry/substrate')
    return atoms


def prepare(bundle):
    data,x0,direction,soft = inputs(bundle)
    points = []
    for t in POLICY['t_A']:
        x = x0+t*direction; atoms = checked_geometry(data,x)
        points.append(dict(id=len(points),t_A=t,coordinates_full_A=x.tolist(),
            quadratic_delta_eV_cell=[soft['soft_gradient_eV_A']*t+.5*k*t*t
                for k in soft['first_direction_Rayleigh_eV_A2']]))
    return dict(schema_version=1,purpose='four_training_energy_force_checks_not_independent_test',
        policy=deepcopy(POLICY),points=points,soft_source_sha256=SOFT_SHA,
        old_config_sha256=CONFIG_SHA,complete_native_audit_sha256=AUDIT_SHA,
        protocol_sha256=sha(bundle/'protocol.md'),worker_text_sha256=study.text_sha(Path(__file__)),
        wrapper_text_sha256=study.text_sha(bundle/'run.slurm'),
        physical_input_sha256=study.CONTRACT,binary_sha256=study.BINARY_SHA,
        original_independent_admission_refusal_preserved=True,
        physical_plus005_generated_or_read=False)


def validate(config,bundle):
    if json.dumps(config,sort_keys=True,allow_nan=False) != json.dumps(prepare(bundle),sort_keys=True,allow_nan=False):
        raise ValueError('frozen four-point/source/physics/budget contract differs')
    return inputs(bundle)


def queue(own=None):
    raw = subprocess.run(['squeue','-h','-u','iai806','-p','hfacnormal01','-o','%i|%100j|%T|%C|%D'],
        check=True,capture_output=True,text=True).stdout
    committed = []
    for line in raw.splitlines():
        f = [x.strip() for x in line.split('|')]
        if len(f) != 5: raise ValueError('unparsed queue')
        if f[1].startswith('hfo2-G3-'):
            if f[3:] != ['32','1']: raise ValueError('study resource layout differs')
            committed.append(dict(job_id=f[0],name=f[1],state=f[2],CPUs=32,nodes=1))
    ready = len(committed)<2 if own is None else len(committed)<=2 and any(
        r['job_id']==own and r['name']==POLICY['job_name'] for r in committed)
    return dict(ready=ready,committed=committed,own_job=own,scheduler_writes=0)


def native_context(config,bundle):
    data,anchor,direction,soft = validate(config,bundle)
    if sha(READER) != READER_SHA: raise ValueError('original native reader differs')
    s = importlib.util.spec_from_file_location('E256_original_native_reader',READER)
    reader = importlib.util.module_from_spec(s); s.loader.exec_module(reader)
    w,_,original = reader.load_worker()  # pins380 runtime/binary/six assets
    data = original['preserving']; record = data['origin_record']
    raw = Path(record['validation']['raw_directory'])
    old,_ = reader.native_record(w,data,anchor,raw.parent/'result.json',raw.parent.parent)
    if old != record: raise ValueError('native cached origin differs')
    return reader,w,data,anchor,direction,soft


def preflight(config,bundle,workdir):
    if workdir.exists() or workdir.resolve()!=EXEC/'inert_decks' or bundle.resolve()!=EXEC/'bundle':
        raise ValueError('fresh exact external namespace required')
    _,_,data,_,_,_ = native_context(config,bundle)
    workdir.mkdir(); decks=[]
    for point in config['points']:
        atoms=checked_geometry(data,np.asarray(point['coordinates_full_A']))
        directory=workdir/f"point_{point['id']:03d}"; directory.mkdir()
        for name in study.CONTRACT:
            if name in ('INPUT','KPT'): shutil.copyfile(study.ASSETS/name,directory/name)
            else: (directory/name).symlink_to(study.ASSETS/name)
        write_structure(directory/'STRU',atoms); geometry_roundtrip(directory/'STRU',atoms)
        hashes={name:sha(directory/name) for name in (*study.CONTRACT,'STRU')}
        if any(hashes[n]!=h for n,h in study.CONTRACT.items()): raise ValueError('inert inputs drift')
        dist=atoms.get_all_distances(mic=True); np.fill_diagonal(dist,np.inf)
        decks.append(dict(id=point['id'],t_A=point['t_A'],input_sha256=hashes,
            minimum_distance_A=float(dist.min()),volume_A3=float(atoms.get_volume())))
    return dict(status='four_native_inert_decks_passed_no_DFT',decks=decks,
        config_sha256=sha(bundle/'config.json'),admission=queue(),new_DFT=0,
        original_independent_admission_refusal_preserved=True)


def run(config,bundle,workdir):
    job=os.environ.get('SLURM_JOB_ID')
    required=dict(RUN_DFT='1',SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',
        SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='32',SLURM_JOB_NAME=POLICY['job_name'],
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    if not job or any(os.environ.get(k)!=v for k,v in required.items()):
        raise ValueError('authorized one-node32MPI allocation required')
    if workdir.exists() or workdir.resolve()!=EXEC/'probes' or bundle.resolve()!=EXEC/'bundle':
        raise ValueError('fresh exact output namespace required; no rerun')
    accepted=json.loads((EXEC/'accepted.json').read_bytes())
    if accepted['job_id']!=job or accepted['config_sha256']!=sha(bundle/'config.json') or not queue(job)['ready']:
        raise ValueError('accepted handle/config/committed-allocation gate failed')
    reader,w,data,_,direction,_=native_context(config,bundle)
    workdir.mkdir()
    evaluator=faces.canary.CalculatorModeEvaluator(data['chart'],
        study.make_factory(parameters={'source_directory':str(study.ASSETS)},command=study.COMMAND),
        workdir=workdir/'evaluations',calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
        minimum_distance_A=1.6,result_validator=study.result_validator,
        validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',max_evaluation_attempts=4)
    for p in config['points']:
        x=np.asarray(p['coordinates_full_A']); checked_geometry(data,x)
        e,g=evaluator(x)
        print(json.dumps(dict(id=p['id'],t_A=p['t_A'],energy_eV_cell=e,
            gradient_along_direction_eV_A=float(direction@g))),flush=True)
    paths=sorted((workdir/'evaluations').glob('eval-*/result.json'))
    if len(paths)!=4 or evaluator.n_evaluation_attempts!=4 or evaluator.n_new_evaluations!=4:
        raise ValueError('exactly four original fresh attempts required')
    rows=[];sources=[]
    for p,path in zip(config['points'],paths):
        row,source=reader.native_record(w,data,np.asarray(p['coordinates_full_A']),path,workdir/'evaluations')
        rows.append(row); sources.append(source)
    study.exclusive_json(workdir/'completed.json',dict(status='four_native_training_statics_complete_not_branch_or_TS',
        job_id=job,config_sha256=sha(bundle/'config.json'),new_native_SCFs=4,new_center_SCFs=0,
        records=rows,native_sources=sources,original_independent_admission_refusal_preserved=True,
        independent_points_computed=0,physical_plus005_generated_or_read=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('prepare','preflight','run'))
    p.add_argument('--bundle',type=Path,required=True);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--config-sha256');p.add_argument('--workdir',type=Path);p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.mode=='prepare': study.exclusive_json(a.config,prepare(a.bundle))
    else:
        config=pinned(a.config,a.config_sha256)
        if a.mode=='preflight':
            proof=preflight(config,a.bundle,a.workdir);study.exclusive_json(a.output,proof)
            print(json.dumps(dict(status=proof['status'],decks=len(proof['decks']),new_DFT=0)))
        else: run(config,a.bundle,a.workdir)
