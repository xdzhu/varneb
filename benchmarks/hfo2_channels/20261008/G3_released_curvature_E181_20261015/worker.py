"""Explicitly amended, finite34D stencils at two existing released Q=0 points.

External immutable bundle over the original E161 runtime. No retries,
optimizer, domain repair, source overwrite or independent-test access.
"""
from __future__ import annotations
import argparse, hashlib, json, os, subprocess, sys
from copy import deepcopy
from pathlib import Path
sys.dont_write_bytecode=True
if (Path(__file__).resolve().parents[4]/'vcneb').is_dir():
    sys.path.insert(0,str(Path(__file__).resolve().parents[4]))
import numpy as np
from scripts import hfo2_G3_conditional_faces as faces
from vcneb import assemble_joint_directional_curvature

study=faces.study
OLD=study.RUNS/'20261015-varneb-G3-faces-E161-r1/source'
RUNTIME_SHA='d07bb011c9e313afe23ca4e78207356ff98adfc43ca231bdd047e6d5d0f1be20'
POLICY=dict(channels=['T','preserving'],retained_Q_A=[0.,0.],released_dimension=34,
    step_A=[.01,.02],sign_order=[1,-1],new_attempts_per_center=136,new_attempts_total=272,
    additional_to_original_raw_peak_attempts=296,additional_verification_base_geometries=2,
    new_Q_training_nodes_per_face=0,training_committed=dict(T=34,preserving=35),
    independent_face_points_untouched_each=8,physical_strain=0.,pressure_GPa=0.,
    maximum_internal_radius_A=.5,minimum_distance_A=1.6,ordinary_NEB_fmax_eV_A=.10,
    base_released_gradient_tolerance_eV_A=.03,job_name='hfo2-G3-E181-release-Hrr',
    partition='hfacnormal01',nodes=1,MPI_ranks=32,allocation_hours=12,
    maximum_committed_study_allocations=2,automatic_retries=0,climb=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def safe_probe(data,full):
    """Reject only; the literal anchor is never reconstructed or moved inward."""
    x=np.asarray(full,dtype=float);b,a=data['basis'],data['axes']
    if (x.shape!=(40,) or not np.isfinite(x).all() or x[-1]!=0.
            or np.linalg.norm(b.T@x)>POLICY['maximum_internal_radius_A']
            or not np.allclose(b@(b.T@x),x,atol=2e-12,rtol=0)
            or not np.allclose(a.T@(b.T@x),0.,atol=2e-12,rtol=0)):
        raise ValueError('probe violates original fixed-Q/control/span/radius; no repair')
    atoms=data['chart'].to_atoms(x)
    if (not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell.array).all()
            or np.linalg.det(atoms.cell.array)<=0 or np.linalg.cond(atoms.cell.array)>1e12
            or not np.array_equal(atoms.cell.array[:2],data['chart'].reference_atoms.cell.array[:2])):
        raise ValueError('invalid probe cell or changed substrate')
    distances=atoms.get_all_distances(mic=True);np.fill_diagonal(distances,np.inf)
    if distances.min()<POLICY['minimum_distance_A']:
        raise ValueError('probe minimum distance; no repair')
    return atoms


def make_channel(data,r,*,declared=None):
    b,a=data['basis'],data['axes'];r=np.asarray(r,dtype=float)
    if (r.shape!=(36,34) or not np.isfinite(r).all()
            or not np.allclose(r.T@r,np.eye(34),atol=2e-12,rtol=0)
            or not np.allclose(a.T@r,0.,atol=2e-12,rtol=0)):
        raise ValueError('complete orthonormal fixed-Q complement required')
    # Exact native bytes in this coordinate convention, not B @ an inverse.
    base=np.asarray(data['origin_record']['coordinates']);directions=b@r
    if declared is not None:
        supplied=np.asarray(declared['literal_directions_40_by_34'])
        if (supplied.shape!=(40,34) or not np.isfinite(supplied).all()
                or not np.allclose(supplied,directions,atol=2e-12,rtol=0)
                or not np.allclose(supplied.T@supplied,np.eye(34),atol=2e-12,rtol=0)
                or np.any(supplied[-1])):
            raise ValueError('declared literal full directions differ')
        directions=supplied
    safe_probe(data,base)
    residual=float(np.linalg.norm(directions.T@data['origin_record']['gradient']))
    if residual>POLICY['base_released_gradient_tolerance_eV_A']:
        raise ValueError('original base is not a released-gradient candidate')
    if declared is not None:
        supplied=declared['base_released_gradient_norm_eV_A']
        if (isinstance(supplied,bool) or not np.isfinite(supplied)
                or not np.isclose(supplied,residual,atol=2e-12,rtol=0)):
            raise ValueError('base residual differs')
        residual=supplied
    probes=[]
    for h in POLICY['step_A']:
        for i in range(34):
            for sign in POLICY['sign_order']:
                x=base+directions[:,i]*(sign*h);safe_probe(data,x)
                probes.append(dict(id=len(probes),step_A=h,direction=i,sign=sign,coordinates_full_A=x.tolist()))
    if len({np.asarray(p['coordinates_full_A']).tobytes() for p in probes})!=136:
        raise ValueError('136unique fresh paired probes required')
    return dict(released_basis_36_by_34=r.tolist(),literal_directions_40_by_34=directions.tolist(),
        literal_native_anchor_40_A=base.tolist(),
        native_anchor_result=data['origin_record'],base_released_gradient_norm_eV_A=residual,probes=probes)


def prepare(bundle,*,declared=None):
    channels={}
    for c in POLICY['channels']:
        data=faces.load_inputs(c)
        r=(np.linalg.svd(data['axes'].T,full_matrices=True)[2][2:].T if declared is None
           else declared['channels'][c]['released_basis_36_by_34'])
        channels[c]=make_channel(data,r,declared=None if declared is None else declared['channels'][c])
    return dict(schema_version=1,purpose='post_release_transverse_derivative_validation_not_TS_or_forecast',
        policy=deepcopy(POLICY),channels=channels,original_input_artifact_pins=deepcopy(faces.PINS),
        physical_input_sha256=deepcopy(study.CONTRACT),binary_sha256=study.BINARY_SHA,
        original_runtime_manifest_sha256=RUNTIME_SHA,explicit_budget_amendment_sha256=sha(bundle/'protocol.md'),
        worker_text_sha256=study.text_sha(Path(__file__)),wrapper_text_sha256=study.text_sha(bundle/'run.slurm'),
        physical_reserved_geometry_generated_or_read=False,independent_prediction_or_stability_certified=False)


def validate(config,bundle):
    expected=prepare(bundle,declared=config)
    # Only already-frozen bases are used to regenerate probes: LAPACK signs
    # are not changed, coordinate tolerances do not weaken exact probe bytes.
    if json.dumps(config,sort_keys=True,allow_nan=False)!=json.dumps(expected,sort_keys=True,allow_nan=False):
        raise ValueError('exact stencil/source/physics/amended-budget contract differs')
    return {c:faces.load_inputs(c) for c in POLICY['channels']}


def scheduler_snapshot(own=None):
    result=subprocess.run(['squeue','-h','-u','iai806','-p','hfacnormal01','-o','%i|%100j|%T|%C|%D'],
        check=True,capture_output=True,text=True).stdout
    committed=[]
    for line in result.splitlines():
        f=[v.strip() for v in line.split('|')]
        if len(f)!=5:raise ValueError('unparsed queue; no admission')
        if f[1].startswith('hfo2-G3-'):
            if f[3:]!=['32','1']:raise ValueError('study resources differ')
            committed.append(dict(job_id=f[0],name=f[1],state=f[2],CPUs=32,nodes=1))
    ready=len(committed)<2 if own is None else (len(committed)<=2 and
        any(r['job_id']==own and r['name']==POLICY['job_name'] for r in committed))
    return dict(ready=ready,committed=committed,own_job=own,scheduler_writes=0)


def preflight(config,bundle):
    data=validate(config,bundle)
    if study.ROOT.resolve()!=OLD.resolve():raise ValueError('original immutable E161 runtime required')
    files=faces.canary.pinned(OLD/'runtime_files.json',RUNTIME_SHA)
    if len(files)!=380 or any(sha(OLD/n)!=h for n,h in files.items()):raise ValueError('old380runtime drift')
    for c in POLICY['channels']:
        faces.native_preflight(faces.prepare(c),c)  # exact original and released native E/F/stress
    return data,dict(status='native_two_released_base_preflight_no_DFT',original_runtime_files=380,
        released_base_native_records_reparsed=2,inert_probe_geometries_checked=272,
        declared_additional_derivative_budget=272,admission=scheduler_snapshot(),new_DFT=0,
        physical_reserved_geometry_generated_or_read=False)


def blocks_from_gradients(directions,gradients):
    gradients=np.asarray(gradients,dtype=float)
    if gradients.shape!=(136,40) or not np.isfinite(gradients).all():
        raise ValueError('complete136finite full gradients required; no partial matrix')
    blocks=[]
    for j,h in enumerate(POLICY['step_A']):
        g=gradients[j*68:(j+1)*68]
        block=assemble_joint_directional_curvature(directions,g[::2],g[1::2],step_A=h)
        blocks.append(dict(step_A=h,complete_pairs=34,full_hessian_action_eV_A2=block.full_hessian_action.tolist(),
            raw_projected_eV_A2=block.raw_projected.tolist(),
            symmetric_projected_eigenvalues_eV_A2=np.linalg.eigvalsh(block.symmetric_projected).tolist(),
            reciprocity_relative_defect=block.reciprocity_relative_defect,
            symmetrization_operator_change_eV_A2=block.symmetrization_operator_change_eV_A2))
    return blocks


def run(config,bundle,workdir):
    if (os.environ.get('RUN_DFT')!='1' or os.environ.get('SLURM_JOB_PARTITION')!='hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES')!='1' or os.environ.get('SLURM_NTASKS')!='1'
            or os.environ.get('SLURM_CPUS_PER_TASK')!='32' or not os.environ.get('SLURM_JOB_ID')
            or os.environ.get('SLURM_JOB_NAME')!=POLICY['job_name']
            or any(os.environ.get(n)!='1' for n in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'))):
        raise ValueError('authorized1node32MPI/singlethread allocation required')
    if (workdir.resolve()!=bundle.resolve().parent/'released_curvature' or workdir.exists()
            or not workdir.resolve().is_relative_to(study.RUNS)):
        raise ValueError('fresh declared output beside external bundle required')
    data,proof=preflight(config,bundle)
    if not scheduler_snapshot(os.environ['SLURM_JOB_ID'])['ready']:raise ValueError('committed allocation gate failed')
    workdir.mkdir();study.exclusive_json(workdir/'launch_preflight.json',proof)
    totals=[]
    for c in POLICY['channels']:
        d=data[c];spec=config['channels'][c];directory=workdir/c;directory.mkdir()
        study.exclusive_json(directory/'reused_native_base.json',spec['native_anchor_result'])
        evaluator=faces.canary.CalculatorModeEvaluator(d['chart'],
            study.make_factory(parameters={'source_directory':str(study.ASSETS)},command=study.COMMAND),
            workdir=directory/'evaluations',calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
            minimum_distance_A=1.6,result_validator=study.result_validator,
            validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1',max_evaluation_attempts=136)
        energies=[];gradients=[]
        for p in spec['probes']:
            x=np.asarray(p['coordinates_full_A']);safe_probe(d,x)
            e,g=evaluator(x);energies.append(e);gradients.append(g)
            print(json.dumps(dict(channel=c,completed=p['id']+1,total=136,step_A=p['step_A'],
                direction=p['direction'],sign=p['sign'],energy_eV_cell=e)),flush=True)
        records=[json.loads(p.read_bytes()) for p in sorted((directory/'evaluations').glob('eval-*/result.json'))]
        if (len(records)!=136 or evaluator.n_evaluation_attempts!=136 or evaluator.n_new_evaluations!=136
                or len(list((directory/'evaluations').glob('attempt-*.json')))!=136):
            raise ValueError('all136fresh native attempts required')
        directions=np.asarray(spec['literal_directions_40_by_34'])
        blocks=blocks_from_gradients(directions,gradients)
        total_seconds=sum(r['validation']['elapsed_seconds'] for r in records)
        result=dict(status='complete_post_release_raw_transverse_curvature_not_stability_certificate',
            channel=c,new_native_SCFs=136,new_center_SCFs=0,energies_eV_cell=energies,
            full_gradients_eV_A=np.asarray(gradients).tolist(),blocks=blocks,
            two_step_raw_block_operator_difference_eV_A2=float(np.linalg.norm(
                np.asarray(blocks[0]['raw_projected_eV_A2'])-blocks[1]['raw_projected_eV_A2'],ord=2)),
            base_released_gradient_norm_eV_A=spec['base_released_gradient_norm_eV_A'],
            sum_native_call_seconds=total_seconds,call_time_proxy_core_hours=total_seconds*32/3600,
            physical_stability_TS_or_independent_prediction_certified=False)
        study.exclusive_json(directory/'curvature.json',result)
        totals.append(dict(channel=c,new_native_SCFs=136,call_time_proxy_core_hours=total_seconds*32/3600))
    study.exclusive_json(workdir/'completed.json',dict(status='complete_both_released_stencils',
        channels=totals,new_native_SCFs=272,new_center_SCFs=0,job_id=os.environ['SLURM_JOB_ID'],
        allocation_core_hours=None,all_other_training_and_preparation_costs_additional=True))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('prepare','preflight','run'));p.add_argument('--bundle',type=Path,required=True)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--config-sha256')
    p.add_argument('--output',type=Path);p.add_argument('--workdir',type=Path);a=p.parse_args()
    if a.mode=='prepare':study.exclusive_json(a.config,prepare(a.bundle));return
    if sha(a.config)!=a.config_sha256:raise ValueError('exact frozen config hash required')
    config=json.loads(a.config.read_bytes())
    if a.mode=='preflight':
        _,proof=preflight(config,a.bundle)
        if a.output is None:raise ValueError('fresh preflight report required')
        study.exclusive_json(a.output,proof);print(json.dumps(proof))
    else:
        if a.workdir is None:raise ValueError('fresh run output required')
        run(config,a.bundle,a.workdir)


if __name__=='__main__':main()
