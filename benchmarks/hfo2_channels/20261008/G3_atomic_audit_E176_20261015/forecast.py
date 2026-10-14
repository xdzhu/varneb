"""Original-center B3 expectations before reading atomic-treatment labels.

All ten points are registered training treatments, not the independent eight
points/face or physical+.005. No parameter fit, favorable-step selection or
new axes are introduced by these formal quadratic evaluations.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT))
from vcneb.quadratic_reduction import condition_quadratic_energy
BASE=ROOT/'benchmarks/hfo2_channels/20261008'
CASE=BASE/'G3_atomic_release_E175_20261015'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def signed_treatment_changes(frozen,atomic,joint,*,origins):
    """Signed EA-EF and EJ-EA; finite local candidates need not be nested minima."""
    if atomic is None:
        return dict(status='atomic_candidate_or_origin_missing',absolute=None,shape=None)
    values=np.asarray([frozen,atomic,joint],dtype=float)
    if not np.isfinite(values).all():raise ValueError('finite same-contract treatment energies required')
    f,a,j=values
    absolute=dict(atomic_minus_frozen_eV=float(a-f),joint_minus_atomic_eV=float(j-a),
        joint_minus_frozen_eV=float(j-f))
    if any(v is None for v in origins):
        return dict(status='origin_missing_absolute_changes_retained',absolute=absolute,shape=None)
    f0,a0,j0=np.asarray(origins,dtype=float)
    if not np.isfinite([f0,a0,j0]).all():raise ValueError('finite same-contract origins required')
    shape=dict(atomic_minus_frozen_eV=float((a-a0)-(f-f0)),
        joint_minus_atomic_eV=float((j-j0)-(a-a0)),
        joint_minus_frozen_eV=float((j-j0)-(f-f0)))
    for block in (absolute,shape):
        if not np.isclose(block['atomic_minus_frozen_eV']+block['joint_minus_atomic_eV'],
                          block['joint_minus_frozen_eV'],atol=2e-11,rtol=0):
            raise ValueError('signed pointwise treatment identity differs')
    return dict(status='signed_candidate_difference_not_global_minimum_or_unique_mode_energy',
        absolute=absolute,shape=shape)


def forecast():
    if sha(CASE/'config.json')!='5bd401c6ca4044ac654f2b19a22371b20c17d3736d207ff7044336d9ed463f38':
        raise ValueError('exact E175 fixed-Q atomic-release configuration required')
    spec=importlib.util.spec_from_file_location('original_atomic_release_for_expectations',CASE/'worker.py')
    w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
    config=json.loads((CASE/'config.json').read_bytes());data=w.validate(config,CASE)
    rows=[]
    for channel in w.POLICY['channels']:
        d=data[channel];a=d['axes'];r=np.asarray(config['channels'][channel]['atomic_release_basis_36_by_31'])
        joint=np.linalg.svd(a.T,full_matrices=True)[2][2:].T;g=d['basis'].T@d['center_gradient']
        for block in d['curvature']['blocks']:
            raw=np.asarray(block['raw_projected_eV_A2'])[:36,:36];h=(raw+raw.T)/2
            atom=condition_quadratic_energy(h,g,d['center_energy'],a,r,stability_floor=0.)
            full=condition_quadratic_energy(h,g,d['center_energy'],a,joint,stability_floor=0.)
            origins=[0.,atom.reference_energy_change,full.reference_energy_change];points=[]
            for p in config['channels'][channel]['points']:
                q=np.asarray(p['q_A']);x=a@q
                f=float(g@x+.5*x@h@x);ea=atom.energy_change(q);ej=full.energy_change(q)
                points.append(dict(id=p['id'],q_A=p['q_A'],
                    model_frozen_energy_change_eV_cell=f,model_atomic_energy_change_eV_cell=ea,
                    model_joint_energy_change_eV_cell=ej,
                    model_atomic_released_coordinates=r.shape[1],
                    **signed_treatment_changes(f,ea,ej,origins=origins)))
            rows.append(dict(channel=channel,probe_step_A=block['step_A'],points=points,
                formal_atomic_origin_release_meV_cell=1000*origins[1],
                formal_joint_origin_release_meV_cell=1000*origins[2],
                minimum_raw_atomic_release_curvature_eV_A2=atom.curvature.minimum_eliminated_curvature))
    return dict(status='unfitted_original_quadratic_B3_training_expectations_not_independent_forecast',
        rows=rows,configuration_sha256=sha(CASE/'config.json'),
        all_two_centers_two_steps_five_points=20,raw_center_origin_offsets_retained=True,
        fitted_coefficients=0,new_DFT=0,atomic_treatment_result_labels_read=False,
        physical_reserved_geometry_generated_or_read=False,independent_slots_used=0,
        post_relaxation_stability_or_physical_error_or_prediction_advantage_certified=False,
        source_text_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_text(encoding='utf-8').encode()).hexdigest()
            for p in (Path(__file__),ROOT/'vcneb/quadratic_reduction.py',CASE/'worker.py')})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('fresh expectation output required')
    result=forecast()
    with a.output.open('x',encoding='utf-8',newline='\n') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps([dict(channel=r['channel'],probe_step_A=r['probe_step_A'],
        atomic_origin_meV=r['formal_atomic_origin_release_meV_cell'],
        joint_origin_meV=r['formal_joint_origin_release_meV_cell']) for r in result['rows']]))
