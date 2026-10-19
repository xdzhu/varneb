"""New finite continuation namespace; preserves the original +1% F and controls."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import tempfile
import numpy as np
from diagnose import CASE, OLD, module
from scripts.audit_hfo2_static_replica import sha256

POLICY=dict(strain=.01, retained=['Qx','Qz','Qy','polar_proxy_z'],
    original_F_job='28812952', original_local_job='28948425',
    starts=['minus','plus'], original_norm_guard_A=.5, maximum_coordinate_norm_A=1.,
    start='last exact native evaluated trial of each original signed A arm',
    atomic_attempts=32, joint_attempts=16, direction_probes=2, maximum_new_calls_per_sign=50,
    gradient_tolerance_eV_A=.03, initial_trust_radius_A=.05, maximum_iterations=200,
    minimum_distance_A=1.6, maximum_deformation=.25, known_direction_probe_step_A=.01,
    optimizer='safeguarded_bfgs', numerical_initializer='unchanged E272 branch/on-site matrix',
    new_physical_Hessian=False, physical_inputs_changed=False, old_results_or_labels_replaced=False,
    original_protocol_success_claimed=False, stationary_TS_or_full_stability_certified=False,
    nodes_per_sign=1, MPI_ranks=32, partition='hfacnormal01', walltime='04:00:00',
    maximum_study_allocations=2, automatic_retry_or_budget_extension=False)


def prepare(output):
    output=Path(output)
    if output.exists(): raise FileExistsError(output)
    diagnosis=json.loads((CASE/'original_failure_diagnosis/diagnosis.json').read_bytes())
    scale=json.loads((CASE/'registered_connection_scale.json').read_bytes())
    if (diagnosis['local_native_calls_reparsed']!=24 or diagnosis['new_DFT_calls']
            or not scale['registration_fixed_by_initial_and_verified_at_final_endpoint']
            or not .5<scale['projected_same_observable_atomic_norm_A']<.9):
        raise ValueError('original failed calls and registered geometric rationale required')
    reader=module(OLD/'analyze_prefix.py','E276_prepare_reader')
    terminal=module(OLD/'analyze_terminal.py','E276_prepare_archive')
    chart,data,e0,g0=reader.setup(2)
    output.mkdir(parents=True)
    shutil.copytree(OLD/'prepared/seeds/2',output/'seed')
    with tempfile.TemporaryDirectory(prefix='E276-prepare-') as td:
        root,_=terminal.restore(Path(td))
        starts=[]
        for arm in diagnosis['arms']:
            sign=arm['sign'];source=root/arm['last_native_relative_path']
            if sha256(source)!=arm['last_result_sha256']:raise ValueError('the original last native trial changed')
            original=source.parent
            target=output/'starts'/sign
            shutil.copytree(original,target)
            starts.append(dict(sign=sign,source_result_sha256=sha256(target/'result.json'),
                coordinate_norm_A=float(np.linalg.norm(arm['last_coordinates'])),
                original_status='attempt_budget_exhausted_not_a_candidate'))
    cfg=dict(policy=POLICY,physical_input_sha256=data['physical_input_sha256'],
        diagnosis_sha256=sha256(CASE/'original_failure_diagnosis/diagnosis.json'),
        registered_scale_sha256=sha256(CASE/'registered_connection_scale.json'),
        original_F_energy_eV_cell=e0, original_F_gradient=g0.tolist(),starts=starts,
        new_DFT_calls=0,execution_admitted=False)
    with (output/'config.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(cfg,f,indent=2,allow_nan=False);f.write('\n')
    return cfg


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();r=prepare(args.output);print(json.dumps(r['policy']))
