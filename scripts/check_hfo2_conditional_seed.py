"""Quadratic implementation replay from E155, not a DFT speed benchmark.

The spectral plane here belongs only to an exact quadratic surrogate. It does
not select/freeze a physical G3 surface, model, branch or predictive domain.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from vcneb import ModePlane, relax_orthogonal_at_q

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/hfo2_channels/20261008'
PINS={'T':'1cc4269a8f1437acddb94c9a7c264f52bae888669db2bcc36521982d280666ea',
      'preserving':'9bfeb3c8a7eca7533246ffd61e6e43df17262fb8ef274dfa87965c1301a10efa'}
MANIFEST_SHA='5f047f54ca83cb2dceb8cb71687eb9206a85c04205cdf26431565bc8a9b98a5a'


def check():
    path=BASE/'G3_curvature_E150_20261014/manifest.json'
    if hashlib.sha256(path.read_bytes()).hexdigest()!=MANIFEST_SHA:
        raise ValueError('original measured-chart manifest differs')
    manifest=json.loads(path.read_bytes());rows=[]
    for channel,pin in PINS.items():
        path=BASE/'G3_first_amplitude_E155_20261014'/(channel+'_native.json')
        if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
            raise ValueError('original native first-step audit differs')
        report=json.loads(path.read_bytes());block=report['blocks'][0]
        if not block['full_stencil_at_this_step'] or block['step_A']!=.01:
            raise ValueError('complete registered first amplitude required')
        raw=np.asarray(block['raw_projected_eV_A2'])[:36,:36];h=(raw+raw.T)/2
        b=np.asarray(manifest['centers'][channel]['directions'])[:,:36]
        g=b.T@np.asarray(report['raw_center_gradient_eV_A'])
        values,vectors=np.linalg.eigh(h)
        if np.count_nonzero(values<0)!=2 or np.any(values==0):
            raise ValueError('this quadratic replay requires the observed two-negative spectrum')
        plane=ModePlane(np.zeros(36),vectors,np.ones(36),np.eye(36)[:,:2],
            ('quadratic spectral1','quadratic spectral2'),'local length','E155-quadratic-replay-only')
        stable=vectors[:,2:]
        oracle=-stable@((stable.T@g)/values[2:])
        oracle_energy=float(.5*oracle@h@oracle+g@oracle)
        for label,seed in (('identity',None),('exact_quadratic_hessian',h)):
            called=[]
            def evaluate(x):
                called.append(x.copy())
                return float(.5*x@h@x+g@x),h@x+g
            try:
                point=relax_orthogonal_at_q(plane,[0.,0.],evaluate,
                    optimizer='safeguarded_bfgs',trial_validator=lambda x: bool(np.isfinite(x).all()),
                    initial_hessian=seed,initial_trust_radius=.05,
                    gradient_tolerance=1e-8,max_iterations=800)
                if point.total_n_evaluations!=len(called):
                    raise ValueError('callback accounting differs')
                row=dict(status='converged',callback_evaluations=len(called),
                    coordinate_oracle_error_A=float(np.linalg.norm(point.coordinates-oracle)),
                    energy_oracle_error_eV=float(abs(point.energy-oracle_energy)),
                    final_orthogonal_gradient_norm_eV_A=point.orthogonal_gradient_norm)
            except RuntimeError as exc:
                row=dict(status='unconverged',reason=str(exc),callback_evaluations=len(called))
            rows.append(dict(channel=channel,optimizer_seed=label,**row))
    return dict(status='exact_quadratic_implementation_replay_not_material_acceleration',rows=rows,
        native_audit_input_sha256=PINS,manifest_sha256=MANIFEST_SHA,
        common_gradient_tolerance_eV_A=1e-8,common_initial_trust_radius_A=.05,
        quadratic_spectral_coordinates_are_not_registered_physical_axes=True,
        new_DFT_or_scheduler_operations=0,physical_curvature_or_TS_certified=False,
        physical_conditional_axes_domains_or_forecasts_selected=False,
        physical_reserved_geometry_generated_or_read=False,
        DFT_speedup_or_novel_algorithm_claimed=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError('fresh quadratic replay report required')
    result=check();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(status=result['status'],rows=result['rows'])))


if __name__=='__main__':main()
