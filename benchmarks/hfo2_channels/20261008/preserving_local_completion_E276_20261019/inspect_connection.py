"""Geometric scale of the existing +1% free band in the old local controls."""
from pathlib import Path
import argparse
import importlib.util
import json
import tempfile
import numpy as np
from ase.io import read as read_atoms
from diagnose import CASE, OLD, REPO, module
from scripts.audit_hfo2_static_replica import sha256
from vcneb import validate_path_geometry


def inspect():
    prior=module(OLD/'analyze_prefix.py','E276_original_chart')
    chart,data,_,_=prior.setup(2)
    connected=CASE.with_name('preserving_strain_E275_20261019')
    read=module(connected/'analyze_terminal.py','E276_connected_native')
    with tempfile.TemporaryDirectory(prefix='E276-connected-scale-') as td:
        band,_=read.read_terminal('strain_0100',Path(td))
        images,obs=read.observation_images(Path(td)/'native_terminal/observation')
        peak=band['sampled_peak_image']
        # One registered per-atom integer offset for the entire polygon,
        # fixed by the common initial endpoint, never peak-wise MIC matching.
        original=REPO/Path(data['observation'])
        oldobs=json.loads(original.read_bytes())
        oldpeak=int(np.argmax(oldobs['relative_enthalpy_meV_fu']))
        oldcenter=read_atoms(original.parent/f'POSCAR_{oldpeak:02d}')
        native_offset=chart.reference_atoms.get_scaled_positions(wrap=False)-oldcenter.get_scaled_positions(wrap=False)
        if not np.allclose(native_offset,np.rint(native_offset),atol=1e-10,rtol=0):
            raise ValueError('original native center needs an explicit integer lift')
        oldinitial=read_atoms(original.parent/'POSCAR_00')
        offset=oldinitial.get_scaled_positions(wrap=False)+np.rint(native_offset)-images[0].get_scaled_positions(wrap=False)
        if not np.allclose(offset,np.rint(offset),atol=1e-10,rtol=0):
            raise ValueError('the exact ordered initial well does not establish an integer registration')
        offset=np.rint(offset)
        oldfinal=read_atoms(original.parent/'POSCAR_08')
        if not np.allclose(oldfinal.get_scaled_positions(wrap=False)+np.rint(native_offset),
                images[-1].get_scaled_positions(wrap=False)+offset,atol=1e-9,rtol=0):
            raise ValueError('one fixed whole-polygon lift must also match the final ordered well')
        registered=images[peak].copy()
        registered.set_scaled_positions(registered.get_scaled_positions(wrap=False)+offset)
        x=chart.from_atoms(registered)
    q=np.asarray(data['orthonormal_39_by_4']);t=np.zeros((39,3))
    t[:36]=np.tile(np.eye(3),(12,1))/np.sqrt(12)
    atomic=x-q@(q.T@x)-t@(t.T@x);atomic[36:]=0.
    metrics=validate_path_geometry([chart.reference_atoms,chart.to_atoms(atomic)],
        cell_scale=chart._joint.cell_scale_A,minimum_distance=1.6,maximum_deformation=.25)
    return dict(source_job='28951227',source_peak_image=peak,
        source_archive_sha256=sha256(connected/'completed/strain_0100/native_terminal.tar.gz'),
        original_local_F_job='28812952', original_four_observables=True,
        source_to_original_chart_norm_A=float(np.linalg.norm(x)),
        source_retained_pattern_shift_A=(np.asarray(data['covectors_39_by_4']).T@x).tolist(),
        projected_same_observable_atomic_norm_A=float(np.linalg.norm(atomic)),
        projected_geometry=metrics, projected_atomic_coordinates=atomic.tolist(),
        projected_energy_or_gradient=None, projected_point_was_not_evaluated=True,
        same_atom_order_frame_and_continuous_lift=True,new_DFT_calls=0,
        whole_polygon_integer_offset=offset.astype(int).tolist(),
        registration_fixed_by_initial_and_verified_at_final_endpoint=True,
        evidence_not_a_local_candidate_or_stability_certificate=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();r=inspect()
    with args.output.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in r.items() if k!='projected_atomic_coordinates'}))
