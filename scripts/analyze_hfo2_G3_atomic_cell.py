"""Atomic/open-cell signatures of the existing complete local quadratics.

Offline arithmetic only: no new DFT, mode selection, physical holdout or
post-relaxation stability claim. Rayleigh terms are curvatures, not barrier
energy contributions. The prescribed substrate control is not eliminated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np

from scripts import summarize_hfo2_G3_curvature as measured

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'benchmarks/hfo2_channels/20261008'
ARCHIVE_SHA = '693fa998fc29448536384e6d69eda6f3a8d67d59eb033f6d4f2cb921d8fe341c'
PINS = {'T': 'a182e4c386889880c77c2069e939cd65665724d98d3838e251a077204cad1c3b',
        'preserving': '5d137fd0e676ad4fb7e87076d76b53e618ab9de25bc86228e02d4bcf0285a778'}


def raw_atomic_cell_signature(hessian, internal_basis, full_gradient):
    """Restrict a full internal quadratic to declared atomic/cell subspaces.

    Zero is an explicitly algebraic sign threshold, not a physical stability
    floor. Positive Hcc permits formal quadratic cell stationarity including
    its nonzero linear term, not certification of the nonlinear released block.
    Negative-subspace traces do not depend on eigenvector signs or rotations
    within that subspace. No individual vector is assigned phonon identity.
    """
    values = (hessian, internal_basis, full_gradient)
    if any(np.iscomplexobj(x) for x in values):
        raise ValueError('real declared Hessian/basis/gradient required')
    h,b,g = [np.array(x,dtype=float,copy=True) for x in values]
    if (h.shape != (36,36) or b.shape != (40,36) or g.shape != (40,)
            or not all(np.isfinite(x).all() for x in (h,b,g))
            or not np.allclose(h,h.T,atol=2e-12,rtol=0)
            or not np.allclose(b.T@b,np.eye(36),atol=2e-12,rtol=0) or np.any(b[-1])):
        raise ValueError('complete symmetric36D block, orthonormal internal embedding and excluded control required')
    _,s,vh = np.linalg.svd(b[:36],full_matrices=False)
    if not np.allclose(s,np.r_[np.ones(33),np.zeros(3)],atol=2e-12,rtol=0):
        raise ValueError('33 atomic/3 open-cell internal partition required')
    a,c = vh[:33].T,vh[33:].T
    if (not np.allclose(b[36:39]@a,0,atol=2e-12,rtol=0)
            or not np.allclose(b[:36]@c,0,atol=2e-12,rtol=0)):
        raise ValueError('atomic/cell subspaces cannot mix physical coordinate rows')
    haa,hcc,hac = a.T@h@a,c.T@h@c,a.T@h@c
    wa,wc = np.linalg.eigvalsh(haa),np.linalg.eigvalsh(hcc)
    w,v = np.linalg.eigh(h)
    neg = v[:,w<0]
    na,nc = a@(a.T@neg),c@(c.T@neg)
    atomic_trace=float(np.trace(na.T@h@na))
    cell_trace=float(np.trace(nc.T@h@nc))
    cross_trace=float(2*np.trace(na.T@h@nc))
    result=dict(atomic_dimension=33,open_cell_dimension=3,internal_dimension=36,
        internal_eigenvalues_eV_A2=w.tolist(),atomic_fixed_cell_eigenvalues_eV_A2=wa.tolist(),
        open_cell_eigenvalues_eV_A2=wc.tolist(),
        internal_raw_negative_count=int(np.count_nonzero(w<0)),
        atomic_raw_negative_count=int(np.count_nonzero(wa<0)),
        open_cell_raw_nonpositive_count=int(np.count_nonzero(wc<=0)),
        atomic_cell_coupling_operator_norm_eV_A2=float(np.linalg.norm(hac,ord=2)),
        negative_subspace=dict(dimension=neg.shape[1],
            open_cell_projector_trace=float(np.linalg.norm(nc,ord='fro')**2),
            atomic_Rayleigh_trace_eV_A2=atomic_trace,open_cell_Rayleigh_trace_eV_A2=cell_trace,
            cross_Rayleigh_trace_eV_A2=cross_trace,sum_negative_eigenvalues_eV_A2=float(w[w<0].sum())),
        formal_cell_stationary_quadratic=None,
        algebraic_zero_signs_not_physical_resolution=True,
        phonon_identity_stationarity_branch_stability_or_barrier_decomposition_certified=False)
    if not np.isclose(atomic_trace+cell_trace+cross_trace,w[w<0].sum(),atol=2e-12,rtol=0):
        raise ValueError('Rayleigh trace does not reconstruct measured negative spectrum')
    if np.all(wc>0):
        # No clipping/pseudoinverse, and no physical uncertainty claim from Hcc>0.
        np.linalg.cholesky(hcc)
        gi=b.T@g;ga,gc=a.T@gi,c.T@gi
        offset=-np.linalg.solve(hcc,gc)
        response=-np.linalg.solve(hcc,hac.T)
        schur=haa+hac@response
        sw=np.linalg.eigvalsh((schur+schur.T)/2)
        result['formal_cell_stationary_quadratic']=dict(
            atomic_curvature_eigenvalues_eV_A2=sw.tolist(),
            atomic_raw_negative_count=int(np.count_nonzero(sw<0)),
            open_cell_gradient_norm_eV_A=float(np.linalg.norm(gc)),
            nonzero_gradient_cell_offset_norm_A=float(np.linalg.norm(offset)),
            constant_energy_offset_eV_cell=float(.5*gc@offset),
            effective_atomic_gradient_norm_eV_A=float(np.linalg.norm(ga+hac@offset)),
            original_atomic_gradient_norm_eV_A=float(np.linalg.norm(ga)),
            eliminated_gradient_assumed_zero=False,
            raw_quadratic_only_not_DFT_prediction_or_measured_error_bound=True)
    return result


def analyze():
    path=BASE/'G3_curvature_E150_20261014/manifest.json'
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=measured.MANIFEST_SHA256:
        raise ValueError('original curvature manifest bytes differ')
    manifest=json.loads(raw)
    path=BASE/'G3_complete_curvature_E159_20261014/native_reports.zip'
    if hashlib.sha256(path.read_bytes()).hexdigest()!=ARCHIVE_SHA:
        raise ValueError('original complete native archive bytes differ')
    reports={}
    with zipfile.ZipFile(path) as archive:
        for channel,pin in PINS.items():
            raw=archive.read(channel+'_native.json')
            if hashlib.sha256(raw).hexdigest()!=pin:
                raise ValueError('original native report bytes differ')
            reports[channel]=json.loads(raw)
    measured.summarize(manifest,reports)
    rows=[]
    for channel in PINS:
        report=reports[channel]
        if report['completed_records']!=148 or len(report['blocks'])!=2:
            raise ValueError('complete two-step internal Hessians required')
        b=np.asarray(manifest['centers'][channel]['directions'])[:,:36]
        for block in report['blocks']:
            raw=np.asarray(block['raw_projected_eV_A2'])[:36,:36]
            rows.append(dict(channel=channel,step_A=block['step_A'],
                **raw_atomic_cell_signature((raw+raw.T)/2,b,report['raw_center_gradient_eV_A'])))
    return dict(status='offline_atomic_cell_raw_quadratic_signatures_not_stationary_prediction',rows=rows,
        units='joint length metric; curvature eV/Angstrom^2; not mass-weighted phonon frequencies',
        cell_scale_A=manifest['cell_scale_A'],internal_control_excluded=True,
        input_sha256=dict(manifest=measured.MANIFEST_SHA256,native_archive=ARCHIVE_SHA,**PINS),
        new_DFT_or_scheduler_operations=0,new_raw_SCF_reparse=False,
        source_modes_forces_ASR_or_calculator_parameters_changed=False,
        conditional_axes_final_representatives_or_forecasts_selected=False,
        physical_reserved_geometry_generated_or_read=False,
        limitations=[
            'Original nonstationary sampled centers, not post-relaxation candidate Hessians.',
            'Finite-step signs, step spread and subspace weights do not bound physical error.',
            'Atomic/cell weights and magnitudes use the fixed declared strain length metric.',
            'Rayleigh traces are directional curvature identities, not path barrier energy partitions.',
            'Formal cell stationarity retains the measured nonzero gradient; it is not a new DFT result.',
            'Local block coupling does not by itself explain the external-strain barrier slope.',
            'No extra full path, Hessian center, surface node, held-out label or new algorithm.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh offline analysis target required')
    result=analyze()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(status=result['status'],center_step_rows=len(result['rows']),new_DFT=0)))


if __name__=='__main__':main()
