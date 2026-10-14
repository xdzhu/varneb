"""Zero-DFT strong-reference transport at the two frozen E150 centers.

Author arrays remain private inputs; only scalar identities/diagnostics are
exported. Frozen native center records are reused, not newly raw-SCF audited.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts import hfo2_G3_local_curvature as study
from scripts.analyze_hfo2_clamped_reference import REGISTRATION_SHA256, pinned_json
from scripts.analyze_hfo2_cmma_path_reference import (MASS, prepare_frame,
                                                     lowest_complete_doublet_indices)
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from vcneb import transport_atomic_reference_span
from vcneb.mode_subspaces import real_mode_subspace
from vcneb.phonons import identify_acoustic_modes
from vcneb.qe_modes import read_qe_gamma_modes

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/hfo2_channels/20261008'
MANIFEST_SHA='5f047f54ca83cb2dceb8cb71687eb9206a85c04205cdf26431565bc8a9b98a5a'


def verify(source_root):
    registered=pinned_json(BASE/'cmma_path_mapping/analysis.json',REGISTRATION_SHA256)
    pins={n:r['sha256'] for n,r in registered['author_source_audit']['sources'].items()}
    pins['Tetragonal.vasp']=registered['additional_T_anchor']['SHA256']
    if any(study.sha256(source_root/n)!=h for n,h in pins.items()):
        raise ValueError('original private author-source bytes differ')
    t_path=BASE/'reference_variants/T.vasp';gamma_path=BASE/'gamma_analysis/T_d0.01.npz'
    if (study.sha256(t_path)!=registered['production_T_SHA256']
            or study.sha256(gamma_path)!=registered['production_T_Gamma_SHA256']):
        raise ValueError('original T and Gamma preparation changed')
    t=read(t_path,format='vasp');_,patterns,_=rotated_t_triplet(t)
    triplet=real_mode_subspace(patterns.reshape(3,36).T*np.repeat(np.sqrt(MASS),3)[:,None],expected_rank=3)
    bases=[('T_geometric_triplet',triplet.basis_columns)]
    with np.load(gamma_path,allow_pickle=False) as gamma:
        if not np.array_equal(gamma['masses_amu'],MASS):raise ValueError('ordered source masses changed')
        acoustic,_=identify_acoustic_modes(gamma['eigenvectors_mass_weighted'],MASS)
        optical=np.setdiff1d(np.arange(36),acoustic)
        low=lowest_complete_doublet_indices(gamma['frequencies_thz'],optical)
        for label,indices in [('T_two_complete_optical_doublets',low),('T_all_optical',optical)]:
            if len(indices) not in (4,33):raise ValueError('complete optical source inventory differs')
            span=real_mode_subspace(gamma['eigenvectors_mass_weighted'][:,indices],expected_rank=len(indices))
            bases.append((label,span.basis_columns))
    frames=registered['selected_frame_registrations']
    if [f['id'] for f in frames]!=[103,106,114,127]:raise ValueError('all four original frames required')
    modes=read_qe_gamma_modes(source_root/'matdyn.modes',MASS)
    cmma=read(source_root/'Cmma.vasp',format='vasp')
    for frame in frames:
        _,basis,_=prepare_frame(modes,cmma,t,registered['author_source_audit']['QE_to_published_POSCAR'],frame)
        bases.append((f"Cmma_frame_{frame['id']}",basis))
    path=BASE/'G3_curvature_E150_20261014/manifest.json'
    m=pinned_json(path,MANIFEST_SHA);p=np.eye(3)[[2,0,1]]
    rows=[]
    for channel in ('T','preserving'):
        chart,_=study.validate(m,channel);record=m['centers'][channel]
        # Only the prescribed substrate rows equal T. The third row is open,
        # not an additional T-cell constraint or a reference to be repaired.
        if not np.allclose(np.asarray(record['substrate_reference_cell_A'])[:2],
                           (p@t.cell.array@p.T)[:2],atol=1e-10,rtol=0):
            raise ValueError('original fixed G2 cell/frame registration differs')
        measured=np.asarray(record['directions']);force=np.asarray(record['raw_source']['forces_eV_A'])
        g=-force.ravel()  # atomic contraction only; no invented cell/control gradient
        old_cell=p.T@chart.reference_cell@p;old_force=force@p
        old_dual=(-old_force@np.linalg.solve(t.cell.array,old_cell).T/np.sqrt(MASS)[:,None]).ravel()
        for label,source_basis in bases:
            result=transport_atomic_reference_span(t,source_basis,MASS,chart,source_to_target_cell_basis=p)
            identity=np.max(np.abs(result.raw_mapped_columns-result.local_basis@result.source_to_local
                                   -result.excluded_translation_columns))
            direct=source_basis.T@old_dual;transported=result.raw_mapped_columns[:36].T@g
            split=result.source_to_local.T@(result.local_basis[:36].T@g)+result.excluded_translation_columns[:36].T@g
            rows.append(dict(channel=channel,reference=label,rank=result.rank,
                expected_rank=source_basis.shape[1],
                reconstruction_max_defect=float(identity),
                independent_source_local_atomic_work_max_defect=float(np.max(np.abs(direct-transported))),
                internal_plus_translation_work_max_defect=float(np.max(np.abs(split-transported))),
                measured_span_projection_defect=float(np.linalg.norm(result.local_basis-measured@(measured.T@result.local_basis))),
                local_translation_overlap_norm=float(np.linalg.norm(chart.rigid_translation_directions().T@result.local_basis)),
                excluded_translation_column_Frobenius_norm=float(np.linalg.norm(result.excluded_translation_columns)),
                all_cell_and_control_components_exact_zero=bool(not np.any(result.local_basis[36:]))))
    return dict(status='registered_reference_transport_scalar_audit_not_model_or_curvature_prediction',rows=rows,
        source_author_sha256=pins,registration_sha256=REGISTRATION_SHA256,curvature_manifest_sha256=MANIFEST_SHA,
        work_covector_defect_unit='eV/(sqrt(amu)*Angstrom)',
        reconstruction_defect_unit='1/sqrt(amu)',
        production_T_SHA256=study.sha256(t_path),production_Gamma_SHA256=study.sha256(gamma_path),
        private_raw_or_transformed_author_arrays_exported=False,source_modes_or_forces_ASR_corrected=False,
        analysis_translation_gauge_restriction=True,source_frequencies_used_as_PBE_curvature=False,
        new_DFT_or_scheduler_writes=0,new_raw_SCF_reaudit=False,conditional_basis_or_forecast_selected=False,
        physical_reserved_geometry_generated_or_read=False,
        code_sha256={n:study.sha256(ROOT/n) for n in ('vcneb/reference_transport.py','scripts/verify_hfo2_reference_transport.py')})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise FileExistsError('fresh scalar audit target required')
    if args.output.resolve().is_relative_to(args.source_root.resolve()):
        raise ValueError('author source remains read-only; use a separate output directory')
    result=verify(args.source_root);args.output.parent.mkdir(parents=True,exist_ok=True)
    study.exclusive_json(args.output,result)
    print(json.dumps(dict(status=result['status'],reference_center_pairs=len(result['rows']),new_DFT=0)))


if __name__=='__main__':main()
