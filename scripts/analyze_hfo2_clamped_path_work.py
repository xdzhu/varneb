"""Audit signed native path work in every preregistered T/Cmma reference span.

Only already admitted ordinary-converged G2 training receipts are accepted.
No new DFT, author-code execution, model/reference selection or holdout access.
Scalar work only is exported; author mode arrays are not redistributed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
from ase.io import read

from scripts.analyze_hfo2_clamped_external_work import ADMITTED_TRAINING_AUDITS
from scripts.analyze_hfo2_clamped_reference import analyze as reference_analysis, terminal_images, original_chart_images, pinned_json
from scripts.analyze_hfo2_cmma_path_reference import MASS, prepare_frame
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.audit_hfo2_static_replica import sha256
from vcneb import ClampedModeCoordinates, clamped_plane_vcneb_boundary
from vcneb.continuous_projection import continuous_reference_coordinates
from vcneb.mode_subspaces import real_mode_subspace
from vcneb.path_work import sampled_path_work
from vcneb.qe_modes import read_qe_gamma_modes
from vcneb.strain_work import configuration_work


def partition(source_basis, masses, cell_dofs):
    """Same translation-free span restriction as the existing gradient audit."""
    basis = np.asarray(source_basis, dtype=float)
    roots = np.sqrt(masses)
    translations = np.zeros((3*len(masses), 3))
    for axis in range(3):
        translations[axis::3, axis] = roots/np.linalg.norm(roots)
    overlap = translations.T@basis
    restricted = basis-translations@overlap
    retained, singular, _ = np.linalg.svd(restricted, full_matrices=False)
    if np.count_nonzero(singular > 1e-10) != basis.shape[1]:
        raise ValueError("registered rank lost in the translation-free work chart")
    whole, _, _ = np.linalg.svd(np.column_stack((translations, retained)), full_matrices=True)
    omitted = whole[:,3+basis.shape[1]:]
    def padded(a):return np.pad(a, ((0,cell_dofs),(0,0)))
    cell = np.zeros((3*len(masses)+cell_dofs,cell_dofs))
    cell[-cell_dofs:] = np.eye(cell_dofs)
    return dict(retained_atomic=padded(retained), omitted_atomic=padded(omitted),
                rigid_translation=padded(translations), open_cell=cell), dict(
                    retained_rank=basis.shape[1], omitted_atomic_rank=omitted.shape[1],
                    source_translation_overlap_Frobenius_norm=float(np.linalg.norm(overlap)),
                    restricted_basis_min_singular_value=float(singular.min()),
                    source_modes_or_forces_ASR_corrected=False)


def path_arrays(images, original, t, rotation, report):
    """Dual common-T atoms plus native permitted-cell coordinates, no mixing."""
    reference_chart = continuous_reference_coordinates(original,t)
    atomic_x = (reference_chart['displacements_A']*np.sqrt(MASS)[None,:,None]).reshape(len(images),-1)
    boundary = clamped_plane_vcneb_boundary(12,np.array(report['mechanical_boundary']['reference_cell_A']),allow_tilt=True)
    chart = ClampedModeCoordinates(images[0],boundary=boundary,cell_scale_A=report['cell_scale_A'])
    native_y = np.array([chart.from_atoms(image) for image in images])
    native_g, energies, atom_g = [], [], []
    for native, old, y in zip(images,original,native_y):
        energy, gradient = chart.energy_and_gradient(y,native,pressure_eV_per_A3=0.)
        energies.append(energy); native_g.append(gradient)
        force = native.get_forces(apply_constraint=False)@rotation
        atom_g.append((-force@np.linalg.solve(t.cell.array,old.cell.array).T/np.sqrt(MASS)[:,None]).ravel())
    coordinates = np.column_stack((atomic_x,native_y[:,36:]))
    gradients = np.column_stack((atom_g,np.array(native_g)[:,36:]))
    direct_segments = []
    for left,right in zip(images[:-1],images[1:]):
        dh=right.cell.array-left.cell.array
        ds=right.get_scaled_positions(wrap=False)-left.get_scaled_positions(wrap=False)
        contractions=[configuration_work(image.cell.array,image.get_stress(voigt=False),dh,
            forces=image.get_forces(apply_constraint=False),fractional_direction=ds) for image in (left,right)]
        direct_segments.append({key:.5*(contractions[0][key]+contractions[1][key])
                                for key in ('cell_eV_per_control','atomic_eV_per_control','total_eV_per_control')})
    return coordinates,gradients,np.array(energies),direct_segments,reference_chart,chart


def analyze(case, source_root, registration, variants, gamma, *, terminal_audit_sha256,
            clamped_seed_manifest, clamped_seed_manifest_sha256):
    if terminal_audit_sha256 not in ADMITTED_TRAINING_AUDITS:
        raise ValueError("admitted completed training audit required before label access")
    coverage=reference_analysis(case,source_root,registration,variants,gamma,
        terminal_audit_sha256=terminal_audit_sha256,clamped_seed_manifest=clamped_seed_manifest,
        clamped_seed_manifest_sha256=clamped_seed_manifest_sha256)
    images,report,residual,receipt=terminal_images(case,terminal_audit_sha256)
    registered=pinned_json(registration,coverage['original_registration_sha256'])
    seeds=pinned_json(clamped_seed_manifest,clamped_seed_manifest_sha256)
    t=read(variants/'T.vasp',format='vasp')
    original=original_chart_images(images,seeds['orientation'],t)
    rotation=np.asarray(seeds['orientation']['cartesian_rotation_rows'])
    _,patterns,_=rotated_t_triplet(t)
    bases=[('T_geometric_triplet',real_mode_subspace(patterns.reshape(3,36).T*
        np.repeat(np.sqrt(MASS),3)[:,None],expected_rank=3).basis_columns)]
    with np.load(gamma,allow_pickle=False) as stored:
        bases.append(('T_two_lowest_complete_optical_doublets',stored['eigenvectors_mass_weighted'][:,
            coverage['T_lowest_two_complete_doublets_rank4_indices']]))
    modes=read_qe_gamma_modes(source_root/'matdyn.modes',MASS)
    cmma=read(source_root/'Cmma.vasp',format='vasp')
    for frame in registered['selected_frame_registrations']:
        _,basis,_=prepare_frame(modes,cmma,t,registered['author_source_audit']['QE_to_published_POSCAR'],frame)
        bases.append((f"Cmma_frame_{frame['id']}",basis))
    x,g,energy,direct,coordinate_audit,chart=path_arrays(images,original,t,rotation,report)
    spans={}
    for name,basis in bases:
        pieces,diagnostic=partition(basis,MASS,chart.cell_dofs)
        work=sampled_path_work(x,g,energy,pieces)
        independent=np.array([row['total_eV_per_control'] for row in direct])
        defect=float(np.max(np.abs(np.array(work['segment_full_work'])-independent)))
        cell_defect=float(np.max(np.abs(np.array(work['segment_partition_work']['open_cell'])-
            [row['cell_eV_per_control'] for row in direct])))
        if max(defect,cell_defect)>1e-10:raise ValueError('native direct and dual-coordinate work disagree')
        spans[name]=dict(**diagnostic,**work,independent_native_work_max_defect_eV_cell=defect,
                        independent_cell_work_max_defect_eV_cell=cell_defect)
    return dict(status='sampled_signed_path_work_not_mode_energies_or_forecast',
        script_sha256=sha256(Path(__file__)),core_source_sha256=sha256(Path(__file__).resolve().parents[1]/'vcneb/path_work.py'),
        source_job_id=receipt['job_id'],terminal_step=receipt['terminal_step'],terminal_audit_sha256=terminal_audit_sha256,
        strain_training_condition=coverage['strain_training_condition'],ordinary_fmax_eV_A=residual['fmax_eV_A'],
        physical_contract='P0/100Ry/full10auDZP/PBE/Gamma2x2x2; original six inputs and clamped tilt-open boundary',
        formula_units=4,energy_and_work_unit='eV per Hf4O8 cell; multiply by250 for meV/HfO2',
        original_registration_sha256=coverage['original_registration_sha256'],author_source_sha256=coverage['author_source_sha256'],
        clamped_seed_manifest_sha256=clamped_seed_manifest_sha256,
        fixed_T_integer_gauge=coordinate_audit['fixed_integer_lattice_shifts_by_atom'].tolist(),
        native_cell_chart_reference='own matched initial PO+ cell, 3 original permitted open deformations',
        atomic_chart='same original T sqrt(amu)*Angstrom coordinates and their dual covectors; no mode-cell mixing',
        absolute_energy_eV_cell=energy.tolist(),sampled_peak_image=int(np.argmax(energy)),
        independent_native_segment_work=direct,representations=spans,
        all_E023_Cmma_frames_retained=True,portable_native_frames_reparsed=9,
        complete_six_physical_bytes_rechecked_locally=False,new_DFT_calls=0,scheduler_mutations=0,
        G3_selection_or_holdout_access=False,production_geometry_source_or_physical_inputs_changed=False,
        limitations=['not separate state-function mode energies, mode coupling coefficients or model accuracy',
                     'coarse quadrature residual is not a measured DFT/barrier uncertainty bound',
                     'reference-dependent signed work, including omitted modes and measured open-cell response',
                     'no spring/NEB-force integration, continuous TS, conditional branch or unobserved mechanism claim'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('case','source-root','registration','variants','gamma','clamped-seed-manifest','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--terminal-audit-sha256',required=True)
    parser.add_argument('--clamped-seed-manifest-sha256',required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('fresh path-work output required')
    result=analyze(args.case,args.source_root,args.registration,args.variants,args.gamma,
        terminal_audit_sha256=args.terminal_audit_sha256,clamped_seed_manifest=args.clamped_seed_manifest,
        clamped_seed_manifest_sha256=args.clamped_seed_manifest_sha256)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(job=result['source_job_id'],references=list(result['representations']),new_DFT_calls=0)))


if __name__=='__main__':main()
