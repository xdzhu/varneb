"""Endpoint-gradient work between the two already observed training conditions.

The explicit connection pairs equal image indices in one ordered, unwrapped
chart. It is not a relaxed path in strain, a unique energy decomposition,
an envelope-theorem barrier derivative, or an independent prediction.
No intermediate geometry, held-out label or calculator call is generated.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from scripts.analyze_hfo2_clamped_reference import terminal_images
from scripts.audit_hfo2_clamped_network import audit as network_audit, IDENTITIES
from scripts.audit_hfo2_manuscript_work_cohort import BASE, CASES
from scripts.audit_hfo2_static_replica import sha256
from vcneb import BiaxialClampedModeCoordinates, clamped_plane_vcneb_boundary
from vcneb.strain_work import configuration_work


def cross_pair(chart, left, right):
    """Same-chart y and dual g; only the two supplied endpoints are evaluated."""
    y0,y1=(chart.from_atoms(a) for a in (left,right))
    values=[chart.energy_and_gradient(y,a) for y,a in ((y0,left),(y1,right))]
    delta=y1-y0; mean=.5*(values[0][1]+values[1][1]); n=3*chart.n_atoms
    translation=chart.rigid_translation_directions()[:n]/np.sqrt(chart.n_atoms)
    rigid=translation@(translation.T@delta[:n])
    parts=dict(internal_atomic=float(mean[:n]@(delta[:n]-rigid)),
        rigid_translation=float(mean[:n]@rigid),
        open_cell=float(mean[n:-1]@delta[n:-1]),
        prescribed_substrate=float(mean[-1]*delta[-1]))
    full=float(mean@delta)
    dh=right.cell.array-left.cell.array
    ds=right.get_scaled_positions(wrap=False)-left.get_scaled_positions(wrap=False)
    direct=[configuration_work(a.cell.array,a.get_stress(voigt=False,apply_constraint=False),dh,
        forces=a.get_forces(apply_constraint=False),fractional_direction=ds) for a in (left,right)]
    independent={k:.5*(direct[0][k]+direct[1][k]) for k in direct[0]}
    defects=[abs(full-independent['total_eV_per_control']),abs(sum(parts.values())-full),
        abs(parts['internal_atomic']+parts['rigid_translation']-independent['atomic_eV_per_control']),
        abs(parts['open_cell']+parts['prescribed_substrate']-independent['cell_eV_per_control'])]
    if max(defects)>1e-10:
        raise ValueError('independent native and dual-chart work contractions disagree')
    return dict(raw_energy_change_eV_cell=values[1][0]-values[0][0],
        endpoint_trapezoid_work_eV_cell=full,parts_eV_cell=parts,
        energy_minus_endpoint_work_eV_cell=values[1][0]-values[0][0]-full,
        independent_native_work=independent,maximum_identity_defect_eV_cell=max(defects),
        control_parameter_values=[float(y0[-1]/chart.strain_scale_A),float(y1[-1]/chart.strain_scale_A)],
        unwrapped_atomic_coordinate_change_norm_A=float(np.linalg.norm(delta[:n])),
        open_coordinate_change_norm_A=float(np.linalg.norm(delta[n:-1])),
        maximum_fractional_lift_change=float(np.max(np.abs(ds))),
        cross_condition_lift_change_above_half_cell=bool(np.max(np.abs(ds))>.5),
        geometry_rematched_wrapped_or_projected=False)


def analyze():
    network=network_audit()  # only the exact eight already admitted training terminals
    groups={}; provenance=[]; identities=[]
    for label,case,pin,_,_ in CASES:
        images,obs,_,receipt=terminal_images(BASE/case,pin)
        candidate,strain=IDENTITIES[receipt['job_id']]
        if Counter(images[0].get_chemical_symbols())!={'Hf':4,'O':8}:
            raise ValueError('same four-formula-unit training composition required')
        groups[candidate,strain]=(images,obs)
        provenance.append(dict(label=label,case=case,terminal_audit_sha256=pin,source_job_id=receipt['job_id']))
        identities.extend((r['raw_log_sha256'],r['input_sha256']['STRU']) for r in obs['raw_image_evaluations'])
    if len(groups)!=8:
        raise ValueError('exact eight unique admitted training sources required')
    reference,obs=groups['T',0.]
    boundary=clamped_plane_vcneb_boundary(12,np.array(obs['mechanical_boundary']['reference_cell_A']),allow_tilt=True)
    chart=BiaxialClampedModeCoordinates(reference[0],boundary=boundary,
        cell_scale_A=obs['cell_scale_A'],anchor_strain=0.,strain_scale_A=obs['cell_scale_A'])
    rows=[]
    for candidate in ('M','T','preserving','reversing'):
        zero,_=groups[candidate,0.]; plus,_=groups[candidate,.01]
        pairs=[dict(image_index=i,**cross_pair(chart,a,b)) for i,(a,b) in enumerate(zip(zero,plus))]
        for row in pairs:
            np.testing.assert_allclose(row['control_parameter_values'],[0.,.01],atol=1e-12,rtol=0)
            row['relative_PO_parts_meV_fu']={k:(v-pairs[0]['parts_eV_cell'][k])*250
                                           for k,v in row['parts_eV_cell'].items()}
            row['relative_PO_energy_change_meV_fu']=(row['raw_energy_change_eV_cell']-pairs[0]['raw_energy_change_eV_cell'])*250
            row['relative_PO_work_meV_fu']=(row['endpoint_trapezoid_work_eV_cell']-pairs[0]['endpoint_trapezoid_work_eV_cell'])*250
            row['relative_PO_energy_minus_work_meV_fu']=(row['energy_minus_endpoint_work_eV_cell']-pairs[0]['energy_minus_endpoint_work_eV_cell'])*250
        energies=[np.array([a.get_potential_energy() for a in chain]) for chain in (zero,plus)]
        peaks=[int(np.argmax(e)) for e in energies]
        observed=(max(energies[1])-energies[1][0]-max(energies[0])+energies[0][0])*250
        matched=pairs[peaks[0]] if peaks[0]==peaks[1] else None
        if matched is not None and abs(matched['relative_PO_energy_change_meV_fu']-observed)>1e-8:
            raise ValueError('matched maxima do not reproduce observed sampled peak change')
        rows.append(dict(candidate=candidate,sampled_peak_indices=peaks,
            matched_peak_pair=matched,observed_sampled_barrier_change_meV_fu=float(observed),pairs=pairs))
    # Native initial cache identities are already exact in network_audit.
    # Separately require the same unwrapped geometric lift, permitting only
    # binary64 serialization roundoff, never an integer image remapping.
    start_defects=[]
    for candidate in ('M','T','preserving','reversing'):
        for strain in (0.,.01):
            candidate_initial=groups[candidate,strain][0][0]
            reference_initial=groups['T',strain][0][0]
            defect=max(float(np.max(np.abs(candidate_initial.positions-reference_initial.positions))),
                float(np.max(np.abs(candidate_initial.cell.array-reference_initial.cell.array))))
            if defect>1e-12:
                raise ValueError('common cached initial state has a different unwrapped lift')
            start_defects.append(defect)
    starts=[r['pairs'][0]['endpoint_trapezoid_work_eV_cell'] for r in rows]
    if max(starts)-min(starts)>1e-12:
        raise ValueError('common initial endpoint work disagrees beyond algebraic roundoff')
    return dict(status='observed_two_condition_endpoint_work_not_stationary_barrier_response',
        script_sha256=sha256(Path(__file__)),physical_contract='original P0/E0/PBE100Ry/full10auDZP/Gamma2x2x2',
        source_receipts=provenance,source_frame_records_used=72,distinct_native_log_STRU_identities=len(set(identities)),
        maximum_initial_ordered_lift_serialization_difference_A=max(start_defects),
        maximum_common_initial_work_roundoff_eV_cell=max(starts)-min(starts),
        algebraic_identity_tolerances_are_not_DFT_error_bounds=True,
        raw_reader_passes=dict(network_admission=72,cross_strain_pairing=72),matched_image_pairs=36,
        mechanical_contract=chart.mechanical_contract,energy_unit='eV/Hf4O8 cell;250 converts to meV/HfO2',
        pressure_GPa=0.,reference='shared zero-strain PO+ ordered continuous lift',
        connection='straight in common fractional-atom/open-row/prescribed-strain coordinates; endpoints only',
        native_network_response=network['response'],candidates=rows,
        maximum_identity_defect_eV_cell=max(p['maximum_identity_defect_eV_cell'] for r in rows for p in r['pairs']),
        new_DFT_calls=0,scheduler_mutations=0,running_source_or_physical_inputs_changed=False,
        complete_six_physical_bytes_rechecked_locally=False,intermediate_geometry_generated=False,
        physical_holdout_generated_or_read=False,new_G3_mode_face_or_forecast_selected=False,
        barrier_error_bound_stationarity_or_TS_certified=False,
        limitations=['equal image index is a declared correspondence, not a canonical strain-dependent MEP parameterization',
            'endpoint quadrature of a nonstationary connection; not unique mode energies or an equilibrium barrier slope',
            'quadrature residual is not a DFT/barrier uncertainty bound or causal mode-strain coupling coefficient',
            'full Hessians, conditional branches, strong references and independent forecasts remain separate tests'])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('fresh analysis output required')
    result=analyze();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps([dict(candidate=r['candidate'],change=r['observed_sampled_barrier_change_meV_fu'],
        parts=r['matched_peak_pair']['relative_PO_parts_meV_fu'] if r['matched_peak_pair'] else None,
        residual=r['matched_peak_pair']['relative_PO_energy_minus_work_meV_fu'] if r['matched_peak_pair'] else None)
        for r in result['candidates']]))


if __name__=='__main__':main()
