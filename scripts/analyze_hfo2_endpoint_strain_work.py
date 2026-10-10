"""Audit ten seen training endpoints and well-only B1 response implications.

No DFT, path-barrier selection, target-condition geometry or forecast freeze.
Stress/energy trapezoids are diagnostics, not measured error bounds.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.stress import voigt_6_to_full_3x3_stress

from scripts.prepare_hfo2_clamped_chains import endpoint_evidence, CHANNEL_FINAL
from scripts.audit_hfo2_static_replica import audited_results, sha256
from vcneb.strain_work import configuration_work


PHASES = ('PO_plus','T','M','PO_minus_T_preserving','PO_minus_T_reversing')
CONDITIONS = ('strain_0000','strain_p0100')
STEP = .01
FU = 4


def endpoint_root(case, condition, phase):
    if condition == 'strain_0000' and phase == 'PO_plus':
        return case/'clamped_PO_continuation_E046_20261009/completed_HF/endpoint'
    return case/'clamped_endpoint_matrix_20261009'/condition/phase/'completed_HF/endpoint'


def analyse(case):
    case = Path(case)
    endpoints, geometries = {}, {}
    for phase in PHASES:
        endpoints[phase], geometries[phase] = [], []
        for condition in CONDITIONS:
            _, atoms, boundary, manifest, provenance = endpoint_evidence(case,condition,phase)
            if manifest['strain'] != (0. if condition == CONDITIONS[0] else STEP):
                raise ValueError('only the two registered seen training conditions are used')
            root = endpoint_root(case,condition,phase)
            terminal = sorted((root/'calculator/image_0000').glob('scf_*'))[-1]
            raw = audited_results(terminal)
            summary = json.loads((root/'endpoint_relax_summary.json').read_text())
            sigma = voigt_6_to_full_3x3_stress(raw['stress'])
            dh = np.zeros((3,3)); dh[:2] = atoms.cell.array[:2]/(1+manifest['strain'])
            partial = configuration_work(atoms.cell.array,sigma,dh)['total_eV_per_control']/FU
            endpoints[phase].append(dict(condition=condition,strain=manifest['strain'],
                energy_eV_cell=float(raw['energy']),volume_A3=atoms.get_volume(),
                imposed_plane_partial_eV_fu_per_strain=partial,
                max_atomic_force_eV_A=float(np.linalg.norm(raw['forces'],axis=1).max()),
                open_traction_kbar=summary['open_traction_norm_kbar'],
                common_substrate_rows_A=atoms.cell.array[:2].tolist(),
                source_audit=provenance,terminal_export_relative=str(terminal.relative_to(case)).replace('\\','/'),
                terminal_call_audit_sha256=sha256(terminal/'call_audit.json'),
                terminal_raw_log_sha256=sha256(terminal/'OUT.ABACUS/running_scf.log'),
                full_pseudo_or_basis_bytes_rechecked_offline=False))
            geometries[phase].append((atoms,raw,sigma))
    # Check the same two substrate rows across all five representations.
    for i in range(2):
        plane = np.asarray(endpoints['PO_plus'][i]['common_substrate_rows_A'])
        if any(not np.allclose(plane,e[i]['common_substrate_rows_A'],rtol=0,atol=1e-10) for e in endpoints.values()):
            raise ValueError('endpoint ensemble mismatch')
    if not np.allclose(np.asarray(endpoints['PO_plus'][1]['common_substrate_rows_A']),
                       1.01*np.asarray(endpoints['PO_plus'][0]['common_substrate_rows_A']),rtol=0,atol=1e-10):
        raise ValueError('registered single biaxial intervention changed')
    pairs = []
    for phase in PHASES:
        lo,hi = endpoints[phase]
        shift = (hi['energy_eV_cell']-lo['energy_eV_cell'])/FU
        a0,a1 = [v[0] for v in geometries[phase]]
        dh = (a1.cell.array-a0.cell.array)/STEP
        ds = (a1.get_scaled_positions(wrap=False)-a0.get_scaled_positions(wrap=False))/STEP
        chord_terms = []
        for a,raw,sigma in geometries[phase]:
            full = configuration_work(a.cell.array,sigma,dh,forces=raw['forces'],fractional_direction=ds)
            release = dh.copy();release[:2] = 0
            open_cell = configuration_work(a.cell.array,sigma,release)['cell_eV_per_control']/FU
            chord_terms.append(dict(atomic_eV_fu_per_strain=full['atomic_eV_per_control']/FU,
                released_cell_eV_fu_per_strain=open_cell,total_eV_fu_per_strain=full['total_eV_per_control']/FU))
        partial_trap = STEP*(lo['imposed_plane_partial_eV_fu_per_strain']+hi['imposed_plane_partial_eV_fu_per_strain'])/2
        chord_trap = STEP*sum(c['total_eV_fu_per_strain'] for c in chord_terms)/2
        pairs.append(dict(phase_representation=phase,energy_shift_meV_fu=1000*shift,
            imposed_plane_trapezoid_meV_fu=1000*partial_trap,
            imposed_plane_trapezoid_defect_meV_fu=1000*(shift-partial_trap),
            configuration_chord_endpoint_terms=chord_terms,
            configuration_chord_trapezoid_meV_fu=1000*chord_trap,
            configuration_chord_defect_meV_fu=1000*(shift-chord_trap),
            finite_difference_secant_eV_fu_per_strain=shift/STEP,
            unmeasured_interior_work_or_quadrature_error_bound=None))
    shifts = {p['phase_representation']:p['energy_shift_meV_fu'] for p in pairs}
    initial_shift = shifts['PO_plus']
    controls = [dict(channel=channel,final_representation=phase,
                    B1_fixed_bottleneck_shift_meV_fu=-initial_shift,
                    B1_final_following_bottleneck_shift_meV_fu=shifts[phase]-initial_shift,
                    interpretation='seen-training endpoint-only null implication; absolute barrier requires a matched audited path')
                for channel,phase in CHANNEL_FINAL.items()]
    return dict(status='seen_training_endpoint_response_NOT_barrier_forecast',endpoints=endpoints,
        pairs=pairs,B1_training_response_implications=controls,formula_units=FU,pressure_eV_A3=0.,
        mechanical_family='same_substrate_tilt_open',training_strains=[0.,STEP],
        new_DFT_calls=0,physical_inputs_changed=False,holdout_generated_or_read=False,
        B0_or_B2_to_B5_predictions_computed=False,forecast_batch_frozen=False,
        full_relaxed_branch_derivative_or_Hessian_stability_certified=False,
        gradient_consistency_or_barrier_error_bound_established=False,
        limitations=['Stress is a local affine partial; released variables have nonzero residuals',
            'Endpoint chord uses the existing ordered seed/history lift without permutations; it is not an optimized intermediate branch',
            'Two-endpoint quadrature defects combine curvature and numerical/relaxation effects, not certified uncertainties',
            'T at+1% is the metric-lowered T descendant, not symmetry-restored P4_2/nmc',
            'Four endpoint-only B1 increments neither choose a channel nor prove H1/H2; full matched barriers remain required'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',type=Path,default=Path('benchmarks/hfo2_channels/20261008'))
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('new analysis file required; preserve earlier data')
    result=analyse(a.case)
    result['source_sha256']={ 'driver':sha256(Path(__file__)),
        'work':sha256(Path(__file__).resolve().parents[1]/'vcneb/strain_work.py')}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(status=result['status'],pairs=result['pairs'],new_DFT_calls=0)))


if __name__ == '__main__':main()
