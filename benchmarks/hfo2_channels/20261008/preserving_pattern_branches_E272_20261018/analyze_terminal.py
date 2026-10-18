"""Local native replay of the whole finite signed batch, including failures."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tarfile
import tempfile

REPOSITORY = Path(__file__).resolve().parents[4]
sys.path.insert(0,str(REPOSITORY))
import numpy as np
from ase.io import read
from scripts.analyze_hfo2_local_accommodation import displacement_parts
from scripts.audit_hfo2_static_replica import sha256
from scripts.plot_hfo2_mechanism import original_chart_images, path_coordinates, rotated_t_triplet

CASE = Path(__file__).resolve().parent


def candidate_motion(chart, coordinates):
    base=CASE.parent
    config=json.loads((base/'matched_pattern_release_E270_20261018/prepared/config.json').read_bytes())
    t_path=base/'reference_variants/T.vasp'
    frame=base/'clamped_endpoint_seeds/clamped_seed_manifest.json'
    if sha256(t_path)!=config['T_reference_sha256'] or sha256(frame)!=config['orientation_sha256']:
        raise ValueError('the common original-T geometric metric changed')
    t=read(t_path); parent,patterns,_=rotated_t_triplet(t)
    orientation=json.loads(frame.read_bytes())['orientation']
    atoms=[chart.reference_atoms,chart.to_atoms(coordinates)]
    images=original_chart_images(atoms,orientation,t)
    q,polar,physical,_=path_coordinates(images,parent,patterns)
    defect=max(float(np.max(np.abs(q[1]-q[0]))),float(abs(polar[1,2]-polar[0,2])))
    if defect>1e-12: raise ValueError('the displayed physical four observables changed')
    parts=displacement_parts(physical['displacements_A'][0],physical['displacements_A'][1],patterns)
    f=[np.linalg.solve(t.cell.array,a.cell.array).T for a in images]
    return dict(atomic_motion_F_to_candidate=parts,measured_four_pattern_defect_A=defect,
                structural_polar_proxy_F_A=polar[0].tolist(),
                structural_polar_proxy_candidate_A=polar[1].tolist(),
                structural_polar_proxy_change_A=(polar[1]-polar[0]).tolist(),
                structural_proxy_is_not_electronic_polarization=True,
                open_y_stretch_change_percentage_points=100*(f[1][1,1]-f[0][1,1]))


def restore(destination):
    archive = CASE/'native_terminal.tar.gz'
    receipt = json.loads((CASE/'native_terminal_receipt.json').read_bytes())
    if (sha256(archive) != receipt['archive_sha256'] or receipt['job_id'] != '28948425'
            or not receipt['all_attempts_exported']):
        raise ValueError('complete original native archive required')
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            p=Path(member.name)
            if (p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0] != 'native_terminal'
                    or not (member.isfile() or member.isdir())
                    or not (destination/p).resolve().is_relative_to(destination.resolve())):
                raise ValueError('unsafe terminal native member')
        tar.extractall(destination)
    root=destination/'native_terminal'; files=json.loads((root/'files.json').read_bytes())
    if len(files) != receipt['exported_files']: raise ValueError('native file count differs')
    for name,digest in files.items():
        path=root/name
        if not path.resolve().is_relative_to(root.resolve()) or sha256(path) != digest:
            raise ValueError('native terminal bytes changed')
    return root,receipt


def analyze():
    spec=importlib.util.spec_from_file_location('E272_independent_terminal_reader',CASE/'analyze_prefix.py')
    reader=importlib.util.module_from_spec(spec); spec.loader.exec_module(reader)
    rows=[]; count=0
    with tempfile.TemporaryDirectory(prefix='varneb_P_terminal_') as tmp:
        root,receipt=restore(Path(tmp)); audit=json.loads((root/'audit.json').read_bytes())
        if audit['job_id'] != '28948425' or audit['bundle_pin'] != reader.BUNDLE_PIN:
            raise ValueError('original whole-batch handle and source required')
        for index in (0,2):
            chart,data,e0,_=reader.setup(index)
            q=np.asarray(data['orthonormal_39_by_4']); c=np.asarray(data['covectors_39_by_4'])
            t=np.zeros((39,3)); t[:36]=np.tile(np.eye(3),(12,1))/np.sqrt(12)
            for sign in ('minus','plus'):
                folder=root/str(index)/sign; summary=json.loads((folder/'treatments.json').read_bytes())
                records=[]
                for path in sorted(folder.glob('evaluations/eval-*/result.json')):
                    native=json.loads(path.read_bytes()); x=np.asarray(native['coordinates'])
                    _,energy,g,_=reader.replay_native(folder,chart,native['coordinates'])
                    rg=g[:36]-q[:36]@(q[:36].T@g[:36])-t[:36]@(t[:36].T@g[:36])
                    records.append(dict(coordinates=native['coordinates'],energy_eV_cell=energy,gradient=g,
                        coordinate_norm_A=float(np.linalg.norm(x)),atomic_released_gradient_eV_A=float(np.linalg.norm(rg))))
                if len(records) != summary['total_new_attempts']:
                    raise ValueError('all complete native attempts required')
                count+=len(records); arms=[]
                for arm,row in zip(('A','J'),summary['treatments']):
                    candidate=row['candidate']; energy=residual=None; motion=None
                    if candidate is not None:
                        matches=[r for r in records if r['coordinates']==candidate['coordinates']]
                        if len(matches)!=1: raise ValueError('one exact local candidate source required')
                        r=matches[0]; g=r['gradient']; x=np.asarray(r['coordinates']); n=36 if arm=='A' else 39
                        residual=float(np.linalg.norm(g[:n]-q[:n]@(q[:n].T@g[:n])-t[:n]@(t[:n].T@g[:n])))
                        energy=r['energy_eV_cell']
                        if (energy != candidate['energy_eV_cell'] or residual>.03 or np.max(np.abs(c.T@x))>1e-10
                                or not np.isclose(residual,candidate['released_gradient_norm_eV_A'],atol=2e-12,rtol=0)
                                or arm=='A' and np.any(x[36:])):
                            raise ValueError('candidate, held observables or released gradient differs')
                        motion=candidate_motion(chart,candidate['coordinates'])
                    arms.append(dict(arm=arm,status=row['status'],new_attempts=row['new_attempts'],
                        energy_eV_cell=energy,dE_minus_F_meV_fu=None if energy is None else 250*(energy-e0),
                        released_gradient_eV_A=residual,actual_geometry=motion))
                check=summary['known_direction_check']; secants=None
                if summary['treatments'][0]['candidate'] is None:
                    if check['status']!='unavailable_A_candidate' or check['new_attempts']!=0:
                        raise ValueError('an unavailable A cannot support a direction check')
                else:
                    points=[next(r for r in records if r['coordinates']==x) for x in check['coordinates']]
                    prep=json.loads((CASE/'prepared/seeds'/str(index)/'branch_preparation.json').read_bytes())
                    direction=np.asarray(prep['signed_starts'][1])/.1
                    xa=np.asarray(summary['treatments'][0]['candidate']['coordinates'])
                    for sign_probe,point in zip((-1,1),points):
                        if not np.allclose(point['coordinates'],xa+sign_probe*.01*direction,atol=2e-12,rtol=0):
                            raise ValueError('the registered checks must stay at A, not J')
                    kg=float((points[1]['gradient']-points[0]['gradient'])@direction/.02)
                    ke=float((points[1]['energy_eV_cell']-2*arms[0]['energy_eV_cell']+points[0]['energy_eV_cell'])/.01**2)
                    if (check['new_attempts']!=2 or check['step_A']!=.01 or check['whole_stability_certified']
                            or not np.isclose(kg,check['gradient_secant_eV_A2'],atol=2e-12,rtol=0)
                            or not np.isclose(ke,check['energy_secant_eV_A2'],atol=2e-12,rtol=0)):
                        raise ValueError('separate native directional secants differ')
                    secants=dict(gradient_eV_A2=kg,energy_eV_A2=ke,at='A',full_stability_certified=False)
                rows.append(dict(index=index,strain_percent=100*data['substrate_strain'],sign=sign,arms=arms,
                    directional_secants_at_A=secants,
                    total_new_attempts=len(records),last_attempt_coordinate_norm_A=records[-1]['coordinate_norm_A'],
                    last_attempt_atomic_released_gradient_eV_A=records[-1]['atomic_released_gradient_eV_A'],
                    last_attempt_energy_change_meV_fu=250*(records[-1]['energy_eV_cell']-e0),
                    last_attempt_is_not_automatically_a_candidate=True))
        if count != receipt['new_SCF_calls'] or count != audit['new_native_SCF_calls']:
            raise ValueError('whole terminal attempt accounting differs')
    return dict(status='all_original_signed_native_attempts_replayed',rows=rows,local_native_SCFs_reparsed=count,
        costs=dict(new_SCF_calls=count,native_proxy_core_hours=receipt['native_call_proxy_core_hours'],
                   allocated_core_hours=receipt['allocation_core_hours'],all_earlier_costs_additional=True),
        archive_sha256=receipt['archive_sha256'],new_DFT_calls=0,old_labels_and_scores_unchanged=True,
        full_orbital_pseudopotential_bytes_rechecked_locally=False,
        stationary_TS_or_full_stability_or_corrected_barrier_certified=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); result=analyze(); args.output.mkdir(parents=True,exist_ok=True)
    target=args.output/'analysis.json'
    if target.exists():
        previous=json.loads(target.read_bytes())
        if previous['archive_sha256']!=result['archive_sha256'] or previous['status']!=result['status']:
            raise ValueError('refusing to replace a different native-cohort analysis')
    with target.open('w',encoding='utf-8',newline='\n') as handle:
        json.dump(result,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(result))
