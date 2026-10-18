"""Ordered endpoint oxygen roles; a geometric check, not a path-family label."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

REPOSITORY=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(REPOSITORY))
from ase.io import read
from ase.neighborlist import neighbor_list
from scripts.audit_hfo2_static_replica import sha256

CASE=Path(__file__).resolve().parent
BASE=CASE.parent
CUTOFFS=(2.4,2.5,2.6,2.7)


def roles(atoms):
    if atoms.get_chemical_symbols()!=['Hf']*4+['O']*8 or not atoms.pbc.all():
        raise ValueError('ordered periodic Hf4O8 required; no atom rematching')
    first,second,distance=neighbor_list('ijd',atoms,4.3)
    rows=[]
    for oxygen in range(4,12):
        ds=sorted(float(d) for i,j,d in zip(first,second,distance) if i==oxygen and j<4)
        if len(ds)<5:raise ValueError('five periodic Hf neighbors required for the distance-gap check')
        counts=[sum(d<cut for d in ds) for cut in CUTOFFS]
        rows.append(dict(ordered_atom_index_zero_based=oxygen,coordination_by_cutoff=counts,
            cutoff_robust_coordination=counts[0] if len(set(counts))==1 else None,
            third_fourth_fifth_Hf_distances_A=ds[2:5]))
    return rows


def analyze():
    register=BASE/'midpoint_reverse_E269_20261018/six_neighborhoods.json'
    config=json.loads((BASE/'matched_pattern_release_E270_20261018/prepared/config.json').read_bytes())
    if sha256(register)!=config['source_registration_sha256']:
        raise ValueError('original six-neighborhood source register changed')
    sources=[r for r in json.loads(register.read_bytes())['rows'] if r['channel']=='preserving']
    if [r['substrate_strain'] for r in sources]!=[0.,.005,.01]:
        raise ValueError('the three original preserving sources required')
    result=[]
    for source in sources:
        path=REPOSITORY/source['observation']
        if sha256(path)!=source['observation_sha256']:
            raise ValueError('original evaluated observation changed')
        observation=json.loads(path.read_bytes());ends=[]
        for i in (0,8):
            poscar=path.parent/f'POSCAR_{i:02d}'
            if sha256(poscar)!=observation['raw_image_evaluations'][i]['snapshot_POSCAR_sha256']:
                raise ValueError('ordered endpoint bytes changed')
            ends.append(dict(image=i,POSCAR_sha256=sha256(poscar),oxygen=roles(read(poscar,format='vasp'))))
        pairs=[(a['cutoff_robust_coordination'],b['cutoff_robust_coordination'])
               for a,b in zip(ends[0]['oxygen'],ends[1]['oxygen'])]
        exchanged=None if any(None in pair for pair in pairs) else all(set(pair)=={3,4} for pair in pairs)
        result.append(dict(substrate_strain=source['substrate_strain'],source_job=source['source_job'],
            source_observation_sha256=source['observation_sha256'],endpoints=ends,
            all_eight_ordered_oxygen_3_4_roles_exchanged=exchanged))
    return dict(status='original_preserving_endpoint_roles_checked',rows=result,
        geometric_Hf_O_cutoffs_A=list(CUTOFFS),periodic_neighbor_search_radius_A=4.3,
        atom_permutations_or_DFT_parameter_changes=False,new_DFT_calls=0,
        native_EFS_reparsed_in_this_geometry_only_analysis=False,
        literature_C_N_S_T_family_assigned=False,periodic_Hf_plane_crossing_classified=False,
        phonon_irrep_or_electronic_polarization_inferred=False,
        interpretation='preserving names a registered reference distortion, not retained O coordination roles')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=analyze()
    with args.output.open('x',encoding='utf-8',newline='\n') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(dict(status=result['status'],all_three_role_checks=[r['all_eight_ordered_oxygen_3_4_roles_exchanged'] for r in result['rows']],new_DFT_calls=0)))
