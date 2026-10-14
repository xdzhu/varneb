"""Case-specific ordered-chart geometry diagnostic; no calculator calls."""
from pathlib import Path
import json
import hashlib
import sys
import zipfile

repo=Path(__file__).resolve().parents[4];sys.path.insert(0,str(repo))
import numpy as np
from ase.geometry import find_mic
from scripts import hfo2_G3_local_curvature as study

base=repo/'benchmarks/hfo2_channels/20261008'
manifest=json.loads((base/'G3_curvature_E150_20261014/manifest.json').read_bytes())
config=json.loads((base/'G3_complete_curvature_E159_20261014/preserving_config.json').read_bytes())
artifact=base/'G3_conditional_canary_E160_20261014/preserving_worker_outputs.zip'
assert hashlib.sha256(artifact.read_bytes()).hexdigest()=='c1955f42e1d9142235971a83eae0fc52c18b4c65521ec5e731e3ed1d2d4c28f9'
with zipfile.ZipFile(artifact) as archive:
    summary=json.loads(archive.read('preserving/canary.json'))
    point=summary['arms'][1]['point']
    rows=[json.loads(archive.read(n)) for n in archive.namelist() if n.startswith('preserving/curvature_seed/eval-') and n.endswith('/result.json')]
    row=next(r for r in rows if np.array_equal(r['coordinates'],point['coordinates_full_A']))
tchart=study.chart_for(manifest['centers']['T']);pchart=study.chart_for(manifest['centers']['preserving'])
basis=np.asarray(manifest['centers']['preserving']['directions'])[:,:36]
axes=np.asarray(config['internal_axes_36_by_2'])
x=np.asarray(point['coordinates_full_A']);g=basis.T@row['gradient'];qg=axes.T@g
atoms=pchart.to_atoms(x);source=pchart.reference_atoms;target=tchart.reference_atoms
other=base/'G3_conditional_canary_E160_20261014/T_worker_outputs.zip'
assert hashlib.sha256(other.read_bytes()).hexdigest()=='82fc05455912488051f52da2396606cbc060fb9be1c06ab42e07487768b2a44f'
with zipfile.ZipFile(other) as archive:
    tpoint=json.loads(archive.read('T/canary.json'))['arms'][1]['point']
t_relaxed=tchart.to_atoms(tpoint['coordinates_full_A'])
comparisons=[]
bt=np.asarray(manifest['centers']['T']['directions'])[:,:36]
for target_label,target_item in [('original_T_center',target),('conditional_T_candidate',t_relaxed)]:
    for label,item in [('P_original_center',source),('P_conditional_candidate',atoms)]:
        common=tchart.from_atoms(item)-tchart.from_atoms(target_item)
        dr=(item.get_scaled_positions(wrap=False)-target_item.get_scaled_positions(wrap=False))@target_item.cell.array
        mic,length=find_mic(dr,target_item.cell.array,pbc=target_item.pbc)
        comparisons.append(dict(structure=label,reference=target_label,
            same_atom_order_no_rotation_or_rematching=True,
            full_unwrapped_common_chart_norm_A=float(np.linalg.norm(common)),
            translation_free_common_chart_norm_A=float(np.linalg.norm(bt.T@common)),
            affine_referenced_ordered_MIC_atomic_norm_A=float(np.linalg.norm(mic)),
            affine_referenced_ordered_MIC_max_atom_A=float(np.max(length)),
            cell_difference_Frobenius_norm_A=float(np.linalg.norm(item.cell.array-target_item.cell.array))))
print(json.dumps(dict(channel='preserving',energy_change_eV_cell=point['energy_eV_cell']-manifest['centers']['preserving']['raw_source']['energy_eV_cell'],
    energy_difference_from_original_T_center_eV_cell=point['energy_eV_cell']-manifest['centers']['T']['raw_source']['energy_eV_cell'],
    energy_difference_from_T_conditional_candidate_eV_cell=point['energy_eV_cell']-tpoint['energy_eV_cell'],
    internal_displacement_norm_A=float(np.linalg.norm(x)),
    atomic_chart_displacement_norm_A=float(np.linalg.norm(x[:36])),
    open_cell_chart_displacement_norm_A=float(np.linalg.norm(x[36:39])),
    retained_gradient_components_eV_A=qg.tolist(),retained_gradient_norm_eV_A=float(np.linalg.norm(qg)),
    full_internal_gradient_norm_eV_A=float(np.linalg.norm(g)),
    released_gradient_norm_eV_A=float(np.linalg.norm(g-axes@qg)),comparisons=comparisons,
    same_phase_branch_TS_or_barrier_reduction_certified=False,new_DFT=0),indent=2))
