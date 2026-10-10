"""Actual clamped pilot observations are not complete training networks."""
import argparse
import hashlib
import json
from pathlib import Path

from scripts.freeze_hfo2_prediction_controls import freeze


def run(root, output):
    output=Path(output)
    if output.exists():
        raise FileExistsError('refusing an existing negative-gate receipt')
    prior=Path(root)/'benchmarks/hfo2_channels/20261008/clamped_G2_E054_20261010/observations'
    paths=[prior/c/'step_0010/observation.json' for c in ('PO_flip_T_pattern_preserving','PO_to_M')]
    records=[json.loads(p.read_text()) for p in paths]
    if any(r['ordinary_residual_pass'] or r['replayed_fmax_eV_A']<=.10 for r in records):
        raise ValueError('registered first pilot observations no longer match this negative-gate case')
    prediction_namespace=output.parent/'prediction_must_not_exist'
    if prediction_namespace.exists():
        raise FileExistsError('a prior prediction namespace must not be overwritten')
    # Both sources are actual zero-strain pilot observations, deliberately
    # passed to the network API to check type/coverage refusal. Neither is
    # represented as a complete0/+1% network or a held-out label.
    try:
        freeze(*paths,prediction_namespace,holdout_path_labels_unread=True)
    except ValueError as error:
        reason=str(error)
    else:
        raise ValueError('unfinished observations were incorrectly accepted as training')
    if prediction_namespace.exists():
        raise ValueError('negative gate emitted predictions')
    result=dict(status='actual_pilot_observations_rejected_as_complete_training',reason=reason,
        input_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        ordinary_residual_pass=[r['ordinary_residual_pass'] for r in records],
        fmax_eV_A=[r['replayed_fmax_eV_A'] for r in records],
        complete_G2_training_network_provided=False,B2_to_B5_material_curvature_provided=False,
        material_predictions_emitted=False,held_out_structure_generated_or_labels_read=False,
        new_DFT_calls=0,source_files_modified=False,
        limitation='This is an actual incomplete-input/type refusal, not a material-prediction performance test')
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.repo_root.resolve(),args.output)
