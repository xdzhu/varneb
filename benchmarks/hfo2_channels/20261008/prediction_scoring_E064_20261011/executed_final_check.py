"""Bounded synthetic forecast-scoring experiment; no material/holdout reads."""
import argparse
import json
from pathlib import Path

import numpy as np

from vcneb.prediction_evaluation import (
    AuditedBarrierReference, FrozenChannelForecast, evaluate_frozen_panel,
    shared_reference_gain_interval,
)


def check():
    rng = np.random.default_rng(61463)
    maximum_grid_extremum_distance = 0.
    for _ in range(48):
        lower, upper = sorted(rng.uniform(-.5, .5, 2))
        candidate = rng.uniform(-.6, .6)
        baselines = rng.uniform(-.6, .6, 4)
        lo, hi = shared_reference_gain_interval(candidate, baselines, (lower, upper))
        grid = np.linspace(lower, upper, 10001)
        gain = np.min(np.abs(baselines[:,None]-grid), axis=0)-np.abs(candidate-grid)
        tolerance = 2*(upper-lower)/10000+1e-14
        if not (lo-1e-14 <= gain.min() <= lo+tolerance and hi-tolerance <= gain.max() <= hi+1e-14):
            raise AssertionError('piecewise extrema and independent dense grid disagree')
        maximum_grid_extremum_distance = max(maximum_grid_extremum_distance,
            abs(lo-gain.min()), abs(hi-gain.max()))
    # An intentionally non-grid-aligned cusp makes the independent-grid
    # approximation visible; random examples often have flat extremal regions.
    cusp_candidate = .2000037
    cusp_baselines = [.1000037, .3000037]
    cusp_bounds = (.12, .28)
    cusp_range = shared_reference_gain_interval(cusp_candidate, cusp_baselines, cusp_bounds)
    cusp_grid = np.linspace(*cusp_bounds, 10001)
    cusp_gain = np.min(np.abs(np.array(cusp_baselines)[:,None]-cusp_grid), axis=0)-abs(cusp_candidate-cusp_grid)
    cusp_distance = cusp_range[1]-float(cusp_gain.max())
    if not (0 < cusp_distance <= 2*(cusp_bounds[1]-cusp_bounds[0])/10000+1e-14):
        raise AssertionError('independent grid must approach, not falsely resolve, the narrow cusp')
    references = [AuditedBarrierReference('synthetic_switch', .2, .005, 'a'*64, 'b'*64),
                  AuditedBarrierReference('synthetic_escape', .4, .005, 'a'*64, 'c'*64)]
    predictions = {'B5': [.2, .4], 'B0': [.201, .399],
                   'B2_T': [.24, .44], 'B2_Cmma_gauge1': [.23, .43],
                   'B2_Cmma_gauge2': [.24, .44], 'B2_Cmma_gauge3': [.25, .45],
                   'B2_Cmma_gauge4': [.26, .46], 'B3_atomic': [.23, .43],
                   'B4_joint': [.2, .4]}
    forecasts = [FrozenChannelForecast(model, ref.case_id, value, None, 'a'*64, format(i+1,'x')*64)
        for i, (model, values) in enumerate(predictions.items())
        for ref, value in zip(references, values)]
    report = evaluate_frozen_panel(forecasts, references, candidate_model='B5',
                                 baseline_models=list(predictions)[1:])
    if report['full_panel_strict_advantage']:
        raise AssertionError('same joint prediction and strong B0 must defeat a manufactured advantage')
    return dict(status='synthetic_scoring_check_passed_not_material_evidence',
        random_shared_reference_intervals_checked=48, independent_grid_points_per_interval=10001,
        maximum_dense_grid_extremum_distance_eV_fu=maximum_grid_extremum_distance,
        non_grid_aligned_cusp_checked=True,
        cusp_grid_extremum_distance_eV_fu=cusp_distance,
        grid_distance_is_discretization_not_DFT_uncertainty=True,
        equal_joint_and_candidate_predictions_and_accurate_B0_reject_advantage=True,
        example_complete_panel=report, new_DFT_calls=0, new_job_submissions=0,
        HfO2_training_or_holdout_read=False, B2_to_B5_material_models_frozen=False,
        independent_material_prediction_advantage_proven=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = check()
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'example_complete_panel'}))


if __name__ == '__main__':
    main()
