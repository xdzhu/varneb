"""Training-only starts for the remaining eight preserving tests; no DFT.

This is a revised algorithm-defined candidate target, not a retrospective
success of the original E185/E188 conditional-plane forecast. Neither an
initial sign nor a small gradient certifies a stable conditional branch.
"""
from itertools import product
import numpy as np

TRAINING_Q = [(0., 0.), (-.1, 0.), (0., -.1), (0., .1), (.1, 0.)]


def targets():
    indices = sorted((p for p in product((-1, 0, 1), repeat=2) if p != (0, 0)),
                     key=lambda p: (abs(p[0])+abs(p[1]), p))
    return [dict(id=i+1, half_grid_index=list(p), q_A=[v*.05 for v in p])
            for i, p in enumerate(indices)]


def training_only_starts(rows, axes, valid):
    """Nearest admissible completed gradient candidate, integer Q tie then ID.

    Energy is not a selection feature. No independent result supplies a
    parent; failed parents are not repaired. ``valid`` is the unchanged
    physical geometry guard applied before any proposed calculator call.
    The external caller must first verify the complete native audit/hash.
    """
    axes = np.asarray(axes, dtype=float)
    if (axes.shape != (36, 2) or not np.isfinite(axes).all()
            or not np.allclose(axes.T@axes, np.eye(2), atol=2e-12, rtol=0)
            or not callable(valid)):
        raise ValueError('original orthonormal two axes and geometry guard required')
    if [row['id'] for row in rows] != list(range(10)):
        raise ValueError('all ten ordered training outcomes required')
    parents = []
    for i, row in enumerate(rows):
        if (row['Q_id'] != i//2 or row['q_A'] != list(TRAINING_Q[i//2])
                or row['start_sign'] != (1 if i % 2 == 0 else -1)):
            raise ValueError('literal original training Q/sign identity required')
        label = row['label']
        if label is None:
            if row['status'] != 'attempt_budget_exhausted_no_retry':
                raise ValueError('explicit original capped outcome required')
            continue
        x = np.asarray(label['coordinates_internal_A'], dtype=float)
        norm = label['released_gradient_norm_eV_A']
        directional = label['zeta_gradient_eV_A']
        if (x.shape != (36,) or not np.isfinite(x).all()
                or not np.allclose(axes.T@x, row['q_A'], atol=2e-12, rtol=0)
                or not np.isfinite(norm) or norm < 0 or not np.isfinite(directional)):
            raise ValueError('finite audited coordinates and physical dual gradients required')
        passes = norm <= .03 and abs(directional) <= 1e-4
        if (row['status'] == 'gradient_candidate_stability_unchecked') != passes or label['gradient_candidate'] != passes:
            raise ValueError('both original gradient stopping tests required')
        if passes:
            parents.append((row, x))
    result = []
    for target in targets():
        half = target['half_grid_index']
        ordered = sorted(parents, key=lambda pair: (
            sum((round(q/.05)-h)**2 for q, h in zip(pair[0]['q_A'], half)), pair[0]['id']))
        rejected = []
        for parent, x in ordered:
            start = x + axes@(np.asarray(target['q_A'])-axes.T@x)
            if valid(start):
                result.append(dict(**target, training_parent_id=parent['id'],
                    starting_internal_coordinates_A=start.tolist(),
                    rejected_parent_ids=rejected,
                    start_status='frozen_training_only_start_not_stability_certificate'))
                break
            rejected.append(parent['id'])
        else:
            result.append(dict(**target, training_parent_id=None,
                starting_internal_coordinates_A=None, rejected_parent_ids=rejected,
                start_status='no_original_guard_admissible_training_start_no_repair'))
    return result
