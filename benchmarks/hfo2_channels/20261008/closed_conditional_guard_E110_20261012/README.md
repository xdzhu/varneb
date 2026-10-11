# E110: geometry guard on fully prescribed conditional points

This is a local software-safety regression, not a new DFT result, a
convergence benchmark or a validated mode–strain prediction.

## Reproduced defect and bounded correction

`relax_orthogonal_at_q(..., optimizer="safeguarded_bfgs", trial_validator=...)`
previously bypassed the geometry validator when no orthogonal direction
remained. A fully prescribed point could therefore call the energy
evaluator even when the supplied guard rejected its geometry. This occurs
both when the two mode coordinates span the whole chart and when an
additional explicit frozen direction closes the remaining subspace.

Four parametrized regression cases first failed against the unpatched
implementation. They check both constructions, rejection before any energy
call, and guard-before-energy ordering for an accepted candidate. Synthetic
callbacks count calls; no first-principles calculation is used.

The correction applies the existing `valid(z)` check to the fully
prescribed branch before its static evaluation. For one supplied start,
an accepted candidate receives one evaluation and a rejected candidate
receives none. Existing multiple-start branch bookkeeping is not changed
or coalesced. The default L-BFGS-B path and nonempty orthogonal optimization
are unchanged. An empty projected gradient does not certify stability of
unmeasured or constrained physical directions.

## Reproduction and evidence

```sh
python -m pytest -q tests/test_mode_surface.py \
 tests/test_conditional_evidence.py tests/test_conditional_coverage.py \
 tests/test_clamped_mode_evaluator.py tests/test_mode_evaluator.py \
 tests/test_clamped_joint_curvature.py
```

The corrected focused run passes 116 tests. The receipt records the
initial four failing cases, this passing run, and the separately completed
clean-archive full regression, including exact file/tree/archive hashes.
The archived code is imported from the archive itself, not the dirty
working tree. One prematurely started full-suite attempt was explicitly
stopped because archive extraction completion had not yet been observed;
it is not counted as validation evidence. The final retry uses a distinct
JUnit filename and starts only after all 5766 tracked files are present
and the extraction process has ended.

See `validation.json` and the retained JUnit artifacts for actual results.
The bundle-local `.gitattributes` disables newline conversion only for
these raw JUnit records, preserving their measured SHA256 across checkouts.
The archive and large temporary outputs remain outside the repository on
E:/TEMP. This evidence bundle contains no pseudopotential or orbital files.

## Scientific and execution boundaries

The two existing HfO₂ training jobs run the original immutable E055
runtime; this local correction is not copied over their running source.
No electronic settings, optimizer options, material references, training
selection rules or registered budgets change. ABACUS/PBE/100 Ry/full
10 au DZP/Gamma 2×2×2/P0/E0 and ordinary 0.10 eV/Å/no CI remain fixed.
This turn initiates no DFT or scheduler mutation, although the existing
authorized production calculations continue to do DFT work.

The full G2 training matrix, subsequent model/curvature gates and
prospective independent holdouts remain incomplete. No G3 selection or
physical +0.5% holdout label is accessed. The scientific goal remains
active; this guard repair does not make the manuscript submission-ready.
