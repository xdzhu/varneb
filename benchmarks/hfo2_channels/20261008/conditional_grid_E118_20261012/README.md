# E118: guarded conditional grids and all-start callback accounting

Two actual wrapper defects are reproduced in before_fix.xml: grid sampling
cannot accept the existing single-point safeguarded optimizer/geometry guard;
and a four-point, three-start prescribed grid executes12callbacks but reports
only4because it counts selected starts alone. This is a software integration
and accounting fix, not a new physical method or HfO2 response result.

The grid forwards optimizer, trial validator, trust radius and frozen-start
policy unchanged. Legacy L-BFGS-B defaults remain. Its n_evaluations reports
all callback calls, including unconverged starts. The single-point selected
count is preserved; total_n_evaluations sums complete outcomes, while legacy
multi-start records without them remain unknown. Backend exact-cache hits
are callbacks, not new SCFs; actual source/cost ledgers remain necessary.

Tests cover pre-evaluation rejection of a fully prescribed point, nonlinear
signed competing starts, invalid-option rejection before any callback, failed
start costs, and the existing clamped joint ASE/EMT evaluator/cache flow.
The nonlinear toy and Cu/EMT tests are implementation evidence, not HfO2
conditional surfaces, stationary TSs, DFT error bounds or independent forecasts.
No material Q coordinates, axes, branches, tolerances or new probes are selected.

```sh
python -m pytest tests/test_conditional_grid_guards_and_cost.py tests/test_mode_surface.py tests/test_clamped_mode_evaluator.py
```

Initial focused tests pass62. The independent full suite now passes2518,
with2skips,0failures/errors and314existing/dependency warnings. Skips are the
unavailable optional ase.io.abacus reader and licensed VASP POTCAR omitted
from clean distribution; neither is disguised as a pass. Actual JUnit time
is486.829s. All5892tracked files are present and changed code/test/doc sources
compare normalized-equal before tests. Only evidence metadata and ledgers
are added afterward. No fresh local full suite or HFpytest is claimed.

Validation.json pins source, archive and raw test evidence. Two rejected
patch-context checks perform no writes before the successful patch. No original
SCF results, manuscript numbers or plots change. The two registered HfO2 chains keep
their frozen remote runtime and original100Ry/full10auDZP input contract;
the local library change is not deployed over a running source checkout.
G2 must finish before G3 material selection or physical holdout access.
Native05:40:19CST corroborates both RUNNING32/hfacnormal01/stderr0, steps
24/.159496 and4/.127876. G2 remains6/8. No cancel/resubmit, new HfO2 DFT call,
deployment, install, monitor, CI, release or submission. User edits remain
untouched/unstaged; goal active and scientifically incomplete.
