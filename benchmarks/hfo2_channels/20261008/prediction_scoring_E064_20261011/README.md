# E064: complete-panel, shared-reference forecast scoring

This is a bounded synthetic implementation experiment for the missing
prospective B0/B1/strong-reference B2–B5 scoring stage. It is **not** a
HfO2 held-out prediction, model freeze, new theory, or proven material
advantage. The actual G2 results and later G3/forecast/holdout gates remain
required. Production inputs, jobs and immutable runtime are untouched.

The module `vcneb.prediction_evaluation` preserves abstentions/signed raw
values and requires a complete registered model×case panel. It scores the
same deterministic reference interval for all models, not independently
fluctuating references; all pairwise and strongest-available baseline
comparisons are kept. Four Cmma example IDs remain four distinct models.
Missing strong baselines do not count as bad accuracy. Identical B4/B5 and
an already accurate B0 prevent a manufactured advantage claim.

See [the interface and limits](../../../../docs/PREDICTION_EVALUATION.md).
The strict full-panel flag is only arithmetic: upstream material audit,
independent-label timing and equal-information/cost still have to be proved.
The1e-12eV/f.u. computational tolerance is not DFT error or a new H1 δ.

`math_check.json` and `executed_initial_check.py` retain the first48random
interval check. `math_check_refined.json` additionally checks an adversarial
non-grid-aligned cusp: exact breakpoint evaluation and an independent10001-
point grid differ by7.4e-6eV/f.u., within the Lipschitz grid-spacing bound.
This distance is numerical discretization, **not** material uncertainty.
The synthetic panel contains nine models and two cases; all values, T/Cmma
labels and hash placeholders are illustrative, not actual material models.

Related local regression covers the new evaluator, registered B0/B1,
nested B2–B5 algebra, harmonic partition null and preparation-cost accounting.
The first and final local suites each have120passes/0failures/0errors/0skips.
The suite with actual historical-HF replay has122passes locally and122in
an independent Git-archive cwd, no failures/errors/skips.14new source/
synthetic-evidence files were byte-identical; later changes are docs/receipts,
not scientific source or observations. The isolated HF numeric check really
ran at02:10:05CST on Python3.10.9/NumPy1.26.4, with source/result hashes in
`HF_receipt.json`. Historical source copies are pinned, while current
numerical replay is checked separately to avoid freezing editable code bytes.
The HF wrapper loads
only the two pinned NumPy files in an explicitly isolated namespace: it
does not claim a full HF pytest or package-installation test.

```sh
python -m scripts.check_prediction_evaluation --output /new/scoring.json
python -m pytest tests/test_prediction_evaluation.py tests/test_prediction_controls.py \
  tests/test_nested_response.py tests/test_response_partition_invariance.py \
  tests/test_response_cost.py -q
```
