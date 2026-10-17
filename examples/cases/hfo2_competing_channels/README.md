# HfO2: replay the competing-channel paper results offline

From a source checkout with Python3.10 or later:

```sh
python -m examples.hfo2_paper_results --output hfo2_paper_tables
```

This command needs only Python's standard library. It does **not** start
ABACUS, connect to HF, access private author arrays or run a new calculation.
Choose a fresh output directory; existing files are never overwritten.
The archived reports are included in the repository, not the installed
library's minimal wheel. Run this example from the repository root.

The output contains:

- `channel_comparison.csv`: all four ordinary-converged channels at0/+1%
  substrate strain, including forward/reverse sampled peaks, reaction energy,
  force residual and original audit identities; energies are meV/HfO2.
- `T_local_model_errors.csv`: all19 frozen model arms on the original eight
  independent local points, retaining coverage and both raw and origin-aligned
  errors. This T label means the PO-to-T bottleneck neighborhood, not the T well.
- `T_local_all_predictions.csv`: all152 model/point records, including missing
  predictions and abstention reasons, with original predictions in eV/cell.
- `replay.json`: the same data and explicit source hashes, cost and claim limits.

The example verifies immutable report bytes and independently recomputes
profile maxima, forward-minus-reverse consistency, common initial energies,
error units and coverage. It does not reparse remote raw ABACUS logs; those
native audits remain linked through the frozen source reports. All original
Hessian/training preparation costs remain additional to the18 test SCFs.

## What these results say

Within the four examined channels, +1% strain lowers the high reversing-flip
sampled peak but raises the lower preserving-flip candidate. Looking only at
the former would give the wrong impression of the examined switching trend.
The T local tests also show why raw and origin-aligned errors must both be
reported: calibration changes the comparison with simple interpolation.
No independent quantitative advantage over the strongest full reference is
established by this export.

These are0K, P=0, same-substrate partial-clamping ABACUS/PBE results under the
archived100Ry/full10auDZP/original KPT contract. The ordinary0.10eV/A target
is not a certified barrier-error bound or a TS certificate. The two selected
physical+0.5% paths and the revised preserving-local tests are not included
before their own prospective protocols, computations and audits are complete.
This example must not be mistaken for a reproduction of finite-temperature
switching, device lifetime or a globally exhaustive channel search.

For new calculations use the [main quick start](../../../README.md) and
[fixed-substrate API](../../../docs/EPITAXIAL_BOUNDARY.md), with your reviewed
calculator and physical settings. The archived experiment's detailed source
index is in the [HfO2 blueprint](../../../docs/HFO2_PAPER_BLUEPRINT_2026-10-17.md);
those research identifiers are not required to run this offline example.
