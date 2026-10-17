# One bounded zero-DFT model experiment and prospective P8 forecasts

2026-10-18. The [registration](protocol.md) was committed as4afab21c before
the fit, after all E260 training data and before every remaining P8 label.
This is not blind model design or a new DFT experiment. It reuses all81
native energy/full-gradient observations, including capped and rejected
trials as point data, **not** as conditional minima. Physical inputs and
every original refusal/forecast remain unchanged.

## Identifiability is not predictive superiority

At each of the two original curvature steps, fit the unmeasured Q/Q block
of a full-space harmonic reference (3coefficients). Compare with the exact
registered retained-direction nonlinear residual (109coefficients total).
Both keep the same source point, native energy/gradient, measured columns
and original positive R33 block. The raw negative R34 direction is retained,
not clipped or silently eliminated. Every fit and all five leave-one-Q-family-
out fits are full rank at the registered SVD cutoff.

Energies below are maximum absolute errors in meV/HfO2, not physical-error bounds.

| Model | Step/A | In-sample,81points | Central-family held out | Worst other held-out family |
|---|---:|---:|---:|---:|
| Full-space harmonic |0.01|3.6167|3.6167|0.0715|
| Retained nonlinear residual |0.01|0.1769|4.0836|0.1896|
| Full-space harmonic |0.02|3.5844|3.5844|0.0752|
| Retained nonlinear residual |0.02|0.1686|4.4101|0.1892|

The nonlinear family fits training energies better but does not establish
better transfer between these training families. Its maximum training full-
gradient error is still0.2124/0.2130eV/A, versus0.6152/0.6160for the harmonic
reference. Do not use energy-only fit quality to certify a force field or
stable conditional sheet. No new degree, weight, source step or training
cohort is selected in response to this result.

These group exclusions also change support: excluding the central family
leaves only|zeta|<=0.1000A in training but tests up to0.3688A. Excluding an
axial extreme removes that extreme Q from the remaining coordinate box.
Thus this is a shifted-support training diagnostic, not the independent
P8 interpolation test or proof of failure at its small-amplitude starts.
Both limitations and all five folds are retained in [training_fit.json](training_fit.json).
The full T-optical+cell re-expression of a complete harmonic model is the
same full-space reference, not a separate adaptive-method advantage.

## Prospective local forecasts, not a stable-PES certificate

The [P8 protocol](P8_PROTOCOL.md) retains the exact reserved eight Q positions,
the E264 training-only starts and the original twelve-call/dual-gradient
optimizer definition. All four families/steps are frozen without choosing
the winner from P8. Surrogate rollouts provide31 numerical predictions and
one explicit budget abstention (harmonic0.02, target8). Surrogate-oracle calls
are not DFT and reveal no independent label. The largest harmonic/nonlinear
differences on jointly predicted points are only of order0.01meV/HfO2; they
must not be advertised as a physically resolved improvement in advance.

[P8_forecasts.json](P8_forecasts.json) has SHA256
`d1d368372537a3d57d6a0cf38716ac8ea67210dbcac463d0dd2296a13e7b2046`.
It also preserves the exact five original standard controls and all14 T/Cmma
reference/step records. Agreement against this explicitly revised numerical
candidate target cannot retrospectively validate their old stable-plane
assumptions. Negative local energy changes are retained; these are **not**
negative activation barriers. Both raw and origin-aligned predictions use
the same fixed measured origin and no test-label rebiasing.

The numerical domain is only the exact eight prospective candidate targets,
their frozen starts, original0.5A total-radius/1.6A distance guard and same
physical contract. The **stable conditional predictive domain remains empty**;
there is no certified interpolation across the failed center. The original
2D refusal, failed candidates0/1/6 and unverified stability remain visible.
Local candidate-energy validation does not replace physical+0.005 path
validation, certify TS identity or demonstrate a JCTC method increment.

Actual P8 DFT requires source/config hash pinning and the original native
preflight before submission. It remains capped at96new attempts, two32MPI
allocations maximum, no retry/rollover or parameter change. T8 are already
seen. Physical+0.005 geometry/labels and private arrays remain sealed; the
complete physical B0--B5 barrier forecast is not frozen by this local panel.
All296+272derivatives,81trainingSCFs/85.6allocation core-hours and all other
preparations remain additional, not a zero-cost claim for the method.

```sh
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python fit.py --output NEW_TRAINING_REPORT.json
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python freeze_P8.py --output NEW_LOCAL_FORECASTS.json
```

These commands are run from this directory in a source checkout and never
invoke DFT/SSH. Existing outputs are refused. Code checks include analytic
gradients, coefficient recovery, rank refusal, preserved raw signs, full
old/new panel coverage and actual frozen-artifact replay; they are not new
independent material results.

Final verification:23related checks pass locally and in an independent
full-source Git export, three modules/eight new checks, no failures/errors/
skips; not fullsuite. Both actual fit and forecast CLIs were independently
rerun, not just their unit fixtures: output bytes match9b301edb... and
d1d36837...exactly. [validation.json](validation.json) records full data/XML/
tree/archive hashes and the externally pinned prospective forecast identity.
The local P8 panel is frozen; actual DFT deployment/preflight is the next
step, not more model tuning. The full physical-strain forecast and paper
goal remain incomplete.
