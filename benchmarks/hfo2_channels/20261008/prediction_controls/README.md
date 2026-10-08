# E026: preregistered simple controls become executable

This delivery implements B0 and both B1 controls from the existing prediction
protocol, with training-only representative selection and fresh hashed
prediction output. It supplies no new scientific novelty, DFT material result,
strain matrix, target error or model-advantage claim. B2--B5 material forecasts
and independent held-out labels are still missing.

The [API/CLI contract](../../../../docs/PREDICTION_CONTROLS.md) explains units,
endpoint input cost, unavailable features, unphysical abstention and why a
propagated training interval is not a prediction-error guarantee.

## Actual material gate check

Use the existing37-image free-cell G1 `network_update_20261009/analysis.json`,
SHA2567bd2dbf7ae54e13765a18a8a96e2c9ebb6f922f57fa226a4d499a35dd0eace1d.
The numeric core rejects it as incomplete; the HfO2 freezer rejects it as
the wrong mechanical/training family before output creation. No real material
forecast file is produced. This check does **not** constitute a prospective
prediction test or successful G1/G2/G3 closure.

Synthetic mathematical tests cover formula normalization/sign, training-only
selection when the upper condition changes the lowest channel, overlap/tie
rules, absence of features, total visible feature cost, no mutation, source
hashes, incomplete/contradictory training, mixed Hamiltonian/boundary, target
label fields, extrapolation, negative barriers and endpoint floors. Synthetic
fixtures are explicitly not strained HfO2 calculations.

```powershell
python -m pytest tests/test_prediction_controls.py tests/test_channel_competition.py -q
```

The original100Ry/full10auDZP/Gamma2x2x2/0.10/noCI contract and production
source are unchanged. G1 chains continue with the same two-active-chain cap;
no G2 or holdout job is submitted by this implementation. Clean staged-tree
regression, source hashes and an HF no-DFT negative-gate replay are recorded
in `validation_delivery.json` after verification.
