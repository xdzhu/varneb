# E254: complete released-center curvature; preserving admission refused

2026-10-16. The original job28868827 ended normally on node9/hfacnormal01,
one node/32 ranks,31120seconds,exit0:0. The unchanged E190 complete reader
reparsed **all272 native SCFs**, not just producer summaries. It checked
the original380-file runtime, binary, six physical assets, literal ordered
coordinates/reservations, native energy/full forces/stress and genuine MPI32.
The E252 first138 native source identities and its T result are unchanged.
No calculation was restarted or submitted to obtain this result.

## The scientific result

The centers are the existing fixed-Q0 local candidates associated with the
T escape and preserving flip paths, **not the endpoint phases**. Both meet
the previously specified released-gradient target0.03eV/A, but only T meets
the necessary raw-positive transverse-curvature condition at both steps.

| Center | Step/A | Minimum eigenvalue of symmetric raw R34/eV A^-2 | Negative count |
|---|---:|---:|---:|
| T escape candidate |0.01|0.164064528478|0|
| T escape candidate |0.02|0.157595936408|0|
| Preserving flip candidate |0.01|-0.022077212864|1|
| Preserving flip candidate |0.02|-0.024089732924|1|

Released-gradient norms are0.026299954459/0.019591084781eV/A(T/preserving).
The preserving candidate therefore **fails the preregistered necessary
sign condition** despite passing its gradient criterion. Its omitted
coordinates cannot simply be treated as a stable eliminated block. We do
not clip the eigenvalue, choose one step, use a pseudoinverse, relabel the
point a minimum, or start the independent tests on that assumption.

Both original nonsymmetric blocks, all full40-by34 gradient actions and
all136 energies/full gradients per center remain available. Their two-step
operator differences are0.054468926314/0.313107773422eV/A^2. These differences
and reciprocity/symmetrization diagnostics are **not total physical-error
bounds**. A raw negative sign is not a resolution-certified instability,
an index-one TS, a phonon frequency, or evidence that a physical branch
does not exist. T's positive signs do not certify its entire25-node sheet.

## Consequence for the registered tests

The E185/E188 forecasts,16 initial decks and16 reserved independent targets
remain unchanged and uncomputed. The physical+0.005 path holdout remains
sealed. This result does not complete final B0--B5 model/domain freezing or
demonstrate predictive superiority. In particular, stable B5=B4 and a full
T-basis re-expression is not an independent new algorithm.

Any later investigation of the preserving soft direction must retain this
first refusal, separate revised training from the original forecast and
respect the original modeling commitments/budgets and geometry guards.
There is no implicit permission to enlarge the domain, change starts in
the reserved test, tighten/retune electronic inputs, add a material or
rerun a whole path. The immediate result is a local validity limit, not a
successful independent accuracy test.

## Cost and evidence

Completed-call proxies are137.839362458148/137.425202945758core-hours;
sum275.264565403906. The final allocation cost is276.622222222222core-hours.
The producer's original completed.json keeps its null final-allocation
field; the complete native reader supplies actual sacct separately.
Original296 derivatives, training paths/faces, atomic treatments and other
preparation remain additional. This is not a speedup measurement.

- native_audit.json:a389154749a07109a3a91499343b86f06a2c185fd310611570637d71193e3381.
- preserving_curvature.json:588e6527936e2ceb8e94d8785293ae39fd36f561c0583004d0709a9fe1a561f2.
- completed.json:70d3e665c86ed3c73fb2f560ec5ba6902b594b40db47bcfd763e64fa21cfa9bb.
- T raw summary remains in../G3_T_released_complete_E252_20261015/T_curvature.json,
  SHA25666f7cdc9648afbfc3b61fc63a3835ec5868350c61aa921e56feea578a8498767.
- Unchanged native reader:../G3_released_audit_E190_20261015/audit.py,
  SHA25651d0452968394ff03d26efface7705ee31c7679215cb2204f6c69e13a5391d4f.

```text
python -m pytest -q tests/test_hfo2_complete_released_curvature.py
```

This offline replay checks complete coverage, all native source identities,
both full reconstructed actions/spectra, preserved negative signs, costs
and the failed necessary sign gate. It is not a second in-place native-log
audit or a physical stability/TS/error certificate. Original
ABACUS/PBE/100Ry/full10auDZP/Gamma2x2x2/P0/E0,ordinaryNEB0.10, guards and
all forecast/test budgets remain unchanged. No CI/install/release/monitor.

Related checks:74 local and74 independent, five modules/eleven new checks,
zero failures/errors/skips; not a full repository run. The previous E205
complete regression remains3510passed/two explicit skips. Tested source
tree d1b9e200ac599950c89da7247736c651b6a70ae8, archive SHA256
0b774c7311c3b7b6b46022c22168751d427dc52919fa888d79a61cdeabf9ab12.
Raw local/independent XML and validation.json preserve this verification.
Only test code is new; production physics/runtime and the manuscript are
unchanged. Scratch and independent exports are on C:, not E:; the five
pre-existing user file hashes are preserved. The full goal remains active
and incomplete: this failed local admission is not a substitute for its
required independent prediction evidence.
