# E143: preserve conditional-grid branch evidence

The grid sampler currently retains selected coordinates and all-start cost,
but discards the per-point records holding alternative/unconverged starts.
The bounded implementation experiment is to preserve those records and bind
them into prospective finite-probe snapshots, without changing optimizers,
selection, energy/gradient evaluation or interpolation/stability gates.

Use only synthetic fixtures and related existing regressions here. Material
G2 remains6/8; no G3 centre/mode is selected, no physical+0.005 holdout is
generated or read, no extra DFT or scheduler writes, no live-source deployment.
100Ry/full10auDZP and the two original32MPI allocations stay unchanged.
Tests, archive replay and actual results follow execution, not this plan.

## Actual implementation and verification

The missing-grid-record regression first fails as expected on the old code.
The signed analytic fixture subsequently uses exact stationary wells to avoid
unrelated finite optimizer tolerance; this is not a before/after speed test.
One doc patch context rejects before any write; actual state is checked and
the exact context corrected. New sampler results now own complete, read-only
point/start final diagnostics. Every record must match its selected grid node
and all-start cost; missing or inconsistent records fail. No extra callback
is made. All-start final outcomes are not start geometries, full trajectories,
raw-SCF audits, stability or global-minimum certificates. An entirely failed
point still raises rather than generating a grid; raw backend failures stay
in their separate audit ledger.

New complete records bind the prospective training snapshot, including an
unselected outcome; tampering/removal is rejected. Four-corner patch screens
slice records correctly. Old empty-record snapshots retain their exact hashes;
old unavailable provenance is never invented. Twenty-eight new checks plus
related existing tests give132local/132independent passes,0fail/error/skip.
Full independently extracted tracked-source regression:2600pass/2skip/0fail,
314existing dependency/deprecation warnings,646.275s JUnit. Skips are optional
ase.io.abacus and the nonredistributed licensed VASP POTCAR; no installation.

Tested tree dbc9e92d59b3abb5dae164a238050411812655b4, archive SHA
aa42b04ea4b79b26541bb36aedc676c4493850f0db5236a95dc8c73b0d798bdc.
All6126tracked files exist before independent tests;7source texts normalized
equal and2original JUnits byte-equal. Only metadata/README/final copied JUnits
and ledger follow that tested snapshot; core implementation/tests unchanged.
Measured toy grid:4nodes/4points/12start outcomes,24actual and reported total
callbacks versus8selected-only callbacks. These are not material SCF/core-hours.

11:26:23CST both original jobs RUNNING32/hfacnormal01/stderr0: preserving
28812952 step8/.113387, reversing28764704 step28/.139222. G2 remains6/8,
no terminal/TS/holdout promotion, stop on rebound or hot source deployment.
The paper/PDF and original protected user changes remain unchanged. This
code-provenance milestone is not completion of the full research objective.
