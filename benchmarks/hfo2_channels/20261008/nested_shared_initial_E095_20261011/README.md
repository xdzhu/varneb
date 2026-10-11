# E095: nested bottleneck controls with one visible initial-well feature

This connects the existing B2--B5 measured-column controls to the existing
shared initial-energy response. It is an implementation delivery, **not**
a new HfO2 material forecast, accuracy claim or physical calculation matrix.

All four arms receive one immutable endpoint feature, the same bottleneck
columns and declared domain. Missing columns, incompatible atomic ablation,
wrong stationary/release assumptions and local-domain failures keep explicit
unavailable outputs. The same feature and preparation-cost reference remain
visible even for unavailable arms. No initial Hessian, derivative, index or
additional quadratic offset is invented. Actual dataset source IDs still
have to deduplicate shared preparation cost externally; four references to
one feature are not four sets of endpoint calculations.

The original physical-block oracle checks frozen, atomic and joint solves,
both rotated and unrotated charts, three parameters, affine anchor offsets,
signed negative gaps, nonlinear measured initial increments, missing columns,
B5 instability promotion, omitted gradients, clamped reactions, domain and
reference guards. Valid joint B4/B5 responses agree as required by the existing
harmonic null. These are synthetic implementation residuals, not material
prediction errors or proof of a connected DFT transition state.

27 new tests. The first targeted run reports178pass/1fail: an assertion
incorrectly expects the word `unmeasured` in the inherited reason string.
The corrected assertion checks the exact inherited reason plus its existing
unmeasured-curvature status; no production code or physical criterion changes
to accommodate that failure. The expanded local related run reports215pass,
0fail/0skip. Independent archive related tests also pass215/0fail/0skip;
the actual full archive suite passes2264/0fail/2existing skips,314warnings,
direct exit0 and355.117JUnitseconds. Core source bytes match the archive;
four new/updated text artifacts undergo Git LF/CRLF normalization with equal
decoded content. Both hashes are retained rather than claiming byte equality.
Only delivery/verification metadata changes after that test archive.
Clean-archive results and source identities are recorded in
[validation.json](validation.json), not inferred from the local tests.

At21:16:03CST both registered32CPU production children remainRUNNING:
28727368 (+1% PO-to-M) is at newstep7/.386236eV/Angstrom;
28727382 (zero-strain preserving flip) is at newstep8/.123621eV/Angstrom.
Both stderr files are empty, with no failure/terminal summary. The two
original +1% flip starters remainJobHeldUser. G2 ordinary completion is still
4/8. The source route itself calls no DFT or scheduler and reads no held-out
label. Registered workers continue their existing calculations.

The final native read at21:24:54CST confirms both stillRUNNING32CPU/stderr0:
M remainsnewstep7/.386236, preserving reachesnewstep9/.120476. These are
current optimizer residuals, not terminal SCF audits or convergence claims.

No extra initial/bottleneck Hessian, G3 choice, holdout, CI, install, release,
publication submission or production-input change is authorized by this
delivery. The full material/strong-reference model batch and prospective
freeze are still incomplete. The original five dirty user files remain
outside the delivery. See the [API documentation](../../../../docs/NESTED_RESPONSE_CONTROLS.md#apply-all-four-arms-with-one-shared-initial-well-energy-feature).

```sh
python -m pytest -q tests/test_nested_shared_initial_response.py tests/test_nested_response.py tests/test_shared_initial_response.py tests/test_stationary_gap.py tests/test_stationary_branch.py tests/test_response_cost.py tests/test_response_partition_invariance.py tests/test_nested_response_benchmark.py
```
