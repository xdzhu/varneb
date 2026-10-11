# E088: shared initial-energy feature without an extra initial Hessian

This closes an implementation gap in the already registered optional endpoint
feature route. It is not a HfO2 forecast, independent prediction advantage,
new harmonic theorem, saddle certificate or completed G2/G3 result. The
observed G2 training remains four ordinary-converged candidates out of eight.

The public `stationary_gap_with_shared_initial_energy` pairs a bottleneck's
controlled stationary branch with one externally audited, immutable initial
energy-change record. Every comparison can reuse the same feature; no initial
Hessian is manufactured or requested. Keep its preparation cost and give the
same visible input to B0--B5. The numeric function itself invokes no calculator,
but visible endpoint preparation is not zero-DFT or blind path prediction.

It checks matched physical/boundary/energy/normalization declarations, P=0,
endpoint audit identities, recorded energy increments against raw totals,
positive endpoint preparation calls, index-one represented bottleneck and
both stationary-anchor/target domains. It preserves the bottleneck's anchor
correction, omitted movable gradients, clamped reactions and signed gap.
Measured initial anchor energies already contain their actual relaxation;
another initial quadratic correction would double-count it. Two endpoint
energies do not determine derivatives/curvature, so these remain absent.
Hash declarations do not audit endpoint geometry, phase or DFT uncertainty.

The new tests use synthetic quadratics, independent block solves, rotated
charts, distinct initial/bottleneck scales and three control conditions. They
also test a nonlinear measured initial increment, common large-energy shifts,
unchanged shared preparation cost, frozen omissions, negative gaps, malformed
metadata, incompatible contracts and index/domain/reference rejection.
There is no material/HfO2 local-Hessian accuracy claim.

## Reproduction

Use an output namespace that does not overwrite these receipts:

```powershell
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:PYTHONDONTWRITEBYTECODE = '1'
python -m pytest -q tests/test_shared_initial_response.py tests/test_stationary_gap.py tests/test_stationary_branch.py tests/test_biaxial_curvature.py tests/test_nested_response.py tests/test_response_cost.py --junitxml=/new/shared-initial-tests.xml
```

The first local focused run actually completed with 198 passes, zero failures,
errors or skips in 2.069 JUnit seconds, including 39 new tests. The actual
clean Git archive is independently tested; its final result is retained in
the delivery receipt rather than inferred from the local pass.

The first combined code/doc patch failed context verification before any
files were created or changed; the corrected patch then applied. This was
not a DFT/native production failure. No endpoint/G3/holdout data was generated
or read, no material basis or optional-feature choice was frozen, and no
extra independent chain, Hessian centre, scheduler write or software install
occurred. The original two live production allocations and fixed 100Ry/
full10-auDZP physical contract remain unchanged. Manuscript/PDF are not
rewritten as if this calculator-free implementation were a material result.

See [the API and formulas](../../../../docs/STATIONARY_GAP_RESPONSE.md) and
[the prospective protocol](../../../../docs/HFO2_PREDICTION_PROTOCOL.md).

## Final validation

The first clean full run had 2102 pass / 1 fail / 2 skip. The same failure was
reproduced in the actual previous HEAD archive: an old comparison test still
expected only two admitted training receipts after E085 added the third,
byte-pinned reversing receipt. The corrected test keeps **exact** three-hash
equality and verifies the real reversing receipt bytes; it does not relax
admission to arbitrary hashes, change production code or read a holdout.
The corrected local focused run passed 203 tests (2.862 JUnit seconds).
The final actual clean Git archive passed 2103 tests, with 2 explicit skips,
zero failures/errors and 314 existing warnings (353.870 JUnit seconds).
Both failed JUnit receipts and the successful final receipt remain here.
Skips are the unavailable optional `ase.io.abacus` reader and intentionally
undistributed licensed VASP POTCAR, not silently passed backend tests.

The NumPy-only `check_numeric.py` directly completed on the final canonical
source locally and on HF: 12 analytic points, maximum identity residual
2.78e-17 eV/f.u. and maximum cross-environment energy-field difference 2.78e-17.
Seven checker/critical-source hashes match byte-for-byte. NumPy1.23.5/1.26.4
are explicit; these residuals are not DFT uncertainty. HF has no pytest, so
no software was installed and no native full-pytest success is claimed.
The synthetic initial Hessian is used only as an independent known truth;
the new API requires no physical initial Hessian. Existing-output protection
was also executed: expected FileExistsError/exit1, original bytes preserved.
The checker and delivery receipts were added after the full-test archive;
core/test logic was unchanged after full regression.

At 17:41:17 CST, preserving28723655 was step34/.158213 and +1%M28709789 was
step3/.849185, both RUNNING32CPU/stderr0. The two +1% flipping starters remain
held. Recent timing puts the bounded segment ends around18:50--19:05 and
18:55--19:15; neither is a guarantee of ordinary convergence. There were no
scheduler mutations, admission retries, production-source/physical-input
writes, G3 choices, extra physical Hessians or+.005 holdout access.
The 19-page material manuscript/PDF remains unchanged. Full G2/G3, model
freezing and independent prediction evidence are still incomplete.
All hashes and direct execution results are in [validation.json](validation.json).
