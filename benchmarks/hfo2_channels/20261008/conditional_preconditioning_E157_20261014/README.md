# E157: curvature-seeded conditional BFGS, not a DFT speedup result

The two E150 native batches continue with unchanged frozen source and
physical inputs. This milestone implements a missing optional optimization
route for the later bounded conditional relaxations; it does not start those
surfaces or choose their final coordinates, branches or forecast domains.

## Backend-independent implementation

`initial_hessian` is a full physical-coordinate approximate Hessian, projected
into the actual metric-orthonormal open space of the fixed-Q problem. Only
that open restriction must be numerically positive definite. Cholesky solves
initialize BFGS; invalid/nonrepresentable seeds are refused before physical
callbacks, without clipping or pseudoinversion. Fixed unstable directions
do not invalidate the optimizer seed. Default None keeps the identity seed.

Geometry guards, Armijo trials, trust radius, physical callbacks, all-start
records and final gradient gates are unchanged. Grid calls forward the same
seed. The existing conservative CalculatorModeEvaluator attempt budget is
reused, not duplicated or bypassed. A real evaluator/fake analytic-calculator
integration test confirms a one-attempt cap blocks the next factory call.
No new INPUT, cutoff, orbital, SCF, NEB threshold or release certificate.
User-facing convention and limitations: docs/CONDITIONAL_OPTIMIZATION.md.

## Actual quadratic replay of the already measured centers

The original E155 first-step native reports and E150 manifest remain pinned
by their exact hashes. Their measured36Dinternal H/g define exact quadratic
functions, not new ABACUS energies. For this **implementation replay only**,
fix the two negative spectral coordinates and minimize the positive complement.
Those spectral coordinates are not registered as physical G3 sampling axes.

Both arms use the same function, start, fixed Q, .05local-length step radius,
1e-8gradient tolerance and800iteration cap. The tight tolerance is for analytic
arithmetic checks, not a changed material NEB/DFT acceptance threshold.

| Quadratic replay | Identity callback count | Exact-H seed callback count |
|---|---:|---:|
| T | 129 | 5 |
| Preserving | 126 | 7 |

Both reach their independent positive-complement linear-solve oracle. These
counts include final verification callbacks, not SCFs or preparation calls.
Using the exact Hessian of a known quadratic naturally helps; the table is
not a novel algorithm, measured DFT acceleration, physical model superiority
or evidence of a stable finite conditional branch. Do not quote it as a
material speedup. Real DFT canary cost and branch checks are still required.

## Verification scope

Tests cover a nonunit metric, fixed negative directions, independent KKT
oracle, grid forwarding, unchanged None traces, invalid seed refusal before
callbacks, screened overly optimistic steps, conservative calculator caps,
and explicitly labeled measured-H quadratic replay. Original physical sources
and the sealed+.005condition remain untouched. No author-mode arrays added.

An initial whole-workspace pytest process was intentionally interrupted:
the working directory contains other sessions' untracked tests, so it is not
the clean tracked-code regression scope. No final JUnit was produced and
that partial run is not counted as passing or as a material failure. The
verified independent archive excludes those unrelated untracked files;
actual related/full test scope and source hashes are recorded after execution
in validation.json. No whole-suite claim before its terminal report exists.

The first clean archive run exposed two E155 arithmetic-replay tests that
incorrectly required bitwise-identical derived floating-point scalars across
BLAS thread configurations. A separate19-test diagnostic reproduced exactly
those two failures. The measured maximum derived-scalar difference is
3.907985046680551e-14; all nonfloat fields are exact. Compare derived arithmetic
with absolute2e-12/relative0 while retaining exact native-byte SHA pins,
metadata/types, failure guards and overwrite protection. This is an algebraic
test-portability policy, **not** a changed physical uncertainty or DFT gate.
Original reports are retained. A second fresh archive tests the corrected
test contract; no physical file is normalized or replaced to make it pass.

Final clean archive:2904passed/2expected skips/0failure/error,504.106JUnit
seconds, with314existing deprecation/compatibility warnings. The skips are
the optional local ase.io.abacus import and the deliberately unredistributed
licensed VASP POTCAR. New related scope:131local/131independent pass with no
skips. Full JUnit reports are losslessly zipped and their decompressed bytes
verified against the original XML hashes; no test evidence is discarded.
Six implementation/test/manual texts are normalized-equal to the tested
archive; physical source bytes are never normalized. Independent quadratic
counts, statuses and pins match exactly, while residual floats differ by at
most9.273325933063254e-15. The full scalar files are not claimed byte-identical.
