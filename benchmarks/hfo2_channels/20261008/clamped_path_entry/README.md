# E037: an explicit clamped-plane production entry, not G2 results

2026-10-09. E036 checked the endpoint delivery, but the generic material
runner still had only free-cell and fixed-cell choices. Reusing its G1
continuation invocation for G2 would have released the prescribed substrate.
This is an actual missing integration, not another scientific hypothesis.

The packaged runner and public config now select the existing exact
`ClampedPlaneBoundary` only when the reference, tilt choice, common cell scale,
unaligned frame and linear interpolation are explicit. No physical calculator
parameters are inferred or changed. Stress remains required. Raw endpoints
and supplied internal images are rejected before any projection or calculator
symbol loading; candidate checks combine the same plane with the declared
geometry limits and existing per-step cell trust radius. Global mode artifacts
are not silently composed with this boundary. Different boundary provenance
cannot overwrite an existing workdir's preflight.

See [the interface and config fragment](../../../../docs/EPITAXIAL_BOUNDARY.md).
`tests/test_clamped_plane_cli.py` covers tilt-released/normal-only oblique
planes, incompatible raw/resume images, strict options, public prepare/run,
serial/threaded few-step EMT and fixed endpoints evaluated once. The mock
factory confirms backend parameters and command pass through unchanged.
These tests do not measure HfO2 barriers or certify a transition state.

## Tested delivery

Only owned staged source was archived, excluding five pre-existing dirty
user files and unrelated untracked work. The final tested tree is
`6b30eb2a0b25df1e14a156dd4568b9343c6cb1b5`; its22,736,722-byte archive has
SHA256 `675bcce0fdcd3aa33cc5eba6023450ceb1ac1be0bce4b3882b808cc0964e813d`.
Full clean regression: **1197passed,2skipped,0failures/errors**,146.56s;
313existing dependency/deprecation warnings are retained. This is a new full
suite, not E032's historic pass count. An earlier staged version also passed
the suite; final delivery additionally uses the explicit scale in geometry
diagnostics and guards. Both receipts are identified separately.

Initial focused tests exposed two machine-roundoff assertions
(4.4--6.7e-16 A after a fractional-coordinate round trip) and a test-only
wrong seed-manifest key. The assertions now check roundoff at1e-14 A and the
fixture follows `manifest_file`/`seed_file`. No mechanical/SCF/NEB threshold
was relaxed. The focused set passed99tests; after shared-scale diagnostic
propagation, all its tests are covered again by the final clean full suite.

The same final archive actually ran the helper once on HF ICU/ASE3.23.1b1,
NumPy1.26.4, in a fresh namespace. All10 registered unrelaxed seeds pass
the three-image geometry entry, and all10 deliberately incompatible internal
images are rejected before calculator loading. Three open cell directions,
source hashes, reference scale and seed bytes are retained. There are **zero
calculator symbols loaded, zeroDFT calls and zerojob submissions**. This
does not establish an interpolated seed as a converged/material-valid path.

```bash
python benchmarks/hfo2_channels/20261008/clamped_path_entry/check_entry.py \
  --source /path/to/immutable-tested-source \
  --output /path/to/fresh/separate-geometry-probe
```

`preflight_HF.json` is the recorded actual stdout JSON, not a fabricated DFT
result. Module/hash values match the tested extracted source. This code is
not injected into the live R12 production source. E036's older endpoint-only
archive is historical and unchanged, not relabeled as this new path entry.

## Actual computation remains gated

At12:40:41CST the same hfacnormal01/32CPU jobs are still RUNNING:
preserving28319570/node11 has complete log step68/.104515eV/A;
reversing28319571/node26 step30/.220242. Neither passes ordinary0.10 yet.
The production archive remains SHA256df4c12ee... and there are at most two
active study chains. No cancellation/restart, CI, parameter change or early
G2/holdout submission occurred. A rebound is not a stop rule.

G1 still requires the final source/SCF/phase/variant and sampling audit.
Only then may the registered four-step clamped PO+ BFGS canary start. The
finite G2 matrix and G3 matched-boundary response/independent predictions
remain required. This fixes delivery correctness; it is not an innovation,
a material prediction, a JCTC-ready verdict or completion of the full goal.
See [the validation receipt](validation_delivery.json).

E038 follow-up closes the reverse reuse direction too: an existing clamped
directory must not be rerun as free-cell merely by omitting clamped flags.
Both directions pass the new local tests; final clean suite1199pass/2skip,
137.82s. An actual HF attempt to use the earlier geometry-only clamped
directory as free-cell raises the expected `FileExistsError` before factory
loading, and its preflight SHA remains identical. The trace is a successful
safety test, not a production DFT failure. E037's source/archive/reports are
historical and unmodified. See [the additional receipt](boundary_guard_delivery.json).
