# E163 — native completed-node audit before using a conditional surface

Scope: the existing two E161 batches, not a new experiment, full path, material
or Hessian center. The previous goal turn E162 was PROGRESS: the completed
local-response evidence entered the connected manuscript and was pushed in
`caa159da`. The full independent-prediction goal remains active and incomplete.

## What the auditor establishes

The standalone `scripts/audit_hfo2_G3_conditional_faces.py` reads the immutable
E161 runtime. Its source manifest (380 Python files), two exact configurations,
original accepted job handles and cached native origin must match. It does not
replace the live runtime, import a changed working checkout as that runtime,
instantiate a calculator, submit, resume or change electronic inputs.

At observation time it snapshots a consecutive prefix of complete point JSONs.
A not-yet-created or partially written point file stops that snapshot and is
reported as pending, not as a failed calculation. No unfinished node's final
trial is promoted. For each completed node, the auditor checks:

1. Fixed point order and Q, nearest safe successful parent and exact continued start.
2. Full evaluator/reference/mechanical contract and conservative twelve-attempt budget.
3. Reservation/coordinate fingerprints and one-to-one complete native records.
4. Original six physical inputs, true DSIZE=32, SCF convergence, ordered geometry,
   energy, forces, full stress and the derived 40-dimensional gradient.
5. Exact prescribed control and fixed substrate rows, internal embedding and
   the explicit final candidate's common released-gradient norm target.
6. Original-order cache-only optimizer replay: trials, cache hits, geometry
   screens, stopping reason and final coordinates. Missing/reordered data
   refuses replay without invoking DFT.

Budget/nonconvergence holes remain holes. Unsafe continuation has zero native
calls. Inflight failure artifacts, if any, are listed separately rather than
hidden by a successful earlier prefix. Native allocation accounting is a fresh
observation of the original handle; timeout or pending state is not silently
converted into a numerical failure or new submission.

`--mode complete` additionally requires the original job to be COMPLETED/0:0,
all 25 points and exact agreement with `face.json`. A prefix report does not
establish completion, released-block stability, unique/global branches, a
valid contour interpolation domain, a TS or held-out predictive accuracy.
The output retains energies, full candidate coordinates and physical gradients
for subsequent analysis without repeating DFT. Native logs/potential payloads
remain at their original HF paths; local report hashes preserve provenance.

## Verification scope

The new unit tests use explicitly injected toy records for arithmetic,
nonzero-Q replay and corruption rejection. They are not material SCFs.
The native execution and exact evidence hashes are retained in validation.json
and the native prefix reports. No physical +0.005 holdout geometry is generated
or read, no extra study allocation is committed, and no running source is edited.
The original curvature/canary/G2 preparation costs remain additional.

Related sources: [live batch protocol](../G3_conditional_faces_E161_20261015/README.md),
[completed origin audit](../G3_conditional_canary_E160_20261014/README.md),
[auditor](../../../../scripts/audit_hfo2_G3_conditional_faces.py),
[tests](../../../../tests/test_hfo2_G3_face_audit.py).

## Actual native results and reproducibility

Both native audits passed: T has nine audited completed points including its
cached origin, with 29 new SCFs reparsed; preserving has six including its
cached origin, with nineteen new SCFs reparsed. All thirteen new nonzero-Q
candidates replay exactly in the original numerical environment. Both original
handles were RUNNING during accounting observations; these reports do not
include SCFs of the current unfinished node. There were no recorded native
failure artifacts at these observations. No new SCF was performed by the audit.

The source was deployed as one standalone file outside the live study, from
the first tested tree `61c7c2e33cc6871a82434561fb2d16a4ce3b88ce`. Its actual
native byte SHA is `ab041335a93b8eaa98fa5cf8a812db5126931105cf972f586e4d15f8d06742da`.
The old 380-file runtime and original physical inputs remain unchanged.
Local arithmetic tests additionally replay the complete exported native
candidate gradients/coordinates and verify both report byte hashes; this is
not a second local reparse of remote raw SCF logs.

The final fourteen-module related suite passes 255 checks locally and 255 in
an independent final archive (zero failures/errors/skips), including 29 new
checks and the earlier manuscript/optimizer regression modules. These are
related checks, not a new entire-repository regression. The final archive tree
is `46243b2af577f32a99141ec7d9f8c09ec0fa26a5`; deployed auditor text is unchanged.
The native report JSONs, test XMLs and validation record here retain exact
evidence without repackaging or modifying the original SCFs.
