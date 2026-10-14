# E170 — frozen matched material arm and native zero-DFT preflight

Previous E169 was PROGRESS, pushed asf0c9cce1: an opt-in feasible-ball solver
and four analytic controls were tested, not material acceleration. This case
freezes a real same-input comparison at original preserving point13,
Q=(-0.2,-0.1) Angstrom. It is training selected after a known failure, not an
independent forecast or a new Hessian center, phase, strain or full band.

## Matched experiment, not an extension of the live point

Use exactly the original requested start, continuation parent5, first actual
starting SCF, first-step full36D Hessian seed, 0.05-Angstrom initial trust
radius and0.03-eV/Angstrom full released-gradient target. The total internal
domain remains0.5 Angstrom, minimum separation1.6 Angstrom, with unchanged Q,
substrate and all six ABACUS/PBE/100 Ry/full10-au DZP/Gamma2x2x2 inputs.
The released coefficient radius is0.4472135954999579 Angstrom, derived from
the same center/metric and fixed-Q contribution. The original total-radius
guard remains authoritative; the numerical ball is not a mechanical reaction.

The old baseline is the already native-audited twelve-SCF capped trajectory,
not a newly computed arm. A lossless thirteen-file archive contains its point
record and all twelve original result JSONs. Each raw byte hash is linked to
the E168 native audit. Original source/ordered geometry/full gradients,
stress and32MPI metadata remain intact. The archive does not contain or
redistribute PP/orbital files or private author arrays.

`SharedExactCache` returns an old result only at the exact full coordinate
hash under the original physical contract. Gradient arrays are copied. A new
delegate must preserve that contract and fixed-control covectors; its own
conservative attempt budget still applies and exceptions are not retried.
The zero-DFT probe has no delegate: its first missing coordinate raises an
explicit exception instead of inventing E/g or launching a calculator.

The prospective arm reserves one additional preserving training slot,
conservatively counting the separate optimizer attempt at an already sampled
Q. Keeping E166's four-per-face reservation, preserving is30/81 committed,
51uncommitted; T remains29/81,52uncommitted. All8independent points/face and
physical+0.005 remain untouched. The old E16125-node/12-attempt policy is not
rewritten. A future fresh one-node32MPI/one-hour arm may use at most12new
attempts. Existing twelve-SCF/cache/curvature/path preparation is additional;
cached hits are not evidence of a net or same-accuracy speedup.

## Native preparation actually performed

Config SHA:f937a2e31ce30ec9b15852fa6a8e06cb73f4abbd4aa9b7f1f26bd29d561b4fd3.
The independent prepare CLI reproduces its bytes exactly. Original cache ZIP
SHA:02cf4f9b4bd169a775762411b798554c345be1496309b64b7ef6218c5ab99610.

An external fresh HF preflight uses the immutable E161 runtime and separately
pinned new solver/helper files. It rechecks all380 old Python files, existing
native center/origin EFS and the twelve original point13 SCFs. The original
trajectory replays with exactly its old order, guards and cap outcome.
The new solver then reuses exactly7 old cache calls, including the original
first SCF, before proposing a novel coordinate. The proposal remains at
total norm0.49999999999999994 Angstrom, fixed Q/control, minimum distance
2.0505811584 Angstrom, and passes actual STRU roundtrip plus all six original
asset hashes. Only one inert INPUT/KPT/STRU deck and original asset symlinks
are written. No new energy, force or stress is generated or inferred.
Native report SHA:2c22c66c58511c51d768271c17f0b6cadec7dfb6781af5106cfbd0822d5ce7c7.

The failed first cache-transfer draft used Windows path separators for remote
source paths; it created no SCF or scheduler state. Its scratch archive is
retained, and a fresh-r2 transfer uses explicit POSIX paths with all bytes
checked before creating the archive. An initial descriptive coordinate-hash
field name was corrected before freezing the tested final config; its draft
is preserved outside the repository. Neither issue changed material input or
was hidden as a backend failure/retry.

## Current boundary and next execution gate

This entry provides only `prepare` and `probe`, not a DFT run/dispatch mode.
The proposed arm is not submitted or computed. Before material execution,
the fresh runner/allocation/source/cost contract must also be frozen and tested.
The already prepared E166 eight statics have priority. Admission must verify
normal original completion and at most two committed study allocations,
including pending. No third job, physical parameter change, silent old-point
retry or last-trial candidate promotion is authorized by this preflight.

The full goal still needs actual matched material results, native complete
faces, released stability and strong-control independent tests. Native safety
and a genuinely novel safe proposal are not proof of convergence or benefit.

## Reproduction

From the repository root, the `matched` module's prepare/probe CLI accepts
explicit native-report, cache-ZIP, solver-file and config paths; existing
targets are refused. Use the tracked E168 report, this case's archive and
`vcneb/feasible_relaxation.py`. Native probes must use the original immutable
E161 source on PYTHONPATH and pinned new helper files, not replace that source.

```text
python -m pytest -q tests/test_hfo2_G3_feasible_material.py
```

Source/config/test/native receipts and exact provenance are retained in this
case's `validation.json`. API context is in the
[feasible conditional-relaxation manual](../../../../docs/FEASIBLE_CONDITIONAL_RELAXATION.md).

## Verification and resource observation

Before native preflight,233 related checks passed locally and233 in the
independent first archive. Final234 local/234 independent checks include16
new source/cache/physical-contract/probe checks across eleven related modules,
with zero failures/errors/skips. This is not a whole-repository regression
or proof of material benefit. The independently generated complete config
is byte-identical; the deployed helper and solver bytes are unchanged from
the first tested treecb9a8adb. Final tree70ec379d9afafbe44a296a3ac5d3af486c69747f
and exact source/archive/native/test hashes are recorded in `validation.json`.

The already prepared E166 controls' one-shot admission was actually checked.
Both original parents were RUNNING and the two allocations occupied, so
ready=false and zero scheduler writes resulted. The operational admission
driver remains outside live source; its observation and source SHA are kept
in `E166_admission_observation.json`. No third job was submitted, held or
queued. At03:57 the T face had24/25 candidates; preserving had17 processed
and15 candidates. They continue under their original policies; neither is
silently restarted or redefined as a completed stable PES.
