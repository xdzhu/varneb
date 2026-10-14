# E159: complete native two-amplitude curvature and bounded DFT canaries

The two original E150 calculations have completed normally: T job28855614
used4:33:37 on node10 and preserving job28855615 used4:20:24 on node27,
both32CPUs, Slurm COMPLETED/0:0. Each supplied148static SCFs, not148NEB
iterations:37registered directions (36internal plus one separate prescribed
control), two signs and two amplitudes (.01/.02Angstrom).

The original E152 reader, run against the unchanged E150 runtime, reparsed
all296native SCFs. Ordered geometry, electronic convergence, genuine DSIZE32,
energy/forces/stress, gradient replay, original six physical input hashes,
center and binary identities passed; no failure artifacts were found.
This completes the original two-center derivative budget, not a third center,
new material or new full NEB band. Earlier E155prefixes overlap this dataset.

## What the completed data actually establish

Both amplitudes at both nonstationary centers give two raw negative internal
eigenvalues. Their lowest three eigenvalues (eV/Angstrom²) are:

| Center | Amplitude/Angstrom | First | Second | Third |
|---|---:|---:|---:|---:|
| T | .01 | -.27619343 | -.16464809 | .05321565 |
| T | .02 | -.27359648 | -.16197420 | .05013192 |
| preserving | .01 | -1.27181885 | -.44744979 | .17365349 |
| preserving | .02 | -1.28182054 | -.45558391 | .17218041 |

The two-step internal operator changes are .10513942/.14997399eV/Angstrom²
for T/preserving. Lowest-direction absolute overlaps are .99997753/.99999755.
These are raw consistency diagnostics, NOT total physical error bounds.
Removing the projected native tangent still leaves one raw negative direction
at both centers and both amplitudes. A single retained coordinate therefore
does not leave a raw-positive complement for these measured matrices.
This is matrix algebra at these centers, not a novel theorem or an index-two
transition-state certificate; their gradients are nonzero.

The original seven strong references are reapplied without rank loss at both
amplitudes, giving28center/amplitude/reference rows. The complete T optical
rank33reference leaves three open-cell coordinates with raw-positive release
blocks (T minima2.95621/2.95656; preserving3.09080/3.09217eV/Angstrom²).
All smaller registered references retain raw negative release directions.
The full optical Schur retained block still has two raw negative directions.
Do not turn this into a claim that every fixed reference fails, nor compare
rank2adaptive and rank33reference spaces as equal-dimensional methods.
Source LDA frequencies are not substituted for the measured PBE curvature.

## Actual next experiment, not a quadratic speed claim

The tested E158worker is deployed in a fresh immutable runtime from015719b.
Only one training modeling node Q=(0,0) per existing center is authorized here.
The first-amplitude two-negative subspace is frozen for this canary; its
34-dimensional complement is raw-positive at both amplitudes. These are
canary-specific axes, NOT final physical surface/forecast domains or phonons.

Identity and curvature-seeded safeguarded BFGS use separate evaluator
namespaces, equal .03eV/Angstrom released-gradient targets, .05Angstrom
initial step radii, .5Angstrom total internal-domain bounds and1.6Angstrom
minimum distances. Each arm reserves at most12new attempts, including
failures; at most24per center. Only exactly zero internal displacement reuses
the audited center E/g. No end-point or center SCF is repeated.

One node/32genuine MPI per study allocation, hfacnormal01,2hour caps,
OMP/MKL/OpenBLAS1, original ABACUS/PBE/100Ry/full10auDZP/Gamma2³/P0/E0.
Ordinary NEB remains.10eV/Angstrom; no CI, cutoff/input changes, automatic
resubmission or production-source overwrite. At most two committed study
allocations including runnable pending. Submission receipts retain one-shot
intent, accepted handle, held resource validation and release.

Actual jobs28859804(preserving) and28859983(T) were submitted once, held
resources checked, then their original handles released. Both first new SCFs
have also been independently reparsed from native output using the pinned
runtime:115.13s for T and95.46s for preserving, genuine32MPI, original six
hashes/ordered geometry/EFS and gradient replay passed. All378runtime Python
file hashes still match. The receipt includes original intent/accepted/held/
release/queue records, not a reconstructed successful submission narrative.

Native preflight is zero-DFT; a submitted/running job or a gradient-converged
candidate is NOT an accepted stable conditional face, TS, equal-branch speed
benchmark or independent prediction. Canary outcomes must be audited after
completion, with all296preparation SCFs and prior G2 costs reported separately.
Individual completed SCFs may be rejected backtracking trials; their last
gradient is not the last accepted optimizer iterate or a convergence verdict.
No physical+.005geometry/label has been generated or read.

## Evidence and replay

`native_reports.zip` losslessly preserves the exact two complete native audit
reports; SHA256 pins, job receipts and validation identify the original bytes.
It contains our own DFT-derived reports, no pseudopotential/orbital payload or
raw/transformed author reference arrays. `complete_summary.json` comes from
the unchanged `scripts.summarize_hfo2_G3_curvature` CLI. The28row
`reference_complements_complete.json` comes from the unchanged
`scripts.analyze_hfo2_G3_reference_complements` CLI, using the same private
registered author-source bank as E156; only scalar diagnostics/source hashes
are exported. Reproduction of that reference calculation requires authorized
access to that bank, not merely this scalar report.

Extract the two native JSONs into a fresh scratch directory; run the existing
CLIs with the original E150manifest, the exact native pins in validation.json,
and a fresh output path. Tests additionally replay complete-data arithmetic
and canary configurations; tolerances on derived floating-point arithmetic
are not physical accuracy criteria. Native bytes and source identities remain
exact. No new plotting, manuscript claim, CI or release is implied here.

Eight new actual-complete-data checks plus the unchanged related worker,
optimizer, evaluator, stencil and reference tests passed154/154locally and
154/154in the independent tracked snapshot (6538files, import roots checked).
The complete-summary CLI replay is byte-identical. The private-bank28row
reference replay has identical scientific/non-code fields including derived
floats; its full bytes differ only in actual code SHA fields from line endings.
Native/configuration artifact bytes remain exact. No new full-suite result is
claimed; E157's complete regression is historical, not rerun for this dataset.
