# E179 — completed preserving batch, native audit and all-hole replay

The original28862441 job finished normally (COMPLETED0:0,21951s) on
node321/hfacnormal01 with32allocated CPUs. All25registered coordinates
were processed;17including the cached origin passed the released-gradient
target and8did not. This is batch completion, not a fully converged surface.
The complete native report is SHA982d2657; terminal face SHA cb97225d.

## Terminology retained from the manuscript

| Term | Used here | Not inferred |
|---|---|---|
| Gradient-converged candidate | Registered fixed Q and released norm at most.03eV/A | Stable minimum, TS or unique branch |
| Conditional candidate sheet | Continued candidates, with missing points retained | Globally minimized or smooth PES |
| Nonconvergence hole | Original algorithm did not supply a target-passing candidate | Failed SCF or nonexistent physical solution |
| Completed-call proxy core-hours | Sum of measured call seconds times32/3600 | Total campaign cost or speedup |

## Native evidence and costs

The unchanged E163 standalone auditor was run in complete mode against
the original E161 runtime, not against the changed checkout. All380original
Python source identities, configurations, handle, cached origin, six physical
assets, genuine DSIZE32, SCF convergence, E/F/stress, fixed substrate/Q,
reservations and all24nonzero-coordinate optimizer trajectories passed.
Exactly172new SCFs/attempts were reparsed; no backend failure artifact was
observed. The cached origin was not recomputed. Native completed-call proxy
cost193.859072782core-hours; scheduler allocation195.12core-hours. Original
curvature/canary/G2 preparation is additional, not free.

All eight holes were then replayed using the unchanged E168 guard observer,
without another raw-EFS parse or DFT call. `replay_holes.py` first checks the
exact complete proof, original380files and every used point/result byte
hash. It retains original safe-parent selection/start and calls the original
cache-only optimizer. The observer classifies the actual unchanged guard
decision and restores its hook. It does not provide a repaired guard.
Full output is SHA d7c89c51; deployed driver a3c277b5.

| Point | Q1,Q2 (A) | New SCFs | Stop | Radius rejections | Minimum observed trial norm (eV/A) | Last radius (A) |
|---|---|---:|---|---:|---:|---:|
| 13 | -.2,-.1 | 12 | attempt cap | 35 | .033791 | .499885 |
| 14 | -.2,+.1 | 12 | attempt cap | 40 | .034194 | .499944 |
| 19 | +.2,-.1 | 12 | attempt cap | 31 | .031416 | .499793 |
| 20 | +.2,+.1 | 12 | attempt cap | 31 | .031335 | .499711 |
| 21 | -.2,-.2 | 12 | attempt cap | 78 | .051496 | .500000 |
| 22 | -.2,+.2 | 12 | attempt cap | 65 | .050990 | .499997 |
| 23 | +.2,-.2 | 10 | no starting point reached target | 75 | .049984 | .500000 |
| 24 | +.2,+.2 | 12 | attempt cap | 81 | .049875 | .499999 |

The table rounds display values only. All94complete trials across these
holes remain in completed_holes.json; none passes.03eV/A or is promoted.
Every one of436recorded rejections was at the original.5A total-radius
short circuit. Rejected proposals were not subsequently distance-tested;
this is not evidence that they would be physically safe outside the domain.
No source geometry is moved inward or tolerance enlarged.

The observations establish domain contact, not boundary KKT, an unreachable
unconstrained minimum, a unique cause or a physical branch termination.
The earlier12call feasible-step arm atpoint13 also failed and is not
repeated here. Neither larger caps nor a radius change is justified as an
automatic response. Holes stay visible; no interpolation fills them.

## Claim-to-evidence and manuscript integration

| Claim | Evidence | Scope |
|---|---|---|
| Batch completed, but17/25gradient coverage | Complete native audit and exact terminal face | Not full conditional convergence |
| Eight nonconvergences are not SCF crashes | All172native records and original replay | No claim about untested inputs |
| Recorded guard rejections are radius-limited | All8unchanged cache-only observer replays | Not physical nonexistence or stability |

Section3.7 and the draft preface now state completion, holes, actual cost
and the boundary limitation. Canonical terms follow the existing terminology
ledger. All119previous table lines,9figure files,13external citations and
the TeX wrapper are preserved. No new theory or independent-success paragraph
is added. The existing LuaLaTeX build passes:32pages,0overfull/underfull or
missing-character warnings; all32contact pages and the updated page19
were inspected. PDF/build products stay on Cscratch, not Git; render_metadata
binds the actual source/PDF and records the inspected evidence pages.

## Validation and limits

First local64check run63pass/1fail: the new test incorrectly assumed every
hole used12calls. The actual point23used10; no source or optimizer was
modified to match the assumption. Correct the expected actual count to94
and retain first_tests.xml. Complete original native audit and all8diagnostic
replays passed before manuscript interpretation. Final related and
independent verification is recorded in validation.json; not a full suite.
Final scope:66local/66fresh-archive checks across3modules,16new,0failures/
errors/skips. Tested tree88510e04/archive4249b066. Manuscript, unchanged
wrapper, deployed driver and both native report files remain byte-identical
in that extraction. This is independent code/artifact validation, not the
protocol's still-unused eight independent material points.

The first independent archive gave65pass/1fail because Git converted the
unchanged TeX wrapper fromLFtoCRLF. Raw wrapper hashes differed, while its
LF-normalized source was exactly identical. Keep first_independent_tests.xml
and the original archive. Add an explicit LFattribute for that one wrapper,
without modifying its bytes or weakening the raw-hash test; validate a fresh
archive. Native-data/input fingerprints and scientific gates are unchanged.

Raw native logs remain at their immutable HFpaths with hash-bound receipts.
No private author vectors, physical+.005geometry/label, new Hessian center,
material, whole path, DFT call or scheduler write was introduced. Settings
remain100Ry/full10auDZP/.10ordinaryNEB; conditional.03target unchanged.
Training commitments T34/81,P35/81 and8independent slots each unchanged.
The atomic-release28867177remained PENDING dueCPUquota at06:54:13CST;
pending is a verified wait, not a backend failure or duplicate-submit trigger.

Next consume actual atomic-treatment evidence, settle a training-only valid
domain with released stability and strong controls, then freeze independent
predictions under the existing protocol. The completed candidate batch is
not by itself independent prediction evidence or JCTC readiness.
