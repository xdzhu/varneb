# E166 — prepared matched frozen-R controls, eight statics and no new band

Previous goal turn E165 was PROGRESS (7d0c8895): original quadratic shape
disagreement was measured on already audited training candidates. The next
discriminating calculation compares the same Q with R=0 and with the saved
released candidate. Four axial points per original face are fixed before
their new frozen SCF labels are read: (-0.1,0), (0,-0.1), (0,0.1), (0.1,0).

The unique E166 training reservation commits four additional structure/node
slots per face:25+4=29 of81, leaving52 training and all8 independent slots.
It is not backdated into the immutable E161 launch policy. Two original
Q0/R0 SCFs and the two released origins are reused without recalculation.
The matched nonzero-Q candidates come from the hash-pinned E163 native audit.
These controls are seen-training diagnostics, not independent forecasts.

## Fixed scientific and execution contract

R=0 is relative to each original nonstationary sampled peak. Q is a coupled
atomic/open-cell coordinate; frozen R does not mean frozen full cell. Both
substrate rows and the separate physical strain control remain exact. The
eight full40D coordinate arrays, source/decision/reference hashes and budget
are fixed in config.json. Coordinate formula comparisons use algebraic
roundoff tolerance; physical INPUT/KPT/potential/orbital hashes remain exact.

ABACUS/PBE/100 Ry/full10-auDZP/Gamma2x2x2/P0/E0 and the validated mpirun-np32
transport are unchanged. Only static SCF calls occur: no optimizer, line
search, reference-axis change, parameter repair, new Hessian center or NEB.
Each point gets one conservative attempt and a fresh namespace; a real
failure stops the remaining calls. Completed prior points and failure files
are preserved; there is no automatic retry/resubmit. Signed frozen-minus-
released energy differences are never clipped to a favorable sign.

The compound job visits preserving then T, four calls each. Request: one
node/32 CPUs/one hour on hfacnormal01, with32 genuine MPI ranks and OMP/MKL/
OpenBLAS1. New work must wait until an original allocation normally completes.
Admission queries both original accounting handles and all visible G3
allocations, including pending jobs. Missing handles or non-normal parent
states require review; they are not inferred to be numerical failure.
At most two study allocations may be committed. The worker rechecks its own
visible job, identity and resource gate before any calculator call.

## Preparation, inert preview and run

From the repository root, prepare once to a fresh target:

```text
python benchmarks/hfo2_channels/20261008/G3_frozen_R_controls_E166_20261015/worker.py prepare --config a-new-config.json
python -m pytest -q tests/test_hfo2_G3_frozen_R_controls.py
```

On HF, preflight the pinned config and request a fresh external --decks
directory. This writes inert INPUT/KPT/STRU decks and symlinks the existing
four potential/orbital assets, validates original bytes and geometry roundtrip,
and reports current admission readiness without submitting or calculating.
Preview decks cannot be written inside the live E161 study, source, assets
or future calculation directory. Do not run ABACUS from a preview manually.

The run.slurm wrapper remains inert without RUN_DFT=1 and the exact approved
allocation/source/config identities. It never submits a child job. The
driver performs exactly eight single calls, not eight optimizations, and
reparses every complete native SCF's E/F/stress, gradient, reservation and
contract before recording a point. A CLI regression enforces one invocation
of the run function; point-count and per-point-attempt checks prevent extra
calls. No completed-Q0 calculation is repeated.

## Interpretation

The controls distinguish a frozen finite-Q slice from the declared relaxation
protocol under identical inputs. They do not by themselves establish a global
conditional minimum, smooth branch, physical curvature error bound, TS or
whole-path barrier. Nonzero retained/frozen residuals are reported, not treated
as SCF failure or corrected by changing the Hamiltonian. Original independent
points and physical +0.005 strain remain sealed.

See [the discriminating hypothesis](../G3_quadratic_agreement_E165_20261015/README.md),
[original candidate-sheet protocol](../G3_conditional_faces_E161_20261015/README.md)
and [native candidate sources](../G3_face_prefix_audit_E163_20261015/README.md).

## Actual preparation and verification, not a submitted calculation

Config SHA:5ffb1fe0d707aec5d0a6546a44e0ad0002cd0f0280f1ba17ca2d955ac2d4c930.
An independent prepare CLI reproduces its bytes exactly. The unique training
reservation is now prepared; four slots per face are reserved, not completed
DFT labels. Later sampling must subtract this reservation from the historical
E161 remaining56 rather than rewriting that old launch policy.

The first tested tree56f11176cb1725753eb35665cb359c7360693059 was deployed to
the fresh HF root20261015-varneb-G3-frozen-R-E166-r1, separately from the live
E161 study. All383 Python files and14 data inputs matched. Native shell
syntax and the unset/zero RUN_DFT gate passed; the latter exited2 before
entering a calculator. All eight inert native decks passed the real ASE
roundtrip and the original six asset hashes; four asset types are symlinks
to the existing library, not copied into this repository.

The native preflight report is preserved verbatim (SHA3ec6390084497dc3f3a23a34b8e02baaef72b80ffc773cbb2437ea2128329b12).
Its resource observation contained the two original RUNNING allocations and
returned admission_ready=false. No frozen-control job was submitted or run,
no frozen_R calculation directory was created and no new SCF occurred.

A separate exact-coordinate cache check rehashed28 registered original trial
result JSONs at these matched Q points. None matches a proposed frozen full
coordinate vector. This is not a global/nearest-geometry cache search. The
two original zero statics remain reused. Its native report is also retained.

Final related checks pass333 locally and333 in the independent final archive,
including26 new checks and the real native preview/config regressions. Toy
eight-call tests are explicitly not DFT. The final tree is
fc2aa5b7c356480034d3d19757ecb2792c80de5e; deployed worker text is unchanged
from the first tested tree. Exact runtime/source/config/test/preview evidence
is in validation.json. This does not establish completed frozen DFT controls,
released stability, a whole PES, a TS, full-cost acceleration or prediction gain.
