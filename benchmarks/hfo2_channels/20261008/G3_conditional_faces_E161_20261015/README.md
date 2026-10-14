# E161: first5x5training conditional-candidate sheets

This continues the original two-center/two-face G3 budget, following the
complete two-amplitude E159data and actual E160seeded-optimizer canaries.
It does not add a material, Hessian center, full band or final predictive
representative. The two sampled representatives remain exploratory under
the original uncertainty/tie rule; physical+0.005geometry and labels stay sealed.

## Prospective scope, fixed before new face labels

Each existing T/preserving center retains its E159two-negative coupled
internal span, in the declared joint length metric, not named parent phonons.
Freeze both axes at[-.2,-.1,0,.1,.2]Angstrom,25nodes/face including the
already native-audited E160Q0candidate. Reuse that exact point/E/g without
new DFT.24new nodes each reserve at most12attempts,288maximum per face.
Failed/interrupted attempts are charged, not refunded.25of81modeling nodes
are occupied;56further training and all8independent test points remain.
The296curvature SCFs,31canary SCFs and prior G2/center costs remain additional.

Same ABACUS/PBE/100Ry/full10auDZP/Gamma2³/P0/E0 and original six asset hashes,
ordered Hf4O8, same biaxial strain0/first two prescribed cell rows, all three
original open-row directions. Three translations are excluded only by the
declared internal embedding; no ASR or source-force change. Fixed control is
excluded exactly, not silently optimized. Raw40Dgradients including the
substrate reaction are retained for subsequent mode-strain analysis.

The successful canary settings continue unchanged: released-gradient norm
target.03eV/Angstrom (NOT ordinary NEB's.10criterion), safeguarded BFGS with
the original measured first-step Hessian seed, initial radius.05Angstrom,
total internal displacement norm.5Angstrom, minimum distance1.6Angstrom.
No eigenvalue clipping, added electronic accuracy, CI or input repair.

Point order is fixed Manhattan distance from Q0, then integer-index order.
At each new Q choose the nearest safe already gradient-converged candidate;
ties use fixed point ID. Carry only its orthogonal component and prescribe
the new Q exactly. One continuation start per point, not a search for the
global conditional minimum. Invalid starts/optimizer-budget outcomes remain
holes; they are never filled by the last completed SCF. Real native backend
failures stop the remaining face and preserve original artifacts. No retries,
automatic submission/rollover or silent parameter changes.

Each completed point is recorded in a new exclusive per-point JSON and owns
a separate evaluator namespace/reservation budget. A status reader must skip
a temporarily incomplete JSON rather than call it a scientific failure. Only successful points
can seed later ones. Cached source files/runtime remain read-only. Native
source/asset/original-center and cached conditional-origin EFS checks precede
any new SCF. The worker refuses an existing output directory/configuration.

## Resources and interpretation

Two study allocations maximum including runnable pending, each one node,
32allocated CPUs/genuine mpirun32, OMP/MKL/OpenBLAS1, hfacnormal01,12h cap.
The controller and serial SCF calls share that allocation; no third manager
job, nested Slurm submission or endpoint recomputation. Expected cost must
use actual SCF counts and Slurm accounting, not25nodes=25SCFs.

This batch produces a training conditional-candidate branch sheet, NOT a
certified PES, global branch, TS or independent prediction. Released-block
stability, competing starts, interpolation validation and the complete
T/Cmma strong controls remain required. No contours across failed nodes or
outside sampled/validated coverage, no full-path barrier replacement from
a local point and no new claim for Schur elimination or preconditioning.

Sources: scripts/hfo2_G3_conditional_faces.py and
cluster/hf_hfo2_G3_conditional_faces.slurm. Prepare freezes configs without
DFT; preflight reparses native source/cache without DFT; run requires exact
config SHA, reviewed allocation and a fresh sibling of the pinned source.
Deployment/dispatch/tests are recorded only after actual execution below.

## Actual tested deployment and launch

189local/189independent related checks pass,0failure/error/skip, including
16new face tests. The two independent prepare CLIs reproduce both frozen
configs byte-exact. Treee3811a09bd1835f72fa00549a6fd06a200792e6e has all6565
tracked files before tests, actual import roots checked. Code/test/wrapper
texts match modulo recorded line endings; original scientific bytes stay exact.
No new whole-suite result is claimed; E157's full regression is historical.

Fresh native runtime:380Python files and10data inputs verified exactly,
4491264byte bundle SHA7f9e0208d84183cdd5f3beccb5d979d8cba8cf72ed1ed32774ae784721210ee7,
Python-file manifestd07bb011c9e313afe23ca4e78207356ff98adfc43ca231bdd047e6d5d0f1be20.
Help import, native shell syntax and RUN_DFT0exit2guard pass. Both native
preflights audit existing center/conditional-zero EFS without new DFT.

One-shot intent/accepted-handle/held resource validation/release receipts
show actual jobs28862319(T,node13) and28862441(preserving,node321), both
RUNNING at initial observation. No duplicate submission/old-source overwrite.
Both first new Q=(-.1,0)SCFs independently reparsed:116.22s/146.12s,
original six asset hashes, SCF/DSIZE32/EFS/ordered geometry/full gradient,
exact fixed control and unchanged first two cell rows pass. All380runtime
Python hashes still match. This confirms NONZERO grid points, not repeatedQ0.
Read-only raw report and original dispatch receipts are in first_native_SCF_audit.json.

The sheet itself remains in progress. Native success at one new point does
not establish25conditional minima, release stability, valid contours or H1/H2.
No physical+.005geometry/label, final predictive representative/domain or
independent face-test outcome has been generated by this launch.
