# E260: author-approved preserving branch trial, frozen before new labels

2026-10-17. The author explicitly approved the October17 blueprint, a new
preserving third-coordinate protocol, at most10 training candidates and
subsequent independent prediction tests, keeping physical inputs unchanged.
This is the new training experiment, not an overwrite of E254/E185/E188.

Five original coupled Q pairs in Angstrom, in order:
(0,0),(-0.1,0),(0,-0.1),(0,+0.1),(+0.1,0).
At each use two separate original-origin starts, +0.10 and -0.10 Angstrom
along the exact sign-fixed E255 omitted direction. No cross-point or
cross-sign continuation; every start is fixed before new labels. Q stays
fixed, zeta and the remaining permitted34 directions can move. There is
no imposed sign constraint, so both starts can reach the same candidate.
Record each signed final amplitude, full released gradient, directional
gradient, ordered geometry, energy, native records, failure and cost.

Use the unchanged physical length chart, fixed substrate,0 strain/P0/E0,
0.5A total radius and1.6A minimum distance. ABACUS/PBE/100Ry/full10auDZP,
original Gamma2x2x2, binary and six physical-asset hashes remain exact.
Use the original380-file runtime read-only; a new separately pinned
calculator-independent branch optimizer is an external numerical helper.

Numerical method: safeguarded BFGS, initial trust0.05A, maximum200 updates.
The original full released-gradient target0.03eV/A alone would admit a
shallow-direction start without resolving it. Require, additionally,
abs(dE/dzeta)<=1e-4eV/A. This is an explicit branch-optimization stop test,
not a changed DFT/SCF tolerance or ordinaryNEB threshold. No CI.

The positive numerical preconditioner uses the first-step measured
post-release R33 restriction and an uncoupled zeta stiffness obtained
from the four already seen E257 gradient probes: fit the odd directional
gradient a*zeta+b*zeta^3, then k_zeta=a+3*b*(0.1)^2. The same preconditioner
serves all10 starts. Its cross blocks are zero by declared numerical
design; it is not a newly measured physical Hessian, eigenvalue clipping,
or a certificate of stability away from the original point. If either
restriction is not positive, refuse before DFT. Original raw negative
curvatures stay unchanged. Armijo energy tests use actual native E/g.
The seed and diagnostic are frozen as numerical arrays. Their reconstruction
uses the existing2e-12 absolute arithmetic-replay tolerance, with all source,
target, physical and budget identities exact; optimization uses the frozen
seed, not a platform-regenerated array. This is not an energy-error tolerance.
The same arithmetic check applies to the frozen internal direction and
starting-coordinate transformations; Q/sign/amplitude/point identity remain
exact. The run uses frozen internal starts, and preflight writes the same
internal-to-full embedding as the actual evaluator. Original literal full
starts are retained for reference, not regenerated as replacements.

Before any new SCF, r1 preflight refused exact JSON equality for these
matrix-product arrays. Read-only native diagnosis found1382 last-bit
differences, maximum1.0658141036401503e-14, with no source/physical/policy
difference and no inert decks or DFT. Keep r1 and its first config/manifest.
The r2 namespace implements the explicit arithmetic check above, without
changing the ten starts, optimization stop tests or physical settings.

Each candidate has12 fresh attempted SCFs maximum;10 candidates reserve
120 maximum, no retries/rollover/additional starts. Nonconverged candidates
remain nonconverged; backend failure stops its shard and is audited before
any recovery. All attempts count, not just successful or accepted steps.
Training commitments go from T34/P39 to T34/P49, each cap81. Derivative,
training, independent and prior-preparation costs stay separate.

Two disjoint shards:0 IDs[0,2,4,6,8],1 IDs[1,3,5,7,9]. Each allocation
hf/hfacnormal01/exactly1node/32genuineMPI/mpirun/OMP,MKL,OpenBLAS1/4h;
at most two study allocations including held/pending. No endpoint or
whole-chain recomputation. Submit held, verify exact resources, record
accepted handles, then release those IDs once. Never resubmit blindly.

This trial asks whether two starts produce distinct reproducible gradient
candidates and how their local zero-strain energy/geometry differ. Neither
sign is an established branch, a global minimum, a TS or a3D surface.
Energy differences need their actual resolution; near-equal labels are not
declared two distinct phases. Small static lowering alone cannot establish
small mechanical response. Any further mixed-derivative or stationarity
probes require a concrete finite registration within the original budgets.

The remaining preserving8 independent targets stay unread/uncomputed
during this training trial. A revised forecast, domain, controls and test
definition must freeze before those labels; keep original refusal visible.
T8 are seen and never reused as new independent evidence. Physical+0.005
geometry/labels and private author arrays are not accessed. Success of this
trial does not by itself authorize reading them before final freeze.
