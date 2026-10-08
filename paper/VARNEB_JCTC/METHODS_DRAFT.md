# Theory and methods working draft

Status: 2026-10-08; definitions and proposed tests, not an evidence-complete paper.
The research design fixes pure periodic HfO₂, zero temperature and zero electric
field. Numerical sections must ultimately name the effective electronic inputs,
mechanical ensemble, ordered atom correspondence and structural variants.

## 1. From NEB to atomic–cell paths

For an ordinary NEB, internal images represent a discrete path connecting two
fixed endpoints. The physical force perpendicular to the local tangent relaxes
the path, while a tangential spring force controls image spacing. A small NEB
residual does not imply that the physical tangential force at every image is
zero; only a properly refined stationary transition structure can support such
an assertion. Therefore the highest sampled image is a transition-state
candidate, not a transition-state certificate.

Variable-cell paths extend the active coordinates to atomic displacements and
cell deformation. We retain an explicit reference cell, ordered atoms, periodic
image choices and a length scale for strain coordinates. Cell rotations and
rigid translations are not reaction coordinates. This is an implementation of
established variable-cell/solid-state NEB ideas, not a new NEB theory. The
optimizer and mechanical active-space projection act on energy/force/stress
results through an ASE-compatible interface independently of the backend.

Under hydrostatic pressure the objective is H=E+PV. At fixed epitaxial vectors
and zero external normal/shear loading, the objective is E restricted to the
declared movable variables. A nonzero constrained stress is a substrate
reaction, not failed convergence. Every image and both endpoints of a channel
must obey the same substrate matrix. Released and clamped zero-strain cells
are distinct ensembles, even if their nominal strain labels coincide.

## 2. Channel selectivity requires a common initial state

For a channel α from the same initial structure I, write

`B_α(c) = E†_α(c) - E_I(c)`,

where c denotes a mechanical condition and E† is the relevant stationary
bottleneck energy, or explicitly a provisional sampled maximum. The identity

`B_α - B_β = E†_α - E†_β`

shows that the common initial-well energy cancels. A changed PO well depth
alone cannot explain a changed barrier *difference* between two routes from
that PO state. Absolute barriers, phase energy differences and channel
selectivity are therefore reported separately. Formation T→PO and decay
PO→M do not share an initial state and are not compared as competing branching
rates merely by their forward barriers. Reverse barriers obey
`B_forward - B_reverse = E_F - E_I`.

At a smooth, stationary constrained branch, the envelope derivative for a
continuously varied strain is

`dB_α/dε = dE†_α/dε - dE_I/dε`.

The derivative is evaluated through the actual strain-to-cell Jacobian and
stress convention; unequal cell volumes and biaxial factors are retained.
For same-initial-state selectivity its derivative is the difference between
the two bottleneck responses. This established identity motivates an
independently testable prediction, rather than being claimed as a new theorem.
It is not applied across a nonsmooth branch crossing or to a nonstationary
image without qualifications.

## 3. Modes, local coordinates and conditional branches

Parent-symmetry distortions classify variants and signs; local phonon
eigenvectors characterize curvature at a chosen structure. They serve
different purposes. A Γ eigenvector of a tetragonal 12-atom cell is not
automatically the X2− eigenvector of a cubic parent. The parent cell and the
integer embedding/folding map are required when finite-q labels are used.
For Γ force constants of the existing 12-atom model, a 1×1×1 displacement
cell suffices; no complete phonon dispersion is implied.

A frozen two-mode slice holds unselected coordinates fixed. A conditional
branch surface instead minimizes over permitted orthogonal coordinates at
fixed mode values. The global lower envelope additionally requires all
relevant branches. These three objects must not be conflated. A large
off-plane reconstruction residual, unstable transverse direction, or branch
switch is an observable failure of a particular reduction, not a reason to
force the path onto a visually attractive plane.

Let q be retained coordinates and r orthogonal movable coordinates in a
common declared metric. For a local quadratic energy with stable H_rr,
orthogonal relaxation gives

`δr* = -H_rr^(-1)(g_r + H_rq δq)`,

`g_eff = g_q - H_qr H_rr^(-1) g_r`,

`K_eff = H_qq - H_qr H_rr^(-1) H_rq`.

At an orthogonally stationary reference g_r=0. The curvature correction is
positive semidefinite, so stable orthogonal release softens local retained
curvature in this approximation. This does **not** prove that release always
reduces an activation barrier: the initial state can relax more strongly than
the bottleneck, endpoints may change, and global branch identities may change.
With the local expansion referenced at r=0, its relaxation energy is
`ΔE_relax = -g_r^T H_rr^(-1) g_r / 2`. Comparing the initial and bottleneck
relaxation contributions is a proposed mechanism diagnostic.

The implementation refuses to invert H_rr when its minimum eigenvalue does
not exceed a measured resolution floor. Acoustic translations are removed
explicitly; a physical unstable direction is retained or assigned a separate
branch, not hidden by a pseudoinverse. The Schur complement is standard
mathematics. The candidate method contribution is a tested branch-aware
selection/augmentation workflow with predictive benefit, if the material
experiments demonstrate that benefit.

## 4. Verification and falsification design

Fixed-Hamiltonian atomic and strain probes compare finite-energy work with
analytic force/stress work at two step sizes. This checks signs, units,
coordinate Jacobians and numerical resolution without altering cutoff or
orbital sets. A local positive or negative directional curvature is not a
full Hessian index and does not certify a variable-cell saddle.

Reduced-model predictions are frozen before independent DFT holdouts. Compare
endpoint-only, atomic-mode, mode–strain-release and branch-augmented models.
Report incorrect rankings, unresolved barrier differences and necessary
extra directions alongside successes. A DFT-backed contour requires sampled
coverage, branch continuity and holdout accuracy; interpolation is not added
DFT evidence. Efficiency is measured in actual energy/force/stress calls and
core-hours, with identical physical contracts and correctly identified
channels, rather than in optimizer steps alone.

The closest prior works and pending full-text audit are listed in
`../../docs/VARNEB_JCTC_HFO2_RESEARCH_PLAN.md`. No conclusion of a favorable
strain window, novel mode coupling, accelerated convergence or JCTC-level
scientific novelty is made here before those tests.
