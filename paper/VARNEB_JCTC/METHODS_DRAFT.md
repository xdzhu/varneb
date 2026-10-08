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

Periodic lifts are part of the input contract: independent image wrapping can
leave all DFT observables correct while introducing spurious jumps into the
unwrapped tangent/spring metric. Under an explicitly declared short-adjacent-
step convention, integer lattice shifts are registered before computation
and the optimizer rejects an inconsistent supplied lift. No runtime remapping
or structural repair is performed. Half-cell ambiguities require additional
path information. Geometry-contract correctness is an engineering prerequisite,
not evidence of accelerated convergence or a new scientific method.

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

An increased decay-minus-switching barrier difference is only improved
relative selectivity. It can occur even when both absolute barriers decrease.
The stronger proposed decoupling criterion separately requires a significant
switching-barrier reduction without a resolved reduction of the easiest
examined nonpolar-decay barrier. An absolute decay-barrier increase is a
stronger outcome. Report both M and T leakage channels and measured numerical
uncertainty; do not turn relative selectivity into a device-retention claim.

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

## 5. Explicit parent-site gauge for HfO2 variants

An ideal fluorite quarter-site scaffold is embedded in the unmodified
twelve-atom T lattice. This geometric scaffold is not a separately calculated
cubic phase. A declared operation maps fractional column coordinates as
`f' = R f + t`; its species-preserving parent-site permutation is recorded
before applying it to a distorted product. Product atoms are not subsequently
relabelled by nearest-site fitting. The operation must preserve the metric
when the cell is retained. Inversion-related product structures are geometric
switching candidates; a Berry-phase branch check is a separate validation.

The four fcc translations in this cell give orthogonal character projectors
onto Gamma and the three X-wavevector sectors, in the scaffold reciprocal
axes. Their sum reconstructs the displacement after a stated translation
gauge is removed. These q sectors contain multiple irreducible representations:
an X-sector amplitude must not be called X2- without further mode identification.
Likewise, the historical T distortion is only one vector in this space;
its scalar projection is not a general channel label.

For the local T-cell dynamical matrix, generate central atomic displacements
at 0.01 and 0.02 Angstrom with a 1x1x1 supercell and identity primitive matrix.
Phonopy symmetry reduces each family to eight independent force calculations.
The Python workflow consumes Angstrom coordinates and eV/Angstrom forces,
so its force constants have units eV/Angstrom^2; it does not use the native
ABACUS Bohr-coordinate conversion factor. Raw force constants, acoustic-sum
drift and permutation asymmetry are retained. ASR/permutation symmetrization
is reported separately, and two-step disagreement is part of the uncertainty
assessment. Acoustic modes are identified by the rigid-translation subspace
overlap, not assumed to be the first three sorted frequencies. This is an
analytic fixed-cell Gamma matrix without LO-TO/NAC corrections, not a
finite-q dispersion or a joint atom-cell saddle certification.

## 6. Electronic validation of inversion-related endpoints

Geometric parent operations alone do not establish an electronic switching
branch. For the common PO starting well and two ordered inversion products,
we evaluate the native LCAO Berry polarization along the third lattice vector.
The labels PO+ and PO- refer to the registered starting/product wells, not a
choice of the smallest absolute polarization value. The original SCF geometry,
100-Ry cutoff, full10-au DZP orbitals, pseudopotentials and2x2x2 mesh are retained.
Only charge/band-gap output is added, and the resulting energy, forces and
stress must reproduce the original endpoint before further analysis.

On that fixed charge, NSCF meshes2x2x2,2x2x4 and2x2x8 test the longitudinal
quadrature. Their total energies are excluded from activation barriers.
The actual ABACUS3.10.0/f7cb1d3 native output and48 occupied bands are audited:
the96-electron count follows from the actual Hf/O UPF valences12/6. For this
version LCAO NSCF eigenvalues require explicit band output; a completed Berry
log does not itself establish a sampled insulating gap.

For row cell vectors a_i the modern polarization quantum lattice is e*a_i/V.
The spin-paired native implementation with even ionic valences reports an
approximately doubled R3 modulus and uses older SI conversion constants.
Both the native value/modulus and the explicit modern-SI conversion are stored;
the polarization is never divided by two or assigned a physical path branch
by minimizing its magnitude. Endpoint inversion is checked modulo the native
period after unit conversion, while self-inverse classes are additionally
checked using the physical eR/V quantum. This endpoint test does not determine
the absolute spontaneous polarization, the full Cartesian polarization or
the continuously unwrapped change along either optimized switching path.
The transverse mesh is not converged by the longitudinal test.

The source-version output protocol follows the
[official native LCAO example](https://github.com/deepmodeling/abacus-develop/tree/f7cb1d3/examples/berryphase/lcao_PbTiO3)
and [Berry implementation](https://github.com/deepmodeling/abacus-develop/blob/f7cb1d3/source/module_io/berryphase.cpp).
These established Berry/unit conventions are validation requirements, not
claimed methodological innovations of VARNEB.

## 7. One mechanical ensemble for endpoints and bands

For a prescribed substrate, two non-collinear ASE row-cell vectors are held
fixed at every endpoint and image. With H=H0 F^T and their unit normal n,
the exact allowed deformation is deltaF=v n^T. Three components of v permit
normal relaxation and two out-of-plane tilts; restricting v parallel to n
instead defines a different, normal-only ensemble. The present finite strain
study selects the tilt-released ensemble. A general rotated substrate cannot
be represented by independently zeroing Cartesian stress components.

The same orthonormal deformation subspace is applied to the endpoint BFGS
filter and the NEB controller. Original seed cells are checked before general
subspace projection, rather than silently converted to a common substrate.
In-plane reaction stresses are recorded but excluded from the open-gradient
stationarity criterion. Geometry and EMT energy/stress finite differences
verify the implementation; they do not establish material barrier accuracy.

For substrate area A=|H_0 cross H_1|, moving only the third cell vector by
w performs work A w dot ((sigma+P I)n). The open-stress residual is therefore
the norm of (sigma+P I)n in the tilt-released ensemble, or its normal
projection in the normal-only ensemble. Unlike a Cartesian component test,
this criterion is invariant under a rigid change of coordinates. Endpoint
BFGS applies separate atomic-force and open-traction gates, 0.03 eV/Angstrom
and 2 kbar respectively; in-plane reaction stresses remain part of the data
but not the convergence gate. The ordinary NEB threshold remains 0.10
eV/Angstrom and is not a saddle-certification criterion.

The finite strain seeds use one proper cyclic rotation, with new Cartesian
x/y/z=old z/x/y and the same cyclic cell-row relabeling. Ordered atomic
identity is unchanged. Every phase shares the T-referenced substrate at
epsilon=0 or +0.01; +0.005 is reserved for an independent holdout. The zero
T-referenced substrate already changes the PO free-cell in-plane lengths by
approximately +3.12% and -3.61%. All affine starters consequently require
fresh clamped relaxation and independent phase/variant checks; free-cell
evaluations are not reused as their energies or forces. These prepared
geometries are not yet strain-dependent material results.

For an explicit cell derivative dotH, the fixed-fractional-coordinate work is
V(sigma+P I):(dotH^T H^{-T}). This is a branch envelope derivative only when
all open variables are stationary. On a continuous endpoint/saddle branch,
subtracting their work derivatives predicts a barrier strain response.
An ordinary NEB peak need not be a stationary saddle, and fixed-cell Gamma
phonons do not certify the joint open atom/cell Hessian index. Released P=0
and clamped epsilon=0 are compared as different boundary conditions, never
as adjacent points of a single strain derivative. These work and envelope
relations are established tools, not claims of new mathematical theory.

Curvature probes must obey the same permitted deformation subspace as the
endpoint and band. A tilt-released fixed plane admits three rank-one
deformation directions v*n^T, generally nonsymmetric. Replacing them with
the six symmetric free-cell strains would change the substrate and the
mechanical experiment. `ActiveJointCurvatureCoordinates.for_clamped_plane`
uses the exact boundary basis, retaining all atoms and removing three
translations for curvature analysis. For twelve atoms this chart has39
coordinates and36 after translation removal. Zero cell directions recover
the atomic fixed-cell chart. The historical symmetric-strain implementation
and its published-source hash are preserved separately. Geometry and EMT
work checks validate this implementation, not an HfO2 bottleneck Hessian,
its index, or independent material-barrier prediction.

For a declared orthonormal joint-coordinate direction matrix B, paired full
physical gradients measure A_j=[g(y0+h B_j)-g(y0-h B_j)]/(2h), approximating
H B. The projected matrix B^T A retains its raw reciprocity defect before
symmetrization. The transverse response A-B(B^T A) exposes coupling outside
the sampled subspace, but does not measure the complement's self-curvature.
It cannot establish stability of unsampled modes or justify their elimination.
Both directions and steps use the same chart metric and mechanical boundary;
NEB-projected and spring forces are not Hessian inputs. Signed energy checks
at two amplitudes provide an independent gradient-consistency diagnostic.
Step-size spread is not a rigorous total DFT uncertainty. The staging and
assembly interfaces call no calculator and change no physical input; separate
geometry/input/raw-SCF provenance remains required for material use.

## 8. Continuous reference projections and numerical observation audit

For a supplied, already continuous periodic lift, integer representatives are
chosen only at image zero relative to the ordered reference. The same
per-atom integer shift n is applied at every image, giving
u_i=(q_i+n-q_ref)H_ref in a fixed reference metric. Independent nearest-image
folding to the reference is not performed: it can hide a winding or introduce
artificial mode-coordinate discontinuities. Ambiguous initial half-cell
representatives, broken adjacent lifts and changed atomic order are rejected.
No calculator is invoked and no image geometry is modified by this analysis.

Cartesian geometric patterns use unit metric weights; Gamma eigenvectors
use their original phonon mass metric. Rigid translation is removed with
the corresponding weighted mean. Complete degenerate-subspace squared
amplitudes are invariant to internal basis rotations; individual mode
coordinates are not. Cell deformation F=solve(H_ref,H_i)^T and Green strain
(F^T F-I)/2 are recorded separately. This is a declared coordinate convention,
not a new decomposition theorem or a partition of energy among modes.

Every retained observation matches a complete, identical ordered SCF and
binds input, structure and raw-log checksums. Complete numeric E/F/stress
are then attached to nonexecuting SinglePointCalculator objects to reproduce
the ordinary extended-space NEB force. Atomic, cell, physical perpendicular
and spring contributions are reported separately. A residual gate, an
optimizer step limit, a scheduler exit status and full-variable saddle
stationarity are distinct checks; none is substituted for another.

## 9. Common-initial-state channel comparisons and bounded responses

For an explicit set of candidate channels, directional barriers are evaluated
from one ordered periodic initial state, raw initial energy, and mechanical
ensemble. The reverse view of an existing T--PO band supplies PO--T without
an additional path calculation. This change of thermodynamic view does not
relift atoms, change their identities or replace the source optimization
metric. Each image retains its raw energy and volume; at prescribed pressure
the profile is H_i=E_i+P V_i. Both directional maxima include the endpoints,
and B_forward-B_reverse=H_final-H_initial is checked.

Let B_s and B_d denote the minimum barriers within the explicitly declared
switching and decay candidates. If B_j lies in a measured interval[L_j,U_j],
its candidate-set minimum lies in[min_j L_j,min_j U_j]. The lowest candidate
may change within these bounds. Differences use interval subtraction rather
than an assumption of independent Gaussian errors. Missing leakage channels
block the minimum; missing numerical/sampling evidence blocks a bounded
selectivity verdict. The ordinary force threshold is not converted into a
barrier uncertainty. Numeric error and sampling audits retain separate hashes.

For S=B_d-B_s, an increased S establishes relative selectivity only; B_d
can still fall. A stronger, declared-path noninferiority comparison requires
a resolved reduction in B_s and a lower bound onDelta B_d no smaller than
-delta, with delta specified before examining responses. Error overlap
alone does not establish noninferiority. Released and clamped ensembles
are not merged into a single strain derivative. These definitions and
interval rules are conventional analysis, not new physical theory or proof
that no unexamined decay channel exists. The implementation calls no DFT
calculator and does not infer a device lifetime or switching rate.

## 10. Prospective validation scope

Before strained material labels are available, we reserve epsilon=+0.005 as
an in-range unseen condition between 0 and +0.01 within the same clamped-plane
ensemble. The candidate selected for each role uses training data only.
Predicted barriers, uncertainty bounds, branch labels, basis, code and data
hashes are fixed before complete holdout path labels are inspected. Independent
local-surface points remain outside the fitting set. A direct interpolation of
training barriers, endpoint-only responses, frozen modes, atomic-only release,
joint atomic/cell release and branch-aware reduction are compared explicitly.
Abstention on an unstable branch is a coverage result, not an accurate barrier
prediction. Two held-out edges do not certify the full network minimum.

The strict primary switching/decay test preregisters zero allowable decay-barrier
loss, separately from relative selectivity. This is an analysis margin, not an
altered calculator or NEB threshold. Unresolved intervals do not establish
noninferiority. Prediction errors and model improvements are compared against
measured reference resolution; an equally accurate interpolation baseline
precludes a claim of improved predictive precision. Corrections after inspecting
holdout labels are recorded as revised training, not successful original forecasts.

Constrained-mode landscapes and NEB-validated switching predictions are already
present in [Zhou et al. (2022)](https://doi.org/10.1126/sciadv.add5953), including
the described random-start subspace search in its supplement.
[Qi et al. (2025)](https://doi.org/10.1103/PhysRevB.111.134106) also discusses
oxygen-crossing pathways beyond its two main categories. We therefore restrict
claims to the declared candidate set, mechanical boundary and observed prediction
advantage. These sources do not establish the novelty or accuracy of our proposed
joint-cell reduction; those remain questions for the independent material tests.

[Behara and Van der Ven (2022)](https://doi.org/10.1103/PhysRevMaterials.6.054403)
also reports strain-coordinate paths, a two-shuffle landscape, and a switching
intermediate changed by fixing the O cell. Our two registered pattern-sign
candidates are not an enumeration of all its translational variants. Complete
cell fixing and the present partial substrate clamp define different ensembles.

## 11. Auditing a strong published parent reference

To avoid attributing a poor fixed-reference choice to a general failure of
fixed-mode models, the prospective comparison includes the Cmma parent of
[Qi and Rabe (2025)](https://doi.org/10.1103/9759-kp38). Its published input,
structure and displacement files are pinned by author commit and byte hashes.
The author's cell and coordinate frames are connected by an explicit proper
rotation, integer basis change, translation and species-preserving assignment;
an independent primitive-to-conventional reconstruction checks the rounded SI
coordinates. Assignment ambiguity, cell mismatch and unreported strain fitting
are rejected. An inconsistency in the printed structure table is retained,
not silently corrected. These operations do not alter a production image.

QE matdyn flvec stores normalized Cartesian displacements rather than
orthonormal mass-weighted dynamical-matrix eigenvectors. With the author masses
and original row order, we restore the mass metric and retain complex phases,
printed normalization and Gram defects. Rigid translations are identified by
subspace overlap; centering character is measured separately from any irrep
label. No acoustic sum rule or implicit orthogonalization is imposed, and
nonzero-q modes are not reclassified as Gamma modes. Only the first complete
Gamma block is imported.

This audit establishes consistency among representations of the published
reference, not its mapping onto a production path or its predictive accuracy.
Published LDA modes are candidate geometric directions, not PBE curvatures or
energies. Preparation of a real orthonormal comparison subspace, ordered-path
gauge mapping and same-data material benchmarking remain explicit gates.
The [source and numerical audit](../../benchmarks/hfo2_channels/20261008/cmma_reference/README.md)
records these limits. No material setting or running calculation is changed.
