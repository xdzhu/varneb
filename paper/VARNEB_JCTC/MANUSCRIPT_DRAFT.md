# Competing switching and phase-escape channels in hafnia: a mode-strain pathway analysis

Working manuscript, 2026-10-15. This is a connected main-text draft, not a
submission-ready article. The main Results organize eight ordinary-converged clamped training bands
by scientific question; dated free-cell and detailed representation checks
are retained in the appendices. All ten
registered endpoint representations across the two matched-substrate training
conditions have passed their physical screens. The sampled G2 network is complete;
complete local joint-curvature measurements, their atomic/open-cell
restrictions and a finite-coordinate quadratic comparison are available.
The first two conditional batches have finished: 25 T and 17 preserving
coordinates passed the released-gradient target, out of 25 per center;
eight preserving holes remain. An uncertainty-resolved
response interpretation, stable conditional branches and independent
material predictions are still missing. An abstract and a
conclusion asserting those unmeasured results are deliberately withheld.
The [methods working record](METHODS_DRAFT.md),
[dated pilot record](PILOT_RESULTS_2026-10-08.md) and
[prospective protocol](../../docs/HFO2_PREDICTION_PROTOCOL_V2_2026-10-09.md)
remain the detailed sources; this draft does not replace their provenance.

## 1. Introduction

A polar crystal must do more than possess a metastable energy minimum to be
usefully switchable. A structural route to the opposite polarization can
compete with routes that leave the polar phase altogether. Lowering a
polarization-reversal barrier is therefore not, by itself, evidence that the
polar phase becomes more resistant to structural escape. Hafnia provides a
particularly useful setting for separating these questions because several
polar and nonpolar structures can be connected through coupled atomic
distortions and lattice deformation. Our question is whether a common
mechanical intervention can facilitate homogeneous switching without
facilitating the easiest examined route out of the polar well.

The connection between mechanical constraints and hafnia pathways is already
well established. Liu and Hanrahan studied orientation, epitaxial strain and
variable-cell transformation barriers, while Delodovici, Barone and Picozzi
identified coupled lattice distortions and their strain response.
([Liu and Hanrahan, 2019](https://doi.org/10.1103/PhysRevMaterials.3.054404);
[Delodovici et al., 2021](https://doi.org/10.1103/PhysRevMaterials.5.064405))
Behara and Van der Ven connected strain coordinates, structural variants and
switching pathways, including a change of intermediate when the entire
orthorhombic cell is fixed. Lattice-mode analysis by Qi, Singh and Rabe
further demonstrates that reversing polarization does not uniquely specify
the accompanying nonpolar distortions or the switching route.
([Behara and Van der Ven, 2022](https://doi.org/10.1103/PhysRevMaterials.6.054403);
[Qi et al., 2025](https://doi.org/10.1103/PhysRevB.111.134106))
Neither mode enumeration nor a strain-induced pathway change is consequently
a new result to be established by repeating these calculations.
Ma and Liu also examine competing M-to-PO mappings involving oxygen and
mixed hafnium-oxygen shuffles. Their supplemental Figure S2 makes clear that
a phase-pair name or one conserved shuffle does not uniquely identify the
path. Our clamped M candidate is therefore not assigned to a literature
mapping by barrier-height agreement; candidate coverage, ordered mapping
and cell-frame transport must remain explicit.
([Ma and Liu, 2023, publicly archived v3 supplement](https://arxiv.org/pdf/2204.09374v3#page=24))
Constrained-mode landscapes have also predicted switching mechanisms that
were subsequently checked by NEB; independent DFT validation alone is not
our novelty criterion.
([Zhou et al., 2022](https://doi.org/10.1126/sciadv.add5953))

The interpretation of a reduced landscape also depends on its reference and
mechanical boundary. Qi and Rabe organize competing structures using unstable
bands of a Cmma reference, providing a strong alternative to a selected
tetragonal-mode representation.
([Qi and Rabe, 2025](https://doi.org/10.1103/9759-kp38))
More recently, phonon-pair condensation has been used to explain low-energy
domain-wall motifs in hafnia. That interface mechanism cannot be equated
with a homogeneous transformation barrier.
([Lee et al., 2026](https://doi.org/10.1103/63lf-k7zs))
Fan, Zhu and Liu show that Pbcn energetics depend on the functional and
boundary, including strain-dependent switching under a fully prescribed
PO-referenced lattice. The present partially clamped T-referenced ensemble
is different; matching its electronic settings does not erase that distinction.
Their phase-resolved energy-gap approximation also agrees with NEB for the
switching routes examined there. Such a proxy and an initial/final-well response
null are distinct controls; neither is assumed valid for an unverified escape
bottleneck.
([Fan et al., 2025](https://doi.org/10.1038/s41524-025-01647-w))

Strain-dependent mode-amplitude models have also been compared with DFT under
uniaxial and biaxial loading and connected to hafnia's formation landscape.
([Lee, Lee and Yu, 2026](https://doi.org/10.1038/s41535-025-00841-9))
Predicting a strain-induced distortion is therefore not, on its own, the
methodological distinction sought here. The test concerns the competing
barrier responses after subtraction of the same initial well, with their
actual mechanical freedoms and strong controls kept explicit.
Their supplementary Part 2 also contrasts a coupled three-mode interpolation
with restricted relaxation. Relaxation-induced lowering of such a profile
is therefore not presented here as a methodological innovation.
Figure S3 uses a coupled interpolation parameter and cell energies, whereas
Table S1 labels its coefficients per formula unit. Amplitude normalization,
energy units and each starting-state reference must consequently remain
explicit; profile lowering alone does not certify our partially clamped
bottleneck or predict its barrier response.
([Lee et al., Supplementary Part 2 and Figure S3](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41535-025-00841-9/MediaObjects/41535_2025_841_MOESM1_ESM.pdf#page=3))

We therefore formulate a narrower test than discovery of a new coupled mode
or a low-energy intermediate. First, we compare switching and phase escape
from the same polar initial state, separating absolute barriers from their
difference. Second, we ask whether measured mode-strain responses and
explicit transverse-stability checks can predict their changes at an
unseen mechanical condition. The proposed reduction is tested against
strong fixed-reference, endpoint-response and direct barrier-interpolation
controls. Its success requires a resolved advantage in independent material
predictions; visual compactness of a mode plot is not sufficient.

## 2. Theory and computational approach

### 2.1. Joint pathways and a common energy reference

Ordinary NEB relaxes the physical force perpendicular to the path while
springs control image spacing. Atomic-cell extensions incorporate lattice
degrees of freedom with an explicit relative metric.
([Sheppard et al., 2012](https://doi.org/10.1063/1.3684549);
[Qian et al., 2013](https://doi.org/10.1016/j.cpc.2013.04.004))
VARNEB applies these established path ideas through energy, force and stress
results obtained via ASE-compatible calculators. The backend, optimizer and
mechanical active space are independent choices. Ordered atoms, reference
cells and continuous periodic lifts are preserved during optimization and
analysis; an image-wise remapping is not a physical pathway improvement.

For candidate channel $\alpha$ at condition $c$, define its bottleneck barrier
from the same initial well $I$ as

$$
B_\alpha(c)=E^\dagger_\alpha(c)-E_I(c).
$$

Here a stationary bottleneck and a discrete sampled maximum are distinguished;
the latter is denoted by a hat in the pilot interpretation. With a common $I$,
the initial-well energy cancels from $B_\alpha-B_\beta$. A change in relative
channel selectivity thus cannot be explained by shifting that common well
alone. In the declared candidate set, we use

$$
B_s=\min_{\alpha\in\mathcal S}B_\alpha,\qquad
B_d=\min_{\alpha\in\mathcal D}B_\alpha,\qquad S=B_d-B_s.
$$

An increased $S$ is relative selectivity, not improved absolute resistance to
escape. The stronger prospective test requires a resolved reduction of $B_s$
and a lower uncertainty bound on $\Delta B_d$ no smaller than zero. An unresolved
decay response is inconclusive, not evidence of noninferiority. This primary
allowable-loss margin is fixed before boundary-response labels. Both T and M escape candidates
are retained. These are candidate-set statements, not an exhaustive network
or a calculation of device lifetime. At finite prescribed pressure the
objective would be $E+PV$; all hafnia results here instead use $P=0$ and $E$.

### 2.2. Mechanical ensembles and branch-aware reduction

All images and endpoints of one path share the same permitted cell space.
The planned epitaxial experiment fixes two T-referenced in-plane vectors
and releases the third vector, including its two tilt components. Atomic
and cell variables are relaxed in that same space. A clamped reaction stress
is not an open-coordinate convergence failure. Zero applied stress and
zero nominal substrate strain are distinct ensembles, and are not adjacent
points of a single strain derivative.

Adding a constraint orders ideal minimax barriers only when the endpoints
and thermodynamic potential are held common. For allowed continuous paths
$\mathcal G_U(I,F)$ and $\mathcal G_C(I,F)\subseteq\mathcal G_U(I,F)$,

$$
B_U^*(I,F)=\inf_{\gamma\in\mathcal G_U(I,F)}
 \max_{s\in[0,1]}E(\gamma(s))-E(I)\le B_C^*(I,F).
$$

The star denotes the ideal minimax value, not a finite-image NEB estimate.
Neither a local optimizer nor discrete sampling establishes that ordering
for two computed candidate bands. More importantly, newly relaxed endpoints
under a different substrate do not satisfy the common-endpoint premise.
Our free-cell and T-substrate-clamped results therefore are not interpreted
as a nested-path inequality, and the matched 0/+1% clamped comparison still
requires each condition's own initial-well subtraction. A higher absolute
bottleneck energy can accompany a lower activation barrier if the well rises
more. Appendix B gives an elementary counterexample, not hafnia evidence or
a new theoretical contribution. This is why the endpoint-response controls
and the separate absolute escape-barrier test cannot be replaced by a generic
statement that clamping raises barriers.

In the registered generalized-coordinate metric, a prescribed active-space
projector separates the allowed physical gradient from support reactions.
Ordinary NEB residuals, allowed physical gradients and excluded reactions
are reported separately. Only enforced mechanical constraints enter this
projection; an unmeasured but physically movable mode is not relabelled as
a support reaction merely because a reduced model omits it.

A retained coordinate $q$ can specify one mode or a declared combination of
modes. A frozen slice keeps the remaining coordinates $r$ fixed. A conditional
surface instead follows a particular local minimum in the permitted r
space at each $q$; the lower envelope over branches is a third object.
With $q$ and $r$ measured from the same expansion centre and stable $H_{rr}$, local
quadratic release gives

$$
r^*(q)=-H_{rr}^{-1}(g_r+H_{rq}q),\qquad
K_{\mathrm{eff}}=H_{qq}-H_{qr}H_{rr}^{-1}H_{rq}.
$$

The offset from a nonzero $g_r$ is retained. The Schur complement is standard
mathematics, not a new theorem. An unresolved or unstable eliminated block
invalidates this release; a physical unstable direction cannot be removed
by a pseudoinverse. Unmeasured directions remain fixed, not implicitly stable.
At a smooth constrained stationary branch, the corresponding barrier strain
derivative is the difference of bottleneck and initial-state work derivatives,
evaluated with their actual cell Jacobians and volumes. An ordinary NEB
sampled peak alone does not justify that stationary-branch formula.

### 2.3. Controlled stationary response and paired energy gaps

Stable orthogonal release and continuation of a stationary bottleneck are
different operations. After the permitted stable r block has been eliminated,
write the local energy in retained internal coordinates q and a prescribed
external coordinate t as

$$
\begin{aligned}
\mathcal E(q,t)={}&\widetilde E_0+\widetilde g_q^Tq+\widetilde g_t t\\
&+\tfrac12q^TAq+t b^Tq+\tfrac12d t^2.
\end{aligned}
$$

All quantities come from one centre, potential and allowed-space chart.
The t coordinate is controlled by the substrate, not released or counted
as an internal unstable direction. If A is resolved and nonsingular, its
specified stationary branch is

$$
q^*(t)=-A^{-1}(\widetilde g_q+b t),\qquad
C=d-b^TA^{-1}b.
$$

In the length-valued response chart, both $q$ and $t$ have Angstrom units and
$A$, $b$ and $d$ have units of $\mathrm{eV}\,\text{\AA}^{-2}$ at the cell-energy level.
Thus $dq^*/d\epsilon=-L A^{-1}b$ has Angstrom per unit strain, not the
energy-density units of an elastic stiffness. These coordinates are not
the mass-weighted amplitudes used for path projections. An acoustic
continuum relation or a fitted optical-mode slope from a different
reference/boundary does not independently validate this internal mixed
response. Same-contract probes and prospective material tests must do so;
the latter remain outstanding here.

For an internal minimum, A is positive definite and C cannot exceed d.
For an index-one saddle, A is indefinite and this curvature correction is
not sign definite. A bottleneck's controlled curvature can therefore harden
even though elimination of its stable directions softens a retained block.
This follows from standard quadratic stationary response, not a new theorem;
softening of one selected mode alone does not determine a barrier change.
An unresolved inverse is rejected rather than repaired by a pseudoinverse.
The index is established only in the represented internal space; physically
movable but unrepresented gradients remain distinct from clamped reactions.

Each centre a uses its own registered scale
$t_a=L_a(\epsilon-\epsilon_0)$. Let e_a(t) denote its stationary model energy
change relative to its raw centre energy, including all nonstationary anchor
corrections. With the independently recorded raw centre gap
$G_{\mathrm{ref}}=E_{S,\mathrm{ref}}-E_{I,\mathrm{ref}}$, define

$$
D(\epsilon_0)=\frac{G_{\mathrm{ref}}+e_S(0)-e_I(0)}{n_{\mathrm{fu}}}.
$$

The small increments $d_a=e_a'(0)t_a+C_a t_a^2/2$ avoid subtracting large
total energies to obtain response. The paired change and its derivatives are

$$
\delta D=\frac{d_S-d_I}{n_{\mathrm{fu}}},\qquad
D(\epsilon)=D(\epsilon_0)+\delta D,
$$

$$
\frac{dD}{d\epsilon}=\frac{L_S e_S'(t_S)-L_I e_I'(t_I)}{n_{\mathrm{fu}}},\qquad
\frac{d^2D}{d\epsilon^2}=\frac{L_S^2 C_S-L_I^2 C_I}{n_{\mathrm{fu}}}.
$$

Consequently neither equal numerical t values in different charts nor a
bottleneck response with the initial-well term omitted represents this
mechanical intervention. The raw gap, stationary-anchor correction and
parameter response are retained separately. A common energy-zero shift
does not affect the response; negative model gaps are not clipped into
apparent activation barriers.

D equals a physical activation barrier only when I and S are the appropriate
actual stationary initial state and connected bottleneck on the same branch.
A restricted quadratic solution cannot establish those premises, identify
the highest bottleneck of a complete path or exclude lower competing routes.
Adding $\delta D$ to an audited NEB barrier is a separate calibration choice:
its anchor mismatch must be disclosed and checked, not silently discarded.
Caller-declared parameter bounds and displacement radii reject obvious
out-of-domain use but do not prove probe-hull coverage or anharmonic accuracy.
The [stationary](../../docs/CONTROLLED_STATIONARY_BRANCH.md) and
[paired-response](../../docs/STATIONARY_GAP_RESPONSE.md) implementations and
their analytic tests are thus methodological preparation, not new hafnia
barrier data or a successful independent material prediction.

A full initial-well quadratic need not be constructed for every comparison.
When audited target endpoint energies are declared visible inputs, the
bottleneck model can instead be paired with that measured initial-well shift.
Every control receives the same endpoint information and its preparation
cost. This is endpoint-assisted bottleneck prediction, not a zero-DFT or
fully predicted well response. It does not add an unregistered initial-well
Hessian to the finite two-bottleneck probe budget. The exact input choice,
basis, domain and calibration must be frozen before complete target path
labels are read; the present module tests alone do not freeze that batch.

The local ablations now share one dataset of measured physical Hessian
actions $D=HM$, where the columns of M are orthonormal in the declared joint
coordinate metric. A control is unavailable unless its retained, released
and prescribed-control directions lie in this measured span; an unmeasured
complement is unknown, not a stable zero block. Frozen, atomic-only and
joint-release controls retain the same affine gradient. Atomic-only release
is unavailable if the retained coordinates already contain cell motion.
The local stability-aware arm promotes the complete released eigensubspace
at or below the registered resolution floor, then requires a resolved
combined internal block of the specified index. This is training-only
increment-dimension continuation, not an anharmonic branch search. When
joint release is already valid, promotion changes no physical stationary
response and cannot establish an accuracy advantage. Measured coverage,
strong reference registration and material validation remain separate gates.

For the same complete quadratic and admissible internal space, valid joint
stationarity is also invariant under changes of the retained subspace and its
dimension. Direct full-space solves verify this null across five partitions,
including a correctly refused unstable-release case restored by promotion.
This prevents a coordinate-partition gain from being claimed as physical
prediction improvement. Independent nonharmonic, branch or equal-accuracy
cost evidence is needed to establish an actual increment; a surface at fixed
retained coordinates is a different object. The
[partition check](../../benchmarks/hfo2_channels/20261008/partition_null_E057_20261010/README.md)
is not a hafnia material forecast.

### 2.4. Fixed material contract and prospective comparisons

The production model is periodic Hf4O8, twelve atoms and four formula units,
with zero external pressure and field. Calculations use ABACUS/PBE, the
original Hf/O pseudopotentials, full 10-au DZP orbitals generated at 100 Ry,
ecutwfc=100 Ry and the original Gamma-centred 2x2x2 mesh. INPUT, KPT,
pseudopotentials and orbital bytes are bound to each observation by hashes.
No electronic setting is retuned between channels or mechanical conditions.
The unchanged electronic input sets scf_thr=1e-8, scf_nmax=150 and
mixing_beta=0.3. Byte-exact [INPUT](../../benchmarks/hfo2_channels/20261008/manuscript_local_reduction_E162_20261015/INPUT)
and [KPT](../../benchmarks/hfo2_channels/20261008/manuscript_local_reduction_E162_20261015/KPT)
snapshots retain the literal settings rather than inferred defaults.
Ordinary NEB uses 0.10 eV/Angstrom, without climbing images; cached fixed
endpoints are not recalculated at every iteration. The repaired T-PO band
has ten total images, while new channel bands have nine total, seven moving.
A residual pass does not supply a barrier error or a saddle certificate.
The preserving candidate subsequently receives six same-contract static
sampling checks, producing a fifteen-image cached reconstruction without
an additional optimization or endpoint recalculation (Appendix C.3). The M
edge requires a twelve-image ordinary refinement after its sampling audit;
this changes resolution on the same edge, not the declared channel set.

The subsequent fixed-plane comparison uses strain 0 and +1%, with +0.5%
reserved as an unseen in-range condition. Predictions and their visible
inputs are frozen before complete held-out path labels are inspected.
Controls include direct barrier interpolation, two endpoint-response nulls,
strong T/Cmma references and restricted versus released reductions.
Unavailable predictions and unstable branches are reported as abstentions.
Actual SCF calls and allocation core-hours, including preparation, define
cost; faster optimization to a different mechanism is not same-path
acceleration. The [detailed methods](METHODS_DRAFT.md) retain unit, gauge,
stability, selection and uncertainty definitions.

For the zero-strain T-escape and preserving-switching sampled peaks, the
local length-metric chart contains 36 translation-free internal coordinates:
33 atomic coordinates and three components of the open cell row. The
prescribed substrate control is separate. Its length scale is
$L=5.12968067458423$ Angstrom; this chart is not the mass-weighted atomic
projection used in the representation figures. Central differences of
the physical energy gradient use steps 0.01 and 0.02 Angstrom in each of
37 internal-plus-control directions, with two signs per step. The cached
center is reused: each center requires 148 new SCFs, 296 in total. Control
probes use $\delta\epsilon=h/L$, not the reserved physical +0.5% condition.
Raw directional responses are retained; spectra below use the symmetrized
internal block only. Finite-step disagreement and symmetrization changes
are diagnostics, not physical error bounds or a stationary-index certificate.

For a bounded conditional-relaxation comparison, the two negative
eigendirections of the first-step internal matrix define fixed retained
coordinates. The remaining 34 directions are released. Both optimizer
arms start at the same cached center with the same retained coordinates,
released-gradient norm target 0.03 eV/Angstrom, initial trust radius
0.05 Angstrom, internal-displacement bound 0.5 Angstrom and minimum atomic
separation 1.6 Angstrom. Each arm has a conservative limit of twelve new
SCF attempts, including unsuccessful calls; cache hits are free. One arm
uses the default identity initialization. The other uses the measured
Hessian restricted to a metric-orthonormal released basis, $U^T H U$,
as a standard quasi-Newton initialization. Cholesky factorization must
succeed without clipping eigenvalues. This algebraic check does not
certify stability after relaxation. The conditional norm target is not
the ordinary NEB residual threshold. Preparation and all trial calls are
charged separately; a budget-exhausted arm supplies a censored cost, not
a completed denominator for a speedup factor.

### 2.5. Signed path work as an energy-consistency diagnostic

Displacement and instantaneous gradient coverage do not measure the energy
change accumulated along a path. We therefore retain a separate, signed
physical-work check on the sampled coordinate polygon. For successive
coordinates $x_i$ and their conjugate physical energy gradients $g_i$, define

$$
W_i=\frac{g_i+g_{i+1}}{2}\cdot(x_{i+1}-x_i),\qquad
R_i=E_{i+1}-E_i-W_i.
$$

Within each declared coordinate block, orthogonal retained, omitted and
translation projectors partition this trapezoidal work exactly. Atomic
mass coordinates and the original permitted open-cell coordinates retain
their own dual gradients and units; no basis mixes these blocks. The cell
contribution completes the physical contraction. Rotating a basis inside
one span changes no result, but choosing another span changes its work.
Neither spring nor projected NEB forces replace the physical gradient.

These contributions are path-dependent signed work, not state-function
energies assigned to phonon modes. Negative work is retained. The cumulative
residual $\sum_i R_i$ diagnoses sampled energy-gradient consistency, but
does not separate unsampled nonlinearity from force/stress error or bound
the activation-barrier error. In particular, cancellation can make final
endpoint closure small while the bottleneck-region discrepancy is large.
The [implementation and finite-data checks](../../docs/SAMPLED_PATH_WORK.md)
introduce no new calculator evaluations or reduction-model accuracy claim.

### 2.6. Work between the two seen strain conditions

To distinguish imposed control from internal configuration changes, the
cross-strain diagnostic uses a single ordered zero-strain PO+ chart
$y=(u,a,t)$. Atomic displacements u, the three open-row coordinates a and
$t=L_\epsilon\epsilon$ are length-valued; $L_\epsilon=5.12968067458423$
Angstrom. They are not the mass-weighted phonon amplitudes in the projection
figures. The first two cell rows follow the prescribed substrate scaling;
the original third row remains open. Physical gradients retain both internal
forces and substrate reaction stress.

For equal-index native frames at the two training conditions, define

$$
W_i^{\mathrm{cross}}=\tfrac12[g_i(0)+g_i(0.01)]^T
 [y_i(0.01)-y_i(0)],\qquad
w_i=W_i^{\mathrm{cross}}-W_I^{\mathrm{cross}}.
$$

The residual is $\Delta[E_i-E_I]-w_i$. Orthogonal atomic, explicit rigid-
translation, open-row and prescribed-control blocks partition w exactly.
The row-cell/fractional-coordinate connection, ordered lift, energy unit
and initial-pair subtraction are fixed for every channel; no minimum-image
convention (MIC) relabelling
or intermediate DFT geometry is introduced. Equal image index is an explicit
discrete pairing, not a canonical correspondence of continuous MEPs.
The observed peak indices must coincide before this energy contrast is
identified with the difference of sampled maxima.

These endpoint trapezoids are connection-dependent work, not unique mode
energies, a frozen-slice measurement or a stationary barrier derivative.
Their closure residual is not a total DFT or activation-barrier uncertainty
bound. The two-configuration decomposition tests what a stress-only account
omits without assigning a causal microscopic coupling.

### 2.7. Atomic/open-cell restrictions and finite-coordinate model checks

The same measured internal matrix is restricted to the 33-dimensional
translation-free atomic space and the three-dimensional open-cell space
to distinguish their raw instabilities. The prescribed substrate is not
part of either internal instability count. Let $P_a$ project onto the
open-cell space in the fixed length metric and $P_-$ onto the entire
negative eigenspace of the symmetrized internal matrix. We report
$W_a=\operatorname{Tr}(P_aP_-)$ rather than assigning individual vectors
across nearly degenerate eigenvalues. This trace is invariant under
orthogonal re-expression of these spaces, including negative-space basis
rotations, but not arbitrary changes of the cell length scale. It is
neither an energy fraction nor a phonon amplitude. A raw-positive cell
block permits a formal cell-stationary Schur reduction; its nonzero
cell gradient and constant energy offset are retained, as in Section 2.2.
Raw signs do not supply a measured stability-resolution floor.

We separately test whether an unchanged original-center quadratic describes
finite-coordinate DFT candidates. Both finite-difference models use the
same retained plane selected from the first-step matrix; no coefficients
are fitted to the candidate energies. To isolate shape disagreement from
the origin offset, we evaluate

$$
d(q)=[E_{\mathrm{model}}(q)-E_{\mathrm{model}}(0)]
     -[E_{\mathrm{candidate}}(q)-E_{\mathrm{candidate}}(0)].
$$

Subtracting the already seen candidate origin makes this a retrospective
shape diagnostic, not a calibration-free or independent prediction. The
common six-coordinate cohort is $(0,0)$, $(-0.1,0)$, $(0,-0.1)$,
$(0,0.1)$, $(0.1,0)$ and $(-0.2,0)$ Angstrom. The audited source snapshot
contains nine T and six preserving candidates, all retained in the analysis.
Equal numeric coordinates at the two centers do not identify the same
physical structure. Finite released residuals also preclude treating
retained partial gradients as exact conditional-envelope derivatives.

For the separate joint length chart, a fixed orthonormal partition
$x=x_0+V_q q+V_r r$ permits the work check of Section 2.5 on neighboring
conditional candidates. With $g_q=V_q^Tg$ and $g_r=V_r^Tg$, the endpoint
trapezoids obey

$$
\Delta E-W_q=W_r+\mathcal R,\qquad
W_\alpha=\tfrac12(g_{\alpha,i}+g_{\alpha,j})^T\Delta\alpha,
\quad\alpha\in\{q,r\}.
$$

Here $\mathcal R=\Delta E-W_q-W_r$ is the full endpoint quadrature
discrepancy, not an energy-error bound. The Cauchy inequality bounds only
the released endpoint term,
$|W_r|\leq(\|g_{r,i}\|+\|g_{r,j}\|)\|\Delta r\|/2$; it does not bound
unsampled energy variation. We retain this term at finite convergence
rather than setting it to zero by an exact-envelope assumption. Closed
grid loops provide a complementary diagnostic: their raw energy changes
telescope to zero, whereas finite endpoint-gradient quadrature need not.
Neither a nonzero circulation nor a large neighboring displacement alone
establishes nonconservative DFT forces or a branch discontinuity.

The matched frozen-R comparison separates two parts of an origin-aligned
shape discrepancy without refitting either quadratic. For each original
center, R=0 at four axial coordinates $q=(-0.1,0)$, $(0,-0.1)$, $(0,0.1)$
and $(0.1,0)$ Angstrom supplies the frozen DFT counterpart. Both original finite-difference
models use the same retained plane and affine gradient. Writing $d_F$ and
$d_R$ for the frozen and released shape discrepancies, respectively, define

$$
\delta_{\mathrm{rel}}(q)=[E_R(q)-E_F(q)]-[E_R(0)-E_F(0)],\qquad
d_R(q)=d_F(q)+\delta_{\mathrm{rel}}^{\mathrm{model}}(q)
                 -\delta_{\mathrm{rel}}^{\mathrm{DFT}}(q).
$$

Every term remains signed and is evaluated at the same q. The original
raw and released origins are exactly reused. Because Q includes open-cell
motion, R=0 is not a fully fixed-cell experiment. Positive released-minus-
frozen candidate differences would be retained rather than clipped; a
finite-gradient candidate is not assumed to be a global conditional minimum.

## 3. Results: competing responses under one mechanical intervention

### 3.1. Opposite candidate responses do not imply easier switching

All eight registered clamped training bands passed the unchanged ordinary
0.10 eV/Angstrom residual criterion. Each substrate condition used one
identical ordered PO+ starting cache for all four candidates; energies at
different conditions were referred to their own relaxed PO+ wells. Table 1
therefore compares the same four registered correspondences, not selected
snapshots with different initial-state offsets. The full forward/reverse
values and final-well differences are retained in Appendix D, Table A4.

Table 1. Sampled forward peaks under the matched boundary, in meV/HfO2.
The last column is the +1% value minus the zero-strain value, not a derivative.

| Candidate | Zero strain | +1% strain | Change |
|---|---:|---:|---:|
| PO-to-M | 98.458122 | 97.227719 | -1.230403 |
| PO-to-T | 84.416624 | 96.763650 | 12.347027 |
| Preserving flip | 96.828183 | 125.424431 | 28.596248 |
| Reversing flip | 338.968965 | 294.944634 | -44.024331 |

The preserving and reversing switching peaks changed by opposite signs.
The decrease of the reversing peak by 44.024331 meV/HfO2 did not make it the
lower examined switching candidate: it remained above the preserving peak
at both training conditions. The lower sampled switching peak increased by
28.596248, whereas the lower sampled escape peak increased by 12.347027.
The sampled selectivity proxy, escape minus switching peak, consequently
changed from -12.411559 to -28.660781 meV/HfO2. Its decrease of 16.249222 is
not the favorable signature proposed in H1.

This is an observed candidate-set contrast, not an error-certified rejection
of H1 or a global pathway ordering. The +1% T/M escape peaks differ by only
0.464068 meV/HfO2, and ordinary residuals do not bound sampled-barrier error
or certify stationary bottlenecks. Both escape candidates and both ordered
switching candidates remain in the analysis; the decreasing but high
reversing candidate is not used alone to claim improved switching.

Final-well response alone also failed to reproduce these sampled changes.
Table 2 gives the preregistered final-following B1 null increments. The
alternative fixed-absolute-bottleneck null gives -6.806502 meV/HfO2 for
every channel because their initial well is common. For PO-to-T the
final-following increment is +1.785908, compared with the observed +12.347027;
for PO-to-M it is -22.206479, compared with -1.230403. These are seen training
comparisons, not independent predictive wins or error-resolved null rejection.

Table 2. Seen final-well shifts and final-following B1 barrier-response
increments, in meV/HfO2. Neither column is a measured bottleneck slope.

| Final endpoint representation | Observed well-energy shift | B1 final-following barrier increment |
|---|---:|---:|
| T descendant | +8.592411 | +1.785908 |
| M | -15.399977 | -22.206479 |
| PO-, T pattern preserving | +6.841761 | +0.035259 |
| PO-, T pattern reversing | +6.841761 | +0.035259 |

The [complete native training replay](../../benchmarks/hfo2_channels/20261008/clamped_plus1_reversing_terminal_E148_20261014/README.md)
retains 72 frame records and 66 distinct raw-log/geometry identities,
including the repeated common initial caches. They are not independent
statistical replicates. The held-out +0.5% physical condition remains
unseen; this training matrix does not establish H2.

### 3.2. Clamping changes the sampled energetic role without changing the group label

The preserving paths offer a direct check of the earlier low-energy-centre
caveat (Appendix C.2). Reanalysis of the original free-cell step-69 band and
the zero-strain clamped terminal band identifies both central samples as
Pbcn at each fixed 0.001/0.01/0.05 Angstrom tolerance and 1-degree angular
tolerance. No structure is standardized, remapped or relaxed in this analysis.
Nevertheless, their sampled energetic roles differ (Figure 1). The free-cell
centre is -14.838201 meV/HfO2 relative to its own PO+ well and lies below
both adjacent samples. The clamped centre is +96.828183 meV/HfO2 relative
to its different, boundary-relaxed PO+ well and lies above both adjacent
samples. The low central region is therefore not present on this clamped
band, even though the same group label remains. This does not establish
that every possible clamped Pbcn branch has disappeared or that either
centre is a stationary minimum or index-one saddle.

![Same Pbcn label but different sampled energetic role](figures/hfo2_preserving_boundary_role_E102_20261012/hfo2_preserving_boundary_role.pdf)

Figure 1. Preserving-pattern candidates with free cell (blue solid) and
zero-strain clamp (orange dashed), using the unchanged electronic contract.
(a) Native energies relative to each band's own PO+ initial well, not one
common initial configuration. (b) Stretch $F_{xx}$ in the common original T
chart. Hollow squares identify the two Pbcn central samples; they are not
certified TSs or stable intermediates. Nine existing samples per path are
shown with straight connectors, with no smooth-MEP fit, statistical error
band or independent replicate interpretation. The original free-cell
nine-point sampled peak is 32.806023 meV/HfO2; the already measured sampling
supplement raises its union maximum to 38.508448, which is not plotted as
an additional optimized chain. Source CSV, the full tolerance sweeps and
native-audit provenance are provided in the dated analysis case.

The geometries provide an associated boundary diagnostic. In the original
T chart, the free-cell centre has $F_{xx}=1.122022$; clamping fixes $F_{xx}$ and $F_{zz}$
to unity throughout the zero-strain band while allowing the original y
direction and two corresponding tilts to respond. The atomic arrangement
and allowed-cell response are both optimized along each band. Thus this
comparison is neither a frozen-geometry energy difference nor an isolated
strain derivative at one state. The six electronic input hashes are the
same, but the mechanical boundaries, initial wells and original generalized
metrics differ. Each band keeps its native normalized arc; equal abscissae
are not a structural matching operation.

This comparison sharpens the material question: predicting competition
requires tracking a branch's energetic role under a defined mechanical
intervention, not carrying a phase symbol or a geometric projection across
boundaries. It supports testing that question, not a demonstrated advantage
of the proposed response model. Conditional branch stability, certified
bottlenecks and independent error-bounded forecasts remain unmeasured.

### 3.3. Geometry, gradients and cell response require different diagnostics

A compact displacement projection did not imply a compact physical-gradient
representation. At the zero-strain T escape peak the geometric T triplet
covered 99.0957% of squared atomic displacement but only 14.3013% of squared
internal atomic gradient. For the M peak the corresponding fractions were
57.9237% and 9.0282%. All four registered rank-four Cmma frames retained
93.5085--93.5087% of the M gradient and 54.0718--54.0765% of the T gradient.
Low displacement coverage therefore did not disqualify Cmma as a strong
energy-model control. The rank-three triplet and the rank-four T optical
and Cmma spaces were not treated as dimension-matched models.

Figure 2 separates energy, displacement and gradient profiles in the same
original-T mass chart. The twenty-seven records include cached endpoints,
not twenty-seven independent replicates or new SCFs. Cell strains and
physical-gradient norms are kept separate from atomic coverage fractions.
The full eight-band peak coverage table is retained in Appendix D, Table A5.

![Matched-substrate escape and switching profiles with distinct displacement/gradient coverage](figures/hfo2_clamped_TMR_E099_20261011/hfo2_clamped_TMR.pdf)

Figure 2. Three ordinary-converged zero-strain training candidates: PO-to-T
(left), PO-to-M (middle) and reversing-pattern switching (right).
(a--c) Nine P=0 energy samples per chain on one energy scale, relative to
the same cached PO+ state. (d--f) Squared translation-free displacement
coverage in the common T chart. (g--i) Squared dual atomic-gradient coverage,
with different denominators from (d--f). Blue, grey and teal denote the T
geometric triplet (rank3), two complete T optical doublets (rank4) and all
four Cmma registrations (each rank4), respectively. Orange dashed guides
mark sampled peaks, not certified TSs; connectors do not imply continuous
MEP interpolation. The 27 records are descriptive samples, not independent
replicates; registration spread is not a statistical uncertainty interval.

The strain comparison retained the same reference metric, original T optical
doublets, all four Cmma registrations and continuous ordered lifts (Figure 3).
At the T candidate's peak the triplet displacement fraction remained
99.0957%/99.0045% at zero/+1%, while its gradient fraction changed from
14.3013% to 69.7245%. The corresponding Cmma fractions changed from
54.0718--54.0765% to 10.5750--10.5766%. The internal gradient norms were
0.010196606 and 0.038221369 eV/(sqrt(amu) Angstrom); the percentages alone
are neither model errors nor a reference-transferability test.

![Same-reference PO-to-T profiles at the two registered training strains](figures/hfo2_clamped_T_strain_E089_20261011/hfo2_clamped_T_strain.pdf)

Figure 3. Zero (left) and +1% (right) substrate-strain training bands.
(a,b) Energies relative to each matched PO+ endpoint. (c,d) Squared
displacement coverage in the same unstrained-T chart. (e,f) Squared dual
atomic-gradient coverage. All four Cmma registrations are retained, not
an uncertainty band. Guides mark highest sampled images; connectors
are not a fitted MEP or pointwise inter-strain correspondence.

The M pair provides a contrasting strong-reference diagnostic (Figure 4).
At +1% its T-triplet displacement and gradient fractions were 54.9251% and
3.9008%; Cmma retained 91.5873--91.5896% of the gradient. Its internal norm
was 0.035561531 eV/(sqrt(amu) Angstrom). Thus the poor geometric-T gradient
coverage cannot establish superiority over the stronger fixed reference.

![Paired M escape bands in unchanged reference spaces](figures/hfo2_clamped_M_strain_E105_20261012/hfo2_clamped_M_strain.png)

Figure 4. M escape training bands at zero (left) and +1% (right) substrate
strain. (a,b) Native energies relative to each matched PO+ well. (c,d)
Squared translation-free atomic-displacement coverage in the same original
T mass chart. (e,f) Squared atomic-gradient coverage in the corresponding
dual chart. Blue circles denote the rank-three geometric T triplet, grey
squares the two complete rank-four T optical doublets, and green dotted
curves all four preregistered rank-four Cmma frames. Vertical dashed lines
mark discrete highest images, not certified TS locations. Each band has its
own normalized generalized arc; equal abscissae do not match configurations.
The eighteen rows include shared cached endpoints and are not independent
replicates. Cmma-frame spread is not a statistical interval. No smoothing,
fitting, new mode selection or physical input change is applied. Scalar
source data and editable vector exports accompany the
[figure contract](../../benchmarks/hfo2_channels/20261008/clamped_M_strain_reference_E105_20261012/figure_contract.md).

The two ordered flip candidates probed different regions of the same
landscape (Figure 5), not the globally easiest homogeneous switching mechanism.
At zero strain their raw peaks were 96.828183 and 338.968965 meV/HfO2.
The preserving peak is not a certified index-one saddle. At that peak the
triplet covered 82.6301% of displacement but 12.9583% of gradient; the
complete optical rank-four control was below 0.0001% for both, whereas
Cmma retained 52.7203% of the gradient. At the reversing peak the triplet
gradient fraction was 0.0320%, compared with 65.1046% for Cmma. Parent-equivalent
minus endpoints remain different at ordered indices; their phase equivalence
does not establish whole-path equivalence or electronic winding.

![Preserving and reversing ordered flips with unchanged atomic reference spaces](figures/hfo2_clamped_flips_E101_20261011/hfo2_clamped_flips.pdf)

Figure 5. Two ordinary-converged zero-strain ordered flip candidates,
preserving (left) and reversing (right). (a,b) Native sampled energies,
with the same PO+ zero and energy scale. (c,d) Squared displacement
coverage in the common T chart. (e,f) Squared dual atomic-gradient coverage
in that chart. Both T reference ranks and all four Cmma registrations are
retained. Orange guides mark highest sampled images, not certified TSs.
Each column has its own normalized generalized arc; equal abscissae are
not matched configurations. Eighteen frame records share the initial
native cache and are not independent replicates. Straight connectors are
not a smooth-MEP reconstruction, and frame spread is not an uncertainty band.

Eight of eight registered G2 training bands are ordinary-converged. These
descriptive projections neither select a reference winner nor certify
conditional energies, transverse stability or independent prediction.
Appendix D retains every reference rank, all four tied Cmma frames and
their nonstatistical registration ranges. The original source modes and
forces are not modified by ASR; raw translation contributions remain visible.

Direct joint-curvature measurements reveal an additional limitation of a
one-coordinate description at the two sampled peaks. At both finite steps,
each internal matrix has two negative eigenvalues; removing the internal
path tangent still leaves one negative direction (Table A6). For T escape,
the lowest two values are approximately -0.276 and -0.165 eV/Angstrom squared
at the smaller step. For preserving switching they are -1.272 and -0.447.
These are nonstationary, finite-step curvatures in a joint atomic-cell
length metric, not phonon frequencies or certified saddle indices. For
these measured matrices, retaining only one internal direction cannot
leave a positive released complement. This is a property of the measured
quadratic, not a new reduction theorem.

The registered reference spaces must be tested against this same internal
problem rather than displacement coverage alone. The rank-three T triplet,
rank-four optical control and each of the four rank-four Cmma registrations
leave negative released curvature at both centers and steps (Table A7).
However, retaining all 33 T optical directions leaves a positive three-cell-
direction complement, while its retained Schur block still has two negative
directions. Thus the evidence does not show that all fixed references fail.
The ranks and released dimensions differ, and a complete quadratic response
remains invariant under valid basis partitions as described in Section 2.3.

Releasing the remaining coordinates at fixed two-dimensional coordinates
lowers the T-center energy by 0.00065214 eV per twelve-atom cell and the
preserving-center energy by 0.04948493 eV per cell. The curvature-seeded
arms reach released-gradient norms 0.02630 and 0.01959 eV/Angstrom in two
and five new SCFs, respectively; neither identity arm reaches the common
target within twelve attempts (Table A8). This is a bounded local result,
not a whole-NEB acceleration benchmark. The preceding 296 curvature SCFs
and earlier path costs are additional, so no net cost advantage follows
from these call counts alone. Post-relaxation released stability remains
unmeasured.

The two resulting candidate energies differ by only 0.00081345 eV per
cell, but energy proximity does not establish a common mechanism or branch.
They retain a 0.19178-Angstrom translation-free distance in the common
ordered T chart, without atom rematching or rotation. Moreover, each
retained plane is anchored at its own original center. These local
fixed-coordinate results therefore do not replace any sampled G2 barrier,
establish branch merger or certify a transition structure. They motivate
checking conditional equilibration and branch stability before using a
low-dimensional surface for prediction.

### 3.4. Endpoint closure can conceal a bottleneck-region discrepancy

The signed-work diagnostic uses the eight completed clamped training bands,
not earlier unfinished snapshots or the reserved strain condition. All nine raw
E/F/stress frames per band are reparsed, with the same ordered continuous
lift, original T mass chart, permitted open-cell space and electronic
contract. These 72 frame records contain 66 distinct raw-log/geometry
identities: four zero-strain paths share one initial PO+ cache, and the
four +1% paths share another. They are
neither new SCFs nor independent statistical replicates. Every energy
change is relative to that band's own matched PO+ well.

Table 3 compares raw sampled peak energies with the full endpoint-gradient
trapezoid accumulated to the same image. Peak and endpoint residuals use
the same sign, raw energy change minus accumulated work. Units are
meV/HfO2 throughout; the Hf4O8-cell values are multiplied by 250.

| Training band | Raw peak | Work to peak | Residual at peak | Residual at endpoint |
|---|---:|---:|---:|---:|
| PO-to-M, zero strain | 98.458122 | 77.560379 | 20.897743 | 3.006016 |
| PO-to-T, zero strain | 84.416624 | 75.955937 | 8.460687 | 0.728018 |
| Reversing flip, zero strain | 338.968965 | 303.693643 | 35.275323 | 0.115373 |
| PO-to-T, +1% strain | 96.763650 | 87.345393 | 9.418257 | 1.997185 |
| Preserving flip, zero strain | 96.828183 | 82.304235 | 14.523948 | 0.032398 |
| PO-to-M, +1% strain | 97.227719 | 75.922579 | 21.305139 | 2.658882 |
| Preserving flip, +1% strain | 125.424431 | 111.393911 | 14.030520 | 0.000050 |
| Reversing flip, +1% strain | 294.944634 | 259.087641 | 35.856993 | -0.000001 |

The reversing candidate illustrates the distinction most clearly. Its
endpoint closure is only 0.115373 meV/HfO2, whereas the discrepancy at the
highest sample is 35.275323. Small final closure therefore does not certify
energy-resolved consistency through the transition region. Conversely,
the present nine-sample trapezoid cannot determine whether this discrepancy
comes from coarse quadrature, physical force/stress error or both. It is
not a measured barrier uncertainty and does not justify retuning the
electronic inputs or declaring an ordinary-converged band failed.

The coordinate pullback itself passes a distinct implementation check:
direct native Cartesian-force/ASE-stress contractions agree with the
dual atomic/open-cell work within $1.12\times10^{-16}$ eV per cell. This
algebraic agreement does not test DFT accuracy. Under the T geometric
triplet, signed retained work to the respective peaks is -0.111271,
2.666921, 195.447066, 1.367408, -0.165926, 0.011623, 3.753494 and 147.404416 meV/HfO2, in Table 3 order.
Omitted atomic work is 68.756917, 72.789643, 105.539695, 85.773095, 82.500795,
67.611606, 107.167048 and 108.609913; open-cell work is 8.914733, 0.499373, 2.706882, 0.204891, -0.030633,
8.299349, 0.473369 and 3.073313. The explicit rigid
translation terms are also retained, at most $4.63\times10^{-9}$ in these
units at the displayed peaks, rather than removed by a source ASR repair.
These signed terms sum to the full sampled work, not to a unique energetic
allocation among modes.

Strong references remain in the comparison. The two complete T optical
doublets yield retained work 12.359091, 11.900861, 0.782574, 11.007299, 18.977823,
11.650637, 26.369936 and 6.179305
in the same order. Across all four registered Cmma frames, retained-work
ranges are 55.618955 to 55.619265 for zero-strain M;
62.066996 to 62.068212 for zero-strain T;
63.846564 to 63.846572 for reversing;
67.984207 to 67.985028 for +1% T;
47.684643 to 47.684652 for preserving;
55.402398 to 55.403003 for +1% M;
26.439087 to 26.439089 for +1% preserving; and
23.010446 to 23.010450 for +1% reversing. Each representation includes
its full omitted atomic contribution and the same open-cell work; ranks
three and four are not treated as identical-dimensional models. This
reference dependence is not a conditional-energy or forecast-accuracy
ranking. In particular, the reversing triplet's integrated work and its
small instantaneous peak-gradient fraction (Table A5) answer different
questions and need not track one another.

The [four original scalar reports](../../benchmarks/hfo2_channels/20261008/clamped_path_work_E092_20261011/README.md),
[preserving report](../../benchmarks/hfo2_channels/20261008/clamped_preserving_reference_E101_20261011/README.md),
[+1% M report](../../benchmarks/hfo2_channels/20261008/clamped_M_strain_reference_E105_20261012/README.md)
and [eight-band raw-frame replay and cohort check](../../benchmarks/hfo2_channels/20261008/clamped_plus1_reversing_terminal_E148_20261014/README.md)
thus add an energy-gradient consistency test to geometric coverage, without
selecting G3 bottlenecks, a reference winner or unseen labels. Measured
conditional response, stability and the prospectively frozen strong-control
comparison remain necessary for a predictive claim. This diagnostic uses
existing data only and does not open an additional sampling or Hessian budget.

### 3.5. Cross-strain work separates imposed control from configuration change

The second work diagnostic pairs equal image indices between the two
training conditions in the single controlled-coordinate chart of Section 2.6.
All four pairs retained their sampled peak index (M:3; T and both flips:4).
The matched-peak energy changes therefore reproduce Table 1 exactly, while
the signed terms in Table 4 partition the declared endpoint trapezoid.

Table 4. Cross-strain peak-relative work in meV/HfO2. Terms are relative to
the same initial-pair work; the final column is energy change minus full
trapezoidal work, not a barrier-error estimate. The retained rigid-translation
terms are smaller than 2.411e-10 at these peaks.

| Candidate | Observed peak change | Atomic work | Open-row work | Substrate work | Energy minus work |
|---|---:|---:|---:|---:|---:|
| PO-to-M | -1.230403 | 2.720687 | -0.174328 | -3.835247 | 0.058485 |
| PO-to-T | 12.347027 | -0.228947 | 0.220068 | 12.275095 | 0.080811 |
| Preserving flip | 28.596248 | 36.895288 | -0.078916 | -12.858888 | 4.638764 |
| Reversing flip | -44.024331 | -32.332083 | -0.270001 | 0.072814 | -11.495061 |

In this connection, the T candidate's work was dominated by the prescribed
substrate term. For preserving, that term had the opposite sign from the
observed peak change: +36.895288 atomic work and a +4.638764 residual
accompanied the -12.858888 substrate term. The isolated mechanical partial
therefore cannot replace the observed finite path response. It is not a
measurement that a frozen strain slice would lower the barrier: the endpoint
gradients here belong to two different internally optimized configurations.

The reversing residual was -11.495061 meV/HfO2. It is preserved rather than
attributed solely to DFT noise or a particular coupling. No intermediate
cross-strain geometries were calculated. Independent Cartesian force/stress
and dual-chart contractions agreed within 6.662e-16 eV/cell, an algebraic
implementation check rather than physical accuracy evidence. The
[36-pair source analysis](../../benchmarks/hfo2_channels/20261008/cross_strain_work_E151_20261014/README.md)
reuses the same 72 frame records/66 distinct identities, without new SCFs,
held-out inputs or reference selection. The next inference requires measured
local curvatures, stable conditional branches and frozen independent tests.

### 3.6. The extra negative direction is already present in the atomic sector

The additional raw negative direction at each sampled peak is not created
solely by opening the cell. Both atomic restrictions already had two
negative eigenvalues at both finite steps, as did the joint matrices;
their open-cell restrictions had three positive eigenvalues (Table 5).
For these measured quadratics, coupling changed the negative curvatures
without adding another negative direction. This observation is specific
to the two nonstationary centers, not a general rule excluding coupling-induced
instabilities or an uncertainty-resolved stationary index.

Table 5. Atomic and cell restrictions of the original sampled-peak matrices.
Steps are in Angstrom and the three curvature columns in eV/Angstrom
squared. The last column is the dimensionless open-cell trace of the
two-dimensional negative subspace, not a percentage of a barrier energy.
The same length scale L=5.12968067458423 Angstrom is used in every row.

| Center | Step | Minimum atomic | Minimum joint | Minimum cell | Cell trace |
|---|---:|---:|---:|---:|---:|
| T | 0.01 | -0.240169 | -0.276193 | 2.956214 | 0.0205773 |
| T | 0.02 | -0.237806 | -0.273596 | 2.956557 | 0.0206105 |
| Preserving | 0.01 | -1.271818 | -1.271819 | 3.090798 | 1.04022e-06 |
| Preserving | 0.02 | -1.281815 | -1.281821 | 3.092172 | 3.09248e-06 |

The T negative subspace contained more open-cell admixture than the
preserving subspace in this metric. However, the complete atomic-cell
coupling operator norms were approximately 7.72 and 8.10 eV/Angstrom
squared, respectively, at the smaller step. A large total coupling norm
therefore did not identify a strongly cell-mixed unstable direction.
Conversely, small open-cell weight at one sampled configuration cannot
establish insensitivity of a complete barrier to the prescribed substrate.
The substrate response and the permitted open-cell fluctuations are
different quantities. All four restrictions and the nonzero-gradient
formal reductions are retained in the
[atomic/open-cell analysis](../../benchmarks/hfo2_channels/20261008/G3_atomic_cell_signature_E164_20261015/README.md).

### 3.7. A useful curvature seed is not a finite-domain energy predictor

The bounded optimizer result in Section 3.3 did not establish the accuracy
of the original quadratic over finite retained displacements. In the
first-step eigenplane, the retained/released coupling was numerically
zero by construction; its Schur curvature correction was of order
$10^{-28}$ eV/Angstrom squared in norm. Formal release changed the energy
offset but supplied no appreciable new shape correction in that plane.
This algebraic null is not a prediction improvement.

After origin alignment, the original quadratics still disagreed with the
six common-coordinate DFT candidates (Table 6). Both probe steps are
reported, not selected according to agreement. The maximum absolute
shape discrepancies were 8.52--8.59 meV/cell for T and 20.52--20.72
meV/cell for preserving, each at $q=(-0.2,0)$ Angstrom. The corresponding
signed errors were positive for T and negative for preserving. All nine
T candidates and all six preserving candidates in the audited source
snapshot remain in the records; the common cohort is used only to make the table's
coordinate coverage explicit.

Table 6. Unfitted original quadratic versus released DFT candidates.
All three energy columns use meV per twelve-atom cell, not meV/HfO2.
Origin release is relative to each original center; the shape discrepancy
uses each model's and candidate's own seen origin. Probe steps are in
Angstrom. These are seen-training diagnostics, not independent errors or
continuous-path barrier errors.

| Center | Probe step | DFT origin release | Model origin release | Max. shape discrepancy |
|---|---:|---:|---:|---:|
| T | 0.01 | -0.652142 | -0.866958 | 8.524106 |
| T | 0.02 | -0.652142 | -0.874429 | 8.587800 |
| Preserving | 0.01 | -49.484930 | -54.521222 | 20.521839 |
| Preserving | 0.02 | -49.484930 | -54.705828 | 20.720640 |

The candidate-energy axial second differences at 0.1 Angstrom were
(-0.405128, -0.243087) eV/Angstrom squared for T and
(-0.337763, -0.184576) for preserving. They are finite differences of
candidate energies, not infinitesimal Hessians of certified conditional
minima. Nonlinear response, finite-difference effects, branch changes
and residual relaxation can all contribute to the model disagreement;
these observations do not distinguish their individual contributions.
Nor does their size define a physical uncertainty bound. The
[quadratic comparison](../../benchmarks/hfo2_channels/20261008/G3_quadratic_agreement_E165_20261015/README.md)
keeps the affine released gradient and both original probe steps.

These results separate two uses of local curvature: an initialization
that is updated during optimization, and a fixed finite-domain energy
model. The former can assist local relaxation even when the latter
misses the observed shape. A matched frozen-coordinate comparison is
needed before assigning the discrepancy to the relaxation response.
Eight such static controls, four per center, have now completed under
the unchanged electronic contract. Native energy, force, stress, input
and 32-MPI checks passed for all eight. They use the original cached
origins and seen released counterparts, reserve training slots only,
and do not replace the independent tests.

The subsequently completed T candidate sheet exposed a separate
finite-convergence issue. Its 25 points required 118 new SCFs in addition
to the cached released origin. Final released-gradient norms ranged from
0.00653 to 0.02978 eV/Angstrom, and the largest total internal displacement
was 0.46378 Angstrom within the original 0.5-Angstrom domain. Sampled
energies spanned -38.531 to 0 meV/cell relative to that released origin.
These are fixed-coordinate candidates, not an end-to-end MEP or certified
stable conditional minima.

Across the 40 nearest grid edges, the largest absolute discrepancy using
only retained endpoint work was 1.209 meV/cell. On that edge the signed
released term was 1.094 meV/cell, leaving a full quadrature discrepancy
of 0.115 meV/cell. Over all edges the largest absolute full discrepancy
was 1.039 meV/cell. Across 16 elementary closed loops, the largest absolute
retained and full endpoint circulations were 2.000 and 0.578 meV/cell,
respectively; these maxima occurred on different loops. This
[complete-sheet diagnostic](../../benchmarks/hfo2_channels/20261008/G3_sheet_consistency_E172_20261015/README.md)
retains finite released residuals and unsampled quadrature separately.
The values are not DFT, interpolation or barrier-error bounds and do not
certify smoothness or stability. They also do not identify the cause of
the original quadratic's finite-domain disagreement.

The completed preserving batch did not attain the same coordinate
coverage. All 25 coordinates were processed, but only 17, including the
cached origin, reached the common released-gradient target. The eight
remaining coordinates stayed as holes. All 172 new SCFs passed native
energy, force, stress, unchanged-input and genuine 32-MPI checks, and all
24 nonzero-coordinate optimization trajectories replayed exactly. Seven
holes exhausted twelve evaluations; one stopped after ten without reaching
the target. They consumed 94 of the 172 SCFs. The batch required 193.859
completed-call proxy core-hours and 195.120 allocation core-hours, excluding
the original curvature and cached preparation.

Cache-only replay attributed every recorded geometry rejection in these
eight trajectories to the original total-displacement radius, not to a
changed electronic input or a failed SCF. Their last completed trials lay
at radii 0.49971--0.50000 Angstrom within the registered 0.5-Angstrom domain;
none of the completed trials met the 0.03 eV/Angstrom released-gradient
target. This establishes contact with the optimizer's domain boundary,
not nonexistence of a solution, boundary stationarity or termination of a
physical branch. The [complete audit and hole replay](../../benchmarks/hfo2_channels/20261008/G3_preserving_terminal_E179_20261015/README.md)
retain the original decisions and stopping reasons. We do not fill these
holes by interpolation or equate batch completion with a stable surface.

The matched controls locate the larger preserving-channel discrepancy
in the release response, within this four-coordinate training cohort
(Table 7). Frozen shape discrepancies stayed below 0.12 meV/cell for
both original probe steps, while released discrepancies reached
4.67--4.72 meV/cell. The DFT release change relative to its origin
reached 4.607 meV/cell. At $q=(-0.1,0)$ Angstrom for the first-step model,
the signed frozen discrepancy -0.064445 and missed release-shape term
-4.606959 summed to the released discrepancy -4.671404 meV/cell.

Table 7. Matched frozen/released shape discrepancies on four axial
coordinates per original center. All energy columns are meV per
twelve-atom cell; probe steps are Angstrom. Separate column maxima may
occur at different coordinates and are not subtracted to form a
decomposition. All 16 model-point records are retained in the
[matched ablation](../../benchmarks/hfo2_channels/20261008/G3_frozen_release_E173_20261015/README.md).
These are training-only diagnostics, not independent predictions or
physical uncertainty bounds.

| Center | Probe step | Max. frozen shape | Max. released shape | Max. DFT release change |
|---|---:|---:|---:|---:|
| Preserving | 0.01 | 0.075880 | 4.671404 | 4.606959 |
| Preserving | 0.02 | 0.116535 | 4.720758 | 4.606959 |
| T | 0.01 | 0.416264 | 1.327177 | 1.168761 |
| T | 0.02 | 0.403160 | 1.346071 | 1.168761 |

Signed cancellation also matters. At T $q=(0.1,0)$ Angstrom, the
first-step frozen discrepancy -0.416264 and release-shape term
+0.378430 partly cancelled, leaving -0.037833 meV/cell. This near agreement
at one coordinate did not establish general model accuracy. The ablation
distinguishes frozen-shape mismatch from missed release-shape change,
but does not isolate anharmonicity, finite-step effects, finite released
residuals or branch changes as a unique microscopic cause. Nor does it
replace the full T/Cmma strong-control or independent-prediction matrix.
The eight statics used 8.137 completed-call proxy core-hours and
9.076 allocation core-hours; curvature and prior candidate preparation
remain additional.

## 4. Discussion and remaining material tests

The completed training comparison demonstrates why a decreasing barrier
on one examined path is insufficient evidence of easier switching. The
higher reversing candidate decreased under tensile substrate strain, while
the lower preserving candidate increased. The candidate-set selectivity
proxy moved opposite to the proposed favorable signature. This observation
must not be converted into an error-certified physical ordering: the
closely separated strained escape peaks and unmeasured continuous bottlenecks
still limit the interpretation.

The representation and work diagnostics explain what the next predictive
test must distinguish, not yet a causal mechanism. A reference can capture
a configuration well and still miss its gradient; a coarse work integral
can close at the final well but retain a large interior discrepancy.
Cross-strain pairing additionally shows that a prescribed stress partial
need not have the sign of the actual peak shift. These are separate
properties, not interchangeable measures of energy-model accuracy.

The completed local-release test makes this distinction operational:
an ordinary-converged chain and a compact displacement projection do not
guarantee equilibration of the coordinates omitted from a proposed surface.
The substantial preserving-center energy change occurs without changing
the electronic contract. It warrants a conditional-branch test, not a
retrospective change of the whole-path barrier table. Conversely, a
gradient-converged fixed-coordinate candidate is not a conditional minimum
until the released block has been checked. The local tests narrow the
next material question to whether a stable, reproducible branch supports
an independent response prediction; they do not already answer it.

The local atomic/open-cell restrictions and the finite-coordinate test
further delimit the mechanism claim. In these raw matrices the extra
negative direction was already atomic, rather than an instability added
only by releasing the open cell. Yet retaining both negative directions
did not make the original quadratic an accurate finite-domain energy
model. Neither a cell-projector trace nor an optimizer call count can
therefore replace explicit conditional sampling and stability tests.
The frozen/relaxed comparison is designed to separate two sources of
disagreement under the same Hamiltonian, not to repair the result by
changing the cutoff or to fit seen data and relabel them independent.

The completed T sheet adds a practical caution for derivative-based
surface construction: finite-residual retained partial gradients need not
form the gradient field of an exactly relaxed conditional energy. The
signed released-work term must remain visible, while endpoint quadrature
and branch validity still require separate checks. This diagnostic does
not substitute for independent surface energies or transverse stability.

Reference adaptation alone cannot establish a methodological improvement.
For the same fully represented quadratic and fixed-control stationary
problem, valid retained/released basis changes give the same response
(Section 2.3). Any advantage over the full T/Cmma controls must therefore
arise in an explicitly tested finite representation, nonlinear branch,
validity domain or equal-accuracy cost comparison, rather than in renaming
coordinates or dropping a strong reference. The existing analytic checks
enforce this null; they are not independent material forecasts.

The decisive unfinished evidence is a same-boundary conditional response
with measured transverse stability and a frozen comparison at unseen
conditions. Ordinary highest images do not become TSs when used as
curvature centers. A valid local index must exclude prescribed substrate
coordinates and retain unresolved directions as unresolved. Predictions
must then meet the registered strong-control and endpoint-information
rules at measured resolution. Added directions selected after target labels
are revised training, not successful original forecasts. Those material
tests and the intended JCTC methodology claim remain incomplete.

All current conclusions concern periodic twelve-atom, zero-field,
zero-temperature energy pathways under the fixed ABACUS/PBE contract.
They do not establish domain-wall nucleation, rates, retention times,
finite-temperature free energies or functional-independent critical strains.
The archival observations below support reproducibility and define failure
boundaries; they are not substitutes for the missing predictive evidence.

## Appendix A. Reproducibility and editorial completion checks

The archived cases provide frozen numeric observations, input/log/geometry
checksums, continuous ordered lifts, analysis commands and figure CSVs.
Raw production logs remain read-only on HF. Earlier incomplete snapshots
remain dated evidence rather than being overwritten by final values.
No external paper figure or unpublished author dataset is copied into this draft.

Before a submission-ready abstract or conclusion is written, this manuscript
requires: completed error/phase audits and resolved matched G2 responses;
same-boundary conditional-branch and transverse-stability evidence;
prospectively frozen predictions with strong-control comparison and actual
cost accounting; final numeric captions and an openly reproducible release
of the results permitted for redistribution. These are requirements, not
claims that those outcomes have been achieved. Details belong in main-text
appendices rather than a separate supplementary manuscript, as requested.

## Appendix B. Why a constraint need not increase a relaxed-endpoint barrier

Consider the dimensionless illustrative potential

$$
\mathcal V(x,y)=(x^2-1)^2+(1+x^2)y^2 .
$$

The coordinates and energy scale are arbitrary; this is not a fitted hafnia
potential, phonon-coordinate model, new DFT point or independent prediction.
The unrestricted minima are $(x,y)=(\pm1,0)$ with energy zero. Every continuous
path connecting them crosses $x=0$, where $\mathcal V=1+y^2\ge1$, and the path
$y=0$ attains this bound. The ideal unrestricted barrier is therefore one.

With the prescribed constraint $y=1$, the relaxed endpoints instead become
$x=\pm1/\sqrt2$, with energy $7/4$. Their constrained bottleneck at $x=0$
has energy two and negative curvature in the allowed x direction. The
activation barrier is consequently $2-7/4=1/4$, despite the higher absolute
bottleneck energy. The nonzero y gradient is a prescribed-coordinate reaction,
not failure of allowed-coordinate stationarity. The unrestricted endpoints
are absent from the constrained path space, so the common-endpoint inequality
has not been violated. This elementary distinction motivates reporting both
well and bottleneck responses; it supplies no favorable hafnia window or
methodological novelty. The analytic derivatives, stationary-point curvatures
and a continuous family of fixed-y examples are checked in the
[analytic checks](../../tests/test_barrier_constraint_interpretation.py).

## Appendix C. Free-cell candidates and finite sampling checks

These dated free-cell observations establish the source correspondences and
diagnostic limits, not the current state of the clamped training matrix.
The qualified G1 review subsequently permitted G2; Table 1 in the main text
contains the completed matched-boundary comparison. Free-cell energies are
not reused as clamped energies. The following source records and four
figures are retained rather than overwritten by later calculations.

### C.1. The formation band and the shared polar starting well

The repaired T-PO band passes the ordinary residual criterion with
0.0598816 eV/Angstrom. Its ten frozen images reproduce complete same-contract
energies, forces and stresses; its endpoints are unchanged. The discrete
T-to-PO maximum is 33.889 meV/f.u. Reversing the thermodynamic view, without
reoptimizing or remapping the band, gives 115.210 meV/f.u. for PO-to-T.
Their difference equals the endpoint energy change, -81.321 meV/f.u.
This is the first residual-passed candidate edge, not a sampling-resolved
absolute activation barrier or a full-variable TS.
([Numerical source](../../benchmarks/hfo2_channels/20261008/converged_gap/README.md))

![Reference projections and residuals on the ordinary T-PO band](figures/hfo2_T_PO_ordinary_20261008/hfo2_T_PO_ordinary_modes.png)

Figure A1. Audited discrete energy, geometric pattern amplitudes, Green strain,
reference-subspace weights and ordinary residuals along T-to-PO. Connections
join samples, not an independently validated smooth MEP. Reference phonons
are not local bottleneck eigenvectors or mode energy contributions; fixed
endpoints have no NEB residual. Complete captions, CSV and reproduction
instructions are in the [figure bundle](figures/hfo2_T_PO_ordinary_20261008/README.md).

All four candidate observations share the ordered PO+ energy
-9783.249675811956 eV/cell. The following table uses the complete frozen
terminal update audited by 18:12 on 9 October, not rounded live-log values:

| Candidate from PO+ | Source job/step | Residual | Maximum | Pass |
|---|---|---:|---:|---|
| PO-to-T, reverse view | 28300425/6 | 0.059882 | 115.210 | Yes |
| PO-to-M, refined | 28392675/1 | 0.096964 | 82.680 | Yes |
| T-pattern-preserving flip | 28319570/69 | 0.099372 | 32.806 | Yes |
| T-pattern-reversing flip | 28319571/45 | 0.094023 | 392.823 | Yes |

Table A1. Four ordinary-residual-passed candidates, not an error-bounded
activation-barrier ranking. Snapshot steps and times differ; M has twelve
images after the separate sampling refinement. Residuals are in eV/Angstrom and
sampled maxima in meV/f.u.; maxima include the endpoints and use four formula
units per cell. The two flip labels denote registered geometric
operations, not certified distinct MEPs or assignments to published irreps.
The [forty-record terminal audit](../../benchmarks/hfo2_channels/20261008/reversing_peak_sampling_20261009/README.md)
includes cached endpoints and reused records, not forty new or independent SCFs.
The earlier nine-image M and incomplete reversing observations remain dated
evidence in their original archives.

PO-to-M terminated normally at the unchanged force threshold. Its ordered
M endpoint is P2_1/c at all three declared symmetry tolerances. The refined
sampled forward/reverse barriers are 82.680/154.065 meV/f.u., with unchanged
reaction energy -71.385 meV/f.u.; their difference follows from the common
endpoint energies. Its initial nine-image maximum was 71.582 meV/f.u.;
the source's physical tangential force at that sampled peak was -0.383495
eV/Angstrom. This earlier diagnostic prompted the finite sampling check
and subsequent ordinary refinement described below. Neither the old nor
the refined residual pass alone establishes a stationary bottleneck.

### C.2. Nonpolar snapshots do not uniquely identify a switching mechanism

The earlier preserving snapshot at step 14 had a Pbcn central sampled
maximum across the declared 0.001/0.01/0.05 Angstrom symmetry tolerances.
The later profile split into two sampled local peaks. At the ordinary-converged
step 69, images 3 and 5 are approximately 32.806 meV/f.u. above PO+.
Both are Pca2_1 at those tolerances. The central image remains Pbcn but is
now below PO+, -14.838 meV/f.u., with an atomic/cell NEB residual of
0.062019/0.099372 eV/Angstrom. Transverse stability is unmeasured. It is neither the current
highest image nor a certified stable intermediate. A phase label or an
earlier central peak cannot therefore fix the eventual bottleneck location.
If that lower nonpolar branch is stable, it also supplies a potential escape
route from the polar well. M/T comparisons alone could then overstate
protection. This coverage limitation is retained explicitly; the registered
branch/stability tests must resolve its interpretation before a global
escape-resistance claim, without adding a post hoc endpoint matrix.

At the two side peaks, the physical generalized tangential forces are
approximately +0.242486 and -0.242483 eV/Angstrom, despite small atomic/cell
NEB residuals at those images. Those forces distinguish the ordinary
perpendicular-force pass from stationary bottlenecks. The central sampled
minimum and the side maxima do not certify a two-saddle sequence or justify
replacing the shared initial well with a new Pbcn endpoint. Its T-referenced
expansion, F_xx=1.122022, is also a prospective branch-response diagnostic:
the already registered substrate contains that direction. The existing G2
matrix must test whether that central branch persists, without adding a
post hoc phase/path or an additional Hessian to the finite budget.

The terminal reversing centre at step 45 is Pa-3 at all three tolerances;
its earlier step-32 Pbca/Pa-3/Pa-3 sensitivity remains recorded. Its amplitude
in the selected rotated-T pattern triplet is essentially zero, but vanishing
coordinates in a truncated representation do not identify a cubic structure. Its group
label in this twelve-atom cell also does not identify a twenty-four-atom
literature variant or a domain-wall motif. The PO-to-M peak retains its
P1/P2_1 tolerance dependence rather than being standardized to a preferred
label.

![Four ordinary terminal bands and physical tangential-force limitation](figures/hfo2_G1_terminal_20261009/rendered/hfo2_G1_terminal.png)

Figure A2. Complete same-PO+ terminal observations at steps 6/1/69/45, corresponding
to Table A1. (a) Discrete low-energy profiles. (b) The higher, ordinary-converged
reversing candidate is shown separately rather than omitted or clipped.
(c) Ordinary NEB max-vector residuals on the moving images.
(d) The preserving band's physical tangential-force magnitude and atomic/cell
NEB residuals are different diagnostics in the same registered source metric.
The dashed 0.10 line is the ordinary NEB target, not an independently
established stationary-TS tolerance for the physical tangent.
Lines connect calculated images only; endpoints have no NEB residual.
No final energetic hierarchy or TS certification follows from this figure.
The [figure contract](figures/hfo2_G1_terminal_20261009/README.md) and CSV
retain phase checks, original lifts and complete provenance. The
[earlier 6/39/69/32 figure](figures/hfo2_G1_preserving_pass_20261009/README.md)
and [initial dated figure](figures/hfo2_G1_provisional_20261009/README.md) are
retained with their original data and hashes, rather than overwritten.

These are mechanism candidates for further relaxation, not discoveries of
Pbcn-mediated switching or of an antipolar wall. A group symbol at one
uniform snapshot does not map it onto the motif, boundary or energy of a
published interface. Fixed and released lattice results must also remain
separate rather than being pooled under a common phase label.

The preserving band's dominant residual changes from cell rows at step 14
to atoms at steps 25/27, then back to cell rows at steps 48 and 69. Its decreasing
sampled maximum and rebound of the global force norm are compatible with
a changing coupled configuration, not proof of a persistent single-block
cause or of a particular optimizer's benefit. The preserving continuation
ended normally on the unchanged force criterion; the reversing candidate
subsequently passed as well. Whole-path polarization and transverse-stability
audits remain distinct. Four ordinary passes alone neither complete G1 nor
fix the stationary locations needed for local mode-strain predictions.

The forty terminal raw band tables are independently replayed with the
audited 96-electron, 48-occupied-band contract. Minimum sampled indirect
gaps are 4.2520, 4.1317, 4.5893 and 4.0173 eV for T/PO, M, preserving and
reversing, respectively. These Gamma2x2x2 sampled gaps are not a full-Brillouin-zone
insulating certificate or Berry-phase branch evidence. The separate native
Berry protocol retains physical $eR/V$ and spin-paired native $2eR/V$ periods;
actual-cell reduced polarization, not a fixed Cartesian period, is needed
along the changing lattice. Neither modular endpoint opposition nor a
unique lift at finite sampled images establishes an absolute spontaneous
polarization or excludes unsampled winding.

The missing R3 path measurements were then made at all fourteen existing
switching interiors, reusing the three audited PO endpoint properties.
Fourteen output-only SCFs reproduce baseline E/F/stress exactly at the
raw-log precision; this is not a claim of zero total numerical error.
Forty-two fixed-charge NSCFs use the registered $2\times2\times2$,
$2\times2\times4$ and $2\times2\times8$ longitudinal quadratures. Maximum
changes from the second to the third quadrature are 0.00063768 and
0.00083843 C/m$^2$ for the preserving and reversing candidates, below
the predeclared 0.01 C/m$^2$ gate. No NSCF energy is used in a barrier.

![Native R3 classes, conditional sampled lift and measured quadrature sensitivity](figures/hfo2_path_Berry_20261009/rendered/hfo2_path_Berry.png)

Figure A3. Native-path electronic property audit. (a) Reported classes in
actual-cell reduced units $p=P_3/(e|R_3|/V)$. (b) Conditional nearest-sample
lift retaining native period 2 and explicit initial integer gauge 0, not
the smallest-absolute-$P$ convention. (c) Changes from $2\times2\times4$
to $2\times2\times8$; the 10 mC/m$^2$ gate is stated above an axis showing
the smaller measured range.
Eighteen plotted rows comprise fourteen new interiors and cached endpoints,
with common PO+ repeated. Connections join samples, not a continuous-path
certificate. Source data and editable exports accompany the
[figure contract](figures/hfo2_path_Berry_20261009/README.md).

Both sampled lifts start at 1.188154, pass near 2 at their central nonpolar
snapshots and end at 2.811846, giving the same reduced increment 1.623691.
The minimum uncertainty-adjusted half-period margins are 0.7291 and 0.6120;
no sampled link is ambiguous under the declared longitudinal sensitivity.
These data do not justify calling the geometric candidates distinct
polarization-winding sectors. The increment is not quantized, and a unique
finite-sample lift still excludes neither unsampled winding nor transverse
electronic components. It is not an absolute spontaneous-P selection.

Native integer h/min/s timing sums to 3155 s, or 28.0444 reported DFT
core-hours at 32 ranks; the fourteen Slurm allocations sum to 4029 s, or
35.8133 core-hours including launch, copying and audit costs. The raw
fourteen SCF and forty-two NSCF logs/band tables replay identically on HF
and locally. Electronic settings and source geometries remain unchanged.

### C.3. A residual pass does not resolve the sampled maximum

After the preserving band's ordinary pass, its cached physical forces and
stress give positive-to-negative energy derivatives along the original
fractional-coordinate/cell segments 2-to-3 and 5-to-6. Three registered
fractions, 0.25, 0.50 and 0.75, are calculated on each segment. All six
SCFs converge with genuine 32-rank MPI and the original six electronic
input hashes; the two sides are calculated rather than imposed as mirrors.

The new highest sample is 38.508448 meV/f.u., 5.702425 meV/f.u. above the
original nine-image sampled maximum of 32.806023. Mirror-side maxima differ
by only 0.0000174 meV/f.u.; that agreement is not a total barrier uncertainty.
An inserted fifteen-image cached band still replays an ordinary maximum
vector residual of 0.099372458 eV/Angstrom, without optimization or extra
SCFs. Thus an unchanged residual pass can coexist with a materially different
sampled maximum. The preserving curve in Table A1 and Figure A2 remains its
original nine-image observation, not a sampling-resolved activation barrier.

![Actual static sampling checks on the preserving reconstruction](figures/hfo2_G1_sampling_bridge_20261009/hfo2_G1_sampling_bridge.png)

Figure A4. Six same-contract static checks of the preserving step-69 path.
(a) Nine original samples and six added energies on the original normalized
generalized arc. (b) The two tested segments, each with two original endpoints
and three new static fractions. Straight connections are sample guides;
the dotted line is the old sampled maximum. These are neither a smooth
MEP fit nor a stationary-TS certificate. The source contains fifteen records,
not fifteen new SCFs or independent replicates. Editable exports, all raw
SCF audits and source data accompany the [case](../../benchmarks/hfo2_channels/20261008/preserving_sampling_bridge_20261009/README.md).

The six SCF evaluations consume 542.98 seconds, or 4.8265 core-hours at
32 ranks. Allocation wall time is 601 seconds, or 5.3422 core-hours including
launch and audit overhead. The checks validate a sampling concern chosen
from already seen gradients; they are not the prospective material holdout.
Their higher energies refer to the specified linear reconstruction, not
an independently relaxed continuous path or a certified barrier error bound.
Stationarity, transverse stability and matched-boundary predictions remain
separate tests. The ordinary target, electronic inputs and no-CI policy are
unchanged.

The same bounded screen was then applied to the other two already completed
edges. Five additional static SCFs checked two near-peak fractions of T-to-PO
and three fractions of PO-to-M segment3-to-4. The T-to-PO and PO-to-M sources,
and the preserving flip, use an identical ordered periodic PO initial well;
their shared energy reference is not an alignment of different endpoint wells.
The M endpoint has a maximum atomic force of 0.011905 eV/Angstrom and maximum
stress component of 1.915095 kbar, within the registered endpoint screen.

| Source path | Old peak | New static | Increase | Residual |
|---|---:|---:|---:|---:|
| T-to-PO | 115.210164 | 115.221524 | 0.011360 | 0.060037 |
| PO-to-M | 71.582207 | 82.747071 | 11.164864 | 0.123408 |

Table A2. Five real reconstruction SCFs, not stationary transition states or
an error-bounded continuous MEP. Energies are meV/HfO2 relative to common PO;
residuals are eV/Angstrom. Each reconstructed band contains 12 total/10 moving
images; only T-to-PO passes the ordinary 0.10 residual target. Force replay
uses no optimization or further SCFs. The respective
forward/reverse sampled barriers are 33.900331/115.221524 for T-to-PO and
82.747071/154.131671 for PO-to-M. Their differences equal the unchanged
endpoint energy differences, -81.321193 and -71.384600 meV/HfO2.
([Actual five-point audit](../../benchmarks/hfo2_channels/20261008/G1_peak_sampling_20261009/README.md))

Job 28380672 completed normally in 544 s at 32 CPUs. The five actual SCFs
consumed 481.761 s, or 4.282322 core-hours; allocation cost was 4.835556 core-hours.
All five retained the original six physical-input hashes, converged SCFs and
complete E/F/stress. The PO-to-M screen estimated 81.718849 meV/HfO2;
the actual sampled maximum of 82.747071 is a direct check, not independent
validation of a barrier predictor. The increased residual required the bounded
ordinary refinement of the twelve-image PO-to-M band described below.
The prior nine-image source is preserved. No electronic retuning, CI or new
mechanistic channel is introduced; neither a higher static sample nor a
denser plot substitutes for transverse relaxation or stationarity.

The required M refinement subsequently passed at 0.0969640 eV/Angstrom
after one new FIRE iteration from all twelve cached initial SCFs. Its sampled
maximum is 82.679957 meV/f.u. and its fixed endpoints are unchanged. This is
a continuation of the same edge and must not be called an acceleration
benchmark: the source representation and force history differ from the
earlier nine-image trajectory. Allocation cost is 9.884444 core-hours.

The reversing terminal band was checked at two registered near-peak
fractions, 3-to-4 at 0.99 and 4-to-5 at 0.01. Both same-contract SCFs yield
392.822011 meV/f.u., below the existing 392.822905 maximum. The screen's
tiny Hermite-predicted increase is therefore not an observed DFT increase
or evidence of sub-microelectronvolt total precision. An eleven-image
cached reconstruction replays 0.0940228 eV/Angstrom, without another
optimization or SCF. These two measurements cost 1.309117 SCF core-hours
and 1.813333 allocation core-hours; the finite peak-static budget has used
thirteen of fourteen points, with no automatic use of the spare point.
The final sampled values and their sampling/refinement histories are
retained separately rather than pooled into stationary barrier estimates.

### C.4. Reference choice changes the apparent compactness of the path

In the unweighted parent-scaffold chart, the three geometric T patterns
capture 95.54% of the formation band's squared displacement at its sampled
maximum, but only 31.43% at PO. Local compactness does not establish a
global few-mode plane. Complete T-Gamma subspaces are evaluated separately
in their mass metric; their displacement weights are not energy fractions.

The stronger published Cmma reference is registered without using modal
coverage to choose its frame. Four tied registrations are retained. In the
common T-origin mass metric, its four-direction span captures 27.18-27.19%
at the formation maximum and 52.24% at PO, versus 89.69% and 30.42% for
the three T patterns. The different reference and metric explain why those
fractions must not be mixed with the preceding scaffold values. The lower
or higher coverage of one chosen span is not a test of a nonlinear fixed-parent
model. Published LDA directions are used as geometric candidates, not as
our PBE Hessian or energy model.
([Registration and scalar evidence](../../benchmarks/hfo2_channels/20261008/cmma_path_mapping/README.md))

### C.5. A measured mixed response exposes a conditional-model limitation

Eight existing same-centre probes at a nonstationary historical image supply
a restricted atomic-strain quadratic block. Releasing only its measured
atomic direction softens the retained strain curvature by 12.36%/12.43%
at the two amplitudes. Four previously seen longer-axis checks give a
maximum energy residual of 0.04169 meV/cell; this is retrospective consistency,
not an independent barrier forecast or a total uncertainty estimate.

The predicted atomic offset, 0.02165 Angstrom, places the specified release
line outside the probe convex hull. The mixed response is measured, but a
conditional DFT branch is unvalidated. No stationary T Hessian is spliced
into this different centre, and stability of unmeasured directions is not
assumed. This failed coverage test identifies what the local model cannot
yet predict, rather than being concealed by a smooth contour.
([Same-centre data and limits](../../benchmarks/hfo2_channels/20261008/restricted_quadratic/README.md))

## Appendix D. Matched endpoints and complete reference tables

### D.1. Endpoint identities, screens and seen well responses

The qualified G1 candidate-connectivity review permits the registered
mechanical comparison without declaring the pilot maxima certified saddles.
The first G2 endpoint calculations fix the two substrate vectors to the
relaxed T reference at nominal biaxial strain zero and release atomic
coordinates and the third lattice vector, including tilt, at zero applied
pressure. Zero strain therefore describes the T-referenced substrate, not
free-cell conditions or zero strain relative to each phase's own lattice.
All screened endpoints use the original ABACUS/PBE electronic contract.

The PO+ four-step canary retains $Pca2_1$ but reaches its step cap before
the physical convergence screen. An exact-geometry continuation reuses only
its terminal E/F/stress, with a fresh BFGS Hessian, and converges after two
new SCFs. The registered T seed passes after its initial SCF, without an
optimization step. Raw outputs reproduce both endpoint energies, forces
and stress; the substrate vectors match. Symmetry labels agree at all three
registered tolerances, without imposing symmetry on either calculation.

| Endpoint | Maximum atomic force | Open traction norm | Relative energy |
|---|---:|---:|---:|
| PO+ | 0.022285 | 1.928435 | 0.000000 |
| T | 0.000500 | 0.009595 | 12.359975 |
| M | 0.018207 | 0.843275 | -92.540494 |
| PO- preserving seed | 0.014888 | 1.692852 | -0.035259 |

Table A3. Screened endpoints at the common T-referenced, partially clamped
zero-strain substrate. Forces are eV/Angstrom, open traction is kbar and
energy is meV/HfO2 relative to the PO+ endpoint in this same ensemble.
The endpoint targets are 0.03 eV/Angstrom and 2 kbar on open components;
reaction stress in clamped components is not a failed convergence criterion.
Neither the screen nor the phase label certifies a Hessian minimum.
The PO- row names its registered seed; polarity and pattern registration
remain separate from a space-group assignment.
([Actual PO+ continuation](../../benchmarks/hfo2_channels/20261008/clamped_PO_continuation_E046_20261009/README.md);
[actual T endpoint and finite matrix](../../benchmarks/hfo2_channels/20261008/clamped_endpoint_matrix_20261009/README.md))

The T-minus-PO+ well separation decreases from 81.321193 meV/HfO2 in the
free-cell ensemble to 12.359975 meV/HfO2 under this substrate. The PO+
clamping energy cost is 68.961218 meV/HfO2, while the screened T reference
is unchanged. This is an observed boundary-dependent well separation,
not a smooth strain derivative between two points in one ensemble.
In particular it determines neither the switching barrier nor the escape
barrier: their maxima and possible intermediate branches still need matched
path calculations. The M-derived endpoint also passes the same physical
screen after seventeen BFGS steps and retains $P2_1/c$ at the three tolerances,
without restoring its free-cell geometry or imposing its symmetry.
Its energy is 92.540494 meV/HfO2 below PO+ in this same ensemble, which does
not determine the intervening escape maximum. The preserving opposite-polarity
seed converges after seven BFGS steps and retains $Pca2_1$. Its observed
-0.035259 meV/HfO2 offset is displayed without enforcing energy equality or
interpreting it as an intrinsic polar bias or a numerical uncertainty bound.
Electronic-polarity and ordered-pattern checks remain separate requirements
before assigning the matched switching channels. The five registered
zero-strain endpoint representations have now completed their physical screens.

The remaining reversing-minus endpoint independently passes the same screen
after seven BFGS steps and eight SCFs. It does not introduce another physical
phase. Composing the previously registered parent operations fixes a fractional
translation (0.5,0.5,0) in the research frame and a species-preserving atom
permutation. The two relaxed minus endpoints agree under this declared
operation to $1.85\times10^{-11}$ Angstrom, but differ by up to 1.36150
Angstrom at the original ordered indices under periodic boundaries. Their
raw energies differ by $7.28\times10^{-12}$ eV/cell; permuted forces and
stresses covary to $3.0\times10^{-10}$ eV/Angstrom and
$2.38\times10^{-12}$ eV/Angstrom$^3$, respectively. These are observed
pair residuals, not a universal numerical precision estimate.
The dominant geometric pattern has opposite signs in the fixed parent chart
(-0.924790 and +0.924790 Angstrom). That sign is reference-origin dependent,
not an intrinsic phase or electronic-polarization identifier. Endpoint
translation equivalence does not establish whole-path equivalence: the
original ordered correspondences and periodic lifts must be retained.
The strict ordered cache correctly does not identify the two endpoints as
the same input geometry. No production image is remapped by this
[read-only audit](../../benchmarks/hfo2_channels/20261008/clamped_endpoint_matrix_20261009/zero_strain_variant_audit.json).

The two fresh PO+ continuation evaluations use 236.309 seconds, and the
single T evaluation uses 125.001 seconds. The eight preserving-seed and
eighteen M evaluations use 1960.598 and 2554.590 seconds, respectively.
The eight reversing-seed evaluations use a further 1244.435 seconds.
Together these thirty-seven fresh evaluations cost 54.408286 SCF core-hours
at 32 ranks. The preceding five-SCF canary cost is recorded separately.
Initial cached-result reuse avoids an identical SCF; it is not a benchmark
of optimizer acceleration or evidence of the proposed predictive gain.

The first two registered +1% training endpoints also pass the same screens,
each after six BFGS steps and seven fresh SCFs. PO+ retains $Pca2_1$;
the T-derived endpoint instead has Ccce symmetry in the analyzed setting.
This label is already present in its affinely strained seed before DFT:
the substrate includes the original long axis and one short axis, so this
partially clamped strain breaks the equivalence of the two original short
axes. All seven evaluated geometries retain Ccce, and the atomic distortion
remains in the original single geometric T-pattern direction, with amplitude
0.831575 to 0.885753 Angstrom and off-pattern residual below
$1.33\times10^{-14}$ Angstrom. It is therefore recorded as a strained-T
descendant, without enforcing tetragonal symmetry or certifying a bulk phase
minimum. Its well-energy separation from the same-substrate PO+ is
14.145883 meV/HfO2. These are seen training observations, not independent
barrier predictions. The registered +1% M endpoint also passes after sixteen
BFGS steps and seventeen SCFs, retaining $P2_1/c$; its force is 0.007755
eV/Angstrom and open traction is 1.452557 kbar. Its energy lies 114.746973
meV/HfO2 below PO+ in this same ensemble. This deeper well does not specify
the escape barrier: the absolute path maxima, not endpoint depths alone,
determine the proposed selectivity test. The [finite matrix and raw records](../../benchmarks/hfo2_channels/20261008/clamped_endpoint_matrix_20261009/README.md)
retain this distinction.

The two registered +1% minus endpoints have now independently passed the
same screens, each after six BFGS steps and seven fresh SCFs. Both retain
$Pca2_1$ at all three tolerances, with maximum atomic force 0.013275
eV/Angstrom and open traction 0.995807 kbar. Their fourteen evaluations
cost a combined 14.250463 SCF core-hours at 32 ranks. These complete the
ten endpoint representations, not ten distinct phases or certified local
minima. Agreement between symmetry-related endpoint energies is not a
numerical error bar for an intervening barrier.

The seen training endpoints already quantify a necessary well-response
control. From zero to 1% substrate strain, PO+ rises by 6.806502 meV/HfO2,
whereas the screened M representation falls by 15.399977 meV/HfO2.
The M-minus-PO+ well separation therefore changes by -22.206479 meV/HfO2.
This is not a measured change in the M escape barrier. The two preregistered
B1 assumptions make this distinction explicit: fixing the absolute
bottleneck energy gives the same -6.806502 meV/HfO2 barrier-response
increment for every channel; allowing it to follow the final well gives
the channel-dependent increments in Table 2. No anchor barrier, unseen
condition or new mechanism is inferred from these increments.

A calculator-independent work contraction evaluates the local imposed-plane
partial as $V\sigma:[H_\varepsilon^T H^{-T}]$, with ASE row cells, tensile
stress and $H_\varepsilon^{0:2}=H^{0:2}/(1+\varepsilon)$ while the open
third vector is held fixed for this partial. It is not an exact relaxed-branch
envelope derivative at the nonzero screened force and open-traction residuals.
Endpoint stress-work trapezoids differ from the observed well-energy changes
by 0.019--0.347 meV/HfO2. Including atomic and released-cell residual work
along the ordered endpoint configuration chord gives defects of
-0.001--0.204 meV/HfO2. These two-endpoint diagnostics mix unmeasured
curvature and finite-relaxation/numerical effects; they are not independent
sampling-error bounds or proof of a smooth stationary branch. The
[raw-work analysis](../../benchmarks/hfo2_channels/20261008/endpoint_strain_work_E059_20261010/README.md)
retains the signed defects, exact sources and all null increments. Actual
path-maxima response must be compared with these well-only controls before
attributing an improvement to a new mode mechanism or prediction method.

### D.2. Complete sampled forward/reverse values

Table A4. Eight ordinary-converged training bands. Energies are meV/HfO2
under the same P=0 electronic contract; reverse peaks and well differences
retain the endpoint-inclusive convention. These are not certified TS values.

| Candidate | Forward discrete peak | Reverse discrete peak | Final-minus-initial well energy |
|---|---:|---:|---:|
| PO-to-M, common T-plane clamping | 98.458122 | 190.998616 | -92.540494 |
| PO-to-M, +1% substrate strain | 97.227719 | 211.974692 | -114.746973 |
| PO-to-T, common T-plane clamping | 84.416624 | 72.056649 | 12.359975 |
| Reversing-pattern switching, same clamping | 338.968965 | 339.004224 | -0.035259 |
| PO-to-T, +1% substrate strain | 96.763650 | 82.617767 | 14.145883 |
| Preserving-pattern switching, same clamping | 96.828183 | 96.863442 | -0.035259 |
| Preserving-pattern switching, +1% substrate strain | 125.424431 | 125.424431 | -0.000000 |
| Reversing-pattern switching, +1% substrate strain | 294.944634 | 294.944634 | 0.000000 |

### D.3. Complete atomic reference comparison

Table A5. Squared atomic-displacement and atomic-gradient coverage (%)
at each sampled peak in the original T chart. The denominators differ;
cell response remains separate. All four Cmma frames 103/106/114/127 are
retained, and their ranges are not uncertainty intervals or model errors.

| Candidate | Atomic reference space | Rank | Peak displacement | Peak gradient |
|---|---|---:|---:|---:|
| PO-to-T | T geometric triplet | 3 | 99.0957 | 14.3013 |
| PO-to-T | Two lowest complete T optical doublets | 4 | 0.2411 | 20.3963 |
| PO-to-T | Registered Cmma span | 4 | 16.6635–16.6636 | 54.0718–54.0765 |
| PO-to-T, +1% | T geometric triplet | 3 | 99.0045 | 69.7245 |
| PO-to-T, +1% | Two lowest complete T optical doublets | 4 | 0.2822 | 2.0694 |
| PO-to-T, +1% | Registered Cmma span | 4 | 18.2383–18.2386 | 10.5750–10.5766 |
| PO-to-M | T geometric triplet | 3 | 57.9237 | 9.0282 |
| PO-to-M | Two lowest complete T optical doublets | 4 | 12.0962 | 32.7096 |
| PO-to-M | Registered Cmma span | 4 | 36.4746–36.4865 | 93.5085–93.5087 |
| PO-to-M, +1% | T geometric triplet | 3 | 54.9251 | 3.9008 |
| PO-to-M, +1% | Two lowest complete T optical doublets | 4 | 12.9756 | 31.7793 |
| PO-to-M, +1% | Registered Cmma span | 4 | 38.4024–38.4141 | 91.5873–91.5896 |
| Reversing-pattern switching | T geometric triplet | 3 | 26.2915 | 0.0320 |
| Reversing-pattern switching | Two lowest complete T optical doublets | 4 | 16.2196 | 10.0215 |
| Reversing-pattern switching | Registered Cmma span | 4 | 46.5637–46.5639 | 65.1046 |
| Preserving-pattern switching | T geometric triplet | 3 | 82.6301 | 12.9583 |
| Preserving-pattern switching | Two lowest complete T optical doublets | 4 | <0.0001 | <0.0001 |
| Preserving-pattern switching | Registered Cmma span | 4 | 14.8336–14.8337 | 52.7203 |
| Preserving-pattern switching, +1% | T geometric triplet | 3 | 52.1751 | 32.8790 |
| Preserving-pattern switching, +1% | Two lowest complete T optical doublets | 4 | <0.0001 | <0.0001 |
| Preserving-pattern switching, +1% | Registered Cmma span | 4 | 16.9707–16.9708 | 49.4667 |
| Reversing-pattern switching, +1% | T geometric triplet | 3 | 16.1410 | <0.0001 |
| Reversing-pattern switching, +1% | Two lowest complete T optical doublets | 4 | 18.1850 | 22.4841 |
| Reversing-pattern switching, +1% | Registered Cmma span | 4 | 52.8497–52.8499 | 63.5220 |

The ordered cyclic G2 frame change is inverted on geometry-only copies
before transport to the original T chart. Native E/F/stress and ordinary
residuals are replayed in their production frame. Dual gradients use
$g_q=-F[H_T^{-1}H]^T/\sqrt m$, with the same translation-free gauge for the
gradient and each declared reference span. The source translation-overlap
norm 0.0036854 is retained; no ASR repair is made to modes or raw forces.
Published LDA Cmma vectors are geometric controls, not our PBE curvature.

### D.4. Complete local joint-curvature measurements

Table A6. Smallest internal eigenvalues of the symmetrized finite-step
matrices, in eV/Angstrom squared. Steps are in Angstrom. The last column
counts negative eigenvalues after excluding the internal path tangent.
Counts are raw signs, not uncertainty-resolved stationary indices. The
prescribed substrate control and three rigid translations are excluded
from these internal spectra. T and preserving denote the original
sampled zero-strain peak centers, not certified TSs.

| Center | Step | First | Second | Third | Internal negative | Tangent-complement negative |
|---|---:|---:|---:|---:|---:|---:|
| T | 0.01 | -0.276193 | -0.164648 | 0.053216 | 2 | 1 |
| T | 0.02 | -0.273596 | -0.161974 | 0.050132 | 2 | 1 |
| Preserving | 0.01 | -1.271819 | -0.447450 | 0.173653 | 2 | 1 |
| Preserving | 0.02 | -1.281821 | -0.455584 | 0.172180 | 2 | 1 |

Full native responses, gradients, reciprocity and two-step diagnostics are
retained in the [complete local-curvature case](../../benchmarks/hfo2_channels/20261008/G3_complete_curvature_E159_20261014/README.md).
No ASR or source-force adjustment is applied. A positive numerical
complement at a nonstationary center is not post-relaxation stability.

### D.5. Strong references in the joint released space

Table A7. Minimum released eigenvalue in eV/Angstrom squared, followed by
its raw negative count in parentheses. Columns distinguish both centers
and finite steps. All four Cmma registrations are retained, including
numerically tied rows. Each reference is an atomic span in the common
36-dimensional internal problem; cell directions remain released. These
unequal retained ranks are not matched-size competing predictors. The
complete T optical reference is a strong control, not an excluded failure.

| Atomic reference | Retained rank | Released dimension | T, 0.01 | T, 0.02 | Preserving, 0.01 | Preserving, 0.02 |
|---|---:|---:|---:|---:|---:|---:|
| T geometric triplet | 3 | 33 | -0.230562 (2) | -0.227205 (2) | -1.271818 (2) | -1.281820 (2) |
| T complete optical doublets | 4 | 32 | -0.175388 (1) | -0.173389 (1) | -0.246013 (1) | -0.253469 (1) |
| All T optical | 33 | 3 | 2.956214 (0) | 2.956557 (0) | 3.090798 (0) | 3.092172 (0) |
| Cmma 103 | 4 | 32 | -0.153434 (1) | -0.151224 (1) | -0.447452 (2) | -0.455596 (2) |
| Cmma 106 | 4 | 32 | -0.153437 (1) | -0.151227 (1) | -0.447451 (2) | -0.455594 (2) |
| Cmma 114 | 4 | 32 | -0.153437 (1) | -0.151227 (1) | -0.447451 (2) | -0.455594 (2) |
| Cmma 127 | 4 | 32 | -0.153434 (1) | -0.151224 (1) | -0.447452 (2) | -0.455596 (2) |

The [reference-complement records](../../benchmarks/hfo2_channels/20261008/G3_complete_curvature_E159_20261014/reference_complements_complete.json)
retain separate atomic-only and joint release tests. This table uses the
joint blocks; published reference frequencies are not substituted for
PBE curvatures. Raw signs and registration spread do not bound physical
error. The all-optical retained Schur block has two raw negative directions
in each case despite its positive released block.

### D.6. Bounded fixed-coordinate relaxation comparison

Table A8. Actual new SCFs, final released-gradient norms (eV/Angstrom)
and native call seconds at fixed retained coordinates. A dash denotes
no gradient-converged final candidate.

| Center | Initialization | New SCFs | Final released norm | Native call seconds |
|---|---|---:|---:|---:|
| T | Identity, budget exhausted | 12 | - | 1489.02 |
| T | Curvature seed | 2 | 0.02630 | 244.81 |
| Preserving | Identity, budget exhausted | 12 | - | 1118.94 |
| Preserving | Curvature seed | 5 | 0.01959 | 452.40 |

Each arm has the same twelve-attempt limit and 0.03-eV/Angstrom
released-gradient norm target. The last evaluated trial is not substituted
for a missing final candidate. Native seconds exclude cached-center
evaluation and the preceding 296 curvature SCFs; earlier path costs are
also additional. Identity costs are censored, so ratios of these call
counts are not measured speedup factors.

The [native canary records](../../benchmarks/hfo2_channels/20261008/G3_conditional_canary_E160_20261014/README.md)
retain all 31 new SCFs, attempt reservations, allocation times and
cache-only optimizer replays. Final candidate energies are
-9782.63681659 eV/cell (T) and -9782.63600314 eV/cell (preserving), using
four-formula-unit cells. No unseen-strain label, new full path or
post-relaxation curvature is inferred from this comparison.
