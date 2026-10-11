# Competing switching and phase-escape channels in hafnia: a mode-strain pathway analysis

Working manuscript, 2026-10-10. This is a connected main-text draft, not a
submission-ready article. The numerical Results below use explicitly dated
pilot observations. All ten registered endpoint representations across the
two matched-substrate training conditions have passed their physical screens;
G2 channel barriers, local conditional branches and
independent material predictions are still missing. An abstract and a
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
([Fan et al., 2025](https://doi.org/10.1038/s41524-025-01647-w))

Strain-dependent mode-amplitude models have also been compared with DFT under
uniaxial and biaxial loading and connected to hafnia's formation landscape.
([Lee, Lee and Yu, 2026](https://doi.org/10.1038/s41535-025-00841-9))
Predicting a strain-induced distortion is therefore not, on its own, the
methodological distinction sought here. The test concerns the competing
barrier responses after subtraction of the same initial well, with their
actual mechanical freedoms and strong controls kept explicit.

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
Ordinary NEB uses 0.10 eV/Angstrom, without climbing images; cached fixed
endpoints are not recalculated at every iteration. The repaired T-PO band
has ten total images, while new channel bands have nine total, seven moving.
A residual pass does not supply a barrier error or a saddle certificate.
The preserving candidate subsequently receives six same-contract static
sampling checks, producing a fifteen-image cached reconstruction without
an additional optimization or endpoint recalculation (Section 3.3). The M
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

## 3. Present results: audited pilot evidence

### 3.1. The formation band and the shared polar starting well

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

Figure 1. Audited discrete energy, geometric pattern amplitudes, Green strain,
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

Table 1. Four ordinary-residual-passed candidates, not an error-bounded
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

### 3.2. Nonpolar snapshots do not uniquely identify a switching mechanism

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

Figure 2. Complete same-PO+ terminal observations at steps 6/1/69/45, corresponding
to Table 1. (a) Discrete low-energy profiles. (b) The higher, ordinary-converged
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

Figure 3. Native-path electronic property audit. (a) Reported classes in
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

### 3.3. A residual pass does not resolve the sampled maximum

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
sampled maximum. The preserving curve in Table 1 and Figure 2 remains its
original nine-image observation, not a sampling-resolved activation barrier.

![Actual static sampling checks on the preserving reconstruction](figures/hfo2_G1_sampling_bridge_20261009/hfo2_G1_sampling_bridge.png)

Figure 4. Six same-contract static checks of the preserving step-69 path.
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

Table 2. Five real reconstruction SCFs, not stationary transition states or
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

### 3.4. Reference choice changes the apparent compactness of the path

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

### 3.5. A measured mixed response exposes a conditional-model limitation

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

### 3.6. A common substrate narrows the T-to-PO well separation

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

Table 3. Screened endpoints at the common T-referenced, partially clamped
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
the channel-dependent increments in Table 4. No anchor barrier, unseen
condition or new mechanism is inferred from these increments.

Table 4. Observed final-well shifts and endpoint-following B1 null
increments across the two seen training strains. Each quantity is in
meV/HfO2; the switching correspondences remain separately registered paths.

| Final endpoint representation | Observed well-energy shift | B1 final-following barrier increment |
|---|---:|---:|
| T descendant | +8.592411 | +1.785908 |
| M | -15.399977 | -22.206479 |
| PO-, T pattern preserving | +6.841761 | +0.035259 |
| PO-, T pattern reversing | +6.841761 | +0.035259 |

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

The first [matched G2 pair](../../benchmarks/hfo2_channels/20261008/clamped_G2_E053_20261010/README.md)
uses the zero-strain preserving candidate and PO-to-M escape candidate.
Their starting polygons retain the ordered G1 correspondence and explicitly
recorded periodic lift, but release the same permitted atomic and cell
variables at the common substrate. Only the two unchanged clamped wells
reuse hash-pinned raw results; seven internal images require fresh DFT.
The free-cell source energies are not transferred to this new ensemble.
Starting-path construction and endpoint convergence still do not establish
the optimized competing barriers, their sampling errors or the proposed
selectivity response.

The first ten-step segments lower the ordinary residual from 3.541469 to
1.466699 eV/Angstrom for the preserving candidate and from 0.994612 to
0.615276 eV/Angstrom for the M candidate. Neither satisfies the 0.10 target.
All 154 fresh interior SCFs have complete audited E/F/stress under the
original input bytes and genuine 32-rank execution. Their transport-time
cost is 154.451016 allocation-core-hours, compared with 155.822222 actual
scheduler allocation-core-hours. This is a partial pilot cost, excluding
earlier endpoint preparation and subsequent continuation, not the full
cost of a predictor. These unfinished observations cannot supply channel
rankings or test the proposed selective response.

All eight registered training starts now have same-ensemble endpoints and
continuous ordered geometry; the remaining six seeds pass zero-DFT geometry
and exact endpoint-cache preflights. Preparation reuses the screened wells
and transfers only G1 geometry, never its free-cell energies. The additional
starts are not optimized bands or successful material forecasts. The
[pilot audit](../../benchmarks/hfo2_channels/20261008/clamped_G2_E054_20261010/README.md)
and [remaining-start audit](../../benchmarks/hfo2_channels/20261008/clamped_G2_remaining_E056_20261010/README.md)
retain these distinctions before the held-out condition is accessed.

The zero-strain, common-substrate PO-to-M candidate subsequently reaches
the ordinary target after 43 updates across its registered segments.
The final segment terminates at step 13 with a same-boundary replayed
residual of 0.099214088 eV/Angstrom. Table 5 reports its endpoint-inclusive
discrete image peaks, not a continuous saddle-point certification or a
sampling-converged barrier. All 91 fresh SCFs of this final segment pass
the original six-file, native 32-rank and complete E/F/stress audit; their
recorded transport time corresponds to 90.730841 allocation-core-hours.
This cost excludes the earlier pilot, 20-update segment and endpoint
preparation. The [terminal evidence](../../benchmarks/hfo2_channels/20261008/clamped_M_terminal_E068_20261011/README.md)
preserves the nine-image replay and those distinctions. One converged
escape candidate does not establish the full G2 network, its strain response,
the selected bottlenecks, H1 or the independently tested H2 forecasts.

Table 5. First ordinary-converged G2 candidate at zero imposed substrate
strain. Energies are meV/HfO2 under the unchanged P=0 electronic contract.

| Candidate | Forward discrete peak | Reverse discrete peak | Final-minus-initial well energy |
|---|---:|---:|---:|
| PO-to-M, common T-plane clamping | 98.458122 | 190.998616 | -92.540494 |

## 4. Discussion and remaining material tests

The present data establish that a compact local projection, an apparently
nonpolar midpoint and an ordinary residual pass answer different questions.
None alone determines the competition between homogeneous switching and
phase escape. The core proposed result remains the channel-selective
response to one matched mechanical intervention, not a phase-label survey.

The next material result must contain the registered fixed-plane responses
of absolute switching and escape barriers, endpoint energies and branch
identities. A favorable selectivity change without preserved escape resistance
will be reported as the weaker outcome. If an endpoint disappears under
that boundary, it will not be forced into a nominal two-well band. Independent
predictions must subsequently outperform the declared strong controls at
measured resolution; interpolation or a fixed-parent model may perform equally
well. New directions added after reading labels are revised training, not
successful original forecasts. These tests are still incomplete.

All present conclusions concern periodic small-cell, zero-field, zero-temperature
energy pathways under the named electronic contract. They do not establish
domain-wall nucleation, switching rates, retention times, finite-temperature
free energies or functional-independent critical strains. The preceding pilot
observations support the finite investigation, but are not sufficient for
the intended JCTC methodology claim.

## Appendix A. Reproducibility and editorial completion checks

The archived cases provide frozen numeric observations, input/log/geometry
checksums, continuous ordered lifts, analysis commands and figure CSVs.
Raw production logs remain read-only on HF. Earlier incomplete snapshots
remain dated evidence rather than being overwritten by final values.
No external paper figure or unpublished author dataset is copied into this draft.

Before a submission-ready abstract or conclusion is written, this manuscript
requires: completed candidate paths and error/phase audits; matched G2 material
responses; same-boundary conditional-branch and transverse-stability evidence;
prospectively frozen predictions with strong-control comparison and actual
cost accounting; final numeric captions and an openly reproducible release
of the results permitted for redistribution. These are requirements, not
claims that those outcomes have been achieved. Details belong in main-text
appendices rather than a separate supplementary manuscript, as requested.
