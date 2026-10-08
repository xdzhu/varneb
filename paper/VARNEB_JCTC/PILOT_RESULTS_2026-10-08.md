# Audited pilot results — not final JCTC conclusions

All new DFT in this record uses the existing ABACUS/PBE100Ry/full10auDZP/
Gamma2x2x2 contract. No cutoff, SCF tolerance, pseudopotential, orbital,
smearing or electronic k mesh has been changed. Source archive and Slurm
handles are in `../../benchmarks/hfo2_channels/20261008/run_registry.json`.

## Parent gauge and switching candidates

T, PO and the literature M seed retain SG137, SG29 and SG14 over the declared
symprec 1e-4/1e-3/1e-2 Angstrom grid. Applying two explicit inversion operations
with fixed fluorite-parent permutations produces SG29 switching candidates.
Their geometric Gamma distortions reverse. In a normalized rotated T-pattern
triplet the initial T amplitude is (0,0,0.831575) Angstrom and PO is approximately
(-0.972368,0,0). One candidate keeps the major x pattern; the other reverses it.
These geometric patterns have not been relabelled as phonon eigenvectors or
specific X irreps. A continuous polarization branch along the paths remains
pending; the later endpoint electronic inversion check is reported below.

Both candidates were calculated independently in array28246022. Their total
energies are -9783.249675811965 and -9783.249675811956 eV/cell, with maximum
forces 0.000925 eV/Angstrom and stresses 0.130872 kbar. Numerical equivalence
therefore passes within this calculation contract; it is not a barrier or
Berry-polarization result.

## T-cell atomic Gamma force constants

Sixteen independent SCFs supply two central-displacement families at 0.01
and 0.02 Angstrom, each symmetry-reduced to eight points in the unchanged
twelve-atom cell (1x1x1). All hashes, one genuine DSIZE32, SCF convergence and
complete forces/stress were re-audited after the batch finished. Raw force
constants and separately ASR/permutation-symmetrized matrices are archived.

No negative optical Gamma eigenvalues were found at either amplitude. The
lowest optical frequencies are 1.020298 and 0.953709 THz. The maximum sorted
frequency change is 0.066590 THz and maximum symmetrized force-constant element
change is 0.057053 eV/Angstrom^2. Raw ASR drift is 1.04e-5/1.50e-8 eV/Angstrom^2;
the largest symmetrization adjustment is 0.003838/0.001364 eV/Angstrom^2.
The appreciable relative spread in the lowest mode must accompany any soft-mode
interpretation. This demonstrates stability in the measured fixed-cell
conventional-cell Gamma atomic subspace, not full-q or joint atom-cell stability.
No LO-TO/NAC correction has been applied.

## M preparation and the unresolved low-barrier claim

The literature M geometry gives a same-contract static energy of
-9783.531876337049 eV/cell, but force0.136625 eV/Angstrom and stress2.520136 kbar
require optimization before it is accepted as an endpoint. Its raw energy
must not be used as a relaxed decay-channel reference. A fixed z reflection
and one species assignment register the author's entire PO->M seed to our PO+
gauge (RMS0.007456/max0.012612 Angstrom). The assignment is not changed along
the path. Only M has now been optimized: job28251302 completed in19min16s,
11 BFGS steps, with final force0.011905 eV/Angstrom and stress1.915095 kbar.
All12 retained SCFs and the final energy/geometry/force/stress pass fresh audits.
The relaxed structure is P1 at symprec0.0001 Angstrom, but P2_1/c at0.001 and
0.01 Angstrom. This tolerance sensitivity is reported rather than removed by
post-hoc symmetrization. Its energy is -9783.535214211037 eV/cell, or
71.384600 meV/f.u. below the common PO+ well. This is a well-energy difference,
not a decay barrier or lifetime. The cached T/PO endpoints are retained.

The seven-total-image guided historical chain has a large gap between
pattern amplitudes (-0.06759,0,0.84249) and (-0.79388,0,0.05429) Angstrom.
Its discrete0.963925meV/f.u. maximum is therefore not treated as a precision
benchmark. All three independent samples completed in array28251303 and passed
the same-input/full-SCF audit. Their energies above the same T reference are
35.795717, 59.386224 and 34.045090 meV/f.u. at interpolation fractions
0.25, 0.50 and 0.75. Thus the straight reconstruction of the sparse historical
chain misses an internal high-energy region. The sub-meV discrete maximum
must not support a physical low-barrier or acceleration claim. These samples
do not determine the globally optimized MEP barrier: intermediate geometry
relaxation is a distinct calculation. The raw historical data are preserved.

## Fixed-cell atomic reduction: independently tested local prediction

Projecting the two measured Gamma matrices onto the rotated T-pattern triplet
and releasing30 orthogonal, nontranslational atomic directions predicts
x/y curvatures of1.97930/1.97680 eV/Angstrom^2, compared with frozen values
4.81315/4.81195. The predicted softening is58.88/58.92 percent; the z curvature
is unchanged within this reduction. The eliminated-block minimum eigenvalues
are0.10084/0.08814 eV/Angstrom^2, above the observed two-step operator spread
0.06834, but the small margin and condition numbers490/561 require caution.
This spread is an operational gate, not a rigorous total numerical-error bound.

This is an application of the established Schur complement, not a new theorem,
full cell release, barrier model or TS result. Eight independently evaluated
geometries in array28257780 now compare frozen and linearly responded Qx at
both signs of0.05 and0.10 Angstrom. The0.01-Angstrom training-matrix prediction
and10-percent assessment criterion were committed before DFT; no holdout was
used to refit the response. All8 raw SCFs pass the unchanged physical contract
and all4 paired tests pass. Curvatures are in eV/Angstrom^2:

| Direction | Q amplitude, Angstrom | Fixed prediction | DFT energy curvature | DFT force curvature |
|---|---:|---:|---:|---:|
| Frozen x | 0.05 | 4.813147 | 4.854014 | 4.820535 |
| Frozen x | 0.10 | 4.813147 | 4.836618 | 4.846827 |
| Responded x | 0.05 | 1.979301 | 2.022935 | 1.983617 |
| Responded x | 0.10 | 1.979301 | 2.000460 | 2.004147 |

The largest prediction error is2.2045 percent. This supports the *local atomic
harmonic response* prediction in the tested neighborhood, not a global model
of hafnia. Orthogonal force residuals are retained in the point records: a
linearly responded structure is not a fully minimized conditional surface.
Channel ordering, joint mode-strain reduction and a boundary-dependent
barrier prediction remain untested and are required for the proposed JCTC thesis.

Primary computational workflow references: [Phonopy settings and displacement
definitions](https://phonopy.github.io/phonopy/setting-tags.html) and
[ABACUS–Phonopy interface](https://abacus-rtd.readthedocs.io/en/v3.5.1/advanced/interface/phonopy.html).
Our Python-force route uses Angstrom coordinates/eV-Angstrom forces directly,
not the native ABACUS-interface Bohr conversion convention. The public M/path
source is [Ma and Liu, PRL130096801](https://doi.org/10.1103/PhysRevLett.130.096801),
with the exact repository commit and geometry checksums in the registration.

## Electronic inversion check of the two switching endpoints

All three output-enabled SCFs reproduce their original endpoint energies,
forces and stress exactly at the recorded precision, without changing the
100-Ry/full10-au-DZP/2x2x2 physical SCF contract. The native Berry batches are
complete (28267393_0 and28267783_1/2), and each of the nine retained NSCF
records has been freshly re-audited against its raw inputs, unchanged charge,
runtime version, completed32-rank log and sampled eigenvalue table.

At the2x2x8 NSCF quadrature, the modern-SI R3 values are+0.715814325 C/m^2 for
the designated PO+ well and-0.715814325 C/m^2 for each inversion product.
The reported spin-paired period is1.204918090 C/m^2; the physical eR3/V quantum
is0.602459056 C/m^2. The inversion residual is zero at the native seven-decimal
print precision. All three224-to228 modular differences are0.000526896 C/m^2,
below the predeclared0.01-C/m^2 threshold. The PO+ sampled gap at228 is
4.5690863 eV; all sampled-gap gates pass.

These are branch-valued endpoint components, **not** absolute spontaneous
polarizations or a continuous switching-polarization change. Both inversion
products remain possible ordered switching endpoints; electronic reversal
does not establish two topologically inequivalent MEPs. The nine-total-image
starting bands retain the parent-defined atom identities and contain seven
interior images. Their barriers and mechanical-boundary response remain to
be measured before any claim of switching/decay selectivity.

The first two wrapper attempts failed in post-calculation audits (missing
NSCF eigenvalue output, then an overly literal doubled-quantum unit check),
not in ABACUS SCF/Berry evaluation. The source-version protocol was corrected
and completed outputs reused. The entire endpoint study required three SCFs
and ten NSCF executions, including one duplicate222 NSCF; no physical input
was retuned. These engineering corrections are not a scientific novelty claim.

## Pilot periodic-coordinate audit

Inspection of the switching seeds exposed a wrapped final-endpoint jump.
The audit was then extended to the two active pilots,28257778/28257779.
Their supplied coordinates also contained5--10-Angstrom jumps that plotting's
implicit unwrapping hid; the corresponding physical short segments were
approximately0.3--0.6Angstrom. The optimizer directly used the unwrapped input
metric, so the defective lifts contaminated tangents and spring forces.
Both pilots were explicitly cancelled after preserving complete snapshots and
raw SCFs. This decision was based on a verified geometry-contract defect,
not a single force rebound. Neither pilot's old residual is convergence or
acceleration evidence. Separate, input-stage lift registration and a rejecting
production guard were implemented; recovery must use identical ordered
periodic geometries, unchanged physical inputs and a fresh optimizer state.
The eight independent curvature tests and electronic endpoint tests do not
use these path tangents and remain valid.

On the identical nineteen cached E/F/stress records, correcting only the
registered integer lattice lifts changes the gap-repair residual from
1.386789 to0.829145eV/Angstrom and PO--M from0.990800 to0.881649. All energies
are unchanged. Public trajectories, integer shifts and raw E/F/stress permit
this attribution to be replayed without the private cluster. The lower
residual is not an accelerator result: it compares erroneous and corrected
input coordinates, not two valid optimizers. Fresh-state ordinary-NEB
allocations28274895/28275259 continue the corrected paths with unchanged
physical settings; final barriers are still pending.

## Complete-SCF path observations and a representation failure boundary

Seven complete snapshots of the corrected bands were exported without a
new electronic calculation (66 previously evaluated images). Each image
matches an identical ordered periodic SCF, with unchanged physical input and
raw-log checksums. Numeric-only trajectories reproduce the production force
and energy logs. All observations remain above the ordinary0.10-eV/Angstrom
NEB threshold; their discrete peaks are not certified transition states.

At step6, the T--PO forward discrete barrier is36.560meV/f.u.; referenced
instead to the same PO well, the reverse value is117.881meV/f.u. At step10,
the PO--M forward value is95.199meV/f.u., with a reverse value166.584meV/f.u.
These are provisional values, not a conclusion about the lowest competing
channel. A normal Slurm completion of the ten-step PO--M allocation coexists
with the material termination `max_steps_reached`, not convergence.

The T--PO residual rises from0.215938 at step3 to0.344031 at step6 while
the maximum moves from image3 to image1. At the latter image, atomic/cell
maximum vectors are0.344031/0.030582, versus a0.023657 spring contribution
and0.346935 physical perpendicular force (eV/Angstrom). Thus the observed
rebound is primarily unresolved atomic relaxation, not a spring-dominated
artifact; a subsequent log row falls to0.300872. The chain is not terminated
on this rebound, nor are electronic parameters changed to suppress it.

In a continuous, fixed-initial-gauge projection, the rotated-T geometric
triplet captures93.91% of the parent-relative displacement squared norm at
the step6 T--PO peak but only52.84% at the step10 PO--M peak. Fixed endpoints
also lie substantially outside this triplet. This rejects a complete common
three-pattern representation of all channels, while preserving its use for
variant tracking. Captured displacement is not an energy contribution.
The stationary-T Gamma full basis reconstructs the atomic part, but this is
not evidence for local harmonic validity or a joint atom/cell saddle index.
These findings guide subsequent basis selection and independent validation;
they do not yet establish a predictive mode--strain mechanism.

## First complete candidate-set observation, not optimized-channel ranking

Four same-contract bands now share the same ordered periodic PO+ initial
structure and energy: the reverse of terminal T--PO step10, PO--M step10,
and two registered flip candidates atstep1. The analysis reuses37 complete
image evaluations with no new DFT calculation. Their provisional discrete
maxima, referenced toPO+, are116.277,95.199,145.779 and435.169meV/f.u.,
respectively. The first T--PO forward maximum is34.956meV/f.u.; it is not
the value to use as a PO decay barrier. All four observations remain above
the ordinary residual threshold and lack measured sampling/error bounds.
Consequently, no final channel ordering, H1 conclusion, or distinction
between two optimized switching MEPs is inferred from this table.

The same T-geometric triplet captures94.92%,52.84%,54.29%, and approximately
zero of the parent-relative squared displacement at those discrete peaks.
The two initial flip peaks have atomic/cell maximum vectors0.590196/2.074970
and0.192700/1.230217eV/Angstrom, with negligible spring forces at their centers.
These actual source-metric residuals are cell-dominated, unlike the earlier
atomic-dominated decay-path observations. This provides a concrete reason
to retain joint cell/atomic freedoms and test path-adapted retained directions,
rather than assuming a complete shared T-triplet plane. It is a descriptive
early-stage basis limitation, not evidence for the final TS modes, energetic
decomposition or predictive model. The source allocations continued normally
at that stage; their terminal health-segment audit follows below.

## Terminal flip-segment observations and a moving optimization bottleneck

Both first ten-step switching allocations completed with a step-limit
termination, not ordinary convergence. Every one of the 18 final-snapshot
images was matched to an exact completed SCF and re-audited for geometry,
fixed input, and full energy/forces/stress. The ordinary residuals replay as
0.980656799 and0.700278856eV/Angstrom. No new DFT was used for this audit.

Their common-PO+ discrete maxima are now97.201 and410.744meV/f.u., compared
with145.779 and435.169 at step1. The preserving candidate differs from the
earlier PO--M observation by only2.002meV/f.u.; none of these unconverged
curves establishes a reliable optimized ordering, let alone noninferiority
under a changed mechanical boundary. Sampling and measured error gates
remain absent, and the selectivity result is deliberately unset.

Both energy maxima (image4) remain cell-dominated: atomic/cell vector maxima
are0.397045/0.980657 and0.143499/0.691671eV/Angstrom. Yet the reversing
candidate's largest *chain* residual has moved to atomic motion at image2.
Energy maxima and optimization bottlenecks need not coincide. Basis and
accelerator selection must therefore use the actual local atomic--cell
response rather than impose one early residual label on every image.

The T-pattern triplet captures57.39% and approximately zero of the peak
parent-relative squared displacement. These are representation diagnostics,
not an energetic decomposition or local phonon/saddle result. Complete
mode--strain curvature, branch stability and independent channel-response
prediction are still required before drawing the proposed JCTC conclusion.

Fresh-state ordinary continuations28319570/28319571 were queued after both
current decay-path allocations28298794/28300425, preserving two active
study chains. The actual immutable production factory reused all18 initial
SCFs with its DFT entry disabled during verification. Only new moved interior
geometries will require SCFs in production. Input parameters, nine-total/seven-
interior image counts, the0.10 criterion, and ordinary non-CI policy remain
unchanged. Geometry continuation is not an acceleration benchmark or a full
FIRE-state restart. Frozen evidence and submission details are in
`../../benchmarks/hfo2_channels/20261008/switching_continuation/` and
`../../benchmarks/hfo2_channels/20261008/submission_handles_r17.json`.

## Measured mixed response and failure of a closed two-direction slice

An exploratory re-analysis of the eight completed G0 probes retains their
full42-coordinate gradient differences. The fixed original directions are an
atomic chain secant and scaled symmetric xx+yy strain; neither is an identified
phonon, and this free-cell chart is not the G2 fixed-substrate ensemble. No new
DFT calculation or electronic parameter change is involved. At0.01 Angstrom,
the raw projected matrix, in eV/Angstrom squared, is

```text
 4.671228720   -3.318399815
-3.318981707   19.070358556
```

The mixed entries agree closely at0.02 Angstrom(-3.317902953/-3.319421332),
while raw reciprocity defects of4.08e-5 and1.07e-4 remain visible. The full
Hessian-action operator changes by0.117872 eV/Angstrom squared across the
two amplitudes. This finite-step spread is not a rigorous total DFT error bar.

The two positive projected eigenvalues(3.94316,19.79843 eV/Angstrom squared)
coexist with transverse action norms5.86712 and7.09425. Thus positive curvature
in the displayed slice does not establish stability of omitted directions.
Their self-curvatures are unknown; no Schur elimination of that complement is
performed. The center gradient norm0.118695 eV/Angstrom is nonzero, independently
excluding a stationarity claim. These measurements support retaining full
atom--cell response in subsequent basis selection, not a full-space saddle,
conditional surface, independent barrier prediction or new mode-coupling theorem.
Fresh raw-source verification is recorded with the corresponding delivery.

The actual HF re-audit verifies all eight probe geometries, calculator-input
bytes, SCF convergence at32MPI, full energy/force/stress and original raw-log
hashes. Rebuilt full Hessian actions agree exactly with the archive replay.
This reuses the initial measured data and is not a new independent validation
of any barrier prediction.

## First ordinary G1 T--PO band closure

The sampling-repair continuation28300425 reaches the ordinary0.10 eV/Angstrom
criterion at step6, with an exactly replayed maximum vector0.0598816. The
allocation ends normally by force threshold at23:13:32 CST, not a step cap.
All ten final images match complete same-contract raw SCFs, and the original
ordered T/PO endpoint geometries are unchanged. The continuation uses48 new
interior SCFs, zero new endpoint SCFs, and46.10 allocation core-hours. These
costs exclude the earlier seed/health stages and are not an acceleration claim.

The discrete T→PO maximum is33.889 meV/f.u.; from the shared PO+ initial state,
the reverse PO→T maximum is115.210 meV/f.u. Their difference equals the endpoint
energy change-81.321 meV/f.u. The maximum is image3, whereas the largest
remaining residual is atomic motion at image2. This supplies the first
ordinary-residual-passed G1 edge, not a full-variable TS certificate or a
sampling/error-converged absolute barrier. No CI task or parameter retuning
follows this pass. The other candidate channels and G2/G3 remain incomplete.

To use the freed capacity, only the existing preserving-flip job28319570 has
its capacity dependency removed; it is confirmed running onnode11 at23:20.
Reversing flip28319571 still waits for PO--M28298794. Thus the two-active-chain
limit remains intact, without a duplicate submission or production-source
change. Frozen evidence is in
`../../benchmarks/hfo2_channels/20261008/converged_gap/`; the original
submission journal and separate scheduling update both remain available.

## A converged band does not define a complete low-dimensional mode plane

![Audited ordinary T--PO band](figures/hfo2_T_PO_ordinary_20261008/hfo2_T_PO_ordinary_modes.png)

**Figure | Reference-mode coverage along the ordinary T--PO band.**
(a) Discrete energy relative to T at zero external pressure, with ten audited
SCF images and straight connections only. (b) Geometric amplitudes in the
three rotated-T-pattern directions. (c) Diagonal Green strains relative to
the original ordered T cell. (d) Squared displacement norm captured by the
complete pattern triplet in the unweighted parent-scaffold metric.
(e) Optical mass-weighted displacement fractions relative to T, grouped into
complete reference-frequency subspaces. The5.22THz doublet and9.69THz singlet
are selected by their weights at the sampled maximum; the complement includes
all other optical modes. The fraction at T is undefined, not zero. (f) Atomic
and scaled-cell ordinary-NEB residuals at the eight moving images; fixed
endpoints have no NEB residual. The dashed line is the ordinary0.10eV/Angstrom
target. The curves are one deterministic chain, not replicate statistics;
no numerical-uncertainty bars have yet been measured. The reference modes are
not local path phonons, mode populations or energy contributions, and the
sampled maximum is not labelled a certified TS. Source data accompany the figure.

The triplet captures95.54% of the squared displacement norm at image3 but
only31.43% at the PO endpoint. Thus a locally compact geometric description
near the maximum cannot be promoted to a complete global path representation.
In a different metric and reference, the two displayed T-Gamma subspaces account
for62.16% and28.53% of the image3 optical squared displacement norm; the two
lowest-frequency optical doublets together account for only3.54%. These are
descriptive weights at a finite displacement, not evidence that the5.22THz
reference doublet becomes the TS unstable direction. All T optical reference
frequencies are positive; full atomic--cell stationarity and stability of the
actual bottleneck still require separate measurements.

This limitation is specific to the measured T reference and selected subspaces;
it is not evidence that all fixed parent representations fail. A mapped Cmma
control is now required by the [v2 addendum](../../docs/HFO2_PREDICTION_PROTOCOL_V2_2026-10-09.md)
before attributing independent predictive benefit to path-adaptive reduction.
At the time of the initial figure, no numerical result for that control was
available. The dated E023 descriptive registration below extends this record;
no prediction advantage has yet been established.

The largest remaining ordinary residual is at image2 rather than the energy
maximum at image3. This distinguishes the optimization bottleneck from the
sampled energetic bottleneck when selecting future local probes. No additional
SCFs, changed Hamiltonian, eliminated unstable block or prospective prediction
are introduced by this figure. The [contract](FIGURE_CONTRACT_HFO2_T_PO.md)
and [reproduction bundle](figures/hfo2_T_PO_ordinary_20261008/README.md) retain
the reference, source hashes, complete-subspace rule and explicit limits.

## E023: a stronger parent is channel-dependent, not automatically sufficient

The author's T orientation and the common production PO+ geometry register
the published Cmma directions without choosing frames by their modal coverage.
Of192 proper-frame/origin candidates,16 are T-compatible and four remain tied
at the closest PO+ registration. All four are retained; reference rows, not
production images, are permuted. A common T displacement chart transports
the reference directions explicitly, and a measured-error real subspace
representation retains four directions without ASR correction or new eigenpair
claims. The small source translation admixture is reported rather than repaired.

In the common T-origin mass metric, the Cmma four-direction span captures
27.18--27.19% of the image3 squared displacement norm and52.24% at PO,
compared with89.69% and30.42% for the rank-three rotated-T triplet. These
percentages deliberately differ from the earlier fluorite-origin geometric
fractions: a reference origin and metric change must not be concealed.
The Cmma-origin linear plane retains endpoint residuals6.43/5.60sqrt(amu)A;
neither endpoint is thereby contained in a globally sufficient four-mode plane.

The same registration applied to nine earlier unconverged snapshots shows
that the relative coverage depends on the candidate channel. In particular,
the recorded reversing-flip maximum has approximately44.8% Cmma coverage
versus29.2% for the T-triplet, unlike the formation-path maximum. These
snapshots are not final MEPs, and coverage is not energy partition or a
barrier predictor. This supports keeping a stronger fixed reference in the
future comparison; it does not establish universal superiority of either
reference, path adaptation or joint-cell reduction. Nonlinear fixed-reference
relaxation and larger spaces remain viable controls. Source/scalar evidence,
all equivalent frames and limits are in the
[reproduction bundle](../../benchmarks/hfo2_channels/20261008/cmma_path_mapping/README.md).

## E024: measured atom–strain softening, with a release-coverage failure

The eight original G0 probes supply a same-reference2x2 mixed block at the
nonstationary guided-chain image02. Combining it with the stationary T atomic
Hessian would mix geometries and is not done. With all other variables fixed,
releasing only the one measured atomic chain-secant direction lowers the model
strain-coordinate curvature by12.36%/12.43% at the two step sizes. Its predicted
atomic offset0.02165Angstrom and reference energy lowering1.095meV/cell arise
from the nonzero orthogonal gradient, not from a stationary phonon instability.

The short-step matrix reproduces four known longer-axis energy changes to
0.04169meV/cell maximum residual, with projected/full-gradient-action errors
0.004427/0.005040eV/Angstrom. This is a retrospective local quadratic check,
not a blind holdout, a claim of total energy uncertainty, or a barrier forecast.
Two-step consistency and projected positive curvature do not measure the
unselected modes' self-curvature.

More importantly, the predicted release line lies wholly outside the convex
hull of the measured axes over the requested cell-coordinate interval. It
cannot be promoted to a validated conditional DFT branch without actual
released-point checks. The measured coupling is retained, while a conditional
surface, full TS index and improved independent prediction remain unestablished.
The [pilot and unit/sampling contract](../../benchmarks/hfo2_channels/20261008/restricted_quadratic/README.md)
record this success/failure boundary with no new SCF or parameter change.

## E025: distinct nonpolar snapshots in the common-initial-state network

![Provisional common-PO mechanisms](figures/hfo2_G1_provisional_20261009/hfo2_G1_provisional_mechanisms.png)

Two newly frozen complete observations, PO→Mstep20 and the preserving flip
step14, retain their original E/forces/stress and calculation hashes. Together
with the ordinary-converged T--PO band and the retained reversing step10,
they supply37 image records with no new DFT calls. At the shared PO+ energy
reference, the sampled maxima are115.210,82.136,65.279 and410.744meV/f.u.,
respectively. Only the first passes ordinary NEB convergence; the other three
values cannot support a final ranking or a mechanically selective window.

The preserving central sampled maximum is Pbcn at all three declared symmetry
tolerances,0.001/0.01/0.05Angstrom. Its scaled-cell residual remains0.243948eV/A;
it is not certified as a TS or a relaxed intermediate. All its discrete sampled
energies lie below the relaxed T energy,81.321meV/f.u. above PO+, but intervening
unsampled energies have not been bounded. The reversing centre is Pbca across
the same tolerance sweep despite essentially zero rotated-T-triplet amplitude.
Thus neither coordinate zeros nor low-dimensional projection alone identify
a phase. The PO→M maximum's P1/P2_1 tolerance dependence is retained, not tuned
to a preferred label. Original atom order and continuous lifts are unchanged.

Pbcn-mediated switching is already known: Behara–Van der Ven report a stable
Pbcn intermediate below T for a selected variant pair, and a change to a
T intermediate when the full cell is fixed. Their SCAN result is not a
same-parameter benchmark for this PBE snapshot. A shared group symbol does
not establish identical variants, energetics or intermediate stability;
the phase label itself is not our novelty claim.
([Primary accepted manuscript, Sec.III.D](https://link.aps.org/accepted/10.1103/PhysRevMaterials.6.054403))

The figure pairs energies with ordinary residuals, the registered geometric
shuffle and T-referenced Green strain. It exposes the current representation
and optimization limits rather than hiding them with a smooth contour.
The finite next step remains completion of G1, then matched-boundary G2 and
same-centre joint response/independent predictions with strong fixed-parent
controls. No new material, cutoff, G0 sampling or early CI is introduced.
The [frozen evidence and reproduction contract](../../benchmarks/hfo2_channels/20261008/network_update_20261009/README.md)
and [caption/figure QA](figures/hfo2_G1_provisional_20261009/README.md) distinguish
these descriptive observations from the uncompleted JCTC hypothesis test.
