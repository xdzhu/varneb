# E268: prospective +0.5% sampled-path response test

Registered 2026-10-18 after E267 and before generating target structures.
The author has delegated this design decision to the researcher. This is a
separate sampled-peak protocol, not a retroactive relaxation of v1/v2's
error-certified selector. The original selector and its refusal remain intact.

## Question and finite outcome

Can the already selected models predict the response of two ordinary,
fixed-resolution paths to the intermediate mechanical condition? Complete
this test once, report accuracy, coverage and all preparation costs, then
close the research scope and organize the four approved Results sections.
No new training set, stationary-point/Hessian campaign or material is added.

The target is the endpoint-inclusive maximum of a converged **nine-image**
ordinary band, relative to its common PO+ endpoint, in meV/HfO2. It is written
as a sampled peak, not a stationary activation barrier. A .10 eV/A residual
does not supply an energy-error bound. Unresolved physical rankings remain
unresolved; H1's original delta=0 test is not changed.

## Selection and information chronology

1. Pin E148's complete eight-band report. Select the smallest sampled peak at
   zero strain within switching and within decay; ties sort by registered
   channel name. This selects `PO_flip_T_pattern_preserving` and `PO_to_T`.
   These are representatives of the examined sample, not certified global
   minima. Do not replace either after observing the target.
2. Freeze this protocol, B0 values, B1 formulas, model availability and the
   complete reference list before generating any +.005 structure. All G2,
   T8 and revised P8 results are already seen; none is a new blinded test.
3. The only initially visible +.005 information is the three relaxed endpoint
   features PO+, preserving PO- and T: ordered geometry, E/F/stress, identity
   audit and actual cost. They are necessary production endpoints and are
   available equally to every arm. No intermediate target image is evaluated
   or its energy exposed in this stage.
4. Substitute these audited endpoint energies into the predeclared B1
   formulas. Freeze the complete numeric/abstaining forecast panel with source
   hashes before evaluating either target path's interiors. No model fitting,
   reference selection, calibration or domain expansion after endpoint access.
5. Run at most the two registered paths; score both, including capped/failed
   paths and abstentions, without selecting a winning metric after the fact.

The private Cmma source bank is not needed and remains unopened.

## Forecasts and applicability

B0 is the arithmetic midpoint of the corresponding zero/+1% sampled peaks.
With anchor Bhat(0), dI=[E_I(.005)-E_I(0)]/4 and
dF=[E_F(.005)-E_F(0)]/4 (eV/f.u.):

- B1-fixed = Bhat(0) - dI;
- B1-follow = Bhat(0) + dF - dI.

These are zero-fit null models. Preserve signed values; a value below the
endpoint-inclusive lower limit is an abstention, not clipped to zero.
Predict each peak, its change from zero strain, and T-minus-preserving.
No propagated training midpoint is presented as an unseen error bound.

Keep B2 frozen, B3 atomic release, B4 joint release and B5 stability/branch
handling, with the path-adapted representation, T-triplet, two complete T
doublets, all T optical modes and all four Cmma registrations
(103/106/114/127); retain both original curvature steps .01/.02 A.
Their current local fits describe the zero-strain internal chart. They have
not supplied a certified continuation of the entire sampled path with
physical epsilon, and the sampled centres fail the stationary-envelope
prerequisites. Record their physical-path forecasts as unavailable with
these specific reasons. The revised 109-parameter P8 model also has no
validated physical-epsilon transfer. Do not substitute its local-Q accuracy
for a prediction of mechanical response, or count an abstention as accuracy.
The original preserving 2D negative-block refusal remains separately visible.
No representation is declared inferior merely because this transfer is absent.

## Unchanged physics and bounded execution

Hf4O8,4 f.u.; original ABACUS/PBE/100 Ry/full 10-au DZP, original Gamma2^3,
six byte hashes and binary; P=0,E=0. T-substrate in-plane rows scaled by
1.005; the third row including tilt and every atom remain free.

Three fresh endpoint starts use the corresponding audited zero-strain
endpoint with only the imposed plane affinely changed. Preserve atom order
and continuous periodic lift. Reuse existing ClampedPlaneFilter BFGS,
fmax=.03 eV/A, open-traction norm<2 kbar, maxstep=.02 A, at most80 updates
(81 SCFs) per endpoint. No old energy cache is reused after strain changes.
Audit phase identity at the same three tolerances and the registered
geometric mode signs. A changed phase or failed endpoint is a reported
failure, not symmetry-restored or force-frozen into the desired result.

For each band, start from the corresponding converged zero-strain G2 ordered
polygon. Affinely impose the target plane and distribute the two endpoint
corrections linearly along the image index. Do not permute atoms, relabel a
different pathway or optimize a seed against target interior energies.
Use nine total images (seven interior), two exact endpoint caches, FIRE,
k=.2 eV/A2, fmax=.10 eV/A, maxstep=.02 A, maximum cell step=.02,
minimum pair distance1.6 A, maximum deformation .25, four candidate-step
retries, P=0, no CI. At most80 updates per path (567 fresh interior SCFs),
24h per allocation. A single force rebound is not a stopping rule.
No extra peak points or stationary refinements are part of this first test.

Only hf/hfacnormal01, each allocation one node/32 CPUs and each native SCF
mpirun -np32, OMP/MKL/OpenBLAS=1. At most two committed allocations including
pending/held. Endpoint allocations: PO+ then preserving PO- in shard0;
T in shard1, each at most8h. Exact original handles are reconciled after
transport errors; no duplicate submissions or automatic budget resets.
Fresh immutable source/bundle, native input preflight, held resource check,
native E/F/stress audit, every attempted SCF and actual core hours retained.

## Scoring and closure

Report per-path signed/absolute error, maximum error across the two paths,
mean absolute error, response error and pair-difference error. Separate
zero-feature B0 from endpoint-feature B1 in costs; give every model access
to the same endpoint inputs. Show both B1 arms and every unavailable arm.
Retain null physical error bounds; no significance or unique TS claim.
If B0 is as good, or all proposed response models abstain, explicitly withdraw
the physical predictive-advantage claim. The software, local validation and
material response remain results, with journal positioning based on that
evidence. This outcome closes this test; it does not trigger a new campaign.
