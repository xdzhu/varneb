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
specific X irreps, and electronic polarization branches are still pending.

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
