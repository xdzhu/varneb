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
The every-step raw SCF and final phase-identity audits remain the acceptance
gate. The cached T/PO endpoints are retained.

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

## Fixed-cell atomic reduction: a prediction awaiting independent tests

Projecting the two measured Gamma matrices onto the rotated T-pattern triplet
and releasing30 orthogonal, nontranslational atomic directions predicts
x/y curvatures of1.97930/1.97680 eV/Angstrom^2, compared with frozen values
4.81315/4.81195. The predicted softening is58.88/58.92 percent; the z curvature
is unchanged within this reduction. The eliminated-block minimum eigenvalues
are0.10084/0.08814 eV/Angstrom^2, above the observed two-step operator spread
0.06834, but the small margin and condition numbers490/561 require caution.
This spread is an operational gate, not a rigorous total numerical-error bound.

This is an application of the established Schur complement, not a new theorem,
independent prediction success, full cell release, barrier model or TS result.
Eight independently evaluated geometries will compare frozen and linearly
responded Qx at +/-0.05 and +/-0.10 Angstrom. Their prediction and evaluation
criteria will be committed before DFT; no holdout will be used to refit the
initial response. A linearly responded structure is not a fully minimized
conditional surface.

Primary computational workflow references: [Phonopy settings and displacement
definitions](https://phonopy.github.io/phonopy/setting-tags.html) and
[ABACUS–Phonopy interface](https://abacus-rtd.readthedocs.io/en/v3.5.1/advanced/interface/phonopy.html).
Our Python-force route uses Angstrom coordinates/eV-Angstrom forces directly,
not the native ABACUS-interface Bohr conversion convention. The public M/path
source is [Ma and Liu, PRL130096801](https://doi.org/10.1103/PhysRevLett.130.096801),
with the exact repository commit and geometry checksums in the registration.
