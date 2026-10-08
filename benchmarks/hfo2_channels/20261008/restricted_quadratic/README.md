# Same-center nonstationary atom–strain quadratic pilot (E024)

This reuses the eight G0 static calculations of job28240257 at the historical
guided-chain image02, **not** the stationary T structure. No T Gamma matrix,
published LDA curvature or different mechanical reference is merged. The
pressure is zero; the original PBE/100Ry/full10-au DZP/KPT/input hashes remain.
All35 other atomic and five other cell coordinates are clamped, including
unselected rigid translations. Their stability is not inferred or needed for
this *restricted* algebraic release. This is not G2 clamped-substrate training.

## Reproduce

```console
python -m scripts.analyze_hfo2_restricted_quadratic --output /new/path/archive.json
python -m scripts.analyze_hfo2_restricted_quadratic --raw-root /historical/r3/probes --output /new/path/raw.json
```

The raw option reads and rechecks the original center, manifest, all eight
signed structures, raw SCF logs,32MPI and fixed input hashes through the
existing G0 audit. Neither option writes into a calculation directory.
The [local report](local_analysis.json) retains raw projected matrices,
reciprocity defects, full measured Hessian actions and transverse response.
The release uses their **explicitly reported symmetric2x2 projection**, not a
fabricated full42x42 Hessian. The two-step action spread0.117872eV/Angstrom2
is an operational stability-screen floor, not a Gaussian sigma or total DFT
error bound. Reports include all API, analysis and source-log hashes.

## Model and measured result

For E=E0+g.y+y.H.y/2, y=Qq+Rr, a stable measured eliminated block yields

```text
r*(q) = -Hrr^-1 gr - Hrr^-1 Hrq q
geff  = gq - Hqr Hrr^-1 gr
Eeff0-E0 = -gr.T Hrr^-1 gr / 2
Keff  = Hqq - Hqr Hrr^-1 Hrq
```

`vcneb.quadratic_reduction.condition_quadratic_energy` implements this
established harmonic reduction, including the nonstationary offset and energy
lowering missing from a curvature-only treatment. It rejects unstable/unresolved
blocks and overlapping/inconsistent bases, retains omitted variables as fixed,
and returns immutable model arrays. Pure analytic and coordinate-rotation tests
check energy, envelope gradient, selected orthogonal stationarity and no hidden
release. This implementation is not itself a new theory or a TS certificate.

The retained direction is scaled equal xx+yy strain; the eliminated direction
is one atomic chain secant, **not a phonon**. At0.01/0.02Angstrom, the frozen
strain-coordinate curvatures19.0704/18.9556eV/Angstrom2 become
16.7126/16.5986 after this one-direction model release:12.36%/12.43% softening.
The effective zero-coordinate gradients are-0.08172/-0.08170eV/Angstrom,
not the original-0.0098765, because the center is nonstationary in the released
atomic direction. Its predicted atomic offset is0.021648/0.021641Angstrom;
model energy lowering is approximately1.095meV/**12-atom cell**, not per atom
or per HfO2. Four formula units are retained explicitly.

The0.01-step matrix and center energy/gradient reproduce the four previously
computed0.02-axis energy changes with maximum residual0.04169meV/cell;
maximum projected/full-gradient-action errors are0.004427/0.005040eV/Angstrom.
No long-step energy is fitted. These points were already examined in earlier
audits and the stability screen uses both step families, so this is a
**retrospective crosscheck**, not a new blind holdout or barrier prediction.

## Units, sampling and why this is not a released DFT branch

In the registered chart, `z_cell=sqrt(2)*L*epsilon`, L=5.127945044745648Angstrom,
for equal epsilon_xx=epsilon_yy. Gradient converts by sqrt(2)*L and curvature
by its square; eV/Angstrom2 times Angstrom2 gives eV/cell per dimensionless
strain squared. The short-step frozen/restricted-released values are
1002.94/878.94eV/cell in that strain convention. They are local nonstationary
restricted derivatives, **not relaxed elastic moduli or substrate responses**.
The runtime lacks Pint; explicit unit-labelled boundaries and the known
Jacobian are checked without installing packages or changing the environment.
No statistical uncertainty distribution or coverage factor is invented.

The known axes have convex hull `|z_atomic|+|z_cell|<=0.02Angstrom`, **not** a
filled2D square. The shortest L1 distance of the predicted release line over
z_cell in[-0.02,0.02] is0.021648Angstrom: its entire requested segment lies
outside that hull. Even a point inside a convex hull would not establish DFT
branch continuity. No released point has been calculated here. Large measured
transverse action also prevents eliminating the unmeasured complement.

**Decision:** retain the measured coupling and nonstationary model implementation;
reject a conditional-surface/TS/independent-barrier claim from these probes.
Future G3 predictions need same-boundary measurements at the actual selected
bottleneck and independent released-point/branch checks within the existing
finite budget. This pilot does not authorize extra G0 SCFs or skip G1/G2.

The [clean archive report](clean_analysis.json),
[HF raw-source replay](hf_raw_analysis.json) and
[delivery receipt](validation_delivery.json) are separate from the initial
working-tree report. The latter differs in three legacy LF/CRLF metadata hashes
only; all other fields match exactly. This does not relax any raw DFT input
hash.45focused tests pass; clean archive997passed/2skipped in100.64s.
The active production source is not overwritten, and no analysis/test source
changes after that full regression.
