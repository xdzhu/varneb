# Explicit G3 budget amendment: post-release transverse curvature

2026-10-15, E181. Reported before preparing native probes or submitting DFT.
The user has continued the bounded HfO2 research goal. Original v1/v2,
E150's296raw-peak SCFs and their results are preserved, not overwritten.
Already seen: eight G2 training bands; the original two raw curvature
centers; E160's two Q=0 seeded candidates; T25/25 and preserving17/25
gradient candidates with eight radius-boundary holes. B3 job28867177 is
pending with zero labels at registration. Physical+.005 and eight independent
face points per face remain sealed. This amendment is exploratory, not blind.

## Why and precisely what changes

The original raw-peak Hessian cannot certify stability after nonlinear
release. A few optimizer secants cannot determine a34-dimensional block.
Measure the entire fixed-Q, fixed-control released block at **only the two
existing E160 curvature-seeded Q=0 geometries**. No center search or choice
based on the new signs is allowed. These are two *additional verification
stencil geometries*, not the original raw geometries: the budget addition is
explicit. It does not authorize a third channel, new reference phonon
matrix, full path, material, orientation, supercell or strain condition.

Use the original40-coordinate chart, internal B40x36 and retained A36x2.
R36x34 spans null(A.T); the literal R, literal native anchor and all probes
are frozen in config.json before DFT. x=x_native+B R(sign*h*e_i),
h=.01,.02Angstrom, i=0..33, sign=(+,-). The center's exact native E/F/stress
is reused, never recomputed.34*2*2=136new SCF attempts per center,
**272additional attempts total**, no automatic retry/refund/rollover.
The original296, canary, Gamma, faces118/172, B3 and G2/preparation costs
remain additional. This is a derivative budget, not12optimizer calls per
surface node; it consumes no additional Q sampling/independent slots.
Training commitments remain T34/81, preserving35/81.

All probes must pass the original .5Angstrom total internal radius,
1.6Angstrom minimum distance, finite positive cell, zero physical strain,
unchanged two Q observables, ordered atoms and substrate before any DFT.
An invalid probe stops preflight: no inward move, axis change or domain repair.
Keep original PBE/100Ry/full10auDZP/Gamma2^3/P0/E0, the exact six files and
binary, full forces and substrate reaction stress. No ASR/force repair.

## Resources and stopping

One sequential two-center job: hf/hfacnormal01,1node/32CPUs,
genuine mpirun -np32, OMP/MKL/OpenBLAS1,12h cap. At most two committed study
allocations including held/pending;28867177 counts. No other account jobs
are cancelled. Admission once: write intent, sbatch --hold, validate exact
allocation, release same handle. An ambiguous intent must be reconciled,
never replayed. Fresh immutable external bundle; E161 runtime stays read-only.
At recent126s/call,272calls suggest about9.5h **after allocation starts**;
queue delay is unknown and this is not a promised completion time.
Failed/partial results remain charged and do not form a complete matrix.

## Decision and limitations

Retain both raw40x34 actions and nonsymmetric34x34 matrices; report reciprocity,
explicit symmetrization change, both-step spectra, their difference, and the
actual base released-gradient residual. No pseudoinverse, clipping, enforced
positive signs, refit or best-step selection. A negative sign is evidence
against the local stable-branch interpretation at that sampling resolution.
Positive raw signs at both steps are useful local evidence, **not a measured
total-error bound, exact stationary minimum, global branch or TS certificate**.
These two points do not certify the rest of either5x5sheet or fill its holes.
Independent material predictions and strong native-chart comparisons remain
required; no JCTC readiness or adaptive accuracy advantage is predeclared.
