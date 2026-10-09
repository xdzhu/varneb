# G1 reversing-channel peak sampling (E043)

Preregistered 2026-10-09 after job 28319571 completed normally. Source:
`../terminal_G1_update_20261009_17/reversing_step45`, observation SHA256
`23c1f7e3799ab0f1b388bd2791418c400f26eccf5679a3353c05dff326ef3f08`.
The nine-image ordinary band replays 0.0940228456663 eV/Å and has a sampled
common-PO maximum of 392.822905197136 meV/HfO2. The highest image is not a
certified stationary transition state.

## Finite experiment

The force-and-stress physical work derivative changes sign in segments 3→4
and 4→5. Hermite screening suggests maxima near λ=.9922934 and .0077066,
respectively, at 392.829955636375 meV/HfO2. These are screening estimates,
not DFT results or independent predictions. Calculate exactly two frozen
linear-reconstruction points: (3→4,.99) and (4→5,.01).

Retain the existing unwrapped fractional lift, ordered Hf4O8, linear cell,
P=0, common PO energy and all six original electronic/basis input hashes.
No relaxation, MIC remapping, endpoint recalculation, CI, new material,
new independent channel or automatic band restart. One hfacnormal01
allocation, one node, 32 true MPI ranks via mpirun, single-thread libraries,
one-hour cap; the two SCFs run serially. With the previous eleven peak
statics, this uses thirteen of the fourteen-SCF G1 sampling cap; the
remaining allowance is not automatically spent.

After completion audit the actual rank probe, charge convergence, complete
energy/forces/stress, six physical input bytes on HF, raw log/call/geometry
hashes and reconstructed eleven-image ordinary residual. If the latter
fails, record it before deciding on any separately bounded refinement.
Two statics do not bound the continuous-path error or certify a TS.

Full G1 still requires variant, mechanism and polarization-branch audits.
In particular the lower Pbcn midpoint of the preserving path must not be
quietly excluded from the nonpolar-escape discussion. No G2/canary or
holdout DFT is submitted by this experiment.

## Actual result

Job **28410062** COMPLETED0:0, 18:08:45–18:12:09 CST, node61, 32 CPUs,
204 s / 1.813333 allocation core-hours. Both static SCFs passed real
32-rank, electronic convergence, complete E/F/stress, original six-byte
input, geometry and raw-log/call audits. New energies relative PO are
392.822010684540 and 392.822010697728 meV/HfO2; neither exceeds the old
392.822905197136 peak. The tiny Hermite-predicted increase is **not** an
observed DFT increase and is not promoted to a numerical precision claim.

The inserted eleven-image cached band replays 0.0940228456663 eV/Å, so no
additional band optimization or SCF is needed for this residual check.
SCF times 75.1612/72.1145 s, total147.2757 s /1.309117 core-hours; this is
distinct from allocation cost. G1 peak statics now13/14, no automatic use
of the spare point. Continuous-MEP error and TS certification remain open.

`analysis_local.json` independently replays the HF report. `completed_HF/`
contains exact prepared geometry and raw logs/inputs/results, with repeated
UPF/orbitals and unused charge-density restart excluded from export (the
originals remain on HF). Full six-byte checks were actually performed on HF.

The separate `terminal_network_local.json` reuses all40 complete ordinary
images: T/PO10, refinedM12, preserving9 and reversing9. Three fixed symmetry
tolerances (.001/.01/.05 Å, angle1°) give Pa-3 for the reversing center and
Pbcn for the preserving center. The latter is−14.838201 meV/HfO2 below PO;
its stability is still unmeasured. Zero T-triplet projection at the former
does not identify it as T or a specific local phonon mode.

`raw_gap_audit_HF/` contains40 copied original `istate.info` files plus
source/log hashes and the exact read-only HF executor. Local replay checks
every byte hash and48-occupied-band table. Minimum sampled gaps in eV:
T/PO4.252037, M4.131676, preserving4.589282, reversing4.017321. These are
Gamma2³ sampled gaps, **not** full-BZ insulation or Berry branch evidence.
All40 raw SCFs' E/F/stress, ordered geometry and six physical bytes were
freshly rechecked on HF; this analysis adds zero DFT.

Reproduce from a repository checkout:

```sh
python -m scripts.prepare_hfo2_reversing_peak_sampling analyze --root benchmarks/hfo2_channels/20261008/reversing_peak_sampling_20261009/completed_HF/prepared --output /new/output/analysis.json
python -m benchmarks.hfo2_channels.20261008.reversing_peak_sampling_20261009.check_terminal_network --specification benchmarks/hfo2_channels/20261008/reversing_peak_sampling_20261009/network_specification.json --output /new/output/network.json
python -m pytest -q tests/test_hfo2_reversing_peak_sampling.py tests/test_hfo2_terminal_network_audit.py
```

## Remaining G1 gate, before G2

Endpoint electronic inversion is already measured; the two full switching
paths' branch-continuous Berry evidence is still missing. The next bounded
experiment should reuse all three audited PO endpoint charges/NSCFs and
evaluate only the14 existing interior geometries (7 per original terminal
band). First check the actual R3 axes; a one-component audit is invalid if
rotation/shear makes it unable to track the polarization component.

At most14 output-only SCFs plus42 fixed-charge NSCF quadratures (2×2×2,
2×2×4,2×2×8 along the already validated R3 recipe), two allocations at most
concurrently,32MPI each. This is a property audit, not fourteen new peak
samples or channels. Preserve baseline physics, add only the already
validated charge/bandgap outputs; reproduce baseline energy/forces/stress
within the existing endpoint tolerances before accepting any Berry value.
NSCF energies must never enter the barrier table. No source is overwritten.

Record the physical and native spin-paired quanta, reduced polarization in
the actual cell and quadrature sensitivity; do not silently halve the native
value, select a smallest-|P| branch or force continuity through ambiguous
links. A sampled insulating gap alone does not close this gate. Unresolved
axes/branch/quadrature require an explicit inconclusive result and reviewed
next action, not automatic additional points or premature G2 submission.

This is the concrete missing measurement, not a new reference-Hessian,
extra endpoint, new material or expansion of the13-channel production cap.
The lower nonpolar configuration remains in the scientific discussion;
G3's existing branch/stability tests, not a symmetry label, must decide its
interpretation. Global escape resistance is not established by M/T alone.

`prepared_HF/` is the retained **failed preflight** generated from the early
CRLF driver (no SCF); only `prepared_HF_LF/` and `completed_HF/` use the
accepted LF runtime. The initial local POSCAR/manifest bytes also require
explicit Git renormalization after their byte-preserving attributes are set;
the delivered snapshot test checks their declared hashes, not just numeric
geometry. This does not change any input, geometry or DFT result.
