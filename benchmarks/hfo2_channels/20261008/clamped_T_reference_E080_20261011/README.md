# E080: paired T/M escape geometry and native atomic gradients

2026-10-11. Zero-DFT training analysis of the already ordinary-converged
PO-to-T job28692775, terminal step10/.093589766eV/Angstrom, compared with
the existing PO-to-M job28722320, step13/.099214088. Fixed electronic input,
100Ry/full10auDZP, PBE/P0/Gamma2x2x2, mechanical substrate/open-cell chart,
original atom order/lifts, T reference, ranks and all four E023 Cmma frames
are unchanged. The full finite G2/G3/control/prediction objective is incomplete.

## What the actual data establish

| Peak atomic representation | T displacement (%) | T gradient (%) | M gradient (%) |
|---|---:|---:|---:|
| T geometric triplet, rank3 |99.0957|14.3013|9.0282|
| Two complete lowest T optical doublets, rank4 |0.2411|20.3963|32.7096|
| All registered Cmma frames, rank4 |16.6635--16.6636|54.0718--54.0765|93.5085--93.5087|

Fractions have different denominators: translation-free mass-weighted
displacement versus the dual atomic energy covector in the original T chart.
They are neither energy partitions nor model/force accuracy scores. T peak4
has internal gradient norm0.010196606eV/(sqrt(amu)*Angstrom); Cmma omitted
norm0.006909923--0.006910274. Absolute norms, translation diagnostics and all
nine records remain in gradient.json and the figure CSV. M peak is image3.
Sampled peaks are not certified continuous saddles. A near-zero force can
produce a sensitive percentage; a large fraction does not certify accuracy.

Geometrical representation does not by itself establish energetic response,
and gradient coverage is channel dependent. This keeps the Cmma control
serious rather than dismissing it from displacement coverage. It does not
select modes/bottlenecks or show independent predictive superiority. Source
LDA modes remain geometric directions, not new PBE eigenpairs/local TS modes.

## Evidence and reproduction

- displacement.json: original reference analyzer, actual nine portable T SCFs
  reparsed, original fixed inverse G2 orientation on geometry-only copies.
- gradient.json: same unchanged references, proper dual force pullback and
  declared translation-free coordinate restriction; no ASR/source repair.
- Terminal T receipt SHA c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e.
  Native77-SCF/full-six-file/runtime audit is historical E079, not a new
  local full-runtime or licensed-orbital audit. Local replay checks9records.
- M scalar analyses E071/E078 and audit E068 are reused unchanged; do not
  overwrite historical reports. The shared analyzer now admits the actual
  completed T receipt explicitly; unknown receipts still reject before labels.
- Author-source bytes are the same pinned public E023 files, kept outside Git.
  No author executable runs, and no mode arrays are redistributed.

Replay the two original module commands from E071/E078 with --case pointing
to clamped_T_terminal_E079_20261011, the above T receipt SHA, the original
seed/reference/source arguments, and **fresh output paths**. Generate the
paired figure with:

```text
python -m scripts.plot_hfo2_clamped_comparison --output FRESH_FIGURE_DIRECTORY
```

The plotter checks fixed report/content/native-receipt pins, individual native
snapshot hashes, identical cached PO+ SCF bytes/EFS and ordered initial geometry.
Reference filenames/comments differ, not the actual boundary. Initial POSCAR
last-digit serialization differs by4.44e-16Angstrom; it is recorded, not repaired.
Each path has its own original normalized arc, not a shared absolute coordinate.
The zero-displacement T endpoint has undefined coverage: flagged in scalar data,
omitted in rendering, never invented as zero physical coverage. All Cmma frames
remain; no confidence band, curve fitting or MEP interpolation is performed.

The figure contract was written before plotting. The initial local figure
had a legend overlapping panel(d)'s label; it was moved to a recoverable
E:/TEMP proof directory and corrected, not published as final. Structural
tests now check legend/axes/panel-label nonoverlap. validation.json records
actual local and clean-archive tests, native-analysis replay and PDF inspection.

No DFT launch, Slurm write, new chain, parameter/runtime/optimizer change,
G3 selection, +0.005 holdout access, extra Gamma/Hessian or model freeze.

Final verification:63related tests pass both locally and in the actual clean
Git-archive cwd, direct exits0, no errors/failures/skips. Actual clean T
gradient replay is byte-identical. The displacement replay differs only in
its executed-script byte hash (existing CRLF/LF checkout policy); normalized
source text and every other JSON field are identical, recorded explicitly
in validation.json. Units-skill static audit has0findings, a heuristic rather
than a physical uncertainty bound. The18-page paper compiles0; actual pages
15--18 and final PNG are inspected, not a full final-submission audit.
