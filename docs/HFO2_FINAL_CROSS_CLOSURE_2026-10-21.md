# HfO2 final two-cross closure

User approved the finite blueprint on 2026-10-21. Active goal: finish the
two missing strain-pattern crosses, interpret either outcome, integrate the
manuscript/main figure/case, compile and show the reading PDF, then push.
No further material, strain, model or Hessian matrix is authorized.

## Fixed scientific argument and terminology

Question: how does internal atomic accommodation under partial clamping
change homogeneous reversal relative to escape from the same polar well?
Existing whole-path evidence reduces the sampled 0-to-1% preserving penalty
from 28.596 to 8.533 meV/HfO2. Local fixed-pattern controls explain part,
but not all, of that change. The final experiment separates imposed strain
from changing the four retained geometric observables.

- Preserving reversal: preserves the signed reference shuffle, not every
  oxygen atom's coordination identity.
- q0/q1: four original-T geometric observables at the original zero/+1%
  sampled peaks; not the earlier two local joint coordinates Q1/Q2.
- J: joint release of the 29 remaining atomic and three permitted cell
  directions, with three translations gauge-fixed.
- Transverse proxy: mean O-minus-Hf displacement, not electronic polarization
  or an independently identified phonon eigenmode.
- U: local candidate energy above the same-strain initial well, not a
  newly calculated path barrier or a certified saddle.

## Finite execution

Protocol/case: `benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021`.
Two existing minus-sign J diagonal candidates are reused. New conditions:
zero strain with q1 and +1% with q0. Hold all four values exactly.
Electronic inputs unchanged: original 100 Ry/full10auDZP/PBE/six assets/
Gamma2x2x2/zero pressure and field. At most 40 SCFs per cross, 80 total;
two one-node32MPI hfacnormal01 allocations, 2h each. No new A treatment,
second sign, midpoint, whole path or automatic numerical-budget extension.

Source `1a4925da` was pushed before allocation release. Bundle pin:
`ecb498b40ec3f19b593bfca198a754d1f2cdc62deb389fcc1dbb6b0d85ea2683`.
Remote namespace:
`/public/home/iai806/abacus/agent-runs/20261021-varneb-strain-pattern-E277-r1`.
Zero-SCF native preflight passed both diagonal caches and cross geometries.
Original held jobs: 28995326 (zero_q1), 28995327 (plus1_q0).
Both original jobs started at 12:52:27 CST and completed normally at
13:12:30 / 13:01:30. They used 10 / 4 fresh SCFs, 15.52 allocation
core-hours total, and reached released gradients 0.015790 / 0.023149 eV/A.
The 114-file terminal native archive was exported and all fourteen new
E/F/stress results independently replayed locally, together with both
reused diagonal candidates. No additional DFT is needed or authorized.

## Completed physical closure

U(0,q0)=83.211101; U(0,q1)=114.553122;
U(1%,q0)=92.398077; U(1%,q1)=107.767116 meV/HfO2.
Same-strain pattern contrasts are 31.342021 / 15.369039, giving
I=-15.972982 meV/HfO2. Fixed-q strain responses are +9.186976 / -6.786007.
The common initial-well shift cancels in I; the q1 relative-energy decrease
alone does not imply absolute local-energy lowering.

At +1%, the q0 cross lies only 0.699140 above the existing free-chain
sampled peak, with transverse proxy 0.120 A versus the free peak's 0.123 A.
The native q1 J is 16.068179 above that peak, with proxy 0.206 A.
This links the retained pattern to the transverse and open-cell response;
it is not a unique single-mode energy partition or another complete barrier.
The code/case and main Table 5 / Figure 11 now contain these actual results.
The finite calculation is closed. All 53 scope tests passed; the 65-page
reading PDF compiled successfully with the existing LuaLaTeX, without
actual typography/reference warnings. Table 5 (page 23), Figure 11 (24),
Abstract/Discussion/Conclusions and adjacent pages were visually inspected.
The compiled PDF and a concise Chinese reading guide are ready for the
author; no new calculation follows this final integration.

## Paper closure and stop rule

Use research-paper, English/JCTC framing. Evidence-first update order:
Results -> Discussion/Conclusions -> Introduction/Abstract. Keep the four
approved Results sections and the existing evidence tables/figures, with
old negative prediction results identified as completed rather than an
unfinished promise. Do not ask the author to reconfirm the approved scope.

The final cross figure uses the existing Python/matplotlib workflow, source
data only, normal-weight (a) panel labels, top/right axes and large fonts.
Four panels link local energy, pattern contrast, transverse displacement
and open-axis response. No mock data or fitted surface is required.

Accept terminal results whether positive, negative or inconclusive; no
automatic expansion. Then update the concise Chinese reading guide and
reproducible case, run relevant native/constraint/figure/manuscript tests,
compile with the existing LaTeX installation, visually inspect changed
pages and present the PDF. Push the completed evidence and source.
Preserve the five unrelated pre-existing user edits in the working tree.
