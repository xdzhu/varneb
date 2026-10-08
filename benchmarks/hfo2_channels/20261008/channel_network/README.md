# Four-channel observations with a common PO+ initial state

Frozen evidence on 2026-10-08, not final competing barriers. There are 37
already evaluated image records and zero new DFT calls in this analysis.
All use ABACUS/PBE/LCAO100Ry/full10auDZP, Hf4O8/4fu, original2x2x2Gamma,
P=0/E=0, ordinaryNEB0.10 and spring0.2. No input, optimizer, running source,
geometry or atom mapping is changed. The two candidate labels refer to
registered parent operations, not certified different MEPs or irrep labels.

| Path viewed from PO+ | Source job / step | Source residual (eV/A) | Provisional discrete maximum (meV/fu) |
|---|---|---:|---:|
| PO→T, reverse of existing T→PO chain | 28274895 / 10 | 0.284995243 | 116.276980 |
| PO→M | 28275259 / 10 | 0.293325305 | 95.199138 |
| PO+→PO−, T-pattern preserving candidate | 28288045 / 1 | 2.074969790 | 145.779475 |
| PO+→PO−, T-pattern reversing candidate | 28288063 / 1 | 1.230216579 | 435.168982 |

The T→PO **forward** maximum is34.955787meV/fu, not the116.276980PO→T
value appropriate for this network. Raw initial PO+ energies agree within
1e-8eV/cell and ordered periodic geometries within1e-9A. Candidate endpoints
may have different integer representatives; no atom permutation is performed
to obtain this comparison. Directional barrier identity is retained.

Every row is unconverged and lacks measured sampling/error bounds. Thus
`ready_for_bounded_discrete_comparison=false`; H1 is not evaluated. The
displayed values neither rank optimized MEPs nor prove that an unexamined
decay channel is absent. Initial high maxima can fall as a path relaxes.
Both flipping source allocations were live while exporting their complete
step1 snapshots; this is a read-only observation, not a live-run fork or
restart. Earlier completed source allocations are continued through R14.

## A useful representation diagnostic, not a final mechanism

At the sampled peak, the same rotated-T geometric triplet captures94.92%
of T→PO parent-relative squared displacement,52.84% forPO→M,54.29% for the
preserving flip candidate, and approximately zero for the reversing candidate.
This is a displacement norm fraction, not an energy contribution. Complete
T-Gamma reconstruction is also not proof of local harmonic validity.

In the two step1 flip peaks, atomic/cell maximum extended-force vectors are
0.590196/2.074970 and0.192700/1.230217eV/A respectively; spring maxima are
below1.1e-16eV/A. Their present residuals are dominated by the generalized
cell block, unlike the two earlier atomic-dominated decay-path observations.
Do not freeze the cell or change electronic parameters to suppress these
forces. Relax the existing paths and then select joint mode/strain directions
from actual bottlenecks. No Hessian or TS certification is claimed here.

## Evidence and replay

`network_specification.json` declares four required candidates and their
source directions. `analysis.json` records every source hash, profile,
residual, mode amplitude, strain and blocker. Both flip observations and
the terminal gap observation include POSCARs, frozen SinglePoint trajectories,
and full raw numeric E/F/stress. PO→M reuses the previously published
`../chain_observations/PO_M_lifted_step10` rather than duplicating it.
Raw ABACUS sources remain on HF with checksums in each observation report.

```bash
python -m scripts.analyze_hfo2_channel_network \
  --specification benchmarks/hfo2_channels/20261008/channel_network/network_specification.json \
  --output /a/new/path/analysis.json
```

The generic library contract and bounded selectivity comparison are documented
in `docs/CHANNEL_COMPETITION.md`. Errors must have measured provenance; this
pilot does not invent them from a force residual. Case/algorithm and clean
archive/HF replay validation records accompany each tested milestone.

Final clean archiveb777ad49/70acd20c:845passed,2skipped,78.69s;54focused
tests pass. HF ASE3.23.1b1 replays the identical complete archive: maximum
force/barrier differences0, parent-patternQ difference3.68e-16A, GammaQ
difference1.78e-15sqrt(amu)A, Green strain difference0. Source and module
hashes agree. Seevalidation.json/hf_replay_check.json. Only these validation
receipts and documentation follow the tested tree, not code/numeric changes.
Live state is separately recorded inlive_jobs_2032CST.json; it does not
replace the frozenstep1 evidence used in this observation study.
