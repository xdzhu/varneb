# Acceleration table: claim-to-evidence audit

This record is a submission gate for Table 3 of `varneb_CPC.tex`, not a new
algorithm benchmark. The six rows use **total** image counts (including two
fixed endpoints) and the public stopping criterion `fmax = 0.10 eV/Å`. The
archived BTO comparisons begin at the same serial step-9 chain; the HfO₂
comparisons begin at the same archived initial chain. All reported savings
are case-specific, not a universal speedup guarantee.

| Case and policy | First crossing | Force at crossing (eV/Å) | New image evaluations in table | Provenance |
| --- | ---: | ---: | ---: | --- |
| BTO, 9 total, global FIRE | step 28 | 0.096704 | 205 | `benchmarks/convergence/results/bto_fire_baseline.log`; summary JSON; job 27675981 |
| BTO, 9 total, BlockFIRE 0.02/image | step 10 | 0.086473 | 79 | `bto_blockfire_002.log`; summary JSON; job 27729004 |
| BTO, 9 total, BlockFIRE 0.01/image | step 5 | 0.089923 | 44 | `bto_blockfire_001.log`; summary JSON; job 27729028 |
| BTO, 9 total, ImageScaledFIRE | step 8 | 0.078068 | 65 | `bto_fire_scaled_002.log`; summary JSON; job 27729012 |
| HfO₂, 7 total, global FIRE | step 42 | 0.0942 | 217 | `hfo2_fire_baseline.log`; summary JSON; job 27678924 |
| HfO₂, 7 total, staged FIRE | coarse step 23 + refine step 1 | 0.092367 | 134 | `hfo2_fire_scaled.log` to step 23, `hfo2_fire_staged_refine23.log`; summary JSON; job 27729231 |

The first-crossing steps and forces above are directly visible in the archived
optimizer logs; each summary reports the same calculator contract and final
path energies. The evaluation totals in `benchmarks/convergence/autoresearch.jsonl`
are explicitly labeled **estimated**. They reconstruct one initial evaluation
of each fixed endpoint plus one evaluation of each interior image at step 0
and every accepted update. For example, BTO global FIRE gives
`2 + 7 × (28 + 1) = 205`, and the 0.01/image BlockFIRE gives
`2 + 7 × (5 + 1) = 44`. The staged HfO₂ value counts both coarse and refine
run startups, including their endpoint evaluations. These arithmetic totals
must not be described as an independent per-image electronic-job audit until
the raw worker/calculator records are reconciled. In particular, manager
manifest entries with near-zero elapsed time can represent reused ASE results
even when `cache_misses` is populated; summing that field blindly overcounts.
Read-only inspection of the original hf worker manifests found 6, 11, and 9
substantial-duration seven-image evaluation batches for the three BTO
alternatives, matching their step-0-plus-update counts. The HfO₂ baseline
manifest has 52 substantial-duration five-image batches over its full run;
the manuscript stops counting at its first crossing (step 42), rather than
using all 52. The old BTO global-FIRE baseline lacks this manifest, so its
205 remains a trace-derived estimate. Duration is not itself proof of an
electronic launch; a definitive per-image count still needs raw calculator
records and endpoint reconciliation.

The resulting BTO and HfO₂ evaluation reductions are respectively
`(205−44)/205 = 78.5%` and `(217−134)/217 = 38.2%` after rounding. The HfO₂
staged path differs from the baseline barrier by 0.00547 eV per 12-atom cell;
this is a path-equivalence diagnostic, not a statistical uncertainty. The
negative 17-image analytic transfer, the slower matched global-cap BTO control,
the incomplete SplitFIRE transfer, and the 20-image HfO₂ BlockFIRE plateau
remain in `benchmarks/convergence/autoresearch.md` and are not hidden by this
table.

## Outstanding before a submission-grade quantitative claim

1. Preserve the tiny source optimizer logs (currently ignored by the general
   `*.log` rule) alongside the summaries or export an equivalent immutable
   machine-readable table. Bind each to its SHA-256 and starting-chain hash.
2. Reconcile estimated new image evaluations against raw worker and ABACUS
   outputs through the **first threshold crossing**. Separate no-op ASE calls,
   cache hits, retries, and endpoint statics. Then update the manuscript's
   current wording that calls these counts directly measured from manager logs.
3. Check final-path equivalence with a common image-arc coordinate, not only
   the endpoint-referenced barrier. Retain the negative-transfer controls.
4. Only if those checks pass should the abstract keep a numerical `38–79%`
   evaluation-reduction claim. Until then, call these reconstructed counts.
