# Acceleration table: claim-to-evidence audit

This record is a submission gate for Table 2 of `varneb_CPC.tex`, not a new
algorithm benchmark. The six rows use **total** image counts (including two
fixed endpoints) and the public stopping criterion `fmax = 0.10 eV/Å`. The
archived BTO comparisons begin at the same serial step-9 chain; the HfO₂
comparisons begin at the same archived initial chain. All reported savings
are case-specific, not a universal speedup guarantee. The launch counts below
were independently checked against Slurm job-step accounting on hf on
2026-09-28; they are no longer only optimizer-trace estimates.

| Case and policy | First crossing | Force at crossing (eV/Å) | ABACUS launches to crossing | Provenance |
| --- | ---: | ---: | ---: | --- |
| BTO, 9 total, global FIRE | step 28 | 0.096704 | 205 | `benchmarks/convergence/results/bto_fire_baseline.log`; summary JSON; job 27675981 |
| BTO, 9 total, BlockFIRE 0.02/image | step 10 | 0.086473 | 79 | `bto_blockfire_002.log`; summary JSON; job 27729004 |
| BTO, 9 total, BlockFIRE 0.01/image | step 5 | 0.089923 | 44 | `bto_blockfire_001.log`; summary JSON; job 27729028 |
| BTO, 9 total, ImageScaledFIRE | step 8 | 0.078068 | 65 | `bto_fire_scaled_002.log`; summary JSON; job 27729012 |
| HfO₂, 7 total, global FIRE | step 42 | 0.0942 | 217 | `hfo2_fire_baseline.log`; summary JSON; job 27678924 |
| HfO₂, 7 total, staged FIRE | coarse step 23 + refine step 1 | 0.092367 | 134 | `hfo2_fire_scaled.log` to step 23, `hfo2_fire_staged_refine23.log`; summary JSON; job 27729231 |

The first-crossing steps and forces above are directly visible in the archived
optimizer logs; each summary reports the same calculator contract and final
path energies. The earlier `benchmarks/convergence/autoresearch.jsonl` values
were explicitly labeled **estimated** because manager manifests alone cannot
distinguish new DFT from reused ASE results: some near-zero-duration records
still populate `cache_misses`. Slurm accounting now supplies an independent
external-launch count. The read-only query was `sacct -j <job ids>
--format=JobID,JobName,State,Elapsed,ExitCode -n -P` on hf. For every job
below, numbered steps are contiguous from zero, named `abacus`, and all are
`COMPLETED` with exit code `0:0`:

| Route | Job ID | All ABACUS launches | Launches through first crossing |
| --- | ---: | ---: | ---: |
| BTO global FIRE | 27675981 | 219 | 205 (step 28) |
| BTO BlockFIRE 0.02 | 27729004 | 79 | 79 (step 10) |
| BTO BlockFIRE 0.01 | 27729028 | 44 | 44 (step 5) |
| BTO ImageScaledFIRE | 27729012 | 65 | 65 (step 8) |
| HfO₂ global FIRE | 27678924 | 262 | 217 (step 42) |
| HfO₂ staged coarse | 27729066 | 157 | 122 (step 23) |
| HfO₂ staged refine | 27729231 | 12 | 12 (step 1) |

The optimizer logs contain exactly 31, 11, 6, 9, 52, 31 and 2 printed
step records for these jobs respectively. Counts equal one initial launch
per fixed endpoint plus one launch per interior image at every printed step:
for example `2 + 7 × (28 + 1) = 205` and `2 + 7 × (5 + 1) = 44`.
The HfO₂ staged first crossing counts **both** starts,
`[2 + 5 × (23 + 1)] + [2 + 5 × (1 + 1)] = 134`.
The chronological boundary is visible independently: Slurm
`27675981.204` ended at 14:25:35, matching the BTO `FIRE: 28` log time;
`27675981.205` began afterward. Likewise `27678924.216` ended before
the HfO₂ `FIRE: 42` line (06:45:09), and `.217` began afterward.
For the staged coarse segment, `27729066.121` ended before its `FIRE: 23`
line (01:45:08), and `.122` began afterward. These checks establish the
**number of successful ABACUS process launches through first crossing**,
not a retrospective per-launch audit of every SCF residual or an independent
physical convergence certificate; threshold forces come from the archived
optimizer logs.

The source `.log` files are already Git-tracked despite the general `*.log`
ignore rule. Their SHA-256 values in table order, with the separate staged
coarse/refine logs last, are:

```
8352e567dcddd91a03d7a3f751c68154f4fe677bf465d97e792b62b75b752672
88275fea44a100740f282b353d7bf64a8988ed8be65139f8b4b0fc5fa107c31e
b5ff8d85e7e80f5da7d9b03f104c3d153481ab3d9f3640feed7846ab0681b0b5
f91dacef95a178e431ea0b81fa96cbcc120db7e144e8e3fca5b25ac942a21f98
fd28f31a2f8ed71f7e09d0d48a7f9e5a1922b1c0e6c5f0d50554021377613887
599a62ddce34fc5ed8e9f798cd3ea5034ebc7db567c97c5ce38ab6b402f4e0a3
10c830cb192d6655eb6ea8dc1a8805f2c5ea1ca935477d151e0d1e62cb744d4d
```

The BTO baseline resume source has 90 frames; its last nine and the
nine-frame alternative-start `initial-vcneb.traj` have exactly equal atomic
species, Cartesian coordinates and cells (maximum difference `0.0 Å`). The
alternative-start file SHA-256 is
`bd1d768bca52f7e29ed83f2d0e44aa49f463a0dfafc8fc5269db186e5d6438c1`.
Both HfO₂ baseline and accelerated variants use the baseline's original
`initial-vcneb.traj` (SHA-256
`a5cff9b6a5d51070fbe72602a2b2072c8ab44b7fb13e5bd716a549b57b6c7f0c`);
the staged refine resumes its own coarse step-23 snapshot (SHA-256
`dcaeb05ba327f0ee30f6745e1146d0b70fe4efdc662b70bcc27dd5259851d889`).

The resulting BTO and HfO₂ ABACUS-launch reductions are respectively
`(205−44)/205 = 78.5%` and `(217−134)/217 = 38.2%` after rounding. The
HfO₂ *final-run summary* barrier difference is 0.00547 eV per 12-atom cell;
it must not be confused with the separately archived first-crossing comparison
below. The
negative 17-image analytic transfer, the slower matched global-cap BTO control,
the incomplete SplitFIRE transfer, and the 20-image HfO₂ BlockFIRE plateau
remain in `benchmarks/convergence/autoresearch.md` and are not hidden by this
table.

## First-crossing whole-chain comparison

Six small evaluated chain snapshots are archived in
`benchmarks/convergence/first_crossing_chains_20260928/`; their SHA-256
values are recorded in the two JSON reports below, and the job/step labels
above identify their original hf runs. `scripts/audit_acceleration_path_equivalence.py` checks ordered
endpoints, species and the common atom–strain metric, then piecewise-linearly
compares each accelerated chain with the baseline at a shared normalized
joint-arc coordinate. The snapshots are at the **first** threshold crossing,
not the later end of the historical baseline job. The raw pressure is zero
for both cases, so the stored energies already equal enthalpies.

| Case and variant | Max. common-arc profile difference (meV/f.u.) | Max. common-arc joint-geometry separation (Å) |
| --- | ---: | ---: |
| BTO BlockFIRE 0.02 | 2.739 | 0.02374 |
| BTO BlockFIRE 0.01 | 1.324 | 0.01260 |
| BTO ImageScaledFIRE | 1.854 | 0.01229 |
| HfO₂ staged FIRE | 4.839 | 0.07946 |

The BTO paths are monotonic, so their identical endpoint-defined maxima do
**not** mean their profiles are identical. For HfO₂, the first-crossing
maximum-relative-enthalpy difference is 1.172 meV/HfO₂, whereas the larger
4.839 meV/HfO₂ discrepancy occurs elsewhere along the aligned profile.
These are deterministic, post-hoc diagnostics and were not a predeclared
path-equivalence acceptance gate. They support a bounded comparison on the
tested chains but not a claim of an identical MEP, saddle basin, or transferable
material-independent speedup. The source-bound numerical reports are
`benchmarks/convergence/bto_first_crossing_common_arc_20260928.json` and
`hfo2_first_crossing_common_arc_20260928.json`; the local audit has two
synthetic tests, including a check that re-spacing the same path is not
misclassified as a new geometry.

## Outstanding before a submission-grade quantitative claim

1. Preserve the above Slurm step accounting and source logs while remote
   accounting retains the records. The number of ABACUS **launches** is now
   independently verified, but their individual raw SCF histories were not
   archived; do not call these retrospectively SCF-audited evaluations.
2. The first-crossing common-arc comparison is complete, but no equivalence
   tolerance was predeclared. Retain the measured discrepancies, negative
   transfers and this limitation; do not recast a post-hoc diagnostic as a
   passed equivalence gate.
3. The abstract may report the measured `38–79%` reduction in **ABACUS
   launches to the shared force threshold on these cases**, but not a universal
   speedup, identical MEP, or matched-saddle guarantee.
