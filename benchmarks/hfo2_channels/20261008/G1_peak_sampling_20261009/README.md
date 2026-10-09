# E041: five bounded peak sampling checks

Preregistered on 2026-10-09, after E040. These are **single static SCFs of
linear reconstruction geometries**, not optimized continuous MEPs, stationary
TS certification, independent predictions, or extra production channels.

Frozen ordinary-converged sources: T→PO job28300425/step6 (10 total images),
PO→M job28298794/step39 (9 total images). Source observation digests, exact
ordered structures and original six-file physical contract are verified.
Common PO energy −9783.249675811956 eV/Hf4O8; four formula units, P=E=0.

Exactly five statics:

- T→PO segment2→3: λ=0.99; segment3→4: λ=0.02.
- PO→M segment3→4: λ=0.25, 0.50, 0.75.

The cached energy/work derivatives suggest a peak around81.718849meV/HfO2
inside the PO→M segment, compared with71.582207 at the old sampled peak.
T→PO screens suggest <0.01meV changes. These numbers select points and are
**not new DFT results or error bounds**. Every interpolation preserves the
existing unwrapped fractional lift, ordered Hf4O8 and cell; no MIC/remapping,
alignment, geometry relaxation or electronic parameter changes.

One hfacnormal01 allocation, one node,32MPI, five serial SCFs,1h cap; at most
one other study chain can be active. No CI, automatic restart or G2 submission.
G1 peak-sampling cap14 SCFs =6 completed preserving +5 this batch +at most3
for the reversing path after its completion. Original13-new-chain limit and
G1→G2→G3 order unchanged. No repeated endpoint evaluation.

Prepare without DFT from a clean repository archive:

```sh
python -m scripts.prepare_hfo2_G1_peak_sampling prepare --repository . --root /fresh/prepared
```

Do not replace existing namespaces or source archives. The runtime refuses
changed source, recipe, physical contract, geometry, repeated point or summary.
Large unchanged UPF/orbital payloads stay on HF; result exports retain their
verified hashes and full raw SCF/input structure records, not duplicate binaries.
The completed result and actual cost must be recorded before drawing conclusions.

## Actual completed checks

Job28380672: COMPLETED0:0,2026-10-09 14:47:13--14:56:17 CST,544 s,
node5,32 CPUs. Real32-rank SCFs, identical six-file physical bytes,
converged charge and complete E/F/stress checked on HF. Result archiveSHA256:
`84321c58fda89b57ba629bd426f5ed8fb9233870dd664529e6f115728f656f3b`.
All raw logs, INPUT/KPT/STRU, geometry, rank probe and recorded input hashes
are in`completed_HF/`; repeated UPF/orbitals remain on HF with verified bytes.

| Source path | New samples (meV/HfO2 relative PO) | New minus old maximum | Inserted ordinary residual |
|---|---|---:|---:|
| T→PO | 115.213825834;115.221524095 | +0.011359934 | 0.060036699,pass |
| PO→M | 82.747070887;82.255037106;68.495859025 | +11.164863720 | 0.123408121,not passed |

Both inserted bands contain12 total/10 moving images. Replays use0 optimizer
steps/extraSCFs. The T→PO forward barrier is33.900331348 and reverse/common-PO
barrier115.221524095. PO→M forward/reverse sampling values are82.747070887 /
154.131670658; endpoint difference−71.384599770 is unchanged. These are
**reconstruction samples**, not relaxed continuousMEP barriers or TS evidence.
In particular the M band cannot be called ordinary-converged after insertion.

Actual SCF sum481.761171 s/4.282322 coreh; allocation544 s/4.835556 coreh.
The Hermite hint81.718849 was lower than the measured sample82.747071;
the difference is not a calibrated barrier error. Offline raw-log, EFS,
recorded six hashes, actual INPUT/KPT bytes and production-writer geometry
checks pass. Local/HF q tails differ≤2.22e-16, cell matrices exactly identical.
Common-well endpoint audit is recorded separately and does not prove stability.

```sh
python -m scripts.analyze_hfo2_G1_peak_sampling --completed benchmarks/hfo2_channels/20261008/G1_peak_sampling_20261009/completed_HF --output /fresh/analysis.json
```

The necessary M refinement preserves all12 input geometries and exact raw SCF
caches, ordinary0.10, FIRE/maxstep0.02/k0.2/P0, fixed endpoints and physical
inputs. Fresh optimizer,20-step/8h cap, at most200 moving-geometry SCFs;
not full optimizer state restoration or acceleration evidence. At most2 live
study calculations. Only submitted after clean-source and actualHF cache/CLI
preflight. No additional scientific edge, new G0 probes, CI or earlyG2.
