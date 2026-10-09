# E039: preserving flip reaches the unchanged ordinary NEB threshold

Frozen complete snapshots exported at 12:55 CST on 2026-10-09, not a final
channel ranking or a certified stationary-TS dataset. All numeric analysis
uses the unchanged source-bound physical audit, with **zero new DFT calls**.
Earlier dated observations and figures remain unmodified.

## Actual source states and common-well comparison

Slurm `sacct -X` at 13:08:27 CST reports preserving job **28319570 COMPLETED
0:0**, ending 12:50:13 after 13:32:57, and reversing **28319571 RUNNING**,
32 CPUs on node26/hfacnormal01. The original preserving terminal summary is
retained, not recreated from a rounded live log. It reports `converged`,
`termination=force_threshold` and 0.13122409286734182 eV/cell sampled barrier.

| Same-PO+ candidate | Frozen step | Replayed fmax (eV/A) | Sampled maximum (meV/HfO2) | Ordinary pass |
|---|---:|---:|---:|---|
| PO to T, reused reverse view | 6 | 0.059881616 | 115.210164 | yes |
| PO to M, reused terminal observation | 39 | 0.097904490 | 71.582207 | yes |
| T-pattern-preserving flip | 69 | 0.099372458 | 32.806023 | yes |
| T-pattern-reversing flip | 32 | 0.200001508 | 394.066333 | no |

Common ordered PO+ energy is -9783.249675811956 eV/Hf4O8 cell, four formula
units. ABACUS/PBE/100 Ry/full 10-au DZP, original INPUT/KPT/UPF/orbital hashes,
Gamma2x2x2, P=0/E=0 and ordinary NEB 0.10 eV/A remain unchanged; no CI.
This export includes 18 records (four cached endpoints), plus 19 reused T--PO/M
records in the network. The 37 records are not 37 new SCFs or replicates.

## New physical interpretation, and what remains unproven

The preserving band now has two sampled peaks at images3/5, approximately
32.806023/32.805958 meV/f.u., both Pca2_1 across registered 0.001/0.01/0.05 A
symmetry tolerances. Its central image4 is Pbcn at all three tolerances and
**-14.838201 meV/f.u.** relative to PO+. Atomic/cell NEB residuals there are
0.062019/0.099372 eV/A; transverse stability is unmeasured. Do not certify a
stable new phase or add a Pbcn endpoint/path because a sampled centre is lower.

At the two peaks the physical generalized tangential forces are
+0.242486/-0.242483 eV/A. They are not stationary TSs even though their
atomic/cell NEB residuals are below the ordinary target. No physical
tangential-force convergence tolerance or CI has been introduced to alter
the ordinary residual criterion. Negative-index and sampling audits remain
separate measurements, not inferred from the connected energy polyline.

The central T-chart F_xx=1.122022 expands a direction within the registered
G2 substrate. The existing matrix can test whether that boundary changes or
eliminates the central branch, not merely shifts a single sampled peak. This
is a prospective interpretation of already seen data, not a frozen independent
prediction or an additional phase/Hessian budget.

The reversing centre at step32 is Pbca/Pa-3/Pa-3 across the same tolerances;
its selected T-pattern triplet is nearly absent. Preserve the sensitivity
sweep rather than assign a desired phase/irrep or a local phonon to a
truncated projection. The dominant whole-chain residual is atomic, image1.
At 13:07:59 a later rounded live row was step33/0.193306; that row does not
replace the complete frozen step32 evidence.

## Reproduce without a calculator

```sh
python -m scripts.analyze_hfo2_chain_observations \
  --root benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255 \
  --variants benchmarks/hfo2_channels/20261008/reference_variants \
  --gamma benchmarks/hfo2_channels/20261008/gamma_analysis/T_d0.01.npz \
  --output /new/path/analysis.json
python -m scripts.analyze_hfo2_network_progress \
  --specification benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/network_specification.json \
  --output /new/path/network.json
```

The exporter checks single DSIZE32, SCF convergence, complete finite E/F/stress,
all six input hashes, ordered periodic geometry, raw logs and complete snapshot
hashes, plus numeric replay against the optimizer row. Both actual HF and local
reports are retained. The original exporter status is deliberately only
`complete_observation_not_final_result`, even for an ordinary-converged band.
See [the new figure](../../../../paper/VARNEB_JCTC/figures/hfo2_G1_preserving_pass_20261009/README.md)
and [delivery receipt](validation_delivery.json). Full G1 phase/polarization/
sampling gates, G2/G3 responses and independent prediction advantage are
still unmeasured; no premature matched-boundary or holdout job is submitted.
