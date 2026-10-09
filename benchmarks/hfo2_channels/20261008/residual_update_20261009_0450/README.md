# E028: complete-frame residual attribution, not final barriers

This dated audit freezes three complete observations from the two existing
G1 production jobs. The export reads the live runs, but writes only to a new
diagnostic namespace. It does not restart an optimizer, change a Hamiltonian,
submit a job or perform another SCF. Every image retains its numeric
energy/forces/stress, original input/log/geometry hashes and exact force replay.

The common initial state is the same ordered PO+ Hf4O8 cell. There are four
formula units, P=0 and external E=0. ABACUS/PBE,100 Ry, full10-au DZP, the
original Gamma2x2x2 KPT and UPF/orbital/input hashes remain unchanged. Ordinary
NEB uses0.10 eV/Angstrom without climbing images; endpoints are cached.

## What was measured

All residuals below are **maximum three-vector norms in the original NEB
generalized metric**, not raw Hellmann--Feynman atomic forces or kbar stress.
The cell rows use the unchanged production cell scale. A sampled energy
maximum is not a stationary transition state or a sampling-converged barrier.

| Candidate | Job/complete step/time CST | Atomic-block residual (eV/A) | Scaled-cell-block residual (eV/A) | Dominant image | Sampled maximum above PO+ (meV/f.u.) |
|---|---|---:|---:|---:|---:|
| PO to M | 28298794/30/04:35:40 | 0.123663133 | 0.063466527 | 3 | 76.131749 |
| T-pattern-preserving flip | 28319570/25/04:16:23 | 0.134176879 | 0.124173502 | 4 | 54.238773 |
| T-pattern-preserving flip | 28319570/27/04:38:38 | 0.140771896 | 0.103344822 | 4 | 51.683480 |

- The preserving residual rises between steps25 and27 while its sampled
  maximum falls by2.555292 meV/f.u. Its dominant block is now atomic, whereas
  [the step14 observation](../network_update_20261009/README.md) was dominated
  by scaled-cell residual0.243947874 eV/A. The cell residual continues falling
  between25 and27. Consequently, the current rebound cannot simply be called
  persistent cell-block dominance.
- The preserving **global** maximum spring contribution is only0.000509103
  and0.000608650 eV/A at25 and27. A spring-spacing term alone does not account
  for the0.134--0.141 eV/A residual. These observations do not identify the
  optimizer's dynamical cause or prove a particular acceleration is required.
  Continue the existing chain; do not stop because of one rebound.
- PO to M is also now atomic-residual dominated. Its sampled maximum at
  image3 has physical tangential generalized force-0.365998838 eV/A, distinct
  from the projected ordinary-NEB residual. Thus even a later ordinary-NEB
  pass would not by itself certify that sampled image as a stationary TS.
- At the preserving maxima, the three rotated-T geometric patterns capture
  65.017% and64.510% of the translation-free squared displacement norm. This
  is descriptive coverage, not an energy fraction, local phonon spectrum,
  polarization or an independent prediction.

All three observations remain above the registered0.10 criterion. They do
not establish a final channel ranking, the network minimum, a stable Pbcn
intermediate, a TS certificate or G1 completion. No new symmetry classification
is inferred from these residuals or pattern amplitudes.

## Evidence and reproduction

The27 image records include6 reused endpoint records and different iterations
of the same channels; they are not27 new DFT calls or independent replicates.
Each observation contains9 POSCARs, `evaluated_chain.traj` and
`observation.json`. Numeric single-point calculators suffice for force replay;
proprietary UPF/orbital files remain on HF and are not redistributed.

`analysis_HF.json` and `analysis_local.json` were actually generated from the
same frozen evidence and unchanged analysis code. Their1923 floating fields
differ by at most2.842171e-14;365 nonfloating fields match exactly. Export
validates SCF completion, finite E/F/stress, unchanged physical input hashes,
ordered periodic geometry and optimizer force/energy replay before writing.
The original production archive checksum remains unchanged.

From the repository root, choose a new output file:

```sh
python -m scripts.analyze_hfo2_chain_observations \
  --root benchmarks/hfo2_channels/20261008/residual_update_20261009_0450 \
  --variants benchmarks/hfo2_channels/20261008/reference_variants \
  --gamma benchmarks/hfo2_channels/20261008/gamma_analysis/T_d0.01.npz \
  --output E:/TEMP/hfo2-residual-replay-new.json
python -m pytest -q tests/test_hfo2_chain_observation.py tests/test_hfo2_observation_analysis.py tests/test_hfo2_network_update.py
```

Those30 existing focused tests passed in29.39s; no implementation was changed
for this audit. The staged-tree archive preserves all37 evidence files byte
for byte, excludes user changes and reproduces `analysis_local.json` exactly;
its30 focused tests also pass in28.26s. Later changes are receipt/journal
metadata only, not source, tests or numeric inputs. Source/evidence hashes and limitations are retained in
[the delivery receipt](validation_delivery.json). This diagnostic adds no G2
or holdout tasks and does not enlarge the finite research matrix. Finish the
existing G1 chains before using them to train, rank or open the next gate.
