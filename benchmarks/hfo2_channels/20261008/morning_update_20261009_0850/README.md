# E029: second ordinary-converged channel and a changing switching profile

This is a dated G1 material update, not the matched-boundary H1/H2 test.
PO to M finished normally at06:29:53CST, Slurm28298794/exit0:0. Its immutable
terminal summary declares `converged=true`, `termination=force_threshold` and
fixed endpoints cached once. Exact raw-SCF force replay independently yields
0.0979044896eV/A at step39, below the unchanged ordinary0.10 criterion.
No CI, additional endpoint optimization or calculator retuning was used.

The common PO+ reference is-9783.249675811956eV per Hf4O8 cell, four formula
units, P=0/E=0. Original ABACUS/PBE100Ry/full10-auDZP and KPT/UPF/orbital/input
hashes remain unchanged. Export checks one genuine DSIZE=32 SCF per raw log,
electronic convergence, finite E/F/stress, exact ordered geometry and hashes.
All production directories and their original source archive remain unmodified.

| Source viewed from common PO+ | Complete step | Ordinary residual (eV/A) | Sampled maximum (meV/f.u.) | Ordinary pass |
|---|---:|---:|---:|---|
| PO to T, unchanged reverse view | 6 | 0.059881616 | 115.210164 | yes |
| PO to M | 39 | 0.097904490 | 71.582207 | yes |
| T-pattern-preserving flip | 48 | 0.204908776 | 39.069483 | no |
| T-pattern-reversing flip | 12 of continuation | 0.455054161 | 397.428755 | no |

These sampled numbers do not certify a final network minimum or a selective
mechanical window. Both switching chains remain unconverged; no new G2 or
holdout jobs are submitted. The ordinary criterion is not tightened.

## What the new complete frames resolve

PO to M has sampled forward/reverse maxima71.582207/142.966807meV/f.u. and
reaction energy-71.384600meV/f.u.; their difference agrees by construction.
M's ordered endpoint is P2_1/c at all declared symmetry tolerances
0.001/0.01/0.05Angstrom, with angle tolerance1degree. Its sampled maximum
at image3 has P1/P2_1/P2_1 tolerance-dependent groups and physical tangential
generalized force-0.383495004eV/A. It is **not** a stationary TS certificate
or an energy/sampling uncertainty bound, despite ordinary NEB convergence.

The preserving flip changes shape. At step48, sampled local maxima are at
images3 and5,39.069477/39.069483meV/f.u., with Pca2_1 at all three tolerances.
The central image4 is now lower,17.440944meV/f.u., and still Pbcn. Its
scaled-cell residual0.204908776eV/A remains the global largest. Therefore the
old central Pbcn peak must not be reused as the current highest image or
preselected as a stationary TS. A lower central snapshot is not yet a stable
intermediate: its gradient and transverse stability are unverified.

Compared with steps25/27, the dominant block has changed back from atoms to
cell rows. The atomic maximum at48 is0.140143245eV/A, cell maximum0.204908776,
and global spring maximum0.003588479. Central stress_xx is-0.007706589eV/A3,
versus-0.003890548 at27, while its first cell length grows5.486711→5.572298A.
The changing configuration/coupling invalidates a fixed single-block cause
for all iterations. It does not prove a launcher error, a cell-only optimizer
overshoot, or the benefit of a particular acceleration. Keep the current
bounded continuation; do not stop it merely because of this rebound.

The reversing continuation is atomic-dominated at image2 and still improving
in the inspected log. Its central Pbca structure is a12-atom observation;
the shared group label is not an identification with a24-atom literature
variant or domain wall. No wrapping, atom remapping or symmetry enforcement
is applied to any production image.

## Reproduction and historical integrity

This folder retains27 newly frozen image records, including6 cached endpoint
records, plus the terminal M summary. The four-channel report reuses the10
T--PO records,37 records in total; none are new SCFs or independent replicates.
The old E025 figure/report and its source-bound replay remain unchanged.
`analyze_hfo2_network_progress.py` reuses the same strict physical analyzer and
changes only progress-schema metadata and its obsolete one-pass caption.
It binds both source hashes; ordinary passes never automatically certify G1,
stationary TSs, sampling errors or independent predictions.

From the repository root, choose a new output path:

```sh
python -m scripts.analyze_hfo2_network_progress \
  --specification benchmarks/hfo2_channels/20261008/morning_update_20261009_0850/network_specification.json \
  --output E:/TEMP/hfo2-morning-progress-replay-new.json
python -m pytest -q tests/test_hfo2_network_progress.py tests/test_hfo2_network_update.py tests/test_hfo2_chain_observation.py tests/test_hfo2_observation_analysis.py
```

The focused suite passes35 tests/58.42s, including the old figure's exact
historical replay. A clean staged archive, excluding user changes, passes
1078 tests with2 skips and313 warnings in126.51s. Its network replay is
byte-identical to the local report. Simple three-frame HF/local analysis
matches1923 floating fields within2.842171e-14 and365 nonfloating fields exactly.
The same tested archive also runs the network CLI in a fresh isolated HF
directory:697 floating fields agree within2.842171e-14 and882 nonfloating
physical/source fields match exactly. Different package versions and38
structure-warning fields are retained in both original reports, not treated
as physical agreement or silently normalized. No remote pytest or installation
is claimed. Exact hashes and runtime metadata are in
[the delivery receipt](validation_delivery.json); unit tests do not replace
material convergence, branch certification or the preregistered predictions.

The09:08:42CST Slurm check confirms both switching jobs still Running.
Their later log values are step50/0.199713 and step13/0.435211eV/A. These
rounded logs are not substituted for the frozen complete step48/12 records.
