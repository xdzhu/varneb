# E025: dated common-PO channel observations

This is a **provisional material observation**, not a final competing-barrier
benchmark. Two newly frozen complete snapshots are combined with the converged
T--PO source and the retained reversing-flip snapshot. No DFT calculation,
job submission, production-source edit or parameter change was made to export,
analyse or plot this bundle.

All candidates start from the same ordered PO+ structure and energy
(-9783.249675811956 eV per 12-atom cell). P=0 and external E=0; four HfO2 formula
units per cell. ABACUS/PBE,100 Ry, full10-au DZP, original2x2x2 Gamma mesh,
UPF/orbital/input hashes and ordinary0.10 eV/Angstrom threshold are unchanged.
The source calculations use32 MPI; fixed endpoints are cached, not recomputed.

| Candidate, viewed from PO+ | Source job/complete step | Replayed ordinary residual (eV/A) | Sampled maximum above PO+ (meV/f.u.) | Group at sampled maximum, all three tolerances |
|---|---|---:|---:|---|
| PO→T | 28300425/6, reverse view | 0.059881616 | 115.210164 | Pc |
| PO→M | 28298794/20 | 0.156520577 | 82.136131 | P1 at0.001 A; P2_1 at0.01/0.05 A |
| T-pattern-preserving flip | 28319570/14 | 0.243947874 | 65.278562 | Pbcn |
| T-pattern-reversing flip | 28288063/10, retained seed | 0.700278856 | 410.743641 | Pbca |

Only the first source passes the ordinary residual criterion. The different
snapshot times are disclosed; these numbers do **not** establish a final
ranking, a minimum over all possible channels or a transition-state certificate.
No numerical/sampling error bounds have yet been measured. Straight lines
joining discrete energies do not certify the energies between samples.

## Measured mechanism observations

- The preserving snapshot's entire **sampled** profile lies below the relaxed
  T energy,81.321193 meV/f.u. above PO+. Its central sampled maximum is Pbcn
  and still has a0.243947874 eV/A scaled-cell residual. A Pbcn group at a
  maximum is not a stable Pbcn intermediate or a certified TS.
- The reversing centre is Pbca even though the three rotated-T geometric
  pattern amplitudes are essentially zero (captured fraction below1e-20).
  Vanishing coordinates in a truncated representation do not identify a
  cubic structure. The full structure must be checked independently.
- The original T reference is P4_2/nmc at each declared tolerance. Space-group
  results are reported for0.001/0.01/0.05 Angstrom and a fixed1-degree angle
  tolerance. No symmetry standardisation, wrapping, atom permutation,
  oxidation-state guessing or relaxation is applied to production images.
- The coordinates are geometric patterns, not polarization, phonon populations
  or energy contributions. The Green-strain norm is a geometry descriptor,
  not an elastic energy. Reversing the T→PO thermodynamic view preserves its
  original atom lifts and NEB force evaluation; arcs are source-specific.

## Provenance and reproduction

`PO_M_step20/` and `PO_flip_preserving_step14/` contain18 newly frozen image
records, including cached endpoints. The other19 records are reused from
`../converged_gap/gap_converged_step06/` and
`../chain_observations/PO_minus_T_reversing_step10/`. The resulting37 records
are neither37 new SCFs nor independent replicates. Each record retains
E/forces/stress, original input/log/geometry hashes and force/energy replay.
The production archive remains the cc536a2b source with SHA256
df4c12eef03d00d02c3ae357e3f47143b6500b96b7bd12f479f2e2083abd5220.

`analysis.json` is the common-initial-state network audit. The final
`material_observations_v3.json` adds an immutable pymatgen structure audit,
all three pattern amplitudes, translation-free displacement and Green strain.
The v3 suffix retains its dated generation identity; earlier development
reports are not included. Runtime versions and captured symmetry-library
warnings are explicit metadata. A replay may emit different counts of the
same cached deprecation warnings; both lists are retained separately. All
physical, source-hash and structural fields remain strict replay checks.

From the repository root, choose **new, nonexistent** output paths:

```powershell
python -m scripts.analyze_hfo2_channel_network --specification benchmarks/hfo2_channels/20261008/network_update_20261009/network_specification.json --output E:/TEMP/network-replay.json
python -m scripts.analyze_hfo2_network_update --specification benchmarks/hfo2_channels/20261008/network_update_20261009/network_specification.json --output E:/TEMP/material-replay.json
python -m scripts.plot_hfo2_network_update --specification benchmarks/hfo2_channels/20261008/network_update_20261009/network_specification.json --material-report E:/TEMP/material-replay.json --output E:/TEMP/hfo2-network-figure-replay
python -m pytest tests/test_hfo2_network_update.py tests/test_hfo2_channel_network.py -q
```

The [figure contract](FIGURE_CONTRACT.md) and
[figure/source-data bundle](../../../../paper/VARNEB_JCTC/figures/hfo2_G1_provisional_20261009/README.md)
disclose the omitted high-energy profile in panel(a); its complete energies
are retained in CSV and the table above. A clean staged-tree regression and
same-archive HF replay are recorded separately in `validation_delivery.json`.

## Consequence for the finite research plan

Finish the existing G1 chains before ranking channels or opening the G2
boundary matrix. After G1, compare the registered mechanical conditions;
then measure the selected bottleneck's joint atomic--cell response and use
the strong fixed-parent, endpoint and interpolation controls in the existing
prediction protocol. These observations select questions, not prospective
prediction successes. No additional G0 probes, decorative dense surfaces,
new materials or early CI are justified by this bundle.
