# Exact-SCF observations, not converged barriers

Seven complete Hf4O8/P=0 observations contain 66 already-computed image
evaluations, with full numeric energy/forces/stress and raw-source checksums.
No new DFT was used for exporting or analyzing them. Original production
source and geometry were left untouched. At the archived steps, all ordinary
NEB residuals still exceed 0.10 eV/Angstrom.

Each folder retains copied POSCARs, an evaluated numeric-only trajectory,
and `observation.json`. `analysis.json` binds those observations, the original
T/Gamma references and both analysis-source hashes. The source SCFs remain
on hf; the public trajectory is sufficient for numerical force replay.
Proprietary pseudopotential/orbital files are not distributed here.

From the repository root, write to a fresh output file:

```sh
python -m scripts.analyze_hfo2_chain_observations \
  --root benchmarks/hfo2_channels/20261008/chain_observations \
  --variants benchmarks/hfo2_channels/20261008/reference_variants \
  --gamma benchmarks/hfo2_channels/20261008/gamma_analysis/T_d0.01.npz \
  --output observation_replay.json
```

The model coordinates use a fixed reference metric and one initial integer
gauge throughout the path, not separate nearest-reference wrapping at every
image. Parent patterns are geometric; T Gamma eigenvectors are stationary-T
reference modes, not local phonons or saddle eigenvectors. Neither set
partitions energy. Degenerate mode amplitudes depend on their basis; complete
subspace norms do not. Cell strain is recorded separately.

The current PO--M peak has only 52.84% of its parent-relative displacement
squared norm in the rotated-T triplet. This is evidence against treating the
same triplet plane as a complete landscape for every competing channel, not
proof of a different topology or failure of a specific irreducible mode.

See `docs/HFO2_CONTINUOUS_CHAIN_OBSERVATIONS_2026-10-08.md` for force attribution,
provisional forward/reverse energy baselines and the next decision gate.

The clean archived source passes791 tests (two optional tests skipped).
`hf_replay_check.json` independently compares the same seven observations on
HF ASE3.23.1b1 with the local result: maximum force difference5.55e-16eV/A,
Gamma coordinate difference1.78e-15sqrt(amu)A, and zero Green-strain difference.
The archive checksum and per-observation errors are retained in that record.
The initial archive correctly failed when globally ignored trajectories were
missing; these seven numeric trajectories are now explicitly versioned.
Cross-platform observation ordering uses an explicit casefold key.
