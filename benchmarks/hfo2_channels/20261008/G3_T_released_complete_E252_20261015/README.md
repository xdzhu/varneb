# E252: first complete released-Q0 transverse stencil, T escape candidate

2026-10-15. This consumes the first completed half of the already submitted
E181 experiment; it does not add a calculation or a sampling point. The
registered point is the existing E160 fixed-Q0 conditional-release candidate
associated with the T escape path, **not the tetragonal endpoint**.

## Measured result and its scope

All 136 T probes, covering 34 released directions, both signs and both
original amplitudes, passed the unchanged E190 native reader. Each SCF's
original energy, full forces/stress, geometry, six input assets and genuine
32-rank MPI identity were checked on HF. The original 380-file runtime and
ABACUS binary were also checked. The first 13 native source identities are
identical to the earlier E190 prefix; no result or provenance path was edited.

The existing released-gradient norm is **0.026299954459 eV/A**, below the
original 0.03 conditional criterion but not zero. In the original fixed
coordinate metric, the measured R34 block has the following diagnostics:

| Probe step (A) | Minimum eigenvalue of symmetric part (eV/A^2) | Negative eigenvalues | Relative reciprocity defect | Explicit symmetrization change, operator norm (eV/A^2) |
|---|---:|---:|---:|---:|
| 0.01 | 0.164064528478 | 0 | 0.000168928039 | 0.004153228791 |
| 0.02 | 0.157595936408 | 0 | 0.000385318890 | 0.010960780748 |

Both original nonsymmetric matrices and all 40-by-34 gradient actions are
retained, not repaired or clipped. The raw-block difference between the two
steps has operator norm **0.054468926314 eV/A^2**. Neither that difference nor
the reciprocity diagnostics is a total physical-error bound. These are joint
coordinate curvatures, not phonon frequencies.

The T point therefore meets the preregistered **necessary local raw-sign
condition** at both steps. This does not certify exact stationarity, all
other T-grid points, global branch uniqueness, the full conditional minimum,
an index-one transition state, or independent prediction accuracy. No step
was selected after seeing its sign.

## Partial batch, full T cohort

The native report is intentionally `--mode prefix`: the whole two-center
Slurm job was still RUNNING. At its 14:09:24 CST completion, it had reparsed
**136 T + 2 preserving SCFs**, with the third preserving reservation pending.
The preserving block is absent, not filled with zeros or inferred positive.
The allocation accounting is a running snapshot, not a completed-job cost.

The complete B3 prerequisite was already audited in E194. The preserving
released-curvature audit and final candidate forecasts still have to be
completed before the independent tests. The eight reserved slots per face
and the physical +0.5% condition are not opened by this review. This is not
the full B0--B5 experiment or evidence of a new predictive advantage.

## Cost and unchanged contract

T used its original **136** new derivative SCFs, with no center SCF and no
retry. Their completed-call proxy is **137.839362458148 core-hours**. The
parent's running accounting snapshot was 140.24 core-hours for the allocation,
not a T-only or final two-center total. The original 296 raw-peak derivatives,
G2 training, faces, B3 and all preparation remain additional. No speedup or
DFT saving is inferred from this check.

The original ABACUS/PBE/100 Ry/full 10 au DZP/Gamma 2x2x2/P0/E0 inputs,
zero-strain clamp, exact literal coordinates, 0.5 A radius and 1.6 A minimum
distance are unchanged. Ordinary NEB remains 0.10 eV/A. The same job28868827
continues serially on hf/hfacnormal01 with one node/32 ranks; it was not
resubmitted, cancelled or modified. No CI, installation, release or new
monitor was started. Scratch and independent exports are on C:, not E:.

## Provenance and reproduction

- `native_audit.json`: SHA256
  `7f2dc871917734f3c2f3826c58d70239ef02545a29363b8b9fac8d7902bb68d4`.
- `T_curvature.json`: original producer summary, SHA256
  `66f7cdc9648afbfc3b61fc63a3835ec5868350c61aa921e56feea578a8498767`.
- Original E181 `config.json`: SHA256
  `71a2e1b8a46b2ef903e4023d329109c5e74975aaf5bab54f4d4638d15c1add39`.
- Unchanged E190 native reader: SHA256
  `51d0452968394ff03d26efface7705ee31c7679215cb2204f6c69e13a5391d4f`.

The real native audit ran on HF against the original paths and produced
`20261015-varneb-G3-released-audit-E190-r1/T_complete_prefix_20261015_140628.json`.
The repository files preserve those raw bytes. An offline checkout can
recompute both actions, projected matrices, all 34 eigenvalues per step,
reciprocity, step disagreement and cost from the audited full gradients:

```text
python -m pytest -q tests/test_hfo2_T_released_curvature.py
```

Offline mathematical replay is not another in-place native SCF-log audit.
The existing native checker also accepts this actual complete T summary and
rejects mutated counts, gradients, blocks or a missing amplitude. No private
author mode arrays or pseudopotential/orbital binaries are redistributed.

Both the working-tree run and the independent source snapshot pass 63 related
tests across four modules, eight new; neither is a new full repository run.
The tested tree is `b93dd0c4d1e6e8873dbab209d09f6417d6745212`, archive SHA256
`4a56f306b7c05913da05e8b3a5447daf47c32b08472550085453365d48d1675e`.
Raw XML and `validation.json` retain the checks and source identities. The
previous E205 complete regression remains separately recorded, including
its two explicit skips. The manuscript and five user edits are unchanged.

A later metadata-only observation at14:22:54 CST showed T136/preserving10,
no failure records and the same running allocation. It does not upgrade
the **two** preserving records in this native audit to ten audited records.
The rest of the preserving matrix and final pre-label forecasts remain
unfinished; no independent test has been started by this review.
