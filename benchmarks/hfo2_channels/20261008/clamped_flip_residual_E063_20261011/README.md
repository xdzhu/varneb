# E063: an atomic plateau with a growing open-cell residual

This is a finite zero-DFT diagnosis of complete steps 13/14/15 of job
28661019 (`PO_flip_T_pattern_preserving`). It is not a new optimizer trial,
a final barrier, a G3 bottleneck selection, or a held-out prediction.

## Measured result

All force values below are in eV/Å in the unchanged generalized coordinate
metric (`cell_scale = 5.12968067458423 Å`). The atomic block is conjugate to
fractional positions mapped through the reference cell, not raw Cartesian
atomic force. The open-cell block is projected onto the three released
deformation directions; fixed-substrate reaction stress is excluded.

| Complete step | NEB fmax / atomic-block max | Open-cell-block max | Perpendicular-only fmax at the same geometry | Spring-only fmax at the same geometry |
| --- | ---: | ---: | ---: | ---: |
| 13 | 0.221735900 | 0.101104650 | 0.220034151 | 0.005340569 |
| 14 | 0.220937676 | 0.158791935 | 0.219311661 | 0.005212512 |
| 15 | 0.220818025 | 0.204300388 | 0.219289136 | 0.005116883 |

The limiting vector in all three frames is image 5, atom index 8 (zero-based,
O in the fixed Hf4O8 order). At step 15, the signed projections of the
perpendicular and spring components on this vector are respectively
0.219286409 and 0.001531616; their sum is the observed 0.220818025.
Component maxima from different vectors must **not** be added.

The snapshot evidence rules out image-spacing/spring residual as the
dominant explanation of this short plateau. It does not show that the cell
is converged: its maximum residual grows and reaches 0.204300388. Removing
the spring at a frozen geometry would still leave ~0.219, above ordinary
0.10. This is not evidence that an actually reoptimized zero-spring path
would behave the same, nor does it establish a Hessian condition number,
an electronic-convergence defect, or a proven acceleration strategy.

Decision: leave the running FIRE segment unchanged through its registered
convergence/step/time boundary. Review its complete terminal chain before
choosing any same-chain continuation. A single rebound or these three
plateau frames are not a stopping criterion. Do not tune INPUT, change the
cell metric to make its residual appear smaller, or add a duplicate trial.

## Provenance and replay

Actual HF audit at 2026-10-11 01:48:31 CST used the original immutable E054
runtime. All 334 archived Python files were byte-compared to its registered
archive; the receipt records the scientific-source hashes. The three
nine-image snapshots were matched to exact ordered native SCFs with the
original six physical input bytes, complete energy/forces/stress and true
DSIZE32. External launches were prohibited during export/decomposition.
New DFT calls, job submissions, source/job/input mutations, G3 selections
and holdout reads are all zero. Proprietary PP/orbitals/binary files are
not copied into this portable case.

The observable transport tar is retained outside the repository at
`E:/TEMP/varneb-E063-observables.tar`; SHA256:
`ab1ec333cbe482c883c0cc3696a7e81caed81e28659bcaf2d12afda90684e5d7`.
Its exact HF copy remains in the E063 namespace. `executed_analysis.py`
and `executed_exporter.py` preserve the analysis bytes actually executed;
they are not another production implementation.

Portable replay needs only ASE/NumPy and the exported single-point data:

```bash
python -m scripts.analyze_hfo2_clamped_residual \
  benchmarks/hfo2_channels/20261008/clamped_flip_residual_E063_20261011/observations/step_0015 \
  /a/new/path/residual.json
python -m pytest tests/test_hfo2_clamped_residual.py \
  tests/test_hfo2_clamped_observation.py tests/test_hfo2_G2_M_handoff.py -q
```

The output must not exist. Native-raw replay tests check portable
INPUT/KPT/STRU/log/audit bytes, ordered geometry, unit-converted E/F/stress,
projected vector sums and all per-image component fields against the HF
results. Full six-file checks are HF evidence, not a claim that proprietary
assets are present locally. The initial tool-only test had 5 passes and
3 not-yet-exported skips; after actual export the related suite had
59 passes with no skips. The same suite in an independent Git-archive
directory also had 59 passes, no failures/errors/skips. All 181 newly
delivered source/material files were byte-identical; ordinarily normalized
root Git metadata and evolving JUnit outputs are not raw-byte attestations.
`validation.json` records the tested tree and archive; later additions are
documentation/receipts only, with no analysis-source or material mutation.
