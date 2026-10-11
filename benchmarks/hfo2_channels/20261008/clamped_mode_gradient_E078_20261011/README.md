# E078: actual gradient coverage is not geometric coverage or energy accuracy

This zero-DFT training analysis reuses all nine portable native SCFs of the
ordinary-converged, zero-strain PO-to-M chain28722320. Its actual audit SHA
is56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c.
The original100Ry/full10auDZP/PBE/P0/2x2x2 contract, ordered geometry,
production runtime, force target0.10 and all E023 reference choices remain
unchanged. No G3 selection, curvature, unseen condition or forecast occurs.

## Why this additional observable matters

E071 showed that the common-T atomic displacement at sampled peak image3
is only36.4746--36.4865% represented by the registered Cmma rank4 span.
This did **not** prove energetic failure. Here the same span represents
93.5085--93.5087% of the local atomic energy-gradient squared norm.
The two fractions have different denominators and answer different questions.
Cmma therefore remains a serious strong control, not a straw-man reference.
No conclusion about conditional energies or independent predictive gain is
drawn from a force projection. The tiny range is over all four retained
registrations, not a statistical/numerical error bar.

| Atomic representation | Rank | Peak retained squared gradient (%) | Peak omitted gradient norm |
|---|---:|---:|---:|
| T geometric triplet |3|9.028221|0.046129817|
| Two complete lowest T optical doublets |4|32.709627|0.039673897|
| All registered Cmma frames |4|93.508466--93.508708|0.012322362--0.012322592|

Norms use eV/(sqrt(amu)*Angstrom), the fixed original-T mass coordinate.
The common internal norm is0.048364688; these are not ordinary eV/Angstrom
max-vector residuals, physical force error bounds or energy fractions.
Rank3/4 results are not a same-dimension accuracy benchmark. All nine images
and the measured norms, ranks and rigid-translation terms are in `analysis.json`.

## Dual transport and an explicitly restricted coordinate gauge

For ASE row cells, a reference displacement u produces current displacement
`u @ inv(H_T) @ H`. The conjugate atomic energy covector is therefore
`g_q = -F @ (inv(H_T) @ H).T / sqrt(m)`, with q=sqrt(m)*u. Treating the
force like a displacement violates virtual work. Native F and geometry
use the prescribed inverse G2 cyclic frame; no per-image rotation fit,
atom rematching or calculator attachment to wrongly framed copies occurs.

Rigid translations are removed as a declared coordinate gauge. The source
basis is restricted to the same gauge, with rank3/4 explicitly required
to survive. This is **not** an ASR repair of source eigenvectors or forces.
The Cmma source translation-overlap norm0.003685405 and minimum restricted
singular value0.999993212 are disclosed; the raw gradient's translation
component is also retained as a diagnostic. The E071 displacement analysis
and its original source bases are untouched. Arbitrary axes within the
restricted span are not named new phonon eigenvectors.

All atomic directions here are mechanically movable. A mode-truncated
complement is not a clamped support reaction; cell derivatives/reactions
are separate and are not silently combined with the mass-gradient norm.
Even a large retained fraction certifies neither stationary branches,
omitted curvatures, negative-mode index, a continuous TS nor barrier accuracy.

## Replay and evidence boundary

```text
python -m scripts.analyze_hfo2_clamped_mode_gradient --case benchmarks/hfo2_channels/20261008/clamped_M_terminal_E068_20261011 --source-root AUTHOR_FILES --registration benchmarks/hfo2_channels/20261008/cmma_path_mapping/analysis.json --variants benchmarks/hfo2_channels/20261008/reference_variants --gamma benchmarks/hfo2_channels/20261008/gamma_analysis/T_d0.01.npz --terminal-audit-sha256 56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c --clamped-seed-manifest benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/clamped_seed_manifest.json --clamped-seed-manifest-sha256 46c19e628b71e0a44a595be5f0ec9a2518e0719db0d1b286b294f5160f4af9a4 --output NEW_REPORT.json
```

AUTHOR_FILES are the same six pinned public files in the original E023/E071
registration, not redistributed arrays. No author executable is run. Source,
seed, registration, portable input/log/geometry/E/F/stress and ordinary-force
audits are replayed; full licensed asset bytes and all91preceding SCFs remain
historical HF evidence, not new local reads. Unknown training receipts reject
before label access. Fresh output is required. No subprocess, scheduler or
calculator launch is needed by this analysis.

Independent analytic pair-energy finite differences, wrong-dual negative
control, oblique frame/atom-order covariance, mass scaling, subspace rotations,
translation admixture, zero-gradient policy and actual native-frame replay
are tested. Validation records distinguish local/clean-delivery tests from
material replay; they are not a new full-project/HF regression or DFT
finite-difference consistency claim. Units-skill static checks are heuristics,
not uncertainty estimates. The full G2/G3/strong-control/independent-forecast
scientific objective remains active and incomplete.

Final validation:68related tests pass in the working tree and the actual
independent Git archive cwd, with direct exit0 and no failures/errors/skips.
The actual replay report is byte-identical to `analysis.json`; its SHA and
recoverable external location are in `validation.json`, avoiding a duplicate
data file in this case. The current17-page paper compiles0 and rendered
pages15/16/17 are inspected; this is not an all-page final submission audit.
