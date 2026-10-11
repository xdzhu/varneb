# E071: fixed-frame mode coverage of the ordinary-converged clamped M chain

2026-10-11. This is a zero-DFT **training-data analysis**, not a new predictor,
G3 bottleneck selection, energy decomposition, local-phonon calculation or TS
certificate. It uses the actual nine terminal frames of job28722320, whose
ordinary clamped residual is0.099214088eV/Angstrom. The full final-segment HF
audit is pinned by SHA56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c.

## Source-level coordinate transport, not structure repair

G2 seed preparation prescribes Cartesian/cell axes `(old z, old x, old y)`
with **identity atom permutation**. Directly using those coordinates with an
unrotated G1 T reference gave a misleading large displacement and tiny
coverage. The rejected trial and the exact executed trial source are retained
as `rejected_untransported_trial.json` and `executed_untransported_trial.py`.
They are not physical findings or the authoritative analysis.

The corrected [analysis.json](analysis.json) uses the pinned original seed
manifest and inverse prescribed proper frame on **geometry-only copies**:
`H_old = R.T @ H_G2 @ R`, `r_old = r_G2 @ R`. Calculator results are not
attached to these copies. Native E/F/stress and ordinary residual are audited
in the original production frame before any coordinate transport. No atom
assignment, image-wise MIC, cell fit, symmetry repair, optimizer step or DFT
input is modified. Unregistered frames and atom permutations are rejected.
Only the original0/+1% training boundaries are admitted; +.005 was not read.

All four E023 Cmma registrations (IDs103,106,114,127), chosen from G1 endpoint
geometry, are kept. The registration itself is unchanged and SHA-pinned.
The original T Gamma reference, author-source bytes and normalization policy
are also unchanged. Authors' raw/transformed mode arrays are not distributed.

## Actual scalar result

Fractions below refer to squared **translation-free atomic displacement in
the same original-T mass metric**, not energy or cell-strain contributions.
The highest sampled image is3, not a certified saddle. The Cmma interval is
the range over four retained registrations, not a statistical error bar.

| Atomic representation | Rank | Image3 fraction (%) | M endpoint fraction (%) |
|---|---:|---:|---:|
| T geometric triplet | 3 | 57.9237 | 24.4400 |
| Two complete lowest T optical doublets | 4 | 12.0962 | 24.2409 |
| Registered full Cmma four-direction span | 4 | 36.4746–36.4865 | 56.3537–56.3667 |

Atomic coverage changes along this real escape candidate. Neither a single
favourable point nor the different ranks justify a stronger-reference failure,
full-path energetic accuracy, nonlinear conditional-surface claim or predictive
advantage. Full Green strains are recorded separately. This does not choose
new modes after reading holdout labels or add curvature/phonon/chain budgets.

## Reproduction and verification scope

```text
python -m scripts.analyze_hfo2_clamped_reference --case benchmarks/hfo2_channels/20261008/clamped_M_terminal_E068_20261011 --source-root AUTHOR_FILES --registration benchmarks/hfo2_channels/20261008/cmma_path_mapping/analysis.json --variants benchmarks/hfo2_channels/20261008/reference_variants --gamma benchmarks/hfo2_channels/20261008/gamma_analysis/T_d0.01.npz --terminal-audit-sha256 56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c --clamped-seed-manifest benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/clamped_seed_manifest.json --clamped-seed-manifest-sha256 46c19e628b71e0a44a595be5f0ec9a2518e0719db0d1b286b294f5160f4af9a4 --output NEW_REPORT.json
```

`AUTHOR_FILES` are the same six public, pinned files described in the
[original registration](../cmma_path_mapping/README.md), outside this repository;
no author executable is run. Output must be new. Local analysis reparses nine
portable native INPUT/KPT/STRU/logs and geometry/E/F/stress; the complete six
physical-file bytes and earlier91SCFs are **historical HF audit evidence**,
not newly checked locally. No scheduler query or mutation occurs in the tool.
Validation and independent-archive reproduction are recorded separately in
`validation.json`; passing mocks alone are not a material-data audit.

The original finite G2/G3, strong-control and independent-prediction objective
remains active and incomplete. No extra DFT, modified physical parameter,
production-source overwrite or holdout access is part of this milestone.
