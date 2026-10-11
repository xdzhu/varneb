# E105: the M escape peak is not its final-well response

2026-10-12. Uses the actual ordinary-converged +1% M continuation28727368,
step22/.098588469684eV/Angstrom, and unchanged zero-strain M28722320.
The historical E104/E093 full native audit covers154fresh SCFs,338runtime
files, the six physical input bytes, actual32MPI and complete E/F/stress.
Each scalar analyzer reparses the nine portable terminal SCFs here; this
is not a second full154-call or licensed-asset audit. No new DFT, submission,
release or dependency write is made by this analysis. Production flips run
independently under the original contract.

The same original T chart/masses/integer lift, rank3 geometric T triplet,
both complete T optical doublets (rank4), and ALL four rank4 Cmma frames
103/106/114/127 remain. No new Gamma calculation, force/mode repair, G2-based
reference choice, G3 bottleneck selection or physical+.005holdout access.
Author arrays stay outside Git and their six original hashes are checked;
no author code is executed or redistributed.

| Highest sampled image quantity | Zero strain | +1% strain |
|---|---:|---:|
| Forward discrete peak, meV/HfO2 |98.458122|97.227719|
| T triplet displacement squared fraction, % |57.9237|54.9251|
| T triplet atomic-gradient squared fraction, % |9.0282|3.9008|
| Cmma atomic-gradient squared fraction, % |93.5085--93.5087|91.5873--91.5896|

The observed forward peak increment is-1.230403416meV/HfO2, versus
-22.206478717meV/HfO2 for the final-minus-initial well increment. The
pre-existing B1 fixed-absolute-peak/final-following nulls give-6.806502 and
-22.206479meV/HfO2, respectively. These are seen-training comparisons, not
independent prediction wins or proof that the small M peak change exceeds
a measured barrier uncertainty. The matched T escape peak rises12.347027
meV/HfO2; both escape routes remain. The two +1%flips are still incomplete,
so the full switching/escape selectivity result is not yet available.

At the +1% sampled peak, the internal atomic-gradient norm is
.035561530520eV/(sqrt(amu) Angstrom). Gradient and displacement fractions
have different denominators and exclude open-cell response. The strong
fixed Cmma controls must not be replaced by the weaker geometric T control
to manufacture a method advantage. Ordinary .10 is not TS authentication,
conditional transverse stability or a barrier-error certificate.

Each profile has its OWN prescribed substrate/relaxed PO+ zero and original
normalized generalized arc. The PO absolute-energy increment is
.027226009508eV/cell; different strains do not share an identical initial
SCF. Equal normalized coordinates do not match structures. Native raw stress
gives the fixed-open-row, own-PO-referenced partial-396.119757meV/HfO2 per
unit strain at the strained peak; this is not a relaxed barrier derivative.

All six signed-work decompositions retain omitted atoms, rigid translation
and open-cell terms. Independent native/dual contractions differ by at most
5.551115123e-17eV/cell. Full work to image3 is75.922579meV/HfO2, with
energy-minus-work21.305139 there and2.658882 at the final image. These coarse
quadrature diagnostics are neither DFT nor barrier uncertainty estimates.

## Reproduce without a calculator

Run from the repository root with the original pinned author files available
outside Git. Use fresh output paths. No calculator or scheduler is called.

```sh
python -m scripts.analyze_hfo2_clamped_reference --case benchmarks/hfo2_channels/20261008/clamped_reversing_release_E104_20261012/plus1_M_terminal --source-root AUTHOR_FILES --registration benchmarks/hfo2_channels/20261008/cmma_path_mapping/analysis.json --variants benchmarks/hfo2_channels/20261008/reference_variants --gamma benchmarks/hfo2_channels/20261008/gamma_analysis/T_d0.01.npz --clamped-seed-manifest benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/clamped_seed_manifest.json --clamped-seed-manifest-sha256 46c19e628b71e0a44a595be5f0ec9a2518e0719db0d1b286b294f5160f4af9a4 --terminal-audit-sha256 b86cdb72eab9e3d8ff2125e16e6fbb2d0954b34d54c0f05099a581c066108389 --output FRESH_DISPLACEMENT.json
```

The same arguments with analyze_hfo2_clamped_mode_gradient or
analyze_hfo2_clamped_path_work reproduce their reports. Raw stress uses:

```sh
python -m scripts.analyze_hfo2_clamped_external_work --case benchmarks/hfo2_channels/20261008/clamped_reversing_release_E104_20261012/plus1_M_terminal --audit-sha256 b86cdb72eab9e3d8ff2125e16e6fbb2d0954b34d54c0f05099a581c066108389 --output FRESH_EXTERNAL_WORK.json
python -m scripts.plot_hfo2_clamped_M_strain --output FRESH_FIGURE_DIRECTORY
python -m pytest -q tests/test_hfo2_clamped_M_strain.py tests/test_hfo2_clamped_T_strain.py tests/test_hfo2_clamped_mode_gradient.py tests/test_hfo2_clamped_path_work.py tests/test_hfo2_clamped_external_work.py
```

The raw-work allowlist now admits precisely six actual completed training
receipts, not unknown/held-out labels. The new plotter is a thin client of
the existing guarded T/M pair builder and renderer. Legacy T default limits
are unchanged; M uses adequate negative-energy margin. Figure9 has18native
records, aligned axes, unbold labels, visible upper/right spines and framed
translucent white legends; all four Cmma curves remain, with no smoothing,
fitting or confidence interval. See the prior figure_contract.md and final
validation.json for actual testing, source replay, PDF and visual review.

Production inputs/optimizers/geometry/source are untouched, as are the five
original dirty user files. No installs, CI, release/PyPI, new monitor or
journal submission. G2 is6/8 ordinary-converged; G3/conditional-stability,
strong-model comparison and prospective independent predictions remain.
The full scientific goal stays active and incomplete.

The initial full working-tree run detected19failures: the old five-receipt
expectation, deliberately pinned historical renderer dependencies after the
opt-in M limit change, and one historical analyzer source-hash comparison.
No DFT or material record failed. The reviewed renderer dependencies and
exact six-receipt assertion are updated; all original scientific pins and
historical figure/report bytes remain. The old/new analyzer bytes are both
explicitly checked, then all scientific values compared. The known historical
biaxial module LF/CRLF variants are independently checked and documented,
never applied to native DFT bytes. An initial README transcription of two
zero-M percentages is corrected from the actual CSV before delivery.
The first full result is retained; the corrected93-test focused run passes.
Independent final full-archive testing and raw scalar/figure replay are
recorded in validation.json rather than implied by that focused count.
