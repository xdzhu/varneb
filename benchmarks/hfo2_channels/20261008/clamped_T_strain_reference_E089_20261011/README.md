# E089: strain response in unchanged reference spaces

2026-10-11. Reuses the actual ordinary-converged +0.01 PO-to-T terminal
job28709788 (E087), latest step5/.0949479822eV/Angstrom. Nine portable native
SCFs are independently reparsed by the original displacement, gradient and
raw-stress analyzers. The historical full42-SCF/six-asset/runtime audit stays
in E087; it is not newly repeated or fully redistributed here. Zero-strain
job28692775/E080 is reused unchanged. No new physical DFT or scheduler writes.

Original T/chart/masses, geometric rank3 triplet, two complete optical
doublets rank4 and all four registered Cmma rank4 frames103/106/114/127 remain.
Author arrays stay outside Git, source bytes audited, no author code runs.
No G2-based reference selection, extra Hessian or +0.005 holdout access.

| Highest sampled image quantity | Zero strain | +0.01 strain |
|---|---:|---:|
| Forward discrete peak, meV/HfO2 |84.416624|96.763650|
| T triplet displacement squared fraction, % |99.0957|99.0045|
| T triplet atomic-gradient squared fraction, % |14.3013|69.7245|
| Two T doublets atomic-gradient fraction, % |20.3963|2.0694|
| All Cmma frames atomic-gradient fraction, % |54.0718--54.0765|10.5750--10.5766|
| Internal atomic-gradient norm, eV/(sqrt(amu) Angstrom) |0.010196606|0.038221369|

This is a seen-training contrast, not independent prediction or model accuracy.
Displacement and covector fractions have different denominators and exclude
cell response; different ranks are disclosed. Percentages at small gradients
are sensitive; all absolute retained/omitted norms remain in scalar data.
The ordinary .10 pass does not certify a continuous saddle or barrier-error
bound. No full channel-ordering or G3 bottleneck claim follows from two T bands.

The two prescribed strains belong to one matched substrate family, not the
same exact physical substrate or initial SCF. Each profile uses its own relaxed
PO+ reference; raw P0 energies are preserved. The measured PO increment is
0.027226009508eV/cell. The peak increment12.347026513meV/HfO2 remains a training
observation; frozen raw-stress partials must not be labelled relaxed barrier
derivatives. Both paths use their own original normalized arc, not pointwise
matched geometry. The unstrained T endpoint fraction is undefined; the strained
endpoint has a nonzero displacement from unstrained T and a defined fraction.

## Reproduction

Use the original E071/E078 analyzer commands with --case pointing to
clamped_plus1_T_terminal_E087_20261011 and its actual terminal receipt SHA
32a6e49a618e6e62632ae71a5e0ed093095254f6858221be2ade3d83cc194a28;
retain the same original registration, variants, Gamma file, external pinned
author files and seed-manifest arguments. Outputs must be fresh paths.
No calculator is called. The external-work analyzer also admits this fourth
completed training receipt explicitly; unknown and held-out receipts reject
before opening labels. The admitted set is exact and byte-pinned in tests.

```text
python -m scripts.plot_hfo2_clamped_T_strain --output FRESH_FIGURE_DIRECTORY
```

The plotter checks frozen scalar/native/receipt pins, same reference/ranks/
integer gauge and six input declarations, prescribed strain and fixed/open
cell contract. It draws all four Cmma registrations, no smoothing, fitting,
confidence interval or invented endpoint coverage. The source CSV contains
18image records including cached endpoints, not18independent replicates.
Figure dimensions183x202mm, aligned axes, unbold labels, visible upper/right
spines and framed translucent white legends follow the stated style.
The contract predates code. Validation and actual rendered review are recorded
in validation.json; preview and editable vector files are in the manuscript
figure directory. Raw production sources and five user-dirty files are untouched.

Two initial plot attempts rejected before exporting: a9.45e-18 unit-normal
serialization difference required an explicit1e-12 comparison, and a missing
channel metadata field was added. No structure or physical input was repaired.
A first test command named a nonexistent test file and ran zero tests; the
corrected actual local run passes60. The failures are recorded, not DFT failures.
Actual clean-archive replay preserves every numerical field; displacement
script and biaxial-curvature source byte hashes differ only under the existing
CRLF/LF checkout policy, with normalized text verified equal. The gradient
report and figure CSV are byte-identical. The raw-work report differs only
in the recorded biaxial source hash, not the measured values. Strict initial
comparison correctly rejected that hash difference; it was inspected and
documented, not silently dropped. Existing local PDFium renders the final
21-page manuscript; no package installation or other-language plotting.
Full G2/strong-control/forecast gates remain incomplete; the scientific goal
is active, not a claim of a submission-ready JCTC package.
