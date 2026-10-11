# E084: finite conditional-surface probe coverage (zero DFT)

This implementation closes the finite-reporting gap identified in E081.
It does not produce an HfO2 surface or select G3 modes/bottlenecks.
The original full-grid interpolation screen and physical calculation contract
are unchanged; the +0.005 material holdout remains unopened.

`local/` contains a training-only synthetic registration written before its
test labels, the 8-of-64 report, a failure/pending report and a verifier receipt.
The synthetic 9x9 node grid has 64 cells; eight exact test values give eight
passing point checks in eight cells, **not whole-grid readiness or a region
error bound**. The negative example keeps a -0.07 eV/cell signed energy error,
0.12 sqrt(amu)*Angstrom coordinate error, one failed calculation and one
pending point. It also rejects a post-registration refit, training-evaluation
reuse and an unregistered observation. Shared audit hashes are declarations
in this toy example, not actual DFT audits or physical-curvature measurements.

Run the commands in [the API guide](../../../../docs/CONDITIONAL_PROBE_COVERAGE.md)
from the repository root, always using a new output namespace. See
`validation.json` for actual tests, archive checks, source hashes and limits.
No scheduler writes, native DFT, model refits, new material matrix, material
contour or scientific predictor advantage are claimed.

Actual local and independent Git-archive focused tests each pass90; adjacent
local tests pass88. Full local regression exits0 with2021passes/one optional
`ase.io.abacus` skip; its scope includes preexisting user-untracked tests.
Five reviewed Python sources are normalized-text identical to the actual
archive, but four differ in CRLF/LF bytes. Both numerical reports and the
registration are byte-identical; verifier receipts differ only in source-byte
hashes. No full-archive or HF regression is claimed. Units audit has no
heuristic findings, not an independent physical error certificate.
