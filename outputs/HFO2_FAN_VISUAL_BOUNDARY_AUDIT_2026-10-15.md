# Fan publication: visual audit and comparison boundary

E177, 2026-10-15. This closes the local-PDF/visual-access gap recorded in
[E070](HFO2_FAN_PHASE_PROXY_AUDIT_2026-10-11.md), not a new material result.
The publisher and author-hosted published PDFs are byte-identical. Main
scientific text pp.1–9 was read; all11pages were rendered as contact sheets,
with pp.4/6/7/8 inspected separately at higher resolution. Source hashes,
access limitations and validation scope are in the
[receipt](HFO2_FAN_VISUAL_BOUNDARY_AUDIT_2026-10-15.provenance.md).

## Source observations

| Published location | Checked observation |
|---|---|
| p.4, Fig.5 | Force parity axes are eV/Angstrom; caption identifies final-training-database configurations. The adjacent MAE values use meV/atom, not force units. |
| pp.6–7, Fig.9 | Each functional has its own relaxed PO strain reference. All lattice lengths are prescribed and uniformly scaled; atoms relax. Panels a/b show E; c/d show Delta E along lambda, in meV/f.u. |
| p.7, Fig.10 and accompanying text | The plotted switching-barrier approximation uses the intermediate T/Pbcn-to-PO gap. Panel b subtracts that functional's zero-strain barrier. |
| p.7, Fig.11 | Strained Pbcn dispersion and a Gamma unstable-mode schematic, not numerical eigenvectors for our cell. |
| p.8, Fig.12/Table2 | Mode-sign domain-wall classification; wall energies are meV/Angstrom^2, not homogeneous activation barriers. |

Source: [Fan, Zhu and Liu, published article](https://www.nature.com/articles/s41524-025-01647-w),
[published PDF](https://www.nature.com/articles/s41524-025-01647-w.pdf).
No curves, mode vectors or coefficient values were digitized/imported.

The publisher HTML fetched normally and contains no supplementary section,
MOESM or ESM attachment reference. Thus **no verified SI attachment was
found on this fetched page**; this does not establish that no supplemental
material exists elsewhere. No attachment URL was guessed.

## Consequences for VARNEB — our interpretation

| Comparison | Allowed conclusion / remaining requirement |
|---|---|
| Fan fully prescribed PO-based lattice vs our partial T-based clamp | Different mechanical ensembles. Equal electronic inputs would not make their strain percentages or barrier thresholds transferable. |
| Intermediate-phase gap vs B1 endpoint-response null | Distinct controls. Use a phase-gap explanation only when phase identity, boundary and energy role are established by existing data; otherwise mark unverified, not failed. |
| Force parity / energetic error | Keep training fit, independent prediction, force error and energy error separate. Fig.5 is not an independent-test label; this audit does not claim their entire study lacks other validation. |
| Lee SI vs present release comparison | [E108's completed Lee SI audit](HFO2_LEE_SI_VISUAL_AUDIT_2026-10-12.md) already establishes restricted relaxation as prior work. Do not repeat that audit or claim release-induced lowering itself as innovation. |

Our proposed contribution remains a same-initial-well switching/escape
comparison with measured freedoms, strong applicable controls and an
independent, prospectively frozen predictive test. It is not established
by relabeling a training reconstruction, compact surface or standard
quadratic elimination as a new physical result.

The existing manuscript introduction already makes the phase-gap/endpoint
and mechanical-boundary distinctions; no duplicate paragraph was added.
No DFT, submission, physical-input change, extra center, holdout access or
new training reservation was performed. Original jobs continue unchanged.

Next material step: audit the original preserving face when complete and
the actual atomic-only release results when available; compare all signed
frozen/atomic/joint changes, including missing candidates and costs.
