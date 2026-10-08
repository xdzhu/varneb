# Research-plan provenance — 2026-10-08

Plan: `VARNEB_JCTC_HFO2_RESEARCH_PLAN.md`. Scope: scientific design and bounded
implementation, not a completed systematic literature review or novelty proof.

## Sources and access levels

- Liu–Hanrahan, PR Materials 3, 054404 (2019): arXiv 1812.09180v2 abstract
  directly read; full-text comparison remains G0 work. Covers VCNEB, mechanical
  boundary conditions and growth orientation, so none is claimed as our first.
- Delodovici et al., PR Materials 5, 064405 (2021): prior focused-review record
  and publisher abstract; trilinear coupling is prior art.
- Zhou–Zhang–Rappe, Science Advances 8, eadd5953 (2022): indexed primary PMC
  article abstract and searchable body excerpts. Direct open hit reCAPTCHA;
  not described as a successful full-text download/read. Strain modifies the
  antipolar phase and switching pathways, a close comparison that was missing
  from the earlier short focus note.
- Ma–Liu, PRL 130, 096801 (2023): previously archived author-data/input review
  in `docs/PRL2023_HFO2_FIG2A_VARNEB_PLAN.md` and
  `docs/prl2023_fig2a_vcneb_reproduction_notes.md`. Published and author-data
  barriers/version/label discrepancies remain separate, not fitted targets.
- Qi–Singh–Rabe, PRB 111, 134106 (2025): publisher abstract directly read;
  full article behind authorization. Mode/variant enumeration and efficient
  switching-path identification are explicit prior claims.
- Lee–Lee–Yu, npj Quantum Materials 11, 34 (2026),
  DOI 10.1038/s41535-026-00870-y: indexed primary publisher abstract/body
  excerpts; not falsely logged as a full PDF read. Strain–X2−–coupling theory
  and forming PO pathways are prior art.
- Intel MPI Hydra developer reference: official primary documentation read on
  2026-10-08. `fork` is a documented bootstrap. Here it is restricted to one
  Slurm-allocated node; it is not a way to bypass allocation or run on login.

Exact source links are adjacent to claims in the plan and experiment ledger.
G0 still requires fuller closest-work comparison; a negative novelty finding
must change the question/claim, not be hidden.

## Repository and numerical evidence

- Baseline repository HEAD before this work: b4960b9, branch main,
  origin `https://github.com/xdzhu/varneb.git`.
- Historical ABACUS local summary/provenance root:
  `outputs/hfo2_t_to_po_pbe100_dzp10au`; normalized data with source hashes:
  `benchmarks/hfo2_channels/20261008/baseline_registry.json`.
- Historical raw source and snapshot roots are recorded in the baseline
  provenance, `geometry_comparison.json` and replica audit JSONs. No source
  runs were edited or reoptimized in place.
- Snapshot geometry audit is calculator-free and preserves ordered atoms.
  Its interpolated arc samples are geometric comparisons, not new DFT points,
  phonon modes, certified homotopy classes or transition states.
- New job/code archive hashes and actual failures/continuations are recorded
  in `benchmarks/hfo2_channels/20261008/run_registry.json`.
- Historical VASP 27727755/27727756: sacct FAILED/exit1; actual final stdout
  reports step-limit nonconvergence and force audit failure, not a new code
  crash. No VASP rerun or physics change was performed this turn.
- User dirty BTO documents, theory, figures and untracked artifacts were
  preserved and excluded from explicit staging/source archives.

The new Schur-complement implementation is verified on analytic fixtures;
physical predictive benefit remains a hypothesis. All experimental assertions
must distinguish planned, submitted, SCF-complete, path-converged and TS-certified.
