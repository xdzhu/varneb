# Published Cmma reference: inert structure and mode audit

2026-10-09, experiment E022. This is preparation for the strong fixed-parent
control in [prediction protocol v2](../../../../docs/HFO2_PREDICTION_PROTOCOL_V2_2026-10-09.md),
not a material prediction result, new DFT run or proof of a method advantage.

## Sources and reproducibility

[Qi–Rabe PRL 135, 046101](https://doi.org/10.1103/9759-kp38), its
[arXiv v2 ancillary SI](https://arxiv.org/src/2412.16792v2/anc/Supplementary_Materials_0707.pdf),
and the authors' [public data repository](https://github.com/yuboqiuab/unstableflatband/tree/a438e4ecf63cddfa68d4ae8ea83d4c19bd238969)
are distinguished from the other Qi–Singh–Rabe paper. The author repository
is pinned at `a438e4ecf63cddfa68d4ae8ea83d4c19bd238969`. The five reviewed raw
files, URLs, Git blob hashes, bytes and SHA256s are listed in [analysis.json](analysis.json).
No explicit redistribution license was visible at the repository root, so raw
structures/eigenvector files and copied figures are not distributed here.

Locally the reviewed originals are at
`E:/TEMP/varneb-cmma-source-20261009-0027`; this host path is not a portable
repository dependency. Obtain the pinned files using the report's public
URLs and keep their original bytes. With existing ASE/pymatgen/numpy/scipy/spglib:

```text
python -m scripts.audit_hfo2_cmma_reference --source-root AUTHOR_FILES --output NEW_REPORT.json
python -m pytest -q tests/test_qe_modes.py tests/test_hfo2_cmma_reference.py
```

Only our reviewed parser/audit runs; no author code, QE, pseudopotential or
DFT calculation is executed. See the [import contract](../../../../docs/QE_GAMMA_MODE_IMPORT.md)
and [delivery receipt](validation_delivery.json).

## Measured reference evidence

- Both source cells are periodic Hf4O8. Space-group number 67 is reproduced
  at `symprec=1e-4/1e-3/1e-2 Å`, angle tolerance 1 degree. The installed
  spglib standard-setting symbol is `Cmme`; the paper's `Cmma` name is retained
  without treating a setting label as a different phase.
- The source README exchanges paper/source y and z. An explicit proper
  frame, cell basis, translation and species-preserving permutation connect
  the QE input to the author POSCAR; maximum site error is `8.292e-7 Å`,
  lattice component error zero. Nothing is scaled, strained or symmetrized
  to force a match.
- An independent reconstruction of the literal six-atom SI Table S2,
  using integer determinant-two conversion to the existing 12-atom reference,
  agrees with the author POSCAR within `1.455e-5 Å` site error and `8.964e-6 Å`
  cell component error. This is an equivalent reference representation,
  **not** an expansion of a production cell or a new phonon calculation.
- Table S1's literal `x=0.05000` remains unchanged. It disagrees with the
  pinned author structure and the independent Table S2 reconstruction, both
  of which support the corresponding `x=0.50000` site. The discrepancy is
  preserved rather than silently corrected in the paper or an input.
- The first Gamma block contains 36 modes. Use author masses Hf 178.49 and
  O 15.9994 amu, not ASE's default O 15.999. QE flvec contains normalized
  displacements; after restoring the mass metric, the measured Gram defect
  is `4.009e-6`. Printed normalization defect is `1.973e-6`; both are retained.
- Three rigid-translation modes are identified by subspace overlap
  `0.9999923`, not by declaring every small negative frequency an instability.
  Four resolved negative reference frequencies are -7.218341, -6.407910,
  -3.251593 and -2.493489 THz. Under the explicit C-centering translation,
  the first two have odd weight above 0.999992 and the next two even weight
  above 0.999992. This character supports folding; it is not an inferred irrep
  label or a claim that a conventional-cell Z mode was computed at Gamma.

## Limits and next gate

Author files are associated by their pinned co-location; underlying flfrc and
DFPT generation were not audited. The co-located PWSCF output identifies 6.3,
but this alone does not prove the version that generated matdyn.modes.
The historical QE 6.3 source and current matdyn docs agree on the vector
convention. Later q blocks are not imported. Complex values are preserved;
there is no QR/ASR/real-phase projection or exact-basis certification.

**Production atom order, lift and path gauge are not yet mapped.** The complete
four-direction strong fixed-parent material benchmark therefore has not passed.
The author LDA geometry/modes are candidates for representation only: all
material energy, force and curvature evidence must remain from our fixed
PBE/ABACUS 100 Ry, full 10-au DZP contract. No reference phonon matrix, strain
condition, new chain or prediction label is added by this audit. Full-variable
TS certification, joint elimination and independent predictive advantage remain
open, and a strong fixed-parent baseline may ultimately perform equally well.
