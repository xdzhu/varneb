# Primary-source comparator read-through — 2026-09-29

This checks the manuscript's four **quantitative external comparators** and
the G-SSNEB methodological boundary against original articles and the
previous page-specific BTO attribution audit. The RSC author PDF was not
re-fetched in this pass (the publisher endpoint returned 403); its Fig. 2
attribution is carried forward from that recorded audit and the publisher's
figure description. It does not
claim to reproduce their DFT protocols or certify every citation in the
bibliography. The manuscript and its plotted source tables were compared at
the current `main` revision; no calculator input or measured energy was
changed to improve agreement.

| Manuscript comparison | Original-article evidence | VARNEB evidence and allowed wording |
| --- | --- | --- |
| GaN B4→B1 and B3→B1, Qian *et al.* | [CPC 184, 2111–2118 (2013), DOI 10.1016/j.cpc.2013.04.004](https://uspex-team.org/static/file/Qian-vcNEB-2013.pdf): PDF pp. 3–6, Figs. 5, 6, and 8. At 45.7 GPa, the rotation-avoiding tetragonal B4→B1 route peaks at image 15, 0.34 eV/GaN; the hexagonal route is 0.39 eV/GaN. At 45.0 GPa the B3→B1 *joined* chain has three approximately 0.57-eV/GaN peaks (images 6, 15, 25) separated by two minima. Its GaN protocol is QE/PW91, ultrasoft potentials, 75 Ry, Γ-centered 8×8×6, 27 **interior** images, 0.03-eV/Å RMS force criterion (0.01 at the saddle). | VARNEB Fig. 4 uses 29 **total** images = 27 interiors at the adopted 45.7 GPa; the five PBE backends share the dominant peak topology but differ in pseudopotential/basis and use a 0.10-eV/Å maximum generalized-force criterion. The 0.3274–0.3385-eV/GaN group is numerically near 0.34, not a matched-protocol reproduction. Separate VASP hexagonal 0.3852 eV/GaN is a nearer transfer; its B3 diagonal 0.9553 eV/GaN is a negative mapping control, not a reproduction of three independently converged segments. Each energy is normalized by its own GaN formula-unit count and own initial endpoint. |
| CdSe RS→WZ, Sheppard *et al.* | [JCP 136, 074103 (2012), DOI 10.1063/1.3684549](https://henkelmanlab.org/pubs/sheppard12_074103.pdf): PDF p. 8, Fig. 11 and Sec. III D. The labeled 2.4-meV/atom DFT maximum is the **initial atom-dominated feature shared** by the black path and early red path; the later red cell-dominated route branches after state (c) and includes atomic rotation. DFT uses VASP/PAW/PW91, 455 eV, an eight-atom cell, and 10×10×10 Monkhorst–Pack sampling. The paper's method treats atomic and cell coordinates jointly with a simulation-cell-aware scaling. | Our PBE cell-mapping path gives 7.22 meV/atom, but neither the electronic contract nor the route/terminal structure is demonstrated identical. Do not call 2.4 meV/atom the full red-branch activation barrier or present 7.22 versus 2.4 as an accuracy score. VARNEB implements the established joint path physics under a documented hydrostatic enthalpy convention; it does not claim the original G-SSNEB formalism or general nonhydrostatic stress treatment. |
| BTO T→C, Akbarian *et al.* | [PCCP 21 (2019), DOI 10.1039/C9CP02955A](https://pubs.rsc.org/en/content/articlelanding/2019/cp/c9cp02955a), author-version pp. 7–8, Fig. 2: approximately 2.1 kcal/mol for a locally restrained PBEsol DFT tetragonal→cubic unit-cell distortion. The paper's NEB computation concerns oxygen-vacancy migration, not the T→C figure. See the page-specific `BTO_LITERATURE_ATTRIBUTION_AUDIT_2026-09-28.md`. | Our ABACUS/PBE seven-total-image T→C path rises monotonically by 0.0871289 eV/BTO. Similar energy scales do not define an activation barrier, an NEB-vs-NEB agreement, or an identical functional/constraint comparison. |
| HfO₂ T→PO, Liu and Hanrahan | [PR Materials 3, 054404 (2019), DOI 10.1103/PhysRevMaterials.3.054404](https://liutheory.westlake.edu.cn/pdf/Liu19p054404.pdf): PDF pp. 2–4, Methods, Fig. 3(a), Table I. The variable-cell **forward** T→PO barrier is 0.032 eV/HfO₂; reverse is 0.108 eV/HfO₂. Forty images, LDA/QE/USPEX, GBRV ultrasoft potentials, 50/250-Ry cutoffs, and a 0.025-eV/Å RMS image-force criterion differ from the separately plotted clamped-cell NEB. | Our 0.1291722-eV barrier belongs to a 12-atom Hf₄O₈ cell, i.e. 32.29 meV per HfO₂ after division by four. Its PBE/ABACUS and path contracts differ from the source. Numerical proximity supports a bounded scale comparison, not a statistical replicate or a clamped-cell claim. |

## Claim decisions

- The current manuscript's Qian and Sheppard descriptions, method qualifiers,
  formula-unit normalization, and caption warnings match these primary-source
  distinctions. The BTO and HfO₂ statements likewise keep distortion versus
  NEB and variable-cell versus clamped-cell results separate. No numerical
  manuscript correction was indicated by this read-through.
- The Fig. 4 Qian curve is approximately digitized from a published graphic;
  its line shape is a qualitative visual comparator, not released original
  electronic-structure data. Barrier bars use the reported 0.34-eV/GaN
  number instead of estimating a new value from pixels.
- This check closes the *four quantitative comparator* attribution pass,
  including the carried-forward page-specific BTO source audit. General-method
  citation coverage, final author-approved wording, CPiP route, and the CPC
  Program Library deposit remain separate submission gates.
