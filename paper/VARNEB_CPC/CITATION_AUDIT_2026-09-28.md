# CPC citation and comparator audit (internal, not a supplement)

On 2026-09-28 all 20 DOI records in `varneb.bib` resolved through the
Crossref works API. DOI, title, year, and author-family sequence were checked
against the committed BibTeX. Five author records needed correction; the
scientific sources below were then checked at the publisher, author-hosted
paper, or indexed original article. This is a metadata and selected comparator
audit, not a claim that every result in the cited papers was reproduced.

| Key | Corrected item | Source |
| --- | --- | --- |
| `Sheppard2012` | Xiao, Chemelewski, and Johnson restored; unrelated Terrell removed. | [Original paper, DOI 10.1063/1.3684549](https://doi.org/10.1063/1.3684549); author-provided local PDF, page 1. |
| `Ghasemi2019` | Authors are Arman Ghasemi, Penghao Xiao, and Wei Gao. | [Original article PDF](https://par.nsf.gov/servlets/purl/10186851). |
| `Larsen2017` | First author's family name is Hjorth Larsen, not Larsen. | [Crossref record for DOI 10.1088/1361-648X/aa680e](https://api.crossref.org/works/10.1088%2F1361-648X%2Faa680e). |
| `Makri2019` | First author is Stela Makri, not Spyros Makri. | [PubMed record and original-article DOI](https://pubmed.ncbi.nlm.nih.gov/30849914/). |
| `GarridoTorres2019` | Hansen, Boes, and Bligaard restored in place of three unrelated names. | [APS journal record](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.156001). |

Comparator interpretation checked against the source pages:

- [Qian et al.](https://uspex-team.org/static/file/Qian-vcNEB-2013.pdf),
  Figs. 5, 6, and 8: GaN tetragonal B4→B1 0.34 eV/GaN at 45.7 GPa;
  hexagonal B4→B1 0.39 eV/GaN; B3→B1 at 45.0 GPa has three peaks around
  0.57 eV/GaN in a joined path. The manuscript does not treat its one B3
  chain as three independently converged segments. The original's Methods
  (PDF page 3) specify QE/PW91, ultrasoft pseudopotentials, 75 Ry, a
  Γ-centered 8×8×6 mesh, 27 interior images, and 0.03-eV/Å RMS image-force
  convergence (0.01 eV/Å at the saddle). VARNEB shares the nominal pressure
  and image count but not that electronic or force-norm contract; 45.7 GPa
  was adopted from the reference, not recomputed as each PBE backend's
  coexistence pressure.
- [Sheppard et al.](https://doi.org/10.1063/1.3684549), Fig. 11: the
  2.4-meV/atom DFT feature is the small initial peak on the atom-dominated
  stage that the later cell-dominated route initially shares. The latter
  subsequently rotates atoms; our 7.22-meV/atom cell mapping has no
  demonstrated matching saddle or identical electronic contract. The
  original explicitly uses PW91; our run uses PBE.
- [Liu and Hanrahan](https://liutheory.westlake.edu.cn/pdf/Liu19p054404.pdf),
  Methods and Table I: 40 images, LDA/QE with GBRV ultrasoft potentials,
  and T→PO forward barrier 0.032 eV/HfO₂. Our 0.1291722-eV value is per
  12-atom Hf₄O₈ cell and becomes 32.29 meV/HfO₂ only after division by four;
  numerical proximity is not a matched-method replicate.

The manuscript's author/affiliation/CRediT and acknowledgement placeholders
refer to VARNEB itself, not these references, and still require author approval.
