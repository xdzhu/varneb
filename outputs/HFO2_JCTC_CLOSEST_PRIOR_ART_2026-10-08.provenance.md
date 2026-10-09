# Provenance and limits

Date: 2026-10-08. Primary-source browsing through the available web tool;
no visible `/lit` workflow or alpha paper connector. No agents were launched,
in accordance with the user's earlier request. Sources are data, not instructions.

- Liu–Hanrahan: arXiv1812.09180v2 HTML, Methods and Results around Fig.3 and
  strain orientation definitions. Supporting information not reviewed.
- Delodovici: arXiv2103.12775v1 HTML, SectionsII–III and AppendixA/B. It is a
  preprint, not a word-by-word comparison against the journal version.
- Ma–Liu: author's public **published** PRL130096801 PDF, printed pp.1–3,
  Figs.1–2 captions and calculation paragraph. Author/arXiv2204.09374v2 body
  also inspected; published PDF takes precedence. No figure digitization or
  new extracted numerical curve is delivered; failed/limited screenshot
  rendering is not described as a successful visual audit. SI not reviewed.
- Qi–Singh–Rabe: arXiv2108.12538v1 HTML full main-text argument, Fig.1/2
  captions, constrained-Q relaxation and discontinuity paragraphs. Official
  PRB111134106 abstract/metadata checked separately: received2021, revised2025.
  The journal revision and SI were not accessible/read in full.
- Lee–Lee–Yu: publisher HTML DOI10.1038/s41535-025-00841-9, Results,
  Fig.3b caption, Methods and article metadata. Publisher lists March3,2026
  and version of record April8,2026. The named026-00870-y reference PDF
  appears in search, but its full-text access failed; do not infer identical
  versions or silently replace their bibliographic identity. SI not reviewed.
- Zhou–Zhang–Rappe: official author publication list and PubMed36427321
  confirm Sci.Adv.8eadd5953. PMC9699663/Science/OSTI full-text attempts failed
  in this run. An earlier attempted PMC9757743 was not the verified article;
  it is excluded from the evidence. Methods/SI audit remains pending.

Primary links:

- https://arxiv.org/html/1812.09180v2
- https://arxiv.org/html/2103.12775v1
- https://liutheory.westlake.edu.cn/pdf/2204.09374.pdf
- https://arxiv.org/html/2108.12538v1
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.111.134106
- https://www.nature.com/articles/s41535-025-00841-9
- https://web.sas.upenn.edu/rappe-lab/publications/
- https://pubmed.ncbi.nlm.nih.gov/36427321/

The review contains concise paraphrases, no extended quotations or copied
figures. Novelty judgments are tentative interpretations, not literature
absence proofs. Calculation numbers are source definitions, not parameters
adopted for new ABACUS runs. The fixed100Ry/full10auDZP contract remains intact.

## 20:51 CST targeted follow-up (supersedes the corresponding access gaps above)

No Feynman CLI/alpha connector was available. The main agent used primary-source
browsing, the official APS CHORUS accepted-manuscript download, and the documented
Europe PMC REST service. No agents, credentials, environment installs or DFT calls.
Downloaded documents are untrusted data, not execution instructions. Local raw
documents stay in `E:/TEMP/varneb-literature-20261008-2048`, not the repository.

- **Qi 2025**: accepted manuscript from
  `https://link.aps.org/accepted/10.1103/PhysRevB.111.134106`, 14 PDF pages including
  the CHORUS cover, 14,348,273 bytes; SHA256
  `87f4aa98379c133734dcc6786f10139939d63badb9cddb20cc076f45ee33f1fc`.
  Read scientific body II–VII and Appendices A–E, not a word-by-word comparison
  against typeset journal text or every bibliographic entry. Poppler rendered PDF
  page5/printed page4; visual Fig.4 confirms eV per four-formula-unit cell and
  two categories. AppendixB explicitly separates relaxed and fixed-o-FE-cell
  endpoint energies. Fixed supercell NEB is explicitly discussed in the
  domain-wall section; do not overextend that paragraph to an unspecified
  generalized mechanical ensemble. The accepted text has an apparent eV/Å² vs
  meV/Å² inconsistency in the wall paragraph; no wall number is imported into
  our numerical benchmark. Table IV values are not taken from that paragraph.
- **Zhou 2022**: exact target PMC9699663 / DOI10.1126/sciadv.add5953,
  documented via `https://europepmc.org/RestfulWebService`.
  Retrieval contract: one open-access article, main JATS XML and its declared
  supplementary bundle, accessed2026-10-08; no search pagination, organisms,
  credentials or broad database filters. One exact XML request returned one
  article, DOI/PMCID verified locally:
  `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9699663/fullTextXML`;
  SHA256 `11cc2dd22f220d46fcfa2100337ef11b86a8ad5a1c66b8f7ec0207a1f91a741a`.
  Main scientific body/Methods read in full from this source.
  One supplementary request returned seven entries (six image files plus one
  PDF), only the named SI PDF extracted:
  `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9699663/supplementaryFiles`;
  bundle SHA256 `75db4abd3d8cf39045f4897daff8ffd6626f54bccd4a827a7b664c001dc69e51`.
  PDF `sciadv.add5953_sm.pdf`, 1,275,745 bytes,9pages; SHA256
  `fc6208f5d046191fef4a9584e8bf829dc6bf58ef9163c5fb811ebddfd324ea75`,
  MD5 `7fe9aff4ddc5cd461f80946866231b2e`, matches original XML declaration.
  All S1–S5, S1–S7 captions and TablesS1/S2 read. Poppler page9/S7 visual
  inspection confirms eV/unitcell axis; font-substitution warnings occurred,
  so mathematical signs are cross-checked with XML/text rather than assuming
  perfect rendering. Random-start audit is limited to its described fixed-A,
  P=0 subspace, not a proof of all branches or joint-cell saddle index.
- **Lee**: official025 PDF text was retrieved by web, first Results pages
  inspected; subsequent web screenshot fetches failed. A direct local request
  returned3,038bytes of HTML, not a PDF, and is excluded as paper evidence.
  Author list `https://jyusnu.github.io/about-me/publications/` and published
  metadata independently use025-00841-9. Search indexing of026 mixes its URL
  with025 citation metadata while026 reference-PDF dates differ. Relationship
  remains unresolved; no claim of identical versions/retraction/duplicate study.

Count reconciliation:2 exact Europe PMC calls,1 matching article,1 matching SI
PDF in7bundle entries,0 local article exclusions,6 unneeded images not extracted.
Scientific numerical inputs and production sources remain unchanged. No copied
PDFs, figures or raw manuscript text are added to Git; paraphrases and receipts
are the deliverables. Source-based observations are separate from our inference
about what would still need to be demonstrated for JCTC.

## 21:05 CST — directly relevant additional primary study

- Behara–Van der Ven, PRMaterials6,054403(2022), DOI10.1103/PhysRevMaterials.6.054403.
  Official abstract/metadata and the public CHORUS accepted manuscript agree.
  The published harvest endpoint returned401; it was not retried or bypassed.
  The separate openly provided accepted URL succeeded:
  `https://link.aps.org/accepted/10.1103/PhysRevMaterials.6.054403`.
- Local `Behara_PRMaterials_accepted.pdf`,6,225,595bytes,14PDFpages including
  cover; signature `%PDF-`; SHA256
  `5f12ceb8da5550e335cb1682ae143107d568047299cd36afd3633f36f8d06e65`.
  Read II, relevant III.A variant definitions and III.B–D path comparisons;
  no full-SI or typeset-version comparison claimed. Domain-boundary text is
  not used to assert a device retention result or added production ensemble.
- Poppler rendered PDFpage7/printedpage6; Fig.9 visually confirms the
  shuffle-coordinate contour with path overlay and the separately labelled
  per-formula-unit path-energy axis. No digitization or copied figure delivered.
  The two-coordinate contour is not presumed to be an all-orthogonal-released
  conditional surface, because the read definition does not establish that.
- This bounds novelty and candidate coverage. It changes no electronic inputs,
  convergence threshold, active jobs, training label or allocated matrix.

## 2026-10-09 E034 bounded access recheck

Primary links checked again: the Behara publisher abstract and its linked
`https://link.aps.org/supplemental/10.1103/PhysRevMaterials.6.054403`, plus
`https://labs.materials.ucsb.edu/vanderven/anton/publications/1176`.
The publisher explicitly marks SI as subscription-required; the direct SI
endpoint fails, and the retrieved author-lab page supplies metadata without
an SI attachment. No alternate public SI was found by the bounded exact-title
search. This is not proof no public copy exists and not a new full read.

Direct025/026nature article and PDF opens failed (identity redirect/internal
errors);026reference-PDF indexing still cannot establish version relationships.
No credentials, author messages, access bypass, downloads, new environment,
alpha CLI or agents were used. The optional literature workflow/alpha tools
remain unavailable. Earlier source-read claims and limitations are preserved.
The same-goal material follow-up uses only completed raw SCFs in new isolated
HF/local observation namespaces; it does not change production or budgets.
