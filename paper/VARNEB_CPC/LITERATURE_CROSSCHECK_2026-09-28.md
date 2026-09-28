# Targeted literature cross-check — 2026-09-28

This is a targeted primary-source audit of three references used for numerical
or methodological comparison, not a claim that the entire bibliography has
been verified.

| Reference | Metadata | Manuscript claim and result |
| --- | --- | --- |
| Akbarian et al., DOI `10.1039/C9CP02955A` | Title, journal, year, and first author match the [publisher record](https://pubs.rsc.org/en/content/articlelanding/2019/cp/c9cp02955a). | **Corrected method mismatch.** The [article](https://pubs.rsc.org/en/content/getauthorversionpdf/C9CP02955A), pp. 7--8, reports 2.1 kcal/mol for a locally restrained PBEsol DFT BTO unit-cell T-to-C distortion. Its NEB calculation is for oxygen-vacancy migration. See `BTO_LITERATURE_ATTRIBUTION_AUDIT_2026-09-28.md`. |
| Qian et al., DOI `10.1016/j.cpc.2013.04.004` | The [published article copy](https://uspex-team.org/static/file/Qian-vcNEB-2013.pdf), p. 1, lists Qian, Dong, Zhou, Tian, Oganov, and Wang. | **Corrected BibTeX author mismatch.** The former entry incorrectly included Colin Glass and omitted four actual coauthors. Title, journal, volume, pages, year, and DOI were already consistent. |
| Liu and Hanrahan, DOI `10.1103/PhysRevMaterials.3.054404` | Title, authors, journal, year and DOI match the [publisher record](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.3.054404). | **Verified.** Table I and Fig. 3a of the [article](https://liutheory.westlake.edu.cn/pdf/Liu19p054404.pdf) give the 0.032-eV/f.u. *forward* T-to-PO barrier for VC-NEB; the reverse is 0.108 eV/f.u. The methods state 40 images, LDA/Quantum ESPRESSO/USPEX, and a 0.025-eV/Å RMS image-force threshold. This is not the separate clamped-cell NEB panel. |

The 0.032-eV/f.u. HfO2 source-data row and manuscript attribution therefore
remain unchanged. The BTO source-data row, manuscript text, generator and
regression check were corrected together; numerical energies and electronic
contracts were not changed. The visible main-paper result remains an
unmatched cross-method comparison, not an accuracy validation against an
identical DFT protocol.
