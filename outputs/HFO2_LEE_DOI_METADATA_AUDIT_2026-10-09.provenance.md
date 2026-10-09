# Lee DOI identity and date correction — E035, 2026-10-09

This bounded registry check supersedes the *current metadata-identity* gap in
the earlier closest-work audits. It does not supersede their full-text/SI
access limits, reread the paper or prove novelty. See [the machine-readable
receipt](HFO2_LEE_DOI_METADATA_AUDIT_2026-10-09.json).

## Observed primary evidence

The [official Crossref REST documentation](https://www.crossref.org/documentation/retrieve-metadata/rest-api/)
defines `/works/{doi}` as an exact-record lookup. At 11:37:08 CST an anonymous,
no-follow GET for [026-00870-y](https://api.crossref.org/works/10.1038/s41535-026-00870-y)
returned HTTP301 with `Location: /works/10.1038/s41535-025-00841-9`.
Following it returned HTTP 200 and DOI 025; direct [025 lookup](https://api.crossref.org/works/10.1038/s41535-025-00841-9)
returned the same 23,901-byte payload/SHA256. The current record is Lee--Lee--Yu,
npj Quantum Materials 11, article 34. Its published/published-online date is
**3 March 2026**; deposit/index timestamps are **8 April 2026**. The [author's
publication list](https://jyusnu.github.io/about-me/publications/) independently
lists article 34 with DOI 025.

Therefore count **one current canonical metadata record**, cite 025, and do not
call April 8 the publication date or infer a version-of-record date from an
index/deposit timestamp. Earlier dated records are retained with this explicit
correction. The main manuscript already cites025 and was not changed.

## Retrieval, failure and interpretation boundaries

Two browser API opens could not retrieve the records; PowerShell anonymous
GETs supplied the documented endpoint fallback. Four bounded native GET calls
covered only these two identifiers: two payload downloads, one final-URL
diagnostic and one no-follow redirect diagnostic. No filters, pagination,
credentials, installs, access bypass, broad search or new agent work. The
database skill has no dedicated Crossref reference; its retrieval contract was
applied using the official REST documentation. Optional alpha tools were not
available and were not simulated.

The first 026 download succeeded but our exact-DOI identity guard rejected its
returned 025 record. That RuntimeException was **not an HTTP failure**; the
explicit 301/Location check explains the mismatch. Both raw JSON payloads stay
outside Git at the paths/SHA256 in the receipt; no abstract, paper or source
figure is redistributed.

The registry `relation` field is empty. Neither that field, equal JSON bytes
nor the redirect establishes the reason for the identifier change, a formal
correction/retraction, or equivalence of previously indexed reference PDFs.
No new Nature full text/SI was read. Scientific claims remain supported only
at their previously recorded access levels, not by metadata. Behara SI and
other remaining source limitations persist.

## Execution boundary

No new DFT, job submission/cancellation, production-source mutation, calculator
retuning, G2/holdout calculation or experiment-budget revision. The two live
G1 jobs checked at 11:28:59 CST remain the same 28319570/28319571 handles; this
receipt is not a later job-status observation. No pytest or recompilation is
claimed for this documentation-only correction. Physical inputs, prediction
protocols, code/tests and the nine-page manuscript remain unchanged. The JSON
receipt, 68 journal records and 32 local Markdown targets validate; the two
prediction protocols and manuscript/PDF retain their previous SHA256 hashes.
