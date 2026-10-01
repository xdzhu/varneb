# CPC source-snapshot candidate

The installable PyPI sdist is intentionally small; it is **not** the complete
source/evidence set behind this manuscript. Use this procedure to prepare a
reviewable candidate from one fixed Git commit. It does not upload to CPC,
Mendeley Data, GitHub Releases, or PyPI.

From the repository root, after choosing the exact reviewed commit:

```powershell
python -m scripts.build_cpc_source_snapshot --revision <full-commit-sha> --output tmp/varneb-cpc-source.zip
```

The command refuses an existing output. It creates the ZIP and an adjacent
`varneb-cpc-source.manifest.json` with the resolved commit, ZIP SHA-256,
entry count, path-inventory SHA-256, and required-file gate. It compares every
ZIP entry name with the pinned Git tree; uncommitted and untracked local files are
not included. The required set covers the installable package, tests, user
manual, runnable toy, manuscript, four main landscape source-data/QA bundles,
selected material audits, and the final GaN chain. The complete tracked tree
is included, not just these gate files. Git export attributes may transform
ordinary text line endings; the raw-byte scientific files have explicit
`.gitattributes` rules and are independently checked by figure regressions.

Common potential/orbital and secret-like filenames cause a preflight refusal.
This is only a narrow automatic gate: a human must still inspect the complete
entry list and licences before distributing the ZIP. Do **not** add VASP
POTCARs, ABACUS orbitals/pseudopotentials, other licensed potential files,
DFT executables, account credentials, or cluster-private data to the public
bundle. Input contracts record identities and hashes; recipients acquire any
licensed computational dependencies separately.

Recommended independent acceptance check after extracting the candidate in
a clean directory with the documented Python/TeX dependencies:

```powershell
python -m pytest -q
python -m scripts.rebuild_cpc_landscapes --output-dir tmp/landscape-rebuild --strict-png
cd paper/VARNEB_CPC
latexmk -pdf -interaction=nonstopmode -halt-on-error varneb_CPC.tex
```

Record the test/skipped counts and rendered-page review against the **same
commit** as the manifest. The strict PNG check is environment-sensitive;
source CSV identity and measured-node counts are the scientific figure gates.
A local build passing these checks is still a candidate, not an approved CPC
Program Library accession or proof that a third-party DFT executable works.
The author controls final metadata, licence review, deposit, and submission.
