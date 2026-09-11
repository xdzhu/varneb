# Provenance: USPEX VCNEB and mode-guided path analysis

Date: 2026-09-12

## Local archive inspected

- Input archive: `H:\ReSearch\VCNEB\USPEX\USPEX_v10.6.tar.gz`
- Extracted copy: `H:\ReSearch\VCNEB\USPEX\extracted_10.6\USPEX_v10.6\`
- Observed contents: `install.sh`, `README`, `MCR_license.txt`, `USPEX_license.txt`, and `USPEX_MATLABruntime.install`
- The last file is a stripped x86-64 Linux ELF self-extracting installer. No readable USPEX source files were found in the archive.
- The archive was not installed or executed; no license acceptance was performed by this analysis.

## Primary web sources

1. Qian et al. VC-NEB paper: https://doi.org/10.1016/j.cpc.2013.04.004
2. Public paper PDF: https://uspex-team.org/static/file/Qian-vcNEB-2013.pdf
3. USPEX VCNEB manual: https://uspex-team.org/online_utilities/uspex_manual_release/EnglishVersion/uspex_manual_english/vcneb.html
4. ABINIT GeoConstraints: https://docs.abinit.org/topics/GeoConstraints/
5. ABINIT relaxation variables: https://docs.abinit.org/variables/rlx/
6. ABINIT transition paths: https://docs.abinit.org/topics/TransPath/
7. Modal NEB analysis repository: https://github.com/kgordiz/modalNEB
8. Modal phonon-ion analysis preprint: https://arxiv.org/abs/2305.01632

## Evidence classification

- `USPEX` algorithm options and workflow: derived from the public Qian paper and official VCNEB manual.
- USPEX package format and license restrictions: derived from the local `README` and `USPEX_license.txt`.
- ABINIT axis/component and linear-combination constraints: derived from the current official ABINIT documentation.
- Mode-guided path and modal projection code: original implementation in `vcneb/modes.py`; regression-tested on `235` and `cu05`.
- No USPEX executable code, decompiled code, or proprietary implementation was copied into this repository.
