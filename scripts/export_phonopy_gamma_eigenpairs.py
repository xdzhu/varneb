"""Export Phonopy-generated, phase-fixed Gamma eigenvectors without running DFT."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

from phonopy import load

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.phonons import phonopy_gamma_eigenpairs, save_phonopy_gamma_eigenpairs


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phonopy-yaml", default="phonopy_disp.yaml")
    parser.add_argument("--force-sets", default="FORCE_SETS")
    parser.add_argument("--output", default="phonopy_gamma_eigenpairs.npz")
    parser.add_argument("--provenance", default="phonopy_gamma_eigenpairs_provenance.json")
    parser.add_argument("--calculator", default="abacus")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    yaml_path = Path(args.phonopy_yaml).resolve()
    force_sets_path = Path(args.force_sets).resolve()
    if not yaml_path.is_file() or not force_sets_path.is_file():
        raise FileNotFoundError("--phonopy-yaml and --force-sets must be existing files")
    phonon = load(
        str(yaml_path),
        force_sets_filename=str(force_sets_path),
        calculator=args.calculator,
        is_compact_fc=False,
    )
    eigenpairs = phonopy_gamma_eigenpairs(phonon)
    output = Path(args.output).resolve()
    save_phonopy_gamma_eigenpairs(output, eigenpairs)
    provenance = {
        "format_version": 1,
        "method": "Phonopy run_qpoints([[0, 0, 0]], with_eigenvectors=True)",
        "eigenvector_convention": "mass-weighted; arbitrary U(1) phase fixed at largest component",
        "frequencies_unit": "THz; negative values denote imaginary harmonic modes",
        "phonopy_yaml": str(yaml_path),
        "phonopy_yaml_sha256": _sha256(yaml_path),
        "force_sets": str(force_sets_path),
        "force_sets_sha256": _sha256(force_sets_path),
        "output": str(output),
        "n_modes": int(len(eigenpairs.frequencies_thz)),
    }
    provenance_path = Path(args.provenance).resolve()
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[DONE] exported {provenance['n_modes']} Phonopy Gamma eigenpairs to {output}")


if __name__ == "__main__":
    main()
