"""Create a VASP BTO static/VCNEB input directory from approved endpoint data.

The script never creates a POTCAR.  It copies an explicitly supplied,
user-licensed file and proves that the copied T/C endpoints match the accepted
ABACUS endpoint identity before any VASP calculation is requested.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile

from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.provenance import compare_endpoint_records, endpoint_structure_record


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    template = ROOT / "examples" / "batio3_vasp_pbe_paw_static"
    parser.add_argument("--initial", required=True, help="accepted tetragonal endpoint structure")
    parser.add_argument("--final", required=True, help="accepted cubic endpoint structure")
    parser.add_argument("--reference-identity", required=True, help="accepted ABACUS endpoint identity JSON")
    parser.add_argument("--potcar", required=True, help="explicit user-licensed Ba/Ti/O PBE PAW POTCAR")
    parser.add_argument("--destination", required=True, help="new VASP case directory")
    parser.add_argument("--incar", default=str(template / "INCAR"))
    parser.add_argument("--kpoints", default=str(template / "KPOINTS"))
    return parser.parse_args()


def setup(
    *,
    initial_path: Path,
    final_path: Path,
    reference_identity: Path,
    potcar: Path,
    destination: Path,
    incar: Path,
    kpoints: Path,
) -> dict:
    for path in (initial_path, final_path, reference_identity, potcar, incar, kpoints):
        if not path.is_file():
            raise FileNotFoundError(path)
    initial, final = read(initial_path), read(final_path)
    reference = json.loads(reference_identity.read_text(encoding="utf-8")).get("endpoint_structures")
    if not isinstance(reference, dict):
        raise ValueError("reference identity lacks endpoint_structures")
    endpoint_records = {"initial": endpoint_structure_record(initial), "final": endpoint_structure_record(final)}
    comparison = compare_endpoint_records(reference, endpoint_records)
    if comparison.get("matches") is not True:
        raise ValueError("supplied VASP endpoints do not match the accepted ABACUS endpoint identity")
    initial_dir, final_dir = destination / "initial", destination / "final"
    if destination.exists():
        raise FileExistsError(f"destination already exists: {destination}")
    initial_dir.mkdir(parents=True)
    final_dir.mkdir(parents=True)
    write(initial_dir / "CONTCAR", initial, format="vasp", direct=True, vasp5=True)
    write(final_dir / "CONTCAR", final, format="vasp", direct=True, vasp5=True)
    shutil.copy2(incar, initial_dir / "INCAR")
    shutil.copy2(kpoints, initial_dir / "KPOINTS")
    shutil.copy2(potcar, initial_dir / "POTCAR")
    manifest = {
        "status": "prepared_no_dft",
        "destination": str(destination.resolve()),
        "endpoint_identity": endpoint_records,
        "endpoint_comparison": comparison,
        "inputs": {
            "initial_source": str(initial_path.resolve()),
            "final_source": str(final_path.resolve()),
            "reference_identity": str(reference_identity.resolve()),
            "potcar_source": str(potcar.resolve()),
            "potcar_sha256": _sha256(potcar),
            "incar_sha256": _sha256(incar),
            "kpoints_sha256": _sha256(kpoints),
        },
    }
    _write_json_atomic(destination / "setup_manifest.json", manifest)
    return manifest


def main() -> int:
    args = parse_args()
    try:
        result = setup(
            initial_path=Path(args.initial).resolve(),
            final_path=Path(args.final).resolve(),
            reference_identity=Path(args.reference_identity).resolve(),
            potcar=Path(args.potcar).resolve(),
            destination=Path(args.destination).resolve(),
            incar=Path(args.incar).resolve(),
            kpoints=Path(args.kpoints).resolve(),
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
