"""Create a calculator-independent VCNEB endpoint identity reference from structures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import tempfile

from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.provenance import endpoint_structure_record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial", required=True, help="ASE-readable initial endpoint")
    parser.add_argument("--final", required=True, help="ASE-readable final endpoint")
    parser.add_argument("--output", default="vcneb_endpoint_identity.json")
    parser.add_argument("--label", default="accepted VCNEB endpoint identity")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    initial_path, final_path = Path(args.initial).resolve(), Path(args.final).resolve()
    initial, final = read(initial_path), read(final_path)
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "format_version": 1,
        "status": "accepted_reference_identity",
        "label": args.label,
        "initial_path": str(initial_path),
        "final_path": str(final_path),
        "endpoint_structures": {
            "initial": endpoint_structure_record(initial),
            "final": endpoint_structure_record(final),
        },
    }
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=output.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(output)
    print(f"[OK] endpoint identity reference={output}")


if __name__ == "__main__":
    main()
