"""Validate a VASP endpoint input before it is reused for VCNEB images."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.vasp import REQUIRED_VCNEB_STATIC_PARAMETERS, prepare_vasp_static_parameters
from vcneb.vasp_contract import canonical_parameters


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, help="directory containing INCAR, KPOINTS and POTCAR")
    parser.add_argument("--output", default=None, help="optional JSON report path")
    return parser.parse_args()


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(json.loads(canonical_parameters(payload)), handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def main() -> None:
    args = parse_args()
    source = Path(args.source).resolve()
    parameters, potcar = prepare_vasp_static_parameters(source)
    payload = {
        "status": "ok",
        "validation_scope": "calculator_free_not_bravais_certification",
        "input_contract_version": 1,
        "source_dir": str(source),
        "potcar": str(potcar.resolve()),
        "required_static_parameters": REQUIRED_VCNEB_STATIC_PARAMETERS,
        "effective_image_parameters": parameters,
    }
    if args.output:
        write_json_atomic(Path(args.output).resolve(), payload)
    print(canonical_parameters(payload))


if __name__ == "__main__":
    main()
