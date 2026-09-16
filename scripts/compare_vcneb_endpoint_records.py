"""Gate a candidate calculator preflight against an ABACUS endpoint identity."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.provenance import compare_endpoint_records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-preflight", required=True, help="Accepted ABACUS vcneb_preflight.json")
    parser.add_argument("--candidate-preflight", required=True, help="QE or VASP vcneb_preflight.json")
    parser.add_argument("--output", default="endpoint_identity_gate.json")
    return parser.parse_args()


def _load(path: str) -> dict:
    source = Path(path).resolve()
    payload = json.loads(source.read_text(encoding="utf-8"))
    records = payload.get("endpoint_structures")
    if not isinstance(records, dict):
        raise ValueError(f"{source} has no endpoint_structures record; regenerate its no-DFT preflight")
    return records


def main() -> None:
    args = parse_args()
    result = compare_endpoint_records(_load(args.reference_preflight), _load(args.candidate_preflight))
    result.update(
        {
            "format_version": 1,
            "reference_preflight": str(Path(args.reference_preflight).resolve()),
            "candidate_preflight": str(Path(args.candidate_preflight).resolve()),
        }
    )
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if not result["matches"]:
        raise SystemExit(f"[BLOCKED] endpoint identities differ; report={output}")
    print(f"[OK] endpoint identities match; report={output}")


if __name__ == "__main__":
    main()
