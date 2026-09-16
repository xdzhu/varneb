"""Refuse a production VCNEB launch unless its endpoint identity gate passed."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gate", required=True, help="endpoint_identity_gate.json from the comparison script")
    args = parser.parse_args()
    path = Path(args.gate).resolve()
    payload = json.loads(path.read_text(encoding="utf-8"))
    endpoints = payload.get("endpoints")
    if payload.get("matches") is not True or not isinstance(endpoints, dict):
        raise SystemExit(f"[BLOCKED] endpoint identity gate did not pass: {path}")
    for label in ("initial", "final"):
        record = endpoints.get(label)
        if not isinstance(record, dict) or record.get("matches") is not True:
            raise SystemExit(f"[BLOCKED] {label} endpoint identity did not pass: {path}")
    print(f"[OK] endpoint identity gate passed: {path}")


if __name__ == "__main__":
    main()
