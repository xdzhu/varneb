"""Backward-compatible script entry point for the packaged VARNEB runner."""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import material_runner as _runner  # noqa: E402

# Keep historical helper imports and test-time calculator injection working.
_load_symbol = _runner._load_symbol
_cell_mask_for_mode = _runner._cell_mask_for_mode
_validate_static_endpoint_identity = _runner._validate_static_endpoint_identity
_climb_enabled = _runner._climb_enabled
parse_args = _runner.parse_args


def main(argv: list[str] | None = None) -> None:
    _runner.main(argv, symbol_loader=_load_symbol)


if __name__ == "__main__":
    main()
