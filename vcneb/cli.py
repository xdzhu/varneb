"""Command-line entry point for VARNEB.

The import package remains ``vcneb`` for source compatibility; ``varneb`` is
the public distribution/project name.  The legacy ``vcneb`` console alias is
kept so existing scripts continue to work.
"""

from __future__ import annotations

import argparse

from .version import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="varneb",
        description="VARNEB: calculator-agnostic variable-cell nudged elastic band tools.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def main() -> int:
    build_parser().parse_args()
    return 0
