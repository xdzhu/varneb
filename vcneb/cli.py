"""Small command-line entry point for the installable VCNEB package."""

from __future__ import annotations

import argparse

from .version import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vcneb",
        description="Calculator-agnostic variable-cell nudged elastic band tools.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def main() -> int:
    build_parser().parse_args()
    return 0
