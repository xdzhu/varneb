"""Derive a reproducible endpoint-displacement mode for a VC-NEB path.

The exported mode is an *endpoint displacement diagnostic*: it is useful for
testing mode-guided and constrained workflows, but it is not a phonon
eigenvector.  The final endpoint is first aligned/mapped with the same
``interpolate_vcneb`` semantics used by the production drivers, so the JSON
mode is expressed in the initial endpoint atom order and cell gauge.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import interpolate_vcneb
from vcneb.core import state_from_atoms


def derive_endpoint_mode(
    initial_path: str | Path,
    final_path: str | Path,
    *,
    mapping: str = "auto",
    mic: bool = True,
    align_translation: bool = True,
) -> dict:
    """Return a JSON-serializable endpoint displacement mode payload."""

    initial = read(initial_path)
    final = read(final_path)
    images = interpolate_vcneb(
        initial,
        final,
        n_images=2,
        align_cells=True,
        mic=mic,
        mapping=None if mapping == "identity" else mapping,
        align_translation=align_translation,
    )
    first, last = images
    reference_cell = first.cell.array
    first_state = state_from_atoms(first, reference_cell)
    last_state = state_from_atoms(last, reference_cell)
    payload = {
        "label": "endpoint-displacement-diagnostic",
        "semantics": "endpoint displacement diagnostic; not a phonon eigenvector",
        "n_atoms": len(first),
        "symbols": first.get_chemical_symbols(),
        "mapping": first.info.get("vcneb_path_metadata", {}).get("mapping"),
        "translation_fractional_final_cell": first.info.get("vcneb_path_metadata", {}).get(
            "translation_fractional_final_cell"
        ),
        "atomic": (
            last_state.q @ reference_cell - first_state.q @ reference_cell
        ).tolist(),
        "cell": (last_state.deform - first_state.deform).tolist(),
        "source": {
            "initial": str(Path(initial_path).resolve()),
            "final": str(Path(final_path).resolve()),
            "mapping": mapping,
            "mic": bool(mic),
            "align_translation": bool(align_translation),
            "cell_alignment": "polar_rotation_removed",
        },
    }
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial", required=True, help="ASE-readable initial endpoint")
    parser.add_argument("--final", required=True, help="ASE-readable final endpoint")
    parser.add_argument("--output", required=True, help="Output JSON mode path")
    parser.add_argument("--mapping", choices=["identity", "auto"], default="auto")
    parser.add_argument("--no-mic", action="store_true")
    parser.add_argument("--no-align-translation", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = derive_endpoint_mode(
        args.initial,
        args.final,
        mapping=args.mapping,
        mic=not args.no_mic,
        align_translation=not args.no_align_translation,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[OK] wrote endpoint displacement mode: {output.resolve()}")


if __name__ == "__main__":
    main()
