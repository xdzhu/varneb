"""Compare calculator-free initial cell paths for an endpoint pair."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import interpolate_vcneb, path_geometry_diagnostics, validate_atom_mapping


FIXTURE = ROOT / "validation" / "hfo2_t_to_po"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial", default=str(FIXTURE / "T_HfO2_12.vasp"))
    parser.add_argument("--final", default=str(FIXTURE / "PO_HfO2_12_mapped.vasp"))
    parser.add_argument("--n-images", type=int, default=7)
    parser.add_argument("--mapping", choices=["identity", "auto"], default="identity")
    parser.add_argument(
        "--output",
        default=str(ROOT / "outputs" / "hfo2_initial_path_comparison.json"),
    )
    parser.add_argument("--trajectory-directory", default=None)
    return parser.parse_args()


def summarize(images) -> dict:
    geometry = path_geometry_diagnostics(images)
    return {
        "volumes_A3": [record["volume_A3"] for record in geometry["images"]],
        "minimum_interatomic_distances_A": [
            record["minimum_interatomic_distance_A"] for record in geometry["images"]
        ],
        "deformation_norms": [
            record["deformation_from_reference_frobenius"] for record in geometry["images"]
        ],
        "minimum_interatomic_distance_over_path_A": min(
            geometry["images"][index]["minimum_interatomic_distance_A"]
            for index in range(len(geometry["images"]))
        ),
        "maximum_deformation_norm": max(
            record["deformation_from_reference_frobenius"] for record in geometry["images"]
        ),
    }


def main() -> None:
    args = parse_args()
    initial = read(args.initial)
    final = read(args.final)
    mapping_value = None if args.mapping == "identity" else "auto"
    mapping_report = validate_atom_mapping(initial, final, mapping_value, mic=True)
    result = {
        "initial": str(Path(args.initial).resolve()),
        "final": str(Path(args.final).resolve()),
        "n_images": args.n_images,
        "mapping_report": mapping_report,
        "paths": {},
    }
    trajectory_directory = None if args.trajectory_directory is None else Path(args.trajectory_directory)
    for strategy in ("linear", "log_strain"):
        images = interpolate_vcneb(
            initial,
            final,
            n_images=args.n_images,
            align_cells=True,
            mic=True,
            cell_interpolation=strategy,
            mapping=mapping_value,
        )
        result["paths"][strategy] = summarize(images)
        if trajectory_directory is not None:
            trajectory_directory.mkdir(parents=True, exist_ok=True)
            write(trajectory_directory / f"initial-{strategy}.traj", images)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for strategy, summary in result["paths"].items():
        print(
            f"{strategy}: min_distance_A={summary['minimum_interatomic_distance_over_path_A']:.6f} "
            f"max_deformation={summary['maximum_deformation_norm']:.6f}"
        )
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()
