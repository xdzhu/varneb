"""Configuration and calculator-free preparation for VARNEB runs."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Mapping

from ase.io import read, write

from .backends import get_backend_spec
from .core import interpolate_vcneb, path_geometry_diagnostics
from .optimizer_registry import get_optimizer_spec
from .provenance import endpoint_structure_record


@dataclass(frozen=True)
class RunConfig:
    """Resolved, calculator-independent run configuration."""

    backend: str
    initial: Path
    final: Path
    workdir: Path
    n_images: int = 7
    fmax_ev_per_angstrom: float = 0.10
    k: float = 0.20
    pressure_gpa: float = 0.0
    cell_interpolation: str = "log_strain"
    mapping: str = "auto"
    mic: bool = True
    align_translation: bool = True
    minimum_distance: float | None = None
    maximum_deformation: float | None = None
    climb: bool = False
    climb_after: int | None = None
    optimizer: str = "FIRE"
    steps: int = 300
    calculator: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        get_backend_spec(self.backend)
        get_optimizer_spec(self.optimizer)
        if self.n_images < 3:
            raise ValueError("n_images must include two endpoints and one interior image")
        if not 0.0 < self.fmax_ev_per_angstrom:
            raise ValueError("fmax_ev_per_angstrom must be positive")
        if self.cell_interpolation not in {"linear", "log_strain"}:
            raise ValueError("cell_interpolation must be 'linear' or 'log_strain'")
        if self.mapping not in {"identity", "auto"}:
            raise ValueError("mapping must be 'identity' or 'auto'")
        if self.steps < 1:
            raise ValueError("steps must be positive")
        if self.minimum_distance is not None and self.minimum_distance <= 0:
            raise ValueError("minimum_distance must be positive when provided")
        if self.maximum_deformation is not None and self.maximum_deformation <= 0:
            raise ValueError("maximum_deformation must be positive when provided")

    @classmethod
    def from_file(cls, path: str | Path) -> "RunConfig":
        source = Path(path).expanduser().resolve()
        data = json.loads(source.read_text(encoding="utf-8"))
        if data.get("schema_version") != 1:
            raise ValueError("config schema_version must be 1")
        base = source.parent
        path_fields = {}
        for key in ("initial", "final", "workdir"):
            value = data.get(key)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"config field {key!r} must be a non-empty path")
            path_fields[key] = (base / value).resolve()
        calculator = data.get("calculator", {})
        if not isinstance(calculator, dict):
            raise ValueError("calculator must be a mapping")
        if not isinstance(calculator.get("parameters", {}), dict):
            raise ValueError("calculator.parameters must be a mapping")
        values = {
            "backend": str(data.get("backend", "")),
            **path_fields,
            "n_images": int(data.get("n_images", 7)),
            "fmax_ev_per_angstrom": float(data.get("fmax_ev_per_angstrom", 0.10)),
            "k": float(data.get("k", 0.20)),
            "pressure_gpa": float(data.get("pressure_gpa", 0.0)),
            "cell_interpolation": str(data.get("cell_interpolation", "log_strain")),
            "mapping": str(data.get("mapping", "auto")),
            "mic": bool(data.get("mic", True)),
            "align_translation": bool(data.get("align_translation", True)),
            "minimum_distance": (
                None if data.get("minimum_distance") is None
                else float(data["minimum_distance"])
            ),
            "maximum_deformation": (
                None if data.get("maximum_deformation") is None
                else float(data["maximum_deformation"])
            ),
            "climb": bool(data.get("climb", False)),
            "climb_after": (
                None if data.get("climb_after") is None else int(data["climb_after"])
            ),
            "optimizer": str(data.get("optimizer", "FIRE")),
            "steps": int(data.get("steps", 300)),
            "calculator": calculator,
        }
        return cls(**values)

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        for key in ("initial", "final", "workdir"):
            result[key] = str(result[key])
        return result


def prepare_run(path: str | Path) -> tuple[RunConfig, Path]:
    """Build and persist a calculator-free initial path from a JSON config."""

    config = RunConfig.from_file(path)
    initial, final = read(config.initial), read(config.final)
    initial_symbols = initial.get_chemical_symbols()
    final_symbols = final.get_chemical_symbols()
    if config.mapping == "identity" and initial_symbols != final_symbols:
        raise ValueError("identity mapping requires identical endpoint atom order/species")
    if sorted(initial_symbols) != sorted(final_symbols):
        raise ValueError("endpoint compositions differ; automatic mapping cannot reconcile them")
    images = interpolate_vcneb(
        initial,
        final,
        config.n_images,
        align_cells=True,
        mic=config.mic,
        cell_interpolation=config.cell_interpolation,
        mapping=None if config.mapping == "identity" else "auto",
        align_translation=config.align_translation,
        minimum_distance=config.minimum_distance,
        maximum_deformation=config.maximum_deformation,
    )
    config.workdir.mkdir(parents=True, exist_ok=True)
    trajectory = config.workdir / "initial-vcneb.traj"
    write(trajectory, images)
    report = {
        "status": "prepared",
        "calculator_attached": False,
        "config": config.to_dict(),
        "endpoint_structures": {
            "initial": endpoint_structure_record(initial),
            "final": endpoint_structure_record(final),
        },
        "initial_path_geometry": path_geometry_diagnostics(
            images,
            minimum_distance=config.minimum_distance,
            maximum_deformation=config.maximum_deformation,
        ),
        "trajectory": str(trajectory),
    }
    report_path = config.workdir / "varneb_preflight.json"
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(report_path)
    return config, report_path


__all__ = ["RunConfig", "prepare_run"]
