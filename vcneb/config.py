"""Configuration and calculator-free preparation for VARNEB runs."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
from typing import Mapping

import numpy as np
from ase.io import read, write

from .backends import get_backend_spec
from .core import interpolate_vcneb, path_geometry_diagnostics
from .optimizer_registry import get_optimizer_spec
from .provenance import endpoint_structure_record


def _text(value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"config field {field_name!r} must be a non-empty string")
    return value


def _integer(value: object, field_name: str, *, nullable: bool = False) -> int | None:
    if value is None and nullable:
        return None
    if type(value) is not int:
        raise ValueError(f"config field {field_name!r} must be an integer")
    return value


def _number(value: object, field_name: str, *, nullable: bool = False) -> float | None:
    if value is None and nullable:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"config field {field_name!r} must be a number")
    try:
        return float(value)
    except OverflowError as exc:
        raise ValueError(f"config field {field_name!r} exceeds the numeric range") from exc


def _boolean(value: object, field_name: str) -> bool:
    if type(value) is not bool:
        raise ValueError(f"config field {field_name!r} must be a boolean")
    return value


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
    cell_mode: str = "full"
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
    image_workers: int = 0
    image_retries: int = 0
    candidate_step_retries: int = 0
    maxstep: float | None = None
    maximum_cell_step: float | None = None
    minimum_endpoint_separation: float | None = None
    endpoint_static_summary: Path | None = None
    calculator: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        get_backend_spec(self.backend)
        get_optimizer_spec(self.optimizer)
        if self.n_images < 3:
            raise ValueError("n_images must include two endpoints and one interior image")
        if not math.isfinite(self.fmax_ev_per_angstrom) or self.fmax_ev_per_angstrom <= 0:
            raise ValueError("fmax_ev_per_angstrom must be finite and positive")
        if not math.isfinite(self.k) or self.k <= 0:
            raise ValueError("k must be finite and positive")
        if not math.isfinite(self.pressure_gpa):
            raise ValueError("pressure_gpa must be finite")
        if self.cell_mode not in {"full", "fixed"}:
            raise ValueError("cell_mode must be 'full' or 'fixed'")
        if self.cell_mode == "fixed" and self.pressure_gpa != 0.0:
            raise ValueError("fixed-cell NEB requires zero pressure_gpa")
        if self.cell_interpolation not in {"linear", "log_strain"}:
            raise ValueError("cell_interpolation must be 'linear' or 'log_strain'")
        if self.mapping not in {"identity", "auto"}:
            raise ValueError("mapping must be 'identity' or 'auto'")
        if self.steps < 1:
            raise ValueError("steps must be positive")
        if any(value < 0 for value in (
            self.image_workers, self.image_retries, self.candidate_step_retries
        )):
            raise ValueError("image worker and retry counts must be non-negative")
        if self.climb_after is not None and self.climb_after < 0:
            raise ValueError("climb_after must be non-negative")
        for name in ("maxstep", "maximum_cell_step", "minimum_endpoint_separation"):
            value = getattr(self, name)
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError(f"{name} must be finite and positive when provided")
        for name in ("minimum_distance", "maximum_deformation"):
            value = getattr(self, name)
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError(f"{name} must be finite and positive when provided")

    @classmethod
    def from_file(cls, path: str | Path) -> "RunConfig":
        source = Path(path).expanduser().resolve()
        data = json.loads(source.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError("VARNEB config must be a JSON object")
        if type(data.get("schema_version")) is not int or data["schema_version"] != 1:
            raise ValueError("config schema_version must be 1")
        unknown = set(data) - set(cls.__dataclass_fields__) - {"schema_version"}
        if unknown:
            raise ValueError(f"unknown config field(s): {', '.join(sorted(unknown))}")
        base = source.parent
        path_fields = {}
        for key in ("initial", "final", "workdir"):
            value = _text(data.get(key), key)
            path_fields[key] = (base / value).resolve()
        calculator = data.get("calculator", {})
        if not isinstance(calculator, dict):
            raise ValueError("calculator must be a mapping")
        extra_calculator = set(calculator) - {
            "kind", "symbol", "parameters", "command", "factory_kwargs"
        }
        if extra_calculator:
            raise ValueError(
                f"unknown calculator field(s): {', '.join(sorted(extra_calculator))}"
            )
        for key in ("kind", "symbol"):
            if calculator.get(key) is not None:
                _text(calculator[key], f"calculator.{key}")
        if calculator.get("command") is not None and not isinstance(calculator["command"], str):
            raise ValueError("calculator.command must be a string")
        if not isinstance(calculator.get("parameters", {}), dict):
            raise ValueError("calculator.parameters must be a mapping")
        if not isinstance(calculator.get("factory_kwargs", {}), dict):
            raise ValueError("calculator.factory_kwargs must be a mapping")
        values = {
            "backend": _text(data.get("backend"), "backend"),
            **path_fields,
            "n_images": _integer(data.get("n_images", 7), "n_images"),
            "fmax_ev_per_angstrom": _number(
                data.get("fmax_ev_per_angstrom", 0.10), "fmax_ev_per_angstrom"
            ),
            "k": _number(data.get("k", 0.20), "k"),
            "pressure_gpa": _number(data.get("pressure_gpa", 0.0), "pressure_gpa"),
            "cell_mode": _text(data.get("cell_mode", "full"), "cell_mode"),
            "cell_interpolation": _text(
                data.get("cell_interpolation", "log_strain"), "cell_interpolation"
            ),
            "mapping": _text(data.get("mapping", "auto"), "mapping"),
            "mic": _boolean(data.get("mic", True), "mic"),
            "align_translation": _boolean(
                data.get("align_translation", True), "align_translation"
            ),
            "minimum_distance": _number(
                data.get("minimum_distance"), "minimum_distance", nullable=True
            ),
            "maximum_deformation": _number(
                data.get("maximum_deformation"), "maximum_deformation", nullable=True
            ),
            "climb": _boolean(data.get("climb", False), "climb"),
            "climb_after": _integer(data.get("climb_after"), "climb_after", nullable=True),
            "optimizer": _text(data.get("optimizer", "FIRE"), "optimizer"),
            "steps": _integer(data.get("steps", 300), "steps"),
            "image_workers": _integer(data.get("image_workers", 0), "image_workers"),
            "image_retries": _integer(data.get("image_retries", 0), "image_retries"),
            "candidate_step_retries": _integer(
                data.get("candidate_step_retries", 0), "candidate_step_retries"
            ),
            "maxstep": _number(data.get("maxstep"), "maxstep", nullable=True),
            "maximum_cell_step": _number(
                data.get("maximum_cell_step"), "maximum_cell_step", nullable=True
            ),
            "minimum_endpoint_separation": _number(
                data.get("minimum_endpoint_separation"),
                "minimum_endpoint_separation", nullable=True
            ),
            "endpoint_static_summary": (
                None if data.get("endpoint_static_summary") is None
                else (base / _text(
                    data["endpoint_static_summary"], "endpoint_static_summary"
                )).resolve()
            ),
            "calculator": calculator,
        }
        return cls(**values)

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        for key in ("initial", "final", "workdir"):
            result[key] = str(result[key])
        if result["endpoint_static_summary"] is not None:
            result["endpoint_static_summary"] = str(result["endpoint_static_summary"])
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
    if config.cell_mode == "fixed" and any(
        not np.allclose(image.cell.array, images[0].cell.array, atol=1e-10, rtol=0.0)
        for image in images
    ):
        raise ValueError("fixed-cell NEB requires identical cells for all images")
    config.workdir.mkdir(parents=True, exist_ok=True)
    trajectory = config.workdir / "initial-vcneb.traj"
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
    if trajectory.exists() or report_path.exists():
        if not trajectory.is_file() or not report_path.is_file():
            raise FileExistsError("prepared workdir is incomplete; refusing to overwrite it")
        previous = json.loads(report_path.read_text(encoding="utf-8"))
        saved_images = read(trajectory, index=":")
        same_chain = len(saved_images) == len(images) and all(
            old.get_chemical_symbols() == new.get_chemical_symbols()
            and np.array_equal(old.pbc, new.pbc)
            and np.allclose(old.cell.array, new.cell.array, atol=1e-12, rtol=0.0)
            and np.allclose(old.positions, new.positions, atol=1e-12, rtol=0.0)
            for old, new in zip(saved_images, images)
        )
        if previous != report or not same_chain:
            raise FileExistsError("prepared path or preflight differs; choose a new workdir")
        return config, report_path
    write(trajectory, images)
    temporary = report_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(report_path)
    return config, report_path


__all__ = ["RunConfig", "prepare_run"]
