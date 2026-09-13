"""Image-level calculator execution backends.

The VC-NEB controller remains a single Python process while an executor can
evaluate independent image calculators concurrently.  In a Slurm allocation,
the calculator command should use ``srun --exclusive`` so each concurrent
call receives an isolated job step and MPI ranks.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from typing import Sequence

import numpy as np
from ase import Atoms


@dataclass(frozen=True)
class ImageEvaluation:
    """Energy, Cartesian forces and full stress for one image."""

    energy: float
    forces: np.ndarray
    stress: np.ndarray

    def validate(self, *, image_index: int, n_atoms: int) -> "ImageEvaluation":
        forces = np.asarray(self.forces, dtype=float)
        stress = np.asarray(self.stress, dtype=float)
        if forces.shape != (n_atoms, 3):
            raise ValueError(
                f"image {image_index} force shape {forces.shape} != {(n_atoms, 3)}"
            )
        if stress.shape != (3, 3):
            raise ValueError(f"image {image_index} stress shape {stress.shape} != (3, 3)")
        if not np.isfinite(float(self.energy)) or not np.all(np.isfinite(forces)) or not np.all(np.isfinite(stress)):
            raise ValueError(f"image {image_index} evaluation contains NaN/Inf")
        return ImageEvaluation(float(self.energy), forces, stress)


class ThreadedCalculatorExecutor:
    """Evaluate independent ASE calculators concurrently.

    This backend intentionally uses threads rather than MPI.  ASE calculator
    calls block while the external ABACUS/VASP process runs, so the Python GIL
    does not serialize the expensive part.  Each image must have a unique
    calculator directory.  The external command is responsible for allocating
    ranks; for Slurm use ``srun --exclusive`` in that command.
    """

    def __init__(
        self,
        max_workers: int,
        *,
        max_retries: int = 0,
        manifest_path: str | Path | None = None,
        cache_dir: str | Path | None = None,
        cache_namespace: str | None = None,
    ):
        if isinstance(max_workers, bool) or int(max_workers) != max_workers or max_workers < 1:
            raise ValueError("max_workers must be a positive integer")
        if isinstance(max_retries, bool) or int(max_retries) != max_retries or max_retries < 0:
            raise ValueError("max_retries must be a non-negative integer")
        self.max_workers = int(max_workers)
        self.max_retries = int(max_retries)
        self.manifest_path = Path(manifest_path) if manifest_path is not None else None
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None
        self.cache_namespace = None if cache_namespace is None else str(cache_namespace)
        if self.cache_dir is not None:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            self._validate_cache_namespace()
        self.last_attempts: dict[int, int] = {}

    def _validate_cache_namespace(self) -> None:
        """Reject accidental reuse of a cache with a different parameter set."""

        assert self.cache_dir is not None
        metadata_path = self.cache_dir / "cache_metadata.json"
        expected = {
            "format_version": 1,
            "namespace": self.cache_namespace,
        }
        if metadata_path.exists():
            try:
                actual = json.loads(metadata_path.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                raise RuntimeError(f"cannot read image cache metadata: {metadata_path}") from exc
            if actual != expected:
                raise ValueError(
                    "image cache metadata does not match the requested namespace; "
                    "use a new cache_dir or the original cache_namespace"
                )
            return
        temporary = metadata_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(expected, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(metadata_path)

    @staticmethod
    def _image_cache_digest(image_index: int, image: Atoms) -> str:
        """Return a byte-stable key for one image state.

        The cache is intentionally exact: even tiny coordinate/cell changes must
        trigger a fresh calculator call.  The caller supplies a namespace through
        ``cache_namespace`` to distinguish calculator parameter sets.
        """

        digest = hashlib.sha256()
        digest.update(str(int(image_index)).encode("ascii"))
        digest.update(b"\0")
        digest.update("\0".join(image.get_chemical_symbols()).encode("utf-8"))
        for value in (
            np.asarray(image.get_positions(), dtype="<f8", order="C"),
            np.asarray(image.cell.array, dtype="<f8", order="C"),
            np.asarray(image.pbc, dtype="<?", order="C"),
        ):
            digest.update(np.asarray(value).tobytes(order="C"))
            digest.update(b"\0")
        return digest.hexdigest()

    def _cache_path(self, image_index: int, image: Atoms) -> tuple[str, Path] | None:
        if self.cache_dir is None:
            return None
        key = self._image_cache_digest(image_index, image)
        return key, self.cache_dir / f"image_{int(image_index):04d}_{key}.npz"

    def _load_cached(self, image_index: int, image: Atoms) -> ImageEvaluation | None:
        cache_path = self._cache_path(image_index, image)
        if cache_path is None:
            return None
        _, path = cache_path
        if not path.exists():
            return None
        try:
            with np.load(path, allow_pickle=False) as data:
                evaluation = ImageEvaluation(
                    float(np.asarray(data["energy"]).reshape(())),
                    np.asarray(data["forces"], dtype=float),
                    np.asarray(data["stress"], dtype=float),
                )
            return evaluation.validate(image_index=image_index, n_atoms=len(image))
        except (OSError, KeyError, ValueError, TypeError, EOFError):
            # A partial/corrupt cache entry is treated as a miss.  The original
            # file is retained for post-mortem inspection and will be replaced
            # atomically after a successful calculator call.
            return None

    def _store_cached(self, image_index: int, image: Atoms, evaluation: ImageEvaluation) -> None:
        cache_path = self._cache_path(image_index, image)
        if cache_path is None:
            return
        _, path = cache_path
        temporary = path.with_suffix(path.suffix + ".tmp")
        with temporary.open("wb") as handle:
            np.savez_compressed(
                handle,
                energy=np.asarray(evaluation.energy, dtype=float),
                forces=np.asarray(evaluation.forces, dtype=float),
                stress=np.asarray(evaluation.stress, dtype=float),
            )
            handle.flush()
        temporary.replace(path)

    def _write_manifest_record(self, record: dict) -> None:
        """Append one durable controller-iteration record to JSONL."""
        if self.manifest_path is None:
            return
        self.manifest_path.parent.mkdir(parents=True, exist_ok=True)
        with self.manifest_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
            handle.flush()

    @staticmethod
    def _evaluate_one(image_index: int, image: Atoms) -> ImageEvaluation:
        try:
            energy = float(image.get_potential_energy())
            forces = np.asarray(image.get_forces(), dtype=float)
            stress = np.asarray(image.get_stress(voigt=False), dtype=float)
            return ImageEvaluation(energy, forces, stress).validate(
                image_index=image_index, n_atoms=len(image)
            )
        except Exception as exc:
            raise RuntimeError(f"parallel calculator evaluation failed for image {image_index}: {exc}") from exc

    def evaluate(
        self,
        images: Sequence[Atoms],
        *,
        indices: Sequence[int] | None = None,
    ) -> list[ImageEvaluation]:
        """Evaluate a batch, optionally identified by chain ``indices``."""
        if not images:
            return []
        if indices is None:
            image_indices = list(range(len(images)))
        else:
            image_indices = [int(index) for index in indices]
            if len(image_indices) != len(images):
                raise ValueError("indices must have the same length as images")
            if len(set(image_indices)) != len(image_indices):
                raise ValueError("indices must be unique")
        self.last_attempts = {}
        started = time.monotonic()
        cached: dict[int, ImageEvaluation] = {}
        missing_images: list[tuple[int, Atoms]] = []
        for image_index, image in zip(image_indices, images):
            result = self._load_cached(image_index, image)
            if result is None:
                missing_images.append((image_index, image))
            else:
                cached[image_index] = result
                self.last_attempts[image_index] = 0
        record = {
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "image_count": len(images),
            "image_indices": image_indices,
            "workers": min(self.max_workers, len(missing_images)),
            "max_retries": self.max_retries,
            "cache_dir": str(self.cache_dir) if self.cache_dir is not None else None,
            "cache_hits": sorted(cached),
            "cache_misses": [image_index for image_index, _ in missing_images],
        }

        def evaluate_with_retry(image_index: int, image: Atoms) -> ImageEvaluation:
            last_error: Exception | None = None
            for attempt in range(self.max_retries + 1):
                try:
                    result = self._evaluate_one(image_index, image)
                    self.last_attempts[image_index] = attempt + 1
                    return result
                except Exception as exc:
                    last_error = exc
                    reset = getattr(getattr(image, "calc", None), "reset", None)
                    if attempt < self.max_retries and callable(reset):
                        reset()
            assert last_error is not None
            self.last_attempts[image_index] = self.max_retries + 1
            raise RuntimeError(
                f"image {image_index} failed after {self.max_retries + 1} evaluation attempt(s): {last_error}"
            ) from last_error

        try:
            if missing_images:
                failure: Exception | None = None
                with ThreadPoolExecutor(
                    max_workers=min(self.max_workers, len(missing_images)),
                    thread_name_prefix="vcneb-image",
                ) as pool:
                    futures = [
                        pool.submit(evaluate_with_retry, image_index, image)
                        for image_index, image in missing_images
                    ]
                    # Collect every future even after one image fails.  A
                    # successful worker may already have completed useful DFT
                    # work; persist that result before propagating the batch
                    # failure so a restart can reuse it from the cache.
                    for (image_index, image), future in zip(missing_images, futures):
                        try:
                            result = future.result()
                        except Exception as exc:
                            if failure is None:
                                failure = exc
                        else:
                            cached[image_index] = result
                            self._store_cached(image_index, image, result)
                if failure is not None:
                    raise failure
        except Exception as exc:
            record["cache_hits"] = sorted(
                index for index in image_indices if self.last_attempts.get(index) == 0
            )
            record["cache_misses"] = [
                index for index in image_indices if self.last_attempts.get(index, 0) > 0
            ]
            record.update(
                status="failed",
                attempts=dict(sorted(self.last_attempts.items())),
                error=repr(exc),
                elapsed_s=time.monotonic() - started,
            )
            self._write_manifest_record(record)
            raise
        record["cache_hits"] = sorted(index for index in image_indices if self.last_attempts.get(index) == 0)
        record["cache_misses"] = [index for index in image_indices if self.last_attempts.get(index, 0) > 0]
        record.update(
            status="ok",
            attempts=dict(sorted(self.last_attempts.items())),
            elapsed_s=time.monotonic() - started,
        )
        self._write_manifest_record(record)
        results = [cached[index] for index in image_indices]
        return results


__all__ = ["ImageEvaluation", "ThreadedCalculatorExecutor"]
