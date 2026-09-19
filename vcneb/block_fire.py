"""Block-adaptive FIRE optimizers for long variable-cell NEB chains.

ASE's FIRE treats a complete NEB band as one vector.  Consequently its
``maxstep`` bounds the norm of *all* image displacements together, and a bad
candidate in one image reduces the time step of the complete band.  The
optimizers in this module keep independent FIRE state and trust radii for
either each image or each atomic/cell block of an image.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase.optimize.optimize import Optimizer

from .step_control import CandidateStepRejected


class BlockFIRE(Optimizer):
    """FIRE with independent state and displacement caps for VCNEB blocks.

    ``block_mode="image"`` assigns one FIRE state to every interior image.
    ``block_mode="atomic-cell"`` further separates its atomic and cell
    generalized coordinates.  Candidate backtracking is local to the image(s)
    named by the validator, so a problematic cell cannot freeze the band.
    """

    def __init__(
        self,
        atoms,
        restart=None,
        logfile="-",
        trajectory=None,
        *,
        block_mode="image",
        dt=0.1,
        dtmax=1.0,
        dtmin=1.0e-4,
        maxstep=0.2,
        cell_maxstep=None,
        Nmin=5,
        finc=1.1,
        fdec=0.5,
        astart=0.1,
        fa=0.99,
        max_candidate_retries=0,
        candidate_retry_factor=0.5,
        candidate_manifest=None,
        **kwargs,
    ):
        if block_mode not in {"image", "atomic-cell"}:
            raise ValueError("block_mode must be 'image' or 'atomic-cell'")
        if not hasattr(atoms, "image_ndofs") or not hasattr(atoms, "n_images"):
            raise ValueError("BlockFIRE requires a VCNEB-like optimizer target")
        if not 0 < float(dtmin) <= float(dt) <= float(dtmax):
            raise ValueError("Require 0 < dtmin <= dt <= dtmax")
        if float(maxstep) <= 0.0 or (cell_maxstep is not None and float(cell_maxstep) <= 0.0):
            raise ValueError("maxstep and cell_maxstep must be positive")
        if not 0 < float(candidate_retry_factor) < 1:
            raise ValueError("candidate_retry_factor must be in (0, 1)")
        if int(max_candidate_retries) != max_candidate_retries or max_candidate_retries < 0:
            raise ValueError("max_candidate_retries must be a nonnegative integer")

        self.block_mode = block_mode
        self.dt0 = float(dt)
        self.dtmax = float(dtmax)
        self.dtmin = float(dtmin)
        self.maxstep = float(maxstep)
        self.cell_maxstep = self.maxstep if cell_maxstep is None else float(cell_maxstep)
        self.Nmin = int(Nmin)
        self.finc = float(finc)
        self.fdec = float(fdec)
        self.astart = float(astart)
        self.fa = float(fa)
        self.max_candidate_retries = int(max_candidate_retries)
        self.candidate_retry_factor = float(candidate_retry_factor)
        self.candidate_manifest = None if candidate_manifest is None else Path(candidate_manifest)
        self.candidate_step_history = []
        self._slices, self._image_for_block, self._caps = self._make_blocks(atoms)
        super().__init__(atoms, restart, logfile, trajectory, **kwargs)

    def _make_blocks(self, chain):
        image_size = int(chain.image_ndofs)
        n_atoms = len(chain.images[0])
        atomic_size = 3 * n_atoms
        slices = []
        image_for_block = []
        caps = []
        for offset, image_index in enumerate(range(1, int(chain.n_images) - 1)):
            start = offset * image_size
            if self.block_mode == "image":
                slices.append(slice(start, start + image_size))
                image_for_block.append(image_index)
                caps.append(self.maxstep)
            else:
                slices.extend(
                    [slice(start, start + atomic_size), slice(start + atomic_size, start + image_size)]
                )
                image_for_block.extend([image_index, image_index])
                caps.extend([self.maxstep, self.cell_maxstep])
        return slices, np.asarray(image_for_block, dtype=int), np.asarray(caps, dtype=float)

    def initialize(self):
        self.vel = None
        self.block_dt = None
        self.block_a = None
        self.block_nsteps = None

    def read(self):
        self.vel, self.block_dt, self.block_a, self.block_nsteps = self.load()

    def _ensure_state(self, ndofs):
        nblocks = len(self._slices)
        if self.vel is None:
            self.vel = np.zeros(ndofs)
            self.block_dt = np.full(nblocks, self.dt0)
            self.block_a = np.full(nblocks, self.astart)
            self.block_nsteps = np.zeros(nblocks, dtype=int)
            return True
        return False

    def _record_candidate(self, event, scales, error=None):
        record = {
            "event": event,
            "optimizer_step": int(self.nsteps),
            "block_mode": self.block_mode,
            "image_scales": {
                str(image): float(min(scales[self._image_for_block == image]))
                for image in sorted(set(self._image_for_block.tolist()))
            },
            "details": [] if error is None else error.details,
            "electronic_job_launched_for_rejected_candidate": False,
        }
        self.candidate_step_history.append(record)
        if self.candidate_manifest is not None:
            self.candidate_manifest.parent.mkdir(parents=True, exist_ok=True)
            with self.candidate_manifest.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(record, allow_nan=False) + "\n")

    @staticmethod
    def _rejected_images(error):
        images = set()
        for item in error.details:
            if isinstance(item, dict) and "image_index" in item:
                images.add(int(item["image_index"]))
        return images

    def _commit_with_local_backtracking(self, original, displacement):
        scales = np.ones(len(self._slices))
        last_error = None
        for attempt in range(self.max_candidate_retries + 1):
            candidate = original.copy()
            for block, section in enumerate(self._slices):
                candidate[section] += scales[block] * displacement[section]
            try:
                self.atoms.set_x(candidate)
            except CandidateStepRejected as error:
                if not np.array_equal(self.atoms.get_x(), original):
                    raise RuntimeError("candidate target violated atomic rejection") from error
                self._record_candidate("rejected", scales, error)
                last_error = error
                if attempt == self.max_candidate_retries:
                    raise
                rejected = self._rejected_images(error)
                affected = np.ones(len(scales), dtype=bool) if not rejected else np.isin(
                    self._image_for_block, list(rejected)
                )
                scales[affected] *= self.candidate_retry_factor
                continue
            shortened = scales < 1.0
            if np.any(shortened):
                # The accepted geometry is already shortened.  Reset only the
                # inconsistent block momentum and retain a recoverable dt floor.
                for block in np.flatnonzero(shortened):
                    section = self._slices[block]
                    self.vel[section] = 0.0
                    self.block_dt[block] = max(
                        self.dtmin,
                        self.block_dt[block] * max(scales[block], self.fdec),
                    )
                    self.block_a[block] = self.astart
                    self.block_nsteps[block] = 0
                self._record_candidate("accepted_backtracked", scales)
            return
        raise last_error  # pragma: no cover - loop always returns or raises

    def step(self, f=None):
        # ASE >=3.27 exposes Optimizer._get_gradient/get_x, whereas the
        # supported cluster ASE 3.23 optimizer directly consumes force arrays.
        # VCNEB's public target API is stable across both environments.
        gradient = np.asarray(self.atoms.get_forces() if f is None else f, dtype=float).ravel()
        fresh = self._ensure_state(gradient.size)

        for block, section in enumerate(self._slices):
            if fresh:
                continue
            force = gradient[section]
            velocity = self.vel[section]
            vf = float(np.vdot(force, velocity))
            force2 = float(np.vdot(force, force))
            if vf > 0.0 and force2 > 0.0:
                speed = np.sqrt(np.vdot(velocity, velocity))
                alpha = self.block_a[block]
                velocity[:] = (1.0 - alpha) * velocity + alpha * force * speed / np.sqrt(force2)
                if self.block_nsteps[block] > self.Nmin:
                    self.block_dt[block] = min(self.block_dt[block] * self.finc, self.dtmax)
                    self.block_a[block] *= self.fa
                self.block_nsteps[block] += 1
            else:
                velocity[:] = 0.0
                self.block_a[block] = self.astart
                self.block_dt[block] = max(self.dtmin, self.block_dt[block] * self.fdec)
                self.block_nsteps[block] = 0

        displacement = np.zeros_like(gradient)
        for block, section in enumerate(self._slices):
            self.vel[section] += self.block_dt[block] * gradient[section]
            move = self.block_dt[block] * self.vel[section]
            norm = float(np.linalg.norm(move))
            cap = self._caps[block]
            if norm > cap:
                move *= cap / norm
            displacement[section] = move

        original = self.atoms.get_x()
        self._commit_with_local_backtracking(original, displacement)
        self.dump((self.vel, self.block_dt, self.block_a, self.block_nsteps))
