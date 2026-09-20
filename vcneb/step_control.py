"""Bounded candidate-only backtracking; never retries failed DFT calculations."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from ase.io import write
from ase.optimize import FIRE


class CandidateStepRejected(ValueError):
    """An uncommitted geometry is unsuitable; no electronic job was launched."""

    def __init__(self, message, *, details=None):
        super().__init__(message)
        self.details = [] if details is None else details
        self.candidate_coordinates = None


class CheckedFIRE(FIRE):
    """FIRE with fixed-direction geometry backtracking and momentum reset.

    The target must reject candidates atomically and attach their generalized
    coordinates to CandidateStepRejected. Only that exception is recoverable.
    On acceptance of a shortened step, velocities are reset and dt reduced;
    forces, spring constants, DFT settings and convergence targets are unchanged.
    Both old ASE ``v`` and new ASE ``vel`` velocity names are supported.
    """

    _state_fields = ("v", "vel", "dt", "a", "Nsteps", "e_last", "r_last", "v_last", "vel_last")

    def __init__(self, atoms, *, max_candidate_retries=8, candidate_retry_factor=0.5,
                 candidate_manifest=None, **kwargs):
        if (isinstance(max_candidate_retries, bool) or int(max_candidate_retries) != max_candidate_retries
                or max_candidate_retries < 0):
            raise ValueError("max_candidate_retries must be a nonnegative integer")
        if max_candidate_retries and not callable(getattr(atoms, "candidate_validator", None)):
            raise ValueError("CheckedFIRE retries require an atomic candidate-validator target")
        if not np.isfinite(candidate_retry_factor) or not 0 < candidate_retry_factor < 1:
            raise ValueError("candidate_retry_factor must be in (0, 1)")
        if kwargs.get("downhill_check", False):
            raise ValueError("geometry backtracking is not combined with FIRE downhill_check")
        self.max_candidate_retries = int(max_candidate_retries)
        self.candidate_retry_factor = float(candidate_retry_factor)
        self.candidate_manifest = None if candidate_manifest is None else Path(candidate_manifest)
        self.candidate_step_history = []
        super().__init__(atoms, **kwargs)

    def _record(self, event, fraction, coordinates, error=None):
        record = {"event": event, "optimizer_step": int(self.nsteps),
                  "proposal_fraction": fraction,
                  "coordinates_sha256": sha256(np.asarray(coordinates, dtype="<f8").tobytes()).hexdigest(),
                  "error": None if error is None else str(error),
                  "details": [] if error is None else error.details,
                  "electronic_job_launched_for_rejected_candidate": False}
        if error is not None and self.candidate_manifest is not None:
            try:
                directory = self.candidate_manifest.parent / "candidate_step_artifacts"
                directory.mkdir(parents=True, exist_ok=True)
                artifact = directory / (
                    f"step_{int(self.nsteps):04d}_event_{len(self.candidate_step_history):04d}.traj"
                )
                preview = self.atoms.candidate_images_from_x(coordinates)
                write(artifact, preview, format="traj")
                record["candidate_trajectory"] = str(artifact)
            except (AttributeError, OSError, ValueError) as artifact_error:
                record["candidate_trajectory_error"] = str(artifact_error)
        self.candidate_step_history.append(record)
        if self.candidate_manifest is not None:
            self.candidate_manifest.parent.mkdir(parents=True, exist_ok=True)
            with self.candidate_manifest.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(record, allow_nan=False) + "\n")

    def step(self, f=None):
        original = self.atoms.get_x().copy()
        saved = {name: deepcopy(self.__dict__[name]) for name in self._state_fields if name in self.__dict__}
        try:
            return super().step(f)
        except CandidateStepRejected as error:
            if not np.array_equal(self.atoms.get_x(), original):
                raise RuntimeError("candidate target violated atomic rejection") from error
            proposed = error.candidate_coordinates
            if proposed is None or np.asarray(proposed).shape != original.shape:
                raise RuntimeError("rejected candidate coordinates were not provided") from error
            delta = np.asarray(proposed) - original
            self._record("rejected", 1.0, proposed, error)
            last_error = error
        for retry in range(1, self.max_candidate_retries + 1):
            fraction = self.candidate_retry_factor ** retry
            candidate = original + fraction * delta
            try:
                self.atoms.set_x(candidate)
            except CandidateStepRejected as error:
                if not np.array_equal(self.atoms.get_x(), original):
                    raise RuntimeError("candidate target violated atomic rejection") from error
                self._record("rejected", fraction, candidate, error)
                last_error = error
                continue
            # Do not retain momentum inconsistent with the shortened move.
            velocity = self.vel if hasattr(self, "vel") else self.v
            velocity[...] = 0
            self.dt = min(float(saved["dt"]), self.dt) * fraction
            self.a = self.astart
            self.Nsteps = 0
            self.dump((velocity, self.dt))
            self._record("accepted_backtracked", fraction, candidate)
            return
        for name in self._state_fields:
            if name in saved:
                self.__dict__[name] = saved[name]
            else:
                self.__dict__.pop(name, None)
        self._record("exhausted", 0.0, original, last_error)
        raise last_error
