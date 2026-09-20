"""Image-count-normalized and staged FIRE controls for VCNEB bands."""

from __future__ import annotations

import numpy as np

from .step_control import CheckedFIRE


class ImageScaledFIRE(CheckedFIRE):
    """Global FIRE whose input ``maxstep`` is interpreted per interior image.

    ASE caps the norm of the full concatenated band.  Multiplying the desired
    per-image displacement scale by ``sqrt(n_interior)`` removes the accidental
    image-count dependence while retaining one global FIRE state.
    """

    def __init__(self, atoms, *, maxstep=0.02, **kwargs):
        n_interior = max(1, int(atoms.n_images) - 2)
        self.per_image_maxstep = float(maxstep)
        if self.per_image_maxstep <= 0.0:
            raise ValueError("maxstep must be positive")
        self.band_scale = float(np.sqrt(n_interior))
        kwargs.setdefault("max_candidate_retries", 0)
        super().__init__(atoms, maxstep=self.per_image_maxstep * self.band_scale, **kwargs)


class StagedFIRE(ImageScaledFIRE):
    """Use image-scaled FIRE far from convergence and a smaller final cap.

    The switch uses the already available VCNEB force and therefore launches
    no additional calculator evaluation.  Momentum is reset at the switch so
    the fine stage does not inherit the coarse-stage overshoot.
    """

    def __init__(
        self,
        atoms,
        *,
        maxstep=0.02,
        switch_fmax=0.12,
        refine_maxstep=0.02,
        **kwargs,
    ):
        if float(switch_fmax) <= 0.0 or float(refine_maxstep) <= 0.0:
            raise ValueError("switch_fmax and refine_maxstep must be positive")
        self.switch_fmax = float(switch_fmax)
        self.refine_maxstep = float(refine_maxstep)
        self.coarse_dt = float(kwargs.get("dt", 0.1))
        self.switched_to_refine = False
        self.switch_history = []
        super().__init__(atoms, maxstep=maxstep, **kwargs)

    def _current_fmax(self, forces):
        gradient_norm = getattr(self.atoms, "gradient_norm", None)
        if callable(gradient_norm):
            return float(gradient_norm(-np.asarray(forces, dtype=float)))
        vectors = np.asarray(forces, dtype=float).reshape(-1, 3)
        return float(np.max(np.linalg.norm(vectors, axis=1)))

    def _switch(self, current_fmax):
        self.maxstep = self.refine_maxstep
        self.dt = self.coarse_dt
        self.a = self.astart
        self.Nsteps = 0
        velocity = getattr(self, "vel", getattr(self, "v", None))
        if velocity is not None:
            velocity[...] = 0.0
        self.switched_to_refine = True
        self.switch_history.append(
            {
                "optimizer_step": int(self.nsteps),
                "fmax_eV_per_A": float(current_fmax),
                "coarse_band_maxstep": self.per_image_maxstep * self.band_scale,
                "refine_band_maxstep": self.refine_maxstep,
                "momentum_reset": True,
                "additional_calculator_evaluations": 0,
            }
        )

    def step(self, f=None):
        forces = self.atoms.get_forces() if f is None else f
        current_fmax = self._current_fmax(forces)
        if not self.switched_to_refine and current_fmax <= self.switch_fmax:
            self._switch(current_fmax)
        # Pass the force array through so checking the stage does not trigger a
        # second electronic-structure evaluation.
        return super().step(forces)
