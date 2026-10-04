#!/bin/bash
# Reproduction outline for the immutable hf experiment; paths are explicit.
# The preparer refuses to overwrite an existing trial case.
set -euo pipefail
root=/public/home/iai806/abacus/agent-runs/20261004-varneb-gan600-ts-certificate
python=/public/home/iai806/.conda/envs/icu/bin/python
cd "${root}/source"
"${python}" -m scripts.prepare_gan_600eV_ts_stationarity_trial \
  --source /public/home/iai806/abacus/agent-runs/gan_600eV_ts_newton_canary_20260928 \
  --hessian-npz "${root}/inputs/joint_hessian.npz" \
  --normal-report "${root}/inputs/normal_strain_audit.json" \
  --output "${root}/case"
# After inspecting case/manifest.json, submit once:
# RUN_DFT=1 CASE_ROOT="${root}/case" sbatch cluster/hf_gan_600eV_ts_stationarity_trial.slurm
