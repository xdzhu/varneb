#!/bin/bash
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "${repo}"

python examples/run_vcneb_convergence.py \
  --output outputs/vcneb_optimizer_screen.json \
  --images 5 9 17 --springs 0.15 --cell-scales 5.0 \
  --optimizers FIRE BlockFIRE SplitFIRE --fmax 0.01 --steps 300
python -m pytest -q

