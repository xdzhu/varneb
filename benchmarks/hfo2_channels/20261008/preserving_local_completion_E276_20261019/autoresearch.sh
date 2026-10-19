#!/bin/bash
# Offline diagnosis only; never launches ABACUS or scheduler jobs.
set -euo pipefail
[[ $# == 1 ]] || { echo "usage: autoresearch.sh FRESH_OUTPUT_DIRECTORY" >&2; exit 2; }
python "$(dirname "$0")/diagnose.py" --output "$1"
