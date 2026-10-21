#!/bin/sh
set -eu
# Calculator-free regression entry point; does not submit jobs.
python -m pytest -q tests/test_hfo2_strain_pattern_cross.py
