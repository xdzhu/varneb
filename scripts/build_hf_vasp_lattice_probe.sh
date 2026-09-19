#!/bin/bash
# Diagnostic only: links a user's licensed object, never modifies VASP.
# These flags match this hf VASP installation, not arbitrary VASP builds.
set -euo pipefail
[[ $# == 2 ]] || { echo "usage: $0 VASP_ROOT NEW_OUTPUT_DIRECTORY" >&2; exit 2; }
vasp_root=$1
output=$2
wrapper=$(dirname "$(readlink -f "$0")")/vasp_lattice_probe.f90
mkdir "$output"  # exclusive; old canaries and binaries remain untouched
ifort -O2 -march=core-avx2 "$wrapper" "$vasp_root/build/std/lattlib.o" -o "$output/vasp_lattice_probe"
{
  echo 'compiler: ifort; flags: -O2 -march=core-avx2'
  sha256sum "$wrapper" "$vasp_root/build/std/lattlib.o" "$vasp_root/bin/vasp_std" "$output/vasp_lattice_probe"
  echo 'Not validated until compared against real VASP on regression structures.'
} > "$output/build_manifest.txt"
