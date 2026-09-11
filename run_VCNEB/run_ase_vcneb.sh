#!/bin/bash
#PBS -N run_VCNEB
#PBS -q gold6248
#PBS -l nodes=1:ppn=40
#PBS -o $PBS_JOBID.log
#PBS -e $PBS_JOBID.err
#PBS -l walltime=500:00:00:00
ulimit -s unlimited
set -euo pipefail

VASP_BIN="${VASP_BIN:-/home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std}"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python)}"

RUN_BASE="${RUN_BASE:-${PBS_O_WORKDIR:-$(pwd)}}"
ROOT_DIR="${ROOT_DIR:-$(cd "$RUN_BASE/.." && pwd -P)}"
INITIAL_STATE_DIR="${INITIAL_STATE_DIR:-$ROOT_DIR/initial_state/relax}"
FINAL_STATE_DIR="${FINAL_STATE_DIR:-$ROOT_DIR/final_state/relax}"
VCNEB_WORKDIR="${VCNEB_WORKDIR:-$RUN_BASE/run}"

N_IMAGES="${N_IMAGES:-7}"
FMAX="${FMAX:-0.05}"
MAX_STEPS="${MAX_STEPS:-300}"
SPRING_K="${SPRING_K:-0.10}"
PRESSURE_GPA="${PRESSURE_GPA:-0.0}"
OPTIMIZER="${OPTIMIZER:-FIRE}"

if [[ -n "${PBS_NODEFILE:-}" && -f "$PBS_NODEFILE" ]]; then
    NP="$(wc -l < "$PBS_NODEFILE")"
else
    NP="${NP:-40}"
fi

run_job() {
    cd "$ROOT_DIR"
    echo "=== ASE VC-NEB start ==="
    echo "ROOT_DIR        : $ROOT_DIR"
    echo "INITIAL_STATE   : $INITIAL_STATE_DIR"
    echo "FINAL_STATE     : $FINAL_STATE_DIR"
    echo "VCNEB_WORKDIR   : $VCNEB_WORKDIR"
    echo "N_IMAGES        : $N_IMAGES"
    echo "NP              : $NP"

    "$PYTHON_BIN" examples/run_vcneb_vasp.py \
        --initial "$INITIAL_STATE_DIR" \
        --final "$FINAL_STATE_DIR" \
        --workdir "$VCNEB_WORKDIR" \
        --n-images "$N_IMAGES" \
        --fmax "$FMAX" \
        --steps "$MAX_STEPS" \
        --k "$SPRING_K" \
        --pressure-gpa "$PRESSURE_GPA" \
        --optimizer "$OPTIMIZER" \
        --vasp-bin "$VASP_BIN" \
        --ncores "$NP"
}

submit_job() {
    if ! command -v qsub >/dev/null 2>&1; then
        echo "[ERROR] qsub not found"
        exit 1
    fi
    qsub -V "$SCRIPT_PATH"
}

main() {
    case "${1:-}" in
        -h|--help)
            echo "Usage: $0 [--local|--submit]"
            exit 0
            ;;
        --submit)
            submit_job
            ;;
        --local|"")
            run_job
            ;;
        *)
            echo "Unknown argument: $1"
            exit 2
            ;;
    esac
}

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
main "$@"
