#!/bin/bash
#PBS -N run_NEB
#PBS -q gold6248
#PBS -l nodes=1:ppn=40
#PBS -o $PBS_JOBID.log
#PBS -e $PBS_JOBID.err
#PBS -l walltime=500:00:00:00
ulimit -s unlimited
set -euo pipefail

# ------------------- 配置 -------------------
# VASP_BIN="${VASP_BIN:-/home/zhuxd/Software/src/vasp/5.4.4_fixC/bin/vasp_std}"
VASP_BIN="${VASP_BIN:-/home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std}"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python)}"

RUN_BASE="${RUN_BASE:-${PBS_O_WORKDIR:-$(pwd)}}"
ROOT_DIR="${ROOT_DIR:-$(cd "$RUN_BASE/.." && pwd -P)}"
INITIAL_STATE_DIR="${INITIAL_STATE_DIR:-$ROOT_DIR/initial_state}"
FINAL_STATE_DIR="${FINAL_STATE_DIR:-$ROOT_DIR/final_state}"
NEB_WORKDIR="${NEB_WORKDIR:-$RUN_BASE/run}"

N_INTERMEDIATE="${N_INTERMEDIATE:-5}"
FMAX="${FMAX:-0.05}"
MAX_STEPS="${MAX_STEPS:-300}"
SPRING_K="${SPRING_K:-0.10}"

# ------------------- 运行函数 -------------------
run_job() {
    echo "=== ASE CI-NEB start ==="
    echo "RUN_BASE        : $RUN_BASE"
    echo "INITIAL_STATE   : $INITIAL_STATE_DIR"
    echo "FINAL_STATE     : $FINAL_STATE_DIR"
    echo "NEB_WORKDIR     : $NEB_WORKDIR"
    echo "N_INTERMEDIATE  : $N_INTERMEDIATE"

    mkdir -p "$NEB_WORKDIR"

    "$PYTHON_BIN" - <<PY
import os
import shutil
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.calculators.vasp import Vasp
from ase.calculators.singlepoint import SinglePointCalculator
from ase.mep import NEB
from ase.optimize import FIRE
from ase.mep import NEBTools

# ------------------- 自动检测核数 -------------------
if "PBS_NODEFILE" in os.environ:
    NP = len(open(os.environ["PBS_NODEFILE"]).readlines())
else:
    try:
        import multiprocessing
        NP = multiprocessing.cpu_count()
    except:
        NP = 1

VASP_BIN = os.environ.get("VASP_BIN", "$VASP_BIN")
WORKDIR = Path(os.environ.get("NEB_WORKDIR", "$NEB_WORKDIR"))
N_INTERMEDIATE = int(os.environ.get("N_INTERMEDIATE", "$N_INTERMEDIATE"))
FMAX = float(os.environ.get("FMAX", "$FMAX"))
MAX_STEPS = int(os.environ.get("MAX_STEPS", "$MAX_STEPS"))
SPRING_K = float(os.environ.get("SPRING_K", "$SPRING_K"))

INITIAL_DIR = Path(os.environ.get("INITIAL_STATE_DIR", "$INITIAL_STATE_DIR")) / "relax"
FINAL_DIR = Path(os.environ.get("FINAL_STATE_DIR", "$FINAL_STATE_DIR")) / "relax"

# ------------------- 读取端点 -------------------
def read_endpoint(d):
    candidates = [d/f for f in ["CONTCAR", "CONTCAR.vasp", "POSCAR"] if (d/f).exists()]
    if candidates:
        atoms = read(candidates[0])
        print(f"[OK] Read {candidates[0]}")
        return atoms
    else:
        raise FileNotFoundError(f"No endpoint structure found in {d}")

def collect_vasp_params(calc):
    params = {}
    for name in [
        "float_params",
        "exp_params",
        "string_params",
        "int_params",
        "bool_params",
        "list_int_params",
        "list_bool_params",
        "list_float_params",
        "special_params",
        "dict_params",
        "input_params",
    ]:
        values = getattr(calc, name, None)
        if values is None:
            continue
        try:
            if len(values) == 0:
                continue
        except TypeError:
            pass
        for key, value in values.items():
            if value is None:
                continue
            if name == "input_params":
                if isinstance(value, dict) and len(value) == 0:
                    continue
                if isinstance(value, (list, tuple)) and len(value) == 0:
                    continue
            params[key] = value
    return params

def read_vasp_inputs(source_dir):
    incar_path = source_dir / "INCAR"
    kpoints_path = source_dir / "KPOINTS"
    potcar_path = source_dir / "POTCAR"
    if not incar_path.exists():
        raise FileNotFoundError(f"INCAR not found: {incar_path}")
    if not kpoints_path.exists():
        raise FileNotFoundError(f"KPOINTS not found: {kpoints_path}")
    if not potcar_path.exists():
        raise FileNotFoundError(f"POTCAR not found: {potcar_path}")

    incar_reader = Vasp()
    if not hasattr(incar_reader, "read_incar"):
        raise AttributeError("ASE Vasp calculator does not provide read_incar().")
    incar_reader.read_incar(str(incar_path))
    incar_params = collect_vasp_params(incar_reader)
    for runtime_key in ("directory", "command", "txt", "label", "atoms"):
        incar_params.pop(runtime_key, None)

    kpoints_reader = Vasp()
    if not hasattr(kpoints_reader, "read_kpoints"):
        raise AttributeError("ASE Vasp calculator does not provide read_kpoints().")
    kpoints_reader.read_kpoints(str(kpoints_path))
    kpoint_params = {
        key: value
        for key, value in collect_vasp_params(kpoints_reader).items()
        if key in {"kpts", "gamma", "reciprocal", "kpts_nintersections"}
    }
    if not kpoint_params:
        raise ValueError(f"Failed to parse KPOINTS with ASE reader: {kpoints_path}")

    return incar_params, kpoint_params, potcar_path

def make_static_vasp(directory, static_params):
    return Vasp(
        directory=str(directory),
        command=f"mpirun -np {NP} {VASP_BIN}",
        txt="vasp.out",
        **static_params,
    )

def snapshot_image(image, require_forces=False):
    energy = image.get_potential_energy()
    results = {"energy": energy}
    if require_forces:
        results["forces"] = image.get_forces()
    snap = image.copy()
    snap.calc = SinglePointCalculator(snap, **results)
    snap.info["energy"] = energy
    return snap

def write_images_with_energy(path, images, mode="w", require_forces=False):
    from ase.io.trajectory import Trajectory

    with Trajectory(str(path), mode) as traj:
        for image in images:
            traj.write(snapshot_image(image, require_forces=require_forces))

initial = read_endpoint(INITIAL_DIR)
final = read_endpoint(FINAL_DIR)

if len(initial) != len(final):
    raise ValueError("Initial and final states have different number of atoms.")
final.set_cell(initial.cell, scale_atoms=False)
final.set_pbc(initial.get_pbc())

# ------------------- 创建 NEB 镜像 -------------------
images = [initial.copy()]
for _ in range(N_INTERMEDIATE):
    images.append(initial.copy())
images.append(final.copy())

neb = NEB(images, climb=True, k=SPRING_K)
neb.interpolate(mic=True)

# ------------------- IDPP 优化 -------------------
try:
    from ase.mep.neb import idpp_interpolate
except ImportError:
    try:
        from ase.mep import idpp_interpolate
    except ImportError:
        idpp_interpolate = None

if idpp_interpolate is not None:
    idpp_interpolate(images, mic=True, steps=100)
    print("[INFO] IDPP interpolation done")
else:
    print("[WARN] ASE version does not support IDPP; using linear interpolation")

# ------------------- 读取 INCAR / KPOINTS 并设置静态自洽参数 -------------------
incar_params, kpoint_params, potcar_path = read_vasp_inputs(INITIAL_DIR)
static_overrides = {
    "ibrion": -1,
    "nsw": 0,
    "isif": 2,
    "isym": 0,
    "lcharg": False,
    "lwave": False,
    "xc": "PBE",
    "pp": "PBE",
}
static_params = dict(incar_params)
static_params.update(kpoint_params)
static_params.update(static_overrides)
print(f"[INFO] Loaded INCAR via ASE reader: {INITIAL_DIR / 'INCAR'}")
print(f"[INFO] Loaded KPOINTS via ASE reader: {INITIAL_DIR / 'KPOINTS'}")
print(f"[INFO] Reference POTCAR found: {potcar_path}")
print(f"[INFO] Static overrides applied: {sorted(static_overrides)}")

for i, img in enumerate(images):
    img_dir = WORKDIR / f"{i:02d}"
    img_dir.mkdir(parents=True, exist_ok=True)
    write(img_dir/"POSCAR.start", img, format="vasp", direct=True, vasp5=True)
    shutil.copy2(potcar_path, img_dir / "POTCAR")
    img.calc = make_static_vasp(img_dir, static_params)

# ------------------- 输出初始路径 traj（确保每张图像有能量） -------------------
initial_traj = WORKDIR / "initial_neb.traj"
write_images_with_energy(initial_traj, images, mode="w")
print(f"[INFO] Initial NEB path written to {initial_traj}")

# ------------------- FIRE 优化 NEB -------------------
neb_traj = WORKDIR / "neb.traj"
if neb_traj.exists():
    neb_traj.unlink()

def append_neb_snapshot():
    write_images_with_energy(neb_traj, images, mode="a", require_forces=True)

append_neb_snapshot()
opt = FIRE(neb, logfile=str(WORKDIR/"neb.opt.log"))
opt.attach(append_neb_snapshot, interval=1)
opt.run(fmax=FMAX, steps=MAX_STEPS)
print(f"[INFO] NEB optimization finished, trajectory: {neb_traj}")

# ------------------- 保存最终结构 -------------------
for i, img in enumerate(images):
    write(WORKDIR/f"{i:02d}"/"POSCAR.final", img, format="vasp", direct=True, vasp5=True)

# ------------------- NEBTools 分析绘图 -------------------
tools = NEBTools(images)
barrier, dE = tools.get_barrier(fit=True)
with open(WORKDIR/"neb_summary.txt","w") as f:
    f.write(f"Forward barrier (fit) = {barrier:.8f} eV\n")
    f.write(f"Reaction energy dE    = {dE:.8f} eV\n")

tools.plot_band()
plt.tight_layout()
plt.savefig(WORKDIR/"neb_barrier.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"[DONE] NEB finished. Barrier={barrier:.4f} eV, dE={dE:.4f} eV")
PY
}

# ------------------- PBS 提交 -------------------
submit_job() {
    if ! command -v qsub >/dev/null 2>&1; then
        echo "[ERROR] qsub not found"
        exit 1
    fi
    qsub -V "$SCRIPT_PATH"
}

# ------------------- 主入口 -------------------
main() {
    case "${1:-}" in
        -h|--help) echo "Usage: $0 [--local|--submit]"; exit 0 ;;
        --submit) submit_job; exit 0 ;;
        --local) run_job; exit 0 ;;
        *) run_job ;;
    esac
}

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
main "$@"
