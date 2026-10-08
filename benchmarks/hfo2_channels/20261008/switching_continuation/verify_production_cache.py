"""Verify an exact restart with the immutable production factory; forbid DFT.

Run from a new verification namespace on HF. This checks the actual old
production calculator, not a replacement mock or the current working tree.
It does not submit jobs, mutate the source SCFs, or restore FIRE momentum.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production-source", type=Path, required=True)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = args.production_source.resolve()
    seed = args.seed.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError("refusing an existing cache-verification namespace")
    sys.path.insert(0, str(source))
    import ase
    import numpy as np
    from ase.io import read
    import examples.hfo2_fixed_input_factory as adapter
    import vcneb

    if (Path(adapter.__file__).resolve() != source / "examples/hfo2_fixed_input_factory.py"
            or Path(vcneb.__file__).resolve() != source / "vcneb/__init__.py"):
        raise ValueError("imports must resolve inside the declared immutable production source")
    manifest = json.loads((seed / "manifest.json").read_text())
    if (manifest["purpose"] != "ordinary_observation_restart"
            or manifest["physical_inputs_changed"] or manifest["old_optimizer_state_reused"]
            or not manifest["all_seed_SCFs_reused"] or manifest["new_DFT_calls"] != 0):
        raise ValueError("exact-cache ordinary geometry continuation required")
    if any(sha256(seed / name) != digest for name, digest in manifest["seed_file_sha256"].items()):
        raise ValueError("restart seed changed")

    def forbidden(*_args, **_kwargs):
        raise RuntimeError("DFT is forbidden during cache verification")

    adapter.FixedHfo2Calculator.calculate = forbidden
    command = "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus"
    factory = adapter.make_seed_cached_factory(
        parameters=json.loads((seed / "factory_parameters.json").read_text()), command=command,
    )
    images = read(seed / "seed.traj", index=":")
    if len(images) != manifest["n_total_images"] or any(image.calc is not None for image in images):
        raise ValueError("expected the complete calculator-free restart seed")
    output.mkdir(parents=True, exist_ok=False)
    for index, (image, raw) in enumerate(zip(images, manifest["raw_image_evaluations"])):
        image.calc = factory(index, image, output / f"image_{index:04d}")
        if (not np.isclose(image.get_potential_energy(), raw["energy_eV_cell"], atol=1e-10, rtol=0)
                or not np.allclose(image.get_forces(), raw["forces_eV_A"], atol=1e-12, rtol=0)
                or not np.allclose(image.get_stress(), raw["stress_ASE_voigt_eV_A3"], atol=1e-12, rtol=0)):
            raise ValueError("actual production cache differs from audited numeric evidence")
    chain = vcneb.VCNEB(images, k=.2, climb=False, pressure=0.)
    fmax = float(np.linalg.norm(chain.get_forces(), axis=1).max())
    if abs(fmax - manifest["restart_fmax_eV_A"]) > 1e-10:
        raise ValueError("actual cached production force replay differs")
    record = {
        "status": "passed_actual_production_factory_without_DFT", "ASE": ase.__version__,
        "n_total_images": len(images), "cached_evaluations": len(images), "new_SCF_calls": 0,
        "fmax_eV_A": fmax, "physical_inputs_changed": False,
        "production_source": str(source),
        "production_factory_sha256": sha256(Path(adapter.__file__)),
        "verification_script_sha256": sha256(Path(__file__)),
        "seed_manifest_sha256": sha256(seed / "manifest.json"),
        "seed_traj_sha256": sha256(seed / "seed.traj"),
    }
    (output / "verification.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record))


if __name__ == "__main__":
    main()
