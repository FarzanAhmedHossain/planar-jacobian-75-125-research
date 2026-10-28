#!/usr/bin/env python3
"""Run the exact selected 10-generator F2 subsystem at Macaulay weight 30 on GPU.

The row-span calculation is over GF(251), so a positive membership result is
only reconnaissance until an exact QQ certificate is reconstructed and checked.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import f2_75_125_211_d19_support12_21_32_38_macaulay_gpu as gpu


TARGET_WEIGHT = 30
MANIFEST = HERE / "f2_75_125_211_d19_r1_support12_support38_manifest_20261028.json"
EXPECTED_INDICES = (6, 7, 8, 9, 10, 11, 14, 15, 16, 17)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def monomial_count(weights: tuple[int, ...], target: int) -> int:
    counts = [0] * (target + 1)
    counts[0] = 1
    for weight in weights:
        for degree in range(weight, target + 1):
            counts[degree] += counts[degree - weight]
    return counts[target]


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    indices = tuple(int(value) for value in manifest["source_generator_indices_zero_based"])
    if indices != EXPECTED_INDICES or len(manifest["source_labels"]) != len(indices):
        raise ValueError("the checked 10-generator source selection changed")

    all_inputs, ring = gpu.capture_compressed_system()
    if len(all_inputs) != 18:
        raise ValueError("the checked 18-generator source system changed")
    selected = [all_inputs[index] for index in indices]

    analysis_path = HERE / "f2_75_125_211_d19_support38_analysis_20261025.json"
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    weights = tuple(int(value) for value in analysis["primitive_positive_integer_weights"])
    degrees = []
    for poly in selected:
        values = {sum(w * e for w, e in zip(weights, monomial))
                  for monomial, _ in poly.terms()}
        if len(values) != 1:
            raise ValueError("selected generator is not weight homogeneous")
        degrees.append(values.pop())

    ncols = monomial_count(weights, TARGET_WEIGHT)
    nrows = sum(monomial_count(weights, TARGET_WEIGHT - degree)
                for degree in degrees if degree <= TARGET_WEIGHT)
    matrix_bytes = nrows * ncols

    platforms = gpu.cl.get_platforms()
    candidates = [(platform, device) for platform in platforms
                  for device in platform.get_devices()
                  if device.type & gpu.cl.device_type.GPU]
    if not candidates:
        raise RuntimeError("no OpenCL GPU device is available")
    platform, device = next(((p, d) for p, d in candidates if "Intel" in d.name), candidates[0])
    preflight = {
        "target_weight": TARGET_WEIGHT,
        "generator_indices_zero_based": list(indices),
        "generator_weight_degrees": degrees,
        "shape": [nrows, ncols],
        "matrix_bytes": matrix_bytes,
        "matrix_GiB": matrix_bytes / (1024 ** 3),
        "device": device.name,
        "device_max_allocation_bytes": int(device.max_mem_alloc_size),
        "device_global_memory_bytes": int(device.global_mem_size),
    }
    print(json.dumps({"stage": "preflight", **preflight}), flush=True)
    if matrix_bytes > int(device.max_mem_alloc_size) or matrix_bytes > int(device.global_mem_size):
        raise MemoryError("preflight says Macaulay matrix exceeds the selected GPU's allocation limits")

    def capture_subset():
        return selected, ring

    gpu.capture_compressed_system = capture_subset
    gpu.OUTPUT_STEM = HERE / "f2_75_125_211_d19_support12_support38_macaulay_gpu_subset10_w"
    result = gpu.run(TARGET_WEIGHT)
    result["scope"] = (
        "GF(251) Macaulay screen for the exact 10-generator necessary subset "
        "(support-12 plus support-38) of the d19 2+1+1, U1=0 system; "
        "not a QQ certificate or branch exclusion"
    )
    result["generator_indices_zero_based"] = list(indices)
    result["generator_source_labels"] = manifest["source_labels"]
    result["selection_manifest"] = MANIFEST.name
    result["selection_manifest_sha256"] = sha256(MANIFEST)
    result["gpu_kernel_script"] = Path(gpu.__file__).name
    result["gpu_kernel_script_sha256"] = sha256(Path(gpu.__file__))
    result["runner_script"] = Path(__file__).name
    result["runner_script_sha256"] = sha256(Path(__file__))
    result["matrix_preflight"] = preflight
    result["warning"] = (
        "Even a positive modular row-span result is only reconnaissance: verify "
        "the exact QQ membership and prove the selected equations are necessary "
        "on the full r!=0 chart before claiming exclusion."
    )
    output = HERE / f"f2_75_125_211_d19_support12_support38_macaulay_gpu_subset10_w{TARGET_WEIGHT}_20261028.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": "subset_report_written", "report": str(output),
                      "rank": result["rank"], "shape": result["macaulay_shape"],
                      "pure_power_membership": result["pure_power_membership"]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
