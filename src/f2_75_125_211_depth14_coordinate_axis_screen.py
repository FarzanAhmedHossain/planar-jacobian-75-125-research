#!/usr/bin/env python3
"""Streaming exact coordinate-axis screen for a large sparse lift checkpoint.

This checks only six coordinate axes of the saved 2+1+1 necessary ideal. It
does not claim a full-branch exclusion. The JSON checkpoint is streamed one
constraint at a time to avoid loading its hundreds of megabytes into RAM.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_sparse_ring_lift_211_roots14_p43zero_d14_terminal_state.json"
OUTPUT = HERE / "f2_75_125_211_depth14_coordinate_axis_screen.json"
TARGETS = tuple(f"p{i}" for i in range(36, 42))
PINNED_ZERO_PARAMETERS = ("p43",)
CHUNK = 1 << 20

sys.path.insert(0, str(HERE))
import f2_75_125_211_depth9_coordinate_axis_screen as poly  # noqa: E402


def source_hash_and_constraints_offset(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    needle = b'"constraints"'
    overlap = b""
    base = 0
    offset = None
    with path.open("rb") as stream:
        while chunk := stream.read(CHUNK):
            digest.update(chunk)
            data = overlap + chunk
            pos = 0
            while True:
                found = data.find(needle, pos)
                if found < 0:
                    break
                suffix = data[found + len(needle) :]
                j = 0
                while j < len(suffix) and suffix[j] in b" \t\r\n":
                    j += 1
                if j < len(suffix) and suffix[j] == ord(":"):
                    j += 1
                    while j < len(suffix) and suffix[j] in b" \t\r\n":
                        j += 1
                    if j < len(suffix) and suffix[j] == ord("["):
                        offset = base - len(overlap) + found + len(needle) + j + 1
                        break
                pos = found + 1
            if offset is not None:
                while chunk := stream.read(CHUNK):
                    digest.update(chunk)
                break
            keep = len(needle) + 32
            overlap = data[-keep:]
            base += len(chunk)
    if offset is None:
        raise ValueError("could not locate the top-level constraints array")
    return digest.hexdigest(), offset


def metadata(path: Path, array_offset: int) -> dict:
    with path.open("rb") as stream:
        prefix = stream.read(min(array_offset, 1 << 20))
    gen_match = re.search(rb'"generators"\s*:\s*(\[[^\]]*\])', prefix)
    branch_match = re.search(rb'"branch"\s*:\s*"([^"]+)"', prefix)
    fmt_match = re.search(rb'"format"\s*:\s*"([^"]+)"', prefix)
    quadratic_match = re.search(rb'"quadratic"\s*:\s*(true|false)', prefix)
    parameter_match = re.search(rb'"nparameters"\s*:\s*(\d+)', prefix)
    if not all((gen_match, branch_match, fmt_match, quadratic_match, parameter_match)):
        raise ValueError("checkpoint metadata is missing or exceeds the header scan")
    return {
        "format": fmt_match.group(1).decode("ascii"),
        "branch": branch_match.group(1).decode("ascii"),
        "quadratic": quadratic_match.group(1) == b"true",
        "generators": json.loads(gen_match.group(1)),
        "nparameters": int(parameter_match.group(1)),
    }


def constraint_items(path: Path, array_offset: int):
    decoder = json.JSONDecoder()
    with path.open("rb") as stream:
        stream.seek(array_offset)
        buffer = ""
        while True:
            buffer = buffer.lstrip()
            if buffer.startswith(","):
                buffer = buffer[1:]
                continue
            if buffer.startswith("]"):
                return
            if not buffer:
                raw = stream.read(CHUNK)
                if not raw:
                    raise ValueError("truncated constraints array")
                buffer = raw.decode("utf-8")
                continue
            try:
                item, end = decoder.raw_decode(buffer)
            except json.JSONDecodeError:
                raw = stream.read(CHUNK)
                if not raw:
                    raise ValueError("truncated constraint item")
                buffer += raw.decode("utf-8")
                continue
            yield item
            buffer = buffer[end:]


def main() -> None:
    positional = [arg for arg in sys.argv[1:] if arg != "--all-parameters"]
    all_parameters = "--all-parameters" in sys.argv[1:]
    source = Path(positional[0]).resolve() if positional else SOURCE
    output = Path(positional[1]).resolve() if len(positional) > 1 else OUTPUT
    digest, array_offset = source_hash_and_constraints_offset(source)
    meta = metadata(source, array_offset)
    if meta["format"] != "f2-sparse-ring-v1" or meta["branch"] != "2+1+1":
        raise ValueError(f"unexpected checkpoint identity: {meta['format']} {meta['branch']}")
    if meta["quadratic"]:
        raise ValueError("this exact screen expects rational coefficients")

    generator_index = {name: i for i, name in enumerate(meta["generators"])}
    targets = (tuple(f"p{i}" for i in range(meta["nparameters"])
                     if f"p{i}" not in PINNED_ZERO_PARAMETERS)
               if all_parameters else TARGETS)
    target_indices = {name: generator_index[name] for name in targets}
    common_gcd = {name: None for name in targets}
    restrictions = {name: 0 for name in targets}
    equations = terms_seen = 0

    for item in constraint_items(source, array_offset):
        if not isinstance(item, list) or len(item) != 2:
            raise ValueError("unexpected constraint encoding")
        tag, encoded_terms = item
        equations += 1
        row = {name: defaultdict(lambda: defaultdict(Fraction)) for name in targets}

        for sparse, values in encoded_terms:
            terms_seen += 1
            if len(values) != 1:
                raise ValueError("non-rational coefficient encoding")
            coefficient = Fraction(values[0])
            if not coefficient:
                continue

            t_degree = 0
            support: dict[str, int] = {}
            auxiliary = set()
            for index, exponent in sparse:
                name = meta["generators"][index]
                if name == "t":
                    t_degree += exponent
                elif name.startswith("x"):
                    auxiliary.add(name)
                elif name.startswith("p"):
                    support[name] = support.get(name, 0) + exponent
                else:
                    auxiliary.add(name)

            # Interpolation variables, if present in a residual, are free
            # coordinates outside this p-axis and cannot be set to zero here.
            if auxiliary:
                raise ValueError(f"residuals still contain auxiliary generators: {sorted(auxiliary)[:8]}")

            if not support:
                # A p-independent term restricts to every parameter axis.
                for target in targets:
                    row[target][t_degree][0] += coefficient
            elif len(support) == 1:
                target, power = next(iter(support.items()))
                if target in row:
                    row[target][t_degree][power] += coefficient

        for target in targets:
            for by_power in row[target].values():
                equation = poly.trim([by_power.get(k, Fraction(0)) for k in range(max(by_power, default=-1) + 1)])
                if not equation:
                    continue
                restrictions[target] += 1
                old = common_gcd[target]
                common_gcd[target] = equation if old is None else poly.gcd_poly(old, equation)

        if equations % 100 == 0:
            print(f"axis_scan equations={equations} sparse_terms={terms_seen}", flush=True)

    results = {}
    for target in targets:
        gcd_value = common_gcd[target]
        if gcd_value is None:
            results[target] = {
                "status": "entire_axis_survives_truncation",
                "nonzero_restrictions": 0,
            }
            continue
        gcd_value = poly.monic(gcd_value)
        nonzero_factor = poly.trim(gcd_value[:])
        while len(nonzero_factor) > 1 and nonzero_factor[0] == 0:
            nonzero_factor = nonzero_factor[1:]
        status = "possible_nonzero_axis_points" if len(nonzero_factor) > 1 else "axis_excluded"
        results[target] = {
            "status": status,
            "nonzero_restrictions": restrictions[target],
            "common_gcd_ascending_coefficients_Q": poly.encode_poly(gcd_value),
            "nonzero_root_factor_ascending_coefficients_Q": poly.encode_poly(nonzero_factor),
        }

    result = {
        "source": source.name,
        "source_sha256": digest,
        "branch": meta["branch"],
        "scope": ("exact restrictions of the saved p43=0 depth-14 necessary ideal to all 84 non-pinned active parameter coordinate axes; not full-branch exclusion"
                  if all_parameters else "exact restrictions of the saved p43=0 depth-14 necessary ideal to six coordinate axes; not full-branch exclusion"),
        "omitted_pinned_zero_axes": list(PINNED_ZERO_PARAMETERS) if all_parameters else [],
        "equations": equations,
        "sparse_terms": terms_seen,
        "axis_count": len(results),
        "axes": results,
        "axis_status_counts": {
            status: sum(entry["status"] == status for entry in results.values())
            for status in sorted({entry["status"] for entry in results.values()})
        },
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
