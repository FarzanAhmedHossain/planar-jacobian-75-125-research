"""Stream the d19 checkpoint and impose exact simple-root divisibilities.

The local analysis gives T | U2,U3 and S | U4,U5 on the U1=0,
2+1+1 branch, where S=t^2-3t+3 and T=(t-1)^3*S^2. The saved P bands
are p[2]=t^3*U2, p[3]=t^3*U3, p[4]=t^3*U4, p[5]=t^2*U5. Since t is
coprime to S*T, the equivalent exact necessary equations are obtained by
reducing p[2],p[3] modulo T and p[4],p[5] modulo S.

Only those four P-band polynomials are decoded. The large constraints, Q
bands, and W bands are streamed past, not retained in memory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path


DEFAULT_SOURCE = (
    "research/f2_75_125_sparse_ring_lift_211_roots19_u1zero_"
    "streamed_state.json"
)
DEFAULT_OUTPUT = (
    "research/f2_75_125_211_simple_root_divisibility_equations_20261025.json"
)
EXPECTED_SOURCE_SHA256 = (
    "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"
)


class ProgressReader:
    def __init__(self, handle, total_bytes):
        self.handle = handle
        self.total_bytes = total_bytes
        self.bytes_read = 0
        self.next_mark = 128 * 1024 * 1024
        self.started = time.monotonic()

    def read(self, size=-1):
        data = self.handle.read(size)
        self.bytes_read += len(data)
        if self.bytes_read >= self.next_mark:
            print(
                f"stream_bytes={self.bytes_read} total_bytes={self.total_bytes} "
                f"elapsed={time.monotonic()-self.started:.1f}s",
                flush=True,
            )
            self.next_mark += 128 * 1024 * 1024
        return data


def read_header(path):
    import ijson

    header = {}
    with path.open("rb") as handle:
        for key in ("format", "branch", "quadratic", "generators"):
            with path.open("rb") as source:
                values = ijson.items(source, key)
                header[key] = next(values)
    return header


def encode_equation(polynomial):
    terms = []
    for powers, coefficient in sorted(polynomial.items()):
        sparse = [[index, int(exponent)]
                  for index, exponent in enumerate(powers) if exponent]
        terms.append([sparse, str(coefficient)])
    return terms


def component_equations(ring, remainder):
    groups = {}
    for powers, coefficient in remainder.items():
        degree = int(powers[0])
        parameter_monomial = (0, *powers[1:])
        groups.setdefault(degree, {})[parameter_monomial] = coefficient
    equations = []
    reconstructed = ring.zero
    for degree in sorted(groups):
        equation = ring.from_prevalidated_dict(groups[degree])
        equations.append({
            "t_degree": degree,
            "term_count": len(equation),
            "terms": encode_equation(equation),
        })
        reconstructed += equation * ring.gens[0]**degree
    if reconstructed != remainder:
        raise AssertionError("remainder coefficient decomposition failed")
    return equations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", default=DEFAULT_SOURCE)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--allow-unverified-source", action="store_true")
    args = parser.parse_args()

    source = Path(args.source)
    output = Path(args.output)
    source_bytes = source.stat().st_size
    if source_bytes != 1_009_425_753:
        raise ValueError(f"unexpected d19 checkpoint size: {source_bytes}")
    if EXPECTED_SOURCE_SHA256 is None and not args.allow_unverified_source:
        raise ValueError("set an expected source digest before running")

    header = read_header(source)
    if header["format"] != "f2-sparse-ring-v1":
        raise ValueError(f"unexpected format: {header['format']}")
    if header["branch"] != "2+1+1" or header["quadratic"]:
        raise ValueError("expected the rational 2+1+1 d19 checkpoint")
    if header["generators"][0] != "t":
        raise ValueError("expected t to be the first generator")

    research_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(research_dir))
    from f2_75_125_flint_mpoly_backend import FlintRing
    from flint import fmpq

    ring = FlintRing(header["generators"])
    t = ring.gens[0]
    S = t**2 - 3*t + 3
    T = (t-1)**3 * S**2
    needed = {2, 3, 4, 5}
    p_bands = {}
    started = time.monotonic()
    with source.open("rb") as handle:
        reader = ProgressReader(handle, source_bytes)
        import ijson
        for index, rows in enumerate(ijson.items(reader, "p.item")):
            if index not in needed:
                continue
            terms = {}
            for sparse, values in rows:
                powers = [0] * len(header["generators"])
                for variable_index, exponent in sparse:
                    powers[int(variable_index)] = int(exponent)
                terms[tuple(powers)] = fmpq(values[0])
            p_bands[index] = ring.from_prevalidated_dict(terms)
            print(
                f"decoded_p{index}_terms={len(p_bands[index])} "
                f"min_t_order={min((powers[0] for powers in terms), default='inf')}",
                flush=True,
            )
            if needed.issubset(p_bands):
                break

    if set(p_bands) != needed:
        raise ValueError(f"checkpoint had only p indices {sorted(p_bands)}")

    modulus_by_band = {2: T, 3: T, 4: S, 5: S}
    results = {}
    for index in sorted(needed):
        band = p_bands[index]
        remainder = band.rem(modulus_by_band[index])
        equations = component_equations(ring, remainder)
        results[f"p{index}_mod_{'T' if index in (2, 3) else 'S'}"] = {
            "band_index": index,
            "modulus_degree": 7 if index in (2, 3) else 2,
            "band_term_count": len(band),
            "remainder_term_count": len(remainder),
            "coefficient_equation_count": len(equations),
            "coefficient_equations": equations,
        }
        print(
            f"p{index}_remainder_terms={len(remainder)} "
            f"coefficient_equations={len(equations)}",
            flush=True,
        )

    report = {
        "scope": "exact necessary simple-root divisibilities for F2 2+1+1 U1=0; not branch exclusion",
        "source": str(source),
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "source_bytes": source_bytes,
        "source_format": header["format"],
        "branch": header["branch"],
        "field": "QQ",
        "generator_count": len(header["generators"]),
        "generator_names": header["generators"],
        "divisors": {
            "S": "t^2-3*t+3",
            "T": "(t-1)^3*(t^2-3*t+3)^2",
        },
        "band_identifications": {
            "p2": "t^3*U2",
            "p3": "t^3*U3",
            "p4": "t^3*U4",
            "p5": "t^2*U5",
        },
        "equations": results,
        "elapsed_seconds": time.monotonic()-started,
        "whole_2+1+1_branch_excluded_percent": 0.0,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(f"output={output} script_sha256={script_hash}", flush=True)


if __name__ == "__main__":
    main()
