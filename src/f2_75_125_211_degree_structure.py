#!/usr/bin/env python3
"""Low-memory exact degree/factor profile for detached F2 constraints.

This is structural reconnaissance only. It never treats a necessary system
as an exclusion certificate. The parser streams one JSON constraint at a time
so the report can be computed alongside the exact recurrence workers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def streamed_array(path: Path):
    decoder = json.JSONDecoder()
    chunk_size = 1 << 20
    with path.open("r", encoding="utf-8") as stream:
        buffer = stream.read(chunk_size)
        pos = 0
        eof = False

        def refill() -> bool:
            nonlocal buffer, pos, eof
            if eof:
                return False
            chunk = stream.read(chunk_size)
            if not chunk:
                eof = True
                return False
            buffer = buffer[pos:] + chunk
            pos = 0
            return True

        def skip_space() -> None:
            nonlocal pos
            while True:
                while pos < len(buffer) and buffer[pos].isspace():
                    pos += 1
                if pos < len(buffer):
                    return
                if not refill():
                    raise ValueError("unexpected EOF in JSON array")

        skip_space()
        if buffer[pos] != "[":
            raise ValueError("detached constraint root is not an array")
        pos += 1
        while True:
            skip_space()
            if buffer[pos] == "]":
                return
            while True:
                try:
                    item, end = decoder.raw_decode(buffer, pos)
                    break
                except json.JSONDecodeError:
                    if not refill():
                        raise
            yield item
            pos = end
            skip_space()
            if buffer[pos] == ",":
                pos += 1
            elif buffer[pos] == "]":
                return
            else:
                raise ValueError("malformed JSON array separator")
            if pos > chunk_size:
                buffer = buffer[pos:]
                pos = 0


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def profile(path: Path, nparameters: int) -> dict:
    constraints = terms = 0
    max_degree_hist: Counter[int] = Counter()
    min_degree_hist: Counter[int] = Counter()
    degree_span_hist: Counter[str] = Counter()
    homogeneous_degrees: Counter[int] = Counter()
    support_hist: Counter[int] = Counter()
    term_degree_hist: Counter[int] = Counter()
    monomial_factor_variables: Counter[str] = Counter()
    common_factor_everywhere = set(range(nparameters))
    low_degree_variable_occurrences: Counter[str] = Counter()
    nonparameter_indices: Counter[int] = Counter()
    samples_nonhomogeneous = []
    samples_factorable = []
    zero_constraints = 0

    for tag, polynomial in streamed_array(path):
        constraints += 1
        if not polynomial:
            zero_constraints += 1
            continue
        degree_counts: Counter[int] = Counter()
        support = set()
        monomial_gcd: dict[int, int] | None = None
        for sparse, _coefficient in polynomial:
            terms += 1
            exponents = {int(i): int(e) for i, e in sparse if e}
            parameter_exponents = {
                i - 1: e for i, e in exponents.items()
                if 1 <= i <= nparameters
            }
            support.update(parameter_exponents)
            degree = sum(parameter_exponents.values())
            degree_counts[degree] += 1
            term_degree_hist[degree] += 1
            for index in exponents:
                if index == 0 or index > nparameters:
                    nonparameter_indices[index] += 1
            if degree <= 2:
                for i in parameter_exponents:
                    low_degree_variable_occurrences[f"p{i}"] += 1
            if monomial_gcd is None:
                monomial_gcd = dict(parameter_exponents)
            else:
                for i in tuple(monomial_gcd):
                    monomial_gcd[i] = min(monomial_gcd[i], parameter_exponents.get(i, 0))
                    if monomial_gcd[i] == 0:
                        del monomial_gcd[i]

        support_hist[len(support)] += 1
        if not degree_counts:
            continue
        lo, hi = min(degree_counts), max(degree_counts)
        min_degree_hist[lo] += 1
        max_degree_hist[hi] += 1
        degree_span_hist[f"{lo}:{hi}"] += 1
        if lo == hi:
            homogeneous_degrees[lo] += 1
        elif len(samples_nonhomogeneous) < 12:
            samples_nonhomogeneous.append({"tag": tag, "min_degree": lo, "max_degree": hi})

        factors = sorted(monomial_gcd or {})
        for i in factors:
            monomial_factor_variables[f"p{i}"] += 1
        common_factor_everywhere.intersection_update(factors)
        if factors and len(samples_factorable) < 12:
            samples_factorable.append({"tag": tag, "common_monomial_factor": {
                f"p{i}": monomial_gcd[i] for i in factors
            }})

    return {
        "source": str(path),
        "source_sha256": digest(path),
        "scope": "exact degree and monomial-factor profile of detached necessary equations; not an exclusion",
        "nparameters": nparameters,
        "constraint_count": constraints,
        "zero_constraint_count": zero_constraints,
        "term_count": terms,
        "parameter_support_size_histogram": dict(sorted(support_hist.items())),
        "term_parameter_degree_histogram": dict(sorted(term_degree_hist.items())),
        "constraint_min_degree_histogram": dict(sorted(min_degree_hist.items())),
        "constraint_max_degree_histogram": dict(sorted(max_degree_hist.items())),
        "constraint_degree_span_histogram": dict(sorted(degree_span_hist.items())),
        "homogeneous_constraint_degree_histogram": dict(sorted(homogeneous_degrees.items())),
        "variables_dividing_constraints_count": dict(sorted(monomial_factor_variables.items())),
        "variables_dividing_every_nonzero_constraint": sorted(f"p{i}" for i in common_factor_everywhere),
        "low_degree_variable_term_occurrences": dict(sorted(low_degree_variable_occurrences.items())),
        "nonparameter_generator_occurrences": dict(sorted(nonparameter_indices.items())),
        "nonhomogeneous_samples": samples_nonhomogeneous,
        "monomial_factor_samples": samples_factorable,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("nparameters", type=int)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    result = profile(args.source, args.nparameters)
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
