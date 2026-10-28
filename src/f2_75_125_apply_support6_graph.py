"""Apply the exact support-six F3 graph to every reduced constraint."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import f2_75_125_quadratic_affine_profile as core
from f2_75_125_graph_support6_linear_eliminate import profile as audit_support6
from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


EXPECTED_ZERO_TAGS = {
    "terminal-alpha-P-distance3:u^2",
    "terminal-alpha-P-distance3:u^3",
    "terminal-alpha-P-distance3:u^4",
    "terminal-alpha-Q-distance3:u^6",
    "terminal-alpha-Q-distance3:u^7",
    "terminal-alpha-Q-distance3:u^8",
    "terminal-alpha-Q-distance6:u^2",
    "terminal-alpha-Q-distance6:u^3",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def e(rational, radical=0):
    return Fraction(rational), Fraction(radical)


def add_polys(left, right):
    result = dict(left)
    for monomial, coefficient in right.items():
        updated = core.add(result.get(monomial, core.ZERO), coefficient)
        if core.iszero(updated):
            result.pop(monomial, None)
        else:
            result[monomial] = updated
    return result


def multiply_polys(left, right):
    result = {}
    for powers_a, coefficient_a in left.items():
        for powers_b, coefficient_b in right.items():
            monomial = tuple(a + b for a, b in zip(powers_a, powers_b))
            value = core.mul(coefficient_a, coefficient_b)
            updated = core.add(result.get(monomial, core.ZERO), value)
            if core.iszero(updated):
                result.pop(monomial, None)
            else:
                result[monomial] = updated
    return result


def encode_coefficient(value):
    rational, radical = value
    if not radical:
        return [str(rational)]
    return [str(radical), str(rational)]


def row_graph(names):
    size = len(names)
    indices = {name: names.index(name) for name in ("p0", "p3", "p5", "p6", "p7", "p8")}
    zero = (0,) * size

    def term(coefficient, **powers_by_name):
        powers = [0] * size
        for name, exponent in powers_by_name.items():
            powers[indices[name]] = exponent
        return tuple(powers), coefficient

    replacement = {
        "p5": dict([
            term(e(-460, -184), p8=1),
            term(e(Fraction(-25, 288), Fraction(11, 288)), p0=3),
            term(e(Fraction(83, 6)), p0=1, p3=1),
        ]),
        "p6": dict([
            term(e(Fraction(195, 2), Fraction(65, 2)), p8=1),
            term(e(Fraction(5, 384), Fraction(-5, 864)), p0=3),
            term(e(Fraction(-145, 48), Fraction(29, 48)), p0=1, p3=1),
        ]),
        "p7": dict([
            term(e(Fraction(-27, 2), Fraction(-27, 10)), p8=1),
            term(e(Fraction(7, 20), Fraction(-7, 60)), p0=1, p3=1),
        ]),
    }
    for name, polynomial in replacement.items():
        if any(not any(monomial) for monomial in polynomial):
            raise AssertionError(f"unexpected constant in {name} graph")
    return indices, replacement, zero


def polynomial_powers(base, maximum, one):
    values = [one]
    for _ in range(maximum):
        values.append(multiply_polys(values[-1], base))
    return values


def apply_one(saved, size, eliminated_indices, powers_by_exponent, product_cache):
    result = {}
    for sparse, coefficient_values in saved:
        powers = [0] * size
        for index, exponent in sparse:
            powers[index] = exponent
        exponents = tuple(powers[index] for index in eliminated_indices)
        if exponents not in product_cache:
            a, b, c = exponents
            expanded = multiply_polys(powers_by_exponent[0][a],
                                      powers_by_exponent[1][b])
            expanded = multiply_polys(expanded, powers_by_exponent[2][c])
            product_cache[exponents] = expanded
        expanded = product_cache[exponents]
        original = list(powers)
        for index in eliminated_indices:
            original[index] = 0
        coefficient = core.from_json(coefficient_values)
        for shift, shift_coefficient in expanded.items():
            monomial = tuple(a + b for a, b in zip(original, shift))
            value = core.mul(coefficient, shift_coefficient)
            updated = core.add(result.get(monomial, core.ZERO), value)
            if core.iszero(updated):
                result.pop(monomial, None)
            else:
                result[monomial] = updated
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("stream_profile", type=Path)
    parser.add_argument("linear_elimination", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    for path in (args.output, args.receipt):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite {path}")

    source_hash = sha256(args.checkpoint)
    source_profile = json.loads(args.stream_profile.read_text(encoding="utf-8"))
    linear_receipt = json.loads(args.linear_elimination.read_text(encoding="utf-8"))
    if source_profile.get("source_sha256") != source_hash:
        raise ValueError("stream profile does not match source constraints")
    if linear_receipt.get("checkpoint_sha256") != source_hash:
        raise ValueError("linear-elimination receipt does not match source constraints")
    replay = audit_support6(args.checkpoint, args.stream_profile)
    for key in ("normalized_row_digest", "rref", "pivot_columns", "rank"):
        if replay.get(key) != linear_receipt.get(key):
            raise ValueError(f"linear-elimination receipt failed exact replay at {key}")
    if replay["rank"] != 3 or replay["pivot_columns"] != [0, 1, 2]:
        raise ValueError("expected the audited rank-three support-six graph")

    text = args.checkpoint.read_text(encoding="utf-8")
    metadata, constraint_pos = top_metadata_and_constraints(text)
    names = metadata["generators"]
    indices, replacement, zero = row_graph(names)
    eliminated_indices = tuple(indices[name] for name in ("p5", "p6", "p7"))
    maxima = [0, 0, 0]
    for _tag, saved in iter_json_array(text, constraint_pos):
        for sparse, _values in saved:
            for j, generator in enumerate(eliminated_indices):
                exponent = next((power for index, power in sparse
                                 if index == generator), 0)
                maxima[j] = max(maxima[j], exponent)
    powers_by_exponent = [
        polynomial_powers(replacement[name], maximum, {zero: core.ONE})
        for name, maximum in zip(("p5", "p6", "p7"), maxima)
    ]
    product_cache = {}

    data = {
        "format": "f2-sparse-ring-v1",
        "branch": metadata["branch"],
        "quadratic": True,
        "generators": names,
        "nparameters": int(metadata["nparameters"]),
        "constraints": [],
        "compaction": {
            "source_checkpoint_sha256": source_hash,
            "component": "F3=0",
            "support_six_row_digest": linear_receipt["normalized_row_digest"],
            "method": "exact graph substitution from the rank-three support-six row reduction",
            "eliminated_parameters": ["p5", "p6", "p7"],
            "scope": "constraint-only necessary system restricted to the verified support-six graph",
        },
    }

    constraint_count = 0
    source_term_count = 0
    output_term_count = 0
    zero_constraints = 0
    zero_tags = set()
    with args.output.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("{")
        first_key = True
        for key, value in data.items():
            if key == "constraints":
                continue
            if not first_key:
                stream.write(",")
            stream.write(json.dumps(key, separators=(",", ":")))
            stream.write(":")
            stream.write(json.dumps(value, separators=(",", ":")))
            first_key = False
        stream.write(',"constraints":[')
        first_constraint = True
        for tag, saved in iter_json_array(text, constraint_pos):
            constraint_count += 1
            source_term_count += len(saved)
            transformed = apply_one(saved, len(names), eliminated_indices,
                                    powers_by_exponent, product_cache)
            output_term_count += len(transformed)
            if not transformed:
                zero_constraints += 1
                zero_tags.add(tag)
            if tag in EXPECTED_ZERO_TAGS and transformed:
                raise AssertionError(f"support-six graph failed to kill {tag}")
            for monomial, _coefficient in transformed.items():
                if any(monomial[index] for index in eliminated_indices):
                    raise AssertionError(f"eliminated parameter reappeared in {tag}")
            encoded = [
                [
                    [[index, exponent] for index, exponent in enumerate(monomial)
                     if exponent],
                    encode_coefficient(coefficient),
                ]
                for monomial, coefficient in sorted(transformed.items())
            ]
            if not first_constraint:
                stream.write(",")
            json.dump([tag, encoded], stream, separators=(",", ":"))
            first_constraint = False
        stream.write("]}\n")

    support_six_zero_tags = zero_tags & EXPECTED_ZERO_TAGS
    if support_six_zero_tags != EXPECTED_ZERO_TAGS:
        raise AssertionError("not all eight support-six source rows vanished")
    output_hash = sha256(args.output)
    receipt = {
        "source_checkpoint": str(args.checkpoint),
        "source_checkpoint_sha256": source_hash,
        "source_stream_profile": str(args.stream_profile),
        "source_linear_elimination": str(args.linear_elimination),
        "component": "F3=0",
        "branch": metadata["branch"],
        "quadratic_field": "Q(sqrt(5))",
        "source_constraint_count": constraint_count,
        "source_constraint_terms": source_term_count,
        "transformed_constraint_terms": output_term_count,
        "zero_constraints_after_graph": zero_constraints,
        "zero_constraint_tags_after_graph": sorted(zero_tags),
        "support_six_rows_vanishing": sorted(support_six_zero_tags),
        "maximum_exponents_eliminated": dict(zip(("p5", "p6", "p7"), maxima)),
        "polynomial_power_cache_entries": len(product_cache),
        "output_checkpoint": str(args.output),
        "output_checkpoint_sha256": output_hash,
        "graph_substitution": {
            "p5": "-(460+184*sqrt(5))*p8-(25/288-11*sqrt(5)/288)*p0^3+(83/6)*p0*p3",
            "p6": "(195/2+65*sqrt(5)/2)*p8+(5/384-5*sqrt(5)/864)*p0^3-(145/48-29*sqrt(5)/48)*p0*p3",
            "p7": "-(27/2+27*sqrt(5)/10)*p8+(7/20-7*sqrt(5)/60)*p0*p3",
            "rref_variables": list(linear_receipt["variables"]),
        },
        "scope": (
            "Exact equivalent graph restriction of all saved F3 necessary constraints. "
            "It is not an exclusion unless the reduced ideal is proved unit."
        ),
    }
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: receipt[key] for key in (
        "component", "source_constraint_count", "source_constraint_terms",
        "transformed_constraint_terms", "zero_constraints_after_graph",
        "support_six_rows_vanishing", "maximum_exponents_eliminated",
        "polynomial_power_cache_entries", "output_checkpoint_sha256", "scope",
    )}, indent=2))


if __name__ == "__main__":
    main()
