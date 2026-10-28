"""Apply the exact distance-four linear graph to every reduced F3 constraint."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import f2_75_125_quadratic_affine_profile as core
from f2_75_125_graph_distance4_linear_eliminate import profile as audit_graph
from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


TARGET_TAGS = {
    *(f"terminal-alpha-P-distance4:u^{i}" for i in range(5)),
    *(f"terminal-alpha-Q-distance4:u^{i}" for i in range(4, 9)),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def e(rational, radical=0):
    return Fraction(rational), Fraction(radical)


def add_poly(left, right):
    result = dict(left)
    for monomial, coefficient in right.items():
        updated = core.add(result.get(monomial, core.ZERO), coefficient)
        if core.iszero(updated):
            result.pop(monomial, None)
        else:
            result[monomial] = updated
    return result


def multiply_poly(left, right):
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


def make_graph(names):
    size = len(names)
    index = {name: names.index(name) for name in
             ("p0", "p3", "p8", "p10", "p11", "p12", "p13", "p14", "p15", "p16")}
    eliminated = tuple(index[f"p{i}"] for i in range(10, 15))
    zero = (0,) * size

    def term(coefficient, **powers_by_name):
        powers = [0] * size
        for name, exponent in powers_by_name.items():
            powers[index[name]] = exponent
        return tuple(powers), coefficient

    # These are the five pivot equations from the independently replayed
    # exact RREF, solved for p10,...,p14. U,V,W,Z stand for the displayed
    # monomials in p0,p3,p8.
    replacement = {
        "p10": dict([
            term(e(Fraction(-3875, 2), Fraction(-1705, 2)), p15=1),
            term(e(-32310, -14360), p16=1),
            term(e(Fraction(-2795, 110592), Fraction(6149, 552960)), p0=4),
            term(e(Fraction(955, 192)), p0=2, p3=1),
            term(e(Fraction(-20, 3), Fraction(-8, 3)), p0=1, p8=1),
            term(e(Fraction(-1825, 8), Fraction(-803, 8)), p3=2),
        ]),
        "p11": dict([
            term(e(Fraction(1995, 2), Fraction(855, 2)), p15=1),
            term(e(15875, 6985), p16=1),
            term(e(Fraction(209, 12288), Fraction(-209, 27648)), p0=4),
            term(e(Fraction(-955, 288), Fraction(191, 288)), p0=2, p3=1),
            term(e(Fraction(17, 8), Fraction(17, 24)), p0=1, p8=1),
            term(e(Fraction(259, 3), 37), p3=2),
        ]),
        "p12": dict([
            term(e(-350, -140), p15=1),
            term(e(Fraction(-10255, 2), Fraction(-4395, 2)), p16=1),
            term(e(Fraction(-715, 110592), Fraction(319, 110592)), p0=4),
            term(e(Fraction(107, 64), Fraction(-107, 192)), p0=2, p3=1),
            term(e(Fraction(-5, 12), Fraction(-1, 12)), p0=1, p8=1),
            term(e(Fraction(-125, 6), Fraction(-25, 3)), p3=2),
        ]),
        "p13": dict([
            term(e(Fraction(165, 2), Fraction(55, 2)), p15=1),
            term(e(1025, 410), p16=1),
            term(e(Fraction(235, 221184), Fraction(-35, 73728)), p0=4),
            term(e(Fraction(-103, 192), Fraction(103, 480)), p0=2, p3=1),
            term(e(Fraction(1, 24)), p0=1, p8=1),
            term(e(3, 1), p3=2),
        ]),
        "p14": dict([
            term(e(Fraction(-25, 2), Fraction(-5, 2)), p15=1),
            term(e(-105, -35), p16=1),
            term(e(Fraction(91, 1152), Fraction(-13, 384)), p0=2, p3=1),
            term(e(Fraction(-5, 24), Fraction(-1, 24)), p3=2),
        ]),
    }
    if tuple(index[name] for name in replacement) != eliminated:
        raise AssertionError("unexpected generator order")
    if any(zero in polynomial for polynomial in replacement.values()):
        raise AssertionError("graph replacement unexpectedly has a constant term")
    return eliminated, replacement, zero


def powers_of(base, maximum, one):
    values = [one]
    for _ in range(maximum):
        values.append(multiply_poly(values[-1], base))
    return values


def apply_one(saved, size, eliminated, powers, cache):
    result = {}
    for sparse, coefficient_values in saved:
        exponents = [0] * len(eliminated)
        original = [0] * size
        for index, exponent in sparse:
            if index in eliminated:
                exponents[eliminated.index(index)] = exponent
            else:
                original[index] = exponent
        key = tuple(exponents)
        if key not in cache:
            expanded = {((0,) * size): core.ONE}
            for i, exponent in enumerate(key):
                if exponent:
                    expanded = multiply_poly(expanded, powers[i][exponent])
            cache[key] = expanded
        coefficient = core.from_json(coefficient_values)
        for shift, shift_coefficient in cache[key].items():
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
    stream_profile = json.loads(args.stream_profile.read_text(encoding="utf-8"))
    linear_receipt = json.loads(args.linear_elimination.read_text(encoding="utf-8"))
    if stream_profile.get("source_sha256") != source_hash:
        raise ValueError("stream profile does not match source constraints")
    if linear_receipt.get("checkpoint_sha256") != source_hash:
        raise ValueError("linear-elimination receipt does not match source constraints")
    replay = audit_graph(args.checkpoint, args.stream_profile)
    for key in ("normalized_row_digest", "rref", "pivot_columns", "rank"):
        if replay.get(key) != linear_receipt.get(key):
            raise ValueError(f"linear-elimination receipt failed replay at {key}")
    if replay["rank"] != 5 or replay["pivot_columns"] != [0, 1, 2, 3, 4]:
        raise ValueError("expected the audited rank-five distance-four graph")

    text = args.checkpoint.read_text(encoding="utf-8")
    metadata, constraint_pos = top_metadata_and_constraints(text)
    names = metadata["generators"]
    eliminated, replacement, zero = make_graph(names)
    maxima = [0] * len(eliminated)
    for _tag, saved in iter_json_array(text, constraint_pos):
        for sparse, _coefficient in saved:
            for index, exponent in sparse:
                if index in eliminated:
                    maxima[eliminated.index(index)] = max(
                        maxima[eliminated.index(index)], exponent)
    powers = [
        powers_of(replacement[name], maximum, {zero: core.ONE})
        for name, maximum in zip(replacement, maxima)
    ]
    cache = {}
    output_data = {
        "format": "f2-sparse-ring-v1",
        "branch": metadata["branch"],
        "quadratic": True,
        "generators": names,
        "nparameters": int(metadata["nparameters"]),
        "compaction": {
            "source_checkpoint_sha256": source_hash,
            "component": "F3=0",
            "linear_graph_digest": linear_receipt["normalized_row_digest"],
            "method": "exact substitution of the rank-five distance-four linear graph",
            "eliminated_parameters": [f"p{i}" for i in range(10, 15)],
            "scope": "constraint-only necessary system restricted to the verified exact graph",
        },
        "constraints": [],
    }
    constraint_count = source_term_count = output_term_count = zero_count = 0
    zero_tags = set()
    with args.output.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps({k: v for k, v in output_data.items() if k != "constraints"},
                                separators=(",", ":"))[:-1])
        stream.write(',"constraints":[')
        first = True
        for tag, saved in iter_json_array(text, constraint_pos):
            constraint_count += 1
            source_term_count += len(saved)
            transformed = apply_one(saved, len(names), eliminated, powers, cache)
            output_term_count += len(transformed)
            if not transformed:
                zero_count += 1
                zero_tags.add(tag)
            if tag in TARGET_TAGS and transformed:
                raise AssertionError(f"distance-four graph failed to kill {tag}")
            if any(any(monomial[index] for index in eliminated)
                   for monomial in transformed):
                raise AssertionError(f"an eliminated generator reappeared in {tag}")
            encoded = [
                [[[i, power] for i, power in enumerate(monomial) if power],
                 encode_coefficient(coefficient)]
                for monomial, coefficient in sorted(transformed.items())
            ]
            if not first:
                stream.write(",")
            json.dump([tag, encoded], stream, separators=(",", ":"))
            first = False
        stream.write("]}\n")
    if not TARGET_TAGS.issubset(zero_tags):
        raise AssertionError("not all ten distance-four rows vanished")
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
        "zero_constraints_after_graph": zero_count,
        "zero_constraint_tags_after_graph": sorted(zero_tags),
        "distance_four_rows_vanishing": sorted(TARGET_TAGS),
        "maximum_exponents_eliminated": {
            f"p{i}": exponent for i, exponent in zip(range(10, 15), maxima)
        },
        "polynomial_power_cache_entries": len(cache),
        "output_checkpoint": str(args.output),
        "output_checkpoint_sha256": output_hash,
        "scope": (
            "Exact necessary-system graph restriction. This is not an exclusion "
            "unless the transformed ideal is shown to be unit and the reduction "
            "is attached to the actual hypothetical pair."
        ),
    }
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: receipt[key] for key in (
        "component", "source_constraint_count", "source_constraint_terms",
        "transformed_constraint_terms", "zero_constraints_after_graph",
        "distance_four_rows_vanishing", "maximum_exponents_eliminated",
        "polynomial_power_cache_entries", "output_checkpoint_sha256", "scope",
    )}, indent=2))


if __name__ == "__main__":
    main()
