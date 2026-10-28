"""Stream exact F3/F4 constraints through the shared p1,p2 graph substitution.

This checks the complete saved necessary-equation list for an immediate
constant or affine unit-ideal obstruction without rewriting the checkpoint's
P/Q/W sequence polynomials.  It is an exact necessary-system audit, not a
branch exclusion unless an exact unit certificate is produced.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import f2_75_125_quadratic_affine_profile as core
from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    top_metadata_and_constraints,
)


ZERO = core.ZERO
ONE = core.ONE


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def neg(value: core.Elt) -> core.Elt:
    return -value[0], -value[1]


def scalar(value: int) -> core.Elt:
    return Fraction(value), Fraction(0)


def scale(value: core.Elt, rational: Fraction) -> core.Elt:
    return value[0] * rational, value[1] * rational


def power_table(value: core.Elt, maximum: int) -> list[core.Elt]:
    result = [ONE]
    for _ in range(maximum):
        result.append(core.mul(result[-1], value))
    return result


def binomial_expansion(exponent: int, alpha: core.Elt, beta: core.Elt):
    """Terms (number of alpha factors, number of beta factors, coefficient)."""
    alpha_powers = power_table(alpha, exponent)
    beta_powers = power_table(beta, exponent)
    return [
        (i, exponent - i,
         scale(core.mul(alpha_powers[i], beta_powers[exponent - i]),
               Fraction(math.comb(exponent, i))))
        for i in range(exponent + 1)
    ]


def graph_coefficients():
    # p2=(p0^2-B*p3)/A; p1=-C*p2-D*p3, in Q(sqrt(5)).
    A = (Fraction(672, 5), Fraction(288, 5))
    B = (Fraction(1728), Fraction(19008, 25))
    C = (Fraction(8), Fraction(8, 5))
    D = (Fraction(459, 10), Fraction(153, 10))
    inv_A = core.inv(A)
    alpha2 = inv_A
    beta2 = core.mul(neg(B), inv_A)
    alpha1 = core.mul(neg(C), inv_A)
    beta1 = core.sub(core.mul(core.mul(C, B), inv_A), D)
    if not core.iszero(core.sub(core.mul(A, alpha2), ONE)):
        raise AssertionError("p2 graph does not solve the normalized g1 relation")
    if not core.iszero(core.add(core.mul(A, beta2), B)):
        raise AssertionError("p2 graph has an incorrect p3 coefficient")
    if not core.iszero(core.add(alpha1, core.mul(C, alpha2))):
        raise AssertionError("p1 graph has an incorrect p0^2 coefficient")
    if not core.iszero(core.add(core.add(beta1, core.mul(C, beta2)), D)):
        raise AssertionError("p1 graph has an incorrect p3 coefficient")
    return alpha1, beta1, alpha2, beta2


def substitute_constraint(saved, generator_count: int, p_index: dict[int, int],
                          alphas_betas):
    alpha1, beta1, alpha2, beta2 = alphas_betas
    p0_gen, p1_gen, p2_gen, p3_gen = (
        p_index[0], p_index[1], p_index[2], p_index[3]
    )
    result: dict[tuple[int, ...], core.Elt] = {}
    for sparse, values in saved:
        powers = [0] * generator_count
        for index, exponent in sparse:
            powers[index] = exponent
        e1, e2 = powers[p1_gen], powers[p2_gen]
        coefficient = core.from_json(values)
        expansion1 = binomial_expansion(e1, alpha1, beta1)
        expansion2 = binomial_expansion(e2, alpha2, beta2)
        powers[p1_gen] = powers[p2_gen] = 0
        for s1, b1, c1 in expansion1:
            for s2, b2, c2 in expansion2:
                monomial = list(powers)
                monomial[p0_gen] += 2 * (s1 + s2)
                monomial[p3_gen] += b1 + b2
                key = tuple(monomial)
                value = core.mul(coefficient, core.mul(c1, c2))
                updated = core.add(result.get(key, ZERO), value)
                if core.iszero(updated):
                    result.pop(key, None)
                else:
                    result[key] = updated
    return result


def encode_coefficient(value: core.Elt) -> list[str]:
    rational, radical = value
    if not radical:
        return [str(rational)]
    # Match AlgebraicField.to_list(): descending powers of sqrt(5).
    return [str(radical), str(rational)]


def profile(source: Path, source_receipt: Path, stream_profile: Path) -> dict:
    source_hash = sha256(source)
    provenance = json.loads(source_receipt.read_text(encoding="utf-8"))
    audited_profile = json.loads(stream_profile.read_text(encoding="utf-8"))
    if audited_profile.get("source_sha256") != source_hash:
        raise ValueError("stream profile does not match the checkpoint bytes")
    component = provenance.get("component")
    if component not in ("F3=0", "F4=0"):
        raise ValueError("expected a verified F3=0 or F4=0 component receipt")

    text = source.read_text(encoding="utf-8")
    metadata, constraint_pos = top_metadata_and_constraints(text)
    if (metadata.get("format") != "f2-sparse-ring-v1"
            or metadata.get("branch") != "2+2"
            or not metadata.get("quadratic")):
        raise ValueError("expected a quadratic-field F2 2+2 checkpoint")
    names = metadata["generators"]
    nparameters = int(metadata["nparameters"])
    p_index = {i: names.index(f"p{i}") for i in range(nparameters)}
    parameter_order = {index: i for i, index in p_index.items()}
    graph = graph_coefficients()
    expected_zero_tags = {
        "terminal-alpha-P-distance2:u^3",
        "terminal-alpha-P-distance2:u^4",
        "terminal-alpha-Q-distance2:u^7",
        "terminal-alpha-Q-distance2:u^8",
        "terminal-alpha-Q-distance2:u^9",
    }
    if component == "F4=0":
        expected_zero_tags |= {
            "terminal-alpha-P-distance4:u^0",
            "terminal-alpha-Q-distance4:u^4",
        }
    constraint_count = 0
    source_terms = 0
    transformed_terms = 0
    zero_constraints = 0
    constants = []
    affine_rows = []
    affine_tags = []
    degree_histogram = Counter()
    support_histogram = Counter()
    digest = hashlib.sha256()
    zero_tags = set()
    reduced_constraints = []

    for tag, saved in iter_json_array(text, constraint_pos):
        constraint_count += 1
        source_terms += len(saved)
        transformed = substitute_constraint(saved, len(names), p_index, graph)
        transformed_terms += len(transformed)
        if not transformed:
            zero_constraints += 1
            if tag in expected_zero_tags:
                zero_tags.add(tag)
        elif tag in expected_zero_tags:
            raise AssertionError(f"verified support-four row {tag} did not vanish")
        digest.update(tag.encode("utf-8"))
        for monomial, coefficient in sorted(transformed.items()):
            digest.update(repr((monomial, coefficient)).encode("utf-8"))
        reduced_constraints.append([
            tag,
            [
                [
                    [[index, exponent] for index, exponent in enumerate(monomial)
                     if exponent],
                    encode_coefficient(coefficient),
                ]
                for monomial, coefficient in sorted(transformed.items())
            ],
        ])

        if len(transformed) == 1:
            (monomial, coefficient), = transformed.items()
            if not any(monomial):
                constants.append({"tag": tag, "value": core.show(coefficient)})

        row = [ZERO] * (nparameters + 1)
        affine = True
        support = set()
        maximum_degree = 0
        for monomial, coefficient in transformed.items():
            degree = 0
            variable = None
            for generator_index, exponent in enumerate(monomial):
                if not exponent:
                    continue
                ordinal = parameter_order.get(generator_index)
                if ordinal is None:
                    affine = False
                    continue
                degree += exponent
                support.add(ordinal)
                if variable is None:
                    variable = ordinal
                elif exponent:
                    variable = -1
            maximum_degree = max(maximum_degree, degree)
            if degree > 1 or variable == -1:
                affine = False
            elif degree == 0:
                row[0] = core.add(row[0], coefficient)
            else:
                row[variable + 1] = core.add(row[variable + 1], coefficient)

        degree_histogram[maximum_degree] += 1
        support_histogram[len(support)] += 1
        if affine:
            affine_rows.append(row)
            affine_tags.append(tag)

    if zero_tags != expected_zero_tags:
        raise AssertionError("not all audited support-four rows were replayed")

    rank, pivots, inconsistency = core.rank_affine(affine_rows)
    affine_certificate = None
    if inconsistency is not None:
        weights = inconsistency["source_row_combination"]
        combined = [ZERO] * (nparameters + 1)
        terms = []
        for index, coefficient in sorted(weights.items()):
            combined = [core.add(value, core.mul(coefficient, term))
                        for value, term in zip(combined, affine_rows[index])]
            terms.append({"tag": affine_tags[index],
                          "coefficient": core.show(coefficient)})
        if combined[0] != ONE or any(not core.iszero(x) for x in combined[1:]):
            raise AssertionError("affine certificate replay did not give 1")
        affine_certificate = {
            "identity": "1 = sum(coefficient_i * substituted_constraint_i)",
            "terms": terms,
            "exact_replay": True,
        }

    return {
        "source_checkpoint": str(source),
        "source_checkpoint_sha256": source_hash,
        "source_receipt": str(source_receipt),
        "stream_profile": str(stream_profile),
        "component": component,
        "branch": metadata["branch"],
        "quadratic_field": "Q(sqrt(5))",
        "constraint_count": constraint_count,
        "source_constraint_terms": source_terms,
        "substituted_constraint_terms": transformed_terms,
        "substituted_constraint_digest": digest.hexdigest(),
        "identically_zero_constraints": zero_constraints,
        "verified_support_four_rows_vanishing": sorted(zero_tags),
        "nonzero_constant_obstructions": constants,
        "maximum_parameter_degree_histogram": dict(sorted(degree_histogram.items())),
        "parameter_support_histogram": dict(sorted(support_histogram.items())),
        "affine_constraint_count": len(affine_rows),
        "affine_rank_over_Qsqrt5": rank,
        "affine_system_inconsistent": affine_certificate is not None,
        "affine_unit_certificate": affine_certificate,
        "affine_pivots": pivots,
        "reduced_checkpoint_data": {
            "format": "f2-sparse-ring-v1",
            "branch": metadata["branch"],
            "quadratic": True,
            "generators": names,
            "nparameters": nparameters,
            "constraints": reduced_constraints,
            "compaction": {
                "source_checkpoint_sha256": source_hash,
                "component": component,
                "basis_relation_digest": digest.hexdigest(),
                "method": "exact substitution of p1,p2 from support-four graph basis",
                "scope": "constraint-only equivalent restriction to the verified graph",
            },
        },
        "substitution": {
            "p1": "-(8+8*sqrt(5)/5)*p2-(459/10+153*sqrt(5)/10)*p3",
            "p2": "(p0^2-(1728+19008*sqrt(5)/25)*p3)/(672/5+288*sqrt(5)/5)",
            "basis_source": "exact replay of the two support-four Groebner generators",
        },
        "scope": (
            "All saved necessary constraints were substituted exactly. A nonzero "
            "constant or replayed affine unit identity excludes this component; "
            "otherwise the report is diagnostic and does not exclude it."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("stream_profile", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reduced-checkpoint", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    if args.reduced_checkpoint and args.reduced_checkpoint.exists():
        raise FileExistsError(f"refusing to overwrite {args.reduced_checkpoint}")
    result = profile(args.checkpoint, args.receipt, args.stream_profile)
    reduced_data = result.pop("reduced_checkpoint_data")
    if args.reduced_checkpoint:
        raw = (json.dumps(reduced_data, separators=(",", ":")) + "\n").encode("utf-8")
        args.reduced_checkpoint.write_bytes(raw)
        result["reduced_checkpoint"] = str(args.reduced_checkpoint)
        result["reduced_checkpoint_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "component", "constraint_count", "source_constraint_terms",
        "substituted_constraint_terms", "identically_zero_constraints",
        "nonzero_constant_obstructions", "maximum_parameter_degree_histogram",
        "parameter_support_histogram", "affine_constraint_count",
        "affine_rank_over_Qsqrt5", "affine_system_inconsistent",
        "affine_unit_certificate", "output",
    ) if key in result}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
