#!/usr/bin/env python3
"""Exact ideal audit on the top-weight boundary of the 2+1+1 slice.

Set the five support variables of weights 1,2,3 to zero and retain the seven
variables of weight 4.  The restriction is a homogeneous ideal in seven
variables generated in degrees two and three.  A zero-dimensional result
with every variable nilpotent would exclude this weighted-projective boundary
stratum, but would not exclude the full 2+1+1 branch.
"""

from __future__ import annotations

import hashlib
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

from f2_75_125_211_degree_structure import streamed_array


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_p43zero_depth13_constraints.jsonfrag"
PROFILE = HERE / "f2_75_125_211_p43zero_depth13_support12_weighted_homogeneity_20261020.json"
RECEIPT = HERE / "f2_75_125_211_support12_topweight_stratum_20261020.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    profile = json.loads(PROFILE.read_text(encoding="utf-8"))
    source_hash = sha256(SOURCE)
    if source_hash != profile["source_sha256"]:
        raise SystemExit("detached source hash does not match weighted-profile report")

    support = tuple(int(name[1:]) for name in profile["support_variables"])
    weights = tuple(int(4 * sp.Rational(value)) for value in profile["positive_weight_vector"])
    if len(support) != 12 or len(weights) != 12:
        raise SystemExit("unexpected support or weight vector")
    high_positions = tuple(i for i, weight in enumerate(weights) if weight == 4)
    low_positions = tuple(i for i, weight in enumerate(weights) if weight < 4)
    high_globals = tuple(support[i] for i in high_positions)
    high_variables = tuple(sp.Symbol(f"p{i}") for i in high_globals)
    high_position = {global_i: local_i for local_i, global_i in enumerate(high_globals)}

    polynomials: list[sp.Expr] = []
    tags: list[str] = []
    weighted_degrees: Counter[int] = Counter()
    ordinary_degrees: Counter[int] = Counter()
    term_count = 0
    for tag, saved in streamed_array(SOURCE):
        active: set[int] = set()
        parsed = []
        terms: defaultdict[tuple[int, ...], sp.Rational] = defaultdict(lambda: sp.S.Zero)
        source_weight_degrees: set[int] = set()
        for sparse, values in saved:
            exponents = {
                index - 1: exponent
                for index, exponent in sparse
                if 1 <= index <= 75 and exponent
            }
            active.update(exponents)
            parsed.append((exponents, values))
        if tuple(sorted(active)) != support:
            continue

        for exponents, values in parsed:
            source_weight_degrees.add(
                sum(weights[support.index(i)] * e for i, e in exponents.items())
            )
            if len(values) != 1:
                raise ValueError("the detached 2+1+1 source should be over QQ")
            if any(exponents.get(support[position], 0) for position in low_positions):
                continue
            vector = [0] * len(high_globals)
            ordinary_degree = 0
            for global_i, exponent in exponents.items():
                if global_i not in high_position:
                    raise ValueError(f"unexpected variable p{global_i} outside the support-12 set")
                vector[high_position[global_i]] = exponent
                ordinary_degree += exponent
            coefficient = sp.Rational(values[0])
            terms[tuple(vector)] += coefficient
            term_count += 1
        if len(source_weight_degrees) != 1:
            raise AssertionError(f"source row {tag} is not positively weighted homogeneous")
        weighted_degree = next(iter(source_weight_degrees))
        restricted = sp.Poly.from_dict(
            {exponent: coefficient for exponent, coefficient in terms.items() if coefficient},
            *high_variables,
            domain=QQ,
        )
        if not restricted:
            continue
        if any(sum(monomial) not in (2, 3) for monomial, _ in restricted.terms()):
            raise AssertionError(f"unexpected top-weight boundary degree in {tag}")
        if any(sum(weights[i] * monomial[j] for j, i in enumerate(high_positions)) != weighted_degree
               for monomial, _ in restricted.terms()):
            raise AssertionError(f"restricted row {tag} lost its claimed weight")
        polynomials.append(restricted.as_expr())
        tags.append(tag)
        weighted_degrees[weighted_degree] += 1
        ordinary_degrees[restricted.total_degree()] += 1

    if not polynomials:
        raise SystemExit("top-weight restriction produced no equations")
    print(
        f"source_sha256={source_hash} selected={len(polynomials)} "
        f"terms_scanned={term_count} weighted={dict(sorted(weighted_degrees.items()))} "
        f"ordinary={dict(sorted(ordinary_degrees.items()))}",
        flush=True,
    )
    print("starting exact QQ grevlex F5B on seven-variable boundary ideal", flush=True)
    began = time.monotonic()
    basis = sp.groebner(polynomials, *high_variables, domain=QQ, order="grevlex", method="f5b")
    elapsed = time.monotonic() - began
    print(
        f"groebner_done basis_size={len(basis.polys)} "
        f"zero_dimensional={basis.is_zero_dimensional} seconds={elapsed:.2f}",
        flush=True,
    )

    basis_expressions = [polynomial.as_expr() for polynomial in basis.polys]
    p6, p7, p8, p9, p10, p11, p12 = high_variables
    l3 = 119 * p10 + 319 * p11 + 726 * p12 + p7 + 8 * p8 + 36 * p9
    cubic_l3 = any(sp.expand(expression - l3**3) == 0 for expression in basis_expressions)
    quadratic_candidates = [
        expression for expression in basis_expressions
        if sp.Poly(expression, *high_variables, domain=QQ).total_degree() == 2
    ]
    if len(quadratic_candidates) != 2 or not cubic_l3:
        raise SystemExit("unexpected boundary basis shape; cannot certify its set-theoretic form")
    q2_candidates = []
    other_quadratic = None
    for expression in quadratic_candidates:
        coefficient = sp.Poly(expression, *high_variables, domain=QQ).coeff_monomial(p6**2)
        if coefficient == 1:
            q2_candidates.append(expression)
        else:
            other_quadratic = expression
    if len(q2_candidates) != 1 or other_quadratic is None:
        raise SystemExit("could not identify the p6^2 quadric and the second quadratic")
    q2 = q2_candidates[0]
    other_factors = sp.factor_list(other_quadratic)[1]
    if not any(sp.expand(factor - l3) == 0 for factor, _power in other_factors):
        raise SystemExit("the remaining quadratic is not divisible by L3")
    l3_substitution = sp.solve(l3, p10)[0]
    q2_on_l3 = sp.expand(q2.subs(p10, l3_substitution))
    if sp.Poly(q2_on_l3, p6).coeff_monomial(p6**2) != 1:
        raise SystemExit("the boundary quadric changed its p6^2 coefficient on L3=0")
    set_theoretic_description = "L3=0 and q2=0"
    boundary_affine_dimension = 5

    nilpotence = {}
    if basis.is_zero_dimensional:
        for variable in high_variables:
            power_found = None
            for power in range(1, 65):
                if basis.reduce(variable**power)[1] == 0:
                    power_found = power
                    break
            nilpotence[str(variable)] = power_found
    result = {
        "source": SOURCE.name,
        "source_sha256": source_hash,
        "weight_profile_sha256": sha256(PROFILE),
        "scope": "exact support-12 subideal restricted to the weight-4 boundary; not a full 2+1+1 exclusion",
        "branch": "2+1+1",
        "coefficient_field": "QQ",
        "monomial_order": "grevlex",
        "algorithm": "SymPy F5B with SYMPY_GROUND_TYPES=flint",
        "boundary_variables": [str(variable) for variable in high_variables],
        "set_to_zero": [f"p{support[i]}" for i in low_positions],
        "selected_generators": len(polynomials),
        "weighted_degree_histogram": dict(sorted(weighted_degrees.items())),
        "ordinary_degree_histogram": dict(sorted(ordinary_degrees.items())),
        "term_count_after_restriction": term_count,
        "basis_size": len(basis.polys),
        "zero_dimensional": basis.is_zero_dimensional,
        "basis": [str(expression) for expression in basis_expressions],
        "basis_factorizations": [str(sp.factor(expression)) for expression in basis_expressions],
        "set_theoretic_zero_locus": set_theoretic_description,
        "L3": str(l3),
        "remaining_quadric_p6_squared_coefficient_on_L3": str(
            sp.Poly(q2_on_l3, *high_variables, domain=QQ).coeff_monomial(p6**2)
        ),
        "affine_dimension_of_boundary_zero_locus": boundary_affine_dimension,
        "variable_nilpotence_powers": nilpotence,
        "all_input_generators_reduce_to_zero": all(
            basis.reduce(polynomial)[1] == 0 for polynomial in polynomials
        ),
        "groebner_seconds": elapsed,
        "warning": "Even an origin-only boundary solution leaves all strata with nonzero lower-weight variables open and does not exclude the affine branch.",
    }
    RECEIPT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"receipt={RECEIPT}", flush=True)


if __name__ == "__main__":
    main()
