"""Exact D10-D11 necessary-subsystem Groebner calculation over QQ.

The input is extracted from the D24 checkpoint by retaining constraints whose
monomial supports lie in generator indices 6, 16, 25..29, 35..40, and 45..49.
This is a necessary subsystem, not the complete D24 ideal.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from time import perf_counter

import sympy as sp
from sympy import Poly, QQ, Rational, groebner, symbols
from sympy.polys.polytools import reduced


INPUT = Path(__file__).with_name("f2_75_125_d24_u1zero_window18_constraints.tsv")
EXPECTED_INPUT_SHA256 = (
    "515046d7e36dd3c1d39972d7071b174a4b8356a6d8339a1e6126b9f36763be01"
)
CHECKPOINT_SHA256 = (
    "2369f36ce13d0d4336053ba0a7a546aa4407d04d77ee544aa4e313dbc88f2c79"
)
OLD_INDICES = (6, 16, 35, 36, 37, 38, 39, 40)
NEW_INDICES = (25, 26, 27, 28, 29, 45, 46, 47, 48, 49)
ALL_INDICES = OLD_INDICES + NEW_INDICES


def parse_polynomial(body: str, variables: tuple[sp.Symbol, ...]) -> Poly:
    positions = {index: position for position, index in enumerate(ALL_INDICES)}
    terms: dict[tuple[int, ...], Rational] = {}
    for item in body.split(";"):
        coefficient_text, monomial_text = item.split(":", 1)
        exponents = [0] * len(ALL_INDICES)
        if monomial_text != "1":
            for factor in monomial_text.split(","):
                index_text, exponent_text = factor.split("^", 1)
                index = int(index_text)
                if index not in positions:
                    raise ValueError(f"unexpected generator index {index}")
                exponents[positions[index]] = int(exponent_text)
        exponent_tuple = tuple(exponents)
        terms[exponent_tuple] = terms.get(exponent_tuple, Rational(0)) + Rational(
            coefficient_text
        )
    return Poly.from_dict(terms, *variables, domain=QQ)


def main() -> None:
    digest = hashlib.sha256(INPUT.read_bytes()).hexdigest()
    if digest != EXPECTED_INPUT_SHA256:
        raise ValueError(f"input digest mismatch: {digest}")

    variables = symbols(" ".join(f"z{index}" for index in ALL_INDICES))
    old_variables = variables[: len(OLD_INDICES)]
    old_set = set(OLD_INDICES)
    old_equations: list[Poly] = []
    new_equations: list[tuple[str, Poly]] = []
    for line in INPUT.read_text(encoding="utf-8").splitlines():
        tag, body = line.split("\t", 1)
        polynomial = parse_polynomial(body, variables)
        active = {
            ALL_INDICES[position]
            for monomial, _coefficient in polynomial.terms()
            for position, exponent in enumerate(monomial)
            if exponent
        }
        if active.issubset(old_set):
            old_terms = {
                monomial[: len(OLD_INDICES)]: coefficient
                for monomial, coefficient in polynomial.terms()
            }
            old_equations.append(
                Poly.from_dict(old_terms, *old_variables, domain=QQ)
            )
        else:
            new_equations.append((tag, polynomial))

    if (len(old_equations), len(new_equations)) != (143, 143):
        raise ValueError(
            f"expected 143 old and 143 added equations; found "
            f"{len(old_equations)} and {len(new_equations)}"
        )
    print(
        f"checkpoint_sha256={CHECKPOINT_SHA256}\n"
        f"input_sha256={digest}\n"
        f"sympy={sp.__version__} old={len(old_equations)} "
        f"added={len(new_equations)} variables={len(variables)}",
        flush=True,
    )

    start = perf_counter()
    old_basis = groebner(
        [polynomial.as_expr() for polynomial in old_equations],
        *old_variables,
        order="grevlex",
        domain=QQ,
        method="f5b",
    )
    print(
        f"old_basis_size={len(old_basis.polys)} "
        f"elapsed_seconds={perf_counter() - start:.3f}",
        flush=True,
    )

    start = perf_counter()
    remainders: list[Poly] = []
    for _tag, polynomial in new_equations:
        _quotients, remainder = reduced(
            polynomial.as_expr(),
            old_basis.polys,
            *variables,
            order="grevlex",
            domain=QQ,
        )
        if remainder != 0:
            remainders.append(Poly(remainder, *variables, domain=QQ))
    print(
        f"nonzero_added_remainders={len(remainders)} "
        f"remainder_terms={sum(len(p.terms()) for p in remainders)} "
        f"elapsed_seconds={perf_counter() - start:.3f}",
        flush=True,
    )

    start = perf_counter()
    basis = groebner(
        [polynomial.as_expr() for polynomial in old_basis.polys]
        + [polynomial.as_expr() for polynomial in remainders],
        *variables,
        order="grevlex",
        domain=QQ,
        method="f5b",
    )
    is_unit = any(polynomial.as_expr() == 1 for polynomial in basis.polys)
    result = [
        f"basis_size={len(basis.polys)} unit_ideal={is_unit} "
        f"zero_dimensional={basis.is_zero_dimensional} "
        f"elapsed_seconds={perf_counter() - start:.3f}",
        *(f"GB{index}={polynomial.as_expr()}" for index, polynomial in enumerate(basis.polys)),
    ]
    output = INPUT.with_suffix(".groebner_result.txt")
    output.write_text("\n".join(result) + "\n", encoding="utf-8")
    print("\n".join(result), flush=True)
    print(f"result_file={output}", flush=True)


if __name__ == "__main__":
    main()
