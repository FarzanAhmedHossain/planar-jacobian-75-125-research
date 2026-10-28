"""Exact QQ Groebner audit of the D24 eight-parameter necessary subsystem.

Input TSV is extracted from the D24 checkpoint by retaining every constraint
whose sparse monomial support is contained in generator indices
{6, 16, 35, 36, 37, 38, 39, 40}.  In the checkpoint's generator order these
are p5, p15, p34, ..., p39.  This subsystem is necessary, but is not the full
138-parameter D24 ideal.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import sympy as sp
from sympy import Poly, QQ, Rational, groebner, symbols


INPUT = Path(__file__).with_name("f2_75_125_d24_u1zero_support8_constraints.tsv")
EXPECTED_INPUT_SHA256 = (
    "5663f9e4a61a0d7518a6fd1b46af93d9ce6404fd42cb5efbbe278ff497548b7d"
)
CHECKPOINT_SHA256 = (
    "2369f36ce13d0d4336053ba0a7a546aa4407d04d77ee544aa4e313dbc88f2c79"
)
GENERATOR_INDICES = (6, 16, 35, 36, 37, 38, 39, 40)
VARIABLE_NAMES = ("p5", "p15", "p34", "p35", "p36", "p37", "p38", "p39")


def read_system() -> tuple[list[str], list[Poly]]:
    digest = hashlib.sha256(INPUT.read_bytes()).hexdigest()
    if digest != EXPECTED_INPUT_SHA256:
        raise ValueError(f"input digest mismatch: {digest}")

    variables = symbols(" ".join(VARIABLE_NAMES))
    rows: list[str] = []
    polynomials: list[Poly] = []
    allowed = set(GENERATOR_INDICES)
    for line in INPUT.read_text(encoding="utf-8").splitlines():
        tag, body = line.split("\t", 1)
        terms: dict[tuple[int, ...], Rational] = {}
        for item in body.split(";"):
            coefficient_text, monomial_text = item.split(":", 1)
            exponents = [0] * len(GENERATOR_INDICES)
            if monomial_text != "1":
                for factor in monomial_text.split(","):
                    index_text, exponent_text = factor.split("^", 1)
                    index = int(index_text)
                    if index not in allowed:
                        raise ValueError(f"unexpected generator index {index} in {tag}")
                    exponents[GENERATOR_INDICES.index(index)] = int(exponent_text)
            key = tuple(exponents)
            terms[key] = terms.get(key, Rational(0)) + Rational(coefficient_text)
        polynomial = Poly.from_dict(terms, *variables, domain=QQ)
        if not polynomial:
            raise ValueError(f"zero constraint unexpectedly retained: {tag}")
        rows.append(tag)
        polynomials.append(polynomial)
    if len(rows) != 143 or len(set(rows)) != 143:
        raise ValueError(f"expected 143 distinct rows, found {len(rows)}")
    return rows, polynomials


def main() -> None:
    rows, polynomials = read_system()
    variables = symbols(" ".join(VARIABLE_NAMES))
    print(f"checkpoint_sha256={CHECKPOINT_SHA256}")
    print(f"input_sha256={EXPECTED_INPUT_SHA256}")
    print(f"sympy={sp.__version__} constraints={len(rows)} variables={VARIABLE_NAMES}")
    print(f"degrees={sorted({p.total_degree() for p in polynomials})}")

    origin_satisfies = all(p.coeff_monomial(1) == 0 for p in polynomials)
    print(f"all_zero_parameter_point_satisfies_subsystem={origin_satisfies}")

    basis = groebner(
        [p.as_expr() for p in polynomials],
        *variables,
        order="grevlex",
        domain=QQ,
        method="f5b",
    )
    is_unit = any(p.as_expr() == 1 for p in basis.polys)
    print(
        f"basis_size={len(basis.polys)} zero_dimensional={basis.is_zero_dimensional} "
        f"unit_ideal={is_unit}"
    )
    for index, polynomial in enumerate(basis.polys):
        print(f"GB{index}={sp.factor(polynomial.as_expr())}")

    # The third reduced-basis element is a nonzero rational multiple of L^3.
    linear = (
        224725 * variables[2]
        + 1639725 * variables[3]
        + 6284625 * variables[4]
        + 15845845 * variables[5]
        + 24712395 * variables[6]
        - 163942748 * variables[7]
    )
    _, remainder = basis.reduce(linear**3)
    print(f"linear_cubic_reduces_to_zero={remainder == 0}")
    print("scope=necessary subsystem only; not a branch exclusion")


if __name__ == "__main__":
    main()
