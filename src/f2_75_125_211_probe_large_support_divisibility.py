"""Probe whether the first post-substitution row above support 17 is redundant."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time

import sympy as sp
from flint import fmpq
from sympy.polys.domains import QQ


ROOT = Path("research")
SOURCE = ROOT / "f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
SUPPORT12_GB = ROOT / "f2_75_125_sparse_ring_lift_211_roots18_u1zero_support12_gb_4gen_weighted_flint.json"
OUTPUT = ROOT / "f2_75_125_211_large_support_divisibility_probe_20261025.json"
SOURCE_SHA256 = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if sha256(SOURCE) != SOURCE_SHA256:
        raise ValueError("d19 checkpoint hash mismatch")

    research = Path(__file__).resolve().parent
    sys.path.insert(0, str(research))
    from f2_75_125_flint_mpoly_backend import FlintRing
    from f2_75_125_211_streaming_checkpoint_scan import JsonStream

    gb_report = json.loads(SUPPORT12_GB.read_text(encoding="utf-8"))
    old_names = gb_report["selected_variables"]
    old_symbols = sp.symbols(" ".join(old_names))
    old_symbol_map = dict(zip(old_names, old_symbols))
    old_basis_poly = sp.Poly(
        sp.sympify(gb_report["basis"][3]["polynomial"],
                   locals=old_symbol_map),
        *old_symbols,
        domain=QQ,
    )

    reader = JsonStream(SOURCE)
    started = time.monotonic()
    result = None
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("expected rational 2+1+1 checkpoint")
        names = metadata["generators"]
        if int(metadata["nparameters"]) != 174 or "p174" not in names:
            raise ValueError("unexpected parameter map")
        ring = FlintRing(names)
        generators = {name: ring.gens[i] for i, name in enumerate(names)}
        p174_index = names.index("p174")
        divisor_terms = {}
        for exponents, coefficient in old_basis_poly.terms():
            full_exponents = [0] * len(names)
            for index, name in enumerate(old_names):
                target = "p174" if name == "r" else name
                full_exponents[names.index(target)] = int(exponents[index])
            divisor_terms[tuple(full_exponents)] = fmpq(
                int(coefficient.p), int(coefficient.q)
            )
        divisor = ring.from_prevalidated_dict(divisor_terms).data

        r = generators["p174"]
        substitutions = {
            "p13": ring.ground_new(fmpq(452, 5))*r,
            "p14": -14*r,
            "p15": r,
            "p34": (
                45*generators["p36"] + fmpq(1654, 5)*generators["p37"]
                + fmpq(7491, 5)*generators["p38"]
                - fmpq(1449178, 25)*generators["p39"]
                + fmpq(612164, 125)*generators["p5"]*r
            ),
            "p35": (
                -10*generators["p36"] - 55*generators["p37"]
                - fmpq(1102, 5)*generators["p38"]
                + fmpq(40222, 5)*generators["p39"]
                - fmpq(251692, 375)*generators["p5"]*r
            ),
        }

        rows_scanned = input_terms = 0
        for tag, sparse_terms in reader.array_values():
            rows_scanned += 1
            terms = {}
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError(f"non-rational coefficient in {tag}")
                exponents = [0] * len(names)
                for index, exponent in sparse:
                    exponents[int(index)] = int(exponent)
                key = tuple(exponents)
                terms[key] = terms.get(key, fmpq(0)) + fmpq(values[0])
            polynomial = ring.from_prevalidated_dict({
                exponents: coefficient
                for exponents, coefficient in terms.items()
                if coefficient
            }).compose(substitutions)
            if not len(polynomial):
                continue

            support = set()
            for exponents, _coefficient in polynomial.items():
                support.update(names[i] for i, exponent in enumerate(exponents)
                               if exponent)
            if len(support) <= 17:
                continue

            quotient, remainder = divmod(polynomial.data, divisor)
            support_names = sorted(
                "r" if name == "p174" else name for name in support
            )
            support_symbols = sp.symbols(" ".join(support_names))
            support_symbol_map = dict(zip(support_names, support_symbols))
            variable_map = {
                "r" if name == "p174" else name: index
                for index, name in enumerate(names)
            }
            row_coefficients = {}
            for exponents, coefficient in polynomial.items():
                monomial = tuple(
                    int(exponents[variable_map[name]]) for name in support_names
                )
                row_coefficients[monomial] = sp.Rational(str(coefficient))
            row_poly = sp.Poly.from_dict(
                row_coefficients, support_symbols, domain=QQ
            )

            divisor_basis = []
            for item in gb_report["basis"]:
                divisor_basis.append(sp.Poly(
                    sp.sympify(item["polynomial"], locals=support_symbol_map),
                    *support_symbols,
                    domain=QQ,
                ))
            old_weights = dict(zip(old_names, gb_report["monomial_weights"]))
            extended_weights = [old_weights.get(name, 1) for name in support_names]

            def extended_weighted_grevlex(monomial):
                return (
                    sum(weight * exponent
                        for weight, exponent in zip(extended_weights, monomial)),
                    tuple(reversed([-exponent for exponent in monomial])),
                )

            _membership_quotients, old_ideal_remainder = sp.reduced(
                row_poly, divisor_basis, *support_symbols,
                order=extended_weighted_grevlex,
            )
            canonical_terms = [
                (tuple(int(x) for x in polynomial.data.monomial(i)),
                 str(polynomial.data.coefficient(i)))
                for i in range(len(polynomial))
            ]
            canonical_terms.sort()
            canonical = json.dumps(canonical_terms, separators=(",", ":"))
            result = {
                "scope": (
                    "Single-row diagnostic only. Exact divisibility by one old "
                    "support-12 Groebner generator is sufficient for redundancy; "
                    "failure of this one-divisor test does not prove nonmembership "
                    "in the full support-12 ideal."
                ),
                "d19_checkpoint_sha256": SOURCE_SHA256,
                "support12_groebner_report_sha256": sha256(SUPPORT12_GB),
                "support12_divisor_basis_index_zero_based": 3,
                "row_tag": tag,
                "row_number": rows_scanned,
                "input_terms_scanned_through_row": input_terms,
                "reduced_support_size": len(support),
                "reduced_support": sorted("r" if n == "p174" else n for n in support),
                "reduced_term_count": len(polynomial),
                "divides_exactly": remainder.is_zero(),
                "quotient_term_count": len(quotient),
                "remainder_term_count": len(remainder),
                "normal_form_mod_full_support12_groebner_basis_is_zero": (
                    old_ideal_remainder.is_zero
                ),
                "normal_form_term_count": len(old_ideal_remainder.terms()),
                "normal_form": (
                    None if old_ideal_remainder.is_zero
                    else str(old_ideal_remainder.as_expr())
                ),
                "reduced_row_sparse_sha256": hashlib.sha256(
                    canonical.encode("utf-8")
                ).hexdigest(),
                "elapsed_seconds": time.monotonic()-started,
            }
            break
        reader.close()
    finally:
        try:
            reader.close()
        except Exception:
            pass

    if result is None:
        raise RuntimeError("no post-substitution row of support greater than 17 found")
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    print(f"report={OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
