"""Extract and factor low-support constraints from an exact F2 checkpoint.

This is a reconnaissance tool for finding a small exact subsystem.  It makes
no exclusion claim: any resulting equations remain necessary conditions, and
the zero parameter point may still lie on their variety.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

from f2_75_125_quadratic_affine_profile_stream import (
    iter_json_array,
    skip_ws,
    top_metadata_and_constraints,
)


def decode_coeff(values):
    root = sp.sqrt(5)
    return sp.Add(*(sp.Rational(value) * root ** (len(values) - 1 - i)
                    for i, value in enumerate(values)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--max-support", type=int, default=4)
    parser.add_argument("--groebner", action="store_true")
    parser.add_argument("--groebner-method", choices=("buchberger", "f5b"),
                        default="buchberger")
    parser.add_argument("--groebner-order", choices=("lex", "grlex", "grevlex"),
                        default="lex")
    parser.add_argument("--force-zero", nargs="*", default=[],
                        help="specialize these parameter generators to zero")
    parser.add_argument("--skip-factor", action="store_true",
                        help="skip factorization and report exact row shapes")
    args = parser.parse_args()

    profile = json.loads(args.receipt.read_text(encoding="utf-8"))
    text = args.checkpoint.read_text(encoding="utf-8")
    metadata, position = top_metadata_and_constraints(text)
    if profile.get("source_sha256"):
        import hashlib
        source_hash = hashlib.sha256(args.checkpoint.read_bytes()).hexdigest()
        if source_hash != profile["source_sha256"]:
            raise ValueError("support profile does not match the checkpoint hash")
    names = metadata["generators"]
    p_index = {
        names.index(f"p{i}"): i
        for i in range(int(metadata["nparameters"]))
    }
    forced = set(args.force_zero)
    unknown_forced = forced - {f"p{i}" for i in range(int(metadata["nparameters"]))}
    if unknown_forced:
        raise ValueError(f"forced parameters are not allocated: {sorted(unknown_forced)}")
    forced_indices = {names.index(name) for name in forced}
    if "constraint_profiles" in profile and not forced:
        wanted = {
            row["tag"] for row in profile["constraint_profiles"]
            if row["parameter_count"] <= args.max_support
        }
    else:
        # The rational-field profiler stores a bounded display sample rather
        # than every support row.  Recover the exact support set by streaming
        # the checkpoint itself, guarded by the profile's source hash.
        wanted = set()
        for tag, saved in iter_json_array(text, position):
            support = set()
            for sparse, _ in saved:
                if any(index in forced_indices and power
                       for index, power in sparse):
                    continue
                support.update(
                    p_index[index]
                    for index, power in sparse
                    if power and index in p_index
                )
            if len(support) <= args.max_support:
                wanted.add(tag)
    if not wanted:
        print("no constraints below support threshold")
        return

    extracted = {}
    for tag, saved in iter_json_array(text, position):
        if tag not in wanted:
            continue
        surviving_terms = [
            (sparse, values) for sparse, values in saved
            if not any(index in forced_indices and power
                       for index, power in sparse)
        ]
        support = sorted({p_index[index] for sparse, _ in surviving_terms
                          for index, power in sparse
                          if power and index in p_index})
        if not surviving_terms:
            continue
        symbols = sp.symbols(" ".join(f"p{i}" for i in support))
        if len(support) == 1:
            symbols = (symbols,)
        local_index = {global_index: i for i, global_index in enumerate(support)}
        terms = defaultdict(lambda: sp.S.Zero)
        for sparse, values in surviving_terms:
            monomial = [0] * len(support)
            for index, power in sparse:
                if index == 0 and power:
                    raise ValueError(f"unexpected t-power in constraint {tag}")
                if index in p_index and power:
                    monomial[local_index[p_index[index]]] = power
            terms[tuple(monomial)] += decode_coeff(values)
        domain = (QQ.algebraic_field(sp.sqrt(5)) if metadata.get("quadratic")
                  else QQ)
        polynomial = sp.Poly.from_dict(dict(terms), *symbols, domain=domain)
        extracted[tag] = (support, symbols, polynomial)

    groups = defaultdict(list)
    print(f"branch={profile['branch']} selected={len(extracted)} max_support={args.max_support} force_zero={sorted(forced)} field={'Q(sqrt(5))' if metadata.get('quadratic') else 'Q'}")
    for tag, (support, symbols, polynomial) in extracted.items():
        groups[tuple(support)].append((tag, symbols, polynomial))
        print(f"\n[{tag}] support={[f'p{i}' for i in support]} terms={len(polynomial.terms())} degree={polynomial.total_degree()}")
        if not args.skip_factor:
            # Factor in the already exact coefficient domain.  Passing an
            # algebraic-coefficient expression back through factor_list can
            # mishandle one-variable rows with scalar algebraic content.
            factored = polynomial.factor_list()
            print(f"content={factored[0]}")
            for factor, exponent in factored[1]:
                print(f"factor^{exponent}={factor}")

    if args.groebner:
        print("\nEXACT GROEBNER BASES BY IDENTICAL SUPPORT")
        for support, rows in groups.items():
            if len(rows) < 2:
                continue
            symbols = rows[0][1]
            expressions = [polynomial.as_expr() for _, _, polynomial in rows]
            if metadata.get("quadratic"):
                basis = sp.groebner(expressions, *symbols,
                                    extension=sp.sqrt(5),
                                    method=args.groebner_method,
                                    order=args.groebner_order)
            else:
                basis = sp.groebner(expressions, *symbols, domain=QQ,
                                    method=args.groebner_method,
                                    order=args.groebner_order)
            print(f"\nsupport={[f'p{i}' for i in support]} rows={len(rows)}")
            print(f"basis={list(basis)}")
            print(f"unit={basis.contains(sp.S.One)}")


if __name__ == "__main__":
    main()
