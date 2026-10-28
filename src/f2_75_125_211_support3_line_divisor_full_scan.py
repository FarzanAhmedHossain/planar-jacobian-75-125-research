"""Propagate the simple-root divisor reduction through all d19 equations.

The exact support-three Groebner basis gives the line
    (p13,p14,p15)=(452*r/5,-14*r,r).
The simple-root U5 divisibility equations, reduced modulo the exact support-10
ideal and restricted to this line, solve for p34 and p35. This script applies
that necessary specialization to every saved d19 equation, one row at a time,
using exact FLINT composition, then checks for constants and univariate-gcd
contradictions. A nonempty/undecided scan is not a branch exclusion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
import sys
import time


SOURCE_HASH = "0f293c437b1cc18f379e3a8d0bb2e98647616e062fa321e91cfdc45ebe67926d"
DEFAULT_SOURCE = Path(
    "research/f2_75_125_sparse_ring_lift_211_roots19_u1zero_streamed_state.json"
)
DEFAULT_REDUCTION = Path(
    "research/f2_75_125_211_simple_root_reduction_membership_20261025.json"
)
DEFAULT_OUTPUT = Path(
    "research/f2_75_125_211_support3_line_divisor_full_scan_20261025.json"
)


class ProgressJsonStream:
    """The existing streaming parser with row-level progress output."""

    def __init__(self, source, parser_class):
        self.reader = parser_class(source)
        self.rows = 0
        self.started = time.monotonic()

    def top_level_constraints(self):
        return self.reader.top_level_constraints()

    def array_values(self):
        for row in self.reader.array_values():
            self.rows += 1
            if self.rows % 100 == 0:
                print(
                    f"rows={self.rows} elapsed={time.monotonic()-self.started:.1f}s",
                    flush=True,
                )
            yield row

    def close(self):
        self.reader.close()


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def find(name, parents):
    parents.setdefault(name, name)
    if parents[name] != name:
        parents[name] = find(parents[name], parents)
    return parents[name]


def union(a, b, parents, ranks):
    ra, rb = find(a, parents), find(b, parents)
    if ra == rb:
        return
    if ranks[ra] < ranks[rb]:
        ra, rb = rb, ra
    parents[rb] = ra
    if ranks[ra] == ranks[rb]:
        ranks[ra] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--reduction", type=Path, default=DEFAULT_REDUCTION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--skip-source-hash", action="store_true")
    args = parser.parse_args()

    if not args.skip_source_hash:
        source_hash = sha256(args.source)
        if source_hash != SOURCE_HASH:
            raise ValueError(f"d19 source hash mismatch: {source_hash}")
    else:
        source_hash = SOURCE_HASH
    reduction = json.loads(args.reduction.read_text(encoding="utf-8"))
    if reduction["d19_source_sha256"] != SOURCE_HASH:
        raise ValueError("reduction receipt belongs to another d19 checkpoint")
    expected_solution = {
        "p34": "45*p36 + 1654*p37/5 + 7491*p38/5 - 1449178*p39/25 + 612164*p5*r/125",
        "p35": "-10*p36 - 55*p37 - 1102*p38/5 + 40222*p39/5 - 251692*p5*r/375",
    }
    actual_solution = reduction["restriction_to_known_support3_line"]["solution"]
    if actual_solution != expected_solution:
        raise ValueError("unexpected exact p34,p35 solution in reduction receipt")

    research = Path(__file__).resolve().parent
    sys.path.insert(0, str(research))
    from f2_75_125_flint_mpoly_backend import FlintRing
    from f2_75_125_211_streaming_checkpoint_scan import (
        JsonStream,
        finish_gcd_state,
        update_gcd_state,
    )
    from flint import fmpq

    parser_stream = ProgressJsonStream(args.source, JsonStream)
    try:
        metadata, _ = parser_stream.top_level_constraints()
        if (metadata.get("format") != "f2-sparse-ring-v1"
                or metadata.get("branch") != "2+1+1"
                or metadata.get("quadratic")):
            raise ValueError("expected rational F2 2+1+1 d19 checkpoint")
        names = metadata["generators"]
        nparameters = int(metadata["nparameters"])
        if nparameters != 174:
            raise ValueError(f"unexpected parameter count: {nparameters}")
        if names[0] != "t" or names[1] != "p0":
            raise ValueError("unexpected checkpoint variable ordering")

        ring = FlintRing(names)
        generators = {name: ring.gens[index]
                      for index, name in enumerate(names)}
        r = generators["p174"]
        r_index = names.index("p174")
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

        support_counts = Counter()
        supports_by_size = defaultdict(Counter)
        occurrences = Counter()
        parents = {}
        ranks = Counter()
        gcd_states = defaultdict(lambda: {
            "count": 0, "gcd": None, "bezout": [], "tags": [], "equations": []
        })
        constant_constraints = []
        input_rows = input_terms = output_terms = zero_rows = nonzero_rows = 0
        max_row_output_terms = 0
        auxiliary_counts = Counter()
        row_hash = hashlib.sha256()

        for tag, sparse_terms in parser_stream.array_values():
            input_rows += 1
            terms = {}
            for sparse, values in sparse_terms:
                input_terms += 1
                if len(values) != 1:
                    raise ValueError(f"non-rational coefficient in row {tag}")
                powers = [0]*len(names)
                for generator, exponent in sparse:
                    exponent = int(exponent)
                    if not exponent:
                        continue
                    if generator < 1 or generator > nparameters:
                        auxiliary_counts[names[generator]] += 1
                        raise ValueError(
                            f"non-parameter generator {names[generator]} in {tag}"
                        )
                    powers[generator] = exponent
                key = tuple(powers)
                coefficient = fmpq(values[0])
                terms[key] = terms.get(key, fmpq(0)) + coefficient
            terms = {powers: coefficient for powers, coefficient in terms.items()
                     if coefficient}
            polynomial = ring.from_prevalidated_dict(terms)
            reduced = polynomial.compose(substitutions)
            if not len(reduced):
                zero_rows += 1
                continue
            nonzero_rows += 1
            output_terms += len(reduced)
            max_row_output_terms = max(max_row_output_terms, len(reduced))

            support = set()
            exponents_by_term = []
            for powers, coefficient in reduced.items():
                monomial_variables = []
                monomial_degree = 0
                for index, exponent in enumerate(powers):
                    if exponent and (1 <= index <= nparameters or index == r_index):
                        name = "r" if index == r_index else names[index]
                        support.add(name)
                        monomial_variables.append((name, int(exponent)))
                        monomial_degree += int(exponent)
                    elif exponent:
                        auxiliary_counts[names[index]] += 1
                        raise ValueError(
                            f"unexpected non-parameter after substitution in {tag}"
                        )
                exponents_by_term.append((tuple(monomial_variables), coefficient,
                                          monomial_degree))
            ordered_support = tuple(sorted(support))
            support_counts[len(ordered_support)] += 1
            supports_by_size[len(ordered_support)][ordered_support] += 1
            for variable in ordered_support:
                occurrences[variable] += 1
                find(variable, parents)
            if ordered_support:
                for variable in ordered_support[1:]:
                    union(ordered_support[0], variable, parents, ranks)

            row_hash.update(tag.encode("utf-8"))
            for powers, coefficient, _ in exponents_by_term:
                row_hash.update(repr(powers).encode("utf-8"))
                row_hash.update(str(coefficient).encode("ascii"))

            if not ordered_support:
                value = sum((coefficient for _, coefficient, _ in exponents_by_term),
                            fmpq(0))
                if value:
                    constant_constraints.append({"tag": tag, "value": str(value)})
                continue
            if len(ordered_support) != 1:
                continue
            variable = ordered_support[0]
            univariate_coefficients = defaultdict(Fraction)
            for powers, coefficient, degree in exponents_by_term:
                univariate_coefficients[degree] += Fraction(str(coefficient))
            equation = [Fraction(0)] * (max(univariate_coefficients, default=0)+1)
            for degree, coefficient in univariate_coefficients.items():
                equation[degree] += coefficient
            if any(equation):
                update_gcd_state(gcd_states[variable], tag, equation)

        components = Counter()
        for variable in list(parents):
            components[find(variable, parents)] += 1
        component_sizes = sorted(components.values(), reverse=True)
        univariate = {name: finish_gcd_state(state)
                      for name, state in gcd_states.items() if state["count"]}
        inconsistent = {name: result for name, result in univariate.items()
                        if result["status"] == "univariate_inconsistent"}
        empty = bool(constant_constraints or inconsistent)
        result = {
            "scope": (
                "Exact propagation of the d19 support-three line and two simple-root "
                "U5 equations through every saved 2+1+1,U1=0 constraint. Empty would "
                "exclude this chart; this is not a 2+1+1 or (75,125) exclusion unless "
                "the necessary-system completeness bridge is separately established."
            ),
            "source": str(args.source),
            "source_sha256": source_hash,
            "reduction_receipt_sha256": sha256(args.reduction),
            "branch": metadata["branch"],
            "nparameters_before_reduction": nparameters,
            "substitution": {
                "p13": "452*r/5",
                "p14": "-14*r",
                "p15": "r",
                "p34": expected_solution["p34"],
                "p35": expected_solution["p35"],
                "r_storage": "reserved p174 slot; absent from the original 174 active parameters",
            },
            "checkpoint_rows_scanned": input_rows,
            "input_sparse_terms": input_terms,
            "nonzero_rows_after_specialization": nonzero_rows,
            "zero_rows_after_specialization": zero_rows,
            "output_sparse_terms": output_terms,
            "maximum_output_terms_in_one_row": max_row_output_terms,
            "specialized_constraint_stream_sha256": row_hash.hexdigest(),
            "support_size_histogram": {
                str(size): count for size, count in sorted(support_counts.items())
            },
            "most_repeated_supports": {
                str(size): [
                    {"variables": list(support), "constraints": count}
                    for support, count in counters.most_common(10)
                ]
                for size, counters in sorted(supports_by_size.items()) if size <= 20
            },
            "incidence_component_sizes": component_sizes,
            "constant_constraints": constant_constraints,
            "univariate_gcds": univariate,
            "inconsistent_univariate_gcds": inconsistent,
            "necessary_constraint_system_empty": empty,
            "whole_2+1+1_branch_excluded_percent": 0.0,
            "elapsed_seconds": time.monotonic()-parser_stream.started,
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
        print(json.dumps({key: result[key] for key in (
            "checkpoint_rows_scanned", "input_sparse_terms",
            "nonzero_rows_after_specialization", "zero_rows_after_specialization",
            "output_sparse_terms", "maximum_output_terms_in_one_row",
            "support_size_histogram", "incidence_component_sizes",
            "constant_constraints", "inconsistent_univariate_gcds",
            "necessary_constraint_system_empty", "elapsed_seconds",
        )}, indent=2), flush=True)
        print(f"report={args.output} specialized_sha256={row_hash.hexdigest()}", flush=True)
    finally:
        parser_stream.close()


if __name__ == "__main__":
    main()
