#!/usr/bin/env python3
"""Low-memory exact incidence and univariate-gcd scan of a sparse F2 checkpoint.

The checkpoint's constraint array is decoded one row at a time.  This avoids
materializing a 500+ MB JSON document while another exact lift is using RAM.
The scan reports necessary-system structure only; a univariate gcd of 1 is an
exact inconsistency certificate for the chart represented by the checkpoint.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

from f2_75_125_211_d13_univariate_gcd import add, encode_poly, mul, scale, trim, xgcd


class JsonStream:
    """Incremental JSON value parser with a bounded buffer between array rows."""

    CHUNK = 1 << 20

    def __init__(self, path: Path):
        self._file = path.open("r", encoding="utf-8")
        self._decoder = json.JSONDecoder()
        self._buffer = ""
        self._position = 0
        self._base_position = 0
        self._ascii_so_far = True
        self._eof = False

    def close(self) -> None:
        self._file.close()

    def _fill(self, size: int = CHUNK) -> bool:
        if self._eof:
            return False
        chunk = self._file.read(size)
        if not chunk:
            self._eof = True
            return False
        self._ascii_so_far = self._ascii_so_far and chunk.isascii()
        self._buffer += chunk
        return True

    def _compact(self, *, force: bool = False) -> None:
        if self._position and (force or self._position >= self.CHUNK):
            self._base_position += self._position
            self._buffer = self._buffer[self._position:]
            self._position = 0

    def absolute_position(self) -> int:
        """Return the absolute character offset, valid as a byte offset for ASCII JSON."""
        self.skip_space()
        if not self._ascii_so_far:
            raise ValueError("raw checkpoint spans require ASCII JSON")
        return self._base_position + self._position

    def seek_absolute(self, position: int) -> None:
        """Seek to an absolute checkpoint offset previously returned above."""
        if not self._ascii_so_far:
            raise ValueError("raw checkpoint spans require ASCII JSON")
        self._file.seek(position)
        self._buffer = ""
        self._position = 0
        self._base_position = position
        self._eof = False

    def skip_polynomial_array(self, *, allowed_indices=None) -> tuple[int, int]:
        """Skip one serialized sparse-polynomial array and return its raw span.

        Each checkpoint term has the shape ``[sparse_powers, coefficients]``.
        Decode and discard those two small arrays individually, keeping neither
        a full polynomial dictionary nor a dense exponent tuple in memory.
        If ``allowed_indices`` is supplied, also certify that every nonzero
        variable in the raw span survives any ring compaction performed by the
        caller.  This matters when the raw JSON span is copied unchanged into a
        later checkpoint.
        """
        start = self.absolute_position()
        self.take("[")
        first = True
        while self.peek() != "]":
            if not first:
                self.take(",")
            first = False
            self.take("[")
            sparse = self.value()  # sparse exponent pairs
            if allowed_indices is not None:
                for index, exponent in sparse:
                    if exponent and index not in allowed_indices:
                        raise ValueError(
                            "raw polynomial span uses a generator omitted by ring compaction"
                        )
            self.take(",")
            self.value()  # exact coefficient encoding
            self.take("]")
        self.take("]")
        end = self.absolute_position()
        if not self._ascii_so_far:
            raise ValueError("raw checkpoint spans require ASCII JSON")
        return start, end

    def skip_space(self) -> None:
        while True:
            end = len(self._buffer)
            while self._position < end and self._buffer[self._position].isspace():
                self._position += 1
            if self._position < len(self._buffer) or not self._fill():
                return
            self._compact(force=True)

    def peek(self) -> str:
        self.skip_space()
        if self._position >= len(self._buffer):
            return ""
        return self._buffer[self._position]

    def take(self, expected: str) -> None:
        actual = self.peek()
        if actual != expected:
            raise ValueError(f"expected {expected!r}, found {actual!r}")
        self._position += 1
        self._compact()

    def value(self):
        self.skip_space()
        fill_size = self.CHUNK
        while True:
            try:
                result, end = self._decoder.raw_decode(self._buffer, self._position)
                self._position = end
                self._compact()
                return result
            except json.JSONDecodeError:
                # A sparse polynomial row can itself be very large. Growing
                # the read chunk geometrically avoids repeatedly reparsing
                # the same incomplete row in fixed 1 MiB increments.
                if not self._fill(fill_size):
                    raise
                fill_size = min(fill_size * 2, 64 << 20)

    def top_level_constraints(self) -> tuple[dict, "JsonStream"]:
        """Consume object metadata through the opening bracket of constraints."""
        self.take("{")
        metadata = {}
        first = True
        while True:
            mark = self.peek()
            if mark == "}":
                raise ValueError("checkpoint has no constraints array")
            if not first:
                self.take(",")
            first = False
            key = self.value()
            if not isinstance(key, str):
                raise ValueError("top-level JSON key is not a string")
            self.take(":")
            if key == "constraints":
                self.take("[")
                return metadata, self
            metadata[key] = self.value()

    def array_values(self):
        """Yield values from an array whose opening bracket was already read."""
        first = True
        while True:
            marker = self.peek()
            if marker == "]":
                self.take("]")
                return
            if not first:
                self.take(",")
            first = False
            yield self.value()


def update_gcd_state(state: dict, tag: str, equation: list[Fraction]) -> None:
    """Incrementally maintain a monic gcd and an exact Bezout witness."""
    equation = trim(equation)
    state["count"] += 1
    if state["gcd"] == [Fraction(1)]:
        return
    if state["gcd"] is None:
        leading = equation[-1]
        state["gcd"] = scale(equation, 1 / leading)
        state["bezout"] = [[1 / leading]]
        state["tags"] = [tag]
        state["equations"] = [equation]
        return
    common, left, right = xgcd(state["gcd"], equation)
    state["bezout"] = [mul(left, coefficient)
                        for coefficient in state["bezout"]]
    state["bezout"].append(right)
    state["tags"].append(tag)
    state["equations"].append(equation)
    state["gcd"] = common


def finish_gcd_state(state: dict) -> dict:
    gcd_poly = state["gcd"]
    if gcd_poly is None:
        return {
            "univariate_constraint_count": 0,
            "status": "no_univariate_constraints",
            "gcd_coefficients_low_to_high": None,
            "bezout_identity_replayed": None,
        }
    identity = [Fraction(0)]
    for coefficient, equation in zip(state["bezout"], state["equations"]):
        identity = add(identity, mul(coefficient, equation))
    if identity != gcd_poly:
        raise AssertionError("univariate Bezout certificate failed exact replay")
    return {
        "univariate_constraint_count": state["count"],
        "certificate_equation_count": len(state["equations"]),
        "certificate_tags": state["tags"],
        "certificate_equations_low_to_high": [
            encode_poly(equation) for equation in state["equations"]
        ],
        "bezout_coefficients_low_to_high": [
            encode_poly(poly) for poly in state["bezout"]
        ],
        "gcd_coefficients_low_to_high": encode_poly(gcd_poly),
        "bezout_identity_replayed": identity == gcd_poly,
        "status": "univariate_inconsistent" if len(gcd_poly) == 1
        else "nonconstant_common_factor",
    }


def scan(source: Path, nparameters: int, source_sha256: str,
         max_listed_support_size: int = 20) -> dict:
    reader = JsonStream(source)
    try:
        metadata, _ = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1":
            raise ValueError("unexpected checkpoint format")
        if metadata.get("branch") != "2+1+1" or metadata.get("quadratic"):
            raise ValueError("expected rational 2+1+1 checkpoint")
        if int(metadata.get("nparameters", -1)) != nparameters:
            raise ValueError("nparameters disagrees with checkpoint metadata")

        size_counts: Counter[int] = Counter()
        supports_by_size: dict[int, Counter[tuple[int, ...]]] = defaultdict(Counter)
        occurrences = [0] * nparameters
        parent = list(range(nparameters))
        ranks = [0] * nparameters
        auxiliary_generators: Counter[int] = Counter()
        gcd_states = {
            variable: {
                "count": 0, "gcd": None, "bezout": [], "tags": [],
                "equations": [],
            }
            for variable in range(nparameters)
        }
        constant_constraints = []
        constraint_count = term_count = sparse_entries = 0

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a: int, b: int) -> None:
            ra, rb = find(a), find(b)
            if ra == rb:
                return
            if ranks[ra] < ranks[rb]:
                ra, rb = rb, ra
            parent[rb] = ra
            if ranks[ra] == ranks[rb]:
                ranks[ra] += 1

        for row in reader.array_values():
            if not isinstance(row, list) or len(row) != 2:
                raise ValueError("malformed constraint row")
            tag, terms = row
            constraint_count += 1
            if constraint_count % 200 == 0:
                print(f"checkpoint scan {constraint_count} rows; "
                      f"terms={term_count}", flush=True)
            support: set[int] = set()
            univariate_possible = True
            univariate_variable = None
            univariate_coefficients: dict[int, Fraction] = {}
            constant_value = Fraction(0)
            has_auxiliary = False

            for sparse, values in terms:
                term_count += 1
                if len(values) != 1:
                    raise ValueError("expected rational checkpoint coefficient")
                coefficient = Fraction(values[0])
                sparse_entries += len(sparse)
                term_variables = set()
                powers_for_term: dict[int, int] = {}
                for generator, exponent in sparse:
                    if not exponent:
                        continue
                    if 1 <= generator <= nparameters:
                        variable = generator - 1
                        term_variables.add(variable)
                        powers_for_term[variable] = exponent
                    else:
                        auxiliary_generators[generator] += 1
                        has_auxiliary = True
                support.update(term_variables)
                if len(support) > 1 or len(term_variables) > 1 or has_auxiliary:
                    univariate_possible = False
                    univariate_coefficients.clear()
                elif term_variables:
                    variable = next(iter(term_variables))
                    if univariate_variable is None:
                        univariate_variable = variable
                        if constant_value:
                            univariate_coefficients[0] = constant_value
                    elif univariate_variable != variable:
                        univariate_possible = False
                        univariate_coefficients.clear()
                    if univariate_possible:
                        exponent = powers_for_term.get(variable, 0)
                        univariate_coefficients[exponent] = (
                            univariate_coefficients.get(exponent, Fraction(0))
                            + coefficient
                        )
                elif univariate_possible and not has_auxiliary:
                    constant_value += coefficient
                    if univariate_variable is not None:
                        univariate_coefficients[0] = (
                            univariate_coefficients.get(0, Fraction(0))
                            + coefficient
                        )

            ordered = tuple(sorted(support))
            size_counts[len(ordered)] += 1
            supports_by_size[len(ordered)][ordered] += 1
            for variable in ordered:
                occurrences[variable] += 1
            if ordered:
                for variable in ordered[1:]:
                    union(ordered[0], variable)

            if not has_auxiliary and len(ordered) == 0 and terms:
                if len(terms) == 1 and constant_value:
                    constant_constraints.append({"tag": tag, "value": str(constant_value)})
            elif (univariate_possible and len(ordered) == 1
                  and univariate_variable is not None):
                polynomial = [Fraction(0)] * (max(univariate_coefficients, default=0) + 1)
                for exponent, coefficient in univariate_coefficients.items():
                    polynomial[exponent] += coefficient
                polynomial = trim(polynomial)
                if any(polynomial):
                    update_gcd_state(gcd_states[univariate_variable], tag, polynomial)

        components: dict[int, list[int]] = defaultdict(list)
        for variable in range(nparameters):
            components[find(variable)].append(variable)
        component_sizes = sorted((len(items) for items in components.values()), reverse=True)
        repeated_supports = {
            str(size): [
                {"variables": list(support), "constraints": count}
                for support, count in counters.most_common(8)
            ]
            for size, counters in sorted(supports_by_size.items())
            if size <= max_listed_support_size
        }
        univariate = {
            f"p{variable}": finish_gcd_state(state)
            for variable, state in gcd_states.items() if state["count"]
        }
        empty = bool(constant_constraints) or any(
            row["status"] == "univariate_inconsistent" for row in univariate.values()
        )
        return {
            "source": str(source),
            "source_sha256": source_sha256,
            "branch": metadata["branch"],
            "nparameters": nparameters,
            "constraints": constraint_count,
            "terms": term_count,
            "sparse_entries": sparse_entries,
            "support_size_histogram": {str(k): v for k, v in sorted(size_counts.items())},
            "distinct_support_sets_by_size": {
                str(k): len(v) for k, v in sorted(supports_by_size.items())
            },
            "most_repeated_supports": repeated_supports,
            "variable_occurrence_counts": {
                f"p{i}": count for i, count in enumerate(occurrences)
            },
            "incidence_component_sizes": component_sizes,
            "auxiliary_generator_term_counts": dict(sorted(auxiliary_generators.items())),
            "constant_constraints": constant_constraints,
            "univariate_gcds": univariate,
            "necessary_constraint_system_empty": empty,
            "scope": (
                "Exact support and univariate gcd scan of the checkpoint's necessary "
                "constraint system. Empty means this chart's saved necessary equations "
                "are inconsistent; branch completeness still requires the root-profile "
                "and F2 normalization bridge audits."
            ),
        }
    finally:
        reader.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("nparameters", type=int)
    parser.add_argument("source_sha256")
    parser.add_argument("report", type=Path)
    parser.add_argument("--max-listed-support-size", type=int, default=20)
    args = parser.parse_args()
    result = scan(args.source, args.nparameters, args.source_sha256,
                  args.max_listed_support_size)
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "source_sha256", "nparameters", "constraints", "terms", "sparse_entries",
        "support_size_histogram", "distinct_support_sets_by_size",
        "incidence_component_sizes", "constant_constraints",
        "univariate_gcds", "necessary_constraint_system_empty",
    )}, indent=2))
    print(f"report={args.report}")


if __name__ == "__main__":
    main()
