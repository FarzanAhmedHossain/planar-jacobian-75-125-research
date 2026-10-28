"""Apply an exact graph substitution to a sparse F2 checkpoint.

The depth-nine F3=0 and F4=0 charts share a two-generator local basis in
Q(sqrt(5)).  This script verifies the basis identities from saved rows, applies
the resulting polynomial graph to every stored polynomial, and writes a
replayable component sidecar.  It is a coordinate reduction, not an exclusion.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ

from f2_75_125_sparse_terminal_q_postcheck import load_sparse


TAG_Q0 = "terminal-alpha-P-distance2:u^3"
TAG_Q1 = "terminal-alpha-P-distance2:u^4"
EXPECTED_F3_TAGS = {
    TAG_Q0,
    TAG_Q1,
    "terminal-alpha-Q-distance2:u^7",
    "terminal-alpha-Q-distance2:u^8",
    "terminal-alpha-Q-distance2:u^9",
}
EXPECTED_F4_TAGS = EXPECTED_F3_TAGS | {
    "terminal-alpha-P-distance4:u^0",
    "terminal-alpha-Q-distance4:u^4",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def digest_state(state: dict) -> str:
    digest = hashlib.sha256()
    for sequence_name in ("p", "q", "w"):
        for polynomial in state[sequence_name]:
            for monomial, coefficient in sorted(polynomial.items()):
                digest.update(repr((monomial, coefficient)).encode())
    for tag, polynomial in state["constraints"]:
        digest.update(tag.encode())
        for monomial, coefficient in sorted(polynomial.items()):
            digest.update(repr((monomial, coefficient)).encode())
    return digest.hexdigest()


def encode(poly, field):
    if poly is None:
        return None
    terms = []
    for powers, coefficient in sorted(poly.items()):
        sparse = [[i, exponent] for i, exponent in enumerate(powers) if exponent]
        values = [coefficient] if field == QQ else coefficient.to_list()
        terms.append([sparse, [str(value) for value in values]])
    return terms


def iter_polynomials(state: dict):
    for name in ("h", "v", "S", "B", "kappa", "kappa2", "kappa3"):
        value = state.get(name)
        if value is not None:
            yield name, value
    for name, value in state["uk"].items():
        yield f"uk[{name}]", value
    for sequence_name in ("p", "q", "w"):
        for i, value in enumerate(state[sequence_name]):
            yield f"{sequence_name}[{i}]", value
    for i, (tag, value) in enumerate(state["constraints"]):
        yield f"constraints[{i}]::{tag}", value


def active_parameters(state: dict) -> list[str]:
    nparameters = int(state["data"]["nparameters"])
    active = set()
    for _name, poly in iter_polynomials(state):
        for powers, _coefficient in poly.items():
            for i in range(nparameters):
                if powers[i + 1]:  # generator 0 is t; p0 starts at index 1.
                    active.add(i)
    return [f"p{i}" for i in sorted(active)]


def coefficient_of(poly, powers):
    for monomial, coefficient in poly.items():
        if monomial == powers:
            return coefficient
    return None


def graph_reduce(source: Path, source_receipt: Path,
                 output: Path, output_receipt: Path) -> dict:
    for path in (output, output_receipt):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite {path}")

    prior = json.loads(source_receipt.read_text(encoding="utf-8"))
    component = prior.get("component")
    if component not in ("F3=0", "F4=0"):
        raise ValueError("this graph reduction is certified for F3=0 or F4=0")
    expected_tags = EXPECTED_F3_TAGS if component == "F3=0" else EXPECTED_F4_TAGS
    source_hash = sha256(source)

    state = load_sparse(source)
    data = state["data"]
    R, field, gens = state["R"], state["field"], state["gens"]
    if (data.get("format") != "f2-sparse-ring-v1" or state["data"].get("branch") != "2+2"
            or not data.get("quadratic")):
        raise ValueError("expected a quadratic-field F2 2+2 sparse checkpoint")
    if prior.get("digest") != digest_state(state):
        raise ValueError("component sidecar digest does not match the source checkpoint")
    if data.get("compaction", {}).get("source_checkpoint_sha256") != prior.get(
        "original_source_checkpoint_sha256"
    ):
        raise ValueError("component sidecar and standardized-source provenance disagree")

    p0, p1, p2, p3 = gens[1:5]
    allowed_indices = {1, 2, 3, 4}
    support_tags = set()
    rows_by_tag = {}
    for tag, poly in state["constraints"]:
        if tag in (TAG_Q0, TAG_Q1):
            rows_by_tag.setdefault(tag, []).append(poly)
        support = set()
        only_p0_to_p3 = True
        for powers, _coefficient in poly.items():
            for i, exponent in enumerate(powers):
                if exponent and i not in allowed_indices:
                    only_p0_to_p3 = False
                if exponent and i in allowed_indices:
                    support.add(i)
        if poly and only_p0_to_p3 and support == allowed_indices:
            support_tags.add(tag)
    if support_tags != expected_tags:
        raise ValueError(
            f"the four-parameter support rows differ from the audited {component} set"
        )
    if (set(rows_by_tag) != {TAG_Q0, TAG_Q1}
            or any(len(rows_by_tag[tag]) != 1 for tag in (TAG_Q0, TAG_Q1))):
        raise ValueError("the two normalized source rows are missing or duplicated")

    p0_squared = [0] * len(R.gens)
    p0_squared[1] = 2
    row0, row1 = rows_by_tag[TAG_Q0][0], rows_by_tag[TAG_Q1][0]
    lead0 = coefficient_of(row0, tuple(p0_squared))
    lead1 = coefficient_of(row1, tuple(p0_squared))
    if lead0 is None or lead1 is None or not lead0 or not lead1:
        raise ValueError("the source rows must have nonzero p0^2 coefficients")
    q0 = row0 * R.ground_new(field.quo(field.one, lead0))
    q1 = row1 * R.ground_new(field.quo(field.one, lead1))
    scalar = lambda value: R.ground_new(field.from_sympy(value))

    A_expr = sp.Rational(672, 5) + sp.Rational(288, 5)*sp.sqrt(5)
    B_expr = 1728 + sp.Rational(19008, 25)*sp.sqrt(5)
    C_expr = 8 + sp.Rational(8, 5)*sp.sqrt(5)
    D_expr = sp.Rational(459, 10) + sp.Rational(153, 10)*sp.sqrt(5)
    lambda_expr = sp.Rational(5779, 204480) + sp.Rational(2333, 204480)*sp.sqrt(5)
    c0_expr = 48 + 48*sp.sqrt(5)
    A, B, C, D = map(scalar, (A_expr, B_expr, C_expr, D_expr))
    lambda_unit, c0 = scalar(lambda_expr), scalar(c0_expr)

    g2 = lambda_unit*(q1-q0)
    g1 = q0-c0*g2
    expected_g2 = p1+C*p2+D*p3
    expected_g1 = p0**2-A*p2-B*p3
    if g2 != expected_g2 or g1 != expected_g1 or not lambda_unit:
        raise AssertionError("the exact two-row Groebner identities failed replay")

    inv_A = R.ground_new(field.quo(field.one, field.from_sympy(A_expr)))
    p2_solution = (p0**2-B*p3)*inv_A
    p1_solution = -C*p2_solution-D*p3
    replacements = {p1: p1_solution, p2: p2_solution}
    if g1.compose(replacements) or g2.compose(replacements):
        raise AssertionError("graph substitution does not solve g1=g2=0")
    selected = [(tag, poly) for tag, poly in state["constraints"]
                if tag in expected_tags]
    if len(selected) != len(expected_tags) or any(
        poly.compose(replacements) for _tag, poly in selected
    ):
        raise AssertionError("not every selected source row vanishes on the graph")

    for name in ("h", "v", "S", "B", "kappa", "kappa2", "kappa3"):
        if state.get(name) is not None:
            state[name] = state[name].compose(replacements)
    state["uk"] = {name: poly.compose(replacements)
                   for name, poly in state["uk"].items()}
    for sequence_name in ("p", "q", "w"):
        state[sequence_name] = [poly.compose(replacements)
                                for poly in state[sequence_name]]
    state["constraints"] = [(tag, poly.compose(replacements))
                            for tag, poly in state["constraints"]]

    if any(poly for tag, poly in state["constraints"] if tag in expected_tags):
        raise AssertionError("a selected local equation survived substitution")
    eliminated_indices = {2, 3}
    for name, poly in iter_polynomials(state):
        for powers, _coefficient in poly.items():
            if any(powers[index] for index in eliminated_indices):
                raise AssertionError(f"eliminated parameter remains in {name}")

    active_after = active_parameters(state)
    reduction = {
        "method": "exact polynomial graph substitution from normalized source rows",
        "component": component,
        "source_checkpoint": str(source),
        "source_checkpoint_sha256": source_hash,
        "source_row_tags": [TAG_Q0, TAG_Q1],
        "complete_small_support_row_tags": sorted(expected_tags),
        "basis_identities": {
            "g2": "lambda*(q1-q0)",
            "g1": "q0-(48+48*sqrt(5))*g2",
            "lambda": "(5779+2333*sqrt(5))/204480",
        },
        "groebner_generators": {
            "g1": "p0^2-(672/5+288*sqrt(5)/5)*p2-(1728+19008*sqrt(5)/25)*p3",
            "g2": "p1+(8+8*sqrt(5)/5)*p2+(459/10+153*sqrt(5)/10)*p3",
        },
        "substitution": {
            "p2": "(p0^2-(1728+19008*sqrt(5)/25)*p3)/(672/5+288*sqrt(5)/5)",
            "p1": "-(8+8*sqrt(5)/5)*p2-(459/10+153*sqrt(5)/10)*p3, with p2 substituted",
        },
        "selected_rows_vanish_exactly": len(selected),
        "eliminated_parameters_absent_from_all_saved_polynomials": ["p1", "p2"],
        "active_parameter_names_after": active_after,
        "active_parameter_count_after": len(active_after),
        "carried_parameter_prefix_preserved": int(data["nparameters"]),
        "scope": (
            f"Exact coordinate reduction of the {component} depth-nine necessary system; "
            "it preserves the saved polynomial zero set and is not an exclusion certificate."
        ),
    }
    compaction = dict(data.get("compaction") or {})
    compaction["exact_graph_reduction"] = reduction
    data["compaction"] = compaction

    for name in ("h", "v", "S", "B", "kappa", "kappa2", "kappa3"):
        if name in data or state.get(name) is not None:
            data[name] = encode(state.get(name), field)
    data["uk"] = {name: encode(poly, field) for name, poly in state["uk"].items()}
    for sequence_name in ("p", "q", "w"):
        data[sequence_name] = [encode(poly, field) for poly in state[sequence_name]]
    data["constraints"] = [[tag, encode(poly, field)]
                           for tag, poly in state["constraints"]]

    digest_after = digest_state(state)
    raw_output = (json.dumps(data, separators=(",", ":")) + "\n").encode("utf-8")
    output.write_bytes(raw_output)
    output_hash = hashlib.sha256(raw_output).hexdigest()
    del state, data
    gc.collect()

    replay = load_sparse(output)
    if digest_state(replay) != digest_after:
        raise AssertionError("saved graph-reduced checkpoint failed digest replay")
    if any(poly for tag, poly in replay["constraints"] if tag in expected_tags):
        raise AssertionError("saved graph-reduced checkpoint lost a selected zero row")
    for name, poly in iter_polynomials(replay):
        for powers, _coefficient in poly.items():
            if any(powers[index] for index in eliminated_indices):
                raise AssertionError(f"eliminated parameter reappeared in saved {name}")

    derived = dict(prior)
    derived["derived_checkpoint"] = str(output)
    derived["derived_checkpoint_sha256"] = output_hash
    derived["digest"] = digest_after
    derived["constraint_count_after"] = len(replay["constraints"])
    derived["constraint_term_count_after"] = sum(
        len(poly) for _tag, poly in replay["constraints"]
    )
    derived["exact_graph_reduction"] = reduction
    derived["reload_digest_verified"] = True
    derived["scope"] = (
        f"Equivalent graph reduction of the retained {component} depth-nine necessary system; "
        "no branch exclusion is claimed."
    )
    output_receipt.write_text(json.dumps(derived, indent=2) + "\n", encoding="utf-8")
    return {
        "source_sha256": source_hash,
        "output_sha256": output_hash,
        "output_bytes": output.stat().st_size,
        "branch": replay["data"]["branch"],
        "component": component,
        "depth": len(replay["p"])-1,
        "carried_parameter_prefix": int(replay["data"]["nparameters"]),
        "active_parameter_count": len(active_after),
        "constraint_count": len(replay["constraints"]),
        "constraint_terms": sum(len(poly) for _tag, poly in replay["constraints"]),
        "small_support_rows_vanishing": len(selected),
        "eliminated_parameters_absent": True,
        "reload_digest_matches": digest_state(replay) == digest_after,
        "scope": f"exact reduction only; {component} exclusion coverage is unchanged",
        "checkpoint": str(output),
        "sidecar": str(output_receipt),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("source_receipt", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("output_receipt", type=Path)
    args = parser.parse_args()
    print(json.dumps(graph_reduce(args.source, args.source_receipt,
                                  args.output, args.output_receipt), indent=2))


if __name__ == "__main__":
    main()
