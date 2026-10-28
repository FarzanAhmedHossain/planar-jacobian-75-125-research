"""Memory-bounded affine audit for a sparse F2 checkpoint.

This parser decodes one saved constraint at a time.  It is intended for the
large distance-nine Q(sqrt(5)) states while another exact continuation is
still using memory.  A contradiction is returned only with an exact,
replayable source-row combination equal to 1.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import f2_75_125_quadratic_affine_profile as core


def skip_ws(text: str, pos: int) -> int:
    while pos < len(text) and text[pos].isspace():
        pos += 1
    return pos


def top_metadata_and_constraints(text: str) -> tuple[dict, int]:
    decoder = json.JSONDecoder()
    pos = skip_ws(text, 0)
    if text[pos] != "{":
        raise ValueError("checkpoint root is not a JSON object")
    pos += 1
    metadata: dict = {}
    keep = {"format", "branch", "quadratic", "generators", "nparameters", "p"}
    while True:
        pos = skip_ws(text, pos)
        if text[pos] == "}":
            raise ValueError("checkpoint has no constraints field")
        key, pos = decoder.raw_decode(text, pos)
        pos = skip_ws(text, pos)
        if text[pos] != ":":
            raise ValueError("malformed top-level JSON object")
        value_pos = skip_ws(text, pos + 1)
        if key == "constraints":
            if text[value_pos] != "[":
                raise ValueError("constraints field is not an array")
            return metadata, value_pos
        value, pos = decoder.raw_decode(text, value_pos)
        if key in keep:
            metadata[key] = value
        pos = skip_ws(text, pos)
        if text[pos] == ",":
            pos += 1
        elif text[pos] != "}":
            raise ValueError("malformed top-level JSON separator")
        else:
            raise ValueError("checkpoint has no constraints field")


def iter_json_array(text: str, pos: int):
    decoder = json.JSONDecoder()
    pos += 1
    while True:
        pos = skip_ws(text, pos)
        if text[pos] == "]":
            return
        item, pos = decoder.raw_decode(text, pos)
        yield item
        pos = skip_ws(text, pos)
        if text[pos] == ",":
            pos += 1
        elif text[pos] == "]":
            return
        else:
            raise ValueError("malformed constraints array")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def profile(source: Path) -> dict:
    text = source.read_text(encoding="utf-8")
    metadata, constraint_pos = top_metadata_and_constraints(text)
    if metadata.get("format") != "f2-sparse-ring-v1":
        raise ValueError("unexpected checkpoint format")
    if not metadata.get("quadratic"):
        raise ValueError("expected a quadratic coefficient field")

    names = metadata["generators"]
    nparameters = int(metadata["nparameters"])
    parameter_names = [f"p{i}" for i in range(nparameters)]
    parameter_indices = {
        names.index(name): i for i, name in enumerate(parameter_names)
    }
    affine_rows: list[list[core.Elt]] = []
    affine_tags: list[str] = []
    term_count = 0
    zero_constraints = 0
    nonlinear_constraints = 0
    nonlinear_reasons: dict[str, int] = {}
    nonlinear_samples: list[dict] = []
    total_degree_histogram: Counter[int] = Counter()
    maximum_degree_histogram: Counter[int] = Counter()
    univariate_constraints: list[dict] = []
    constraint_profiles: list[dict] = []
    constant_obstructions: list[dict] = []
    constraint_count = 0

    for tag, polynomial in iter_json_array(text, constraint_pos):
        constraint_count += 1
        if not polynomial:
            zero_constraints += 1
            continue
        row = [core.ZERO for _ in range(nparameters + 1)]
        affine = True
        reasons: set[str] = set()
        parameter_support: set[int] = set()
        degree_histogram: Counter[int] = Counter()
        max_degree = 0
        for sparse, values in polynomial:
            term_count += 1
            coefficient = core.from_json(values)
            degree = 0
            p_variable = None
            t_power = 0
            other_power = 0
            for generator_index, power in sparse:
                if generator_index == 0:
                    t_power += power
                elif generator_index in parameter_indices:
                    degree += power
                    if power:
                        parameter_support.add(parameter_indices[generator_index])
                    if p_variable is None:
                        p_variable = parameter_indices[generator_index]
                    elif power:
                        p_variable = -1
                else:
                    other_power += power
            degree_histogram[degree] += 1
            max_degree = max(max_degree, degree)
            total_degree_histogram[degree] += 1
            if t_power or other_power or degree > 1:
                affine = False
                if t_power:
                    reasons.add("formal_t_variable")
                if other_power:
                    reasons.add("auxiliary_generator")
                if degree > 1:
                    reasons.add("nonlinear_parameter_monomial")
            elif degree == 0:
                row[0] = core.add(row[0], coefficient)
            elif p_variable is None or p_variable < 0:
                affine = False
                reasons.add("multiple_parameter_monomial")
            else:
                row[p_variable + 1] = core.add(
                    row[p_variable + 1], coefficient
                )

        if affine:
            affine_rows.append(row)
            affine_tags.append(tag)
        else:
            nonlinear_constraints += 1
            for reason in reasons:
                nonlinear_reasons[reason] = nonlinear_reasons.get(reason, 0) + 1
            if len(nonlinear_samples) < 12:
                nonlinear_samples.append({"tag": tag, "reasons": sorted(reasons)})
        if len(polynomial) == 1 and not polynomial[0][0]:
            constant_obstructions.append(
                {"tag": tag, "value": core.show(core.from_json(polynomial[0][1]))}
            )
        maximum_degree_histogram[max_degree] += 1
        degree_data = {str(k): v for k, v in sorted(degree_histogram.items())}
        profile_row = {
            "tag": tag,
            "terms": len(polynomial),
            "parameter_count": len(parameter_support),
            "maximum_parameter_degree": max_degree,
            "parameter_degree_histogram": degree_data,
        }
        constraint_profiles.append(profile_row)
        if len(parameter_support) == 1:
            profile_row["parameter"] = parameter_names[next(iter(parameter_support))]
            univariate_constraints.append(profile_row)

    rank, pivots, inconsistency = core.rank_affine(affine_rows)
    certificate = None
    if inconsistency is not None:
        weights = inconsistency["source_row_combination"]
        combined = [core.ZERO for _ in range(nparameters + 1)]
        terms = []
        for row_index, coefficient in sorted(weights.items()):
            combined = [
                core.add(value, core.mul(coefficient, term))
                for value, term in zip(combined, affine_rows[row_index])
            ]
            terms.append(
                {"tag": affine_tags[row_index], "coefficient": core.show(coefficient)}
            )
        if combined[0] != core.ONE or any(
            not core.iszero(value) for value in combined[1:]
        ):
            raise AssertionError("reconstructed affine combination is not exactly 1")
        certificate = {
            "identity": "1 = sum(coefficient_i * constraint_i)",
            "terms": terms,
            "exact_replay": True,
        }

    return {
        "source_checkpoint": str(source),
        "source_sha256": file_sha256(source),
        "branch": metadata["branch"],
        "quadratic_field": "Q(sqrt(5))",
        "checkpoint_depth": None,
        "carried_parameter_count": nparameters,
        "constraint_count": constraint_count,
        "constraint_term_count": term_count,
        "zero_constraint_count": zero_constraints,
        "nonlinear_constraint_count": nonlinear_constraints,
        "nonlinear_reasons_by_constraint": nonlinear_reasons,
        "nonlinear_constraint_samples": nonlinear_samples,
        "nonzero_constant_obstructions": constant_obstructions,
        "term_total_parameter_degree_histogram": {
            str(k): v for k, v in sorted(total_degree_histogram.items())
        },
        "constraint_maximum_degree_histogram": {
            str(k): v for k, v in sorted(maximum_degree_histogram.items())
        },
        "univariate_constraint_count": len(univariate_constraints),
        "univariate_constraints": univariate_constraints,
        "constraint_profiles": constraint_profiles,
        "affine_constraint_count": len(affine_rows),
        "affine_coefficient_rank_over_Qsqrt5": rank,
        "affine_system_inconsistent": certificate is not None,
        "affine_unit_certificate": certificate,
        "affine_pivots": pivots,
        "scope": (
            "An inconsistent affine subsystem proves the checkpoint ideal is the "
            "unit ideal. Otherwise this exact profile is diagnostic only; nonlinear "
            "constraints remain to be decided."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(f"refusing to overwrite {args.receipt}")
    result = profile(args.checkpoint)
    args.receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    keys = (
        "source_sha256", "branch", "quadratic_field", "checkpoint_depth",
        "carried_parameter_count", "constraint_count", "constraint_term_count",
        "affine_constraint_count", "affine_coefficient_rank_over_Qsqrt5",
        "affine_system_inconsistent", "affine_unit_certificate",
        "term_total_parameter_degree_histogram",
        "constraint_maximum_degree_histogram", "univariate_constraint_count",
        "univariate_constraints",
    )
    print(json.dumps({key: result[key] for key in keys}, indent=2))
    print(f"receipt={args.receipt}")


if __name__ == "__main__":
    main()
