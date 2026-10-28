"""Exact sparse-parameter replay of the coupled F2 edge bands.

All coefficient arithmetic stays in K[t,p_0,...,p_{r-1}], where K is QQ or
QQ(sqrt(5)); no expanded SymPy expression domain is used for parameter
polynomials.  At each affine block a constant pivot submatrix is inverted
over K, and EVERY original coefficient equation is checked after solving.

The fifth block retains the full source support P_10=t^2 U_5, deg U_5<=12.
This includes the source-translation parameter omitted by the centered
fifthbands script.  Necessary support relaxation only; no pair exclusion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from math import comb, factorial
from pathlib import Path
import time

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring as sympy_ring

from f2_75_125_power_relation_pq_nextbands import weight49_parameters
from f2_75_125_211_streaming_checkpoint_scan import JsonStream


POLY_BACKEND = "sympy"
DEFAULT_PARAMETER_RESERVE = 192


def ring(names, field):
    """Select the exact sparse polynomial backend for a QQ recurrence."""
    if POLY_BACKEND == "flint-mpoly":
        if field != QQ:
            raise TypeError("flint-mpoly backend currently supports only QQ")
        from f2_75_125_flint_mpoly_backend import FlintRing
        return FlintRing(names, field)
    return sympy_ring(names, field)


class LazyConstraintPolynomial:
    """A checkpoint constraint kept as an exact source-file span.

    Old compatibility equations are not part of the next affine solve.  Keeping
    their serialized polynomial spans avoids holding the full constraint set
    in FLINT memory; replacements are recorded and applied one row at a time
    when a checkpoint is written or a digest is explicitly requested.
    """

    __slots__ = ("path", "start", "end", "transformations")

    def __init__(self, path, start, end, transformations=()):
        self.path = Path(path)
        self.start = int(start)
        self.end = int(end)
        self.transformations = tuple(transformations)

    def compose_lazily(self, replacements):
        return LazyConstraintPolynomial(
            self.path, self.start, self.end,
            (*self.transformations, dict(replacements)),
        )


def _decode_sparse_checkpoint_polynomial(reader, target_ring, field,
                                         old_names, new_indices, quadratic):
    """Decode one sparse checkpoint polynomial into the requested ring.

    Large checkpoint coefficients can contain millions of terms.  Building
    their full dense-exponent dictionary at once duplicates a substantial
    amount of memory before FLINT takes ownership of the polynomial.  Flush
    bounded dictionaries into exact ring polynomials and accumulate those
    chunks instead.  This preserves every term while limiting the Python-side
    exponent-tuple peak.
    """
    target_indices = getattr(target_ring, "_symbol_index", None)
    result = {}
    assembled = None
    chunk_limit = 65536
    term_count = 0
    next_progress = 1_000_000
    builder = (target_ring.from_prevalidated_dict
               if hasattr(target_ring, "from_prevalidated_dict")
               else target_ring.from_dict)

    def flush_chunk():
        nonlocal result, assembled
        if not result:
            return
        chunk = result
        result = {}
        partial = builder(chunk)
        del chunk
        if assembled is None:
            assembled = partial
        else:
            assembled = assembled + partial

    reader.take("[")
    first_term = True
    while reader.peek() != "]":
        if not first_term:
            reader.take(",")
        first_term = False
        reader.take("[")
        sparse = reader.value()
        reader.take(",")
        values = reader.value()
        reader.take("]")
        powers = [0]*len(target_ring.gens)
        for old_index, exponent in sparse:
            name = old_names[old_index]
            if target_indices is None:
                index = new_indices[name]
            else:
                index = target_indices.get(name)
            if index is None:
                if exponent:
                    raise ValueError("saved polynomial uses an unavailable generator")
                continue
            powers[index] = exponent
        if quadratic:
            ground_coefficients = [field.dom.dtype(value) for value in values]
            coefficient = field.dtype(ground_coefficients, field.mod, field.dom)
        else:
            coefficient = field.dtype(values[0])
        result[tuple(powers)] = coefficient
        term_count += 1
        if len(result) >= chunk_limit:
            flush_chunk()
        if term_count >= next_progress:
            print(f"checkpoint_polynomial_decoded_terms={term_count}",
                  flush=True)
            next_progress += 1_000_000
    reader.take("]")
    if result:
        flush_chunk()
    return target_ring.zero if assembled is None else assembled


class SparseLift:
    def __init__(self, tag, h_expr, v_expr, field,
                 parameter_reserve=DEFAULT_PARAMETER_RESERVE):
        self.tag, self.field = tag, field
        if parameter_reserve < 1:
            raise ValueError("parameter reserve must be positive")
        names = ["t", *(f"p{i}" for i in range(parameter_reserve)),
                 *(f"x{i}" for i in range(73))]
        self.ring, *generators = ring(names, field)
        self.t = generators[0]
        self.parameters = generators[1:1+parameter_reserve]
        self.unknowns = generators[1+parameter_reserve:]
        self.new_start = 1 + len(self.parameters)
        self.affine_ring = (self.ring.prefix(self.new_start)
                            if getattr(self.ring, "backend_name", None) == "flint-mpoly"
                            else None)
        self.nparameters = 0
        self.compaction = None
        self.receipt = []
        self.constraints = []
        self.started = time.monotonic()
        self.h = self.ring.from_expr(h_expr)
        self.v = self.ring.from_expr(v_expr)
        self.u = self.t*self.h
        old, s_expr, f_expr, _ = weight49_parameters(h_expr, v_expr, field, tag)
        replacement = dict(zip(old, sp.symbols(f"p0:{len(old)}")))
        self.S = self.ring.from_expr(s_expr.xreplace(replacement))
        self.w = [self.u**7*self.v, self.ring.from_expr(f_expr.xreplace(replacement))]
        self.nparameters = len(old)
        self._refresh_affine_ring()
        self.uk = {1: self.h**2*self.S}
        self.p = [self.u**3, self.t**3*self.uk[1]]
        self.q = [self.u**5, self.const(sp.Rational(5, 3))*self.t**5*self.h**2*self.uk[1]]
        self.B = None
        self.kappa = self.ring.zero
        self.kappa2 = self.ring.zero
        self.kappa3 = self.ring.zero

    @classmethod
    def from_checkpoint(cls, path, *, verify_saved_terminal_rows=True,
                        parameter_reserve=DEFAULT_PARAMETER_RESERVE):
        """Resume exact state; optionally skip replay of saved terminal-P rows.

        The default rechecks every saved row.  A separate calculation may set
        ``verify_saved_terminal_rows=False`` only when the checkpoint has
        already been independently audited; this skips that duplicate replay,
        not any coefficient constraints stored in the checkpoint.
        """
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        field = QQ.algebraic_field(sp.sqrt(5)) if data["quadratic"] else QQ
        old_names = data["generators"]
        obj = cls.__new__(cls)
        obj.tag, obj.field = data["branch"], field
        if parameter_reserve < 1:
            raise ValueError("parameter reserve must be positive")
        stored_indices = {name: i for i, name in enumerate(old_names)}
        if "t" not in stored_indices:
            raise ValueError("checkpoint has no t generator")
        parameter_names = sorted(
            (name for name in old_names
             if name.startswith("p") and name[1:].isdigit()),
            key=lambda name: int(name[1:]),
        )
        expected_parameters = [f"p{i}" for i in range(len(parameter_names))]
        if parameter_names != expected_parameters:
            raise ValueError("checkpoint parameter generators must be contiguous from p0")
        reserve = max(parameter_reserve, len(parameter_names))
        unknown_names = sorted(
            (name for name in old_names
             if name.startswith("x") and name[1:].isdigit()),
            key=lambda name: int(name[1:]),
        )
        if unknown_names != [f"x{i}" for i in range(73)]:
            raise ValueError("checkpoint must retain unknown generators x0 through x72")
        # Checkpoints are commonly compacted to the parameters that have been
        # used so far. Restore a bounded reserve before the x-unknowns so the
        # next affine block can allocate fresh parameters without changing any
        # saved sparse exponent vector.
        insertion = old_names.index("x0")
        extra_parameters = [f"p{i}" for i in range(len(parameter_names), reserve)]
        ring_names = old_names[:insertion] + extra_parameters + old_names[insertion:]
        ring_indices = {name: i for i, name in enumerate(ring_names)}
        obj.ring, *generators = ring(ring_names, field)
        if ring_names[0] != "t" or ring_names[1:1+reserve] != [f"p{i}" for i in range(reserve)]:
            raise ValueError("checkpoint generators are not ordered as t, p-parameters, x-unknowns")
        new_indices = ring_indices
        obj.t = generators[new_indices["t"]]
        parameter_names = [f"p{i}" for i in range(reserve)]
        obj.parameters = [generators[new_indices[name]] for name in parameter_names]
        obj.unknowns = [generators[new_indices[name]] for name in unknown_names]
        obj.new_start = 1 + len(obj.parameters)
        obj.affine_ring = (obj.ring.prefix(obj.new_start)
                           if getattr(obj.ring, "backend_name", None) == "flint-mpoly"
                           else None)
        obj.compaction = data.get("compaction")

        def decode(terms):
            result = {}
            for sparse, values in terms:
                powers = [0]*len(generators)
                for old_index, exponent in sparse:
                    powers[new_indices[old_names[old_index]]] = exponent
                if data["quadratic"]:
                    # The JSON stores ANP.to_list() in the field's polynomial
                    # basis.  Rebuilding a SymPy expression with sqrt(5) and
                    # then converting it back is equivalent but much slower
                    # for large checkpoints; construct the exact ANP directly.
                    ground_coefficients = [
                        field.dom.dtype(value) for value in values
                    ]
                    coefficient = field.dtype(
                        ground_coefficients, field.mod, field.dom
                    )
                else:
                    coefficient = field.dtype(values[0])
                result[tuple(powers)] = coefficient
            return obj.ring.from_dict(result)

        obj.h, obj.v, obj.S = (decode(data[key]) for key in ("h", "v", "S"))
        obj.u = obj.t * obj.h
        obj.B = None if data["B"] is None else decode(data["B"])
        obj.kappa = decode(data["kappa"])
        obj.kappa2 = decode(data["kappa2"]) if "kappa2" in data else obj.ring.zero
        obj.kappa3 = decode(data["kappa3"]) if "kappa3" in data else obj.ring.zero
        obj.uk = {int(k): decode(v) for k, v in data["uk"].items()}
        obj.p = [decode(v) for v in data["p"]]
        obj.q = [decode(v) for v in data["q"]]
        obj.w = [decode(v) for v in data["w"]]
        obj.nparameters = data["nparameters"]
        obj._refresh_affine_ring()
        obj.receipt = data["receipt"]
        obj.constraints = [(tag, decode(poly))
                           for tag, poly in data.get("constraints", [])]
        obj.started = time.monotonic()
        obj.local_C = obj.local_coefficient(obj.local_slice(0, 6), 6)
        zero_powers = (0,)*len(obj.ring.gens)
        c0 = obj.local_C[zero_powers]
        obj.center_a = -obj.local_coefficient(obj.local_slice(1, 5), 5)/(6*c0)
        obj.center_b = obj.ring.zero
        second = obj.centered_slice(2, 5)
        obj.center_b = -obj.local_coefficient(second, 5)/(6*c0)
        if verify_saved_terminal_rows:
            # Recheck every saved terminal-P slice on resume.  The original
            # 32 forbidden coefficients are distances 0..5; later
            # checkpoints may carry additional terminal slices, so the guard
            # must not silently stop at those first six.
            for distance in range(len(obj.p)):
                bound = (71-5*distance)//12
                if bound >= 0 and obj.centered_slice(distance, bound):
                    raise AssertionError(f"saved P terminal rows fail at distance {distance}")
        print(f"resumed={path} branch={obj.tag} depth={len(obj.p)-1} parameters={obj.nparameters}", flush=True)
        return obj

    @classmethod
    def from_checkpoint_streaming(cls, path, *, verify_saved_terminal_rows=True,
                                  parameter_reserve=DEFAULT_PARAMETER_RESERVE):
        """Resume an exact checkpoint while decoding one sparse row at a time.

        The ordinary JSON loader keeps both the complete 500+ MB JSON object
        and its reconstructed FLINT polynomials alive at once.  The checkpoint
        writer places ``constraints`` before the remaining polynomial state,
        so this reader can build the ring from the small metadata prefix and
        decode each retained equation directly into that ring.  It preserves
        every saved constraint and performs the same terminal-row replay as
        :meth:`from_checkpoint`.
        """
        reader = JsonStream(Path(path))
        try:
            metadata, _ = reader.top_level_constraints()
            if metadata.get("format") != "f2-sparse-ring-v1":
                raise ValueError("unexpected checkpoint format")
            if "branch" not in metadata or "quadratic" not in metadata:
                raise ValueError("checkpoint metadata is incomplete")
            field = QQ.algebraic_field(sp.sqrt(5)) if metadata["quadratic"] else QQ
            old_names = metadata["generators"]
            obj = cls.__new__(cls)
            obj.tag, obj.field = metadata["branch"], field
            if parameter_reserve < 1:
                raise ValueError("parameter reserve must be positive")
            stored_indices = {name: i for i, name in enumerate(old_names)}
            if "t" not in stored_indices:
                raise ValueError("checkpoint has no t generator")
            parameter_names = sorted(
                (name for name in old_names
                 if name.startswith("p") and name[1:].isdigit()),
                key=lambda name: int(name[1:]),
            )
            expected_parameters = [f"p{i}" for i in range(len(parameter_names))]
            if parameter_names != expected_parameters:
                raise ValueError("checkpoint parameter generators must be contiguous from p0")
            active_parameters = metadata["nparameters"]
            if not 0 <= active_parameters <= len(parameter_names):
                raise ValueError("checkpoint active parameter count is invalid")
            reserve = max(parameter_reserve, active_parameters)
            unknown_names = sorted(
                (name for name in old_names
                 if name.startswith("x") and name[1:].isdigit()),
                key=lambda name: int(name[1:]),
            )
            if unknown_names != [f"x{i}" for i in range(73)]:
                raise ValueError("checkpoint must retain unknown generators x0 through x72")
            if old_names != ["t", *parameter_names, *unknown_names]:
                raise ValueError("checkpoint generators are not ordered as t, p-parameters, x-unknowns")
            # A checkpoint may reserve many parameter slots that have never
            # been allocated.  Keep only the requested reserve (and the
            # carried parameters) so the resumed FLINT ring uses shorter
            # exponent vectors.  Sparse terms are remapped by generator name;
            # state polynomials using a dropped slot fail during decoding, and
            # the raw lazy constraints are checked below before being copied.
            ring_names = ["t", *(f"p{i}" for i in range(reserve)),
                          *unknown_names]
            ring_indices = {name: i for i, name in enumerate(ring_names)}
            obj.ring, *generators = ring(ring_names, field)
            if ring_names[0] != "t" or ring_names[1:1+reserve] != [f"p{i}" for i in range(reserve)]:
                raise ValueError("checkpoint generators are not ordered as t, p-parameters, x-unknowns")
            new_indices = ring_indices
            obj.t = generators[new_indices["t"]]
            obj.parameters = [generators[new_indices[f"p{i}"]] for i in range(reserve)]
            obj.unknowns = [generators[new_indices[name]] for name in unknown_names]
            obj.new_start = 1 + len(obj.parameters)
            obj._checkpoint_path = Path(path)
            obj._checkpoint_old_names = tuple(old_names)
            obj._checkpoint_new_indices = new_indices
            obj._checkpoint_quadratic = bool(metadata["quadratic"])
            obj.affine_ring = (obj.ring.prefix(obj.new_start)
                               if getattr(obj.ring, "backend_name", None) == "flint-mpoly"
                               else None)
            obj.compaction = metadata.get("compaction")
            obj.nparameters = metadata["nparameters"]
            obj._refresh_affine_ring()
            obj.receipt = metadata["receipt"]
            obj.constraints = []
            # Saved compatibility constraints involve only t and allocated
            # carried parameters, never the 73 temporary affine unknowns or
            # unused reserve slots.  Decode them in this shorter prefix ring
            # to reduce both the Python exponent-tuple peak and FLINT term
            # storage during large checkpoint resumes.
            obj.constraint_ring = (
                obj.ring.prefix(1+obj.nparameters)
                if getattr(obj.ring, "backend_name", None) == "flint-mpoly"
                else obj.ring
            )

            def decode_polynomial(target_ring=None):
                target_ring = target_ring or obj.ring
                return _decode_sparse_checkpoint_polynomial(
                    reader, target_ring, field, old_names, new_indices,
                    metadata["quadratic"],
                )

            first_constraint = True
            while reader.peek() != "]":
                if not first_constraint:
                    reader.take(",")
                first_constraint = False
                reader.take("[")
                tag = reader.value()
                reader.take(",")
                polynomial_start, polynomial_end = reader.skip_polynomial_array(
                    allowed_indices=frozenset(range(1+active_parameters))
                )
                reader.take("]")
                obj.constraints.append((tag, LazyConstraintPolynomial(
                    path, polynomial_start, polynomial_end,
                )))
            reader.take("]")

            def decode_array():
                reader.take("[")
                result = []
                first = True
                while reader.peek() != "]":
                    if not first:
                        reader.take(",")
                    first = False
                    polynomial = decode_polynomial()
                    result.append(polynomial)
                    print(f"checkpoint_array_polynomials={len(result)} "
                          f"terms={len(polynomial)}", flush=True)
                reader.take("]")
                return result

            def decode_map():
                reader.take("{")
                result = {}
                first = True
                while reader.peek() != "}":
                    if not first:
                        reader.take(",")
                    first = False
                    key = reader.value()
                    reader.take(":")
                    result[int(key)] = decode_polynomial()
                reader.take("}")
                return result

            obj.h = obj.v = obj.S = obj.ring.zero
            obj.B = None
            obj.kappa = obj.kappa2 = obj.kappa3 = obj.ring.zero
            obj.uk = {}
            obj.p = obj.q = obj.w = []
            required = {"h", "v", "S", "kappa", "uk", "p", "q", "w"}
            seen = set()
            reader.take(",")
            first = True
            while reader.peek() != "}":
                if not first:
                    reader.take(",")
                first = False
                key = reader.value()
                reader.take(":")
                if key in {"h", "v", "S", "kappa", "kappa2", "kappa3"}:
                    setattr(obj, key, decode_polynomial())
                elif key == "B":
                    if reader.peek() == "n":
                        obj.B = reader.value()
                    else:
                        obj.B = decode_polynomial()
                elif key == "uk":
                    obj.uk = decode_map()
                elif key in {"p", "q", "w"}:
                    setattr(obj, key, decode_array())
                elif key == "compaction":
                    obj.compaction = reader.value()
                else:
                    raise ValueError(f"unexpected checkpoint field after constraints: {key}")
                seen.add(key)
            reader.take("}")
            if reader.peek():
                raise ValueError("unexpected trailing data after checkpoint object")
            if not required.issubset(seen):
                raise ValueError(f"checkpoint is missing fields: {sorted(required-seen)}")
            obj.u = obj.t * obj.h
            obj.started = time.monotonic()
            obj.local_C = obj.local_coefficient(obj.local_slice(0, 6), 6)
            zero_powers = (0,)*len(obj.ring.gens)
            c0 = obj.local_C[zero_powers]
            obj.center_a = -obj.local_coefficient(obj.local_slice(1, 5), 5)/(6*c0)
            obj.center_b = obj.ring.zero
            second = obj.centered_slice(2, 5)
            obj.center_b = -obj.local_coefficient(second, 5)/(6*c0)
            if verify_saved_terminal_rows:
                for distance in range(len(obj.p)):
                    bound = (71-5*distance)//12
                    if bound >= 0 and obj.centered_slice(distance, bound):
                        raise AssertionError(f"saved P terminal rows fail at distance {distance}")
            print(f"resumed={path} branch={obj.tag} depth={len(obj.p)-1} parameters={obj.nparameters} constraints={len(obj.constraints)} streaming=true", flush=True)
            return obj
        finally:
            reader.close()

    def materialize_constraint_polynomial(self, polynomial):
        """Decode one lazy constraint and apply its queued substitutions."""
        if not isinstance(polynomial, LazyConstraintPolynomial):
            return polynomial
        reader = JsonStream(polynomial.path)
        try:
            reader.seek_absolute(polynomial.start)
            value = _decode_sparse_checkpoint_polynomial(
                reader, self.ring, self.field,
                self._checkpoint_old_names,
                self._checkpoint_new_indices,
                self._checkpoint_quadratic,
            )
        finally:
            reader.close()
        for replacements in polynomial.transformations:
            value = value.compose(replacements)
        return value

    def const(self, value):
        return self.ring.ground_new(self.field.convert(value))

    def polynomial_from_terms(self, terms):
        """Build from validated dense exponents without recopying FLINT input."""
        direct = getattr(self.ring, "from_prevalidated_dict", None)
        return direct(terms) if direct is not None else self.ring.from_dict(terms)

    def polynomial_linear_combination(self, summands):
        """Combine scalar multiples of sparse polynomials exactly.

        The rational FLINT backend performs scalar multiplication and addition
        inside C, avoiding Python iteration over every monomial of the large
        coefficient slices. The SymPy path keeps the one-construction sparse
        accumulation, which avoids repeated normalization there.
        """
        if getattr(self.ring, "backend_name", None) == "flint-mpoly":
            result = None
            for polynomial, scalar in summands:
                if not scalar or not polynomial:
                    continue
                # Free parameters may still be represented by the FLINT
                # ring's lightweight FGen wrapper.  Normalize them to FPoly
                # before using the backend-only scale operation.
                polynomial = self.ring.as_poly(polynomial)
                scaled = polynomial.scale(scalar)
                if result is None:
                    result = scaled
                elif result.shift == scaled.shift:
                    result.data.iadd(scaled.data)
                else:
                    result = result + scaled
            return self.ring.zero if result is None else result

        terms = {}
        for polynomial, scalar in summands:
            if not scalar or not polynomial:
                continue
            if not hasattr(polynomial, "items") and hasattr(polynomial, "_poly"):
                polynomial = polynomial._poly()
            for powers, coefficient in polynomial.items():
                value = coefficient * scalar
                if not value:
                    continue
                previous = terms.get(powers)
                if previous is None:
                    terms[powers] = value
                else:
                    total = previous + value
                    if total:
                        terms[powers] = total
                    else:
                        del terms[powers]
        return self.polynomial_from_terms(terms) if terms else self.ring.zero

    def _refresh_affine_ring(self):
        """Keep the FLINT coefficient subring aligned with active parameters."""
        if getattr(self.ring, "backend_name", None) == "flint-mpoly":
            self.affine_ring = self.ring.prefix(1+self.nparameters)
        else:
            self.affine_ring = None

    def allocate(self, count):
        begin = self.nparameters
        new_count = begin + count
        if new_count > len(self.parameters):
            raise ValueError("parameter reserve exhausted")
        self.nparameters = new_count
        self._refresh_affine_ring()
        return self.parameters[begin:self.nparameters]

    def derivative(self, value):
        return value.diff(self.t)

    def split_affine(self, expression):
        """Extract all rows, including forcing-only rows, directly by monomial.

        Accumulate each t-slice of the forcing part as a sparse term dictionary
        and construct that polynomial once.  Adding one one-term FLINT mpoly at
        a time repeatedly normalizes the entire growing polynomial and becomes
        quadratic in the number of terms on the large distance-16 slice.
        """
        if getattr(self.ring, "backend_name", None) == "flint-mpoly":
            return self._split_affine_flint_projected(expression)
        rows, forcing_terms = {}, {}
        for powers, coefficient in expression.items():
            degree = powers[0]
            new_powers = powers[self.new_start:]
            new_degree = sum(new_powers)
            if new_degree:
                if new_degree != 1 or any(powers[1:self.new_start]):
                    raise AssertionError("new-coefficient matrix is not constant affine")
                j = new_powers.index(1)
                row = rows.setdefault(degree, [self.field.zero]*len(self.unknowns))
                row[j] += coefficient
            else:
                monomial = (0, *powers[1:])
                terms = forcing_terms.setdefault(degree, {})
                if monomial in terms:
                    terms[monomial] += coefficient
                else:
                    terms[monomial] = coefficient
        forcing = {degree: self.polynomial_from_terms(terms)
                   for degree, terms in forcing_terms.items()}
        degrees = sorted(set(rows) | set(forcing))
        return degrees, [rows.get(d, [self.field.zero]*len(self.unknowns)) for d in degrees], [forcing.get(d, self.ring.zero) for d in degrees]

    def _split_affine_flint_projected(self, expression):
        """Split a FLINT affine equation in one pass, flushing each t-slice.

        A monomial either contributes to the forcing term (no fresh unknown)
        or to one constant matrix entry (one fresh unknown to the first power
        and no carried parameter).  Native terms are ordered with t first, so
        each forcing slice can be built and released as soon as its t-degree
        changes.  This avoids retaining a Python dictionary for every forcing
        monomial and avoids a derivative/projection for each fresh unknown.
        """
        affine_ring = self.affine_ring
        if affine_ring is None:
            raise AssertionError("FLINT affine subring was not initialized")
        if affine_ring.symbols != self.ring.symbols[:1+self.nparameters]:
            raise AssertionError("FLINT affine subring is stale for active parameters")
        rows = {}
        forcing = {}
        active_count = self.nparameters
        current_degree = None
        current_terms = {}
        seen_degrees = set()
        for powers, scalar in expression.native_items():
            degree = int(powers[0])
            if current_degree is None:
                current_degree = degree
                seen_degrees.add(degree)
            elif degree != current_degree:
                if degree in seen_degrees:
                    raise AssertionError(
                        "FLINT native terms are not grouped by t-degree")
                if current_terms:
                    small_slice = affine_ring.from_native_dict(current_terms)
                    forcing[current_degree] = small_slice.embed_into(self.ring)
                    current_terms.clear()
                current_degree = degree
                seen_degrees.add(degree)

            new_powers = powers[self.new_start:]
            new_degree = sum(new_powers)
            if new_degree:
                if new_degree != 1 or any(powers[1:self.new_start]):
                    raise AssertionError(
                        "new-coefficient matrix is not constant affine")
                column = new_powers.index(1)
                row = rows.setdefault(
                    degree, [self.field.zero]*len(self.unknowns))
                row[column] += scalar
                continue
            if any(powers[1+active_count:self.new_start]):
                raise AssertionError(
                    "forcing term uses an unallocated carried parameter")
            monomial = (0, *powers[1:1+active_count])
            current_terms[monomial] = scalar

        if current_terms:
            small_slice = affine_ring.from_native_dict(current_terms)
            forcing[current_degree] = small_slice.embed_into(self.ring)
            current_terms.clear()

        degrees = sorted(set(rows) | set(forcing))
        return (
            degrees,
            [rows.get(d, [self.field.zero]*len(self.unknowns))
             for d in degrees],
            [forcing.get(d, self.ring.zero) for d in degrees],
        )

    def _split_affine_flint_terms_reference(self, expressions):
        """Split streamed affine summands without a campaign-sized term dict.

        Each summand is consumed in FLINT's t-lex order.  We retain only one
        t-slice of native terms at a time, convert that slice in the active
        carried-parameter ring, and merge it into the corresponding forcing
        slice.  Fresh x-unknowns contribute only constant matrix entries.
        """
        affine_ring = self.affine_ring
        if affine_ring is None:
            raise AssertionError("FLINT affine subring was not initialized")
        rows, forcing = {}, {}
        active_forcing = {}
        active_count = self.nparameters

        def retain_slice(degree, terms):
            if degree is None or not terms:
                return
            value = affine_ring.from_native_dict(terms)
            if not value:
                return
            previous = forcing.get(degree)
            value = value if previous is None else previous + value
            if value:
                forcing[degree] = value
            else:
                forcing.pop(degree, None)

        for expression in expressions:
            expression = self.ring.as_poly(expression)
            data = expression.data
            shift = expression.shift
            current_degree = None
            current_terms = {}
            seen_degrees = set()
            # Read each FLINT monomial exactly once and use its native tuple
            # directly.  ``native_items`` used to rebuild a second full
            # exponent tuple for every term; at distance 23 this scan can
            # exceed half a million terms in a 300-generator context.
            # Keep this indexed traversal bounded-memory: python-flint's
            # ``terms()`` implementation materializes both ``monoms()`` and
            # ``coeffs()`` before yielding, which is unsafe for these very
            # large, high-dimensional residuals.
            for term_index in range(len(data)):
                powers = data.monomial(term_index)
                scalar = data.coefficient(term_index)
                degree = int(powers[0])+shift
                if current_degree is None:
                    current_degree = degree
                    seen_degrees.add(degree)
                elif degree != current_degree:
                    if degree in seen_degrees:
                        raise AssertionError(
                            "FLINT native terms are not grouped by t-degree")
                    retain_slice(current_degree, current_terms)
                    current_terms = {}
                    current_degree = degree
                    seen_degrees.add(degree)
                new_powers = powers[self.new_start:]
                new_degree = sum(new_powers)
                if new_degree:
                    if new_degree != 1 or any(powers[1:self.new_start]):
                        raise AssertionError(
                            "new-coefficient matrix is not constant affine")
                    column = new_powers.index(1)
                    row = rows.setdefault(
                        degree, [self.field.zero]*len(self.unknowns))
                    row[column] += scalar
                    continue
                if any(powers[1+active_count:self.new_start]):
                    raise AssertionError(
                        "forcing term uses an unallocated carried parameter")
                monomial = (0, *powers[1:1+active_count])
                current_terms[monomial] = scalar
            retain_slice(current_degree, current_terms)

        rows = {degree: row for degree, row in rows.items() if any(row)}
        degrees = sorted(set(rows) | set(forcing))
        return (
            degrees,
            [rows.get(d, [self.field.zero]*len(self.unknowns))
             for d in degrees],
            [forcing[d].embed_into(self.ring) if d in forcing
             else self.ring.zero for d in degrees],
        )

    def _split_affine_flint_terms(self, expressions):
        """Split tagged FLINT summands into forcing and constant affine rows.

        For a summand containing fresh unknowns, its terms linear in a fresh
        variable can also contain carried parameters.  Collect those
        coefficients in the active ring and combine them across all summands
        before checking that the resulting affine matrix is over the ground
        field.  Native FLINT derivatives and context projections extract these
        coefficients without scanning every sparse monomial in Python.  The
        derivative projection is checked exactly to ensure each summand is
        affine in the fresh variables.  The earlier pair of independent
        projections lost coefficient cancellation and rejected valid inputs.
        """
        if self.affine_ring is None:
            raise AssertionError("FLINT affine subring was not initialized")

        rows, forcing = {}, {}
        active_forcing = {}
        matrix_coefficients = {}
        active_count = self.nparameters
        active_ring = self.affine_ring
        if active_ring.symbols != self.ring.symbols[:1+active_count]:
            raise AssertionError("active FLINT ring does not match the parameter prefix")
        fresh_names = self.ring.symbols[self.new_start:]
        if len(fresh_names) != len(self.unknowns):
            raise AssertionError("fresh FLINT generators and unknowns differ")

        def retain_slice(degree, value):
            if degree is None or not value:
                return
            if value.ring is active_ring:
                accumulator = active_forcing
            elif value.ring is not self.ring:
                raise AssertionError("forcing slice has an unexpected FLINT context")
            else:
                accumulator = forcing
            previous = accumulator.get(degree)
            if previous is None:
                if value:
                    accumulator[degree] = value
                return
            if previous.shift != value.shift:
                raise AssertionError("forcing slices at one t-degree have different shifts")
            # python-flint exposes fmpq_mpoly.iadd as a native in-place add.
            # The ordinary `previous + value` keeps both large operands alive
            # while allocating a third full polynomial; this was the d23
            # allocation-failure hot spot.  The accumulator is owned only by
            # this local dictionary, so mutating its FLINT payload is safe.
            previous.data.iadd(value.data)
            if not previous:
                accumulator.pop(degree, None)

        def retain_matrix_slice(degree, column, value):
            if not value:
                return
            key = (degree, column)
            previous = matrix_coefficients.get(key)
            if previous is None:
                matrix_coefficients[key] = value
            else:
                if previous.shift != value.shift:
                    raise AssertionError("matrix slices at one t-degree have different shifts")
                previous.data.iadd(value.data)
                if not previous:
                    matrix_coefficients.pop(key, None)

        for summand_index, (expression, has_fresh_unknown) in enumerate(expressions):
            source_ring = getattr(expression, "ring", None)
            if source_ring is active_ring:
                if has_fresh_unknown:
                    raise AssertionError("active-ring correction tagged as fresh")
                for degree, source_slice in expression.first_variable_coefficients():
                    retain_slice(degree, source_slice)
                continue
            else:
                expression = self.ring.as_poly(expression)
                source_ring = self.ring

            # FLINT uses degree -1 in every variable for zero; skip before
            # applying the reserved-parameter guard.
            if expression.data.is_zero():
                continue
            variable_degrees = expression.data.degrees()
            reserved_start = 1 + active_count
            if any(variable_degrees[reserved_start:self.new_start]):
                offending = [
                    (index, self.ring.symbols[index], int(variable_degrees[index]))
                    for index in range(reserved_start, min(self.new_start, len(variable_degrees)))
                    if variable_degrees[index]
                ]
                raise AssertionError(
                    "affine summand uses an unallocated carried parameter "
                    f"(summand={summand_index}, active_count={active_count}, "
                    f"reserved={offending[:12]}, terms={len(expression)}, "
                    f"ring_width={len(source_ring.symbols)})")
            fresh_indices = [
                index for index in range(self.new_start, len(variable_degrees))
                if variable_degrees[index]
            ]
            if fresh_indices and not has_fresh_unknown:
                raise AssertionError(
                    "known-only summand contains a fresh unknown")

            # Separate the forcing part once, in FLINT, and read its t-slices
            # in the smaller active context.
            forcing_expression = (
                expression if expression.ring is active_ring
                else expression.project_to_ring(active_ring)
            )
            for degree, forcing_slice in forcing_expression.first_variable_coefficients():
                retain_slice(degree, forcing_slice)

            # A fresh coefficient is the derivative in that variable,
            # evaluated at zero in all fresh variables by projection to the
            # active context.  If the derivative changes under that
            # projection, the summand was nonlinear in the fresh unknowns.
            # This exact guard catches both powers and mixed products.
            for variable_index in fresh_indices:
                derivative = expression.data.derivative(variable_index)
                coefficient_data = derivative.project_to_context(active_ring.ctx)
                if coefficient_data.project_to_context(self.ring.ctx) != derivative:
                    raise AssertionError(
                        "new-coefficient matrix is not constant affine")
                coefficient_expression = type(expression)(
                    active_ring, coefficient_data, expression.shift)
                column = variable_index - self.new_start
                for degree, coefficient_slice in coefficient_expression.first_variable_coefficients():
                    retain_matrix_slice(degree, column, coefficient_slice)

        for (degree, column), coefficient in matrix_coefficients.items():
            data = coefficient.data
            if not data.is_zero():
                if len(data) != 1 or any(data.monomial(0)):
                    raise AssertionError(
                        "fresh-variable coefficient still depends on carried "
                        f"parameters at t^{degree}, column {column} "
                        f"(terms={len(data)})"
                    )
                row = rows.setdefault(
                    degree, [self.field.zero]*len(self.unknowns))
                row[column] += data.coefficient(0)

        # Merge active-context forcing only after all summands have been
        # combined.  Embedding each contribution before accumulation can copy
        # the same large coefficient slice many times; one embedding per
        # surviving t-degree preserves exact cancellation and reduces both
        # native work and peak memory.
        for degree, active_value in active_forcing.items():
            embedded = active_value.embed_into(self.ring)
            previous = forcing.get(degree)
            if previous is None:
                if embedded:
                    forcing[degree] = embedded
                continue
            if previous.shift != embedded.shift:
                raise AssertionError("forcing slices at one t-degree have different shifts")
            previous.data.iadd(embedded.data)
            if not previous:
                forcing.pop(degree, None)

        rows = {degree: row for degree, row in rows.items() if any(row)}
        degrees = sorted(set(rows) | set(forcing))
        return (
            degrees,
            [rows.get(d, [self.field.zero]*len(self.unknowns))
             for d in degrees],
            [forcing.get(d, self.ring.zero) for d in degrees],
        )

    def record_univariate_constraints(self, expression, tag,
                                      exact_duplicate_candidates=()):
        """Retain every t-coefficient of a polynomial compatibility equation."""
        if getattr(self.ring, "backend_name", None) == "flint-mpoly":
            # The FLINT context uses lex order with t as the first variable,
            # so equal t-degrees are contiguous in its term iterator.  Build
            # each slice once instead of rescanning/materializing the full
            # large polynomial once for every t-degree.
            source_ring = getattr(expression, "ring", None)
            if source_ring is None:
                expression = self.ring.as_poly(expression)
                source_ring = self.ring
            if source_ring is not self.ring and source_ring is not self.affine_ring:
                raise TypeError("constraint polynomial is in an unexpected FLINT ring")
            current_degree = None
            terms = {}

            def retain_slice(degree, slice_terms):
                if degree is None or not slice_terms:
                    return
                value = source_ring.from_native_dict(slice_terms)
                if source_ring is not self.ring:
                    value = value.embed_into(self.ring)
                if value:
                    if not any(value == old
                               for old in exact_duplicate_candidates):
                        self.constraints.append(
                            (f"{tag}:t^{degree}", value))

            for powers, coefficient in expression.native_items():
                degree = int(powers[0])
                if current_degree is not None and degree != current_degree:
                    retain_slice(current_degree, terms)
                    terms = {}
                current_degree = degree
                terms[(0, *powers[1:])] = coefficient

            retain_slice(current_degree, terms)
            return

        for degree in sorted({powers[0] for powers in expression}):
            value = self.local_coefficient(expression, degree)
            if value:
                if any(value == old for old in exact_duplicate_candidates):
                    continue
                self.constraints.append((f"{tag}:t^{degree}", value))

    def affine_solve(self, equations, nvariables, step, reuse_free=None,
                     streamed=False):
        rows, forcing, row_tags = [], [], []
        for side, expression in enumerate(equations):
            if streamed:
                if getattr(self.ring, "backend_name", None) == "flint-mpoly":
                    split = self._split_affine_flint_terms(expression)
                else:
                    split = self.split_affine(sum(expression, self.ring.zero))
                degrees, local_rows, local_forcing = split
            else:
                degrees, local_rows, local_forcing = self.split_affine(expression)
            rows.extend(row[:nvariables] for row in local_rows)
            forcing.extend(local_forcing)
            row_tags.extend((side, degree) for degree in degrees)
            assert all(not any(row[nvariables:]) for row in local_rows)
        # Entries already belong to K.  from_list would call K(entry), which
        # tries to coerce an existing ANP through QQ in the quadratic branch.
        matrix = DomainMatrix(rows, (len(rows), nvariables), self.field)
        _, pivot_cols = matrix.rref()
        _, pivot_rows = matrix.transpose().rref()
        rank = len(pivot_cols)
        assert len(pivot_rows) == rank
        free_cols = [j for j in range(nvariables) if j not in pivot_cols]
        free_values = (self.allocate(len(free_cols)) if reuse_free is None
                       else [reuse_free[j] for j in free_cols])
        self.last_free = tuple(free_values)
        pivot = matrix.extract(pivot_rows, pivot_cols)
        inverse = pivot.inv().to_list()
        solution = [self.ring.zero]*nvariables
        for j, value in zip(free_cols, free_values):
            solution[j] = value
        pivot_forcing = []
        for i in pivot_rows:
            summands = [(forcing[i], 1)]
            for j in free_cols:
                if rows[i][j]:
                    summands.append((solution[j], rows[i][j]))
            pivot_forcing.append(self.polynomial_linear_combination(summands))
        for j, inverse_row in zip(pivot_cols, inverse):
            summands = [(rhs, coefficient)
                        for coefficient, rhs in zip(inverse_row, pivot_forcing)
                        if coefficient and rhs]
            solution[j] = -self.polynomial_linear_combination(summands)
        residuals = []
        pivot_row_set = set(pivot_rows)
        for row_index, (tag, row, rhs) in enumerate(zip(row_tags, rows, forcing)):
            # These basis rows were used to solve the pivot variables, so their
            # residuals vanish identically. Only dependent rows can impose new
            # compatibility constraints on the carried parameters.
            if row_index in pivot_row_set:
                continue
            summands = [(rhs, 1)]
            for coefficient, entry in zip(row, solution):
                if coefficient and entry:
                    summands.append((entry, coefficient))
            value = self.polynomial_linear_combination(summands)
            if value:
                residuals.append((tag, value))
        record = dict(step=step, rows=len(rows), variables=nvariables, rank=rank,
                      fibre=len(free_cols), carried_parameters=self.nparameters,
                      nonzero_residuals=len(residuals), all_rows_verified=not residuals)
        self.receipt.append(record)
        print(json.dumps(dict(branch=self.tag, **record)), flush=True)
        if residuals:
            for tag, value in residuals[:3]:
                # Converting a high-degree sparse polynomial over an
                # algebraic number field to a SymPy expression can recurse
                # deeply in SymPy's coefficient gcd/simplification path.
                # The exact polynomial is retained below as a constraint;
                # its tag and sparse term count are sufficient diagnostics.
                print(f"residual={tag} terms={len(value)} retained=True", flush=True)
            for tag, value in residuals:
                self.constraints.append((f"{step}:{tag[0]}:t^{tag[1]}", value))
        return solution

    def replace_parameters(self, replacements):
        self.S = self.S.compose(replacements)
        if self.B is not None:
            self.B = self.B.compose(replacements)
        self.kappa = self.kappa.compose(replacements)
        self.kappa2 = self.kappa2.compose(replacements)
        self.kappa3 = self.kappa3.compose(replacements)
        self.uk = {key: value.compose(replacements) for key, value in self.uk.items()}
        self.p = [value.compose(replacements) for value in self.p]
        self.q = [value.compose(replacements) for value in self.q]
        self.w = [value.compose(replacements) for value in self.w]
        replaced_constraints = []
        for tag, value in self.constraints:
            if isinstance(value, LazyConstraintPolynomial):
                value = value.compose_lazily(replacements)
            elif getattr(self.ring, "backend_name", None) == "flint-mpoly":
                if value.ring is not self.ring:
                    value = value.embed_into(self.ring)
                value = value.compose(replacements)
            else:
                value = value.compose(replacements)
            replaced_constraints.append((tag, value))
        self.constraints = replaced_constraints
        for name in ("center_a", "center_b"):
            if hasattr(self, name):
                setattr(self, name, getattr(self, name).compose(replacements))

    def impose_congruence(self, expression, modulus, targets, tag):
        """Solve a necessary congruence in a selected previous affine fibre."""
        remainder = expression.rem(modulus)
        renaming = dict(zip(targets, self.unknowns))
        equation = remainder.compose(renaming)
        previous_constraint_count = len(self.constraints)
        solution = self.affine_solve((equation,), len(targets), tag, reuse_free=targets)
        # affine_solve has already retained its exact compatibility residuals.
        # After the substitutions below, some or all t-slices of the original
        # congruence can be literally identical to those rows.  Keep every
        # nonduplicate slice, but do not carry a second copy of an identical
        # polynomial through later lifts and checkpoint serialization.
        duplicate_candidates = [value for _, value
                               in self.constraints[previous_constraint_count:]]
        replacements = dict(zip(targets, solution))
        # A rank-zero solve on an already-free fibre returns each target
        # parameter unchanged.  Composing that identity map through every
        # carried polynomial is mathematically inert but can be very costly
        # once the sparse lift has hundreds of thousands of terms.
        identity = all(replacements[target] == target for target in targets)
        if identity:
            remainder_after = remainder
        else:
            self.replace_parameters(replacements)
            remainder_after = expression.compose(replacements).rem(modulus)
        if remainder_after:
            self.record_univariate_constraints(
                remainder_after, f"{tag}:postsolve",
                exact_duplicate_candidates=duplicate_candidates,
            )

    def impose_tuned_contact(self, distance):
        """Overapproximate the all-roots-tuned branch, retaining zero cases.

        K6/K7 at tuned roots force h|(3B-S²), h|(27U3-S³).
        K8 then forces rad(h)|(S*Delta/h-81U4).  These are applied as
        necessary congruences; S is never inverted or assumed nonzero here.
        """
        targets = self.last_free
        if distance == 2:
            self.impose_congruence(3*self.B-self.S**2, self.h, targets, "tuned-D")
        elif distance == 3:
            self.impose_congruence(27*self.uk[3]-self.S**3, self.h, targets, "tuned-E")
        elif distance == 4:
            d = 3*self.B-self.S**2
            delta = 27*self.uk[3]+2*self.S**3-9*self.S*self.B
            assert not d.rem(self.h) and not delta.rem(self.h)
            radical = self.h.exquo(self.h.gcd(self.h.diff(self.t)))
            contact = self.S*delta.exquo(self.h)-81*self.uk[4]
            self.impose_congruence(contact, radical, targets, "tuned-K8")
            d = 3*self.B-self.S**2
            delta = 27*self.uk[3]+2*self.S**3-9*self.S*self.B
            n8 = ((self.S*delta-81*self.uk[4]*self.h)**2-d*delta**2
                  +6*self.S*d**2*delta-3*self.S**2*d**3+3*d**4
                  -243*d**2*self.uk[4]*self.h)
            # The omitted U5 term contains h²*Delta, hence is divisible by h³.
            assert not n8.rem(self.h**3)

    def impose_root_values(self, distance):
        """Unconditional K6/K7 value constraints, including deeper roots."""
        # Root-value divisibility enters only at K6/K7 (distances 2 and 3).
        # Avoid recomputing the squarefree radical on later resumed blocks:
        # the sparse multivariate ring's generic gcd routine can recurse on an
        # algebraic coefficient field even though those blocks add no such row.
        if distance not in (2, 3):
            return
        radical = self.h.exquo(self.h.gcd(self.h.diff(self.t)))
        if distance == 2:
            self.impose_congruence(3*self.B-self.S**2, radical, self.last_free,
                                   "K67-root-D")
        else:
            self.impose_congruence(27*self.uk[3]-self.S**3, radical, self.last_free,
                                   "K67-root-E")
            assert not self.sixth_numerator().rem(self.h)

    def impose_deep_rootwise_dd(self, distance):
        """Use the exhaustive DD divisor consequences for the 2+1+1 branch.

        The exact K6--K10 rootwise reduction, together with the terminal
        condition S(1)=0, leaves only the deep/deep profile at the two simple
        roots of R=t^2-3t+3. It gives

            (t-1)^4 R^2 | U2,  (t-1) R | U3,  R | U4.

        These are necessary divisibilities for every surviving solution,
        including U1=0. Solving them on the current affine fibre keeps the
        later exact lift smaller; this is a parameterization aid, not a new
        exclusion. The caller must also enable the exact terminal and K6/K7
        root-value reductions used by the proof of the DD split.
        """
        if self.tag != "2+1+1":
            raise ValueError("the DD divisor reduction is specific to 2+1+1")
        R = self.t**2 - 3*self.t + 3
        expected_h = (self.t-1)**2 * R
        if self.h != expected_h:
            raise AssertionError("unexpected normalized 2+1+1 carrier")
        if distance == 2:
            divisor, expression = (self.t-1)**4 * R**2, self.uk[2]
        elif distance == 3:
            divisor, expression = (self.t-1) * R, self.uk[3]
        elif distance == 4:
            divisor, expression = R, self.uk[4]
        else:
            return
        self.impose_congruence(expression, divisor, self.last_free,
                               f"DD-rootwise-U{distance}")

    def sixth_numerator(self):
        d = 3*self.B-self.S**2
        delta = 27*self.uk[3]+2*self.S**3-9*self.S*self.B
        return delta**2-3*d**3

    def local_coefficient(self, polynomial, degree):
        return self.polynomial_from_terms({(0, *powers[1:]): coefficient
                                           for powers, coefficient in polynomial.items()
                                           if powers[0] == degree})

    def local_slice(self, distance, bound):
        """Taylor-expand X^-15 P(X^5,X^-1(1+u)) in u=t.

        A global band monomial t^i becomes (1+u)^(5i-15+s).
        The exponent is its original nonnegative y exponent, so this
        transformation uses the actual source monomials and no Laurent cut.
        """
        terms = {}
        for powers, coefficient in self.p[distance].items():
            exponent = 5*powers[0]-15+distance
            if exponent < 0:
                raise AssertionError("P monomial is outside original polynomial support")
            for degree in range(min(exponent, bound)+1):
                monomial = (degree, *powers[1:])
                value = coefficient*comb(exponent, degree)
                previous = terms.get(monomial)
                if previous is None:
                    if value:
                        terms[monomial] = value
                else:
                    total = previous+value
                    if total:
                        terms[monomial] = total
                    else:
                        del terms[monomial]
        return self.polynomial_from_terms(terms)

    def local_q_slice(self, distance, bound):
        """Taylor-expand X^-25 Q(X^5,X^-1(1+u)) in the local variable u."""
        terms = {}
        low = -(-(25-distance)//5)
        high = min(25, (125-distance)//5)
        for powers, coefficient in self.q[distance].items():
            if powers[0] < low or powers[0] > high:
                # Those coefficients are separately retained as source-support
                # equations; work with their projection on that exact locus.
                continue
            exponent = 5*powers[0]-25+distance
            if exponent < 0:
                raise AssertionError("Q monomial is outside original polynomial support")
            for degree in range(min(exponent, bound)+1):
                monomial = (degree, *powers[1:])
                value = coefficient*comb(exponent, degree)
                previous = terms.get(monomial)
                if previous is None:
                    if value:
                        terms[monomial] = value
                else:
                    total = previous+value
                    if total:
                        terms[monomial] = total
                    else:
                        del terms[monomial]
        return self.polynomial_from_terms(terms)

    def local_q_slice_derivative(self, distance, bound, order):
        """Return ``d^order/du^order local_q_slice`` through ``u^bound``.

        Construct the differentiated Taylor coefficients directly, rather
        than expanding a larger slice and repeatedly differentiating a large
        sparse multivariate polynomial.
        """
        source_bound = bound+order
        low = -(-(25-distance)//5)
        high = min(25, (125-distance)//5)
        source = self.q[distance]
        # Collect every requested output slice during one source-term pass.
        # The previous implementation copied the full polynomial to a dict
        # and rescanned it once for each of the (bound+1) output degrees.
        slices = {}
        for degree in range(order, source_bound+1):
            slices[degree] = {}
        for powers, coefficient in source.items():
            if powers[0] < low or powers[0] > high:
                continue
            exponent = 5*powers[0]-25+distance
            if exponent < 0:
                raise AssertionError("Q monomial is outside original polynomial support")
            for degree in range(order, source_bound+1):
                if degree > exponent:
                    continue
                output_degree = degree-order
                factor = comb(exponent, degree)*factorial(degree)//factorial(output_degree)
                monomial = (output_degree, *powers[1:])
                value = coefficient*factor
                terms = slices[degree]
                previous = terms.get(monomial)
                if previous is None:
                    if value:
                        terms[monomial] = value
                else:
                    total = previous+value
                    if total:
                        terms[monomial] = total
                    else:
                        del terms[monomial]
        out = self.ring.zero
        for degree in range(order, source_bound+1):
            terms = slices[degree]
            if terms:
                out += self.polynomial_from_terms(terms)
        return out

    def impose_q_source_support(self, distance):
        """Retain every Q band coefficient outside the original rectangle."""
        low = -(-(25-distance)//5)
        high = min(25, (125-distance)//5)
        # Group every forbidden t-slice in one pass.  Calling
        # ``local_coefficient`` once per bad degree rescans the full sparse
        # Q polynomial repeatedly; at deep lifts this polynomial can contain
        # hundreds of thousands of terms.  The grouped dictionaries are
        # exactly the same coefficient slices, with no relaxation.
        bad_slices = {}
        for powers, coefficient in self.q[distance].items():
            degree = powers[0]
            if degree < low or degree > high:
                bad_slices.setdefault(degree, {})[(0, *powers[1:])] = coefficient
        before = len(self.constraints)
        for degree, terms in sorted(bad_slices.items()):
            value = self.polynomial_from_terms(terms)
            if value:
                self.constraints.append((f"Q-source-support-distance{distance}:t^{degree}", value))
        if len(self.constraints) > before:
            print(f"branch={self.tag} Q-source-support-distance={distance} "
                  f"equations={len(self.constraints)-before}", flush=True)

    def centered_slice(self, distance, bound):
        result = self.ring.zero
        modulus = self.t**(bound+1)
        for i in range(distance+1):
            deficit = distance-i
            for n in range(deficit//2+1):
                m = deficit-2*n
                derivative = self.local_slice(i, bound+m+n)
                for _ in range(m+n):
                    derivative = derivative.diff(self.t)
                a_power = self.ring.one if m == 0 else self.center_a**m
                b_power = self.ring.one if n == 0 else self.center_b**n
                result += (derivative*a_power*b_power
                           *self.field.convert(sp.Rational(1, factorial(m)*factorial(n)))).rem(modulus)
        return result.rem(modulus)

    def centered_q_slice(self, distance, bound):
        result = self.ring.zero
        modulus = self.t**(bound+1)
        for i in range(distance+1):
            deficit = distance-i
            for n in range(deficit//2+1):
                m = deficit-2*n
                a_power = self.ring.one if m == 0 else self.center_a**m
                b_power = self.ring.one if n == 0 else self.center_b**n
                # Skip a zero chart summand before expanding its potentially
                # large local-Q derivative and sparse polynomial products.
                if not a_power or not b_power:
                    continue
                derivative = self.local_q_slice_derivative(i, bound, m+n)
                result += (derivative*a_power*b_power
                           *self.field.convert(sp.Rational(1, factorial(m)*factorial(n))))
        return result.rem(modulus)

    def impose_terminal_q(self, distance):
        bound = (119-5*distance)//12
        if bound < 0:
            return
        self.impose_congruence(self.centered_q_slice(distance, bound),
                               self.t**(bound+1), self.last_free,
                               f"terminal-Q-distance{distance}")

    def initialize_terminal_chart(self):
        """Necessary terminal-chain constraint at the distinguished root 1.

        Source theorem and its scope are recorded in
        f2_terminal_chain_local_selectors_20261018.md.  Type-III collapses
        are retained as u -> u+a*w+b*w²; only nonzero scalar units divide.
        """
        targets = tuple(self.parameters[:self.nparameters])
        self.impose_congruence(self.S, self.t-1, targets, "terminal-S(1)")
        self.local_C = self.local_coefficient(self.local_slice(0, 6), 6)
        zero_powers = (0,)*len(self.ring.gens)
        if set(self.local_C) != {zero_powers}:
            raise AssertionError("terminal leading coefficient is not a nonzero scalar")
        c0 = self.local_C[zero_powers]
        self.center_a = -self.local_coefficient(self.local_slice(1, 5), 5)/(6*c0)
        self.center_b = self.ring.zero
        assert not self.centered_slice(0, 5)
        assert not self.centered_slice(1, 5)
        print(f"branch={self.tag} terminal_chart_initialized a_terms={len(self.center_a)}", flush=True)

    def impose_terminal_chart(self, distance):
        targets = self.last_free
        if distance == 2:
            # The b*A0' term starts at u5, so lower coefficients are the
            # exact compatibility equations before determining b.
            self.impose_congruence(self.centered_slice(2, 4), self.t**5,
                                   targets, "terminal-P-distance2")
            zero_powers = (0,)*len(self.ring.gens)
            c0 = self.local_C[zero_powers]
            self.center_b = -self.local_coefficient(self.centered_slice(2, 5), 5)/(6*c0)
            assert not self.centered_slice(2, 5)
            print(f"branch={self.tag} terminal_center_b_terms={len(self.center_b)}", flush=True)
        else:
            bound = (71-5*distance)//12
            if bound >= 0:
                expression = self.centered_slice(distance, bound)
                self.impose_congruence(expression, self.t**(bound+1), targets,
                                       f"terminal-P-distance{distance}")
                assert not self.centered_slice(distance, bound)

    def substitute_affine(self, expression, solution):
        forcing_terms = {}
        variable_terms = [{} for _ in self.unknowns]
        for powers, coefficient in expression.items():
            new_powers = powers[self.new_start:]
            if any(new_powers):
                if sum(new_powers) != 1:
                    raise AssertionError("substitution is not affine")
                j = new_powers.index(1)
                monomial = (*powers[:self.new_start], *((0,)*len(self.unknowns)))
                variable_terms[j][monomial] = coefficient
            else:
                forcing_terms[powers] = coefficient
        out = self.polynomial_from_terms(forcing_terms)
        for terms, value in zip(variable_terms, solution):
            if terms and value:
                out += self.polynomial_from_terms(terms)*value
        return out

    def series_power_coefficients(self, values, exponent, max_distance):
        """Return all coefficients through ``max_distance`` in one pass.

        The former callers recomputed every lower coefficient independently.
        This vector form uses the same exact convolution recurrence while
        avoiding that repeated polynomial multiplication.
        """
        coeffs = [self.ring.one] + [self.ring.zero]*max_distance
        for _ in range(exponent):
            coeffs = [sum((coeffs[j]*values[i-j] for j in range(i+1)), self.ring.zero)
                      for i in range(max_distance+1)]
        return coeffs

    def selected_series_power_coefficients(
        self, values, requests, max_distance, *, polynomial_ring=None
    ):
        """Compute several powers together and retain only needed coefficients.

        ``requests`` maps an exponent to ``(prefix_degree, indices)``.  The
        result for that exponent is a sparse ``degree -> coefficient`` map.
        Work backward from those requests to compute only the coefficient
        degrees needed by later powers.  In particular, the final power no
        longer constructs every coefficient from degree zero through
        ``max_distance`` when only its top requested coefficient is used.
        """
        if not requests:
            return {}

        requested = {}
        needed = {}
        max_exponent = max(requests)
        for exponent, (prefix_degree, indices) in requests.items():
            degrees = set(indices)
            if prefix_degree is not None:
                degrees.update(range(prefix_degree+1))
            if any(degree < 0 or degree > max_distance for degree in degrees):
                raise IndexError("requested power coefficient is outside the computed range")
            requested[exponent] = degrees
            needed.setdefault(exponent, set()).update(degrees)

        # c[e,d] = sum(c[e-1,j] * values[d-j], j=0..d).  Therefore a
        # requested degree d at exponent e needs the complete prefix through d
        # at exponent e-1; propagate these dependencies downward.
        for exponent in range(max_exponent, 1, -1):
            degrees = needed.get(exponent, set())
            if degrees:
                needed.setdefault(exponent-1, set()).update(
                    range(max(degrees)+1)
                )

        # In the Q-power calculation we need one coefficient of Q^3 and a
        # sparse subset of Q^2.  The generic recurrence below retains the
        # entire Q^2 prefix until it has formed Q^3, which can hold many large
        # sparse FLINT polynomials at once.  Stream that one Q^3 coefficient
        # through the symmetric Q^2 coefficients instead, retaining only the
        # requested Q^2 outputs.  The branch is deliberately narrow; all
        # other power requests continue through the general implementation.
        if (max_exponent == 3 and set(requested) == {2, 3}
                and len(requested[3]) == 1 and requested[2]):
            target = next(iter(requested[3]))
            if max(requested[2]) <= target:
                polynomial_ring = polynomial_ring or self.ring
                square_outputs = {}
                cube_target = polynomial_ring.zero
                for degree in range(target+1):
                    square_coefficient = polynomial_ring.zero
                    for index in range((degree+1)//2):
                        square_coefficient += (
                            2*values[index]*values[degree-index]
                        )
                    if degree % 2 == 0:
                        middle = degree//2
                        square_coefficient += values[middle]*values[middle]
                    if degree in requested[2]:
                        square_outputs[degree] = square_coefficient
                    cube_target += square_coefficient*values[target-degree]
                return {
                    2: {degree: square_outputs[degree]
                        for degree in sorted(requested[2])},
                    3: {target: cube_target},
                }

        polynomial_ring = polynomial_ring or self.ring
        selected = {}
        coefficients = {}
        for exponent in range(1, max_exponent+1):
            degrees = needed.get(exponent, set())
            if exponent == 1:
                coefficients = {degree: values[degree] for degree in degrees}
            elif exponent == 2:
                # The coefficient of a square is symmetric:
                #   [z^d] A(z)^2 = 2*sum_{i<d-i} a_i*a_{d-i}
                #                    + (d even ? a_{d/2}^2 : 0).
                # Avoid computing both products in every off-diagonal pair.
                previous = coefficients
                coefficients = {}
                for degree in sorted(degrees):
                    # Accumulate the symmetric products one at a time.  A
                    # materialized list can retain many enormous sparse
                    # FLINT polynomials simultaneously and cause a peak
                    # allocation failure even though the final coefficient
                    # fits in memory.
                    coefficient = polynomial_ring.zero
                    for index in range((degree+1)//2):
                        coefficient += (
                            2*previous[index]*values[degree-index]
                        )
                    if degree % 2 == 0:
                        middle = degree//2
                        coefficient += previous[middle]*values[middle]
                    coefficients[degree] = coefficient
            else:
                previous = coefficients
                coefficients = {
                    degree: sum(
                        (previous[index]*values[degree-index]
                         for index in range(degree+1)),
                        polynomial_ring.zero,
                    )
                    for degree in sorted(degrees)
                }
            if exponent in requested:
                selected[exponent] = {
                    degree: coefficients[degree]
                    for degree in sorted(requested[exponent])
                }
        return selected

    def series_product_coefficients(self, left, right, degrees, *, polynomial_ring=None):
        """Return selected coefficients of the exact product of two series."""
        polynomial_ring = polynomial_ring or self.ring
        selected = {}
        for degree in sorted(set(degrees)):
            first = max(0, degree-(len(right)-1))
            last = min(degree, len(left)-1)
            selected[degree] = sum(
                (left[i]*right[degree-i] for i in range(first, last+1)),
                polynomial_ring.zero,
            )
        return selected

    def series_power_coefficient(self, values, exponent, distance):
        return self.series_power_coefficients(values, exponent, distance)[distance]

    def new_q_legacy(self, distance, p_current):
        """Solve the exact Q^3-P^5 resonance law through distance fifteen.

        The removed faces occur at distances five, ten, and fifteen:
        -3*kappa*P^3*Q, -3*kappa2*P*Q^2, and -3*kappa3*P^4 in
        Q^3-P^5. Their exact coefficients are retained as free scalars.
        """
        pseries = [*self.p, p_current]
        q_known = [*self.q, self.ring.zero]
        # Retain only the specific power coefficients reused by the coupled
        # block, plus the two convolution prefixes needed for the resonances.
        p_requests = {5: (None, (distance,)), 4: (None, (distance,))}
        q_requests = {3: (None, (distance,)), 2: (None, (distance,))}
        if distance >= 5:
            p_requests[3] = (distance-5, ())
            p_requests[2] = (distance-5, ())
        if distance >= 10:
            q_requests[2] = (distance-10, (distance,))
        if distance >= 15:
            p_requests[4] = (None, (distance, distance-15))
        if self.affine_ring is None:
            # SymPy fallback: retain the original full-ring exact path.
            p_powers = self.selected_series_power_coefficients(
                pseries, p_requests, distance
            )
            q_powers = self.selected_series_power_coefficients(
                q_known, q_requests, distance
            )
        else:
            # All previously solved slices depend only on t and carried
            # parameters.  Compute their power coefficients in that smaller
            # ring; only the new p_current slice contains fresh affine
            # unknowns.  In a coefficient of total distance ``distance``,
            # p_current can occur only once, paired with four (or three)
            # copies of p_0.  Restore those exact linear terms after embedding
            # the parameter-only coefficients.
            affine_ring = self.affine_ring
            p_small = [value.project_to_ring(affine_ring) for value in self.p]
            p_small.append(affine_ring.zero)
            q_small = [value.project_to_ring(affine_ring) for value in self.q]
            q_small.append(affine_ring.zero)
            # The generic recurrence computes P^3 through the full distance
            # to obtain P^5 at only its top requested coefficient.  Use
            # P^2, P^4=(P^2)^2, and P^5=P^4*P instead; P^3 is needed only
            # through distance-5 for the resonance correction.  This avoids
            # constructing the expensive, unused upper P^3 coefficients.
            p2_map = self.selected_series_power_coefficients(
                p_small, {2: (distance, ())}, distance,
                polynomial_ring=affine_ring,
            )[2]
            p2_small = [p2_map[i] for i in range(distance+1)]
            p4_map = self.series_product_coefficients(
                p2_small, p2_small, range(distance+1),
                polynomial_ring=affine_ring,
            )
            p4_small = [p4_map[i] for i in range(distance+1)]
            p5_map = self.series_product_coefficients(
                p4_small, p_small, (distance,), polynomial_ring=affine_ring,
            )
            # The P^3*Q resonance correction at distance d uses every
            # coefficient P^3[i] through i=d-5; retain that full prefix.
            p3_degrees = range(distance-5+1) if distance >= 5 else ()
            p3_map = self.series_product_coefficients(
                p2_small, p_small, p3_degrees, polynomial_ring=affine_ring,
            )
            p2_degrees = (range(distance-5+1) if distance >= 5 else (distance,))
            p4_degrees = (distance, distance-15) if distance >= 15 else (distance,)
            p_powers = {
                2: {degree: p2_map[degree].embed_into(self.ring)
                    for degree in p2_degrees},
                4: {degree: p4_map[degree].embed_into(self.ring)
                    for degree in p4_degrees},
                5: {distance: p5_map[distance].embed_into(self.ring)},
            }
            if distance >= 5:
                p_powers[3] = {
                    degree: p3_map[degree].embed_into(self.ring)
                    for degree in p3_degrees
                }
            del p_small, p2_map, p2_small, p4_map, p4_small, p5_map, p3_map
            # Do not keep both the parameter-ring powers and their full-ring
            # embeddings alive while constructing the second family.
            q_small_powers = self.selected_series_power_coefficients(
                q_small, q_requests, distance, polynomial_ring=affine_ring
            )
            del q_small
            q_powers = {
                exponent: {
                    degree: value.embed_into(self.ring)
                    for degree, value in coefficients.items()
                }
                for exponent, coefficients in q_small_powers.items()
            }
            del q_small_powers
            p_powers[5][distance] += 5*self.p[0]**4*p_current
            p_powers[4][distance] += 4*self.p[0]**3*p_current
        p5 = p_powers[5][distance]
        q3 = q_powers[3][distance]
        correction = self.ring.zero
        if distance >= 5:
            target = distance - 5
            correction = 3*self.kappa*sum(
                (p_powers[3][i]*self.q[target-i] for i in range(target+1)),
                self.ring.zero)
        if distance >= 10:
            target = distance - 10
            correction += 3*self.kappa2*sum(
                (pseries[i]*q_powers[2][target-i] for i in range(target+1)),
                self.ring.zero)
        if distance >= 15:
            target = distance - 15
            correction += 3*self.kappa3*p_powers[4][target]
        numerator = p5 + correction - q3
        denominator = 3*self.u**10
        quotient, remainder = numerator.div(denominator)
        if remainder:
            self.record_univariate_constraints(remainder,
                                                f"corrected-power-distance{distance}")
            print(f"branch={self.tag} corrected_power_distance={distance} "
                  f"remainder_terms={len(remainder)} "
                  f"retained_constraints={len(self.constraints)}", flush=True)
        # p^5[d] and q^3[d] have served their sole purpose in constructing
        # q_current.  The Jacobian block below only consumes p^4[d], p^3 at
        # d-5 (and d-15), p^2 through d-5, and q^2 at its requested degrees.
        # Drop the two largest unused slices and trim the unused lower prefixes
        # before the coupled Jacobian equations are assembled.
        del p_powers[5]
        del q_powers[3]
        p_powers[4] = {distance: p_powers[4][distance]}
        if distance >= 5:
            p3_degrees = {distance-5}
            if distance >= 15:
                p3_degrees.add(distance-15)
            p_powers[3] = {
                degree: p_powers[3][degree]
                for degree in p3_degrees
            }
        return quotient, p_powers, q_powers

    def new_q(self, distance, p_current):
        """Build q_current and only the Jacobian slices actually consumed.

        In the FLINT path, all old-series powers and the two combined
        Jacobian corrections stay in the active carried-parameter ring.  Only
        q_current, P^4[d], and the combined corrections are embedded in the
        full ring containing the current affine unknowns.
        """
        if self.affine_ring is None:
            quotient, p_powers, q_powers = self.new_q_legacy(
                distance, p_current
            )
            p4_full = p_powers[4][distance]
            p_extra = self.ring.zero
            if distance:
                p_extra = -3*q_powers[2][distance]
                if distance >= 5:
                    p_extra += 3*self.kappa*p_powers[3][distance-5]
                if distance >= 10:
                    target = distance-10
                    p_extra += 6*self.kappa2*sum(
                        (self.p[i]*q_powers[2][target-i]
                         for i in range(target+1)), self.ring.zero)
            q_extra = self.ring.zero
            if distance >= 5:
                target = distance-5
                q_extra -= 9*self.kappa*sum(
                    (p_powers[2][i]*self.q[target-i]
                     for i in range(target+1)), self.ring.zero)
            if distance >= 10:
                q_extra -= 3*self.kappa2*q_powers[2][distance-10]
            if distance >= 15:
                q_extra -= 12*self.kappa3*p_powers[3][distance-15]
            return quotient, p4_full, p_extra, q_extra, self.ring.zero

        self._refresh_affine_ring()
        affine_ring = self.affine_ring
        print(f"branch={self.tag} distance={distance} stage=affine-ring-ready "
              f"parameters={self.nparameters} p_slices={len(self.p)} "
              f"q_slices={len(self.q)}", flush=True)
        p_small = [value.project_to_ring(affine_ring) for value in self.p]
        p_small.append(affine_ring.zero)
        q_small = [value.project_to_ring(affine_ring) for value in self.q]
        q_small.append(affine_ring.zero)
        print(f"branch={self.tag} distance={distance} stage=old-series-projected "
              f"p_slices={len(p_small)-1} q_slices={len(q_small)-1}",
              flush=True)

        kappa = self.kappa.project_to_ring(affine_ring)
        kappa2 = self.kappa2.project_to_ring(affine_ring)
        kappa3 = self.kappa3.project_to_ring(affine_ring)

        print(f"branch={self.tag} distance={distance} stage=p2-prefix-begin",
              flush=True)
        p2_map = self.selected_series_power_coefficients(
            p_small, {2: (distance, ())}, distance,
            polynomial_ring=affine_ring,
        )[2]
        print(f"branch={self.tag} distance={distance} stage=p2-prefix-complete "
              f"coefficients={len(p2_map)} terms="
              f"{sum(len(value) for value in p2_map.values())}", flush=True)
        p2_small = [p2_map[i] for i in range(distance+1)]

        # Only P^4 at distance d and, at a resonance, d-15 are reused.  Stream
        # each P^4 coefficient into the single requested P^5 coefficient
        # instead of retaining the full P^4 series alongside P^2, P^3, and
        # Q-power coefficients.  This is algebraically the same convolution,
        # but substantially reduces the peak sparse-polynomial footprint.
        p5 = affine_ring.zero
        p4_old = affine_ring.zero
        p4_resonance = affine_ring.zero
        for degree in range(distance+1):
            # Consume each off-diagonal product immediately.  Retaining the
            # entire list keeps every large sparse FLINT product alive while
            # the sum is accumulated, which needlessly raises peak memory.
            p4_coefficient = sum(
                (2*p2_small[index]*p2_small[degree-index]
                 for index in range((degree+1)//2)),
                affine_ring.zero,
            )
            if degree % 2 == 0:
                middle = degree//2
                p4_coefficient += p2_small[middle]*p2_small[middle]
            p5 += p4_coefficient*p_small[distance-degree]
            if degree == distance:
                p4_old = p4_coefficient
            if distance >= 15 and degree == distance-15:
                p4_resonance = p4_coefficient
            del p4_coefficient
        print(f"branch={self.tag} distance={distance} "
              "stage=p5-convolution-complete", flush=True)

        # The P^3/Q convolution and the two exceptional P^3 coefficients are
        # the only further uses of P^3.  Accumulate the convolution directly
        # and discard each intermediate coefficient immediately.
        p3_q_correction = affine_ring.zero
        p3_at_distance_minus5 = affine_ring.zero
        p3_at_distance_minus15 = affine_ring.zero
        q_extra = affine_ring.zero
        if distance >= 5:
            target = distance-5
            for degree in range(target+1):
                first = max(0, degree-(len(p_small)-1))
                last = min(degree, len(p2_small)-1)
                p3_coefficient = sum(
                    (p2_small[i]*p_small[degree-i]
                     for i in range(first, last+1)), affine_ring.zero,
                )
                p3_q_correction += p3_coefficient*q_small[target-degree]
                if degree == target:
                    p3_at_distance_minus5 = p3_coefficient
                if distance >= 15 and degree == distance-15:
                    p3_at_distance_minus15 = p3_coefficient
                del p3_coefficient
            print(f"branch={self.tag} distance={distance} "
                  "stage=p3-q-convolution-complete", flush=True)
            q_extra -= 9*kappa*sum(
                (p2_map[i]*q_small[target-i] for i in range(target+1)),
                affine_ring.zero,
            )
        del p2_map, p2_small

        # P-series intermediates are now gone; compute Q powers in a separate
        # phase so their large sparse coefficients do not overlap the full
        # P^2/P^3/P^4 prefix.
        q_requests = {3: (None, (distance,)), 2: (None, (distance,))}
        if distance >= 10:
            q_requests[2] = (distance-10, (distance,))
        q_small_powers = self.selected_series_power_coefficients(
            q_small, q_requests, distance, polynomial_ring=affine_ring
        )
        print(f"branch={self.tag} distance={distance} stage=q-powers-complete "
              f"exponents={','.join(str(e) for e in sorted(q_small_powers))}",
              flush=True)

        correction = 3*kappa*p3_q_correction
        if distance >= 10:
            target = distance-10
            correction += 3*kappa2*sum(
                (p_small[i]*q_small_powers[2][target-i]
                 for i in range(target+1)), affine_ring.zero,
            )
        if distance >= 15:
            correction += 3*kappa3*p4_resonance
        numerator = p5+correction-q_small_powers[3][distance]
        denominator = 3*self.u.project_to_ring(affine_ring)**10
        q_old, remainder = numerator.div(denominator)
        if remainder:
            self.record_univariate_constraints(
                remainder, f"corrected-power-distance{distance}"
            )
            print(f"branch={self.tag} corrected_power_distance={distance} "
                  f"remainder_terms={len(remainder)} "
                  f"retained_constraints={len(self.constraints)}", flush=True)
        del remainder

        # The power-compatibility equation has now been divided and its exact
        # remainder retained.  These objects are not used by the Jacobian
        # corrections below; release them before multiplying the large
        # carried-parameter slices that form p_extra and q_extra.
        del numerator, p5, correction, q_small_powers[3]
        del p3_q_correction, p4_resonance
        print(f"branch={self.tag} distance={distance} "
              "stage=power-remainder-released", flush=True)

        p_extra = affine_ring.zero
        if distance:
            p_extra = -3*q_small_powers[2][distance]
            if distance >= 5:
                p_extra += 3*kappa*p3_at_distance_minus5
            if distance >= 10:
                target = distance-10
                p_extra += 6*kappa2*sum(
                    (p_small[i]*q_small_powers[2][target-i]
                     for i in range(target+1)), affine_ring.zero,
                )
        print(f"branch={self.tag} distance={distance} stage=p-extra-built",
              flush=True)

        if distance >= 10:
            q_extra -= 3*kappa2*q_small_powers[2][distance-10]
        if distance >= 15:
            q_extra -= 12*kappa3*p3_at_distance_minus15

        # q_old and p4_old are known-parameter slices.  Keep them in the
        # smaller active ring throughout the streamed affine solve.  Their
        # exact old-only contributions are folded into p_extra/q_extra; q_old
        # is embedded once later, when the solved q coefficient is stored.
        p_extra -= 6*q_small[0]*q_old
        w0_active = self.w[0].project_to_ring(affine_ring)
        self._active_w0_block = w0_active
        self._active_p_block = p_small[:-1]
        self._active_q_block = q_small[:-1]
        q_extra += ((25-distance)*q_old*self.derivative(w0_active)
                    - 39*self.derivative(q_old)*w0_active
                    - 5*p4_old)
        print(f"branch={self.tag} distance={distance} "
              "stage=jacobian-extras-built", flush=True)

        del p_small, q_small, q_small_powers
        del p3_at_distance_minus5, p3_at_distance_minus15
        del denominator, kappa, kappa2, kappa3, p4_old

        q_current = self.const(sp.Rational(5, 3))*self.u**2*p_current
        p4_full = 4*self.p[0]**3*p_current
        # The correction terms contain only already-solved coefficients.
        # Keep them in the active parameter ring; the streamed FLINT affine
        # splitter consumes them there without a full-ring forcing projection.
        return q_current, p4_full, p_extra, q_extra, q_old

    def block(self, distance):
        t = self.t
        print(f"branch={self.tag} distance={distance} stage=begin parameters={self.nparameters} constraints={len(self.constraints)}",
              flush=True)
        w_weight = 39-distance
        w_lo = -(-w_weight//5)
        w_hi = (300+w_weight)//5
        w_count = w_hi-w_lo+1
        if distance >= 5:
            p_degree = 15-distance
            p_low = -((-p_degree)//5)
            p_high = min(15, (75-distance)//5)
            u_count = p_high-p_low+1
        else:
            p_low = None
            u_count = 8 if distance == 2 else 12
        nw = sum((self.unknowns[j]*t**(w_lo+j) for j in range(w_count)), self.ring.zero)
        nu = sum((self.unknowns[w_count+j]*t**j for j in range(u_count)), self.ring.zero)
        print(f"branch={self.tag} distance={distance} stage=unknown-series w_count={w_count} u_count={u_count}",
              flush=True)
        if distance == 5:
            self.kappa = self.allocate(1)[0]
        if distance == 10:
            self.kappa2 = self.allocate(1)[0]
        if distance == 15:
            self.kappa3 = self.allocate(1)[0]
        p_current = (t**p_low*nu if distance >= 5 else
                     t**3*(self.h*nu if distance == 2 else nu))
        q_current, p4_full, p_extra, q_extra, q_current_old = self.new_q(
            distance, p_current
        )
        print(f"branch={self.tag} distance={distance} stage=series-powers p_current={len(p_current)} q_current={len(q_current)}",
              flush=True)
        ps, qs, ws = [*self.p, p_current], [*self.q, q_current], [*self.w, nw]
        flint_affine = getattr(self.ring, "backend_name", None) == "flint-mpoly"

        if flint_affine:
            active_ring = self.affine_ring
            p_old_active = self._active_p_block
            q_old_active = self._active_q_block
            print(f"branch={self.tag} distance={distance} "
                  "stage=old-jacobian-slices-begin", flush=True)
            # For positive distance, only the endpoints i=0 and i=d contain
            # fresh affine unknowns.  Compute the i=1,...,d-1 known slices
            # wholly in the smaller active ring; the endpoint terms are
            # emitted separately below in the full ring.  Project each old
            # w-slice just in time: retaining all d-1 projected slices at
            # once caused the d24 worker to exhaust FLINT's allocation arena
            # before it entered the streamed affine solver.
            p_known = active_ring.zero
            q_known = active_ring.zero
            for i in range(1, distance):
                w_slice = self.w[distance-i].project_to_ring(active_ring)
                w_derivative = self.derivative(w_slice)
                p_known += (
                    (15-i)*p_old_active[i]*w_derivative
                    -(39-(distance-i))*self.derivative(p_old_active[i])
                    *w_slice
                )
                q_known += (
                    (25-i)*q_old_active[i]*w_derivative
                    -(39-(distance-i))*self.derivative(q_old_active[i])
                    *w_slice
                )
                del w_slice, w_derivative
            print(f"branch={self.tag} distance={distance} "
                  f"stage=old-jacobian-slices-complete "
                  f"p_terms={len(p_known)} q_terms={len(q_known)}",
                  flush=True)
            p0_full, q0_full, w0_full = self.p[0], self.q[0], self.w[0]
            del self._active_p_block, self._active_q_block
            del self._active_w0_block, p_old_active, q_old_active

        def affine_term(value, has_fresh_unknown):
            # The FLINT splitter can skip a full projection for old-only
            # forcing summands.  Keep the symbolic backend's input unchanged.
            return ((value, has_fresh_unknown) if flint_affine else value)

        def p_equation_terms():
            if flint_affine:
                if p_known:
                    yield affine_term(p_known, False)
                if 15-distance:
                    yield affine_term(
                        (15-distance)*p_current*self.derivative(w0_full)
                        -39*self.derivative(p_current)*w0_full,
                        True,
                    )
                yield affine_term(
                    15*p0_full*self.derivative(nw)
                    -(39-distance)*self.derivative(p0_full)*nw,
                    True,
                )
                yield affine_term(p_extra, False)
                if distance:
                    yield affine_term(-6*q0_full*q_current, True)
                else:
                    yield affine_term(-3*q_current**2, True)
                return
            for i in range(distance+1):
                j = distance-i
                if 15-i:
                    yield affine_term(
                        (15-i)*ps[i]*self.derivative(ws[j]),
                        i == distance or j == distance,
                    )
                if 39-j:
                    yield affine_term(
                        -(39-j)*self.derivative(ps[i])*ws[j],
                        i == distance or j == distance,
                    )
            # q_known has a zero at the not-yet-solved distance.  Restore
            # the cross terms involving q_current in [t^d]Q^2.
            if distance:
                yield affine_term(p_extra, False)
                yield affine_term(-6*self.q[0]*q_current, True)
            else:
                yield affine_term(-3*q_current**2, True)

        def q_equation_terms():
            if flint_affine:
                if q_known:
                    yield affine_term(q_known, False)
                yield affine_term(
                    (25-distance)*q_current*self.derivative(w0_full)
                    -39*self.derivative(q_current)*w0_full,
                    True,
                )
                yield affine_term(
                    25*q0_full*self.derivative(nw)
                    -(39-distance)*self.derivative(q0_full)*nw,
                    True,
                )
                yield affine_term(-5*p4_full, True)
                yield affine_term(q_extra, False)
                return
            for i in range(distance+1):
                j = distance-i
                if 25-i:
                    yield affine_term(
                        (25-i)*qs[i]*self.derivative(ws[j]),
                        i == distance or j == distance,
                    )
                if 39-j:
                    yield affine_term(
                        -(39-j)*self.derivative(qs[i])*ws[j],
                        i == distance or j == distance,
                    )
            yield affine_term(-5*p4_full, True)
            yield affine_term(q_extra, False)

        print(f"branch={self.tag} distance={distance} stage=streamed-affine-solve-input",
              flush=True)
        solution = self.affine_solve(
            (p_equation_terms(), q_equation_terms()),
            w_count+u_count, distance, streamed=True,
        )
        print(f"branch={self.tag} distance={distance} stage=affine-solve-output variables={w_count+u_count} free={len(solution)}",
              flush=True)
        self.w.append(self.substitute_affine(nw, solution))
        solved_u = self.substitute_affine(nu, solution)
        if distance == 2:
            self.B = solved_u
            self.uk[2] = self.h*solved_u
        else:
            self.uk[distance] = solved_u
        self.p.append(self.substitute_affine(p_current, solution))
        solved_q_current = self.substitute_affine(q_current, solution)
        if q_current_old:
            solved_q_current += q_current_old.embed_into(self.ring)
        self.q.append(solved_q_current)
        self.impose_q_source_support(distance)
        print(f"branch={self.tag} distance={distance} p_terms={len(self.p[-1])} q_terms={len(self.q[-1])} w_terms={len(self.w[-1])} elapsed={time.monotonic()-self.started:.2f}", flush=True)

    def digest(self):
        result = hashlib.sha256()
        for sequence in (self.p, self.q, self.w):
            for polynomial in sequence:
                for monomial, coefficient in sorted(polynomial.items()):
                    result.update(repr((monomial, coefficient)).encode())
        for tag, polynomial in self.constraints:
            result.update(tag.encode())
            polynomial = self.materialize_constraint_polynomial(polynomial)
            for monomial, coefficient in sorted(polynomial.items()):
                result.update(repr((monomial, coefficient)).encode())
        return result.hexdigest()

    def save_checkpoint(self, path):
        """Write exact terms incrementally; avoid a second copy of the state."""
        hasher = hashlib.sha256()
        buffered = []
        buffered_chars = 0

        def write(text):
            nonlocal buffered_chars
            buffered.append(text)
            buffered_chars += len(text)
            if buffered_chars >= (1 << 16):
                flush()

        def flush():
            nonlocal buffered_chars
            if buffered:
                chunk = "".join(buffered)
                handle.write(chunk)
                hasher.update(chunk.encode("utf-8"))
                buffered.clear()
                buffered_chars = 0

        def write_json(value):
            write(json.dumps(value, separators=(",", ":")))

        def write_polynomial(polynomial):
            if isinstance(polynomial, LazyConstraintPolynomial):
                if not polynomial.transformations:
                    remaining = polynomial.end-polynomial.start
                    with polynomial.path.open("rb") as source:
                        source.seek(polynomial.start)
                        while remaining:
                            chunk = source.read(min(1 << 20, remaining))
                            if not chunk:
                                raise OSError("truncated lazy checkpoint constraint")
                            write(chunk.decode("ascii"))
                            remaining -= len(chunk)
                    return
                polynomial = self.materialize_constraint_polynomial(polynomial)
            if polynomial is None:
                write("null")
                return
            write("[")
            first_term = True
            backend = getattr(getattr(polynomial, "ring", None),
                              "backend_name", None)
            sparse_terms = (backend == "flint-mpoly"
                            and hasattr(polynomial, "sparse_items"))
            if sparse_terms:
                terms = polynomial.sparse_items()
            else:
                terms = polynomial.items()
                if backend != "flint-mpoly":
                    terms = sorted(terms)
            # FLINT's indexed term iterator already follows its canonical
            # context monomial order.  Sorting it in Python first materializes
            # every term tuple at once, which is costly for the multi-million
            # term d23 checkpoint.  The serialized term order may differ from
            # the SymPy path, but the polynomial and its exact coefficients do
            # not; the checkpoint hash is computed over the emitted bytes.
            for powers, coefficient in terms:
                if not first_term:
                    write(",")
                first_term = False
                sparse = (powers if sparse_terms else
                          [[i, int(power)] for i, power in enumerate(powers) if power])
                values = [coefficient] if self.field == QQ else coefficient.to_list()
                write("[")
                write_json(sparse)
                write(",")
                write_json([str(value) for value in values])
                write("]")
            write("]")

        with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
            first_field = True

            def field(name, emit):
                nonlocal first_field
                if not first_field:
                    write(",")
                first_field = False
                write_json(name)
                write(":")
                emit()

            def scalar_field(name, value):
                field(name, lambda: write_json(value))

            write("{")
            scalar_field("format", "f2-sparse-ring-v1")
            scalar_field("branch", self.tag)
            scalar_field("quadratic", self.field != QQ)
            scalar_field("polynomial_backend",
                         getattr(self.ring, "backend_name", "sympy"))
            scalar_field("generators", [str(symbol) for symbol in self.ring.symbols])
            scalar_field("nparameters", self.nparameters)
            scalar_field("receipt", self.receipt)

            def write_constraints():
                write("[")
                for index, (tag, polynomial) in enumerate(self.constraints):
                    if index:
                        write(",")
                    write("[")
                    write_json(tag)
                    write(",")
                    write_polynomial(polynomial)
                    write("]")
                write("]")

            field("constraints", write_constraints)
            for name in ("h", "v", "S", "B", "kappa", "kappa2", "kappa3"):
                field(name, lambda name=name: write_polynomial(getattr(self, name)))

            def write_uk():
                write("{")
                for index, (key, polynomial) in enumerate(self.uk.items()):
                    if index:
                        write(",")
                    write_json(str(key))
                    write(":")
                    write_polynomial(polynomial)
                write("}")

            field("uk", write_uk)
            for name in ("p", "q", "w"):
                def write_array(name=name):
                    write("[")
                    for index, polynomial in enumerate(getattr(self, name)):
                        if index:
                            write(",")
                        write_polynomial(polynomial)
                    write("]")
                field(name, write_array)
            if getattr(self, "compaction", None) is not None:
                scalar_field("compaction", self.compaction)
            write("}\n")
            flush()

        self.last_checkpoint_sha256 = hasher.hexdigest()
        byte_count = Path(path).stat().st_size
        print(f"checkpoint={path} bytes={byte_count} sha256={self.last_checkpoint_sha256}", flush=True)


def verify_fifth_q_formula():
    h, s, b, c3, c4, c5, t, w, k = sp.symbols("h s b c3 c4 c5 t w k")
    delta = s/h*w+b/h**2*w**2+c3/h**3*w**3+c4/h**3*w**4+c5/(t*h**3)*w**5
    direct = t**5*sp.expand(h**5*sum(sp.binomial(sp.Rational(5,3),j)*delta**j for j in range(6))).coeff(w,5)+k*(t*h)**4
    expected = (sp.Rational(5,3)*t**4*h**2*c5+t**5*(sp.Rational(10,9)*(h*s*c4+b*c3)-sp.Rational(5,27)*(s**2*c3+s*b**2)+sp.Rational(20,243)*s**3*b-sp.Rational(7,729)*s**5)+k*(t*h)**4)
    assert sp.cancel(direct-expected) == 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch", choices=("2+2", "2+1+1"), default="2+1+1")
    parser.add_argument("--through", type=int, default=5,
                        help="last exact distance to lift (must be at least 2)")
    parser.add_argument("--backend", choices=("sympy", "flint-mpoly"), default="sympy")
    parser.add_argument("--flint-threads", type=int,
                        default=max(1, os.cpu_count() or 1),
                        help="maximum FLINT worker threads for the flint-mpoly backend")
    parser.add_argument("--parameter-reserve", type=int, default=DEFAULT_PARAMETER_RESERVE,
                        help="number of reserved affine parameters, also appended when resuming")
    parser.add_argument("--resume", help="resume from an exact sparse checkpoint")
    parser.add_argument("--legacy-resume", action="store_true",
                        help="use the all-at-once JSON loader for a resumed checkpoint")
    parser.add_argument("--receipt")
    parser.add_argument("--tuned", action="store_true")
    parser.add_argument("--checkpoint")
    parser.add_argument("--terminal", action="store_true")
    parser.add_argument("--terminal-q", action="store_true")
    parser.add_argument("--root-values", action="store_true")
    parser.add_argument("--deep-rootwise-dd", action="store_true",
                        help="impose the exact 2+1+1 DD divisors on U2--U4")
    parser.add_argument("--skip-polynomial-digest", action="store_true",
                        help="skip the slow in-memory term digest; report the exact checkpoint SHA-256 instead")
    args = parser.parse_args()
    if args.through < 2:
        parser.error("--through must be at least 2")
    if args.flint_threads < 1:
        parser.error("--flint-threads must be at least 1")
    if args.deep_rootwise_dd and (args.branch != "2+1+1" or
                                  not args.terminal or not args.root_values):
        raise ValueError("--deep-rootwise-dd requires 2+1+1, --terminal, and --root-values")
    global POLY_BACKEND
    POLY_BACKEND = args.backend
    if args.backend == "flint-mpoly":
        import flint
        flint.ctx.threads = args.flint_threads
        print(f"flint_threads={flint.ctx.threads}", flush=True)
    verify_fifth_q_formula()
    t = sp.symbols("t")
    if args.branch == "2+1+1":
        r = t**2-3*t+3
        h, v, field = (t-1)**2*r, -sp.Rational(1,3)*t*(t-1)*r, QQ
    else:
        alpha = (3+sp.sqrt(5))/2
        h = (t-1)**2*(t-alpha)**2
        v = t*(t-1)*(t-alpha)*((4*alpha-11)/3*t+3-alpha)
        field = QQ.algebraic_field(sp.sqrt(5))
    resume_loader = (SparseLift.from_checkpoint if args.legacy_resume
                     else SparseLift.from_checkpoint_streaming)
    lift = (resume_loader(args.resume,
                          parameter_reserve=args.parameter_reserve)
            if args.resume else
            SparseLift(args.branch, sp.expand(h), sp.expand(v), field,
                       parameter_reserve=args.parameter_reserve))
    if args.resume and lift.tag != args.branch:
        raise ValueError("--branch disagrees with checkpoint")
    checkpoint_saved = False
    if args.terminal and not args.resume:
        if args.tuned:
            raise ValueError("terminal and tuned are separate branch relaxations")
        lift.initialize_terminal_chart()
    for distance in range(len(lift.p), args.through+1):
        lift.block(distance)
        if args.tuned:
            lift.impose_tuned_contact(distance)
        if args.terminal:
            lift.impose_terminal_chart(distance)
        if args.terminal_q:
            lift.impose_terminal_q(distance)
        if args.root_values:
            lift.impose_root_values(distance)
        if args.deep_rootwise_dd:
            lift.impose_deep_rootwise_dd(distance)
        if args.checkpoint:
            lift.save_checkpoint(args.checkpoint)
            checkpoint_saved = True
    report = dict(branch=args.branch, field=str(field), through=args.through,
                  p10_support="full t^2 times degree<=12", blocks=lift.receipt,
                  digest=None if args.skip_polynomial_digest else lift.digest(),
                  digest_kind=("omitted; use checkpoint_sha256" if args.skip_polynomial_digest
                               else "canonical sparse-polynomial term stream"),
                  checkpoint_sha256=getattr(lift, "last_checkpoint_sha256", None),
                  tuned_contact=args.tuned, terminal_chain=args.terminal,
                  terminal_q=args.terminal_q,
                  root_values=args.root_values,
                  deep_rootwise_dd=args.deep_rootwise_dd,
                  retained_constraint_count=len(lift.constraints),
                  retained_constraint_tags=[tag for tag, _ in lift.constraints],
                  scope="support-relaxed necessary system; not exclusion",
                  elapsed_seconds=time.monotonic()-lift.started)
    if args.receipt:
        Path(args.receipt).write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    if args.checkpoint and not checkpoint_saved:
        lift.save_checkpoint(args.checkpoint)
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
