"""Exact QQ sparse-polynomial adapter backed by FLINT's fmpq_mpoly.

This backend is an optional acceleration layer for the rational-coefficient
F2 recurrence. It deliberately rejects number fields; the Q(sqrt(5)) branches
continue to use SymPy until a separate exact backend is available.
"""

from __future__ import annotations

from typing import Any

import sympy as sp
from flint import fmpq, fmpq_mpoly, fmpq_mpoly_ctx


def _as_fmpq(value: Any) -> fmpq:
    if isinstance(value, fmpq):
        return value
    if isinstance(value, sp.Rational):
        return fmpq(int(value.p), int(value.q))
    try:
        return fmpq(value)
    except (TypeError, ValueError):
        return fmpq(str(value))


class FlintRing:
    """A small duck-typed ring interface used by ``SparseLift`` over QQ."""

    backend_name = "flint-mpoly"

    def __init__(self, names, field=None):
        if field is not None and str(field) != "QQ":
            raise TypeError("FLINT mpoly backend supports only QQ")
        self.symbols = tuple(str(name) for name in names)
        self.ctx = fmpq_mpoly_ctx.get(self.symbols, "lex")
        self.gens = tuple(FGen(self, i) for i in range(len(self.symbols)))
        self._gen_data = self.ctx.gens()
        self._zero_data = fmpq_mpoly(0, ctx=self.ctx)
        self._one_data = fmpq_mpoly(1, ctx=self.ctx)
        self.zero = FPoly(self, self._zero_data)
        self.one = FPoly(self, self._one_data)
        self._sympy_symbols = tuple(sp.Symbol(name) for name in self.symbols)
        self._symbol_index = {name: i for i, name in enumerate(self.symbols)}

    def __iter__(self):
        # Match the return shape of sympy.polys.rings.ring(...).
        yield self
        yield from self.gens

    def _coerce(self, value):
        if isinstance(value, FPoly):
            if value.ring is not self:
                raise TypeError("polynomials belong to different FLINT rings")
            return value.data
        if isinstance(value, FGen):
            if value.ring is not self:
                raise TypeError("generators belong to different FLINT rings")
            return self._gen_data[value.index]
        if isinstance(value, fmpq_mpoly):
            if value.context() != self.ctx:
                raise TypeError("FLINT polynomials belong to different contexts")
            return value
        return fmpq_mpoly(_as_fmpq(value), ctx=self.ctx)

    def from_dict(self, terms):
        clean = {}
        shift = 0
        if terms:
            shift = min(0, min(int(powers[0]) for powers in terms))
        for powers, coefficient in terms.items():
            powers = tuple(int(power) for power in powers)
            if len(powers) != len(self.symbols):
                raise ValueError("monomial has the wrong number of exponents")
            coefficient = _as_fmpq(coefficient)
            if coefficient:
                if shift:
                    powers = (powers[0]-shift, *powers[1:])
                clean[powers] = coefficient
        return FPoly(self, self.ctx.from_dict(clean), shift)

    def from_prevalidated_dict(self, terms):
        """Build from decoder-produced dense exponent tuples without recopying.

        The streaming checkpoint reader already guarantees exact tuple arity,
        integer exponents, and FLINT-compatible rational coefficients. Running
        those entries through ``from_dict`` again duplicated a Python pass over
        every exponent of every saved term, which dominated large checkpoint
        resumes.
        """
        if terms and min(int(powers[0]) for powers in terms) < 0:
            return self.from_dict(terms)
        return FPoly(self, self.ctx.from_dict(terms))

    def from_native_dict(self, terms):
        """Build directly from FLINT's native exponent tuples.

        ``fmpq_mpoly.monomial`` returns tuples of ``fmpz``.  FLINT's own
        ``from_dict`` accepts those entries, so callers that already have
        native tuples can avoid converting every exponent to a Python int.
        """
        if terms and min(int(powers[0]) for powers in terms) < 0:
            return self.from_dict(terms)
        return FPoly(self, self.ctx.from_dict(terms))

    def from_nonnegative_native_dict(self, terms):
        """Build from native terms when the caller already proves nonnegativity.

        This avoids rescanning large dictionaries just to check an invariant
        guaranteed by FLINT's polynomial exponent tuples and by the caller's
        construction.  Use only for terms whose exponents are all known to be
        nonnegative.
        """
        return FPoly(self, self.ctx.from_dict(terms))

    def as_poly(self, value):
        if isinstance(value, FPoly):
            if value.ring is not self:
                raise TypeError("polynomials belong to different FLINT rings")
            return value
        if isinstance(value, FGen):
            if value.ring is not self:
                raise TypeError("generators belong to different FLINT rings")
            return value._poly()
        return FPoly(self, self._coerce(value))

    def prefix(self, length):
        """Return the lexicographic subring on the first ``length`` names."""
        return FlintRing(self.symbols[:length])

    def from_expr(self, expression):
        expression = sp.expand(expression)
        if not expression:
            return self.zero
        polynomial = sp.Poly(expression, *self._sympy_symbols, domain=sp.QQ)
        return self.from_dict(dict(polynomial.terms()))

    def ground_new(self, coefficient):
        return FPoly(self, self._coerce(coefficient))

    def term_new(self, powers, coefficient):
        powers = tuple(int(power) for power in powers)
        return self.from_dict({powers: coefficient})


class FGen:
    __slots__ = ("ring", "index")

    def __init__(self, ring, index):
        self.ring = ring
        self.index = int(index)

    def _poly(self):
        return FPoly(self.ring, self.ring._gen_data[self.index])

    def __hash__(self):
        return hash((id(self.ring), self.index))

    def __eq__(self, other):
        return isinstance(other, FGen) and self.ring is other.ring and self.index == other.index

    def __str__(self):
        return self.ring.symbols[self.index]

    def __repr__(self):
        return str(self)

    def items(self):
        """Expose a generator as the one-term sparse polynomial it is."""
        powers = [0] * len(self.ring.symbols)
        powers[self.index] = 1
        return iter(((tuple(powers), 1),))

    def __add__(self, other):
        return self._poly() + other

    def __radd__(self, other):
        return other + self._poly()

    def __sub__(self, other):
        return self._poly() - other

    def __rsub__(self, other):
        return other - self._poly()

    def __mul__(self, other):
        return self._poly() * other

    def __rmul__(self, other):
        return other * self._poly()

    def __truediv__(self, other):
        return self._poly() / other

    def __pow__(self, exponent):
        exponent = int(exponent)
        if exponent < 0:
            if self.index != 0:
                raise ValueError("negative powers are supported only for t")
            return FPoly(self.ring, self.ring._one_data, exponent)
        return self._poly() ** exponent

    def __neg__(self):
        return -self._poly()

    def compose(self, replacements):
        return self._poly().compose(replacements)


class FPoly:
    __slots__ = ("ring", "data", "shift")

    def __init__(self, ring, data, shift=0):
        self.ring = ring
        self.data = data
        self.shift = int(shift)

    def __len__(self):
        # FLINT exposes the term count directly.  Converting a large
        # multivariate polynomial to a Python dictionary just to count its
        # terms can multiply peak memory in the high-dimensional F2 lift.
        return len(self.data)

    def __iter__(self):
        # python-flint's monoms() materializes every exponent tuple.  Walk
        # FLINT's indexed term storage instead so a large sparse polynomial
        # does not require a second full support list in Python memory.
        data = self.data
        shift = self.shift
        return ((int(data.monomial(i)[0])+shift,
                 *(int(power) for power in data.monomial(i)[1:]))
                for i in range(len(data)))

    def __getitem__(self, monomial):
        powers = tuple(int(power) for power in monomial)
        powers = (powers[0]-self.shift, *powers[1:])
        if powers[0] < 0:
            return 0
        return self.data[powers]

    def __bool__(self):
        return not self.data.is_zero()

    def __repr__(self):
        if self.shift:
            return f"t^{self.shift}*({self.data!r})"
        return repr(self.data)

    def __eq__(self, other):
        try:
            other = self.ring.as_poly(other)
            if self.shift == other.shift:
                return self.data == other.data
            base = min(self.shift, other.shift)
            t = self.ring._gen_data[0]
            left = self.data * (t ** (self.shift-base))
            right = other.data * (t ** (other.shift-base))
            return left == right
        except (TypeError, ValueError):
            return False

    def __add__(self, other):
        other = self.ring.as_poly(other)
        # The sparse series recurrences add very large polynomials with the
        # same Laurent shift.  In that common case, multiplying each operand
        # by t**0 only allocates/copies FLINT polynomials before adding them.
        # Preserve the shared shift and add the native polynomials directly.
        if self.shift == other.shift:
            return FPoly(self.ring, self.data + other.data, self.shift)
        base = min(self.shift, other.shift)
        t = self.ring._gen_data[0]
        left = (self.data if self.shift == base else
                self.data * (t ** (self.shift-base)))
        right = (other.data if other.shift == base else
                 other.data * (t ** (other.shift-base)))
        return FPoly(self.ring, left + right, base)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        other = self.ring.as_poly(other)
        # As in __add__, keep a common Laurent shift outside the native
        # polynomial.  Multiplying both operands by t**0 needlessly copies
        # very large FLINT polynomials in the streamed Jacobian slices.
        if self.shift == other.shift:
            return FPoly(self.ring, self.data - other.data, self.shift)
        base = min(self.shift, other.shift)
        t = self.ring._gen_data[0]
        left = (self.data if self.shift == base else
                self.data * (t ** (self.shift-base)))
        right = (other.data if other.shift == base else
                 other.data * (t ** (other.shift-base)))
        return FPoly(self.ring, left - right, base)

    def __rsub__(self, other):
        return self.ring.as_poly(other) - self

    def __mul__(self, other):
        other = self.ring.as_poly(other)
        return FPoly(self.ring, self.data * other.data,
                     self.shift + other.shift)

    def scale(self, scalar):
        """Scale coefficients directly, without making a constant mpoly."""
        return FPoly(self.ring, self.data * _as_fmpq(scalar), self.shift)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        quotient, remainder = self.div(other)
        if remainder:
            raise ArithmeticError("polynomial division is not exact")
        return quotient

    def __pow__(self, exponent):
        exponent = int(exponent)
        if exponent < 0:
            if self.data != self.ring._one_data:
                raise ValueError("negative powers are supported only for a t monomial")
            return FPoly(self.ring, self.ring._one_data,
                         self.shift * exponent)
        return FPoly(self.ring, self.data ** exponent, self.shift * exponent)

    def __neg__(self):
        return FPoly(self.ring, -self.data, self.shift)

    def items(self):
        # fmpq_mpoly.terms() looks lazy but currently zips two fully
        # materialized lists from monoms() and coeffs().  Iterate indexed
        # FLINT storage so large high-dimensional coefficient polynomials do
        # not create those lists in Python memory.
        data = self.data
        shift = self.shift
        return (((int(data.monomial(i)[0])+shift,
                  *(int(power) for power in data.monomial(i)[1:])),
                 data.coefficient(i)) for i in range(len(data)))

    def native_items(self):
        """Iterate indexed terms with one native monomial read per term.

        Calling ``monomial(i)`` twice rebuilt the full exponent tuple twice.
        The distance-23 affine extraction scans hundreds of thousands of
        terms in a ring with hundreds of generators, so keep the single
        native tuple alive while reading its coefficient.
        """
        data = self.data
        shift = self.shift

        def iter_terms():
            for i in range(len(data)):
                powers = data.monomial(i)
                yield ((int(powers[0])+shift, *powers[1:]),
                       data.coefficient(i))

        return iter_terms()

    def sparse_items(self):
        """Iterate sparse exponent pairs without scanning unused variables.

        FLINT stores each monomial as a dense exponent tuple.  Most checkpoint
        polynomials use only a small subset of the several hundred generators,
        so obtain that subset once from ``degrees()`` and inspect only those
        coordinates for each term.  The Laurent shift applies to generator 0.
        """
        data = self.data
        shift = self.shift
        active = [index for index, degree in enumerate(data.degrees()) if degree]
        if shift and 0 not in active:
            active.insert(0, 0)

        def iter_terms():
            for term_index in range(len(data)):
                powers = data.monomial(term_index)
                sparse = []
                for index in active:
                    power = int(powers[index]) + (shift if index == 0 else 0)
                    if power:
                        sparse.append((index, power))
                yield sparse, data.coefficient(term_index)

        return iter_terms()

    def first_variable_coefficients(self):
        """Yield exact coefficients in the first generator by a single term pass.

        The FLINT affine solver needs the coefficient polynomials of the first
        (Laurent) variable.  Repeated differentiation and evaluation scans the
        full sparse polynomial once for each exponent.  Under the F2 lift this
        can mean hundreds of scans over multi-million-term polynomials.  FLINT
        stores lexicographic terms grouped by the first exponent, so collect
        just one coefficient slice at a time and construct it in the same
        context.  The transient dictionary is bounded by the largest single
        slice rather than the whole polynomial.
        """
        data = self.data
        term_count = len(data)
        if not term_count:
            return

        shift = self.shift
        first = int(data.monomial(0)[0])
        last = int(data.monomial(term_count - 1)[0])
        if first == last:
            coefficient = data.subs({0: fmpq(1)})
            if not coefficient.is_zero():
                yield first + shift, FPoly(self.ring, coefficient)
            return

        current_degree = None
        current_terms = {}

        def emit_slice(degree, terms):
            if degree is None or not terms:
                return None
            return self.ring.from_nonnegative_native_dict(terms)

        for term_index in range(term_count):
            powers = data.monomial(term_index)
            degree = int(powers[0])
            if current_degree is None:
                current_degree = degree
            elif degree != current_degree:
                coefficient = emit_slice(current_degree, current_terms)
                if coefficient:
                    yield current_degree + shift, coefficient
                current_terms = {}
                current_degree = degree
            # FLINT's canonical support has no duplicate monomials; within a
            # fixed first exponent, dropping that exponent remains injective.
            current_terms[(0, *powers[1:])] = data.coefficient(term_index)

        coefficient = emit_slice(current_degree, current_terms)
        if coefficient:
            yield current_degree + shift, coefficient

    def project_to_ring(self, target_ring):
        """Set generators absent from ``target_ring`` to zero in FLINT."""
        return FPoly(target_ring,
                     self.data.project_to_context(target_ring.ctx), self.shift)

    def embed_into(self, target_ring):
        """Embed this polynomial into a larger ring with matching names."""
        # `project_to_context` also handles a same-name superset context: it
        # maps each source generator to its matching target generator and
        # pads the additional target exponents with zero.  This is an exact
        # embedding, while generic `compose` needlessly multiplies each large
        # sparse polynomial through every source generator substitution.
        for name in self.ring.symbols:
            if name not in target_ring._symbol_index:
                raise KeyError(name)
        return FPoly(
            target_ring,
            self.data.project_to_context(target_ring.ctx),
            self.shift,
        )

    def diff(self, variable):
        if isinstance(variable, FGen):
            index = variable.index
        elif isinstance(variable, str):
            index = self.ring._symbol_index[variable]
        else:
            index = int(variable)
        result = FPoly(self.ring, self.data.derivative(index), self.shift)
        if index == 0 and self.shift:
            result += self * (self.shift) / self.ring.gens[0]
        return result

    def _strip_t(self):
        """Factor the common t-power into the Laurent shift."""
        if self.data.is_zero():
            return self.shift, self.data
        # The context uses lex order with t as its first variable, so native
        # FLINT terms are grouped monotonically by t-degree.  The minimum is
        # therefore at one endpoint; scanning every term in a multi-million
        # term numerator here was an avoidable Python-level bottleneck.
        term_count = len(self.data)
        minimum = min(
            int(self.data.monomial(0)[0]),
            int(self.data.monomial(term_count-1)[0]),
        )
        if minimum == 0:
            return self.shift, self.data
        t_power = self.ring._gen_data[0] ** minimum
        quotient, remainder = divmod(self.data, t_power)
        if not remainder.is_zero():
            raise AssertionError("computed t-valuation did not divide the FLINT polynomial")
        return self.shift + minimum, quotient

    def div(self, divisor):
        divisor = self.ring.as_poly(divisor)
        dividend_shift, dividend_data = self._strip_t()
        divisor_shift, divisor_data = divisor._strip_t()
        quotient, remainder = divmod(dividend_data, divisor_data)
        return (FPoly(self.ring, quotient, dividend_shift-divisor_shift),
                FPoly(self.ring, remainder, dividend_shift))

    def rem(self, divisor):
        divisor = self.ring.as_poly(divisor)
        divisor_shift, divisor_data = divisor._strip_t()
        terms = divisor_data.to_dict()
        # The terminal slices repeatedly reduce modulo t^k. Avoid invoking a
        # general multivariate division for this common monomial modulus.
        if len(terms) == 1:
            powers, coefficient = next(iter(terms.items()))
            degree = int(powers[0]) + divisor_shift
            if coefficient == 1 and degree > 0 and not any(powers[1:]):
                kept = {m: c for m, c in self.items() if m[0] < degree}
                return self.ring.from_dict(kept)
        _, remainder = self.div(divisor)
        return remainder

    def exquo(self, divisor):
        return self / divisor

    def gcd(self, other):
        other = self.ring.as_poly(other)
        left_shift, left_data = self._strip_t()
        right_shift, right_data = other._strip_t()
        return FPoly(self.ring, left_data.gcd(right_data),
                     min(left_shift, right_shift))

    def compose(self, replacements):
        args = list(self.ring._gen_data)
        t_was_replaced = False
        for key, value in replacements.items():
            if isinstance(key, FGen):
                index = key.index
            elif isinstance(key, str):
                index = self.ring._symbol_index[key]
            else:
                raise TypeError(f"unsupported FLINT substitution key: {key!r}")
            args[index] = self.ring._coerce(value)
            t_was_replaced |= index == 0
        result = FPoly(self.ring, self.data.compose(*args))
        if self.shift and not t_was_replaced:
            return FPoly(self.ring, result.data, self.shift)
        if self.shift:
            return result * self.ring.as_poly(args[0]) ** self.shift
        return result

