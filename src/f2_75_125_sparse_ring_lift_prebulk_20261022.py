"""Exact sparse-parameter replay of the coupled F2 edge bands.

All coefficient arithmetic stays in K[t,p_0,...,p_159], where K is QQ or
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
from math import comb, factorial
from pathlib import Path
import time

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring as sympy_ring

from f2_75_125_power_relation_pq_nextbands import weight49_parameters


POLY_BACKEND = "sympy"


def ring(names, field):
    """Select the exact sparse polynomial backend for a QQ recurrence."""
    if POLY_BACKEND == "flint-mpoly":
        if field != QQ:
            raise TypeError("flint-mpoly backend currently supports only QQ")
        from f2_75_125_flint_mpoly_backend import FlintRing
        return FlintRing(names, field)
    return sympy_ring(names, field)


class SparseLift:
    def __init__(self, tag, h_expr, v_expr, field):
        self.tag, self.field = tag, field
        names = ["t", *(f"p{i}" for i in range(160)), *(f"x{i}" for i in range(73))]
        self.ring, *generators = ring(names, field)
        self.t = generators[0]
        self.parameters = generators[1:161]
        self.unknowns = generators[161:]
        self.new_start = 1 + len(self.parameters)
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
        self.uk = {1: self.h**2*self.S}
        self.p = [self.u**3, self.t**3*self.uk[1]]
        self.q = [self.u**5, self.const(sp.Rational(5, 3))*self.t**5*self.h**2*self.uk[1]]
        self.B = None
        self.kappa = self.ring.zero
        self.kappa2 = self.ring.zero
        self.kappa3 = self.ring.zero

    @classmethod
    def from_checkpoint(cls, path, *, verify_saved_terminal_rows=True):
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
        obj.ring, *generators = ring(old_names, field)
        new_indices = {name: i for i, name in enumerate(old_names)}
        if "t" not in new_indices:
            raise ValueError("checkpoint has no t generator")
        obj.t = generators[new_indices["t"]]
        parameter_names = sorted(
            (name for name in old_names
             if name.startswith("p") and name[1:].isdigit()),
            key=lambda name: int(name[1:]),
        )
        expected_parameters = [f"p{i}" for i in range(len(parameter_names))]
        if parameter_names != expected_parameters:
            raise ValueError("checkpoint parameter generators must be contiguous from p0")
        unknown_names = sorted(
            (name for name in old_names
             if name.startswith("x") and name[1:].isdigit()),
            key=lambda name: int(name[1:]),
        )
        if unknown_names != [f"x{i}" for i in range(73)]:
            raise ValueError("checkpoint must retain unknown generators x0 through x72")
        obj.parameters = [generators[new_indices[name]] for name in parameter_names]
        obj.unknowns = [generators[new_indices[name]] for name in unknown_names]
        obj.new_start = 1 + len(obj.parameters)
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

    def const(self, value):
        return self.ring.ground_new(self.field.convert(value))

    def allocate(self, count):
        begin = self.nparameters
        self.nparameters += count
        if self.nparameters > len(self.parameters):
            raise ValueError("parameter reserve exhausted")
        return self.parameters[begin:self.nparameters]

    def derivative(self, value):
        return value.diff(self.t)

    def split_affine(self, expression):
        """Extract all rows, including forcing-only rows, directly by monomial."""
        rows, forcing = {}, {}
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
                forcing[degree] = forcing.get(degree, self.ring.zero) + self.ring.term_new(monomial, coefficient)
        degrees = sorted(set(rows) | set(forcing))
        return degrees, [rows.get(d, [self.field.zero]*len(self.unknowns)) for d in degrees], [forcing.get(d, self.ring.zero) for d in degrees]

    def record_univariate_constraints(self, expression, tag):
        """Retain every t-coefficient of a polynomial compatibility equation."""
        for degree in sorted({powers[0] for powers in expression}):
            value = self.local_coefficient(expression, degree)
            if value:
                self.constraints.append((f"{tag}:t^{degree}", value))

    def affine_solve(self, equations, nvariables, step, reuse_free=None):
        rows, forcing, row_tags = [], [], []
        for side, expression in enumerate(equations):
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
            value = forcing[i]
            for j in free_cols:
                if rows[i][j]:
                    value += solution[j]*rows[i][j]
            pivot_forcing.append(value)
        for j, inverse_row in zip(pivot_cols, inverse):
            value = self.ring.zero
            for coefficient, rhs in zip(inverse_row, pivot_forcing):
                if coefficient and rhs:
                    value -= rhs*coefficient
            solution[j] = value
        residuals = []
        for tag, row, rhs in zip(row_tags, rows, forcing):
            value = rhs
            for coefficient, entry in zip(row, solution):
                if coefficient and entry:
                    value += entry*coefficient
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
        self.constraints = [(tag, value.compose(replacements))
                            for tag, value in self.constraints]
        for name in ("center_a", "center_b"):
            if hasattr(self, name):
                setattr(self, name, getattr(self, name).compose(replacements))

    def impose_congruence(self, expression, modulus, targets, tag):
        """Solve a necessary congruence in a selected previous affine fibre."""
        remainder = expression.rem(modulus)
        renaming = dict(zip(targets, self.unknowns))
        equation = remainder.compose(renaming)
        solution = self.affine_solve((equation,), len(targets), tag, reuse_free=targets)
        replacements = dict(zip(targets, solution))
        self.replace_parameters(replacements)
        remainder_after = expression.compose(replacements).rem(modulus)
        if remainder_after:
            self.record_univariate_constraints(remainder_after, f"{tag}:postsolve")

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

    def sixth_numerator(self):
        d = 3*self.B-self.S**2
        delta = 27*self.uk[3]+2*self.S**3-9*self.S*self.B
        return delta**2-3*d**3

    def local_coefficient(self, polynomial, degree):
        return self.ring.from_dict({(0, *powers[1:]): coefficient
                                    for powers, coefficient in polynomial.items()
                                    if powers[0] == degree})

    def local_slice(self, distance, bound):
        """Taylor-expand X^-15 P(X^5,X^-1(1+u)) in u=t.

        A global band monomial t^i becomes (1+u)^(5i-15+s).
        The exponent is its original nonnegative y exponent, so this
        transformation uses the actual source monomials and no Laurent cut.
        """
        out = self.ring.zero
        for powers, coefficient in self.p[distance].items():
            exponent = 5*powers[0]-15+distance
            if exponent < 0:
                raise AssertionError("P monomial is outside original polynomial support")
            for degree in range(min(exponent, bound)+1):
                monomial = (degree, *powers[1:])
                out += self.ring.term_new(monomial, coefficient*comb(exponent, degree))
        return out

    def local_q_slice(self, distance, bound):
        """Taylor-expand X^-25 Q(X^5,X^-1(1+u)) in the local variable u."""
        out = self.ring.zero
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
                out += self.ring.term_new(monomial, coefficient*comb(exponent, degree))
        return out

    def local_q_slice_derivative(self, distance, bound, order):
        """Return ``d^order/du^order local_q_slice`` through ``u^bound``.

        Construct the differentiated Taylor coefficients directly, rather
        than expanding a larger slice and repeatedly differentiating a large
        sparse multivariate polynomial.
        """
        out = self.ring.zero
        source_bound = bound+order
        low = -(-(25-distance)//5)
        high = min(25, (125-distance)//5)
        for powers, coefficient in self.q[distance].items():
            if powers[0] < low or powers[0] > high:
                continue
            exponent = 5*powers[0]-25+distance
            if exponent < 0:
                raise AssertionError("Q monomial is outside original polynomial support")
            for degree in range(order, min(exponent, source_bound)+1):
                output_degree = degree-order
                factor = comb(exponent, degree)*factorial(degree)//factorial(output_degree)
                monomial = (output_degree, *powers[1:])
                out += self.ring.term_new(monomial, coefficient*factor)
        return out

    def impose_q_source_support(self, distance):
        """Retain every Q band coefficient outside the original rectangle."""
        low = -(-(25-distance)//5)
        high = min(25, (125-distance)//5)
        bad_degrees = sorted({powers[0] for powers in self.q[distance]
                              if powers[0] < low or powers[0] > high})
        before = len(self.constraints)
        for degree in bad_degrees:
            value = self.local_coefficient(self.q[distance], degree)
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
        out = self.ring.zero
        for powers, coefficient in expression.items():
            new_powers = powers[self.new_start:]
            if any(new_powers):
                if sum(new_powers) != 1:
                    raise AssertionError("substitution is not affine")
                j = new_powers.index(1)
                monomial = (*powers[:self.new_start], *((0,)*len(self.unknowns)))
                out += self.ring.term_new(monomial, coefficient)*solution[j]
            else:
                out += self.ring.term_new(powers, coefficient)
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

    def series_power_coefficient(self, values, exponent, distance):
        return self.series_power_coefficients(values, exponent, distance)[distance]

    def new_q(self, distance, p_current):
        """Solve the exact Q^3-P^5 resonance law through distance fifteen.

        The removed faces occur at distances five, ten, and fifteen:
        -3*kappa*P^3*Q, -3*kappa2*P*Q^2, and -3*kappa3*P^4 in
        Q^3-P^5. Their exact coefficients are retained as free scalars.
        """
        pseries = [*self.p, p_current]
        q_known = [*self.q, self.ring.zero]
        p5 = self.series_power_coefficient(pseries, 5, distance)
        q3 = self.series_power_coefficient(q_known, 3, distance)
        correction = self.ring.zero
        if distance >= 5:
            target = distance - 5
            p3 = self.series_power_coefficients(pseries, 3, target)
            correction = 3*self.kappa*sum(
                (p3[i]*self.q[target-i] for i in range(target+1)),
                self.ring.zero)
        if distance >= 10:
            target = distance - 10
            q2 = self.series_power_coefficients(q_known, 2, target)
            correction += 3*self.kappa2*sum(
                (pseries[i]*q2[target-i] for i in range(target+1)),
                self.ring.zero)
        if distance >= 15:
            target = distance - 15
            p4 = self.series_power_coefficients(pseries, 4, target)
            correction += 3*self.kappa3*p4[target]
        numerator = p5 + correction - q3
        denominator = 3*self.u**10
        quotient, remainder = numerator.div(denominator)
        if remainder:
            self.record_univariate_constraints(remainder,
                                                f"corrected-power-distance{distance}")
            print(f"branch={self.tag} corrected_power_distance={distance} "
                  f"remainder_terms={len(remainder)} "
                  f"retained_constraints={len(self.constraints)}", flush=True)
        return quotient

    def block(self, distance):
        t = self.t
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
        if distance == 5:
            self.kappa = self.allocate(1)[0]
        if distance == 10:
            self.kappa2 = self.allocate(1)[0]
        if distance == 15:
            self.kappa3 = self.allocate(1)[0]
        p_current = (t**p_low*nu if distance >= 5 else
                     t**3*(self.h*nu if distance == 2 else nu))
        q_current = self.new_q(distance, p_current)
        ps, qs, ws = [*self.p, p_current], [*self.q, q_current], [*self.w, nw]
        p_equation = self.ring.zero
        q_equation = self.ring.zero
        for i in range(distance+1):
            j = distance-i
            p_equation += (15-i)*ps[i]*self.derivative(ws[j])-(39-j)*self.derivative(ps[i])*ws[j]
            q_equation += (25-i)*qs[i]*self.derivative(ws[j])-(39-j)*self.derivative(qs[i])*ws[j]
        p_equation -= 3*self.series_power_coefficient(qs, 2, distance)
        q_equation -= 5*self.series_power_coefficient(ps, 4, distance)
        if distance >= 5:
            correction_distance = distance - 5
            p_correction = self.series_power_coefficient(ps, 3, correction_distance)
            p2 = self.series_power_coefficients(ps, 2, correction_distance)
            q_correction = sum(
                (p2[i]*qs[correction_distance-i]
                 for i in range(correction_distance + 1)), self.ring.zero)
            p_equation += 3*self.kappa*p_correction
            q_equation -= 9*self.kappa*q_correction
        if distance >= 10:
            correction_distance = distance - 10
            pq_correction = sum(
                (ps[i]*self.series_power_coefficient(qs, 2,
                                                      correction_distance-i)
                 for i in range(correction_distance+1)), self.ring.zero)
            q2_correction = self.series_power_coefficient(
                qs, 2, correction_distance)
            p_equation += 6*self.kappa2*pq_correction
            q_equation -= 3*self.kappa2*q2_correction
        if distance >= 15:
            correction_distance = distance - 15
            q_equation -= 12*self.kappa3*self.series_power_coefficient(
                ps, 3, correction_distance
            )
        solution = self.affine_solve((p_equation, q_equation), w_count+u_count, distance)
        self.w.append(self.substitute_affine(nw, solution))
        solved_u = self.substitute_affine(nu, solution)
        if distance == 2:
            self.B = solved_u
            self.uk[2] = self.h*solved_u
        else:
            self.uk[distance] = solved_u
        self.p.append(self.substitute_affine(p_current, solution))
        self.q.append(self.substitute_affine(q_current, solution))
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
            for monomial, coefficient in sorted(polynomial.items()):
                result.update(repr((monomial, coefficient)).encode())
        return result.hexdigest()

    def save_checkpoint(self, path):
        """Use explicit exact terms; SymPy PolyRing pickling mutates its state."""
        def encode(polynomial):
            terms = []
            for powers, coefficient in sorted(polynomial.items()):
                sparse = [[i, int(power)] for i, power in enumerate(powers) if power]
                values = [coefficient] if self.field == QQ else coefficient.to_list()
                terms.append([sparse, [str(value) for value in values]])
            return terms

        state = dict(format="f2-sparse-ring-v1", branch=self.tag,
                     quadratic=self.field != QQ,
                     polynomial_backend=getattr(self.ring, "backend_name", "sympy"),
                     generators=[str(symbol) for symbol in self.ring.symbols],
                     nparameters=self.nparameters, receipt=self.receipt,
                     constraints=[(tag, encode(poly)) for tag, poly in self.constraints],
                     h=encode(self.h), v=encode(self.v), S=encode(self.S),
                     B=encode(self.B), kappa=encode(self.kappa),
                     kappa2=encode(self.kappa2), kappa3=encode(self.kappa3),
                     uk={key: encode(value) for key, value in self.uk.items()},
                     p=[encode(value) for value in self.p],
                     q=[encode(value) for value in self.q],
                     w=[encode(value) for value in self.w])
        if getattr(self, "compaction", None) is not None:
            state["compaction"] = self.compaction
        text = json.dumps(state, separators=(",", ":"))+"\n"
        Path(path).write_text(text, encoding="utf-8")
        print(f"checkpoint={path} bytes={len(text)} sha256={hashlib.sha256(text.encode()).hexdigest()}", flush=True)


def verify_fifth_q_formula():
    h, s, b, c3, c4, c5, t, w, k = sp.symbols("h s b c3 c4 c5 t w k")
    delta = s/h*w+b/h**2*w**2+c3/h**3*w**3+c4/h**3*w**4+c5/(t*h**3)*w**5
    direct = t**5*sp.expand(h**5*sum(sp.binomial(sp.Rational(5,3),j)*delta**j for j in range(6))).coeff(w,5)+k*(t*h)**4
    expected = (sp.Rational(5,3)*t**4*h**2*c5+t**5*(sp.Rational(10,9)*(h*s*c4+b*c3)-sp.Rational(5,27)*(s**2*c3+s*b**2)+sp.Rational(20,243)*s**3*b-sp.Rational(7,729)*s**5)+k*(t*h)**4)
    assert sp.cancel(direct-expected) == 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch", choices=("2+2", "2+1+1"), default="2+1+1")
    parser.add_argument("--through", type=int, choices=tuple(range(2,19)), default=5)
    parser.add_argument("--backend", choices=("sympy", "flint-mpoly"), default="sympy")
    parser.add_argument("--resume", help="resume from an exact sparse checkpoint")
    parser.add_argument("--receipt")
    parser.add_argument("--tuned", action="store_true")
    parser.add_argument("--checkpoint")
    parser.add_argument("--terminal", action="store_true")
    parser.add_argument("--terminal-q", action="store_true")
    parser.add_argument("--root-values", action="store_true")
    args = parser.parse_args()
    global POLY_BACKEND
    POLY_BACKEND = args.backend
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
    lift = (SparseLift.from_checkpoint(args.resume) if args.resume else
            SparseLift(args.branch, sp.expand(h), sp.expand(v), field))
    if args.resume and lift.tag != args.branch:
        raise ValueError("--branch disagrees with checkpoint")
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
        if args.checkpoint:
            lift.save_checkpoint(args.checkpoint)
    report = dict(branch=args.branch, field=str(field), through=args.through,
                  p10_support="full t^2 times degree<=12", blocks=lift.receipt,
                  digest=lift.digest(), tuned_contact=args.tuned, terminal_chain=args.terminal,
                  terminal_q=args.terminal_q,
                  root_values=args.root_values,
                  retained_constraint_count=len(lift.constraints),
                  retained_constraint_tags=[tag for tag, _ in lift.constraints],
                  scope="support-relaxed necessary system; not exclusion",
                  elapsed_seconds=time.monotonic()-lift.started)
    if args.receipt:
        Path(args.receipt).write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    if args.checkpoint:
        lift.save_checkpoint(args.checkpoint)
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
