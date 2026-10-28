#!/usr/bin/env python3
"""Run the exact d13 F5B calculation with a bounded cached grevlex key."""

from functools import lru_cache

from sympy.polys.orderings import ReversedGradedLexOrder

# SymPy repeatedly recomputes the same pure-Python grevlex key in F5B. Cache
# up to 250k monomials in this separate process; the uncached reference job
# remains available for comparison.
ReversedGradedLexOrder.__call__ = lru_cache(maxsize=250_000)(
    ReversedGradedLexOrder.__call__
)

from f2_75_125_211_d13_first_block_groebner import main

if __name__ == "__main__":
    main()
