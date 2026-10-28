"""Run the missing exact 2+2 weight-45/55 power-relation lift.

This extends the already audited 2+2 family through weights 45 and 55.  It is
still a support-relaxed necessary subsystem, not a whole-branch exclusion.
The other root of alpha^2 - 3*alpha + 1 is its Galois conjugate.
"""

import sympy as sp
from sympy.polys.domains import QQ
from pathlib import Path

from f2_75_125_power_relation_fifthbands import (
    audit_fifth_block,
    verify_new_band_formulas,
)


def main() -> None:
    verify_new_band_formulas()
    t = sp.symbols("t")
    sqrt5 = sp.sqrt(5)
    alpha = (3 + sqrt5) / 2
    h = sp.expand((t - 1) ** 2 * (t - alpha) ** 2)
    slope = (4 * alpha - 11) / 3
    intercept = 3 - alpha
    v = sp.expand(t * (t - 1) * (t - alpha) * (slope * t + intercept))
    U = t * h
    selector = sp.expand(5 * U * sp.diff(v, t) - 4 * sp.diff(U, t) * v - U)
    if selector != 0:
        raise AssertionError("the normalized 2+2 edge selector failed")
    residuals = audit_fifth_block(
        "2+2",
        h,
        v,
        QQ.algebraic_field(sqrt5),
        Path(__file__).with_name("f2_75_125_22_first_two_blocks_46_56.pkl"),
    )
    print(
        f"scope=one 2+2 Galois representative through corrected weights 45/55; "
        f"residual_count={len(residuals)}; not a degree-pair exclusion",
        flush=True,
    )


if __name__ == "__main__":
    main()
