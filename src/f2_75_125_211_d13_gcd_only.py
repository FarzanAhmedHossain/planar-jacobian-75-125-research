#!/usr/bin/env python3
"""Compute the exact common polynomial gcd of the d13 five-variable inputs."""

from __future__ import annotations

import contextlib
import io
import json
import sys
from pathlib import Path

import f2_75_125_211_d13_first_block_groebner as driver


report_path = Path(sys.argv[2])
# Reuse the driver's exact specialization, support filtering, and rational gcd
# calculation, but deliberately do not run Groebner after the gcd is computed.
driver.ring_groebner = lambda *_args, **_kwargs: None
with contextlib.redirect_stdout(io.StringIO()):
    driver.main()

report = json.loads(report_path.read_text(encoding="utf-8"))
report.update({
    "groebner_method": None,
    "monomial_order": None,
    "groebner_basis_length": None,
    "groebner_basis": None,
    "unit_ideal": None,
    "zero_dimensional": None,
    "scope": (
        "Only the exact common gcd of the selected equations was computed. "
        "No Groebner basis or ideal inconsistency was computed."
    ),
})
report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print(json.dumps({
    "input_equations": report["input_equations"],
    "unique_equations_up_to_nonzero_rational_scale": report[
        "unique_equations_up_to_nonzero_rational_scale"
    ],
    "common_gcd": report["common_gcd"],
    "scope": report["scope"],
}, indent=2))
print(f"report={report_path}")
