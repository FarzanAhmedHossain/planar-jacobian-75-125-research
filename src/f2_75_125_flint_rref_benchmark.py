"""Exact DomainMatrix benchmark mirroring the F2 Q(sqrt(5)) lifts.

Set F2_BENCH_FLINT=1 to load the local python-flint package and select
SymPy's FLINT ground types before importing SymPy. The checksum is over the
exact reduced matrix, so the two modes can be compared for identical output.
"""

from __future__ import annotations

import hashlib
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
if os.environ.get("F2_BENCH_FLINT") == "1":
    sys.path.insert(0, str(HERE / "f2_75_125_flint_deps"))
    os.environ["SYMPY_GROUND_TYPES"] = "flint"

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix


def main() -> None:
    root5 = sp.sqrt(5)
    field = QQ.algebraic_field(root5)
    rows = []
    nrows, ncols = 48, 32
    for i in range(nrows):
        row = []
        for j in range(ncols):
            if (17 * i + 13 * j + i * j) % 9 < 3:
                a = ((17 * i + 13 * j + i * j) % 19) - 9
                b = ((7 * i + 11 * j + i * j) % 13) - 6
                den = ((i + 3 * j) % 7) + 1
                value = sp.Rational(a, den) + sp.Rational(b, den + 1) * root5
            else:
                value = sp.S.Zero
            row.append(field.from_sympy(value))
        rows.append(row)
    matrix = DomainMatrix(rows, (nrows, ncols), field)
    start = time.perf_counter()
    reduced, pivots = matrix.rref()
    elapsed = time.perf_counter() - start
    canonical = [
        [str(field.to_sympy(entry)) for entry in row]
        for row in reduced.to_list()
    ]
    payload = repr((canonical, pivots)).encode()
    print(f"ground_types={os.environ.get('SYMPY_GROUND_TYPES', 'default')}")
    print(f"qq_dtype={QQ.dtype}")
    print(f"field_dtype={field.dtype}")
    print(f"shape={nrows}x{ncols} rank={len(pivots)}")
    print(f"seconds={elapsed:.6f}")
    print(f"sha256={hashlib.sha256(payload).hexdigest()}")


if __name__ == "__main__":
    main()
