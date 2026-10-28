#!/usr/bin/env python3
"""Exact degree reduction from the K6-K10 root-jet split, F2 2+1+1.

This translates the certified root-order alternatives into divisibility
conditions on U1,...,U4 for h=(t-1)^2 R, with R squarefree quadratic.  It is a
necessary parameter reduction only; it does not eliminate any root-profile
stratum or the full 2+1+1 branch.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "f2_75_125_211_rootwise_divisor_reduction_20261021.json"
SOURCES = [
    HERE / "f2_75_125_k6_local_jet_certificate.py",
    HERE / "f2_75_125_k7_double_root_certificate.py",
    HERE / "f2_75_125_k8_local_jet_certificate.py",
    HERE / "f2_75_125_211_k10_double_root_profile_exclusion_20261021.py",
]
K10_RECEIPT = HERE / "f2_75_125_211_k10_double_root_profile_exclusion_20261021.json"
K10_RECEIPT_SHA256 = "2895a21bf58a01282e02ee6205a9cf9e56698b38e534793618bfd606240632e7"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def coeff_count(degree_bound: int) -> int:
    return max(0, degree_bound + 1)


def main() -> None:
    k10_receipt = json.loads(K10_RECEIPT.read_text(encoding="utf-8"))
    assert sha256(K10_RECEIPT) == K10_RECEIPT_SHA256
    assert k10_receipt["status"] == "EXACT_DOUBLE_ROOT_PROFILE_EXCLUDED_AT_K10"
    # R has two distinct roots alpha,beta, neither equal to the double root 1.
    # A root is "deep" when ord(U1)>=3 and ord(U2)>=2; otherwise the exhaustive
    # K6/K7 alternative is the tuned (2,1) jet. K8 adds U4(root)=0 at deep
    # simple roots. At the double root, K8 excludes ord(U2)=3 with ord(U3)=1,
    # and the exact K10 certificate excludes ord(U2)=3 with ord(U3)>=2. Since
    # K7 forces U3(1)=0, every surviving profile has ord(U2)>=4.
    profiles = []
    for deep_count, names in ((0, ["TT"]), (1, ["DT", "TD"]), (2, ["DD"])):
        d = deep_count
        bounds = {
            "U1": {"factor_degree": 9 + d, "quotient_degree": 2 - d},
            "U2": {"factor_degree": 6 + d, "quotient_degree": 5 - d},
            "U3": {"factor_degree": 1 + d, "quotient_degree": 10 - d},
            "U4": {"factor_degree": d, "quotient_degree": 11 - d},
        }
        counts = {name: coeff_count(row["quotient_degree"])
                  for name, row in bounds.items()}
        assert bounds["U1"]["factor_degree"] + bounds["U1"]["quotient_degree"] == 11
        assert bounds["U2"]["factor_degree"] + bounds["U2"]["quotient_degree"] == 11
        assert bounds["U3"]["factor_degree"] + bounds["U3"]["quotient_degree"] == 11
        assert bounds["U4"]["factor_degree"] + bounds["U4"]["quotient_degree"] == 11
        profiles.append({
            "profile_names": names,
            "deep_simple_root_count": d,
            "deep_root_factor": "D(t)=product of t-alpha over the deep simple roots",
            "necessary_forms": {
                "U1": "(t-1)^5 * R(t)^2 * D(t) * A(t)",
                "U2": "(t-1)^4 * R(t) * D(t) * B(t)",
                "U3": "(t-1) * D(t) * C(t)",
                "U4": "D(t) * E(t)",
            },
            "quotient_degree_bounds": {
                name: row["quotient_degree"] for name, row in bounds.items()
            },
            "free_quotient_coefficient_counts": counts,
            "total_U1_to_U4_coefficients": sum(counts.values()),
        })

    assert [p["free_quotient_coefficient_counts"] for p in profiles] == [
        {"U1": 3, "U2": 6, "U3": 11, "U4": 12},
        {"U1": 2, "U2": 5, "U3": 10, "U4": 11},
        {"U1": 1, "U2": 4, "U3": 9, "U4": 10},
    ]

    result = {
        "branch": "2+1+1",
        "carrier": "h=(t-1)^2 R, deg(R)=2, R squarefree, R(1)!=0",
        "root_jet_dependencies": {
            path.name: sha256(path) for path in SOURCES
        },
        "K10_profile_exclusion_receipt": {
            "path": K10_RECEIPT.name,
            "sha256": K10_RECEIPT_SHA256,
            "status": k10_receipt["status"],
        },
        "local_conditions_used": [
            "double root: ord(U1)>=5, ord(U2)>=4, U3(1)=0; K8 and K10 jointly exclude every ord(U2)=3 profile",
            "each simple root: either tuned (ord(U1)=2, ord(U2)=1) or deep (ord(U1)>=3, ord(U2)>=2, U3(alpha)=0)",
            "K8 at each deep simple root: U4(alpha)=0",
        ],
        "degree_caps": {"U1": 11, "U2": 11, "U3": 11, "U4": 11},
        "profile_count": sum(len(p["profile_names"]) for p in profiles),
        "profiles": profiles,
        "scope": "exact factor-degree arithmetic from K6-K10 root-jet conditions; necessary local parameterization only, not a whole-branch exclusion",
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "profiles": [
            (p["profile_names"], p["free_quotient_coefficient_counts"],
             p["total_U1_to_U4_coefficients"])
            for p in profiles
        ],
        "source_hashes": result["root_jet_dependencies"],
        "receipt": OUTPUT.name,
        "scope": result["scope"],
    }, indent=2))


if __name__ == "__main__":
    main()
