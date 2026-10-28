#!/usr/bin/env python3
"""Continue the exact U1=0 simple-root profile screen from K35 to K40."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import f2_75_125_211_zero_u1_simple_root_k25_k35_screen as prior_screen


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "f2_75_125_211_zero_u1_simple_root_k25_k35_screen.json"
OUT = HERE / "f2_75_125_211_zero_u1_simple_root_k36_k40_screen.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    prior = json.loads(SOURCE.read_text(encoding="utf-8"))
    profiles = [tuple(row) for row in prior["survivor_profiles"]]
    expected = prior["profiles_surviving_K35"]
    if len(profiles) != expected or expected != 18_876:
        raise AssertionError("K35 survivor receipt count is inconsistent")

    by_band = {}
    examples = {}
    survivors = profiles
    for n in range(36, 41):
        main, corrections = prior_screen.bands_for(n)
        survivors, excluded, sample = prior_screen.screen_band(
            survivors, main, corrections
        )
        by_band[str(n)] = {"newly_excluded": excluded, "survivors": len(survivors)}
        examples[str(n)] = sample
        print(f"U1=0 simple-root K{n}: excluded={excluded} survivors={len(survivors)}", flush=True)

    result = {
        "scope": "2+1+1 U1=0 simple-root profiles, exact first-integral valuation screen K36-K40",
        "field": "exact rational generalized-binomial coefficients and integer valuations",
        "input_receipt": SOURCE.name,
        "input_sha256": sha256(SOURCE),
        "initial_profile_count": prior["initial_profile_count"],
        "profiles_entering_K36": len(profiles),
        "new_exclusions_by_band": by_band,
        "K9_K40_excluded": prior["K9_K35_excluded"] + sum(
            entry["newly_excluded"] for entry in by_band.values()
        ),
        "K9_K40_excluded_percent": round(
            100
            * (
                prior["K9_K35_excluded"]
                + sum(entry["newly_excluded"] for entry in by_band.values())
            )
            / prior["initial_profile_count"],
            6,
        ),
        "profiles_surviving_K40": len(survivors),
        "survivor_profiles": [list(row) for row in survivors],
        "examples_by_band": examples,
        "interpretation": (
            "This is only a finite local order-profile screen. Surviving profiles "
            "represent necessary valuation patterns, not polynomial coefficient solutions."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in result.items() if key != "survivor_profiles"}, indent=2), flush=True)
    print(f"receipt={OUT} sha256={sha256(OUT)}", flush=True)


if __name__ == "__main__":
    main()
