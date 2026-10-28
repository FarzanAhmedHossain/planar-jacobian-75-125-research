"""Run the hash-bound d19 streaming extractor on the support-38 family."""

from pathlib import Path

from f2_75_125_211_extract_d19_support32_family import main as extract_family


if __name__ == "__main__":
    extract_family(
        target_size=38,
        output=Path("research/f2_75_125_211_d19_support38_family_input_20261025.json"),
    )
