"""Analyze the exact support-38 family extracted from the d19 checkpoint."""

from pathlib import Path

from f2_75_125_211_analyze_d19_support32_family import main as analyze_family


if __name__ == "__main__":
    analyze_family(
        input_path=Path("research/f2_75_125_211_d19_support38_family_input_20261025.json"),
        output_path=Path("research/f2_75_125_211_d19_support38_analysis_20261025.json"),
    )
