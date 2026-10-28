#!/usr/bin/env python3
"""Small streaming audit of p79 specialization/projection in checkpoints."""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


TARGET = {"p1", "p7", "p15", "p16", "p17"}


def main() -> None:
    path = Path(sys.argv[1])
    reader = JsonStream(path)
    try:
        metadata, _ = reader.top_level_constraints()
        generators = metadata["generators"]
        p79 = generators.index("p79")
        target_indices = {generators.index(name): name for name in TARGET}
        hist: Counter[str] = Counter()
        p79_terms = 0
        first_rows = []
        constraints = 0
        for tag, terms in reader.array_values():
            constraints += 1
            if constraints <= 4:
                term_preview = []
                for sparse, values in terms[:2]:
                    term_preview.append({
                        "sparse": [(generators[i], e) for i, e in sparse],
                        "values": values,
                    })
                first_rows.append({"tag": tag, "preview": term_preview})
            support = set()
            for sparse, _values in terms:
                if any(i == p79 and e for i, e in sparse):
                    p79_terms += 1
                    continue
                if any(i not in target_indices and e for i, e in sparse):
                    continue
                support.update(target_indices[i] for i, e in sparse if e)
            hist["".join(sorted(support)) or "<empty>"] += 1
        print(json.dumps({
            "path": str(path),
            "nparameters": metadata.get("nparameters"),
            "ngenerators": len(generators),
            "p79_index": p79,
            "target_indices": {name: generators.index(name) for name in sorted(TARGET)},
            "constraints": constraints,
            "terms_with_p79": p79_terms,
            "projected_support_histogram": dict(hist),
            "first_rows": first_rows,
        }, indent=2))
    finally:
        reader.close()


if __name__ == "__main__":
    main()
