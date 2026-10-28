"""Stream exact D24 U1=0 constraints into nested sparse support windows.

The 6.9 GB checkpoint stores 1,999 compatibility polynomials before its much
larger recurrence state.  This reader consumes only that prefix, one sparse
term at a time, and emits constraints supported entirely on the D10, D11, and
D12 parameter windows.  It never materializes the checkpoint or a full
constraint polynomial in RAM.

The extracted systems are necessary subsystems only.  A proper ideal does not
exclude the F2 branch; only a unit-ideal certificate would do so.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
import time
from pathlib import Path

from f2_75_125_211_streaming_checkpoint_scan import JsonStream


ROOT = Path(r"D:\codex-f2-75-125")
CHECKPOINT = ROOT / (
    "f2_75_125_sparse_ring_lift_211_roots24_u1zero_"
    "support3line_streamw_t1_r230_20261028_memrelease_state.json"
)
RECEIPT = CHECKPOINT.with_name(CHECKPOINT.name.replace("_state.json", "_receipt.json"))
EXPECTED_CHECKPOINT_SHA256 = (
    "2369f36ce13d0d4336053ba0a7a546aa4407d04d77ee544aa4e313dbc88f2c79"
)
EXPECTED_CHECKPOINT_BYTES = 6_939_873_971

WINDOWS = {
    "support8": frozenset((6, 16, *range(35, 41))),
    "window18": frozenset((6, 16, *range(25, 30), *range(35, 41), *range(45, 50))),
    "window23": frozenset((6, 16, *range(25, 30), *range(35, 41),
                           *range(45, 50), *range(54, 59))),
}
EXPECTED_EXISTING = {
    "support8": (
        Path(__file__).with_name("f2_75_125_d24_u1zero_support8_constraints.tsv"),
        "5663f9e4a61a0d7518a6fd1b46af93d9ce6404fd42cb5efbbe278ff497548b7d",
        143,
    ),
    "window18": (
        Path(__file__).with_name("f2_75_125_d24_u1zero_window18_constraints.tsv"),
        "515046d7e36dd3c1d39972d7071b174a4b8356a6d8339a1e6126b9f36763be01",
        286,
    ),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def serialize_term(sparse: list, values: list) -> tuple[str, frozenset[int]]:
    if len(values) != 1:
        raise ValueError("the pinned U1=0 checkpoint must have rational coefficients")
    factors = []
    active = set()
    for index, exponent in sparse:
        index, exponent = int(index), int(exponent)
        if exponent:
            active.add(index)
            factors.append(f"{index}^{exponent}")
    monomial = ",".join(factors) if factors else "1"
    return f"{values[0]}:{monomial}", frozenset(active)


def extract() -> dict:
    if CHECKPOINT.stat().st_size != EXPECTED_CHECKPOINT_BYTES:
        raise ValueError("checkpoint byte size changed; refusing stale-source extraction")
    receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
    if receipt.get("checkpoint_sha256") != EXPECTED_CHECKPOINT_SHA256:
        raise ValueError("receipt does not bind the expected checkpoint digest")
    if receipt.get("scope") != "support-relaxed necessary system; not exclusion":
        raise ValueError("unexpected checkpoint receipt scope")

    outputs = {
        name: Path(__file__).with_name(
            f"f2_75_125_d24_u1zero_{name}_verified_reextract.tsv"
        )
        for name in WINDOWS
    }
    outputs["window23"] = Path(__file__).with_name(
        "f2_75_125_d24_u1zero_d12_window23_verified_constraints.tsv"
    )
    for path in outputs.values():
        if path.exists():
            raise FileExistsError(f"refusing to overwrite {path}")

    started = time.monotonic()
    counts = {name: 0 for name in WINDOWS}
    terms = {name: 0 for name in WINDOWS}
    max_window = WINDOWS["window23"]
    reader = JsonStream(CHECKPOINT)
    handles = {}
    try:
        metadata, reader = reader.top_level_constraints()
        if metadata.get("format") != "f2-sparse-ring-v1" or metadata.get("branch") != "2+1+1":
            raise ValueError("unexpected checkpoint metadata")
        if metadata.get("quadratic") is not False:
            raise ValueError("the U1=0 source is expected over QQ")
        if metadata.get("nparameters") != 224:
            raise ValueError("unexpected active parameter count")
        if metadata.get("generators", [None])[6] != "p5":
            raise ValueError("generator-index mapping changed")
        for name, path in outputs.items():
            handles[name] = path.open("w", encoding="utf-8", newline="\r\n")

        first_constraint = True
        scanned = 0
        while reader.peek() != "]":
            if not first_constraint:
                reader.take(",")
            first_constraint = False
            reader.take("[")
            tag = reader.value()
            reader.take(",")
            reader.take("[")

            active: set[int] = set()
            keep_terms = True
            term_count = 0
            with tempfile.SpooledTemporaryFile(
                max_size=8 << 20, mode="w+t", encoding="utf-8", newline="",
                dir=Path(__file__).parent,
            ) as spool:
                first_term = True
                while reader.peek() != "]":
                    if not first_term:
                        reader.take(",")
                    first_term = False
                    reader.take("[")
                    sparse = reader.value()
                    reader.take(",")
                    values = reader.value()
                    reader.take("]")
                    serialized, term_support = serialize_term(sparse, values)
                    active.update(term_support)
                    term_count += 1
                    if keep_terms and active.issubset(max_window):
                        if term_count > 1:
                            spool.write(";")
                        spool.write(serialized)
                    elif keep_terms:
                        keep_terms = False
                        spool.seek(0)
                        spool.truncate(0)

                reader.take("]")
                reader.take("]")
                if keep_terms:
                    for name, support in WINDOWS.items():
                        if active.issubset(support):
                            handle = handles[name]
                            handle.write(f"{tag}\t")
                            spool.seek(0)
                            shutil.copyfileobj(spool, handle, length=1 << 20)
                            handle.write("\n")
                            counts[name] += 1
                            terms[name] += term_count

            scanned += 1
            if scanned % 250 == 0:
                print(
                    f"scanned={scanned} counts={counts} "
                    f"checkpoint_char={reader.absolute_position()} "
                    f"elapsed={time.monotonic()-started:.1f}s",
                    flush=True,
                )

        reader.take("]")
        if scanned != receipt.get("retained_constraint_count"):
            raise AssertionError(
                f"receipt expected {receipt.get('retained_constraint_count')} constraints, got {scanned}"
            )
    finally:
        reader.close()
        for handle in handles.values():
            handle.close()

    report = {
        "checkpoint": str(CHECKPOINT),
        "checkpoint_bytes": CHECKPOINT.stat().st_size,
        "checkpoint_sha256_from_bound_receipt": EXPECTED_CHECKPOINT_SHA256,
        "receipt_scope": receipt["scope"],
        "constraints_scanned": scanned,
        "windows": {},
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }
    for name, path in outputs.items():
        report["windows"][name] = {
            "indices": sorted(WINDOWS[name]),
            "constraints": counts[name],
            "terms": terms[name],
            "path": str(path),
            "sha256": sha256(path),
            "scope": "necessary support subsystem; not an exclusion",
        }

    for name, (existing, expected_hash, expected_count) in EXPECTED_EXISTING.items():
        extracted = outputs[name]
        match = (
            sha256(existing) == expected_hash
            and counts[name] == expected_count
            and sha256(extracted) == expected_hash
        )
        report["windows"][name]["replay_matches_pinned_artifact"] = match
        report["windows"][name]["pinned_path"] = str(existing)
        report["windows"][name]["pinned_sha256"] = expected_hash
        if not match:
            raise AssertionError(f"{name} extraction failed pinned-artifact replay")

    result = Path(__file__).with_name("f2_75_125_d24_u1zero_d12_window23_extraction.json")
    result.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == "__main__":
    extract()
