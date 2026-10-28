"""Bounded exact SymPy F5B attempt on the saved reduced alpha-local ideal."""

from __future__ import annotations

import ctypes
import hashlib
import json
import subprocess
import sys
import time
from ctypes import wintypes
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parent
IDEAL_RECEIPT = ROOT / "f2_75_125_depth6_alpha_q4s2_aux_ideal_agent_20261018.json"
RESULT = ROOT / "f2_75_125_depth6_alpha_q4s2_aux_gb_agent_20261018.json"
ERROR_LOG = ROOT / "f2_75_125_depth6_alpha_q4s2_aux_gb_agent_20261018.stderr.txt"
TIMEOUT_SECONDS = 120
MEMORY_CAP_BYTES = 1024**3


class IO_COUNTERS(ctypes.Structure):
    _fields_ = [(name, ctypes.c_ulonglong) for name in (
        "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
        "ReadTransferCount", "WriteTransferCount", "OtherTransferCount",
    )]


class BASIC_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [
        ("PerProcessUserTimeLimit", ctypes.c_longlong),
        ("PerJobUserTimeLimit", ctypes.c_longlong),
        ("LimitFlags", wintypes.DWORD),
        ("MinimumWorkingSetSize", ctypes.c_size_t),
        ("MaximumWorkingSetSize", ctypes.c_size_t),
        ("ActiveProcessLimit", wintypes.DWORD),
        ("Affinity", ctypes.c_size_t),
        ("PriorityClass", wintypes.DWORD),
        ("SchedulingClass", wintypes.DWORD),
    ]


class EXTENDED_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [
        ("BasicLimitInformation", BASIC_LIMIT_INFORMATION),
        ("IoInfo", IO_COUNTERS),
        ("ProcessMemoryLimit", ctypes.c_size_t),
        ("JobMemoryLimit", ctypes.c_size_t),
        ("PeakProcessMemoryUsed", ctypes.c_size_t),
        ("PeakJobMemoryUsed", ctypes.c_size_t),
    ]


def compute() -> None:
    receipt = json.loads(IDEAL_RECEIPT.read_text(encoding="utf-8"))
    reduction = receipt["constant_linear_reduction"]
    variable_names = reduction["post_substitution_variables"]
    symbols = {name: sp.Symbol(name) for name in variable_names}
    locals_ = {"sqrt": sp.sqrt, **symbols}
    field = QQ.algebraic_field(sp.sqrt(5))
    generator_records = reduction["post_substitution_generators"]
    expressions = [sp.sympify(row["polynomial"], locals=locals_) for row in generator_records]

    # Recheck the exact saved term counts before starting F5B.
    polynomials = [sp.Poly(expr, *[symbols[n] for n in variable_names], domain=field)
                   for expr in expressions]
    saved_term_counts = [row["terms"] for row in generator_records]
    if [len(poly.terms()) for poly in polynomials] != saved_term_counts:
        raise AssertionError("saved reduced generator strings failed exact term-count replay")

    basis = sp.groebner(
        expressions, *[symbols[n] for n in variable_names], order="grevlex",
        domain=field, method="f5b",
    )
    basis_text = "\n".join(str(poly.as_expr()) for poly in basis.polys)
    result = {
        "status": "completed",
        "input_ideal_receipt": str(IDEAL_RECEIPT.relative_to(ROOT.parent)),
        "receipt_sha256": hashlib.sha256(IDEAL_RECEIPT.read_bytes()).hexdigest(),
        "input_checkpoint_sha256": receipt["input_state_sha256"],
        "input_generator_set_sha256": reduction["post_substitution_generator_set_sha256"],
        "cas": f"SymPy {sp.__version__}",
        "method": "sp.groebner(..., order='grevlex', method='f5b', domain=QQ(sqrt(5)))",
        "variable_count": len(variable_names),
        "generator_count": len(expressions),
        "input_term_count": sum(saved_term_counts),
        "basis_size": len(basis.polys),
        "basis_term_count": sum(len(poly.terms()) for poly in basis.polys),
        "is_unit": len(basis.polys) == 1 and basis.polys[0].as_expr() == 1,
        "basis_sha256": hashlib.sha256(basis_text.encode()).hexdigest(),
    }
    RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


def create_limited_job() -> int:
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.CreateJobObjectW.restype = wintypes.HANDLE
    job = kernel32.CreateJobObjectW(None, None)
    if not job:
        raise ctypes.WinError(ctypes.get_last_error())

    info = EXTENDED_LIMIT_INFORMATION()
    # JOB_OBJECT_LIMIT_PROCESS_MEMORY | JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
    info.BasicLimitInformation.LimitFlags = 0x00000100 | 0x00002000
    info.ProcessMemoryLimit = MEMORY_CAP_BYTES
    kernel32.SetInformationJobObject.argtypes = [
        wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD,
    ]
    kernel32.SetInformationJobObject.restype = wintypes.BOOL
    ok = kernel32.SetInformationJobObject(
        job, 9, ctypes.byref(info), ctypes.sizeof(info),
    )
    if not ok:
        raise ctypes.WinError(ctypes.get_last_error())
    return job


def supervise() -> None:
    if RESULT.exists():
        RESULT.unlink()
    ideal = json.loads(IDEAL_RECEIPT.read_text(encoding="utf-8"))
    reduction = ideal["constant_linear_reduction"]
    base_result = {
        "input_ideal_receipt": str(IDEAL_RECEIPT.relative_to(ROOT.parent)),
        "receipt_sha256": hashlib.sha256(IDEAL_RECEIPT.read_bytes()).hexdigest(),
        "input_checkpoint_sha256": ideal["input_state_sha256"],
        "input_generator_set_sha256": reduction["post_substitution_generator_set_sha256"],
        "cas": f"SymPy {sp.__version__}",
        "method": "F5B Groebner basis over Q(sqrt(5)) with grevlex order",
        "generator_count": reduction["post_substitution_equation_count"],
        "variable_count": reduction["post_substitution_variable_count"],
        "input_term_count": reduction["post_substitution_term_count"],
        "input_max_generator_terms": reduction["post_substitution_max_polynomial_terms"],
    }
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
    kernel32.AssignProcessToJobObject.restype = wintypes.BOOL
    kernel32.TerminateJobObject.argtypes = [wintypes.HANDLE, wintypes.UINT]
    kernel32.TerminateJobObject.restype = wintypes.BOOL
    job = create_limited_job()
    start = time.monotonic()
    with ERROR_LOG.open("w", encoding="utf-8") as err:
        process = subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()), "--compute"],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=err,
        )
        if not kernel32.AssignProcessToJobObject(job, wintypes.HANDLE(process._handle)):
            process.kill()
            process.wait()
            raise ctypes.WinError(ctypes.get_last_error())

        timed_out = False
        while process.poll() is None:
            if time.monotonic() - start >= TIMEOUT_SECONDS:
                timed_out = True
                kernel32.TerminateJobObject(job, 124)
                break
            time.sleep(1)
        exit_code = process.wait()
    kernel32.CloseHandle(wintypes.HANDLE(job))

    if RESULT.exists():
        result = json.loads(RESULT.read_text(encoding="utf-8"))
    else:
        result = base_result
    result.update({
        "timeout_seconds": TIMEOUT_SECONDS,
        "memory_cap_bytes": MEMORY_CAP_BYTES,
        "memory_cap_enforced_by_windows_job_object": True,
        "timed_out": timed_out,
        "child_exit_code": exit_code,
        "stderr_log": str(ERROR_LOG.relative_to(ROOT.parent)),
    })
    if timed_out:
        result["status"] = "timed_out"
    elif exit_code != 0 and result.get("status") != "completed":
        result["status"] = "failed"
    RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    if "--compute" in sys.argv:
        compute()
    else:
        supervise()
