#!/usr/bin/env python3
"""Run the exact F2 sparse lift with FLINT's requested thread limit."""

from __future__ import annotations

import os
import runpy
import sys
from pathlib import Path

import flint

THREADS = min(8, os.cpu_count() or 1)
flint.ctx.threads = THREADS

target = Path(__file__).with_name("f2_75_125_sparse_ring_lift_stream_resume.py")
sys.argv = [str(target), *sys.argv[1:]]
print(f"flint_threads={THREADS} target={target.resolve()}", flush=True)
runpy.run_path(str(target), run_name="__main__")
