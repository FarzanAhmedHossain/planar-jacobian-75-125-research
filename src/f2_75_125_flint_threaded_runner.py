"""Run the exact sparse F2 lift with an explicit python-flint thread limit.

Usage (from the workspace root)::

    python research/f2_75_125_flint_threaded_runner.py --flint-threads 7 \
        --branch 2+1+1 --through 16 ...

The remaining arguments are forwarded unchanged to
``f2_75_125_sparse_ring_lift.py``. Set ``PYTHONPATH`` so python-flint can be
imported from the research dependency bundle before invoking this wrapper.
"""

from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--flint-threads", type=int, required=True)
    parser.add_argument(
        "--stream-resume", action="store_true",
        help="use the streaming checkpoint loader and writer",
    )
    options, forwarded = parser.parse_known_args()
    if options.flint_threads < 1:
        parser.error("--flint-threads must be positive")

    from flint import ctx

    ctx.threads = options.flint_threads
    target_name = ("f2_75_125_sparse_ring_lift_stream_resume.py"
                   if options.stream_resume
                   else "f2_75_125_sparse_ring_lift.py")
    target = Path(__file__).with_name(target_name)
    sys.argv = [str(target), *forwarded]
    print(f"flint_threads={ctx.threads} target={target}", flush=True)
    runpy.run_path(str(target), run_name="__main__")


if __name__ == "__main__":
    main()
