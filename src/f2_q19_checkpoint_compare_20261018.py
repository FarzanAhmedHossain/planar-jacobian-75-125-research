"""Read-only exact Q19 optimization audit on a saved depth-six 2+2 state.

Uses its own checkpoint decoder and terminal coefficient construction.
Does not import or mutate the primary sparse-lift engine.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path
import time

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--checkpoint",default=str(Path(__file__).with_name(
        "f2_75_125_sparse_ring_lift_22_terminal_roots6_q_state.json")))
    ap.add_argument("--receipt",default=str(Path(__file__).with_suffix(".json")))
    ap.add_argument("--skip-legacy",action="store_true")
    args=ap.parse_args()
    started=time.monotonic()
    raw=Path(args.checkpoint).read_bytes()
    data=json.loads(raw)
    assert data["branch"]=="2+2" and data["quadratic"]
    assert len(data["p"])==len(data["q"])==7
    K=QQ.algebraic_field(sp.sqrt(5))
    nparams=data["nparameters"]
    R,*gens=ring(["t",*(f"p{i}" for i in range(nparams))],K)
    t=gens[0]
    index={str(symbol):i for i,symbol in enumerate(R.symbols)}
    old_names=data["generators"]
    @lru_cache(None)
    def coefficient(values):
        return K.new([QQ.convert(sp.Rational(value)) for value in values])
    def decode(encoded):
        terms={}
        for sparse,values in encoded:
            powers=[0]*len(gens)
            for old_i,power in sparse:
                name=old_names[old_i]
                assert name in index, f"unexpected active generator {name}"
                powers[index[name]]=power
            terms[tuple(powers)]=coefficient(tuple(values))
        return R.from_dict(terms)
    def const(value):
        return R.ground_new(K.convert(value))
    h,S,B,kappa=(decode(data[key]) for key in ("h","S","B","kappa"))
    p=[decode(value) for value in data["p"]]
    q=[decode(value) for value in data["q"]]
    print(f"decoded_depth6 checkpoint_bytes={len(raw)} seconds={time.monotonic()-started:.3f}",flush=True)
    assert p[0]==(t*h)**3 and q[0]==(t*h)**5
    assert p[1]==t**3*h**2*S and p[2]==t**3*h*B
    U3=p[3].exquo(t**3)
    U4=p[4].exquo(t**3)
    U5=p[5].exquo(t**2)
    U6=p[6].exquo(t**2)
    d=3*B-S**2
    delta=27*U3+2*S**3-9*S*B
    N=delta**2-3*d**3
    J,remainder=N.div(h)
    remainder_degrees=sorted({powers[0] for powers in remainder})
    assert not remainder
    q19=(const(sp.Rational(5,3))*(t*h)**2*p[6]
          +const(sp.Rational(10,9))*t**2*h*S*p[5]
          +t**2*(const(sp.Rational(10,9))*B-const(sp.Rational(5,27))*S**2)*p[4]
          +const(sp.Rational(5,6561))*t**5*J
          +const(sp.Rational(4,3))*kappa*t**4*h**3*S)
    saved_difference=q19-q[6]
    assert not saved_difference
    assert all(4<=powers[0]<=23 for powers in q19)
    assert all(2<=powers[0]<=13 for powers in p[6])
    print(f"optimized_Q19_matches_saved terms={len(q19)} h_remainder=0 seconds={time.monotonic()-started:.3f}",flush=True)

    # Independently reconstruct all terminal P coefficients. The local
    # variable u is represented temporarily by t in this separate ring;
    # every coefficient in a,b is checked to be t-independent.
    def tc(poly,degree):
        return R.from_dict({(0,*powers[1:]):c for powers,c in poly.items()
                            if powers[0]==degree})
    local=[]
    for s,poly in enumerate(p):
        out=R.zero
        for powers,c in poly.items():
            exponent=5*powers[0]-15+s
            assert exponent>=0
            for j in range(min(exponent,9)+1):
                out+=R.term_new((j,*powers[1:]),c*comb(exponent,j))
        local.append(out)
    C=tc(local[0],6)
    assert C and all(not any(powers) for powers in C)
    scalar=C[(0,)*len(gens)]
    a=-tc(local[1],5)/(6*scalar)
    r2=local[2]+a*local[1].diff(t)+a**2*local[0].diff(t).diff(t)/2
    b=-tc(r2,5)/(6*scalar)
    assert all(powers[0]==0 for poly in (a,b) for powers in poly)
    terminal=[]
    for s in range(7):
        bound=(71-5*s)//12
        expression=R.zero
        for i in range(s+1):
            for nn in range((s-i)//2+1):
                mm=s-i-2*nn
                derivative=local[i]
                for _ in range(mm+nn):
                    derivative=derivative.diff(t)
                expression+=derivative*a**mm*b**nn*K.convert(
                    sp.Rational(1,factorial(mm)*factorial(nn)))
        rows=[tc(expression,j) for j in range(bound+1)]
        assert all(not row for row in rows), f"terminal P distance {s} failed"
        terminal.append(dict(distance=s,positions=len(rows),nonzero=0))
    assert sum(record["positions"] for record in terminal[:6])==32
    print(f"terminal_P_32_and_depth6_rows_all_zero seconds={time.monotonic()-started:.3f}",flush=True)

    # Reconstruct the legacy corrected-power numerator by explicit sixth
    # convolutions. Comparing it to 3M^10*q19 proves that its exact quotient
    # is q19 and its remainder is zero, without repeating the slow division.
    legacy_status="not_run"
    if not args.skip_legacy:
        p0,p1,p2,p3,p4,p5,p6=p
        q0,q1,q2,q3,q4,q5=q[:6]
        p5c=(5*p0**4*p6+20*p0**3*(p1*p5+p2*p4)+10*p0**3*p3**2
              +30*p0**2*p1**2*p4+60*p0**2*p1*p2*p3+10*p0**2*p2**3
              +20*p0*p1**3*p3+30*p0*p1**2*p2**2+5*p1**4*p2)
        q3c=(6*q0*(q1*q5+q2*q4)+3*q0*q3**2+3*q1**2*q4
              +6*q1*q2*q3+q2**3)
        correction=3*kappa*(3*p0**2*p1*q0+p0**3*q1)
        numerator=p5c+correction-q3c
        legacy_residual=numerator-3*(t*h)**10*q19
        assert not legacy_residual
        legacy_status="exact numerator=3*(t*h)^10*optimized_Q19; quotient agrees and remainder zero"
        print(f"legacy_quotient_identity_exact seconds={time.monotonic()-started:.3f}",flush=True)
    payload=dict(
        checkpoint=str(Path(args.checkpoint).resolve()),checkpoint_sha256=sha256(raw).hexdigest(),
        checkpoint_bytes=len(raw),branch=data["branch"],depth=6,
        allocated_parameters=nparams,saved_constraints=len(data.get("constraints",[])),
        pole_numerator_terms=len(N),pole_quotient_terms=len(J),
        h_remainder_coefficient_positions_checked=4,h_remainder_nonzero_degrees=remainder_degrees,
        optimized_Q19_terms=len(q19),saved_Q19_terms=len(q[6]),
        optimized_vs_saved_difference_terms=len(saved_difference),
        legacy_quotient_comparison=legacy_status,
        P_terminal_positions_through5=32,P_terminal_positions_through6=36,
        P_terminal_nonzero_positions=0,P_terminal_checks=terminal,
        elapsed_seconds=round(time.monotonic()-started,3),
        checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        status="PASS" if not args.skip_legacy else "PASS_optimized_and_saved_only",
        scope="read-only exact optimization audit on this checkpoint; no new branch exclusion")
    Path(args.receipt).write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2),flush=True)


if __name__=="__main__":
    main()
