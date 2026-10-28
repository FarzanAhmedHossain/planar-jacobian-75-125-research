# D24 (2+1+1, U_1=0): eight-parameter necessary subsystem

## Scope

This audit extracts every retained D24 constraint whose sparse support is contained in checkpoint generator indices

```text
{6, 16, 35, 36, 37, 38, 39, 40}
```

The checkpoint generator order identifies these with
`(p5, p15, p34, p35, p36, p37, p38, p39)`. The extracted 143 equations are a necessary subsystem of the D24 model; they are not the full D24 ideal and cannot by themselves establish branch exclusion.

## Exact computation

- D24 checkpoint SHA-256: `2369f36ce13d0d4336053ba0a7a546aa4407d04d77ee544aa4e313dbc88f2c79`.
- Extracted TSV SHA-256: `5663f9e4a61a0d7518a6fd1b46af93d9ce6404fd42cb5efbbe278ff497548b7d`.
- The exact replay is [f2_75_125_d24_u1zero_support8_groebner.py](../../src/f2_75_125_d24_u1zero_support8_groebner.py); it checks the input digest and computes a grevlex Gröbner basis over ℚ.
- SymPy 1.14.0 returned a three-polynomial basis of degrees (4,3,3), not the unit ideal; the ideal is not zero-dimensional. The origin satisfies all 143 equations.

One basis element is a nonzero rational multiple of

\[
L^3,\qquad L=224725p_{34}+1639725p_{35}+6284625p_{36}+15845845p_{37}+24712395p_{38}-163942748p_{39}.
\]

Thus any solution of this subsystem satisfies (L=0). Another basis element is a multiple of (L M), so it adds no condition after (L=0); the remaining basis element is quartic. This gives a compact necessary condition but no contradiction.

## Branch status

Branchwide excluded remains **0%** for (2+1+1, U_1=0). The calculation neither closes nor rules out the D24 branch; further constraints involving additional parameters must be coupled to this reduced subsystem.
