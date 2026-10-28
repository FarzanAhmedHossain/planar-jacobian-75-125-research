# K6–K10 rootwise divisor reduction for F2 (75,125), 2+1+1

This note converts the exact K6, K7, K8 and K10 local valuation certificates into
degree-bounded polynomial parameterizations. It is an exact necessary
reduction; the four resulting profile strata have not been eliminated.

Let
\[
h(t)=(t-1)^2R(t),\qquad \deg R=2,
\]
with (R) squarefree and (R(1)\ne0). Denote its two simple roots by
\(\alpha,\beta\). The K6/K7 certificates give the exhaustive simple-root
alternatives:

* **Tuned** at a simple root: \(\operatorname{ord}(U_1)=2\) and
  \(\operatorname{ord}(U_2)=1\), with the K6/K7 normalized residue equations.
* **Deep** at a simple root: \(\operatorname{ord}(U_1)\ge3\),
  \(\operatorname{ord}(U_2)\ge2\), and \(U_3(\alpha)=0\). K8 further forces
  \(U_4(\alpha)=0\).

At the distinguished double root, the same certificates give
\(\operatorname{ord}_1(U_1)\ge5\) and \(U_3(1)=0\). K8 excludes the
subcase \(\operatorname{ord}_1(U_2)=3,\ \operatorname{ord}_1(U_3)=1\). The
exact K10 certificate excludes all remaining local data with
\(\operatorname{ord}_1(U_2)=3,\ \operatorname{ord}_1(U_3)\ge2\), for arbitrary
regular slices \(U_4,U_5\). Since the K7 root-order split already gives
\(\operatorname{ord}_1(U_2)\ge3\), every surviving profile therefore has
\(\operatorname{ord}_1(U_2)\ge4\).

Choose a subset (S\subseteq\{\alpha,\beta\}) of deep simple roots and put
\(D_S(t)=\prod_{\gamma\in S}(t-\gamma)\), (d=|S|\). Combining those local
orders with \(\deg U_i\le11\) yields
\[
\begin{aligned}
U_1&=(t-1)^5R^2D_S A,&\deg A&\le2-d,\\
U_2&=(t-1)^4RD_S B,&\deg B&\le5-d,\\
U_3&=(t-1)D_S C,&\deg C&\le10-d,\\
U_4&=D_S E,&\deg E&\le11-d.
\end{aligned}
\]
Thus the four root profiles have these quotient coefficient counts:

| Profile | (U_1) | (U_2) | (U_3) | (U_4) | Total |
|---|---:|---:|---:|---:|---:|
| TT | 3 | 6 | 11 | 12 | 32 |
| DT | 2 | 5 | 10 | 11 | 28 |
| TD | 2 | 5 | 10 | 11 | 28 |
| DD | 1 | 4 | 9 | 10 | 24 |

The tuned-root residue equations and all later band equations still constrain
these quotient coefficients; the open conditions distinguishing T from D
must also be retained. In particular, even DD only reduces (U_1) to a
single scalar multiple of \((t-1)^5R^3\); it does not show that scalar or the
other quotients vanish.

The arithmetic replay is
`f2_75_125_211_rootwise_divisor_reduction_20261021.py`; its receipt pins the
K6–K10 checker source hashes and the exact K10 exclusion receipt. This local
reduction has not yet been attached to the complete exact coefficient ideal,
so full-branch coverage remains 0%.
