# Exact first-five-band local split for the (2+1+1) branch

This note records a characteristic-zero necessary-condition reduction for the
quartic edge factor
\[
h(t)=(t-1)^2R(t),
\]
where (R) is squarefree of degree two and (R(1)\ne0). It treats the
distinguished double root at (t=1) and the two simple roots of (R). It is
not a whole-branch exclusion.

## First divisors

The third band is polynomial exactly when
\[
h^4\mid 18h^3U_1U_2-U_1^3.
\]
At a root of (h) of multiplicity (e), if (k=\operatorname{ord}(U_1)
<4e/3), then (U_1^3) has order (3k<4e) and strictly smaller order than
the other summand. Thus
\[
\operatorname{ord}_\alpha(U_1)\ge\left\lceil\frac{4e}{3}\right\rceil.
\]
At the two simple roots this gives order at least two. The attached terminal
root condition at the double root gives order at least five. With the degree
bound \(\deg U_1\le11\), every candidate therefore has
\[
U_1=(t-1)^5R^2B,\qquad \deg B\le2.
\]
The fourth-band condition is
\[
h^7\mid54h^6U_1U_3+27h^6U_2^2-9h^3U_1^2U_2+U_1^4.
\]
If (U_2) does not vanish at a root of (h) of multiplicity (e\in\{1,2\}),
the (27h^6U_2^2) term has order (6e<7e), strictly below the other terms.
So (U_2) vanishes at all three distinct roots and
\[
U_2=(t-1)RC,\qquad \deg C\le8,
\]
using \(\deg U_2\le11\).

## Fifth-band split at the double root

The fifth-band condition is
\[
h^{10}\mid N_5,
\]
where
\[
N_5=810h^9(U_1U_4+U_2U_3)
-135h^6(U_1^2U_3+U_1U_2^2)
+60h^3U_1^3U_2-7U_1^5.
\]
At either simple root, the already forced orders
\(\operatorname{ord}(U_1)\ge2\), \(\operatorname{ord}(U_2)\ge1\) make every
summand have order at least ten, so the condition is automatic there.

At the double root write (u=t-1), (k=\operatorname{ord}(U_1)\ge5),
\(q=\operatorname{ord}(U_2)\ge1), and (r=\operatorname{ord}(U_3)\ge0).
The six summand orders, in the order shown in (N_5), are
\[
18+k+\operatorname{ord}(U_4),\quad18+q+r,\quad12+2k+r,\quad
12+k+2q,\quad6+3k+q,\quad5k.
\]
The order threshold is 20. This gives exactly these local alternatives:

1. (q\ge2): every summand has order at least 20.
2. (q=1, k\ge6, r\ge1): every summand has order at least 20.
3. (q=1, k=5, r=0): the two order-19 terms must cancel.

The other possibilities are impossible by a unique order-19 term:
\(q=1,k=5,r\ge1\) leaves \(-135h^6U_1U_2^2\), and
\(q=1,k\ge6,r=0\) leaves \(810h^9U_2U_3\).

For the tied case write
\[
h=cu^2+O(u^3),\quad U_1=Au^5+O(u^6),\quad
U_2=Du+O(u^2),\quad U_3=E+O(u),
\]
with (cD\ne0). The order-19 coefficient is exactly
\[
135c^6D(6c^3E-AD),
\]
so (AD=6c^3E). The parameterizations above give
\(A=c^2B(1)\) and \(D=cC(1)\), hence the tied condition becomes
\[
B(1)C(1)=6U_3(1).
\]

Equivalently, the remaining local parameter space lies in the union
\[
\begin{array}{ll}
\text{(i)} & C(1)=0,\\
\text{(ii)}& C(1)\ne0,\ B(1)=0,\ U_3(1)=0,\\
\text{(iii)}& C(1)\ne0,\ B(1)\ne0,\ B(1)C(1)=6U_3(1).
\end{array}
\]
This is an exhaustive local split under the first five bands and the terminal
double-root condition. It has not been substituted into the complete
2+1+1 coefficient ideal, and none of these three loci has been excluded.

The later sixth-band local certificate already strengthens the double-root
contact to \(\operatorname{ord}(U_2)\ge2\). Re-running
[f2_75_125_k6_local_jet_certificate.py](../../src/f2_75_125_k6_local_jet_certificate.py)
verifies that the unique order-\(-5\) term is
\(-5U_2^3/(81h^4)\), which rules out \(q=1\). Thus only stratum (i),
\(C(1)=0\), remains in the current deeper-band model. This is a local jet
restriction already used by that model, not a whole-branch exclusion.

## Exact replay and scope

The SymPy replay checks the order-19 coefficient, all six-term valuation
censuses for the simple and double roots, and the translation of the tied
condition into (B,C,U_3). The script is
[f2_75_125_211_first5_band_root_split_20261021.py](../../src/f2_75_125_211_first5_band_root_split_20261021.py)
and the latest receipt is
[f2_75_125_211_first5_band_root_split_20261021_v2.json](../../data/json/f2_75_125_211_first5_band_root_split_20261021_v2.json).

The 2+1+1 whole-branch exclusion remains at **0%**. This local refinement
does not change the 2+2 status or its paused computations.
