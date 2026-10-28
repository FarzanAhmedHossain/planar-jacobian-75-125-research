# Transport of the MLT D=75 principal edge to F2 bands

Date: 2026-10-18. This is a coordinate derivation for one Makar-Limanov–Trakhtenberg (MLT) D=75 Newton-resolution row. It does not exclude (75,125). The row is a conditional resolution branch unless its applicability is separately established for a given mate.

## Source row and its relation to h

The primary source is Makar-Limanov–Trakhtenberg, [*Properties of a Jacobian mate*, MPIM 2024-33](https://archive.mpim-bonn.mpg.de/5148/1/mpim-preprint_2024-33.pdf). The table preface says it records vertices before and after resolution and the principal edge before resolution (PDF text lines 1466–1478). Its D=75 row at PDF text lines 1977–1987 gives

\[
v_1=3(7/5,2)=(21/5,6),\quad v'_2=(4/5,1),\quad
\lambda_0=5/3,\quad
\phi_0=x(xy^5-r_1)^2q_0(xy^5),\quad \deg q_0=2.
\]

In the F2 normalization, \(t=xy^5\) and \(P_{15}=a t^3h(t)^3\), so the leading factor \(\phi_0\) is \(x h(t)\), up to a nonzero scalar. Thus this row selects a repeated root \(r_1\) of \(h\) and has \(h(t)=(t-r_1)^2q_0(t)\), with \(\deg q_0=2\).

An allowed nonzero diagonal rescaling of \(x,y\) sends \(t\) to a scalar multiple of itself, so normalize \(r_1=1\). The two quartic patterns are then

- \(2+2:\ h=(t-1)^2(t-\alpha)^2\), with \(\alpha\ne1\);
- \(2+1+1:\ h=(t-1)^2(t-\alpha)(t-\beta)\), with \(1,\alpha,\beta\) distinct.

The selected root is 1. In the 2+2 pattern this calculation does not automatically apply the same row at the other double root \(\alpha\).

## Endpoint line in the ramified F2 chart

For the selected root take \(y=x^{-1/5}+Y\), set \(v=x^{1/5}Y\), and then put \(x=X^5\), \(Y=X^{-1}u\), \(w=X^{-1}\). The root normalization gives \(t=xy^5=(1+v)^5\). Write
\[
P=\sum_{s\ge0}p_{15-s}(t)z^{15-s},\qquad z=y^{-1},
\]
and define the normalized band function
\[
B_s(t)=p_{15-s}(t)/(a t^3)
\]
near \(t=1\). The exact contribution from band \(s\), after the first Newton shift, is
\[
a x^{(15-s)/5}(1+v)^s B_s((1+v)^5).
\]
For the normalized bands used in the current F2 model, \(B_0=h^3\), \(B_s=U_s\) for \(1\le s\le4\), and \(B_5=U_5/t\) after the recorded \(x\)-translation gauge. In particular the chart slices are \(A_0(u)=h((1+u)^5)^3\), \(A_s(u)=(1+u)^sU_s((1+u)^5)\) for \(1\le s\le4\), and \(A_5(u)=U_5((1+u)^5)\). For later bands retain the definition \(B_s=p_{15-s}/(a t^3)\), which is regular near \(t=1\) even if it is Laurent at \(t=0\).
Consequently, if
\[
C_{s,j}=[v^j]\big((1+v)^s B_s((1+v)^5)\big),
\]
then \(C_{s,j}\ne0\) gives a support point \((\mu,j)=((15-s+j)/5,j)\) in the shifted \((x,Y)\)-chart.

The MLT endpoints determine
\[
25\mu-17j=3,
\]
since both \((21/5,6)\) and \((4/5,1)\) satisfy it. The supporting side relevant to the F2 terminal edge is \(25\mu-17j\le3\). Substituting \(\mu=(15-s+j)/5\) gives
\[
25\mu-17j=75-5s-12j\le3
\quad\Longleftrightarrow\quad
5s+12j\ge72.
\]
Thus every \(C_{s,j}\) with \(5s+12j<72\) must vanish. The endpoint coefficients at \((s,j)=(0,6)\) and \((12,1)\) are nonzero. The line and endpoints are direct table data; identifying the inward side with \(5s+12j\ge72\) uses the F2 complete-chain chart, whose terminal edge is recorded in the linked GGV audit.

In the no-extra-shear chart, define \(A_s(u)=(1+u)^sB_s((1+u)^5)\). The lower endpoint gives \(A_{12}(0)=B_{12}(1)=0\), while the nonzero boundary coefficient gives
\[
[u]A_{12}=12B_{12}(1)+5B'_{12}(1)=5B'_{12}(1)\ne0.
\]

The least permitted u-order in that chart is

| Band s | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Least j | 6 | 6 | 6 | 5 | 5 | 4 | 4 | 4 | 3 | 3 | 2 | 2 | 1 | 1 | 1 |

This is the same terminal normal and endpoint as the F2 complete-chain calculation in research/f2_terminal_chain_local_selectors_20261018.md: the terminal P weight is 3, the terminal Q weight is 5, and
\[
\lambda_0=(25-17)/3-1=5/3.
\]
The MLT row therefore gives a primary-source description of that same terminal edge. The paper's formula for \(\lambda_0\) is at PDF text lines 507–515.

## Coordinate changes must be retained

The preceding table is for the resolved chart. To pull it back to the existing normalized F2 bands, let
\[
F(u,w)=\sum_i A_i(u)w^i,\qquad A_i(u)=(1+u)^iB_i((1+u)^5),
\]
after removing the common nonzero scalar. The F2 chain permits earlier type-III source shears \(u\mapsto u+aw+bw^2\). Their exact centered s-th slice is
\[
T_s(u)=[w^s]F(u+aw+bw^2,w)
=\sum_{i=0}^{s}\ \sum_{m+2n=s-i}
\frac{a^m b^n}{m!n!}A_i^{(m+n)}(u).
\]
The transported edge constraint is
\[
[u^j]T_s=0\quad(5s+12j<72),\qquad [u]T_{12}\ne0.
\]
This retains the cases \(a=0\), \(b=0\), and \(ab\ne0\). Until the equivalence of the specific MLT chart with a gauge in which these shears vanish is proved, they must not be silently set to zero. The separate source gauge \(x\mapsto x+c\) also must be transported on lower slices; its compatibility with the terminal edge ideal is checked in research/f2_terminal_centering_identities_20261018.md.

For the first two bands, \(A_0(u)=h((1+u)^5)^3\) has order 6 and the exact centered slices begin
\[
T_1=A_1+aA'_0,\qquad
T_2=A_2+aA'_1+bA'_0+(a^2/2)A''_0.
\]
The low coefficients of these slices must vanish according to \(5s+12j<72\). If one additionally proves \(a=b=0\), the row gives raw orders \(\operatorname{ord}_{t=1}U_1,\operatorname{ord}_{t=1}U_2\ge6\), \(\operatorname{ord}_{t=1}U_3,\operatorname{ord}_{t=1}U_4\ge5\), and \(\operatorname{ord}_{t=1}U_5\ge4\). With the full F2 chain retained, the raw implications are weaker and the tangent relations already recorded in the terminal-chain notes apply.

## Effect on the two quartic branches

With \(a,b\) retained, this is the same terminal half-space already imposed by the GGHV-derived F2 model. It supplies no additional independent equation beyond that coefficient system, and the endpoint alone excludes neither quartic pattern.

Even in the special chart \(a=b=0\), the row's orders fit within the degree bounds. Combining the selected-root orders with the existing local double-root/simple-root divisibilities gives

- \(2+2:\ U_1=(t-1)^6(t-\alpha)^4S_1\), \(\deg S_1\le1\); \(U_2=(t-1)^6(t-\alpha)^2S_2\), \(\deg S_2\le3\).
- \(2+1+1:\ U_1=(t-1)^6(t-\alpha)^2(t-\beta)^2S_1\), \(\deg S_1\le1\); \(U_2=(t-1)^6(t-\alpha)(t-\beta)S_2\), \(\deg S_2\le3\).

These are only degree-compatible divisibility forms, not solutions of the complete Keller band system. They show that the endpoint does not force the current \(U_i\) degree bounds to fail. The no-shear s=12 condition \(B_{12}(1)=0,\ B'_{12}(1)\ne0\) is also compatible with the allowed P3 slice.

## Scope warning

MLT's D=75 table lists multiple alternative resolution rows. The \(\lambda_0=5/3\) row is one branch of that enumeration; its equations are not automatically universal for every hypothetical F2 mate unless this row is shown to occur (or the other rows are eliminated). In the branch where it does occur, the principal edge gives an exact necessary half-space, but that half-space is already present in the GGHV complete-chain model. This audit adds a primary-source cross-check and a coordinate transport, not a new global exclusion.

## Status of the steps

* **Source-checked:** the D=75 row data and table convention, and MLT's \(\lambda_0\) formula cited above.
* **Exact algebra:** the \(B_s\) band expansion, the endpoint-line substitution, the ceiling orders, and the Taylor formula for \(T_s\).
* **Cross-source identification:** the inward half-space and the free \(a,b\) transport are identified with the GGHV F2 terminal chart. The matching normal/endpoints are exact; identifying the precise MLT resolution coordinates with a gauge where \(a=b=0\) has not been proved here.
* **Conditional scope:** the \(\lambda_0=5/3\) row is one enumerated D=75 branch, not a universal equation for all hypothetical mates. No branch exclusion follows.

