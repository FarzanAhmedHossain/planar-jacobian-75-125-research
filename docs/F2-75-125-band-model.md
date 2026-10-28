# Working band model for the (75,125) degree pair

## Scope

This note records a finite outer support model and band equations for the \((75,125)\) degree pair. A written, minimality-free front-end reduction in the upstream proof bundle carries every hypothetical pair of these degrees to the standard \(F_2,j=1\) row; I inspected that extraction and replayed its custody checker. The finite coefficient system below is therefore exhaustive subject to the cited GGV/Makar-Limanov inputs and the written proof extraction. The band calculations exclude several local double-root strata through K15, audit the first Q-only top-y compatibility block including all its resonances, and use \(Q^3-\lambda P^5\) to reduce the quartic edge carrier to two algebraic root configurations. A Magnus-IV branchwise cube-support calculation forces the pre-generator \(R\) into a 76-monomial edge region and gives 11 finite \((R,Z,Q)\) support models. A new high-y Euclidean expansion reduces \(Q\) to a 708-term tail plus two scalars and forces the top slice of \(R\) to be a pure fifth power in the four highest Magnus branches. These are further necessary conditions only; the pair is not excluded.

## The degree-25 homogeneous carrier and factor-count branches

For a hypothetical pair with total degrees $(75,125)$ and constant
Jacobian, the top homogeneous forms $P_{75}$ and $Q_{125}$ satisfy
$[P_{75},Q_{125}]=0$. In characteristic zero, writing them in an affine
coordinate on $\mathbf P^1$ gives the corresponding one-variable relation
$3fg'-5f'g=0$. Therefore $g^3=c f^5$ for a nonzero scalar $c$, and unique
factorization gives $P_{75}=aH^3$ and $Q_{125}=bH^5$ for a binary
homogeneous form $H$ of degree $25$ and nonzero scalars $a,b$. The pair
$(H^3,H^5)$ itself has Jacobian zero, so it is only the principal
homogeneous part, not a Keller pair by itself.

Over an algebraic closure, every such $H$ factors into linear forms. Counting
distinct factors makes the split exhaustive: one factor $H=cL^{25}$; two
factors $H=cL^eM^{25-e}$, $1\le e\le24$; or at least three distinct
factors with positive multiplicities summing to $25$. In particular,
\[
x^{25}+y^{25}=\prod_{j=0}^{24}(x-\zeta_jy),
\qquad \zeta_j^{25}=-1,
\]
has 25 distinct linear factors in characteristic zero. It lies in the
three-or-more-factor branch, which is excluded by Abhyankar's at-most-two-
points-at-infinity theorem for plane Jacobian pairs. More generally, this
rules out the entire three-or-more-factor stratum. The cited one-point result
does not by itself exclude an individual pair whose leading form has one
factor: Abhyankar's Theorem (19.4) uses a universally quantified implication
from every Jacobian pair having one point at infinity to the Jacobian
conjecture. That scope cannot be replaced by a pointwise conclusion for one
pair. Likewise, Jung--van der Kulk degree divisibility is a consequence for
an already automorphic pair, not a standalone exclusion of a hypothetical
counterexample. Accordingly, this ordinary leading-form audit closes the
three-or-more-factor stratum; the one-factor and exactly-two-factor strata
need separate, scope-matched arguments.

This factor-count split is **not** the same as the active F2 split below.
The F2 edge reduction uses a quartic polynomial $h(t)$ in the edge
coordinate and, after the complete-chain root condition, its only remaining
root types are $2+2$ and $2+1+1$. The quartic $h$ is not the original
degree-25 homogeneous $H$: these are finer residual cases after the
Newton-chain reduction of the surviving exactly-two-factor H stratum. They
must still be excluded using the lower-band coefficient equations.

## Published corner data

GGV's table gives, for \(F_2\),

\[
A_0=(5,20),\quad A'_0=(1,0),\quad (m,n)=(j+2,2j+3).
\]

At \(j=1\), this is \((m,n)=(3,5)\), with degrees \(3(5+20)=75\) and \(5(5+20)=125\). The two edge endpoints in the scaled Newton diagrams are therefore

\[
P: (3,0)\longrightarrow(15,60),
\qquad
Q: (5,0)\longrightarrow(25,100).
\]

The edge direction is \((1,5)\). Put

\[
t=xy^5,\qquad z=y^{-1}.
\]

Then \(x=tz^5\), \(y=z^{-1}\), and

\[
x^iy^j=t^i z^{5i-j}.
\]

Thus the displayed \(P\)- and \(Q\)-edges are the \(z^{15}\)- and \(z^{25}\)-bands, respectively:

\[
P_{15}(t)=t^3f(t),\quad \deg f=12;
\qquad
Q_{25}(t)=t^5g(t),\quad \deg g=20.
\]

The endpoint conditions give nonzero constant and leading coefficients for \(f\) and \(g\).

## The top band equation

For band decompositions \(P=\sum_k p_k(t)z^k\) and \(Q=\sum_l q_l(t)z^l\), direct use of the chain rule gives

\[
[P,Q]_{x,y}
=\sum_{k,l}\bigl(kp_kq_l'-lp_k'q_l\bigr)z^{k+l-4}.
\]

Hence, if \([P,Q]=1\), every coefficient band other than \(z^0\) vanishes. If the support envelope really gives \(k\le15\) for \(P\) and \(l\le25\) for \(Q\), then \(z^{36}\) is the unique top output band, contributed only by \((k,l)=(15,25)\). Its coefficient is

\[
15P_{15}Q_{25}'-25P_{15}'Q_{25}
=5t^8\bigl(3fg'-5f'g\bigr).
\]

Therefore

\[
3fg'-5f'g=0.
\]

Equivalently, \((g^3/f^5)'=0\), so \(g^3=c f^5\) for \(c\ne0\). Unique factorization over \(\overline K[t]\) then gives

\[
f=a h^3,\qquad g=b h^5,\qquad \deg h=4,
\]

for nonzero constants \(a,b\). Thus the common edge is a power pair:

\[
P_{15}=a\,C^3,\qquad Q_{25}=b\,C^5,
\qquad C=xh(xy^5).
\]

This is a concrete rigidity reduction: the two edge polynomials, initially 13 and 21 coefficients, are controlled by one quartic \(h\) and two scalars. It does not yet constrain the lower bands enough to rule out a Keller pair.

## The first adjacent band still has a large kernel

Write the next bands as

\[
P_{14}=a\,t^3U(t),\quad \deg U\le11;
\qquad
Q_{24}=b\,t^5V(t),\quad \deg V\le19.
\]

Only \((15,24)\) and \((14,25)\) contribute to output band \(z^{35}\). Its equation is

\[
15P_{15}Q_{24}'-24P_{15}'Q_{24}
+14P_{14}Q_{25}'-25P_{14}'Q_{25}=0.
\]

After substituting \(P_{15}=a\,t^3h^3\) and \(Q_{25}=b\,t^5h^5\), and dividing by the nonzero scalar \(ab\), this becomes

\[
t^7h^7\,D\!\left(3\frac{V}{h^5}-5\frac{U}{h^3}\right)=0,
\qquad
D(w):=hw+th'w+5thw'.
\]

For \(w\in\overline K(t)\), \(D(w)=0\) implies \((thw^5)'=0\), hence \(thw^5\) is constant. Since \(h(0)\ne0\), the valuation at \(t=0\) rules out a nonzero rational \(w\): \(1+5\operatorname{ord}_0(w)\) cannot be \(0\). Thus \(w=0\), and the entire \(z^{35}\) condition reduces to

\[
3V=5h^2U.
\]

It leaves all 12 coefficients of \(U\) free; \(V\) is then determined and automatically has degree at most \(19\). So the first neighboring band does **not** pin down \(h\) or eliminate the new band variables. The next useful constraints must come from lower output bands, where these free parameters interact with further unknown bands.

## The next band gives a root-multiplicity condition

Write the following bands as \(P_{13}=a\,t^3W\), \(Q_{23}=b\,t^5Z\), with \(\deg W\le11\), \(\deg Z\le19\). Output band \(z^{34}\) receives exactly the pairs \((15,23),(14,24),(13,25)\). After using \(3V=5h^2U\), the coefficient equation is

\[
2h^7D_{5/2}\!\left(
3\frac{Z}{h^5}-5\frac{W}{h^3}-\frac53\frac{U^2}{h^6}
\right)=0,
\qquad
D_{5/2}(w):=hw+th'w+\frac52thw'.
\]

If \(D_{5/2}(w)=0\) for rational \(w\), then \((t^2h^2w^5)'=0\). Since \(h(0)\ne0\), the order at \(t=0\) rules out a nonzero rational solution. Therefore the band equation is equivalent to

\[
9hZ-15h^3W=5U^2.
\]

In particular, \(h\mid U^2\). This condition is also sufficient for this one band equation: if \(h\mid U^2\), choose \(W=0\) and \(Z=5U^2/(9h)\), which has degree at most \(18\). Thus the first two lower-band equations yield an exact, local divisibility test rather than an exclusion.

Factoring \(h\) over \(\overline K\) turns that test into a small case split. If a root of \(h\) has multiplicity \(e\), then \(U\) must vanish there to order at least \(\lceil e/2\rceil\). In the squarefree quartic case, \(h\mid U\); for repeated-root patterns the required vanishing is weaker. This makes the factor-pattern idea more relevant for the quartic edge factor \(h\) than for the degree-25 homogeneous leading carrier.

For a quartic \(h\), the first divisibility test splits into exactly these multiplicity patterns:

| Root multiplicities of \(h\) | Necessary divisibility of \(U\) |
| --- | --- |
| \(4\) | \((t-\alpha)^2\mid U\) |
| \(3+1\) | \((t-\alpha)^2(t-\beta)\mid U\) |
| \(2+2\) | \((t-\alpha)(t-\beta)\mid U\) |
| \(2+1+1\) | \((t-\alpha)(t-\beta)(t-\gamma)\mid U\) |
| \(1+1+1+1\) | \(h\mid U\) |

Here the displayed roots are distinct where their labels differ, and none is \(0\) because \(h(0)\ne0\). These are necessary conditions from one band, not yet exclusions of any pattern.

The GGV complete-chain data give the missing attachment of a repeated root to every standardized pair. In the published \(F_2\) row, \(A_1=(7/5,2)\); Theorem 2.20's type-II.b successor formula says there is a root \(\lambda\ne0\) of the first edge polynomial with multiplicity \(m_\lambda\), and
\[
A_1=(1,0)+\frac{m_\lambda}{3}(1/5,1).
\]
Thus \(m_\lambda=6\). In the ramified edge coordinate \(z=x^{1/5}y\), that edge polynomial is \(a h(z^5)^3\). Since \(z\mapsto z^5\) is unramified at \(\lambda\ne0\), the corresponding root \(\tau=\lambda^5\) of \(h(t)\) has multiplicity exactly two. A nonzero diagonal rescaling lets us set \(\tau=1\). Consequently every standardized pair has
\[
h(t)=(t-1)^2R(t),\qquad \deg R=2,\quad R(0)R(1)\ne0.
\]
This uses the published complete-chain successor condition, not a conjectural carrier ansatz. The quartic factor split is therefore only two cases:

| Branch | Quartic edge factor | Consequence of \(h\mid U^2\) |
| --- | --- | --- |
| Two factors, \(2+2\) | \((t-1)^2(t-\alpha)^2\) | \((t-1)(t-\alpha)\mid U\) |
| Three factors, \(2+1+1\) | \((t-1)^2(t-\alpha)(t-\beta)\) | \((t-1)(t-\alpha)(t-\beta)\mid U\) |

Here \(\alpha\ne1\), and in the second row \(\beta\) is distinct from both \(1\) and \(\alpha\); all roots are nonzero. The remaining lower-band equations still have to exclude both branches.

There is a useful rootwise refinement of the first complete-chain step. In the \(F_2\) case, the type-II.b root cutoff from GGV Proposition 5.16 / GGHV Proposition 2.5(4) is
\[
m_\lambda\ge \frac{m}{l}\frac{a\rho+bl\sigma}{\rho+\sigma}
=3\cdot\frac{5\cdot5+20\cdot(-1)}{5-1}=\frac{15}{4}.
\]
For \(h=(t-1)^2(t-\alpha)^2\), every root of the edge polynomial \(h(z^5)^3\) above either double root has multiplicity \(6\), so roots above both \(t=1\) and \(t=\alpha\) meet this cutoff. Applying the root-shear proposition at either choice gives the same immediate successor corner \(A_1=(7/5,2)\), because the successor depends on the multiplicity \(m_\lambda=6\). The downstream terminal chart is valid rootwise as well. After either shear, Proposition 5.18 gives an \((3,5)\)-pair in \(L^{(5)}\) with that regular corner. The continuation in the proof of GGHV Theorem 2.20 can be run from the chosen shear: remove any type-III corners by Proposition 5.17 and continue through type-II children until the terminal type-I corner. Each resulting chain is admissible. In the GGHV table, for initial edge \((A_0,A'_0)=((5,20),(1,0))\) and \((m,n)=(3,5)\), the only admissible row is \(F_2\), whose terminal corner is \(A_1=(7/5,2)\). Thus every qualifying root choice has the same terminal support normal, although the local centering coefficients are root-specific. For \(2+1+1\), only the distinguished double root qualifies; its simple roots have multiplicity 3, below the \(15/4\) cutoff.

For the second \(2+2\) chart, choose \(\theta^5=\alpha\) and normalize the ramified local coordinate as \(z=\theta(1+u)\). A monomial from band distance \(d\), with coefficient index \(i\), contributes a common nonzero scalar times
\[
\alpha^i(1+u)^{5i-15+d}\quad(P),\qquad
\alpha^i(1+u)^{5i-25+d}\quad(Q).
\]
Therefore the second-root terminal rows lie over \(\mathbb Q(\sqrt5)\), despite \(\theta\) itself having degree 10 over \(\mathbb Q\). The independent exact pointwise replay in `f2_75_125_alpha_terminal_chart_probe_20261018.py` verifies that one saved depth-six witness satisfying the \(t=1\) chart fails the \(t=\alpha\) chart already at forbidden position P, distance 1, \(u^4\). More generally, that first forbidden coefficient at distance 1 is a nonzero scalar multiple of \(S(\alpha)\), so the second chart requires \(S(\alpha)=0\). In the distance-seven gauged checkpoint, \(S\) depends only on \(p_1,p_2\), and exact evaluation gives
\[
S(\alpha)=\frac{-671-231\sqrt5}{758}p_1+
\frac{-2255-913\sqrt5}{1516}p_2.
\]
Together with \(S(1)=0\), an exact \(72\times64\) weight-49 kernel calculation and a rank certificate modulo 59 show the kernel dimension drops from 3 to 1: the modular ranks before adding root values, after \(S(1)=0\), and after both root values are 61, 62, and 63. Equivalently, \(p_2=((\sqrt5-5)/5)p_1\) on this necessary slice. This removes one free weight-49 direction; it is not a branch exclusion. Replays and receipts are `f2_75_125_dual_root_S_slice_20261018.py`, `f2_75_125_dual_root_weight49_rank_20261018.py`, `f2_75_125_dual_root_S_slice_receipt_20261018.json`, and `f2_75_125_dual_root_weight49_rank_receipt_20261018.json`.

## The near-edge recurrence

For \(1\le s\le4\), write the bands \(15-s\) and \(25-s\) as

\[
P_{15-s}=a\,t^3U_s,\qquad Q_{25-s}=b\,t^5V_s,
\]

where \(\deg U_s\le11\), \(\deg V_s\le19\), and set \(U_0=h^3,V_0=h^5\). The equation in output band \(z^{36-s}\), after dividing by \(ab\,t^7\), has the form

\[
s h^7D_{5/s}\!\left(3\frac{V_s}{h^5}-5\frac{U_s}{h^3}\right)+R_s=0,
\]

where

\[
D_{5/s}(w)=hw+th'w+\frac5s\,thw'
\]

and the already-known interaction term is

\[
R_s=
\sum_{a_0=1}^{s-1}
\left[
(-5a_0+3b_0)U_{a_0}V_{b_0}
+t\bigl((15-a_0)U_{a_0}V_{b_0}'-(25-b_0)U_{a_0}'V_{b_0}\bigr)
\right],
\qquad b_0=s-a_0.
\]

For \(s=1\), \(R_s=0\). For \(s=2\), the formula recovers the \(h\mid U_1^2\) condition above. The operator has the exact kernel test

\[
\bigl(t^s h^s w^5\bigr)'
=s\,t^{s-1}h^{s-1}w^4D_{5/s}(w).
\]

Because \(h(0)\ne0\), a rational solution of \(D_{5/s}(w)=0\) can be nonzero only when \(5\mid s\). Thus the first four near-edge operators are injective on rational functions; the first formal resonance is \(s=5\). At \(s=5\) the monomial support bounds also change (the lowest allowed \(t\)-powers drop), so the recurrence must be re-derived with the correct band supports there. This identifies a natural finite block to study next: compute the \(s=3,4\) compatibility conditions exactly, then treat the \(s=5\) support transition separately.

## The first four \(Q\)-bands are forced by a binomial series

Let \(w\) mark the distance from the top edge and set

\[
F(t,w)=\sum_{s=0}^4 U_s(t)w^s,\qquad
G(t,w)=\sum_{s=0}^4 V_s(t)w^s,
\quad U_0=h^3,\quad V_0=h^5.
\]

Up to the nonzero scalars \(a,b\), the near-edge parts are \(P=x^3F(t,w)\) and \(Q=x^5G(t,w)\), with \(t=xy^5\) and \(w=y\) in this truncated band chart. Direct differentiation shows that the coefficients of the nonconstant output bands are the coefficients of

\[
\Phi(F,G)
=w(-3FG_w+5GF_w)
+t(25F_tG-15FG_t)
+tw(-F_tG_w+F_wG_t).
\]

The series \(G=F^{5/3}\), expanded in \(\overline K(t)[[w]]\), makes \(\Phi\) vanish identically. For each \(1\le s\le4\), the difference between any other solution's \(V_s\) and this binomial coefficient lies in the kernel of \(D_{5/s}\). That kernel is zero, so the first four band equations force

\[
V_s=[w^s]F(t,w)^{5/3}\qquad (1\le s\le4).
\]

The obstruction is polynomiality of these rational binomial coefficients. The first three are

\[
\begin{aligned}
V_1={}&\frac53h^2U_1,\\
V_2={}&\frac53h^2U_2+\frac59\frac{U_1^2}{h},\\
V_3={}&\frac53h^2U_3+\frac{10}{9}\frac{U_1U_2}{h}
-\frac5{81}\frac{U_1^3}{h^4}.
\end{aligned}
\]

Thus \(V_3\in\overline K[t]\) exactly when

\[
h^4\mid 18h^3U_1U_2-U_1^3.
\]

If a root \(\alpha\) of \(h\) has multiplicity \(e\), this condition is equivalent to

\[
\operatorname{ord}_{\alpha}(U_1)\ge\left\lceil\frac{4e}{3}\right\rceil.
\]

Indeed, below that threshold the \(U_1^3\) term has order \(3\operatorname{ord}_\alpha(U_1)<4e\) and strictly lower order than the other numerator term, so it cannot cancel; at or above the threshold both terms are divisible by \((t-\alpha)^{4e}\). This sharpens the \(h\mid U_1^2\) condition from the previous band.

The fourth coefficient gives the next explicit compatibility condition:

\[
V_4=\frac53h^2U_4+
\frac{5}{243h^7}
\left(54h^6U_1U_3+27h^6U_2^2-9h^3U_1^2U_2+U_1^4\right).
\]

Hence the parenthesized numerator must be divisible by \(h^7\). This is now the first coupled condition on \(U_2,U_3\) after the root-contact restrictions on \(U_1\). It is a finite algebraic problem for a quartic \(h\), not yet an exclusion.

## The resonant fifth descent

The support bounds change exactly here. On bands \(10\) and \(20\),

\[
P_{10}=a t^2U_5,\qquad Q_{20}=b t^4V_5,
\qquad \deg U_5\le12,\quad \deg V_5\le20.
\]

Normalize the band series by \(a t^3\) and \(b t^5\). Thus their fifth coefficients are \(U_5/t\) and \(V_5/t\). With

\[
F=h^3+U_1w+U_2w^2+U_3w^3+U_4w^4+(U_5/t)w^5+O(w^6),
\]

the fifth-band equation has the binomial solution \([w^5]F^{5/3}\). Its homogeneous difference is governed by

\[
D_1(r)=hr+th'r+thr'=(thr)',
\]

whose rational kernel is \(r=\lambda/(th)\). Hence the actual coefficient is the binomial one plus this one-dimensional resonance. Expanding the coefficient gives the exact relation

\[
V_5=\frac53h^2U_5+\frac{\lambda}{3}h^4+
\frac{t}{729h^{10}}N_5,
\]

where

\[
N_5=810h^9(U_1U_4+U_2U_3)
-135h^6(U_1^2U_3+U_1U_2^2)
+60h^3U_1^3U_2-7U_1^5.
\]

Because \(h(0)\ne0\), \(\gcd(t,h)=1\). Since \(V_5,U_5,h\) are polynomials, the fifth band is solvable exactly when

\[
h^{10}\mid N_5.
\]

The resonance scalar \(\lambda\) does not affect this divisibility condition. This is the first coupled compatibility involving \(U_1,U_2,U_3,U_4\) after the earlier binomial constraints.

It also sharpens the local multiplicity bound. If \(\alpha\) is a root of \(h\) of multiplicity \(e\), and \(r=\operatorname{ord}_\alpha(U_1)<3e/2\), then \(-7U_1^5\) is the unique lowest-order term of \(N_5\) at \(\alpha\): every other summand has strictly larger order. This contradicts \(h^{10}\mid N_5\). Therefore

\[
\operatorname{ord}_\alpha(U_1)\ge\left\lceil\frac{3e}{2}\right\rceil.
\]

For the two factor patterns allowed by the GGV successor corner, the resulting minimum root factors of \(U_1\) are:

| Root multiplicities of \(h\) | Necessary factor of \(U_1\) from descent five |
| --- | --- |
| \(2+2\) | \((t-\alpha)^3(t-\beta)^3\) |
| \(2+1+1\) | \((t-\alpha)^3(t-\beta)^2(t-\gamma)^2\) |

The first four bands and the fifth descent together sharpen the contact at every double root further. Let \(u=t-\alpha\), where \(\operatorname{ord}_\alpha h=2\), and suppose \(\operatorname{ord}_\alpha U_1=3\). Write
\[
h=c u^2+O(u^3),\qquad U_1=A u^3+O(u^4),\qquad U_2=D+O(u),
\]
with \(cA\ne0\). Divisibility of the fourth binomial coefficient by \(h^7\) forces \(D\ne0\) and its order-12 coefficient gives
\[
27x^2-9x+1=0,\qquad x:=c^3D/A^2.
\]
The order-15 coefficient of \(N_5\) must also vanish because \(h^{10}\mid N_5\); it gives
\[
-135x^2+60x-7=0.
\]
The incompatibility has the exact Bézout certificate
\[
7=(70-225x)(27x^2-9x+1)+(9-45x)(-135x^2+60x-7).
\]
Therefore both equations cannot hold over any characteristic-zero field. Thus
\[
\operatorname{ord}_\alpha U_1\ge4
\]
at every double root of \(h\). This is an exact local exclusion of the contact-order-3 subbranch, not of either full factorization branch. Consequently a nonzero \(U_1\) has degree at least \(8\) on both allowed patterns: \(4+4\) for \(2+2\), or \(4+2+2\) for \(2+1+1\). Since \(\deg U_1\le11\), at most three residual degrees of freedom remain in its root factor. The remaining band equations must still eliminate these branches (and separately handle \(U_1=0\)).

This also makes the fourth-band condition especially simple. For the allowed root multiplicities, the fourth numerator
\[
N_4=54h^6U_1U_3+27h^6U_2^2-9h^3U_1^2U_2+U_1^4
\]
is divisible by \(h^7\) if and only if \(\operatorname{rad}(h)\mid U_2\). At a root of multiplicity \(e\in\{1,2\}\), if \(U_2\) does not vanish, the term \(27h^6U_2^2\) has order \(6e<7e\), while the other terms have order at least \(7e\); it cannot cancel. If \(U_2\) does vanish, all four terms have order at least \(7e\), using \(\operatorname{ord}(U_1)\ge2\) for \(e=1\) and \(\ge4\) for \(e=2\). Thus the initial bands reduce to the compact parameterization
\[
\begin{array}{c|c|c}
\text{root type of }h & U_1 & U_2\\ \hline
2+2 &(t-1)^4(t-\alpha)^4S_1 &(t-1)(t-\alpha)S_2\\
2+1+1 &(t-1)^4(t-\alpha)^2(t-\beta)^2S_1 &(t-1)(t-\alpha)(t-\beta)S_2
\end{array}
\]
with \(\deg S_1\le3\), \(\deg S_2\le9\) in the first row and \(\deg S_2\le8\) in the second. Here \(U_1=0\) is included by taking \(S_1=0\). These are exact consequences of the current band equations and the attached GGV double-root condition; they leave a small, explicit local algebra problem rather than an exclusion.

At a double root \(\alpha\), the fifth-band condition gives one more useful jet split. Write \(u=t-\alpha\), \(h=cu^2+O(u^3)\), and let \(r=\operatorname{ord}_\alpha U_1\ge4\), \(q=\operatorname{ord}_\alpha U_2\ge1\). If \(q=1\), then \(r=4\) is impossible: the unique order-18 term of \(N_5\) is \(-135h^6U_1U_2^2\). If \(r=5\), divisibility by \(h^{10}\) requires
\[
U_3(\alpha)=\frac{A D}{6c^3},\qquad
U_1=A u^5+O(u^6),\quad U_2=D u+O(u^2),
\]
from cancellation of the order-19 terms. If \(r\ge6\) (including \(U_1=0\)), it instead forces \(U_3(\alpha)=0\), since \(810h^9U_2U_3\) would otherwise be the unique order-19 term. For a simple root of \(h\), the already established orders \(\operatorname{ord} U_1\ge2\), \(\operatorname{ord}U_2\ge1\) make the fifth-band divisibility automatic locally. These are branch conditions for the next exact elimination, not a contradiction of the remaining global system.

## The sixth band raises the $U_2$ contact at double roots

The next band rules out the $q=1$ alternative above, without requiring a case split among $r=4$, $r=5$, and $r\ge6$. Work in the formal band variable $w$, and write
\[
F(t,w)=h^3+\sum_{s\ge1}F_s(t)w^s,\qquad
G(t,w)=h^5+\sum_{s\ge1}G_s(t)w^s,
\]
for the near-edge bands $0\le s\le9$, with $F_s=U_s,G_s=V_s$ for $s\le4$, and $F_s=U_s/t,G_s=V_s/t$ for $5\le s\le9$. Put $G=F^{5/3}R$, $S=th$, and $\rho=U_1/h^3$. The full Keller equation becomes
\[
w\left(3+t\frac{F_t}{F}\right)R_w
+t\left(15-w\frac{F_w}{F}\right)R_t
=\frac{w^{36}}{t^7F^{8/3}}.
\]
Thus the homogeneous recurrence displayed above is valid through the
coefficient of $w^{35}$; the constant Jacobian first appears at band
distance $36$.
The first four bands give $R_1=\cdots=R_4=0$. At the fifth and sixth powers of $w$, respectively, the coefficient equations give
\[
R_5=\frac{\kappa}{S},\qquad
6\left(3\frac{S'}h\right)R_6+5t\rho'R_5+15tR_6'-t\rho R_5'=0,
\]
so
\[
R_6=-\frac{\kappa\rho}{3S}.
\]
The only homogeneous ambiguity is $C S^{-6/5}$; it vanishes because $R_6$ is rational in $t$, whereas $S=th$ has a simple zero at $t=0$ ($h(0)\ne0$), making $S^{-6/5}$ have nonintegral valuation there. Consequently
\[
G_6-[w^6]F^{5/3}=\frac{4\kappa hU_1}{3t},\qquad
V_6=t[w^6]F^{5/3}+\frac{4\kappa hU_1}{3}.
\]

Now let $\alpha\ne0$ be a double root of $h$, and suppose $q=\operatorname{ord}_\alpha U_2=1$. We already know $r=\operatorname{ord}_\alpha U_1\ge4$; the remaining $U_3,\ldots,U_6$ are polynomials and have nonnegative order. In
\[
[w^6]F^{5/3}=h^5\sum_{\sum j m_j=6}
\binom{5/3}{k}\frac{k!}{\prod_jm_j!}
\prod_{j=1}^6\left(\frac{F_j}{h^3}\right)^{m_j},
\qquad k=\sum_jm_j,
\]
the unique term of least order is the partition (6=2+2+2):
\[
-\frac5{81}\frac{U_2^3}{h^4},\qquad
\operatorname{ord}_\alpha=-5.
\]
Every other partition has order greater than $-5$. The correction $4\kappa hU_1/(3t)$ is regular at $\alpha$, so it cannot cancel this pole; this contradicts polynomiality of $V_6$. The exact partition and valuation census is replayed by [f2_75_125_k6_local_jet_certificate.py](../src/f2_75_125_k6_local_jet_certificate.py). Thus the sixth band strengthens the local condition to
\[
\operatorname{ord}_\alpha U_2\ge2
\]
at every double root of $h$.

For the two allowed root types, this updates the necessary factorization of $U_2$ to
\[
\begin{array}{c|c|c}
\text{root type of }h&\text{necessary factor of }U_2&\text{residual degree bound}\\\hline
2+2&(t-1)^2(t-\alpha)^2&\deg S_2\le7\\
2+1+1&(t-1)^2(t-\alpha)(t-\beta)&\deg S_2\le7
\end{array}
\]
where $U_2$ is that factor times $S_2$. This removes all fifth-band $q=1$ subbranches, including the tuned $r=5$ case, but it is still only a local necessary condition and does not exclude either global factorization branch.

## The seventh band classifies the remaining double-root jets

The sixth-band formula extends one more step. Set \(\sigma=U_2/h^3\). The coefficient of $w^7$ in the ratio equation is
\[
7\left(3\frac{S'}h\right)R_7+6t\rho'R_6
+5t(\sigma'-\rho\rho')R_5+15tR_7'
-t\rho R_6'+t(\rho^2-2\sigma)R_5'=0.
\]
Substitution of the already determined $R_5,R_6$ gives
\[
R_7=\frac{\kappa}{S}\left(-\frac{\sigma}{3}+\frac{2\rho^2}{9}\right).
\]
The homogeneous solution would be $C S^{-7/5}$, which cannot be rational at the simple zero $t=0$ of $S=th$; hence $C=0$. Extracting the seventh coefficient of $G=F^{5/3}R$ gives
\[
G_7-[w^7]F^{5/3}
=\frac{\kappa}{t}\left(\frac43hU_2+\frac{2U_1^2}{9h^2}\right),
\]
or, since $G_7=V_7/t$,
\[
V_7=t[w^7]F^{5/3}+\kappa\left(\frac43hU_2+\frac{2U_1^2}{9h^2}\right).
\]
At a double root, the correction is regular under the already known bounds \(\operatorname{ord}U_1\ge4\), \(\operatorname{ord}U_2\ge2\); at a simple root it is regular under \(\operatorname{ord}U_1\ge2\), \(\operatorname{ord}U_2\ge1\). Thus polynomiality of $V_7$ forces the binomial coefficient \([w^7]F^{5/3}\) to have no pole there.

This yields a sharper exhaustive local split at any double root \(\alpha\). Write $u=t-\alpha$, $h=cu^2+O(u^3)$, $U_1=A u^r+O(u^{r+1})$, and $U_2=B u^q+O(u^{q+1})$, with $c\ne0$. The case $q=1$ was already eliminated by the sixth band. For $r=4,q=2$, put
\[
X=A/c^2,\qquad Y=B/c,\qquad Z=U_3(\alpha).
\]
The order-(-2) term of the sixth coefficient and the order-(-4) term of the seventh coefficient give
\[
\begin{aligned}
K_6={}&729Z^2-81Y^3-486XYZ+162X^2Y^2+108X^3Z-63X^4Y+7X^6=0,\\
K_7={}&729Y^2Z+729XZ^2-324XY^3-972X^2YZ
+378X^3Y^2+189X^4Z-126X^5Y+13X^7=0.
\end{aligned}
\]
These are weighted homogeneous of weights $(1,2,3)$ on $(X,Y,Z)$. Since $X\ne0$, scale to $X=1$. The exact resultant in $Z$ is
\[
\operatorname{Res}_Z(K_6(1,Y,Z),K_7(1,Y,Z))
=-3^{15}(3Y-1)^7.
\]
Therefore $Y=1/3$, and both specialized equations become $(27Z-1)^2=0$. In the unscaled coordinates this forces
\[
Y=X^2/3,\qquad Z=X^3/27.
\]

The other contact patterns are excluded by two short leading-term comparisons. If $r=4,q\ge3$, the sixth and seventh residues, after scaling $X=1$, are respectively
\[
729Z^2+108Z+7=0,\qquad 729Z^2+189Z+13=0;
\]
their difference gives $Z=-2/27$, which leaves $3=0$ in the first equation. If $r\ge5,q=2$, the sixth residue gives $9Z^2-Y^3=0$, while the seventh gives $Y^2Z=0$; $q=2$ means $Y\ne0$, so these are incompatible. Finally, if $r\ge5,q\ge3$, the unique lowest term of the sixth coefficient is $5U_3(\alpha)^2/(9h)$, forcing $U_3(\alpha)=0$. The exhaustive weighted-partition censuses, normalized equations, and exact resultant are replayed by [f2_75_125_k7_double_root_certificate.py](../src/f2_75_125_k7_double_root_certificate.py).

Thus every double root lies either on the exceptional tuned jet $r=4,q=2$, with the two displayed coefficient relations, or on a deeper-contact branch $r\ge5,q\ge3$, where $U_3(\alpha)=0$.

The identical classification holds at a simple root, with the minimal pair of orders $(r,q)=(2,1)$. Indeed, for $e=1$ the K6/K7 valuation minima and normalized residue equations are the same as the double-root calculation after dividing all root orders and pole orders by two. Hence at a simple root the mixed cases $r=2,q\ge2$ and $r\ge3,q=1$ are impossible; the surviving local alternatives are $r=2,q=1$ with $Y=X^2/3$, $U_3(\alpha)=X^3/27$, or $r\ge3,q\ge2$ with $U_3(\alpha)=0$. Here again $X=(U_1/h^2)(\alpha)$ and $Y=(U_2/h)(\alpha)$ in the tuned case.

In particular, for every root of $h$ of multiplicity $e\in\{1,2\}$, all possibilities left by these two bands lie in one of the patterns
\[
(r,q)=(2e,e)\text{ with }Y=X^2/3,\ U_3(\alpha)=X^3/27,
\qquad\text{or}\qquad
r\ge2e+1,\ q\ge e+1,\ U_3(\alpha)=0.
\]
These are still necessary local patterns, not a contradiction of the global factorization branches or the full coefficient system.

Because these alternatives hold at every root of $h$, they have a concise global form. The root orders imply
\[
U_1=h^2A,\qquad U_2=hB,
\qquad \deg A\le3,\quad \deg B\le7,
\]
using \(\deg U_1,\deg U_2\le11\). At a tuned root, $B=A^2/3$ and $U_3=A^3/27$ after evaluation; at a deeper root all three values $A,B,U_3$ vanish. Therefore
\[
\operatorname{rad}(h)\mid 3B-A^2,
\qquad
\operatorname{rad}(h)\mid 27U_3-A^3.
\]
This packages all K6/K7 root conditions into two exact divisibilities and gives a smaller global parameterization for either quartic root type. It also shows that the fifth-band numerator condition is then automatic:
\[
\frac{N_5}{h^{10}}=
810hAU_4+810BU_3-135A^2U_3-135hAB^2+60A^3B-7A^5,
\]
which is polynomial. So the next new information must come from at least the eighth band or from other rows of the full outer system.

There is also a compact exact form of the K6/K7 polynomiality conditions, including their higher root multiplicities. After substituting $U_1=h^2A$, $U_2=hB$, write $C=U_3$. The terms with negative powers of $h$ in the two binomial coefficients are
\[
[w^6]F^{5/3}=\frac{W_6}{h}+\text{terms regular at the roots of }h,
\qquad
[w^7]F^{5/3}=\frac{W_7}{h^2}+\frac{W_{7,1}}h+\text{terms regular at the roots of }h,
\]
where
\[
\begin{aligned}
W_6={}&\frac59C^2-\frac5{81}B^3-\frac{10}{27}ABC+\frac{10}{81}A^2B^2
+\frac{20}{243}A^3C-\frac{35}{729}A^4B+\frac{35}{6561}A^6,\\
W_7={}&-\frac5{27}B^2C-\frac5{27}AC^2+\frac{20}{243}AB^3+\frac{20}{81}A^2BC
-\frac{70}{729}A^3B^2-\frac{35}{729}A^4C+\frac{70}{2187}A^5B-\frac{65}{19683}A^7,\\
W_{7,1}={}&\frac{10}{9}CU_4-\frac{10}{27}ABU_4+\frac{20}{243}A^3U_4.
\end{aligned}
\]
The omitted terms may contain (t^{-1}), because the fifth and later bands are divided by (t); (t) is a unit at every root of (h), and the factor (t) in (V_6,V_7) clears those denominators. Therefore polynomiality requires \(h\mid\widehat W_6\) and \(h^2\mid\widehat W_7+h\widehat W_{7,1}\). Set $D=3B-A^2$, $E=27C-A^3$, and \(\Delta=E-3AD\), and rescale the pole numerators by \(\widehat W_6=(6561/5)W_6\), \(\widehat W_7=(-19683/5)W_7\), and \(\widehat W_{7,1}=(-19683/5)W_{7,1}\). Direct substitution simplifies these to
\[
\widehat W_6=\Delta^2-3D^3,
\qquad
\widehat W_7=A\widehat W_6+3D^2\Delta,
\qquad
\widehat W_{7,1}=-162U_4\Delta.
\]
This explains the leading root relations above: both $D$ and $E$ vanish at every root. Hence \(\operatorname{rad}(h)\mid D,E\), so \(\Delta\) also vanishes at every root. It follows that \(\widehat W_6=\Delta^2-3D^3\) has order at least two at every root and therefore the full K6 condition \(h\mid\widehat W_6\) follows. At a simple root, the three terms in \(\widehat W_7+h\widehat W_{7,1}\) have order at least two, so the full K7 condition is automatic there. At a tuned double root, the full divisibility \(h^2\mid\widehat W_7+h\widehat W_{7,1}\) further forces
\[
\operatorname{ord}_\alpha(D)\ge2,
\qquad
\operatorname{ord}_\alpha(E)\ge2.
\]
At a deeper double root, write $A=a_1u+\cdots$, $D=d_1u+\cdots$, and \(\Delta=e_1u+\cdots\), with $h=cu^2+\cdots$. The order-three coefficient of the same divisibility is
\[
e_1\left(a_1e_1+3d_1^2-162cU_4(\alpha)\right)=0.
\]
The exact pole-part expansion and the change of variables $(A,B,C)\mapsto(A,D,E)$ are checked in [f2_75_125_k7_double_root_certificate.py](../src/f2_75_125_k7_double_root_certificate.py). These conditions retain possible branches; they do not yet force the full system to be inconsistent.

## Eighth-band checks on the local branches

The K8 recurrence can also be solved explicitly. Write \(F=h^3 f\), \(f=1+\rho w+\sigma w^2+\tau w^3+\cdots\), and set $D=t\,d/dt$, $a=D\log S$, and $L_j=[w^j]\log f$. The coefficient recurrence gives
\[
R_8=\frac{\kappa}{S}\left(-\frac{14}{81}\rho^3+\frac49\rho\sigma-\frac13\tau\right).
\]
Its homogeneous ambiguity $C S^{-8/5}$ is again excluded by rationality at $t=0$. Coefficient extraction then yields
\[
G_8-[w^8]F^{5/3}
=\frac{4\kappa}{81t}\left(27hU_3+\frac{9U_1U_2}{h^2}-\frac{U_1^3}{h^5}\right),
\]
so the correction in the polynomial band $V_8$ is
\[
V_8=t[w^8]F^{5/3}+\frac{4\kappa}{81}\left(27hU_3+\frac{9U_1U_2}{h^2}-\frac{U_1^3}{h^5}\right).
\]
Under either K6/K7 root pattern this correction is regular at every root of $h$. Thus any pole of the binomial coefficient must vanish on its own.

There is an exact global form of this K8 pole condition. Put $C=U_3$ and retain
\[
D=3B-A^2,\qquad E=27C-A^3,\qquad \Delta=E-3AD.
\]
After $U_1=h^2A$, $U_2=hB$, the resonance correction simplifies to
\[
\frac{4\kappa}{81}h(27C+9AB-A^3),
\]
so it contributes no pole along $h=0$. The singular part of the binomial coefficient is
\[
[w^8]F^{5/3}=\frac{N_8}{h^3}+\text{terms regular at every root of }h,
\]
where the exact numerator is
\[
\boxed{\displaystyle
N_8=\frac5{59049}\left[
\bigl(A\Delta-81U_4h\bigr)^2
-D\Delta^2+6AD^2\Delta-3A^2D^3+3D^4
-243D^2U_4h+486\Delta U_5h^2
\right].}
\]
The numerator is obtained independently from all weighted partitions of $8$ in the binomial expansion and is then checked against this compressed expression in [f2_75_125_k8_local_jet_certificate.py](../src/f2_75_125_k8_local_jet_certificate.py). Since $h(0)\ne0$, polynomiality of $V_8$ forces the exact divisibility $h^3\mid N_8$.

This expression sharpens the root-jet split. Write $J_8=A\Delta-81U_4h$ and let $e=\operatorname{ord}_\alpha h\in\{1,2\}$. At a tuned root, the K6/K7 conditions give $\operatorname{ord}_\alpha D,\operatorname{ord}_\alpha\Delta\ge e$ (the double-root statement uses the full K7 divisibility). Every term after $J_8^2$ then has order at least $3e$. Thus $h^3\mid N_8$ forces
\[
\operatorname{ord}_\alpha(J_8)\ge\left\lceil\frac{3e}{2}\right\rceil,
\]
namely order at least two at a simple tuned root and at least three at a double tuned root. At a deeper simple root, $\operatorname{ord}_\alpha A,\operatorname{ord}_\alpha D,\operatorname{ord}_\alpha\Delta\ge1$. If $U_4(\alpha)\ne0$, then $J_8$ has order one, whose square cannot be canceled by the remainder, which has order at least three. Hence $U_4(\alpha)=0$. At a deeper double root, the only possible order-three term in the remainder is the leading part of $-D\Delta^2$; $J_8^2$ has order at least four. Therefore its order-three coefficient must vanish, forcing
\[
\operatorname{ord}_\alpha D\ge2
\quad\hbox{or}\quad
\operatorname{ord}_\alpha\Delta\ge2.
\]
In particular, the earlier excluded subcase $q=3$, $\operatorname{ord}_\alpha(U_3)=1$ is recovered: it has both orders equal to one. These are exact necessary constraints on the surviving branches, not yet an exclusion of the tuned branches or of the full degree pair.

An explicit all-deep assignment shows that these local divisibility tests, by themselves, cannot finish the exclusion. For any quartic $h$ with $h(0)\ne0$, set
\[
U_1=0,\qquad U_2=h^2,\qquad U_3=U_4=h,\qquad U_5=0.
\]
At a root of multiplicity $e$, this gives $\operatorname{ord}(U_2)=2e\ge e+1$ and $\operatorname{ord}(U_3)=e>0$, so it lies on the deeper-contact alternative. Direct substitution into the displayed K4--K8 conditions gives
\[
\begin{aligned}
N_4&=27h^{10},&
N_5/h^{10}&=810h^2,\\
\widehat W_6&=729h^2-81h^3,&
\widehat W_7+h\widehat W_{7,1}&=-3645h^3,\\
N_8&=\frac5{59049}\left(-2187h^3+4617h^4\right).
\end{aligned}
\]
Here $h^7\mid N_4$, $N_5$ is $h^{10}$ times the displayed polynomial quotient, $h\mid\widehat W_6$, $h^2\mid\widehat W_7+h\widehat W_{7,1}$, and $h^3\mid N_8$. The exact rational checker also verifies this substitution. It is a survivor of the local pole and divisibility conditions for either quartic root type; it is **not** a solution of the full band equations or a constructed Keller pair. This pins down why the next attack must use cross-band equations, degree bounds, or a global compatibility beyond the local root tests.

The exact K8 partition census gives these additional local facts:

* At a tuned root, $r=2e,q=e$ and $Y=X^2/3$, $U_3(\alpha)=X^3/27$. The lowest possible pole residue of \([w^8]F^{5/3}\) cancels exactly. This leading K8 check therefore leaves the tuned stratum alive; higher Laurent coefficients may still constrain it.
* At a deeper simple root, $r\ge3,q\ge2$, $U_3(\alpha)=0$, the unique lowest term is \(5U_4(\alpha)^2/(9h)\). Polynomiality forces $U_4(\alpha)=0$.
* At a deeper double root, $r\ge5,q\ge3$, $U_3(\alpha)=0$, the subcase $q=3,\operatorname{ord}_\alpha(U_3)=1$ has a unique uncancellable order-$-3$ term, \(-5U_2U_3^2/(27h^4)\). Hence this subcase is impossible; a surviving deeper double root must have $q\ge4$ or \(\operatorname{ord}_\alpha(U_3)\ge2\). Further K8 residues in those surviving subcases have not yet been eliminated.

The exact recurrence coefficient arithmetic and these finite valuation censuses are replayed by [f2_75_125_k8_local_jet_certificate.py](../src/f2_75_125_k8_local_jet_certificate.py). This gives local pruning for the deeper-contact branches, while confirming that the same leading-residue strategy does not remove the tuned branches.

### Ninth-band pruning on the $U_1=0$ deeper-double-root branch

There is a further exact condition at $w$-degree nine, still below the
inhomogeneous source at degree 36. Let $\alpha\ne0$ be a double root of $h$,
put $u=t-\alpha$, and assume the deeper branch $U_1=0$. Write
\[
q=\operatorname{ord}_\alpha U_2,\quad
r=\operatorname{ord}_\alpha U_3,\quad
s=\operatorname{ord}_\alpha U_4,\quad
v=\operatorname{ord}_\alpha(U_5/t).
\]
The K7/K8 conditions give $q\ge3$, $r\ge1$, $s\ge1$, and when $q=3$ they
require $r\ge2$. For
\[
F=h^3+U_2w^2+U_3w^3+U_4w^4+(U_5/t)w^5+\cdots,
\]
the degree-nine coefficient of the binomial part is the sum over the eight
partitions
\[
9,\ 7+2,\ 6+3,\ 5+4,\ 5+2+2,\ 4+3+2,\ 3+3+3,\ 3+2+2+2.
\]
Writing $F_j=[w^j]F$, the exact coefficient is
\[
\begin{aligned}
[w^9]F^{5/3}={}&\frac53h^2F_9
+\frac{10}{9h}(F_7F_2+F_6F_3+F_5F_4)
-\frac5{27h^4}F_5F_2^2
-\frac{10}{27h^4}F_4F_3F_2\\
&-\frac5{81h^4}F_3^3
+\frac{20}{243h^7}F_3F_2^3.
\end{aligned}
\]
At a double root, the potentially singular contributions from
$5+4$, $5+2+2$, $4+3+2$, $3+3+3$, and $3+2+2+2$ have
orders respectively
\[
s+v-2,\quad v+2q-8,\quad s+r+q-8,\quad 3r-8,\quad r+3q-14
\]
The homogeneous first-integral expression derived below shows that the only
resonance reaching degree nine is the K5 mode. Its contribution is
\[
\frac{\kappa_1}{t}[w^4]F^{4/3}
=\frac{\kappa_1}{t}\left(\frac43hU_4+\frac29\frac{U_2^2}{h^2}\right).
\]
The coefficients of
the two terms that can tie on the surviving $q=3$ strata are exactly
\[
\frac{20}{243}\frac{U_3U_2^3}{h^7}
\quad\hbox{and}\quad
-\frac5{27}\frac{(U_5/t)U_2^2}{h^4}.
\]

Since $q\ge3$, $s\ge1$, and $t$ is a unit at $\alpha$, this correction has
order at least two. It cannot cancel any negative-order term below. The exact
root-order census therefore yields these necessary K9 conditions:

* If $q=3$, the case $r=2$ is impossible. If $r=3$, then $v=0$; if $r=4$,
  then $v=1$. In either tied case, writing
  $h=cu^2+\cdots$, $U_2=du^3+\cdots$, $U_3=eu^r+\cdots$, and
  $U_5/t=gu^{r-3}+\cdots$, cancellation forces
  \[
  4de=9c^3g.
  \]
  If $r\ge5$ (including $U_3=0$), then $v\ge2$.
* If $q\ge4$, the cases $r=1,2$ are impossible. For $r\ge3$, one must have
  $s+v\ge2$.

The coefficient identities and all bounded local-order branches through the
support caps are replayed exactly by
[f2_75_125_k9_zero_first_slice.py](../src/f2_75_125_k9_zero_first_slice.py).
The checker verifies 19,992 order strata. These conditions further prune the
$U_1=0$ deeper-double-root branch, but leave the tied coefficient loci and
the higher-order cases alive; they do not exclude either quartic root pattern
or the degree pair.

## First integral for the homogeneous band block

The homogeneous part of the ratio PDE has a closed characteristic invariant.
It organizes the top band block without solving each nonresonant recurrence
separately. Write
\[
F=h^3f,\qquad S=th,\qquad
a=D_t\log S,\quad L=D_t\log f,\quad M=w\,\frac{f_w}{f},
\qquad D_t=t\frac{d}{dt}.
\]
The homogeneous operator is
\[
\mathcal L(R):=w(3a+L)R_w+(15-M)D_tR.
\]
The actual Keller equation has $\mathcal L(R)=w^{36}/(t^7F^{8/3})$,
so the top block through $w^{35}$ satisfies $\mathcal L(R)=0$.
Define
\[
I=\frac{w^5}{S f^{1/3}}.
\]
The exact derivative identities $wI_w=I(5-M/3)$ and
$D_tI=I(-a-L/3)$ give
\[
\frac{\mathcal L(I)}I
=(3a+L)(5-M/3)+(15-M)(-a-L/3)=0.
\]
Therefore every formal series $R=\Phi(I)$, with $\Phi(0)=1$, solves the
homogeneous equation. Conversely, the coefficient recurrence at $w^s$ has
homogeneous part
\[
\left(D_t+\frac{s}{5}a\right)R_s=0,
\]
whose solutions are scalar multiples of $S^{-s/5}$. Since $S$ has order one
at $t=0$ and $R_s\in K(t)$, a homogeneous ambiguity is rational exactly
when $5\mid s$. At each such resonance there is one new scalar. Thus the
actual top block is organized by
\[
\Phi(z)=1+\kappa_1z+\kappa_2z^2+\cdots,\qquad
R\equiv\Phi\!\left(\frac{w^5}{S f^{1/3}}\right)\pmod{w^{36}}.
\]
If $C=F^{1/3}=h f^{1/3}$, this becomes
\[
\boxed{\displaystyle
G\equiv C^5\Phi\!\left(\frac{w^5}{tC}\right)
=\sum_{r\ge0}\kappa_r\,t^{-r}w^{5r}C^{5-r}
\pmod{w^{36}},\qquad \kappa_0=1.}
\]
The K6--K8 corrections are the first coefficients of the $\kappa_1$ mode.
At band $10$, a second rational resonance $\kappa_2$ is available; formulas
that retain only $\kappa_1$ cannot be extrapolated past that band.

The outer support envelope gives exact $t=0$ pole-depth bounds for the
coefficients of these series:
\[
\operatorname{pole}_0(F_s)\le
\min\!\left(\left\lfloor\frac{s}{5}\right\rfloor,3\right),\qquad
\operatorname{pole}_0(G_s)\le
\min\!\left(\left\lfloor\frac{s}{5}\right\rfloor,5\right).
\]
Here $\operatorname{pole}_0$ is the allowed pole order at $t=0$; the
[support audit](../src/f2_75_125_support_audit.py)
checks these bounds for every band. At a resonance $s=5r$, the new $\kappa_r$
mode begins in $G_s$ with $t^{-r}h^{5-r}$. This fits the $Q$ support for
$r\le5$. The first mode whose own pole exceeds the support cap is $\kappa_6$
at band distance $30$, where its $t^{-6}$ term must cancel against
contributions from earlier modes and the $C^5$ term. That is a concrete
support-compatibility equation in the homogeneous top block.

### Tenth-band pruning on the $U_1=0$ double-root branch

At band distance ten, the first integral gives exactly
\[
G_{10}=[w^{10}]F^{5/3}
+\frac{\kappa_1}{t}[w^5]F^{4/3}
+\frac{\kappa_2}{t^2}h^3.
\]
Assume $U_1=0$ and let $\alpha\ne0$ be a double root of $h$. Put
\[
q=\operatorname{ord}_\alpha U_2,\quad
r=\operatorname{ord}_\alpha U_3,\quad
s=\operatorname{ord}_\alpha U_4,\quad
v=\operatorname{ord}_\alpha(U_5/t).
\]
The surviving K9 strata have either $q=3$ with $(r,v)=(3,0)$,
$(4,1)$, or $r\ge5,v\ge2$, or they have $q\ge4,r\ge3,s+v\ge2$.
Since $t$ is a unit at $\alpha$,
\[
[w^5]F^{4/3}=\frac43hF_5+\frac49\frac{F_2F_3}{h^2}
\]
has order at least two on every such stratum; the $\kappa_2$ term has order
six. Therefore polynomiality of $G_{10}$ forces the binomial coefficient
\([w^{10}]F^{5/3}\) itself to be regular at $\alpha$.

The exact weighted-partition census has twelve partitions of ten using
parts $2,\ldots,10$. For $q=3$, if $s\ge1$, its unique lowest term has
order $-5$ and is
\[
-\frac7{729}\frac{U_2^5}{h^{10}},
\]
which is impossible. Hence $s=0$. At $s=0$, the three terms of order
$-5$ are exactly
\[
-\frac7{729}\frac{U_2^5}{h^{10}}
+\frac{20}{243}\frac{U_2^3U_4}{h^7}
-\frac5{27}\frac{U_2U_4^2}{h^4}.
\]
Write $h=cu^2+\cdots$, $U_2=du^3+\cdots$, and $U_4=a+\cdots$ at
$u=t-\alpha$, with $c,d\ne0$. Cancellation of the $u^{-5}$ coefficient
requires
\[
135c^6a^2-60c^3d^2a+7d^4=0.
\]
Equivalently, $z=c^3a/d^2$ satisfies
\[
135z^2-60z+7=0,
\]
whose discriminant is $-180$ and whose two roots over $\mathbb C$ are
$(10\pm i\sqrt5)/45$. Thus K10 forces $U_4(\alpha)\ne0$ and pins its
normalized leading value to two possibilities. It does not close the branch;
the K9 relation $4de=9c^3g$ on the $r=3,4$ subbranches involves other
leading coefficients and is not contradicted by this quadratic.

For $q\ge4$, K9 plus K10 imply two further restrictions. If $v=0$, K9
requires $s\ge2$, but the unique lowest K10 pole is then the nonzero
$(5+5)$ term $(5/9)(U_5/t)^2/h$, of order $-2$; therefore $v\ne0$.
If $q=4$ and $s=0$ or $1$, K9 requires respectively $v\ge2$ or $v\ge1$.
The unique lowest K10 term is then $(2+4+4)$, namely
$-(5/27)U_2U_4^2/h^4$, of order $-4$ or $-2$. Hence $q=4$ forces
$s\ge2$. These are robust unique-pole obstructions: the terms involved use
only already-known slices, while all $F_6,\ldots,F_{10}$ were allowed their
worst-case order zero in the exhaustive census.

The exact coefficient and order checks are replayed by
[f2_75_125_k10_zero_first_slice.py](../src/f2_75_125_k10_zero_first_slice.py).
K10 advances local pruning on the $U_1=0$ branch, but the $q=3$ quadratic
has solutions over the algebraic closure and higher $q$ strata remain. It
does not exclude either quartic root pattern or the degree pair.

### Eleventh-band elimination of the two K9 ties

The next coefficient improves the $q=3$ double-root branch. At a double
root $\alpha$, put $u=t-\alpha$ and use the K10 conclusions
\[
h=cu^2+\cdots,\qquad U_2=du^3+\cdots,\qquad
U_4=a+\cdots,\qquad z=\frac{c^3a}{d^2},\qquad
135z^2-60z+7=0.
\]
K10 has already forced $s=\operatorname{ord}_\alpha U_4=0$. Write
$U_3=eu^r+\cdots$ and $F_5=U_5/t=gu^v+\cdots$ when these slices are
nonzero. K9 leaves $r=3,v=0$ with $4de=9c^3g$, $r=4,v=1$ with the same
relation, or $r\ge5,v\ge2$.

Since $U_1=0$, the first-integral formula gives
\[
G_{11}=[w^{11}]F^{5/3}+\frac{\kappa_1}{t}[w^6]F^{4/3}.
\]
On the $q=3,s=0$ stratum the resonant correction has order at least $-1$.
The K9 tied cases have binomial poles of order $-5$ and $-4$, respectively,
so their leading residues cannot be changed by that correction. Exact
enumeration of all fourteen weighted partitions of $11$ shows that in both
cases the entire lowest-order residue comes from
\[
-\frac5{27}\frac{U_3U_4^2}{h^4}
-\frac{10}{27}\frac{U_2U_4F_5}{h^4}
+\frac{20}{81}\frac{U_2^2U_3U_4}{h^7}
+\frac{20}{243}\frac{U_2^3F_5}{h^7}
-\frac{35}{729}\frac{U_2^4U_3}{h^{10}}.
\]
Normalize $\gamma=c^3g/(de)$. After factoring out the nonzero
$d^4e/c^{10}$, this residue is
\[
A(z)+\gamma B(z),\qquad
A(z)=-\frac5{27}z^2+\frac{20}{81}z-\frac{35}{729},\quad
B(z)=-\frac{10}{27}z+\frac{20}{243}.
\]
K9 fixes $\gamma=4/9$ on both tied cases. The K10 equation then reduces the
residue to
\[
\frac{-405z^2+180z-25}{2187}=-\frac4{2187}\ne0.
\]
This contradiction eliminates both local strata $(r,v)=(3,0)$ and $(4,1)$.

The same K11 census further sharpens the $r\ge5$ alternatives. If $r=5$,
then $v\ge3$ leaves a unique order-$-3$ residue $A(z)$; K10 reduces it to
$4(30z-7)/729$, which is nonzero at either K10 root. Thus $r=5$ forces
$v=2$. At $v=2$, K11 instead imposes the solvable tie
\[
\gamma=\frac{2(30z-7)}{15(9z-2)}.
\]
If $r\ge6$ (including $U_3=0$), $v=2$ leaves the unique order-$-3$
residue $B(z)$, which is nonzero at the K10 roots; hence $v=2$ is excluded.
Through K11, the remaining $q=3$ cases have $r=5,v=2$, or
$r\ge6,v\ge3$ (including zero-slice cases compatible with the same order
bound). The next band tests both residual strata.

The exact partition coefficients, local order minima, and residue
reductions are replayed by
[f2_75_125_k11_q3_double_root.py](../src/f2_75_125_k11_q3_double_root.py).
The checker uses rational arithmetic throughout. It certifies this local
K11 pruning only; it does not eliminate the $q\ge4$ double-root strata, the
simple-root strata, or the degree pair.

### Twelfth-band exclusion of the full $q=3$ double-root branch

At distance twelve the first-integral formula is
\[
G_{12}=[w^{12}]F^{5/3}
+\frac{\kappa_1}{t}[w^7]F^{4/3}
+\frac{\kappa_2}{t^2}[w^2]F.
\]
On the K11 survivors $q=3,s=0$, namely $r=5,v=2$ or $r\ge6,v\ge3$,
exact enumeration of all twenty-one weighted partitions of $12$ gives a
unique lowest order $-8$ residue. Its four terms involve only $U_2$ and
$U_4$:
\[
-\frac5{81}\frac{U_4^3}{h^4}
+\frac{10}{81}\frac{U_2^2U_4^2}{h^7}
-\frac{35}{729}\frac{U_2^4U_4}{h^{10}}
+\frac{35}{6561}\frac{U_2^6}{h^{13}}.
\]
The $\kappa_1$ term has order at least one on every residual stratum; the
$\kappa_2$ term has order three. Neither can cancel this pole. With
$h=cu^2+\cdots$, $U_2=du^3+\cdots$, and $U_4=a+\cdots$, the residue is
\[
\frac{a^3}{c^4}\left(-\frac5{81}+\frac{10}{81z}
-\frac{35}{729z^2}+\frac{35}{6561z^3}\right),
\qquad z=\frac{c^3a}{d^2}.
\]
After clearing $6561z^3$, the numerator is
$-405z^3+810z^2-315z+35$. Modulo the K10 relation
$135z^2-60z+7=0$, this reduces to $\frac73(1-6z)$. The K10 polynomial
does not vanish at $z=1/6$ (its value is $3/4$), and $a,c,d$ are nonzero.
So the K12 residue cannot vanish. This excludes every $q=3$ stratum at the
double root on the $U_1=0$ branch.

The exact partition census, correction-order bounds, and polynomial
reduction are replayed by
[f2_75_125_k12_q3_double_root.py](../src/f2_75_125_k12_q3_double_root.py).
It checks 78 bounded/zero-slice order strata. The same exact K12 census
also closes the remaining low-$s$ part of the $q\ge5$ branch: if $s=0$,
the unique order-$-8$ term is $-(5/81)U_4^3/h^4$; if $s=1$, the unique
order-$-5$ term is the same. For these cases the $\kappa_1$ correction has
order at least $-1$ or $0$, respectively, and the $\kappa_2$ term has
order at least five. Thus $q\ge5$ forces $s\ge2$. Since K10 already gave
$q=4\Rightarrow s\ge2$, the entire surviving double-root branch now has
\[
q=\operatorname{ord}_\alpha U_2\ge4,\qquad
s=\operatorname{ord}_\alpha U_4\ge2.
\]
This still leaves the $q\ge4,s\ge2$ double-root branch and the simple-root
branch; neither the quartic carrier nor the degree pair is excluded.

### Three-band exclusion for $q=4,s=2,v=1,r\ge4$

One part of the remaining double-root branch closes by combining K11--K13.
Let $\alpha$ be a double root, $u=t-\alpha$, and assume
\[
q=\operatorname{ord}_\alpha U_2=4,\quad
s=\operatorname{ord}_\alpha U_4=2,\quad
v=\operatorname{ord}_\alpha(U_5/t)=1,\quad
r=\operatorname{ord}_\alpha U_3\ge4.
\]
Here $r=\infty$ means $U_3=0$. Writing
\[
h=cu^2+\cdots,\quad U_2=du^4+\cdots,\quad U_4=au^2+\cdots,
\quad U_5/t=gu+\cdots,\quad U_6=f+\cdots,
\]
the K11, K12, and K13 binomial parts have lowest pole orders $-1,-2,-3$.
Their exact lowest terms, respectively, are
\[
\frac{10}{9}\frac{F_5F_6}{h}
-\frac{10}{27}\frac{U_2U_4F_5}{h^4}
+\frac{20}{243}\frac{U_2^3F_5}{h^7},
\]
\[
\frac59\frac{F_6^2}{h}
-\frac5{81}\frac{U_4^3}{h^4}
-\frac5{27}\frac{U_2F_5^2}{h^4}
-\frac{10}{27}\frac{U_2U_4F_6}{h^4}
+\frac{10}{81}\frac{U_2^2U_4^2}{h^7}
+\frac{20}{243}\frac{U_2^3F_6}{h^7}
-\frac{35}{729}\frac{U_2^4U_4}{h^{10}}
+\frac{35}{6561}\frac{U_2^6}{h^{13}},
\]
and
\[
-\frac5{27}\frac{U_4^2F_5}{h^4}
-\frac{10}{27}\frac{U_2F_5F_6}{h^4}
+\frac{20}{81}\frac{U_2^2U_4F_5}{h^7}
-\frac{35}{729}\frac{U_2^4F_5}{h^{10}}.
\]
The respective $\kappa_1$ resonance corrections have orders at least
$2,1,0$, so none can affect those residues; the $\kappa_2$ modes are
$U_1=0$, $U_2$, and $U_3$, also of higher order.

Set $z=c^3a/d^2$ and $\phi=c^6f/d^3$. K11 reduces to
\[
\phi=\frac{9z-2}{27}.
\]
After substituting this into K13, its residue is a nonzero scalar multiple
of $-(3z-1)^2$, so $z=1/3$. Finally, set $W=c^9g^2/d^5$. The K12 residue,
after K11, is a nonzero scalar multiple of
\[
81W+(3z-1)^3.
\]
At $z=1/3$ this forces $W=0$, hence $g=0$, contradicting $v=1$.
Therefore the local stratum $q=4,s=2,v=1,r\ge4$ is empty.

The partition census and all three exact residue reductions are replayed by
[f2_75_125_k11_k13_q4_v1_double_root.py](../src/f2_75_125_k11_k13_q4_v1_double_root.py).
The remaining $r=3$ case is excluded at K15 below. The still-open cases here
include $q=4,s\ge3$, $q\ge5,s\ge2$, and the simple-root branch.

### K15 excludes the remaining $q=4,s=2,v=1,r=3$ jet

Now take the same double-root normalization, but set
\[
q=4,\qquad s=2,\qquad v=1,\qquad r=3.
\]
Write
\[
h=cu^2+\cdots,\quad U_2=du^4+\cdots,\quad U_3=eu^3+\cdots,
\quad U_4=au^2+\cdots,\quad F_5:=U_5/t=gu+\cdots,
\quad U_6=f+\cdots,
\]
where $c,d,e,a,g$ are nonzero. Use the four nonzero normalized parameters
\[
z=\frac{c^3a}{d^2},\qquad
\gamma=\frac{c^3g}{de},\qquad
\epsilon=\frac{e^2c^3}{d^3},\qquad
\phi=\frac{c^6f}{d^3}.
\]
In particular, $\epsilon\ne0$ on this order stratum.

The K15 first-integral formula is
\[
G_{15}=[w^{15}]F^{5/3}
+\frac{\kappa_1}{t}[w^{10}]F^{4/3}
+\frac{\kappa_2}{t^2}[w^5]F
+\frac{\kappa_3}{t^3}h^2.
\]
At the double root the binomial term has lowest order $-5$, contributed by
exactly 20 weighted partitions of 15 using parts $2,\ldots,6$. The three
resonance corrections have orders at least $-2$, $1$, and $4$, respectively;
$t$ is a unit there. Thus they cannot alter the order-$-5$ residue. After
dividing by the nonzero common factor $d^6e/c^{16}$, its equation is
\[
\begin{aligned}
E_{15}:={}&-189\epsilon^2-1215\epsilon\gamma^3
+4860\epsilon\gamma^2+4860\epsilon\gamma z-5670\epsilon\gamma
+1620\epsilon\phi-3780\epsilon z+2100\epsilon\\
&-7290\gamma\phi z+4860\gamma\phi+4860\gamma z^2
-3780\gamma z+630\gamma-3645\phi^2\\
&+9720\phi z-3780\phi+1620z^3-5670z^2+3150z-455=0.
\end{aligned}
\]
The earlier K11--K14 lowest-residue equations are $E_{11}=E_{12}=E_{13}
=E_{14}=0$ in the same four parameters. A degree-three exact cofactor
identity gives
\[
\epsilon^3=Q_{11}E_{11}+Q_{12}E_{12}+Q_{13}E_{13}+Q_{14}E_{14}+Q_{15}E_{15},
\]
where all five rational polynomials $Q_i$ are stored as sparse exact
coefficients in
[the certificate JSON](../results/f2_75_125_k15_q4_v1_r3_certificate.json).
The checker regenerates the 20-partition residue, confirms that the K11--K14
residues agree with their displayed normalizations, and expands this identity
over $\mathbb Q$ exactly:
[f2_75_125_k15_q4_v1_r3_certificate.py](../src/f2_75_125_k15_q4_v1_r3_certificate.py).
Because $\epsilon\ne0$, the identity contradicts the five band equations.
This excludes the last $r=3$ subbranch. Combined with the preceding K11--K13
result and K10's exclusion of $v=0$, the stratum $q=4,s=2$ requires finite
$v\ge2$ (or $F_5=0$, which is outside this notation's finite-order case).

### K11--K15 exclude the $r=3$, $v\ge2$ continuation

Now take $q=4$, $s=2$, and $r=3$, but assume only
$\operatorname{ord}_\alpha F_5\ge2$ (including $F_5=0$). The worst case for
the order census is $v=2$. Write
\[
h=cu^2+\cdots,\quad U_2=du^4+\cdots,\quad U_3=eu^3+\cdots,
\quad U_4=au^2+\cdots,\quad F_6=f+\cdots,
\]
with $c,d,e,a\ne0$, and set
\[
z=\frac{c^3a}{d^2},\qquad
\epsilon=\frac{e^2c^3}{d^3},\qquad
\phi=\frac{c^6f}{d^3}.
\]
The K11--K15 binomial residues have orders $-1,-2,-3,-4,-5$ and come from
5, 11, 7, 13, and 11 lowest-order partitions, respectively. None contains
part 5, so the leading coefficient of $F_5$ does not enter. The respective
$(\kappa_1,\kappa_2,\kappa_3)$ correction order lower bounds are
$(2,\infty,\infty)$, $(1,4,\infty)$, $(0,3,\infty)$,
$(-1,2,\infty)$, and $(-2,2,4)$; every correction is strictly above its
binomial residue. A degree-two identity over $\mathbb Q[z,\epsilon,\phi]$
then gives
\[
\epsilon=Q_{11}A_{11}+Q_{12}A_{12}+Q_{13}A_{13}+Q_{14}A_{14}
+Q_{15}A_{15},
\]
where $A_{11},\ldots,A_{15}$ are the normalized lowest-residue equations
and the exact sparse cofactors are in
[the certificate JSON](../results/f2_75_125_k15_q4_v2_r3_certificate.json).
The replay script regenerates all five residues, checks the correction
orders, and verifies the identity identically over $\mathbb Q$:
[f2_75_125_k15_q4_v2_r3_certificate.py](../src/f2_75_125_k15_q4_v2_r3_certificate.py).
Since $\epsilon\ne0$, this excludes the entire $r=3$, $v\ge2$ continuation,
including $F_5=0$. Together with K10 for $v=0$ and the previous K15 result
for $v=1$, the local $q=4,s=2,r=3$ branch is now closed. This does not close
$q=4,s=2,r\ge4$, $q=4,s\ge3$, $q\ge5,s\ge2$, or the simple-root branch.

At a simple root of the $2+1+1$ pattern, the K10 local census for the
$U_1=0$ branch is as follows. Let $\alpha$ be a simple root of
$h$, set $u=t-\alpha$, and write
\[
q=\operatorname{ord}_\alpha U_2\ge2,\quad
r=\operatorname{ord}_\alpha U_3\ge1,\quad
s=\operatorname{ord}_\alpha U_4\ge1,\quad
v=\operatorname{ord}_\alpha(U_5/t)\ge0.
\]
The bounds $q\ge2$, $r\ge1$ come from the deeper K6/K7 branch; K8 gives
$s\ge1$. At a simple root the $\kappa_1$ correction in $G_{10}$ has order
at least one, and the $\kappa_2$ term has order three. The only potentially
singular terms of $[w^{10}]F^{5/3}$ are therefore
\[
\frac{10}{81}\frac{U_2^2U_3^2}{h^7}
-\frac{10}{27}\frac{U_2U_3F_5}{h^4}
-\frac5{27}\frac{U_3^2U_4}{h^4}
+\frac59\frac{F_5^2}{h},
\qquad F_5=U_5/t.
\]
Their orders are respectively $-7+2q+2r$, $-4+q+r+v$,
$-4+2r+s$, and $-1+2v$. Each is at least $-1$, so K10 polynomiality
requires only the cancellation of whichever terms attain order $-1$.
The exact exhaustive census of 22,176 bounded order strata gives these
surviving local alternatives:

| Orders at the simple root | Required K10 condition |
| --- | --- |
| $q=2,r=1,s=1,v\ge1$ | $2d^2-3c^3f=0$ |
| $q=2,r=1,v=0,s\ge2$ | $2d^2e^2-6c^3deg+9c^6g^2=0$ |
| $q=2,r=1,s=1,v=0$ | $2d^2e^2-6c^3deg-3c^3e^2f+9c^6g^2=0$ |
| $q\ge3,r=1,s=1,v=0$ | $e^2f-3c^3g^2=0$ |
| $q\ge3,r=1,s\ge2,v\ge1$ | No K10 pole at this root |
| $r\ge2,v\ge1$ | No K10 pole at this root |

Here $h=cu+\cdots$, $U_2=du^q+\cdots$, $U_3=eu^r+\cdots$,
$U_4=fu^s+\cdots$, and $F_5=gu^v+\cdots$ whenever the relevant order is
finite; coefficients appearing in a row are nonzero. Any unlisted
combination has a unique nonzero order-$-1$ term and is impossible. In
particular, if $v=0$ and $r\ge2$, K10 excludes the stratum. The four tie
equations above still have solutions over the algebraic closure, so this is
branch pruning, not an exclusion of the $2+1+1$ type. The checker also
replays every exact partition coefficient and verifies that no other
partition can contribute a pole on the stated order domain.

### Simple-root K9--K10 cross-band divisibilities on $U_1=0$

At each simple root of $S(t)=t^2-3t+3$, the exact K9 pole rules out
$r=\operatorname{ord}(U_3)=1$: the unique leading term is
$-(5/81)U_3^3/h^4$ of order $-1$. Hence $\operatorname{ord}(U_3)\ge2$.
With this bound, K10 rules out $v=\operatorname{ord}(U_5/t)=0$: the unique
leading pole is $(5/9)(U_5/t)^2/h$. The K9 and K10 correction terms have
strictly larger local orders, so neither can cancel the indicated pole.
Since the roots of $S$ are distinct and simple, these facts imply
\[
S^2\mid U_3,\qquad S\mid (U_5/t).
\]
The earlier simple-root bounds give $S^2\mid U_2$ and $S\mid U_4$; the
distinguished double-root conditions give
$\operatorname{ord}_{t=1}U_2,\operatorname{ord}_{t=1}U_3\ge3$. Therefore,
for $T=(t-1)^3S^2$,
\[
T\mid U_2,\qquad T\mid U_3,\qquad S\mid U_4,
\qquad S\mid(U_5/t),
\]
and $\deg(U_2/T),\deg(U_3/T)\le4$ because $\deg U_2,\deg U_3\le11$.
This is an exact necessary divisibility reduction for the d19 coefficient
system, not a $2+1+1$ branch exclusion. The bounded simple-root profile
screen adds no further eliminations beyond these order implications through
K125.

### Exact $t=0$ face at distances 30--35

The first forbidden-pole block can be read particularly cleanly after the
ramified substitution $t=\tau^5$, $w=\tau v$. A monomial $t^e w^s$ has
$\tau$-order $5e+s$. Thus, for $s=30,\ldots,34$, the only possible pole
deeper than the $Q$ cap of five is $t^{-6}$, detected by
\[
 [\tau^{s-30}v^s]\,G(\tau^5,\tau v)=0.
\]
At $s=35$, both $t^{-7}$ and $t^{-6}$ exceed the cap, giving the two
conditions $[v^{35}]G(0,v)=0$ and
$[\tau^5v^{35}]G(\tau^5,\tau v)=0$. There are seven forbidden-pole
coefficients in this block, not just the six depth-six coefficients.

The first-integral formula becomes
\[
G(\tau^5,\tau v)=\sum_{r=0}^7\kappa_rv^{5r}
  F(\tau^5,\tau v)^{(5-r)/3}\pmod{w^{36}}.
\]
These equations do not by themselves contradict the support envelope. For
an explicit exact face witness, take either allowed quartic root pattern
\[
h_{22}(t)=(t^2-1)^2,
\qquad
h_{211}(t)=-\frac54(t-1)^2(t-2)(t+\tfrac25).
\]
Both have $h(0)=1$, $h'(0)=0$, and a double root at $1$; the first has
root type $2+2$, the second $2+1+1$. Set
\[
F(t,w)=h(t)^3+t^{-1}w^5.
\]
Then $F(\tau^5,\tau v)=1+\xi+O(\tau^{10})$, where $\xi=v^5$.
Normalize the removable shear modes $\kappa_2=\kappa_5=0$, and set
$\kappa_1=\kappa_3=\kappa_4=0$. The two remaining over-cap face
coefficients vanish for
\[
\kappa_6=-\binom{5/3}{6}=-\frac{35}{6561},
\qquad
\kappa_7=-\binom{5/3}{7}
  -\kappa_6\binom{-1/3}{1}=\frac{10}{6561}.
\]
Indeed, these choices make the $\xi^6$ and $\xi^7$ coefficients of
\[
(1+\xi)^{5/3}+\kappa_6\xi^6(1+\xi)^{-1/3}
  +\kappa_7\xi^7(1+\xi)^{-2/3}
\]
zero. The $\tau^1,\ldots,\tau^5$ face coefficients also vanish, so all
seven forbidden-pole equations at distances 30--35 pass in both root
patterns. The arithmetic and the two quartic factorizations are replayed by
[f2_75_125_t0_overflow_witness.py](../src/f2_75_125_t0_overflow_witness.py).

This is intentionally only a $t=0$ face witness. It does not establish that
the resulting $Q$ coefficients have no poles at the nonzero roots of $h$,
does not solve the inhomogeneous distance-36 equation, and is not a Keller
pair. Its consequence for the search is narrower but decisive: the first
over-cap resonance block cannot exclude either GGV root pattern by itself.
The next useful step has to couple this face to root polynomiality or to the
distance-36 source equation.

### First face of the inhomogeneous distance-36 equation

The same ramified chart gives a compact necessary equation at the first
source band. Write $z=v^5$ and expand
\[
\widehat F(\tau,v):=F(\tau^5,\tau v)
 =A(z)+\tau vB(z)+O(\tau^2),\qquad
\widehat G(\tau,v):=G(\tau^5,\tau v)
 =G_0(z)+\tau vC(z)+O(\tau^2).
\]
The exact support envelope gives
\[
\deg A\le3,\quad \deg B\le3,\quad
\deg G_0\le5,\quad \deg C\le5,
\]
with $A(0)=h(0)^3\ne0$. Let $R=G/F^{5/3}$, and set
$D(z)=15A(z)-5zA'(z)$ and $E(z)=zG_0'(z)-5G_0(z)$. Substituting
$t=\tau^5,w=\tau v$ into the full ratio PDE and taking its $\tau^1$
coefficient gives the exact polynomial identity
\[
\boxed{D(z)C(z)+5B(z)E(z)=5z^7.}
\]
Here $\deg D\le2$ and $\deg E\le4$: the nominal top terms cancel in
$15A-5zA'$ and $zG_0'-5G_0$. The coefficient of $z^7$ is
\[
a_2c_5-b_3g_4=1,
\]
where $a_2,b_3,c_5,g_4$ are the indicated coefficients of $A,B,C,G_0$.
Under the chart $P=x^3F$, $Q=x^5G$, these are respectively the linear
coefficients of $x$ in $P$, $y$ in $P$, $y$ in $Q$, and $x$ in $Q$; this
equation is exactly the determinant-one condition on the linear part. So
the first source face recovers a required Jacobian normalization rather
than an obstruction.

There is an exact face jet compatible with both quartic root patterns:
\[
A=1,\quad B=z^3,\quad G_0=1-z^4,\quad C=\frac53z^3.
\]
It satisfies $DC+5BE=5z^7$. Its homogeneous first-integral parameters are
$\kappa_0=1,\kappa_4=-1$ on the $\tau^0$ face. The predicted $\tau^1$
homogeneous series is $\frac53z^3-\frac13z^7$; the latter term is at
distance 36, where the source equation replaces the homogeneous recurrence.
The chosen support-compatible face coefficient $C$ agrees through distance
35, and the exact distance-36 identity supplies the source term. If lifted,
this face would have linear parts $P_{\rm lin}=y$ and $Q_{\rm lin}=-x$,
whose determinant is $1$. The polynomial identity and
this substitution are checked by
[f2_75_125_distance36_face.py](../src/f2_75_125_distance36_face.py).

This entire first-face identity is just the restriction of the original
Jacobian equation to $y=0$. In fact, set
\[
p_0(x)=x^3A(1/x),\quad p_1(x)=x^3B(1/x),\quad
q_0(x)=x^5G_0(1/x),\quad q_1(x)=x^5C(1/x).
\]
Then $DC+5BE=5z^7$ is equivalent to
\[
p_0'(x)q_1(x)-p_1(x)q_0'(x)=1.
\]
Since $\deg p_0=3$ and $\deg q_0=5$, the derivatives have degrees two
and four. Whenever $p_0'$ and $q_0'$ are coprime, the extended Euclidean
algorithm supplies $q_1,p_1$ with degree bounds $\deg q_1<4$ and
$\deg p_1<2$, which are within the F2 support limits. Thus this face is
generically solvable; it cannot be the rigidity block by itself. The
worked face has $p_0=x^3$, $p_1=1$, $q_0=x^5-x$, and
$q_1=5x^2/3$, giving the identity exactly.

### The second inhomogeneous $\tau$ face also survives

Write the next coefficients as
\[
\widehat F=A(z)+\tau vB(z)+\tau^2v^2F_2(z)+O(\tau^3),\qquad
\widehat G=G_0(z)+\tau vC(z)+\tau^2v^2G_2(z)+O(\tau^3),
\]
and let $R=\widehat G/\widehat F^{5/3}$. The change of variables
$(t,w)=(\tau^5,\tau v)$ transforms the exact ratio PDE into
\[
\left(15-\frac{v\widehat F_v}{\widehat F}\right)R_\tau
+v\frac{\widehat F_\tau}{\widehat F}R_v
=5v^{36}\widehat F^{-8/3}.
\]
For the compatible first face above, set $A=1$, $B=z^3$,
$G_0=1-z^4$, and $C=5z^3/3$. If
\[
r_1=C-\frac53G_0B,\qquad
r_2=G_2-\frac53G_0F_2-\frac53CB+\frac{20}{9}G_0B^2,
\]
the coefficient of $\tau$ in this PDE, after removing the common factor
$v^2$, is
\[
30r_2-16z^3r_1+B(r_1+5zr_1')
+5z(2F_2-B^2)G_0'=-\frac{40}{3}z^{10}.
\]
The exact choice
\[
F_2=\frac53z^2,\qquad G_2=\frac{25}{9}z^2
\]
satisfies this identity; it also has $B(0)=0$, the visible $U_1=0$
face condition. The corresponding axis two-jet is
\[
P=x^3+y+\frac53xy^2+O(y^3),\qquad
Q=x^5-x+\frac53x^2y+\frac{25}{9}x^3y^2+O(y^3),
\]
and its Jacobian is $1+O(y^2)$. The exact rational replay is
[f2_75_125_distance36_tau2_face.js](../src/f2_75_125_distance36_tau2_face.js).

This verifies compatibility through the second $\tau$ face at $t=0$ only.
It neither lifts the data at the nonzero roots of $h$ nor closes the later
faces under the finite-support bounds. The obstruction search must couple
those cross-root conditions to the complete coefficient system; the
branchwide exclusion remains open.

The next coefficient uses $F_3(z)$ and $G_3(z)$. Their degree caps are not
extra ansatz assumptions: they follow directly from the exact outer support
envelopes. A supported monomial $x^iy^j$ in $P=x^3F(t,w)$ has
$t$-exponent $i-3$ and $w$-exponent $j-5(i-3)$, so its $\tau$-order after
$t=\tau^5,w=\tau v$ is exactly $j$. Hence the $\tau^3$ face comes only from
the $j=3$ support slice. The $P$ edge bound $5i-j\le15$ gives
$0\le i\le3$, producing $z$-degrees $3-i=0,1,2,3$; for $Q=x^5G$,
$5i-j\le25$ gives $0\le i\le5$, producing degrees $5-i=0,1,\ldots,5$.
Thus $\deg F_3\le3$ and $\deg G_3\le5$ are exact support-derived bounds.
The replay in
[f2_75_125_distance36_tau3_face.py](../src/f2_75_125_distance36_tau3_face.py)
now enumerates these slices before solving the exact $\mathbb Q$ face system.
Its compatible two-parameter witness is therefore valid for this entire
$\tau^3$ coefficient space, while still being only a $t=0$ jet.

For the particularly small witness $a=b=0$, the next $\tau^3$ coefficient
of the transformed ratio PDE also has an exact support-compatible solution.
At $\tau^4$, the same projection of the $j=4$ slices gives
$\deg F_4\le3$ and $\deg G_4\le5$; the exact linear system has 8 equations
in 10 coefficients, rank 8, and a two-parameter affine solution. The replay
is [f2_75_125_distance36_tau3_coefficient.py](../src/f2_75_125_distance36_tau3_coefficient.py).
This extends one local jet by another source-face coefficient and confirms
that the $t=0$ slice alone is not currently producing an obstruction; it
still gives no lift at the nonzero roots of $h$.

In the original coordinates, the $\kappa_r$ kernel mode adds
$\kappa_rP^{(5-r)/3}$ to $Q$. The cases $r=2$ and $r=5$ are the polynomial
target shears $\kappa_2P$ and $\kappa_5$; they can be normalized away while
preserving the Keller condition and degree pair. Since $S_P\subseteq S_Q$,
these shears also preserve the exhaustive outer support envelope. The other modes have
fractional or negative powers and cannot be removed by a polynomial target
shear. The first-integral identity is checked in
[f2_75_125_k8_local_jet_certificate.py](../src/f2_75_125_k8_local_jet_certificate.py).
The remaining task is to impose the exact Laurent support and degree
restrictions, then join this homogeneous block to the inhomogeneous row at
distance $36$.

## Front-end reduction and finite exhaustive outer envelope

### Normalizing the nonzero Jacobian constant

The standard-pair theorems are stated for a Jacobian in \(\overline K^\times\), whereas the coefficient system below uses \([P,Q]=1\). After standardization, write \([P,Q]=\delta\ne0\) and replace \(Q\) by \(\delta^{-1}Q\). Bilinearity gives \([P,\delta^{-1}Q]=1\). This target scaling preserves the degree and support of \(Q\), all Newton corners and edge directions, and the subalgebra \(\overline K[P,Q]\); it also preserves whether the map is an automorphism. Therefore the normalized pair is still a standard counterexample with the same F2 data, and setting the band-system right-hand side to \(1\) loses no cases.

The starting point is not an assumption on the literal coordinates. The upstream Lemma A extraction proves, for an arbitrary Keller pair of degrees \((75,125)\), an automorphic standardization over \(\overline K\) to a standard \((3,5)\)-pair with transformed degree gcd \(g\le25\). Its key degree-transport identity is \(25=M\bar a+N\bar b\), while the standard corner has \(g=\bar a+\bar b\). GGV5's complete-chain theorem and its census through \(g\le25\) leave only \(F_2,j=1\), with \(A_0=(5,20)\) and \(g=25\). Thus the degrees of the standardized pair are again exactly \((75,125)\), and \(M=N=1\) in the transport identity. This deduction uses the published normal-form and GGV chain inputs; the upstream bundle is a written proof extraction, not a proof-kernel formalization.

Here is the support-transfer argument, with the two geometric steps kept separate. Before the final standardization, Lemma A's subrectangularizer gives \(\operatorname{Supp}(\bar P)\subseteq[0,a]\times[0,b]\), with corner \((a,b)=m(\bar a,\bar b)\). Its positive-weight argument proves that \((n\bar a,n\bar b)\) is the unique maximizer of every positive weight on \(\operatorname{Supp}(\bar Q)\). Choosing weights \((T,1)\) and \((1,T)\), with \(T\) larger than the opposite corner coordinate, rules out any support exponent beyond either coordinate; hence \(\operatorname{Supp}(\bar Q)\subseteq[0,n\bar a]\times[0,n\bar b]\). The F2 census gives \((\bar a,\bar b)=(5,20)\), while the degree-transport identity is \(25=5M+20N\) with positive integers \(M,N\), forcing \(M=N=1\). The final standardizing map in GGV1 Proposition 5.20 is \(x\mapsto x,\ y\mapsto y+\lambda\); it preserves these rectangular bounds because it only lowers the \(y\)-exponent.

For the present problem over \(\mathbb C\), the subrectangularization input also has a direct published source, so this step need not rest on an unverified reading of van den Essen Corollary 10.2.21. Makar-Limanov's 2025 paper puts an arbitrary complex counterexample into a reduced coordinate form by polynomial automorphisms; it then states that for every positive weight the leading form of the first component is the single upper-right monomial, and concludes that its Newton polygon lies in the rectangle with that corner (pp. 304–305). If the original Jacobian is \(c\ne0\), replace \(Q\) by \(c^{-1}Q\) first; this leaves \(\mathbb C[P,Q]\) unchanged. This supplies the needed rectangle and corner for \(P\) over \(\mathbb C\). The bound for \(Q\), the standard-pair reduction, and the F2 edge half-spaces still come from the separate GGV arguments above; Makar-Limanov's result alone does not establish the 2,162-variable model.

The stricter half-space bounds below are not transported through that translation. In fact, \(x^{15}y^{60}\) becomes \(x^{15}(y+\lambda)^{60}\), which contains \(x^{15}\) when \(\lambda\ne0\); the latter exponent has \(5i-j=75\), outside the proposed \(P\)-half-space \(5i-j\le15\). Instead, apply GGHV Theorem 2.20 to the already standardized pair. Its complete chain starts from that standard pair, and its initial regular corner is a corner of that pair itself. The \(F_2,j=1\) table row gives \(A'_0=(1,0)\), \(A_0=(5,20)\), and \(m=3,n=5\). Thus the \(P\)-edge has vertices \(mA'_0=(3,0)\) and \(mA_0=(15,60)\); its supporting line is \(5i-j=15\), and all of \(\operatorname{Supp}(P)\) lies on the polygon side \(5i-j\le15\). The corresponding \(Q\)-edge is scaled by \(n/m\), with vertices \((5,0)\) and \((25,100)\), so all of \(\operatorname{Supp}(Q)\) lies in \(5i-j\le25\). This uses the F2 corner as a property of the final standard pair, not as an edge property claimed to survive the shear.

The extraction's subrectangular witness and its positive-weight support argument give

\[
\operatorname{Supp}(P)\subseteq[0,15]\times[0,60],\qquad
\operatorname{Supp}(Q)\subseteq[0,25]\times[0,100].
\]

After the translation, the rectangular support bounds are therefore

\[
\operatorname{Supp}(P)\subseteq[0,15]\times[0,60],\qquad
\operatorname{Supp}(Q)\subseteq[0,25]\times[0,100].
\]

The post-standardization \(F_2\) corner calculation above gives the additional half-space bounds. Taken together, these imply that every monomial of the standardized pair lies in

\[
\begin{aligned}
S_P&=\{(i,j)\in\mathbb Z^2:0\le i\le15,\ 0\le j\le60,\ 5i-j\le15\},\\
S_Q&=\{(i,j)\in\mathbb Z^2:0\le i\le25,\ 0\le j\le100,\ 5i-j\le25\}.
\end{aligned}
\]

This is an exhaustive **outer** support model for every hypothetical \((75,125)\) pair after the front-end automorphism: it may contain monomials that cannot coexist in an actual Keller pair, but it contains every monomial allowed by the standard rectangles and the forced corner edge. Therefore, proving its coefficient system inconsistent would suffice for exclusion. It avoids needing to attach the narrower approximate-root/Puiseux carrier model; the tradeoff is a much larger coefficient system.

The arithmetic and the shear distinction are replayed by [f2_75_125_edge_bridge_audit.py](../src/f2_75_125_edge_bridge_audit.py). It checks the corner normal, endpoint support values, exact support counts and band ranges, and the explicit failure of the translation to preserve the stricter half-space. This checker verifies the finite arithmetic conditional on the cited standard-pair and F2-corner theorems; it does not replace those theorem inputs.

Set \(t=xy^5\), \(z=y^{-1}\), so \(x^iy^j=t^iz^{5i-j}\). The band index is \(k=5i-j\). For the two components the exact ranges are

\[
-60\le k\le15\quad(P),\qquad -100\le l\le25\quad(Q).
\]

Within a band, the possible \(t\)-exponents are

\[
\begin{aligned}
i_P^-(k)&=\max(0,\lceil k/5\rceil),&
i_P^+(k)&=\min(15,\lfloor(k+60)/5\rfloor),\\
i_Q^-(l)&=\max(0,\lceil l/5\rceil),&
i_Q^+(l)&=\min(25,\lfloor(l+100)/5\rfloor).
\end{aligned}
\]

Thus \(P_k(t)=\sum_{i=i_P^-(k)}^{i_P^+(k)}a_{i,5i-k}t^i\), and similarly for \(Q_l\). The support counts are independently checkable from

\[
|S_P|=4\cdot61+\sum_{i=4}^{15}(76-5i)=586,
\qquad
|S_Q|=6\cdot101+\sum_{i=6}^{25}(126-5i)=1576.
\]

In particular, \(P_{15}=t^3f(t)\) with \(\deg f=12\), and \(Q_{25}=t^5g(t)\) with \(\deg g=20\); the two endpoints of each forced edge give nonzero constant and leading coefficients of \(f,g\). The small pure-Python enumerator [f2_75_125_support_audit.py](../src/f2_75_125_support_audit.py) checks the lattice counts, endpoint inclusion, and band intervals without external dependencies.

For \(P=\sum_k p_k(t)z^k\) and \(Q=\sum_l q_l(t)z^l\), the Jacobian equation becomes the finite family

\[
E_r(t):=\sum_{k+l=r+4}\bigl(kP_kQ_l'-lP_k'Q_l\bigr)
=\begin{cases}1,&r=0,\\0,&r\ne0.\end{cases}
\]

This gives \(586+1576=2162\) coefficient variables before endpoint nonvanishing and normalizations. A raw Gröbner run in all variables is not the useful first calculation; eliminate bands from the top, beginning with the quartic edge factorization derived above.

### Exact sparse-system replay

The exact sparse system is mechanically enumerable. [f2_75_125_exact_band_system.py](../src/f2_75_125_exact_band_system.py) loops over every P-support point and every Q-support point. For monomials `x^i y^j` and `x^u y^v`, it checks termwise the coefficient identity `i*v - j*u = (5*i-j)*u - (5*u-v)*i` and places the result in output row `(z-band, t-degree) = (5*i-j+5*u-v-4, i+u-1)`. The replay gives **3,790 nonzero coefficient rows** and **912,713 nonzero bilinear terms** from the **923,536** possible P-Q monomial pairs; the other 10,823 pairs have zero bracket coefficient. Every nonzero row has at least two terms, so the raw system has no immediate one-term contradiction. The unique constant-Jacobian right-hand side is row `(0,0)`. The top output band `z^36` contains 31 rows, at t-degrees 8 through 38, with 268 nonzero terms. The possible t-degree-39 corner interaction vanishes exactly because `15*100 - 60*25 = 0`. The deterministic bilinear-term stream has SHA-256 `2e36049f836f07e455fa076279a3fa41eef70d9f30ddc48ce2d03965ac19117e`.

This replay establishes the coefficient transcription once the support envelope is accepted; it does not prove the cited standardization or solve the system. Its value is that all later band elimination can be audited against one fixed, exact equation stream rather than a manually selected subset.

## The first 60 y-adic slices are flexible

The same finite support system can be read in the original \(y\)-adic direction. Write

\[
P=\sum_{j=0}^{60}p_j(x)y^j,\qquad Q=\sum_{j=0}^{100}q_j(x)y^j.
\]

The outer support gives, for \(1\le j\le60\),

\[
\deg p_j\le d_j:=3+\left\lfloor\frac j5\right\rfloor,
\qquad
\deg q_j\le d_j+2.
\]

The \(y^n\)-coefficient of \([P,Q]=1\) is

\[
(n+1)(p_0'q_{n+1}-p_{n+1}q_0')+T_n=\delta_{n0},
\]

where

\[
T_n=\sum_{i=1}^{n}(n-i+1)p_i'q_{n-i+1}
-\sum_{i=0}^{n-1}(i+1)p_{i+1}q_{n-i}'.
\]

For an actual pair, the \(n=0\) equation implies \(\gcd(p_0',q_0')=1\); the forced endpoints give \(\deg p_0'=2\) and \(\deg q_0'=4\). Fix all earlier layers and set \(d=d_{n+1}\). The new layer enters through

\[
L_d(u,v)=p_0'v-uq_0',\qquad \deg u\le d,\quad\deg v\le d+2.
\]

This map is onto the polynomials of degree at most \(d+4\). Indeed, its kernel consists exactly of \((u,v)=(p_0'r,q_0'r)\) with \(\deg r\le d-2\), so its kernel dimension is \(d-1\); the domain and target dimensions are \(2d+4\) and \(d+5\), respectively. Also \(\deg T_n\le d+4\), by the floor-sum inequality \(\lfloor a/5\rfloor+\lfloor b/5\rfloor\le\lfloor(a+b)/5\rfloor\). Thus every non-edge slice \(0\le n<60\) can be solved once the earlier layers are fixed.

There is one compatibility to check when \(j=n+1=5k\): the top coefficients of \(p_j\) and \(q_j\) are already specified by the edge polynomials \(P_{15}=t^3f(t)\), \(Q_{25}=t^5g(t)\). The coefficient of \(x^{d+4}y^{5k-1}\) is a top-output coefficient, so support forces every contributing monomial pair to lie on both top edges. Its vanishing is precisely the corresponding coefficient of the already imposed edge identity \(3fg'-5f'g=0\). The kernel projects onto the one-dimensional space of changes in the two leading coefficients that preserve this top-output equation, so the prescribed compatible pair can be attained. Such an edge layer therefore retains \(d-2\) free parameters; a non-edge layer retains \(d-1\).

Summed over layers \(j=1,\ldots,60\), this gives 462 syzygy parameters before pinning the edge coefficients and 450 after the 12 edge layers are pinned. The exact positive-control script [f2_75_125_early_layer_lift.py](../src/f2_75_125_early_layer_lift.py) constructs a rational truncation through output \(y^{59}\), using \(p_0=x^3+x\), \(q_0=x^5+x\), and \(h(t)=1+t^4\) on the top edge. This is a partial-truncation witness, not a full Keller pair. It confirms that the first 60 y-adic slices alone cannot provide the exclusion; the next obstruction must use the \(Q\)-only layers \(60\le n\le99\), the closed tail \(100\le n\le159\), or additional global constraints.

## The upper y-edge creates a second carrier

The coefficient of \(y^{159}\) in \([P,Q]\) is the polynomial

\[
100p_{60}'q_{100}-60p_{60}q_{100}'=0.
\]

The forced corner monomials imply \(\deg p_{60}=15\) and \(\deg q_{100}=25\). Dividing the displayed relation by \(20\) and taking a logarithmic derivative gives

\[
\left(q_{100}^3/p_{60}^5\right)'=0.
\]

Over \(\overline K[x]\), unique factorization therefore gives

\[
p_{60}=aR(x)^3,\qquad q_{100}=bR(x)^5,\qquad \deg R=5,
\]

for \(a,b\ne0\). This degree-five \(y\)-edge carrier is distinct from the quartic \(h(t)\) in the \((5,-1)\) Newton-edge chart. Their common corner coefficients are linked by the top-edge equation, but no equality between \(R\) and \(h\) is assumed. The subrectangle bounds do not force these top-\(y\) slices to be monomials: lower-\(x\) terms at \(y=60\) or \(y=100\) lie on a boundary face and are not ruled out by uniqueness of maximizers for strictly positive weights. Thus the full degree-five \(R\) branch must be retained.

### Exact valuation pruning for the pure fifth-power upper-y carrier

Here \(R(x)\) is the univariate carrier of the top y-slices. It is distinct
from the ordinary total-degree carrier \(H\) and from the quartic Newton-edge
polynomial \(h(t)\).

There is a branch-specific consequence of the first four descending upper-y
equations. Suppose the degree-five carrier in the top y-slices is a pure
fifth power, \(R(x)=\lambda(x-\alpha)^5\). The source translation
\(x\mapsto x+\alpha\) preserves the F2 outer support envelope, and nonzero
rescaling of \(P,Q\) lets us normalize the top slices to
\(p_{60}=x^{15}\), \(q_{100}=x^{25}\). Write
\(A_s=p_{60-s}\), \(B_s=q_{100-s}\). The fixed part of the depth-\(s\)
equation is

\[
C_s=\sum_{i=1}^{s-1}\bigl((100-s+i)A_i'B_{s-i}-(60-i)A_iB_{s-i}'\bigr),
\]

while its new-slice terms have x-degree at least 14 from \(B_s\) and at
least 24 from \(A_s\). The depth-1 equation gives
\(B_1=\frac53x^{10}A_1\). If \(m=\operatorname{ord}_0(A_1)<\infty\),
the leading term of \(C_2\) has degree \(2m+9\) and coefficient
\(\frac{50}{3}(4m-59)c^2\), where \(c=[x^m]A_1\ne0\). Thus depths 2 and 3
rule out \(m=0,1,2,3,4,5,6\): for \(m=0,1,2\), \(C_2\) has a nonzero
term below degree 14; for \(3\le m\le6\), solving the forced \(B_2\)
coefficient gives a nonzero \(C_3\) term of degree \(3m-6<14\).

The borderline \(m=7\) survives depth 3 but not depth 4. The successive
forced coefficients are
\[
 [x^9]B_2=\frac59c^2,\qquad
 [x]B_3=-\frac5{81}c^3,
\]
so the fixed depth-4 residual has coefficient
\([x^7]C_4=-\frac{3100}{81}c^4\ne0\). Neither new depth-4 slice reaches
x-degree 7, and the other fixed summands start at degrees at least 8 and
16. Therefore this branch must satisfy
\[
\operatorname{ord}_0(p_{59})\ge8\quad\text{or}\quad p_{59}=0.
\]
The exact coefficient arithmetic is replayed by
[f2_75_125_pure_y_carrier_valuation.py](../src/f2_75_125_pure_y_carrier_valuation.py).
This forces the eight coefficients of \(x^0,\ldots,x^7\) in \(p_{59}\) to
vanish in one upper-y carrier stratum. The next two paragraphs exclude the
order-8 and order-9 cases; orders at least 10 and \(p_{59}=0\) still remain.

The \(m=8\) stratum is excluded one row later. Put
\(c=[x^8]p_{59}\ne0\) and \(u=[x]p_{58}\). At depth 4, the \(x^9\)
residual row forces \([x^0]p_{58}=0\), while the \(x^{11}\) row gives
\[
27u^2-9c^2u+c^4=0.
\]
At depth 5, the coefficient of \(x^4\) in the fixed residual is
\[
 [x^4]C_5=\frac{25}{3}c\,u(c^2-18u).
\]
The new depth-5 slices start at \(x\)-degrees 14 and 24, so they cannot
alter this row. The two necessary equations have the explicit certificate
\[
(90u+7c^2)(27u^2-9c^2u+c^4)
+(135u-27c^2)u(c^2-18u)=7c^6.
\]
Since \(c\ne0\), the right side is nonzero, contradicting both equations.
Thus the \(m=8\) valuation branch is excluded. The checker verifies the
depth-4 selector, the depth-5 residual with generic higher slice
coefficients, and this polynomial identity. The selector is the same
quadratic as the normalized F2 root-jet equation, but the polynomial alone
does not identify the parameters or establish that carrier crosswalk.

The \(m=9\) stratum is excluded by one more exact local certificate. Put
\(c=[x^9]p_{59}\ne0\) and \(u=[x^3]p_{58}\). The depth-4 rows of degrees
9, 11, and 13 force the coefficients of \(x^0,x^1,x^2\) in \(p_{58}\) to
vanish. Its degree-15 row then gives the same selector
\(F=27u^2-9c^2u+c^4=0\). At depth 5, the degree-9 fixed residual is
\[
[x^9]C_5=\frac{575}{81}c\,u(c^2-18u)+\frac{4025}{243}cF.
\]
The new depth-5 slices start at degrees 14 and 24, so this coefficient must
vanish. With \(G=u(c^2-18u)\), the identity
\[
(90u+7c^2)F+(135u-27c^2)G=7c^6
\]
again contradicts \(c\ne0\). The exact checker keeps the unforced higher
slice coefficients symbolic. Consequently the pure fifth-power upper-y
carrier branch now requires \(p_{59}=0\) or
\(\operatorname{ord}_x(p_{59})\ge10\); those cases remain open.

There is an exact operator description of every subsequent upper-\(y\) slice. At depth \(s\in\{1,\ldots,60\}\), let

\[
u=p_{60-s},\qquad v=q_{100-s},\qquad
D_P(s)=3+\left\lfloor\frac{60-s}{5}\right\rfloor,\quad
D_Q(s)=5+\left\lfloor\frac{100-s}{5}\right\rfloor.
\]

Then \(D_Q(s)=D_P(s)+10\). The part of the \(y^{159-s}\) coefficient involving the two new layers is

\[
\mathcal L_s(u,v)=100u'q_{100}-(60-s)u q_{100}'
+(100-s)p_{60}'v-60p_{60}v'.
\]

After substituting \(p_{60}=aR^3,\ q_{100}=bR^5\), it factors exactly as

\[
\mathcal L_s(u,v)
=R^2\mathcal D_s\!\left(3av-5bR^2u\right),\qquad
\mathcal D_s(w)=(100-s)R'w-20Rw'.
\]

The map \((u,v)\mapsto 3av-5bR^2u\) is onto the polynomials of degree at most \(D_Q(s)\). Hence the image of the new-layer operator is exactly

\[
R^2\,\mathcal D_s\bigl(\overline K[x]_{\le D_Q(s)}\bigr).
\]

This gives a compact necessary condition for each slice: after moving the contributions from already-fixed shallower layers to the other side, the resulting polynomial must be divisible by \(R^2\), and the quotient must lie in the image of \(\mathcal D_s\). The rank is

\[
\operatorname{rank}\mathcal L_s=D_Q(s)+1-\epsilon_s,
\]

where \(\epsilon_s=1\) exactly when a nonzero polynomial \(w\) of degree at most \(D_Q(s)\) satisfies \(w'/w=((100-s)/20)R'/R\); otherwise \(\epsilon_s=0\). In root-multiplicity terms, for a degree-five \(R\) with at least two distinct roots this occurs only at \(s=20,40,60\). For \(R=c(x-\alpha)^5\), it occurs at every \(s\) divisible by \(4\) between \(4\) and \(60\). Thus the full coefficient-space cokernel has dimension \(14+\epsilon_s\).

### The first Q-only equation in the descending top-y block

At descending depth \(s=61\), the output is the coefficient of \(y^{98}\). The P slice \(p_{60-s}\) is absent, while \(q_{100-s}=q_{39}\) is still allowed and satisfies \(\deg q_{39}\le12\). The terms already fixed at depths \(1,\ldots,60\) are

\[
C_{61}=\sum_{i=1}^{60}\left((39+i)p_{60-i}'q_{39+i}-(60-i)p_{60-i}q_{39+i}'\right).
\]

The exact \(y^{98}\) equation is

\[
C_{61}+R^2\mathcal D_{61}(3a q_{39})=0,
\qquad \mathcal D_{61}(w)=39R'w-20Rw'.
\]

For a monomial \(q=x^d\), \(0\le d\le12\), the new-layer map has leading term

\[
R^2\mathcal D_{61}(3a x^d)=3a r_5^3(195-20d)x^{d+14}+\text{lower powers of }x,
\]

where \(r_5\ne0\) is the leading coefficient of \(R\). Since \(195-20d\ne0\) for every integer \(d\), these thirteen images have distinct degrees \(14,\ldots,26\). Thus the operator on \(q_{39}\) is injective of rank \(13\), inside the \(27\)-dimensional target of polynomials of degree at most \(26\). Solvability therefore imposes exactly fourteen independent linear compatibility conditions on \(C_{61}\); this is the first top-down slice where the P-side no longer supplies a new variable. They can be extracted without an unproved pivot: eliminate degrees \(26,25,\ldots,14\) in that order using the known nonzero leading coefficients above; the resulting remainder has degree at most \(13\), and its fourteen coefficients must vanish. It is a concrete residual target, not yet a contradiction. The support bound, slice degree, uniform leading-coefficient rank, and exact image-remainder controls are replayed by [f2_75_125_q_only_first_slice.py](../src/f2_75_125_q_only_first_slice.py).

The first upper-\(y\) slice is especially rigid. At depth \(s=1\), the equation is \(R^2\mathcal D_1(w_1)=0\), with \(w_1=3a q_{99}-5bR^2p_{59}\). The differential equation \(\mathcal D_1(w_1)=0\) would require \(w_1'/w_1=(99/20)R'/R\). Since every root multiplicity of a degree-five \(R\) is at most five, no nonzero polynomial \(w_1\) can satisfy this. Hence

\[
q_{99}=\frac{5b}{3a}R^2p_{59}.
\]

At depth \(s=2\), the middle-layer contribution is

\[
C_2=99p_{59}'q_{99}-59p_{59}q_{99}'
=\frac{5b}{3a}\left(40R^2p_{59}p_{59}'-118RR'p_{59}^2\right).
\]

The new-layer image is divisible by \(R^2\), so the depth-two equation implies \(R^2\mid C_2\) and

\[
\mathcal D_2(w_2)=-C_2/R^2,\qquad w_2=3a q_{98}-5bR^2p_{58}.
\]

Let \(\alpha\) be a root of \(R\) of multiplicity \(m\le5\), and let \(r=\operatorname{ord}_\alpha(p_{59})\). The leading coefficient of \(C_2\) at \(\alpha\) is proportional to \(40r-118m\), which never vanishes for these integer multiplicities. If \(w_2\) has order \(\ell\), the leading coefficient of \(\mathcal D_2(w_2)\) is proportional to \(98m-20\ell\), also nonzero for \(m\le5\). Matching orders gives \(\ell=2r-m\), so every root must satisfy

\[
\operatorname{ord}_\alpha(p_{59})\ge\left\lceil\frac m2\right\rceil.
\]

Equivalently, the polynomial \(\prod_\alpha(x-\alpha)^{\lceil m_\alpha/2\rceil}\) divides \(p_{59}\), with \(m_\alpha=\operatorname{ord}_\alpha(R)\); \(p_{59}=0\) is included. The finite leading-coefficient census is replayed by [f2_75_125_y_top_edge_jet.py](../src/f2_75_125_y_top_edge_jet.py).

The depth-three equation strengthens this jet bound. Put \(p_1=p_{59}\), \(p_2=p_{58}\), \(q_1=q_{99}\), and \(q_2=q_{98}\), and write \(w_2=3a q_2-5bR^2p_2\). After using \(q_1=(5b/3a)R^2p_1\), its already-fixed-layer contribution is exactly

\[
C_3=\frac{5b}{3a}\left(40R^2(p_1p_2)'-234RR'p_1p_2\right)
+\frac{98p_1'w_2-59p_1w_2'}{3a}.
\]

The new-layer equation is \(R^2\mathcal D_3(w_3)=-C_3\), where \(w_3=3a q_{97}-5bR^2p_{57}\). At a root of multiplicity \(m\), every polynomial input \(w_3\) gives \(\operatorname{ord}_\alpha(\mathcal D_3(w_3))\ge m-1\), including resonant cases where its leading coefficient vanishes. Therefore \(\operatorname{ord}_\alpha(C_3)\ge3m-1\).

If \(p_1\) has order \(r\) and \(p_2\ne0\) has order \(t\), the two possible leading orders in \(C_3\) are \(A=3r-m-1\) and \(B=2m+r+t-1\). Their leading coefficients are nonzero: they are proportional to \(59m-20r\) and \(40(r+t)-234m\), respectively. If \(A<B\), the required order \(3m-1\) forces \(r\ge\lceil4m/3\rceil\). If \(B<A\), the inequality itself gives \(r>3m/2\), which implies the same bound. If \(A=B\), then \(r=(3m+t)/2\ge3m/2\), again implying it. When \(p_2=0\), only the \(A\)-term remains and yields the same bound. Thus any nonzero \(p_{59}\) must satisfy

\[
\operatorname{ord}_\alpha(p_{59})\ge\left\lceil\frac{4m_\alpha}{3}\right\rceil\quad\text{at every root of }R.
\]

Equivalently, \(\prod_\alpha(x-\alpha)^{\lceil4m_\alpha/3\rceil}\mid p_{59}\). Since \(\deg R=5\), this forced divisor has degree between 7 and 10, depending on the multiplicity partition of 5; it does not by itself contradict \(\deg p_{59}\le14\). The script checks the exact identities and all bounded integer multiplicity/order cases.

The fourth slice adds a further leading-jet test. Set \(p_3=p_{57}\) and \(w_3=3a q_{97}-5bR^2p_{57}\). The depth-four already-fixed contribution is

\[
\begin{aligned}
C_4={}&\frac{5b}{3a}\left(40R^2(p_1p_3)'-232RR'p_1p_3
+40R^2p_2p_2'-116RR'p_2^2\right)\\
&+\frac{97p_1'w_3-59p_1w_3'+98p_2'w_2-58p_2w_2'}{3a},
\end{aligned}
\]

and \(R^2\mathcal D_4(w_4)=-C_4\), with \(w_4=3a q_{96}-5bR^2p_{56}\), so again \(\operatorname{ord}_\alpha(C_4)\ge3m-1\). At the minimum \(r=\lceil4m/3\rceil\) from the depth-three bound, when its \(A\)-term is uniquely leading, \(w_3\) has order \(3r-4m\). The \(p_1w_3\) term in \(C_4\) is then a unique nonzero term below order \(3m-1\), except in two tied depth-three cases. The exact order census gives:

- For \(m=3,r=4\), the unique low term rules out equality; hence \(r\ge5\).
- For \(m=5,r=7\), it rules out equality; hence \(r\ge8\).
- For \(m=2,r=3\) and \(m=4,r=6\), if \(\operatorname{ord}_\alpha(p_{58})\ge1\) (or \(p_{58}=0\)), the depth-three \(A\)-term is uniquely leading and the same fourth-slice contradiction applies. The only remaining possibility at this stage is \(\operatorname{ord}_\alpha(p_{58})=0\), where the depth-three leading orders tie. Cancelling them forces the leading coefficient \(e\) of \(p_{58}\), after normalizing the local parameter so the initial term of \(R\) is \(u^m\) and writing the leading coefficient of \(p_{59}\) as \(d\), to satisfy \(e=d^2/(18a)\). For these two multiplicities, \(97m-20\ell\ne0\) for every integer \(\ell\), so the depth-three operator preserves the order of its leading jet. It follows that the \(p_1w_3\) term at depth four has order at least \(2m\). After dividing through by the nonzero derivative unit of the local parameter, the coefficient of order \(2m-1\) is exactly \(-116m\kappa e(d^2/(3a)+e)\ne0\), where \(\kappa=5b/(3a)\); all other terms have higher order. This also contradicts the required order \(3m-1\).

Consequently a nonzero \(p_{59}\) must satisfy the stronger root conditions

\[
\operatorname{ord}_\alpha(p_{59})\ge
\begin{cases}
2,&m_\alpha=1,\\
\lceil4m_\alpha/3\rceil+1,&2\le m_\alpha\le5.
\end{cases}
\]

This raises the minimum forced divisor degree for \(p_{59}\) to 8 across the degree-five multiplicity partitions, but still leaves room under \(\deg p_{59}\le14\). The replay script checks the exact depth-three and depth-four operator identities and the bounded local order cases, including the two tied branches. These are local pruning conditions, not an exclusion of the pair.

### The separate \(p_{59}=0\) branch

The root-jet bound does not cover \(p_{59}=0\), so keep that branch explicit. At depths 1, 2, and 3, all prior-layer contributions vanish successively. The homogeneous equations force \(w_1=w_2=w_3=0\): a nonzero kernel element of \(\mathcal D_s\) would have root order \((100-s)m/20\), which is nonintegral for every \(m\in\{1,\ldots,5\}\) and \(s=1,2,3\). Thus \(q_{99}=0\), \(q_{98}=\kappa R^2p_{58}\), and \(q_{97}=\kappa R^2p_{57}\), with \(\kappa=5b/(3a)\).

At depth 4, writing \(p_2=p_{58}\), its fixed-layer term is

\[
C_4=\kappa\left(40R^2p_2p_2'-116RR'p_2^2\right).
\]

At a root of multiplicity \(m\), if \(p_2\ne0\) has order \(t\), the leading coefficient is proportional to \(40t-116m\), never zero for \(m\le5\). Since the new-layer image has order at least \(3m-1\), this forces \(t\ge\lceil m/2\rceil\). If \(p_2=0\), then \(w_4\in\ker\mathcal D_4\). Factor multiplicities show that this kernel is zero unless \(R\) is a pure fifth power; for \(R=c(x-\alpha)^5\), it is one-dimensional, generated by \((x-\alpha)^{24}\). That exceptional branch must be retained.

At depth 5, if both \(p_{58}\) and \(p_{57}\) are nonzero with root orders \(t,z\), then

\[
C_5=\kappa\left(40R^2(p_{58}p_{57})'-230RR'p_{58}p_{57}\right).
\]

For \(t+z<m\), its leading coefficient cannot vanish and its order is too low to lie in the new-layer image. Hence \(t+z\ge m\) at every root. Depth 6 sharpens this when \(p_{58}\ne0\). Put \(p_4=p_{56}\) and \(w_4=3a q_{96}-5bR^2p_{56}\); the fixed-layer term is

\[
\begin{aligned}
C_6={}&\kappa\left(40R^2(p_{58}p_{56})'-228RR'p_{58}p_{56}
+40R^2p_{57}p_{57}'-114RR'p_{57}^2\right)\\
&+\frac{96p_{58}'w_4-58p_{58}w_4'}{3a}.
\end{aligned}
\]

If \(p_{58}\) has order \(t\), depth 4 forces \(w_4\) to have order \(2t-m\); the relevant leading coefficient is \(116m-40t\ne0\) for \(m\le5\), \(t\le14\). The only indicial resonance of \(\mathcal D_4\) at a root with \(m\le5\) is \((m,\ell)=(5,24)\), above \(2t-m\le23\), so it cannot change this leading order. The last term of \(C_6\) then has order \(3t-m-1\) with nonzero coefficient. If \(p_{57}\ne0\), write \(z=\operatorname{ord}_\alpha(p_{57})\). The depth-5 condition gives \(t+z\ge m\). At depth 6, the three possible leading orders are \(A=3t-m-1\) from the \(w_4\)-term, \(B=2m+t+v-1\) from \(p_{58}p_{56}\) (where \(v=\operatorname{ord}_\alpha(p_{56})\)), and \(C=2m+2z-1\) from \(p_{57}^2\). Their coefficients are nonzero for \(1\le m\le5\) and the support caps. The exact bounded order census shows that a term below the required image order \(3m-1\) can survive only in the tied stratum \(t=m,z=0\); there the depth-4 equation gives \(w_4=F u^m+\cdots\), with \(F=\kappa e^2/c\) if \(R=c u^m+\cdots\) and \(p_{58}=e u^m+\cdots\), while \(p_{57}=d+\cdots\). Cancellation of the depth-6 coefficient of order \(2m-1\) forces \(e^3=9a c^3d^2\). The next slice removes this exception. Set \(p_5=p_{55}\), \(w_5=3a q_{95}-5bR^2p_{55}\); the exact depth-seven fixed term is

\[
\begin{aligned}
C_7={}&\kappa\left(40R^2\bigl((p_{58}p_{55})'+(p_{57}p_{56})'\bigr)
-226RR'(p_{58}p_{55}+p_{57}p_{56})\right)\\
&+\frac{95p_{58}'w_5-58p_{58}w_5'+96p_{57}'w_4-57p_{57}w_4'}{3a}.
\end{aligned}
\]

At the tied stratum, depth 5 forces \(w_5=G+\cdots\) with \(G=2\kappa ed/c\). The only term of order \(m-1\) in \(C_7\) is then

\[
\frac{m}{3a}(95eG-57dF)u^{m-1}
=\frac{133m\kappa e^2d}{3ac}u^{m-1},
\]

which is nonzero. All other fixed terms have order greater than \(m-1\), whereas every new-layer term \(R^2\mathcal D_7(w_7)\) has order at least \(3m-1\). This contradiction excludes the tied stratum. Thus, when both \(p_{58}\) and \(p_{57}\) are nonzero, the local conditions sharpen to \(t\ge\lceil4m/3\rceil\) and \(z\ge\lceil m/2\rceil\). If \(p_{57}=0\), depth 6 itself forces \(t\ge\lceil4m/3\rceil\). When \(p_{58}=0\) and \(p_{57}\ne0\), the depth-6 diagonal term gives \(\operatorname{ord}_\alpha(p_{57})\ge\lceil m/2\rceil\), since \(40z-114m\ne0\) for \(m\le5\). The homogeneous \(w_4\) resonance does not enter that equation because it is multiplied by \(p_{58}=0\). The exact identities and local order census are replayed by [f2_75_125_zero_first_slice.py](../src/f2_75_125_zero_first_slice.py).

This opens the zero-slice branch into separate cases—\(p_{58}\ne0\), \(p_{58}=0,p_{57}\ne0\), and a deeper zero tail with a possible pure-power resonance. The new divisibility bounds still fit the available slice degrees, so none is yet excluded.

There is also a first coupling condition in the \(p_{58}=0,p_{57}\ne0\) branch. For a degree-five carrier, \(\mathcal D_5\) has no polynomial kernel: a kernel element would need root orders \(19m_\alpha/4\) at every root, which would require every multiplicity in the partition of 5 to be divisible by 4. Thus \(w_5=0\) after the zero depth-five fixed term. If \(R\) has at least two distinct roots, \(\mathcal D_4\) also has no kernel, so \(w_4=0\). At depth 7, if \(p_{56}\ne0\), the fixed term reduces to

\[
C_7=\kappa\left(40R^2(p_{57}p_{56})'-226RR'p_{57}p_{56}\right).
\]

At a root of multiplicity \(m\le4\), putting \(z=\operatorname{ord}(p_{57})\) and \(v=\operatorname{ord}(p_{56})\), the leading coefficient for \(z+v<m\) is \(40(z+v)-226m\ne0\), and the term is below the required image order. Hence \(z+v\ge m\). If \(p_{56}=0\) but \(p_{55}\ne0\), the depth-eight cross term is

\[
C_8=\kappa\left(40R^2(p_{57}p_{55})'-224RR'p_{57}p_{55}\right),
\]

and the same argument, using \(40s-224m\ne0\) for \(m\le4\), gives \(\operatorname{ord}(p_{57})+\operatorname{ord}(p_{55})\ge m\). These are coupled divisibility conditions, not a contradiction: later slices or global degree constraints are still needed. The pure fifth-power carrier remains separate because its possible \(w_4\) kernel can enter depth 7. The symbolic depth-seven/eight identities and the finite resonance checks are included in [f2_75_125_zero_first_slice.py](../src/f2_75_125_zero_first_slice.py).

There is a useful extension to the first-nonzero-slice split. Let \(P_s=p_{60-s}\), \(Q_s=q_{100-s}\), and \(w_s=3aQ_s-5bR^2P_s\). Suppose \(P_k\ne0\) is the first nonzero \(P\)-slice and the homogeneous operators \(\mathcal D_j\) have no kernel for \(j\le k\). Then \(w_j=0\) for \(j<k\), and \(w_k=0\). At depth \(2k\), every fixed-layer term vanishes except the diagonal pair \((P_k,Q_k)\), so

\[
C_{2k}=\kappa\left(20R^2(P_k^2)'-(120-2k)RR'P_k^2\right),\qquad \kappa=\frac{5b}{3a}.
\]

For a carrier with at least two distinct roots, the first homogeneous-kernel depth is 20. Thus when \(2\le k\le19\) and \(P_k\) is the first nonzero slice, the order condition at \(2k\) forces \(\operatorname{ord}_\alpha(P_k)\ge\lceil m_\alpha/2\rceil\), except when the leading coefficient vanishes. Within the actual support caps, the exceptional \((k,m,t)\)—depth, root multiplicity, and slice order—are exactly

\[
(4,5,14),\ (5,4,11),\ (8,5,13),\ (10,2,5),\ (10,4,10),\ (12,5,12),\ (15,4,9),\ (16,5,11).
\]

These are resonance candidates for the *leading coefficient test*, not established solutions; higher jets can still eliminate them. For a pure fifth-power carrier, the first homogeneous-kernel depth is 4, but the same half-multiplicity condition has no leading-coefficient exception for first nonzero depths \(k=2,3\). At \(k\ge4\), the resonant \(w_4\) branch must be carried explicitly. The general diagonal identity and the finite exceptional census are included in [f2_75_125_zero_first_slice.py](../src/f2_75_125_zero_first_slice.py).

The exact rank replay [f2_75_125_y_top_edge_rank_audit.py](../src/f2_75_125_y_top_edge_rank_audit.py) checks all 60 depths in the squarefree example \(R=x^5+1\) and the pure-power example \(R=x^5\). The factorization and rank formula above are uniform in \(R\); the two examples serve as exact arithmetic controls. The recurrence still has to satisfy all subsequent compatibility conditions, so the calculations here are local reductions through depth four for \(p_{59}\ne0\), the depth-seven tie exclusion and selected depth-eight cross condition on \(p_{59}=0\), not an exclusion of the degree pair.

### Resonances across the entire Q-only tail

For descending depths \(61\le s\le100\), only the Q-slice \(q_{100-s}\) is new. Put \(u=100-s\) and \(D=5+\lfloor u/5\rfloor\), so \(\deg q_{100-s}\le D\). The new-layer image is \(R^2\mathcal D_s(3a q_{100-s})\), and its target degree is at most \(D+14\). A nonzero kernel element must satisfy

\[
q' / q = \frac{u}{20}\,R'/R,
\qquad\text{hence}\qquad q=cR^{u/20}.
\]

If the root multiplicities of \(R\) are \(m_1,\ldots,m_r\), this is polynomial exactly when \(20\mid u m_i\) for every \(i\); its degree \(u/4\) then lies within the support cap. The complete partition census is therefore:

| Root multiplicities of \(R\) | Resonant depths \(61\le s\le100\) |
| --- | --- |
| \(5\) | \(64,68,72,76,80,84,88,92,96,100\) |
| Any other partition of \(5\) | \(80,100\) |

Thus the one-slice operator has rank \(D+1\) and cokernel dimension \(14\) off resonance; on each listed resonant depth it has rank \(D\) and cokernel dimension \(15\). The \(s=100\) kernel is the constant term of \(Q\), which the Jacobian cannot detect. This census makes every exceptional pivot explicit throughout the forty-slice Q-only block, but the residual compatibility equations still need to be imposed. The 40 support-target caps and all seven multiplicity partitions are replayed by [f2_75_125_q_only_resonance_audit.py](../src/f2_75_125_q_only_resonance_audit.py).

## A zero-Jacobian control for the two F2 carriers

The quartic edge carrier and the degree-five upper-\(y\) carrier can coexist with essentially independent coefficients, subject only to their shared corner. Let

\[
h(t)=(t-1)^2(s_0+s_1t+s_2t^2),\qquad
R(x)=s_2x^5+r_4x^4+r_3x^3+r_2x^2+r_1x+r_0,
\]

with \(s_0s_2\ne0\), and define

\[
A(x,y)=x\,h(xy^5)+y^{20}\bigl(R(x)-s_2x^5\bigr).
\]

Then \(\operatorname{Supp}(A)\subseteq\{0\le i\le5,\ 0\le j\le20,\ 5i-j\le5\}\). Its \((5,-1)\)-edge is exactly \(x h(xy^5)\), while its \(y^{20}\)-slice is exactly \(R(x)\). Therefore \(P=A^3,\ Q=A^5\) lie in the exhaustive \((75,125)\) outer support envelopes, have edge faces \((xh)^3,(xh)^5\), and upper-\(y\) faces \(R^3,R^5\). The exact chain-rule identity gives \([P,Q]=0\).

Thus the zero-right-hand-side band equations admit a family with independently variable \(h\) and \(R\) (apart from the common leading coefficient). This rules out an exclusion based only on those homogeneous face rows. It does not satisfy the Keller equation \([P,Q]=1\); the inhomogeneous low-degree rows must supply the missing obstruction. The support sums, face identities, and Jacobian identity are replayed by [f2_75_125_dual_carrier_zero_jacobian_control.py](../src/f2_75_125_dual_carrier_zero_jacobian_control.py).

## Relation to the proposed \(H=x^{25}+y^{25}\)

There are two different objects that can be called a common carrier:

- The ordinary total-degree leading forms of a pair of degrees \((75,125)\) are necessarily proportional to \(H^3,H^5\) for some homogeneous \(H\) of degree \(25\).
- In this \(F_2\) edge chart, the edge carrier is the nonhomogeneous polynomial \(C=xh(xy^5)\), with \(h\) quartic. Its highest total-degree term is proportional to \(x^5y^{20}\).

Thus \(H=x^{25}+y^{25}\) is not the leading carrier in the standardized \(F_2\) chart; there the homogeneous carrier is the two-factor monomial \(x^5y^{20}\). The factor-pattern strategy could instead be applied to the quartic edge polynomial \(h(t)\). The front-end reduction recorded below supplies the link between the generic top-form picture and this chart, subject to its cited inputs and proof extraction.

There is a stronger conditional consequence of the front-end reduction recorded above. In the standardized F2 representative, the GGV rectangle for \(A_0=(5,20)\) has a unique vertex maximizing ordinary total degree, namely \((15,60)\) for \(P\) and \((25,100)\) for \(Q\). Hence that representative's top homogeneous pieces are single monomials proportional to \((x^5y^{20})^3\) and \((x^5y^{20})^5\). The degree-transport identity gives \(g_0=M\bar a+N\bar b=25\); since the census forces \((\bar a,\bar b)=(5,20)\) and \(g=25\), positivity forces \(M=N=1\). Thus the final normalization is affine and preserves the multiplicities of the factors in this representative. Subject to the cited GGV inputs and written reduction, every hypothetical \((75,125)\) pair can therefore be moved to a representative with carrier multiplicity type \(5+20\). This is a conditional normal-form reduction, not a classification of the ordinary top form in the original coordinates: the preceding polynomial automorphisms need not preserve ordinary degree forms or their factor counts.

### What the fixed Fermat-type form does and does not cover

For any hypothetical pair with ordinary top homogeneous pieces \(P_{75}\) and \(Q_{125}\), the degree-198 part of \([P,Q]=1\) is \([P_{75},Q_{125}]=0\). In characteristic zero this implies

\[
P_{75}=aH^3,\qquad Q_{125}=bH^5
\]

for a homogeneous degree-25 form \(H\) and nonzero scalars \(a,b\). One way to see this is to write the forms in the chart \(x\ne0\); the bracket equation becomes \(3FG'-5F'G=0\), so \(G^3/F^5\) is constant, and unique factorization gives the powers.

Over \(\mathbb C\),

\[
x^{25}+y^{25}=\prod_{\zeta^{25}=-1}(x-\zeta y)
\]

has 25 distinct linear factors. It lies in the at-least-three-factor stratum,
which Abhyankar's at-most-two-points-at-infinity theorem excludes for a plane
Jacobian pair. It is not a representative for the whole stratum: general
degree-25 binary forms have root multiplicity partitions and cross-ratio data
not equivalent to this special form under a linear change of variables.

Before imposing the stronger \(F_2\) reduction, the generic top-form condition has the exhaustive factor split

\[
H=c\prod_{i=1}^{r} L_i^{m_i},\qquad \sum_i m_i=25,
\]

with distinct linear forms \(L_i\). The cases \(r=1\), \(r=2\), and \(r\ge3\) are exhaustive. Abhyankar's two-point theorem excludes \(r\ge3\). The one-factor stratum is not excluded by the universal one-point criterion: that criterion quantifies over every Jacobian pair and is equivalent to the Jacobian conjecture; it does not establish automorphy of a single pair. The conditional \(F_2\) reduction standardizes the surviving carrier to multiplicity type \(5+20\), subject to its stated source hypotheses. The F2 quartic edge patterns \(2+2\) and \(2+1+1\) are residual branches within that normalized route, and still require a contradiction from the lower coefficients. A necessary local root condition alone is not an exclusion unless it contradicts every allowed point of its branch or is combined with the remaining equations.

The finite degree-transport, corner, and factor-count arithmetic is replayed by [f2_75_125_top_form_branch_audit.py](../src/f2_75_125_top_form_branch_audit.py). It confirms the two-factor \(5+20\) consequence and that \(x^{25}+y^{25}\) has a different factor count; it is conditional on the standardization/F2 theorem inputs and does not exclude the pair.

Finally, the exact pair \((H^3,H^5)\) has Jacobian zero because its two entries are algebraically dependent. The candidate is a pair \((P,Q)\) whose *top homogeneous parts* are proportional to these powers and whose lower terms would still have to make \([P,Q]=1\).

### A total-degree face audit

The same support envelope can be sliced by ordinary total degree. For deficit
\(0\le d\le15\), its entire homogeneous slices are

\[
P_{75-d}=\sum_{k=0}^{\lfloor5d/6\rfloor}p_{d,k}x^{15-d+k}y^{60-k},\qquad
Q_{125-d}=\sum_{k=0}^{\lfloor5d/6\rfloor}q_{d,k}x^{25-d+k}y^{100-k}.
\]

After scaling the nonzero top coefficients to 1 (which preserves the vanishing
of every positive-degree Jacobian row), the equation of total output degree
\(198-d\) has, at the monomial indexed by \(k\), the form

\[
C_{d,k}+5(4d-5k)(3q_{d,k}-5p_{d,k})=0,
\]

where \(C_{d,k}\) is determined entirely by the previously chosen slices.
Away from \(5k=4d\), this row solves for \(q_{d,k}\) after choosing
\(p_{d,k}\). At the resonances \((d,k)=(5,4),(10,8),(15,12)\), the
new coefficients disappear from the row; exact symbolic substitution gives
\(C_{d,k}=0\) identically in each case. There are also no nonzero output
coefficients beyond the support cap \(k\le\lfloor5d/6\rfloor\).

The checker [f2_75_125_total_degree_frontier.py](../src/f2_75_125_total_degree_frontier.py) enumerates every support slice and replays these equations exactly through deficit 15. It checks 108 supported rows and 27 out-of-support rows; the resulting high-degree truncation still has 111 free coefficients. Thus the ordinary total-degree face supplies no early rigidity block through output degree 183. This is a flexibility result only: the lower-degree rows and their coupling to the Newton-edge and upper-\(y\) recurrences remain to be imposed.

The exact recurrence has now been extended to deficit 20. It verifies 186
supported coefficient rows and the out-of-support rows, including the new
resonance \((d,k)=(20,16)\); all four resonance compatibility polynomials
through \(d=20\) vanish identically. The truncated family therefore still
has 175 free coefficients through output degree 178. This adds no rigidity,
but it clears the first resonance beyond the earlier 111-parameter block.
The extension is replayed by
[f2_75_125_total_degree_resonances.py](../src/f2_75_125_total_degree_resonances.py);
an attempted continuation through deficit 30 had no captured terminal output
when its execution session disappeared, so no result beyond deficit 20 is
counted here.

Terminology: the disjoint exhaustive split is by the number \(r\) of distinct linear factors. The “pure-power” case above means \(r=1\), namely \(H=L^{25}\). If “pure power” instead means \(H=G^e\) for any \(e>1\), it overlaps the other cases: the surviving F2 form satisfies \(x^5y^{20}=(xy^4)^5\) while having two distinct linear factors. State the cases by \(r\) and the multiplicity partition to avoid that overlap.

### A power relation gives an exact edge-root selector

The edge identity \(P_{15}=aU^3\), \(Q_{25}=bU^5\), with \(U=th\), also supplies a useful derived polynomial. Set

\[
W=Q^3-\lambda P^5,\qquad \lambda=b^3/a^5.
\]

The weight-75 faces cancel, so \(\omega(W)\le74\), where \(\omega(x^iy^j)=5i-j\). Since \([P,Q]=1\),

\[
[P,W]=3Q^2.
\]

The bracket lowers this weight by at least four: \(\omega([P,W])\le15+\omega(W)-4\). The right side has top weight 50, hence \(\omega(W)\ge39\).

If the top band of \(W\) has weight \(s>39\), write it as \(z^s g_s(t)\). Its bracket with \(P_{15}\) must vanish at weight \(s+11>50\), giving

\[
5U g_s'-sU'g_s=0.
\]

Because \(U=th\) has a simple root at \(t=0\), a nonzero polynomial solution requires \(5\mid s\). For every possible weight above 39, a power product cancels that face:

| Weight \(s\) | Product with that top weight | Weight of its bracket correction \([P,P^iQ^j]\) |
| ---: | --- | ---: |
| 40 | \(PQ\) | 15 |
| 45 | \(P^3\) | 0 (the bracket is zero) |
| 50 | \(Q^2\) | 25 |
| 55 | \(P^2Q\) | 30 |
| 60 | \(P^4\) | 0 (the bracket is zero) |
| 65 | \(PQ^2\) | 40 |
| 70 | \(P^3Q\) | 45 |

Subtracting suitable multiples successively removes all bands above 39. Every nonzero bracket correction has weight below 50, so the weight-50 equation remains unchanged. The resulting corrected polynomial \(W^*\) has top band exactly 39. This correction step is important: \(W^*\) need not equal \(Q^3-\lambda P^5\), but it still has the same bracket equation at weight 50.

Write that band as \(z^{39}g(t)\). Equating weight 50 in \([P,W^*]=3Q^2+\text{terms of weight below }50\) gives

\[
5Ug'-39U'g=\mu U^8,\qquad \mu=b^2/a\ne0.
\]

At a root of \(U\) of multiplicity \(e=1\) or 2, if \(g\) has order \(\ell\), the leading coefficient on the left is proportional to \(5\ell-39e\), which is nonzero. If \(\ell<7e+1\), its order \(e+\ell-1\) is below the right side's order \(8e\), impossible. Therefore

\[
\operatorname{ord}_{\tau}(g)\ge7e+1
\]

at every root \(\tau\) of \(U\). At infinity, with \(\deg U=5\), the leading term has degree \(\deg(g)+4\) and coefficient proportional to \(\deg(g)-39\), unless \(\deg(g)=39\) and that term cancels. Matching the right side's degree 40 leaves \(\deg(g)\in\{36,39\}\). The root-order lower bounds total 38 in the \(2+2\) case and 39 in the \(2+1+1\) case, forcing \(\deg(g)=39\) in either case.

It follows that \(g=U^7V\), where \(\deg V=4\) and the squarefree part of \(U\) divides \(V\). After rescaling \(V\), its equation is

\[
5UV'-4U'V=U.
\]

This gives a finite residual calculation in the two root patterns already forced by the F2 successor data:

* **Pattern \(2+2\).** Set \(h=(t-1)^2(t-\alpha)^2\), with \(\alpha\ne0,1\). Then \(V=t(t-1)(t-\alpha)(mt+n)\), with \(m\ne0\). The coefficients of \(t^7\) and then \(t^6\) force
  \[
  -3(\alpha+1)m-5n=0,\qquad \alpha^2-3\alpha+1=0.
  \]
  The latter quadratic has two distinct nonzero roots, neither equal to 1. For either root, \(m=(4\alpha-11)/3\), \(n=3-\alpha\) satisfies the full polynomial equation.
* **Pattern \(2+1+1\).** Set \(h=(t-1)^2(t-\alpha)(t-\beta)\). Now \(V=c\,t(t-1)(t-\alpha)(t-\beta)\), with \(c\ne0\). The coefficients of \(t^7,t^6\) force \(\alpha+\beta=3\) and \(\alpha\beta=3\), so the two remaining roots satisfy \(t^2-3t+3=0\). They are distinct and avoid 0 and 1; \(c=-1/3\) satisfies the full equation.

Thus the power relation collapses the quartic edge carrier to two explicit algebraic configurations, but does **not** exclude either one. The exact checker [f2_75_125_power_relation_edge_certificate.py](../src/f2_75_125_power_relation_edge_certificate.py) replays the semigroup cancellations, root-order and degree alternatives, and coefficient reductions. It does not prove the support/standardization inputs or produce a completion to a full Keller pair. The result is a necessary edge condition and a smaller target for the subsequent bands, not an exclusion of \((75,125)\).

### The next power-relation band is a small exact linear problem

The next output band, of weight 49 in \([P,W^*]\), couples the edge selector to the first adjacent bands. Normalize the nonzero top-edge scalars and divide by the nonzero Jacobian scalar. Write

\[
P_{15}=U^3,\quad Q_{25}=U^5,\quad P_{14}=A=t^3U_1,
\quad Q_{24}=q_{24}=\frac53t^5h^2U_1,
\quad W^*_{39}=g=U^7V,\quad W^*_{38}=f.
\]

Here \(V\) is the normalized quartic solution in the previous subsection. Since the power-product corrections to \(W\) have bracket weight at most 45, they do not enter this row. The band equation is

\[
15U^3f'-114U^2U'f+14Ag'-39A'g-6U^5q_{24}=0.
\]

The outer rectangle bounds \(W^*\) by \(0\le i\le75,\ 0\le j\le300\). On weight 38 this allows exactly the relaxed polynomial range \(f\in\operatorname{span}\{t^8,\ldots,t^{67}\}\). The fifth through fifteenth band conditions already give \(U_1=h^2S\), \(\deg S\le3\). Thus this row is a homogeneous linear system in 60 coefficients of \(f\) and four coefficients of \(S\), 64 unknowns altogether.

The exact ranks over the relevant coefficient fields are:

| Carrier pattern | Field | Rank | Kernel dimension |
| --- | --- | ---: | ---: |
| \(2+2\), either root \(\alpha=(3\pm\sqrt5)/2\) | \(\mathbb Q(\sqrt5)\) | 61 | 3 |
| \(2+1+1\), \(h=(t-1)^2(t^2-3t+3)\) | \(\mathbb Q\) | 62 | 2 |

When \(U_1=0\), the row reduces to \(5Uf'-38U'f=0\). Its only polynomial solution is \(f=0\), since a nonzero solution would be proportional to \(U^{38/5}\), which has nonintegral order at the simple root \(t=0\). Therefore projection of the displayed kernel to the four coefficients of \(S\) is injective: this row leaves a 3-dimensional allowed \(S\)-space in either \(2+2\) conjugate and a 2-dimensional space in the \(2+1+1\) case. It prunes the first adjacent band but does not force \(U_1=0\) or exclude either carrier pattern. Later rows and the \(U_1=0\) branch must still be handled.

The replay [f2_75_125_power_relation_second_band.py](../src/f2_75_125_power_relation_second_band.py) constructs the full relaxed coefficient matrix and computes these ranks exactly. The reduction uses only this one band plus previously established necessary divisibilities; it is not a Gröbner computation of the full Jacobian system and does not claim an exclusion.

### Weight 48 adds no further condition in the relaxed model

At the next output weight, 48, the bracket pairs are exactly \((15,37),(14,38),(13,39)\). Put \(U_1=h^2S\), \(U_2=hT\), with \(\deg S\le3\), \(\deg T\le7\), and write

\[
A=t^3U_1,\quad A_2=t^3U_2,\quad q_{24}=\frac53t^5h^2U_1,
\quad q_{23}=\frac53t^5h^2U_2+\frac59t^5h^3S^2,
\quad W^*_{39}=g=U^7V,\quad W^*_{38}=f,\quad W^*_{37}=k.
\]

No correction term from the power products contributes here because their brackets have weight at most 45. The row is

\[
15U^3k'-111U^2U'k
+14Af'-38A'f
+13A_2g'-39A_2'g
=6U^5q_{23}+3q_{24}^2.
\]

Allowing the full support-permitted ranges \(k\in\operatorname{span}\{t^8,\ldots,t^{67}\}\) and \(T\in K[t]_{\le7}\), the 72 coefficient equations are linear in their 68 new coefficients. Their exact left-kernel compatibility conditions vanish identically on the smaller weight-49 parameter space in both branches:

| Carrier | Field | Weight-49 parameter dimension | Weight-48 rank | Left-kernel dimension | Nonzero residuals |
| --- | --- | ---: | ---: | ---: | ---: |
| \(2+2\), \(\alpha=(3+\sqrt5)/2\) | \(\mathbb Q(\sqrt5)\) | 3 | 61 | 11 | 0 |
| \(2+1+1\) | \(\mathbb Q\) | 2 | 62 | 10 | 0 |

The other \(2+2\) root is the Galois conjugate of the displayed calculation. Therefore this relaxed weight-48 row allows a solution for every parameter left by weight 49; it does not even shrink those parameter spaces. The linear solution fibres have dimensions 7 and 6, respectively. This shows why further edge rows, or coupling with the original \([P,Q]\) bands, are needed before a unit-ideal certificate can appear. The exact replay is [f2_75_125_power_relation_weight48_audit.py](../src/f2_75_125_power_relation_weight48_audit.py); it checks all rank and zero-residual counts and explicitly labels the system as a relaxation.

### The companion \([Q,W^*]\) band is row-space redundant

The same corrected power relation has the exact identity

\[
[Q,W^*]=5\lambda J P^4+\text{the brackets of the power-product corrections}.
\]

Each correction bracket has weight at most 55, so it cannot enter weight 59. At that weight only \([Q_{25},W^*_{38}]\) and \([Q_{24},W^*_{39}]\) contribute. With the same normalization as above, the necessary row is

\[
5U^4(5Uf'-38U'f)+24q_{24}g'-39q_{24}'g-20U^9A=0.
\]

I stacked this row's exact coefficient matrix with the weight-49 \([P,W^*]\) matrix, using the same 64 coefficients of \(S\) and \(f\). In the \(2+2\) carrier over \(\mathbb Q(\sqrt5)\), the separate ranks are 61 and the stacked rank is 61; in the \(2+1+1\) carrier over \(\mathbb Q\), they are 62 and 62. Since the stacked rank equals each separate rank, the two row spaces coincide in each branch: this companion band adds no linear condition to the relaxed system. The exact replay is [f2_75_125_power_relation_cross_commutator.py](../src/f2_75_125_power_relation_cross_commutator.py). It is a negative result about this proposed rigidity route, conditional on the same edge and support model; it does not exclude either branch.

### One deeper coupled band pair also has no obstruction

I extended the exact calculation by one more level, imposing weight 47 of
\([P,W^*]\) together with weight 57 of \([Q,W^*]\), after solving the
weight-48/58 block. This is a linear system in 72 new coefficients: 60 for
the weight-38 slice and 12 for the next shared power coefficient. In the
\(2+2\) branch over \(\mathbb Q(\sqrt5)\), the matrix rank is 61 and all 91
left-kernel compatibility polynomials vanish identically; in the \(2+1+1\)
branch over \(\mathbb Q\), the rank is 62 and all 90 compatibility
polynomials vanish identically. The resulting affine fibres have dimensions
11 and 10, respectively. Thus this pair of rows also extends every point
remaining from the weight-49 edge kernel; it supplies no polynomial
obstruction in this relaxed model. The exact replay is
[f2_75_125_power_relation_pq_thirdbands.py](../src/f2_75_125_power_relation_pq_thirdbands.py).
This is a verified negative result for the chosen band route, not an
exclusion of either carrier or the degree pair.

The repaired exact \(2+2\) replay has now independently confirmed the
47/57 block: its 152-row, 72-variable matrix has rank \(61\), all 91
left-kernel compatibility residuals vanish on the full preceding affine
family, and the affine solve verifies. Thus its 11-dimensional new-slice
fibre extends every point left by 48/58; this block also adds no restriction.

### The next 46/56 coefficient matrix has a 13-dimensional fibre

The new-slice part of weight 46 in \([P,W^*]\) and weight 56 in
\([Q,W^*]\) is a parameter-independent homogeneous linear system in the 61
coefficients of \(W^*_{35}\) and the 12 coefficients of \(U_4\). An exact
matrix-only extraction gives 72 nonzero rows from the first bracket and 80
from the second. In both the \(2+2\) branch over \(\mathbb Q(\sqrt5)\) and
the \(2+1+1\) branch over \(\mathbb Q\), the combined 152-row matrix has rank
60, hence a 13-dimensional solution fibre. This rank calculation does not
check compatibility with the forcing from earlier bands; forcing-only rows
must also be retained. The independent matrix replay is
[f2_75_125_power_relation_fourthband_matrix.py](../src/f2_75_125_power_relation_fourthband_matrix.py),
and the direct band-convolution formula audit is included in
[f2_75_125_power_relation_fourthbands.py](../src/f2_75_125_power_relation_fourthbands.py).
The full 152-row forcing compatibility has now been verified in both
branches: each has rank \(60\), left-kernel dimension \(92\), and zero
residuals. Thus each leaves the predicted 13-dimensional affine fibre; the
block does not exclude either branch. The \(2+2\) sparse replay initially
stopped because SymPy tried to rediscover the known field
\(\mathbb Q(\sqrt5)\) from very large exact Flint coefficients and overflowed
internally; this was a software failure, not an algebraic obstruction. A
separate dense exact replay has verified the full \(2+2\) compatibility.
The later fifth-band run reached the \(47/57\) pivot solve, then failed when
SymPy's expression domain (`EX`) tried to test a huge exact coefficient for
zero through a floating logarithm. The parameter-polynomial checker now
avoids constructing an `EX` polynomial: it treats \(\sqrt5\) formally,
clears a parameter-free denominator, and reduces coefficients modulo
\(\alpha^2-5\) over \(\mathbb Q\). Smoke checks with \(2^{10000}\), algebraic
denominators, and zero/nonzero left-kernel residuals pass. A repaired 2+2
fifth-band replay is active. The solver also reuses the already computed
transposed RREF to build the left-kernel basis and carry the rank into
compatibility checking; the exact quadratic-field smoke check verifies this
reconstruction. In the repaired run, the (2+2) 48/58 block has now passed
exactly with rank (61), (91) left-kernel equations, and zero residuals,
leaving a (7)-dimensional affine fibre. The run has advanced to the 47/57
block: its exact pivot profile has rank (61) and the solve is underway, but
no compatibility result for that block is available yet.

The next corrected band pair, weights \(45/55\), includes the weight-70
resonance correction \(-3\kappa P^3Q\). The exact convolution formulas
for these rows and the corresponding binomial coefficient \(q_{20}\) pass
the dedicated formula audit in
[f2_75_125_power_relation_fifthbands.py](../src/f2_75_125_power_relation_fifthbands.py).
The exact \(2+1+1\) calculation at weights \(45/55\) has completed: the
combined 152-row, 72-variable matrix has rank \(63\), its left kernel has
dimension \(89\), and all compatibility residuals vanish. The affine fibre
therefore has dimension \(9\), and the carried family has 32 parameters.
Every point left by the preceding bands extends through these equations, so
this corrected band pair adds no branch restriction and does not exclude the
carrier. The 2+2 calculation through these weights remains to be completed
with the repaired parameter-polynomial checker.

### Magnus IV gives a finite inner-polynomial split

Assuming the standardized F2 outer support envelope, a generic common source
translation makes both Newton polygons exact rectangles
\(\operatorname{Rect}_{15,60}\) and \(\operatorname{Rect}_{25,100}\).
Lee–Li's proved Magnus-IV bound then restricts the northeast vertex of the
inner polynomial to 12 lattice points, including the origin. The origin
would make the first polynomial a nontrivial composition, impossible for a
nonzero constant-Jacobian pair; hence there are 11 possible nonconstant
inner-vertex branches. Their union implies the necessary support condition
\(\operatorname{Supp}(\zeta(P))\subseteq[0,9]\times[0,36]\) after specializing
the polynomial pre-generator recursion back to the original coordinates.
Equivalently, there must be a pre-generator \(R\) with
\(P=R^3+Z\), \(\operatorname{Supp}(R)\subseteq[0,5]\times[0,20]\), and
\(\operatorname{Supp}(Z)\subseteq[0,9]\times[0,36]\). This is a further
necessary condition, not an exclusion. A new exact branchwise coefficient
audit uses \(R^3=P-Z\) outside the union of the F2 support polygon and each
Magnus rectangle. Across all 11 branches, those equations force
\(\operatorname{Supp}(R)\subseteq\{5i-j\le5\}\), a 76-monomial region;
then \(\operatorname{Supp}(Z)\) shrinks to the intersection of its branch
rectangle with the F2 support polygon. The smallest and largest resulting
models have 1,659 and 1,916 coefficient variables before the Jacobian
equations. The translation and specialization argument, exact branch counts,
and source details are recorded in
[magnus_75_125_application.md](progress/magnus_75_125_application.md); the support
deduction is replayed by
[magnus_75_125_cube_support.py](../src/magnus_75_125_cube_support.py).

The top y-slice relation \(q_{100}=b r_{20}^5\) and a Euclidean expansion of
\(Q\) in powers of this pre-generator further give
\(Q=bR^5+cR^4+S\), with \(\deg_yS\le60\). This cuts each raw branch model
from 1,576 Q coefficients to 708 tail coefficients and two scalars. In the
four branches with \(\deg_yZ=24,28,32,36\), the leading coefficient
\(r_{20}(x)\) must be a fifth power of a linear polynomial. The first three
upper-y rows of \([Z,R]\) then force \(x^3\mid r_{19}\) in the first three
of these four branches and \(x^2\mid r_{19}\) in the fourth. The
resonant fourth row strengthens the last three to \(x^3\mid r_{19}\) and
\(x\mid r_{18}\); in the \((7,28)\) branch it also factors as
\((5b_1-a_3^2)(5b_1-2a_3^2)=0\). The \((6,24)\) branch has the separate
split \([x^3]r_{19}\,r_{18}(0)=0\). The congruence
\(R^2\mid1-[Z,S]\) forces \(x^{11-e/4}\mid[y^{60}]S\), and its next
remainder row forces \(x^{9-e/4}\mid[y^{59}]S\) for the first three
branches. The resulting high-branch leaves have 907/906, 944/943, 986, and
1030 variables (including the indicated intersections). The derivations and
branch counts are written in
[magnus_75_125_application.md](progress/magnus_75_125_application.md) and replayed by
[magnus_75_125_q_tail.py](../src/magnus_75_125_q_tail.py) and
[magnus_75_125_commutator_front.py](../src/magnus_75_125_commutator_front.py) and
[magnus_75_125_tail_remainder.py](../src/magnus_75_125_tail_remainder.py) and
[magnus_75_125_tail_remainder_next_slice.py](../src/magnus_75_125_tail_remainder_next_slice.py).
These are still necessary conditions only; the remaining Jacobian equations
and exact certificates are outstanding.

Continuing the quotient recursion one row further, the \(y^{e+57}\) remainder
has fixed \(x\)-coefficients below the first possible contribution from
\(S_{58}\). Write \(r_{19}=x^3(a_3+a_4x)\),
\(r_{18}=\rho_0+\rho_1x+\rho_2x^2+\cdots\), and
\(S_{60}=x^{11-s}(u_0+u_1x+\cdots)\). For \(s=7,8,9\),
\(\rho_0=0\) is already known. The fixed equations are
\[
s=8:\quad u_0(2529a_3^2-15125\rho_1)=0,
\qquad
s=9:\quad
\begin{cases}
u_0(59a_3^2-240\rho_1)=0,\\
55a_3^2u_1+114a_3a_4u_0-235\rho_1u_1-190\rho_2u_0=0.
\end{cases}
\]
There are no such fixed coefficients for \(s=6,7\). The equations retain the
zero cases \(u_0=0\) and \(a_3=0\); no coefficient is divided out. Exact direct
bracket and remainder replay is in
[magnus_75_125_tail_remainder_lowrows.py](../src/magnus_75_125_tail_remainder_lowrows.py).
Resolving the pivots gives \(s=8\): either \(u_0=0\), or
\(u_0\ne0\) and \(\rho_1=2529a_3^2/15125\). For \(s=9\), the
\(u_0\ne0\) branch pins \(\rho_1=59a_3^2/240\) and
\(\rho_2=57a_3a_4/95-133a_3^2u_1/(9120u_0)\); if \(u_0=0\), then
either \(u_1=0\) with no further condition from these rows, or
\(u_1\ne0\) and \(\rho_1=11a_3^2/47\). These extra selectors have not
yet been combined into an exclusion. One quotient row later, after solving
the necessary \(S_{58}\) coefficients, the \(s=8\) branch with \(u_0\ne0\)
forces \(a_3=\rho_1=\sigma_0=0\) and \(\sigma_1=2a_4\rho_2/5\); the
\(s=9\), \(u_0\ne0\) branch forces
\(a_3=\rho_1=\rho_2=\sigma_0=\sigma_1=0\) and pins
\(\sigma_2=3a_4\rho_3/5-4a_4^3/25\). The \(s=9\),
\(u_0=0,u_1\ne0\) branch forces \(a_3=\rho_1=\sigma_0=0\). The
remaining vanishing-\(S_{60}\) leaves carry exact additional equations; none
has been excluded. The full branch equations and replay are in
[magnus_75_125_application.md](progress/magnus_75_125_application.md) and
[magnus_75_125_tail_remainder_depth4.py](../src/magnus_75_125_tail_remainder_depth4.py).

## Translating Helali's certificate pattern to this pair

### Coupling the next pair of bands still leaves every edge parameter

The first companion commutator was redundant, so I checked one level deeper: weight 48 of \([P,W^*]\) together with weight 58 of \([Q,W^*]\). Both equations are linear in the same 60 coefficients of \(W^*_{37}\) and eight coefficients of \(U_2=hB\), once the weight-49 kernel has parametrized \(S\) and \(f=W^*_{38}\). Their forcing terms are at most quadratic in the weight-49 parameters. The weight-58 equation used here is

\[
5U^4(5Uk'-37U'k)+24q_{24}f'-38q_{24}'f
+23q_{23}g'-39q_{23}'g-20U^9A_2-30U^6A^2=0.
\]

The exact combined matrix has 68 new variables and rank 61 in the \(2+2\) branch over \(\mathbb Q(\sqrt5)\), and rank 62 in the \(2+1+1\) branch over \(\mathbb Q\). Every left-kernel forcing polynomial vanishes identically on the full weight-49 kernel in both cases. Thus every point of the 3-dimensional (respectively 2-dimensional) weight-49 parameter space extends through this pair of rows; the new solution fibres have dimensions 7 and 6. The replay is [f2_75_125_power_relation_pq_nextbands.py](../src/f2_75_125_power_relation_pq_nextbands.py). This rules out the hoped-for rigidity at these two adjacent bands in the relaxed edge model; the next useful attack must bring in later bands or equations from the wider support system.

Helali's useful computational lesson is to postpone a large Gröbner-basis computation until the geometry has reduced the unknowns. His released replay splits one configuration into explicit sign branches, treats a vanishing specialization separately, and closes the terminal systems with exact polynomial identities or unit-ideal certificates. The release also states its boundary clearly: the coefficient systems are excluded, while the implication from the published reduction to the degree-pair conclusion remains conditional on that reduction and its transcription ([paper and replay](https://github.com/bilLkarkariy/jc2-72-108-exact-certificates)).

For \((75,125)\), the analogous order of attack is:

1. Keep the exhaustive outer support envelope as the theorem-facing model; do not assume the narrower normalized carrier equations until their attachment to an actual hypothetical pair is proved.
2. Eliminate linear slices first. The upper-edge equations expose \(p_{60}=aR^3\), \(q_{100}=bR^5\), then reduce each new slice to the one-variable operator \(\mathcal D_s\). The first four slices now force the root-jet conditions above.
3. Split by the multiplicity partition of \(R\), and by any pivot that is divided out. Each zero-pivot locus must get its own equations or certificate; saturation alone must not silently drop it.
4. Use only the remaining necessary equations to form a small residual ideal. If that relaxed ideal is the unit ideal, it excludes the full system. Otherwise retain its exact branches and continue to later bands.
5. End every proposed exclusion with a characteristic-zero identity and a separate checker that verifies the identity and its branch hypotheses. Modular ranks can discover a certificate, but by themselves do not prove a characteristic-zero unit ideal.

The present jet conditions remain pruning rules: their forced divisors fit within the available slice degrees. The next concrete step is to substitute the recursively constrained depths 1–60 into the fourteen exact remainder conditions at depth 61, while carrying the separate resonance branches at depths 64–100. Those constraints can then be coupled to the low-y equations before any large Gröbner computation.

## What still has to be proved

1. Independently referee the minimality-free front-end extraction and its cited normal-form inputs. Its finite-witness custody checker has been replayed, but that checker does not prove the GGV/Makar-Limanov arguments.
2. Eliminate the remaining bands of the 2,162-variable outer system without dividing by an unproved nonzero coefficient; branch explicitly at vanishing pivots.
3. For each surviving algebraic branch, produce exact characteristic-zero certificates (for example, a unit-ideal identity or a polynomial identity forcing a parameter to vanish on a nonzero branch) and an independent replay checker.
4. If pursuing the additional normalized carrier system instead of the outer envelope, prove its full source/support attachment. The GGV successor condition now attaches the exact double root of \(h\), but that alone does not establish every normalization or equation used in the narrower carrier system.

  The GGV table and the rectangle/corner consequences above are published inputs with an explicit derivation. The front-end reduction is supported by the written Lemma A extraction and its finite-witness transfer; the accompanying checker validates custody and finite arithmetic, not the theorem proofs. The outer envelope is broader than the narrower normalized carrier model used in exploratory work. It now supplies a concrete finite system for every hypothetical pair after standardization, but the system has not been eliminated and the pair remains open.

## Sources

- GGV, [*Some algorithms related to the Jacobian Conjecture*](https://arxiv.org/html/1708.07936), especially the \(F_2\) family table, Theorem 2.20's type-II.b successor formula, and the \((75,125)\) row.
- Valqui et al., [*On the shape of possible counterexamples to the Jacobian Conjecture*](https://arxiv.org/abs/1401.1784), especially Definition 4.3 and Proposition 4.7 on standard pairs and subrectangular support.
- Makar-Limanov, [*On the shape of a counterexample to the two-dimensional Jacobian conjecture*](https://doi.org/10.55630/serdica.2025.51.299-314), pp. 304–305, for the complex-field Newton-polygon rectangle.
- Exploratory [F2 upstream carrier/support audit](https://github.com/royvanrijn/jacobian-research/blob/main/plane-jc/F2_UPSTREAM_CARRIER_EXTRACTION_PROFILE.md), which records the support envelope and its trust boundary.
- Exploratory [F2 \((75,125)\) derivation audit](https://github.com/royvanrijn/jacobian-research/blob/main/plane-jc/F2_75_125_DERIVATION.md), which labels its result a forced skeleton rather than a normal-form certificate.

## Terminal F2 chart and exact sparse replay

The F2 complete-chain attachment gives a terminal type-I.b corner after at
most the two removable type-III shears of slopes (2/5) and (3/5). In the
normalized ramified chart the remaining source change is
(u\mapsto u+a w+b w^2), and the terminal halfspace requires
([w^s u^j]F=0) when (5s+12j<72), with the analogous (120) bound for
(G). At the distinguished double root this first forces (S(1)=0), hence
excludes the tuned nonzero (S(1)) face. The zero cases remain in the system.
The source derivation and its scope are recorded in
[f2_terminal_chain_local_selectors_20261018.md](progress/f2_terminal_chain_local_selectors_20261018.md).

There are exactly 32 forbidden (P)-positions through band distance five:
the slice counts are (6,6,6,5,5,4). The centering coefficients (a,b)
are fixed by the (u^5w) and (u^5w^2) rows; the other rows are then imposed
as equations. The independent generic Taylor checker verified all 32
positions and the source-translation covariance in
[f2_terminal_centering_identities_20261018.py](../src/f2_terminal_centering_identities_20261018.py).
The exact sparse solver retains the full (P_{10}=t^2U_5) support.

In the (2+2) carrier, the sparse lift through distance five verified every
coupled coefficient row and every terminal (P)-row with zero exact residual.
An independent local (Q)-chart check found all 56 forbidden (Q)-positions
through the same distance vanish identically on that family. Thus these
terminal conditions have not yet excluded the two-factor carrier.

For the (2+1+1) control carrier, a stronger exact lift includes the root
conditions \(\operatorname{rad}(h)\mid(3B-S^2)\) and
\(\operatorname{rad}(h)\mid(27U_3-S^3)\), together with the 32 terminal
rows, through distance seven. All 152-row coupled blocks and all terminal
rows through that distance verify with zero residual. The (Q)-chart check
then verifies 72 forbidden positions through distance seven. The exact
replay receipts are
[the sparse band receipt](../data/json/f2_75_125_sparse_ring_lift_211_terminal_roots7_resume.json),
[the Q-chart receipt](../data/json/f2_75_125_sparse_terminal_q_211_terminal_roots7.json),
and [the corrected power-law audit](../src/f2_75_125_fractional_power_band_audit.py).
These are necessary-condition checks, not an exclusion. The (2+2) lift
with both root conditions is being extended separately; the saved
through-five result without those root conditions is not substituted for it.

Extending only the (2+1+1) control to distance nine exposed 40 coefficient
equations at each of distances eight and nine from polynomiality of the
corrected power relation, and 61 bracket-compatibility equations at each
distance. The original (Q)-support rectangle adds four equations at each
of those distances; the terminal (Q)-halfspace adds seven at each. The
retained system has 224 distinct equations in 18 parameters, all homogeneous
of degree eight or nine, so the origin satisfies this truncated residual
ideal. It is therefore not a unit-ideal certificate. The exact state and
profile are [here](../data/json/f2_75_125_sparse_ring_lift_211_terminal_roots9_full_state.json)
and [here](../data/json/f2_75_125_sparse_ring_lift_211_terminal_roots9_ideal_profile.json).

## Exact Gröbner split of the `2+2` two-root slice

After imposing the exact second-root condition
`p2=((sqrt(5)-5)/5)*p1`, the depth-seven `2+2` necessary system has 123
constraints in 16 active parameters over `Q(sqrt(5))`. Singular 4.4.1 computes
its standard basis with one generator, and that generator is not a unit. Its
factorization is a nonzero rational scalar times a weighted-degree-three
factor and a weighted-degree-four factor, each once. The weights are 1 on
`p1`, 2 on `p7,p8,p9`, 3 on `p16,...,p20`, and 4 on `p27,...,p33`.
The factor product is verified exactly inside Singular. Each of the 123 input generators reduces to the one-element basis, and an explicit lift reconstructs that basis generator from the inputs with zero remainder. This verifies equality of the principal ideal in both directions. This converts the
slice into the two algebraic cases `F3=0` and `F4=0` for the next terminal
rows. It does not exclude either case or the whole `2+2` branch. The full
input, exact run log, and hash receipt are listed in
[f2_75_125_branch_coverage_status_20261018.md](progress/branch-coverage-ledger.md).

Each component now also has an exact global pivot chart: `F3=0` solves for
`p16` with a nonzero constant coefficient in `Q(sqrt(5))`, and `F4=0` solves
for `p27` the same way. Exact substitutions verify both identities. The
receipt and replay script are
[f2_75_125_Salpha_factor_charts_20261018.json](../data/json/f2_75_125_Salpha_factor_charts_20261018.json)
and
[f2_75_125_Salpha_factor_charts_20261018.py](../src/f2_75_125_Salpha_factor_charts_20261018.py).
These charts facilitate extending the necessary equations separately; neither
component is yet excluded.

The two components are affine 15-spaces and meet in an affine 14-space. The
factor multiplicities are both one, so the exact slice ideal is radical. The
intersection pivot remains nonzero after solving the cubic component first;
this is recorded in the chart receipt. This makes the decomposition exhaustive
without any missing pivot-open strata.

## K9 local pruning for nonzero `U1` at the `2+1+1` double root

At the distinguished double root `t=1`, the F2 terminal condition gives
`ord(U1)>=5`; K7/K8 put the nonzero-`U1` deeper branch in
`ord(U2)>=3`, `ord(U3)>=1`, with `ord(U3)>=2` if `ord(U2)=3`. Write
`ord(U1)=k`, `ord(U2)=q`, `ord(U3)=r`, and `ord(U4)=s`. The K9 first-integral
correction is
\[
\frac{\kappa_1}{t}[w^4]F^{4/3}.
\]
Since `t` is a unit at this root, the five weighted partitions of four show
that this correction has order at least two on every permitted deeper
stratum: their lower bounds are respectively `2+s`, `k+r-4`, `2q-4`,
`2k+q-10`, and `4k-16`. It therefore cannot cancel a K9 negative pole.

The exact checker
[`f2_75_125_k9_nonzero_first_slice.py`](../src/f2_75_125_k9_nonzero_first_slice.py)
enumerates all 51,870 supported order profiles of `(U1,...,U5)` in this
stratum, treats `U6,...,U9` at order zero, and certifies only profiles whose
unique least negative term uses `U1,...,U5`. The support bounds
`U1=(t-1)^5R^2B`, `deg B<=2`, and `U2=(t-1)^2RC`, `deg C<=7`, give
`ord(U1)<=7` and `ord(U2)<=9` at `t=1`. The checker finds 16,200 certified
profiles (31.2319% of this finite order-profile set, **not** of the whole
branch).
Because all terms involving `U6,...,U9` already have order at least the
unique minimum in the baseline, raising any later order cannot remove these
certificates. Two sample obstructions are:

* `k=5,q=3,r>=3,s>=1`: the unique lowest term is
  `-(35/729)U1*U2^4/h^10`, of order `-3`.
* `k=5,q=4,r=1,s=0,ord(U5/t)=0`: the unique lowest term is
  `-(5/81)U3^3/h^4`, of order `-5`.

The checker also verifies that the K9 correction has minimum order two on the
entire enumerated stratum. This is new local pruning of the nonzero-`U1`
deeper double-root branch, but other profiles remain and it does not exclude
the `2+1+1` carrier or the degree pair.

## K10 local pruning for nonzero `U1` at the `2+1+1` double root

The first-integral identity gives
\[
G_{10}=[w^{10}]F^{5/3}
+\frac{\kappa_1}{t}[w^5]F^{4/3}
+\frac{\kappa_2}{t^2}h^3.
\]
At `t=1`, the exact seven-partition valuation screen finds
`ord([w^5]F^(4/3))>=1` on every supported nonzero-`U1` profile; the final
resonance term has order six because `ord(h)=2`. Thus neither resonance can
cancel a negative pole in the distance-ten binomial coefficient.

The exact checker
[`f2_75_125_k10_nonzero_first_slice.py`](../src/f2_75_125_k10_nonzero_first_slice.py)
enumerates the same 51,870 supported `(U1,...,U5)` order profiles as K9, with
`U6,...,U10` set to the worst-case order zero. It counts only a unique least
negative term supported on `U1,...,U5`, so increasing later orders cannot
remove the obstruction. K10 has 16,994 such profiles; 3,420 are new beyond
the 16,200 K9 profiles. The K9/K10 union is 19,620 / 51,870 = 37.8253% of
this finite order-profile set, **not** of the `2+1+1` branch.

For example, `ord(U1)=5`, `ord(U2)=3`, `ord(U3)=2`, `ord(U4)=1` gives the
unique K10 term `-(7/729)U2^5/h^10`, of pole order five. This profile is not
removed by the robust K9 screen, while the K10 resonance correction has order
at least one. The new exact screen reduces local profiles only; the whole
`2+1+1` branch remains unexcluded.

## K11--K12 local pruning for nonzero `U1` at the `2+1+1` double root

The next first-integral bands are
\[
G_{11}=[w^{11}]F^{5/3}
+\frac{\kappa_1}{t}[w^6]F^{4/3}
+\frac{\kappa_2}{t^2}U_1,
\qquad
G_{12}=[w^{12}]F^{5/3}
+\frac{\kappa_1}{t}[w^7]F^{4/3}
+\frac{\kappa_2}{t^2}U_2.
\]
The `kappa_2` terms are regular at the distinguished double root. The
`kappa_1` corrections are not uniformly regular: their minimum orders across
the supported profiles are `-2` at K11 and `-3` at K12. The exact screen
therefore counts a main-term pole only when it is unique, uses no later slice
than `U5`, and has order strictly below the full `kappa_1` correction.

The exact-rational checker
[`f2_75_125_k11_k12_nonzero_first_slice.py`](../src/f2_75_125_k11_k12_nonzero_first_slice.py)
also replays the K9/K10 profile totals and requires them to agree with the
saved K10 receipt. On the 51,870 supported profiles, it finds 14,430 robust
K11 poles (514 new after K9/K10) and 20,176 robust K12 poles (1,388 new after
K9/K11). Together with the previous bands this gives 21,522 / 51,870 =
41.4922% of the finite local order-profile set, **not** of the whole branch.

Two new examples are:

* K11, `ord(U1)=6`, `ord(U2)=4`, `ord(U3)=2`, `ord(U4)=0`: the unique main
  term is `-(5/27)U3*U4^2/h^4`, of pole order six, below the K11 correction
  order zero.
* K12, `ord(U1)=5`, `ord(U2)=5`, `ord(U3)=3`, `ord(U4)=1`,
  `ord(U5/t)=0`: the unique main term is `-(5/81)U4^3/h^4`, of pole order
  five, below the K12 correction order zero.

These additional exact local obstructions still do not eliminate the full
`2+1+1` branch.

## K13--K14 local pruning for nonzero `U1` at the `2+1+1` double root

At the next two bands, the first-integral expansions are
\[
G_{13}=[w^{13}]F^{5/3}
+\frac{\kappa_1}{t}[w^8]F^{4/3}
+\frac{\kappa_2}{t^2}U_3,
\qquad
G_{14}=[w^{14}]F^{5/3}
+\frac{\kappa_1}{t}[w^9]F^{4/3}
+\frac{\kappa_2}{t^2}U_4.
\]
The `kappa_2` terms have orders at least one and zero, respectively. The
`kappa_1` terms can have negative order, so the checker compares each
unique main pole against the full correction order profile by profile.

The exact-rational checker
[`f2_75_125_k13_k14_nonzero_first_slice.py`](../src/f2_75_125_k13_k14_nonzero_first_slice.py)
replays the K9--K12 exclusion flags for all 51,870 supported profiles and
checks their saved per-band totals and union before counting new profiles.
It finds 228 new K13 profiles and 80 new K14 profiles. The cumulative union is
21,830 / 51,870 = 42.086% of this finite order-profile set, **not** of the
whole `2+1+1` branch. On the profiles surviving K9--K12, the minimum
`kappa_1` correction orders are `-4` at K13 and `-5` at K14.

Examples are:

* K13, `(ord(U1),ord(U2),ord(U3),ord(U4),ord(U5/t))=(6,5,3,3,1)`: the unique
  term `-(5/27)U3*U5^2/h^4` has pole order three, below correction order zero.
* K14, `(6,6,5,3,1)`: the unique term `-(5/27)U4*U5^2/h^4` also has pole
  order three, below correction order zero.

This remains local profile pruning; the full `2+1+1` branch is not excluded.

## K15--K16 local pruning for nonzero `U1` (2026-10-21)

The next first-integral coefficients are
\[
\begin{aligned}
G_{15}&=[w^{15}]F^{5/3}
 +\frac{\kappa_1}{t}[w^{10}]F^{4/3}
 +\frac{\kappa_2}{t}U_5+\kappa_3h^2,\\
G_{16}&=[w^{16}]F^{5/3}
 +\frac{\kappa_1}{t}[w^{11}]F^{4/3}
 +\frac{\kappa_2}{t}U_6+\kappa_3[w]F^{2/3}.
\end{aligned}
\]
At the distinguished double root, the additional resonance orders are at least
zero at K15 and at least zero at K16 (`ord(U1/h) >= 3` for the last K16
term). The `kappa_1` correction can have a pole, so each main pole is compared
against the full correction order for its profile.

The exact-rational checker
[`f2_75_125_k15_k16_nonzero_first_slice.py`](../src/f2_75_125_k15_k16_nonzero_first_slice.py)
replays K9--K14 profile by profile and checks the saved per-band counts and
union before adding K15/K16. Among 51,870 profiles, K15 adds 1,000 new
profiles; K16 adds none after K9--K15. The cumulative K9--K16 union is
22,830 / 51,870 = 44.0139% of this finite local profile grid, **not** of the
whole `2+1+1` branch. There are 29,040 grid profiles left after this screen.

For example, `(ord(U1),ord(U2),ord(U3),ord(U4),ord(U5/t))=(6,5,4,3,1)`
gives the unique K15 term `-(5/81)U5^3`, with pole order five and full
correction order `-2`. On the profiles reaching K15, the minimum
`kappa_1` correction order is `-7`; among the profiles reaching K16 it is
`-8`. These are local profile results only; the `2+1+1` branch remains open.

## K25–K29 local valuation screen

The exact rational checker [f2_75_125_k25_k29_nonzero_first_slice.py](../src/f2_75_125_k25_k29_nonzero_first_slice.py) continues the finite profile analysis from K24. It evaluates the exact first-integral terms through K29 and conservatively retains the possible order-zero `kappa_5` correction at K25; the next resonance, `kappa_6`, starts at K30. Across all 29,040 profiles surviving K24, no new unique main pole lies strictly below every correction at K25, K26, K27, K28, or K29. The receipt is [f2_75_125_k25_k29_nonzero_first_slice.json](../data/json/f2_75_125_k25_k29_nonzero_first_slice.json).

The finite-grid union therefore remains 22,830 / 51,870 = 44.013881%, with 29,040 profiles still surviving this screen. This gives **0% whole-branch exclusion** for the `2+1+1` branch; it only shows that this particular first-integral valuation screen adds nothing through K29. The global exact coefficient lift is still the active route for this branch.

## Leading-form factor strata for `(75,125)`

For a hypothetical characteristic-zero Keller pair with degrees `(75,125)`, write the ordinary degree forms as `P_75 = a H^3`, `Q_125 = b H^5`, `deg(H)=25`. The number of distinct linear factors of `H` over the algebraic closure gives three exhaustive strata: one, exactly two, or at least three. Abhyankar's two-point theorem excludes the last stratum. The one-factor case is not excluded by the universally quantified one-point criterion: Abhyankar's Theorem (19.4) says that the implication from *every* Jacobian pair having one point at infinity to the Jacobian conjecture is equivalent to the conjecture. That does not establish that a particular pair with a one-factor leading form is automorphic. Jung--van der Kulk divisibility likewise applies after automorphy has been proved. Consequently, this ordinary top-form argument leaves both the one-factor and exactly-two-factor strata open. The exactly-two-factor stratum is separate from the conditional F2 front-end reduction to a standardized `5+20` carrier; the reduction's hypotheses and attachment must be audited separately.

`H=x^25+y^25` has 25 distinct linear factors over `C`, so it lies in the excluded three-or-more-factor stratum, not in the one- or two-factor stratum. The pair `(H^3,H^5)` has zero Jacobian and represents leading forms only; it is not a Keller pair. Any computational exclusion of a factorization stratum must quantify over its entire lower-coefficient family. References: [Abhyankar's two-point theorem and global equivalent formulations](https://mathweb.tifr.res.in/Documents/Publications/Lectures/tifr57.pdf); [Abhyankar's exposition](https://www.indianmathsoc.org/ms/mathstudent2007.pdf).

## K30–K35 valuation screen

The exact local screen adds no profiles at K30–K35 after retaining `kappa_6` and `kappa_7`: 29,040 / 51,870 finite local profiles survive; the cumulative K9–K35 screen remains 44.013881% of the finite grid only. Receipt and SHA-256 are recorded in `f2_75_125_branch_coverage_status_20261018.md`. Whole-branch coverage is unchanged; the exact depth-16 FLINT lift remains the active computation.

## K36–K49 local screen update

The exact-rational valuation screen continues through K49 and adds no local profiles at K36–K45 or K46–K49. It remains 22,830 / 51,870 = 44.013881% of the finite local-profile grid, with 29,040 profiles surviving; this does not imply any whole-branch exclusion. Receipts: `f2_75_125_211_k36_k45_nonzero_first_slice_20261022.json` and `f2_75_125_211_k46_k49_nonzero_first_slice_20261022.json`.

A streamed-target Intel Iris Xe weight-26 Macaulay scan is now running on the active `2+1+1` support-23 slice, alongside the exact depth-16 lift. It has not yet entered device elimination or produced a receipt at the last poll. The scan is over GF(251), so it is only a candidate-discovery pass; no modular rank can be promoted to characteristic-zero branch exclusion without an exact-Q certificate.

## Live follow-up at K50 and GPU weight 26

As recorded in `f2_75_125_branch_coverage_status_20261018.md`, a streamed-target OpenCL Macaulay run at weighted degree 26 is active on the Intel Iris Xe; it has reached rank 150 after scanning 7,692 of 316,694 columns. This is a modular support-23 screen and currently contributes no branch exclusion. The exact local valuation continuation through K50–K54 is being computed with the full new r=10 correction included from K50 onward; the first band is not complete yet. No branch-wide percentage has changed.

## Live poll at W26 rank 200 — 2026-10-23 22:57 +06:00
The streamed-target GF(251) W26 scan has advanced to rank 200/4,304 after 10,747/316,694 columns, but remains modular reconnaissance without a receipt. The exact depth-16 lift and K50–K54 valuation screen are still running; the latter has not completed K50. Physical RAM free is about 730 MiB and commit headroom about 13.3 GiB. No branch-wide exclusion has been added.

## K50 local-profile screen — 2026-10-23
The exact rational first-integral profile screen completed K50 and excluded zero additional profiles. It retains 29,040 of 51,870 finite order profiles; K9–K50 cumulative exclusion remains 22,830/51,870 = 44.013881% of that finite grid only. K51–K54 are still running. Receipt: `f2_75_125_211_k50_k54_nonzero_first_slice_20261023.json`, SHA-256 `8b19115511a9979e68a43e8164fbfbc6fbfe191b69af1cf738b192eb67fc28e4`. This yields 0% whole-branch exclusion for `2+1+1`; `2+2` stays paused at 0%.

## K51 live continuation — 2026-10-23

K51 finished with zero new finite order-profile exclusions; 29,040/51,870 profiles survive, leaving the finite-grid tally at 44.013881% and whole-branch 2+1+1 exclusion at 0%. K52–K54 remain live. The exact depth-16 lift remains live without a receipt; the W26 GPU scan reached rank 550/4,304 after 16,883 columns and remains modular reconnaissance. The separate 2+2 branch remains paused at 0% whole-branch exclusion. The Appelgate–Onishi prime-factor criterion and Magnus/Nagata gcd criteria do not directly apply to (75,125).

## Resource reprioritization — 2026-10-23

The K50–K54 finite-profile worker is paused after K51 (K51 added zero profiles; the K51 receipt is saved and hash-verified). The W26 OpenCL scan was stopped at rank 600/4,304 and 17,828 columns with no receipt; its per-column blocking readback made the modular screen too slow under current memory pressure. The exact depth-16 FLINT lift remains active. Whole-branch exclusion is still 0% in the active 2+1+1 branch and 0% in paused 2+2.

### K52–K125 exact local-profile continuation and strict five-variable GB result — 2026-10-24 15:57 +06:00

The OpenCL integer-DP continuation is complete through K125 for the finite F2 2+1+1 local valuation grid. Exact CPU/OpenCL parity was checked on 1,771,440 profile-band outputs. K65–K125 removes zero new profiles; cumulative pruning remains 22,830/51,870 = 44.013881% of this finite grid, with 29,040 profiles surviving. This is local finite-grid pruning only: whole-branch exclusion remains **0% for 2+1+1** and **0% for paused 2+2**.

The corrected strict-support d13 Groebner calculation is finished. It keeps only complete equations supported in {p1,p7,p15,p16,p17} after the valid p79=0 specialization: 112 equations, 105 distinct up to rational scale. The exact F5B/Grevlex basis has 48 elements, is nonunit, and is positive-dimensional. Therefore this necessary subsystem does not exclude the chart. Report: f2_75_125_sparse_ring_lift_211_roots17_dd_20261024_d13_strict5_gb_cached_gmpy2_ringapi.json, SHA-256 4D544110D025A1F10FC05747669B7E11369B2AC0ED690E1D6EA310428EECF998.

A Rabinowitsch saturation check for the open condition p1 != 0 is now running on that exact basis as a diagnostic; no branch nonvanishing theorem for p1 has been established, so even a forced p1=0 would need a separate branch attachment before it could exclude anything.

## Geometric-degree and generic-fiber topology filter — 2026-10-26

Let delta be the geometric degree (number of points in a generic fiber) of a hypothetical Keller map with component degrees 75 and 125. Nguyen Van Chau's Theorem B, applied after swapping the components so that 125=25*5 is the larger degree, gives delta=5r+3s for nonnegative integers r,s with r+s>=1 and delta>=min(6,5)=5. Zoladek's theorem excludes geometric degree at most 5, so a nonautomorphism has delta>=6; the semigroup condition also rules out delta=7.

Chau's Theorem C gives, in the original component labels, 5*chi_P - 3*chi_Q = 2*delta. The same paper's puncture formula gives, for a nonautomorphism, at least 4 punctures in a generic P-fiber and at least 6 in a generic Q-fiber: the reduced degree factors are 3 and 5. These are branches of the generic affine fibers, not the two distinct points at infinity determined by H.

These conditions do not yet exclude the F2 2+1+1 or 2+2 quartic-root branches. They are compatible, for example, with delta=6=3+3, a rational generic P-fiber with 5 punctures (chi_P=-3), and a rational generic Q-fiber with 11 punctures (chi_Q=-9), since 5*(-3)-3*(-9)=12=2*6. A useful next bridge would express Chau's separative-series data and puncture counts in terms of the exact F2 quartic-root/lower-band coefficients; that attachment has not been proved. No branch coverage percentage changes. Sources: Nguyen Van Chau, https://matwbn.icm.edu.pl/ksiazki/apm/apm71/apm7135.pdf; Zoladek, https://doi.org/10.1016/j.top.2008.04.001.

Clarification: the genus/puncture data above is only a numerical compatibility witness for the listed necessary constraints; it is not asserted to come from any Keller map.

