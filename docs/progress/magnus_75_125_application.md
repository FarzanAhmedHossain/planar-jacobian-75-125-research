# Applegateâ€“Onishiâ€“Nagata and Magnus tools for the \(75,125\) case

## Arithmetic shortcuts do not cover this pair

The plane results associated with Magnus, Applegateâ€“Onishi, and Nagata exclude
several arithmetic classes of degree pairs. They do not directly cover
\((75,125)\):

\[
75=3\cdot5\cdot5,\qquad 125=5\cdot5\cdot5,
\qquad \gcd(75,125)=25.
\]

Each component degree has three prime factors counted with multiplicity, and
the gcd is neither prime nor at most 8. The Applegateâ€“Onishi degree criterion
and the Magnusâ€“Applegateâ€“Onishiâ€“Nagata gcd criterion therefore cannot by
themselves exclude this case. Their structural lemmas and Magnus's homogeneous
component formula remain useful inputs to a degree-specific reduction.

## Magnus IV gives a finite inner-vertex split

The F2 outer model has

\[
\operatorname{Supp}(P)\subseteq[0,15]\times[0,60],\qquad
\operatorname{Supp}(Q)\subseteq[0,25]\times[0,100],
\]

with nonzero upper-right coefficients at \((15,60)\) and \((25,100)\).
To apply Leeâ€“Li's Theorem 5.8 and Corollary 5.10, first make a generic source
translation \(P_{c,d}(x,y)=P(x+c,y+d)\), \(Q_{c,d}(x,y)=Q(x+c,y+d)\).
The top-right coefficients stay nonzero. For \(P\), the translated
coefficient of \(x^{15}y^0\) is a nonzero polynomial in \(d\), because its
leading coefficient is the original coefficient of \(x^{15}y^{60}\). The
coefficient of \(x^0y^{60}\) is likewise a nonzero polynomial in \(c\).
The same statements hold for \(Q\). Thus the same generic \(c,d\) can be
chosen so that
\(N^0(P_{c,d})=\operatorname{Rect}_{15,60}\) and
\(N^0(Q_{c,d})=\operatorname{Rect}_{25,100}\). Scaling the two
components normalizes the top corners and preserves a nonzero constant
Jacobian.

The Magnus-IV hypotheses then hold with
\[
(a,b,m,n)=(3,5,15,60),\qquad bm/a=25,\qquad bn/a=100.
\]
For the inner polynomial \(\zeta(P_{c,d})\), Corollary 5.10 places its
northeast vertex in

\[
0\le y<40,\qquad \frac{53x-15}{13}\le y\le4x.
\]

An exact rational enumeration gives

\[
\begin{aligned}
&(0,0),\ (1,3),\ (1,4),\ (2,7),\ (2,8),\ (3,12),\\
&(4,16),\ (5,20),\ (6,24),\ (7,28),\ (8,32),\ (9,36).
\end{aligned}
\]

The origin case would make \(\zeta\) constant or zero, hence
\(P_{c,d}=R_{c,d}^{3}+c_0\) for the pre-generator \(R_{c,d}\).
That is a nontrivial one-variable composition, which is incompatible with a
nonzero constant Jacobian by the standard composition lemma. Thus only the
other 11 points remain. In particular,

\[
\operatorname{Supp}(\zeta(P_{c,d}))\subseteq[0,9]\times[0,36].
\]

This is a finite branch split: in each branch the northeast vertex is fixed
to one of 11 points, and the support is confined to its corresponding
rectangle.

## Pulling the support bound back to the F2 coordinates

The translation can be removed from the necessary support condition. The
pre-generator \(R_{c,d}\) is constructed from \(P_{c,d}\) by the
triangular coefficient recursion for the cube root: its leading monomial is
\(x^5y^{20}\), and \(P_{c,d}-R_{c,d}^3\) vanishes on
\([10,15]\times[40,60]\). At each step the new coefficient of
\(R_{c,d}\) has coefficient 3 in the defining equation; all earlier terms
are already determined. Therefore every coefficient of \(R_{c,d}\), and
of \(\zeta(P_{c,d})=P_{c,d}-R_{c,d}^3\), is a polynomial in \(c,d\)
and the coefficients of \(P\).

For generic \(c,d\), Magnus IV says every coefficient of
\(\zeta(P_{c,d})\) outside \([0,9]\times[0,36]\) is zero.
Each such coefficient is a polynomial in \(c,d\); since it vanishes on a
dense open set, it vanishes identically. Specializing to \(c=d=0\) gives
the following necessary condition directly in the standardized F2
coordinates:

> There is a polynomial \(R\) supported in \([0,5]\times[0,20]\),
> with normalized \(x^5y^{20}\)-coefficient, such that
> \(P=R^3+Z\) and \(\operatorname{Supp}(Z)\subseteq[0,9]\times[0,36]\).

Substituting this into the Keller equation gives the exact identity
\[
3R^2[R,Q]+[Z,Q]=1,
\]
so \(R^2\) must divide \(1-[Z,Q]\). This is a concrete algebraic condition
to impose inside each of the 11 inner-vertex branches; it has not yet been
used to exclude any branch.

This specialization follows from the coefficient recursion; Theorem 5.8 is
not being applied directly at \(c=d=0\), where the translated polynomial
need not have all rectangle corners nonzero. This consequence uses the
proved Magnus-IV inner-vertex theorem and the explicit dense-open argument
above. It does not require computing the
principal polynomial \(W_F\), its divisor \(\delta\), or the
innermost correction terms. It is still only a necessary condition and does
not yet exclude the pair.

For a branch with northeast vertex \((x,y)\), let \(B_{x,y}=[0,x]\times[0,y]\).
Since \(R^3=P-Z\), the Newton polygon identity gives
\[
3N(R)=N(R^3)\subseteq\operatorname{conv}(S_P\cup B_{x,y}).
\]
This first cuts \(R\)'s 126-term rectangle to 76--96 possible lattice points,
depending on the branch. We then use the actual coefficient equations of
\(R^3=P-Z\) at exponents outside \(S_P\cup B_{x,y}\). After normalizing
\(r_{5,20}=1\), exact convolution shows that their pure-power consequences
force every extra \(R\)-coefficient outside \(5i-j\le5\) to vanish in all
11 branches. For example, the \((4,16)\) branch has the forced coefficient
equation \([x^6y^{12}]R^3=r_{2,4}^3=0\). The larger branches require several
rounds of the same exact single-monomial deduction; products of distinct
unknowns are not treated as individually zero.

Consequently \(R\) has at most the same 76 support points in every branch.
Then \(R^3\) is supported in \(S_P\), so \(Z=P-R^3\) is supported in
\(S_P\cap B_{x,y}\), improving the branch rectangles too. The reduced
variable counts (75 free coefficients of \(R\), the indicated \(Z\)-support,
and 1,576 coefficients of \(Q\)) are:

| Northeast vertex of \(Z\) | Initial \(Z\)-rectangle terms | Forbidden \(R^3\) coefficient rows | Forced-zero \(R\) terms | Remaining \(Z\) terms | Variables \((R,Z,Q)\) |
|---:|---:|---:|---:|---:|---:|
| \((1,3)\) | 8 | 0 | 0 | 8 | 1,659 |
| \((1,4)\) | 10 | 0 | 0 | 10 | 1,661 |
| \((2,7)\) | 24 | 0 | 0 | 24 | 1,675 |
| \((2,8)\) | 27 | 0 | 0 | 27 | 1,678 |
| \((3,12)\) | 52 | 0 | 0 | 52 | 1,703 |
| \((4,16)\) | 85 | 14 | 1 | 80 | 1,731 |
| \((5,20)\) | 126 | 45 | 6 | 111 | 1,762 |
| \((6,24)\) | 175 | 57 | 9 | 145 | 1,796 |
| \((7,28)\) | 232 | 66 | 12 | 182 | 1,833 |
| \((8,32)\) | 297 | 72 | 15 | 222 | 1,873 |
| \((9,36)\) | 370 | 75 | 20 | 265 | 1,916 |

The exact rational hull, all forbidden coefficient positions, and the
iterated pure-power deductions are replayed by
[magnus_75_125_cube_support.py](../../src/magnus_75_125_cube_support.py). This is a
branchwise support reduction, not an exclusion: the full \([P,Q]=1\) system
still has to be imposed on these finite models.

## The Magnus cube-root edge is the F2 quartic edge

This gives an exact coordinate bridge between the F2 quartic `h` and the
Magnus pre-generator `\\mathcal R`; it does not identify `h` with the separate
degree-five upper-y carrier discussed below. Use the F2 coordinates
`t=xy^5`, `z=y^{-1}`, and the edge weight `\\ell=5i-j`. The F2 top edge is
`P_15(t)=a t^3 h(t)^3`, with `h` a monic quartic. The Magnus reduction gives
`P=\\mathcal R^3+Z`, support `\\ell(\\mathcal R)\\le5`, and
`Supp(Z)\\subseteq[0,9]\\times[0,36]` in the same standardized coordinates.

Only products of three `\\ell=5` terms of `\\mathcal R` can contribute to
`\\ell=15` in `\\mathcal R^3`. Write that edge as
`\\mathcal R_5(t)=t g(t)`, where `g` has degree four. A point of the `Z`
rectangle on `\\ell=15` has `i=3,...,9`, so `Z_15(t)` has degree at most 9.
Consequently `t^3(a h^3-g^3)=Z_15`, and `a h^3-g^3` has degree at most 6.
The top corner coefficient gives `a=lc(g)^3`. Choose `\\omega^3=a`; then
`\\omega h` and `g` have the same nonzero leading coefficient. If they were
distinct quartics, the factorization
`(\\omega h)^3-g^3=(\\omega h-g)((\\omega h)^2+(\\omega h)g+g^2)`
would have degree at least 8: the second factor has degree 8 and leading
coefficient three times the square of the common leading coefficient, while
the first has degree at least 0. This contradicts the degree-at-most-6
bound. Hence `g=\\omega h`, and after normalizing the top coefficient,
`\\mathcal R_5(t)=t h(t)`.

Thus the two local factor cases impose explicit coefficient values on the
Magnus pre-generator edge, at `(x,y)=(1,0),(2,5),(3,10),(4,15),(5,20)`:

| F2 quartic | Edge coefficients in increasing powers of `t` |
| --- | --- |
| `2+2`, `\\alpha^2-3\\alpha+1=0` | `3\\alpha-1, 2-8\\alpha, 7\\alpha, -2(\\alpha+1), 1` |
| `2+1+1` | `3, -9, 10, -5, 1` |

The exact support/coefficient replay is
[`f2_75_125_magnus_edge_bridge_20261021.py`](../../src/f2_75_125_magnus_edge_bridge_20261021.py).
All five edge positions survive the 76-point support cut, so this bridge
alone excludes neither branch. It makes the branch equations available
inside the same `\\mathcal R,Z,Q` systems; the next elimination must impose
these edge coefficients together with
`3\\mathcal R^2[\\mathcal R,Q]+[Z,Q]=1`, including every zero-pivot stratum.

## Coupling the edge pins into the Magnus Q-tail support model

The edge bridge also pins the remainders. On the `ell=15` face, `P` and
`\\mathcal R^3` agree, so `Z=P-\\mathcal R^3` vanishes at every supported
point of that face. On the `ell=25` face of `Q`, the `b\\mathcal R^5` term is
the full face and `c\\mathcal R^4` has weight at most 20; therefore the tail
`S` in `Q=b\\mathcal R^5+c\\mathcal R^4+S` has zero coefficient at each
supported `ell=25` point. In the `y<=60` tail envelope these are the 13
points `(i,5i-25)`, `i=5,...,17`.

After normalizing the corner coefficient `r_(5,20)=1`, the quartic edge
bridge fixes the other four supported coefficients of `\\mathcal R` on
`ell=5`. The exact support model consequently has 71 free `\\mathcal R`
coefficients (rather than 75) and 695 free `S` coefficients (rather than
708); the number of `Z` face zeros depends on the Magnus vertex. The exact
replay `magnus_75_125_q_tail.py` verifies these counts and the coefficient
identities.

| Magnus vertex | `Z` support terms | `Z` face zeros | Raw variables | Edge-refined variables | After permitted x-shift |
| --- | ---: | ---: | ---: | ---: | ---: |
| `(1,3)` | 8 | 0 | 793 | 776 | 776 |
| `(1,4)` | 10 | 0 | 795 | 778 | 778 |
| `(2,7)` | 24 | 0 | 809 | 792 | 792 |
| `(2,8)` | 27 | 0 | 812 | 795 | 795 |
| `(3,12)` | 52 | 1 | 837 | 819 | 819 |
| `(4,16)` | 80 | 2 | 865 | 846 | 846 |
| `(5,20)` | 111 | 3 | 896 | 876 | 876 |
| `(6,24)` | 145 | 4 | 930 | 909 | 898 |
| `(7,28)` | 182 | 5 | 967 | 945 | 933 |
| `(8,32)` | 222 | 6 | 1007 | 984 | 971 |
| `(9,36)` | 265 | 7 | 1050 | 1026 | 1012 |

These are exact necessary support cuts and reduced variable counts, not
dimensions of solution sets and not exclusions. The full commutator equations
and `\\mathcal R^2 | 1-[Z,S]` still have to be imposed on every surviving
coefficient stratum.

In the existing `(9,36)` depth-eight recurrence, the pinned coefficient
`r_(4,15)` occurs exactly once in the solved coefficient `[x]S_54`, as
`(2/5)c_2 r_(4,15)`. The edge bridge therefore makes this contribution
`-2c_2` on the `2+1+1` quartic branch, and
`-(4/5)(alpha+1)c_2` on the `2+2` branch, where
`alpha^2-3alpha+1=0`. The current retained fixed rows do not involve
`r_(4,15)`, so this only specializes the solved `S_54` slice; it must be
carried into the next remainder row before it can yield an exclusion.

## A high-y Euclidean expansion cuts the Q tail

Write \(R=\sum_{j=0}^{20}r_j(x)y^j\), and set
\(\rho(x)=r_{20}(x)\). The support condition \(5i-j\le5\) and the
normalization \(r_{5,20}=1\) give \(\deg\rho=5\). Since
\(P=R^3+Z\) and \(\deg_y Z\le36\), the top y-slice of \(P\) is
\(p_{60}=\rho^3\). Comparing the coefficient of \(y^{159}\) in
\([P,Q]=1\) gives
\[
  3p_{60}q_{100}'-5p_{60}'q_{100}=0,
\]
so \(q_{100}=b\rho^5\) for a nonzero scalar \(b\).

Over \(K(x)[y]\), divide by the degree-20 polynomial \(R\) repeatedly:
\[
  Q=\sum_{k=0}^{5}A_k(x)R^k+S_0,\qquad \deg_y S_0<20.
\]
The top slice gives \(A_5=b\). Also
\([R^3,Q]=1-[Z,Q]\), whose y-degree is at most 135, and
\([R^3,Q]=3R^2[R,Q]\); hence \(\deg_y[R,Q]\le95\). The identity
\[
  [R,A_k(x)R^k]=-A_k'(x)R^kR_y
\]
shows that a nonzero \(A_4'\) would contribute a noncancellable
y-degree-99 term to \([R,Q]\): terms with \(k\le3\) have degree at most 79,
and \([R,S_0]\) has degree at most 38. Therefore \(A_4'=0\), so
\(A_4=c\in K\). This proves
\[
  Q=bR^5+cR^4+S,\qquad \deg_y S\le60.
\]
The supports of \(R^5\) and \(R^4\) lie inside the original Q envelope.
The remaining tail \(S\) therefore has only the 708 envelope monomials with
\(y\le60\), rather than 1,576 Q coefficients.

There is a further consequence in the four branches whose northeast
vertices are \((6,24),(7,28),(8,32),(9,36)\). Since
\([R^3,Q]=[R^3,S]\) has y-degree at most 119,
\([Z,Q]=1-[R^3,Q]\) also has degree at most 119. The terms
\([Z,R^4]\) and \([Z,S]\) have degrees at most 115 and 95, respectively;
because \(b\ne0\), this forces \(\deg_y[Z,R^5]\le119\), hence
\(\deg_y[Z,R]\le39\). If \(e=\deg_y Z\) is 24, 28, 32, or 36 and
\(z_e(x)=[y^e]Z\ne0\), the coefficient of \(y^{e+99}\) in
\([Z,R^5]\) must vanish:
\[
  5\rho^4(20\rho z_e'-e z_e\rho')=0.
\]
Thus \(z_e^{20}=C\rho^e\) for a nonzero constant \(C\). Since
\(20/\gcd(20,e)=5\) in each of these four cases and \(\deg\rho=5\),
\(\rho\) must be a scalar fifth power of a linear polynomial; specifically,
\(\rho=\lambda(x-\alpha)^5\) and
\(z_e=\mu(x-\alpha)^{e/4}\). This prunes the top-slice multiplicity
pattern, but leaves the corresponding branches alive. The source translation
\(x\mapsto x+\alpha\) preserves all of the support envelopes (it only lowers
the x-exponent of monomials), so in these four branches we may further take
\(\rho=x^5\) and \(z_e=\mu x^{e/4}\). This sets the five lower coefficients
of \(\rho\) and the \(e/4\) lower coefficients in the top y-slice of \(Z\)
to zero.

The same decomposition gives a useful congruence. Since
\[
 [Z,Q]=[Z,S]+(5bR^4+4cR^3)[Z,R]
\]
and \(1-[Z,Q]=3R^2[R,Q]\), it follows that
\[
  R^2\mid 1-[Z,S].
\]
This is a low-y necessary condition on the 708-term tail. It is not yet a
contradiction; the next elimination can use the exact remainder of
\(1-[Z,S]\) on division by \(R^2\), together with the high-y slice
conditions above.

The leading y-slice of this remainder condition already prunes the tail.
Write \(1-[Z,S]=R^2T\). Since \(\deg_y[Z,S]\le e+59\) and
\(\deg_yR^2=40\), the coefficient of \(y^{e+59}\) in \([Z,S]\) must be
divisible by the leading coefficient \(\rho^2=x^{10}\) of \(R^2\). Put
\(S_{60}(x)=[y^{60}]S\). This coefficient is
\[
  60z_e'S_{60}-e z_e S_{60}'
  =4s\mu\sum_{i=0}^{17}(15-i)[x^i]S_{60}\,x^{s+i-1},
  \qquad e=4s,\quad z_e=\mu x^s.
\]
For \(i<11-s\), the displayed x-power is less than 10 and its scalar is
nonzero, so divisibility forces \(x^{11-s}\mid S_{60}\). The four conditions
are \(x^5\mid S_{60}\), \(x^4\mid S_{60}\), \(x^3\mid S_{60}\), and
\(x^2\mid S_{60}\), respectively. This removes another 5, 4, 3, and 2
coefficients from the four tail models.

The next coefficient of the same division also constrains
\(S_{59}(x)=[y^{59}]S\). Writing \(T_{e+19}=[y^{e+19}]T\), the
\(y^{e+58}\) row is
\[
 -[Z,S]_{e+58}=x^{10}T_{e+18}+2x^5r_{19}T_{e+19}.
\]
The second term has x-order at least 8 because \(x^3\mid r_{19}\). The part
of \([Z,S]_{e+58}\) involving \(S_{59}\) is
\(59z_e'S_{59}-e z_e S_{59}'\). Its remaining part involves
\(z_{e-1}\) and \(S_{60}\), and has x-order at least 8 from
\(\operatorname{ord}_0(z_{e-1})\ge s-2\) and
\(\operatorname{ord}_0(S_{60})\ge11-s\). Thus below x-order 8 only the
\(S_{59}\) operator remains; its multiplier on coefficient \([x^i]S_{59}\)
is \(\mu s(59-4i)x^{s+i-1}\), which is nonzero for every supported i.
Consequently \(x^{9-s}\mid S_{59}\) for \(s=6,7,8\); there is no new
condition when \(s=9\). This removes 3, 2, 1, and 0 further coefficients.

The same row also fixes two coefficients immediately above those initial
orders. Write
\[
r_{19}=x^3(a_3+a_4x),\qquad
S_{60}=x^{11-s}(b_0+b_1x+\cdots),\qquad
S_{59}=\sum_i c_i x^i.
\]
The quotient row requires the coefficients of \(x^8\) and \(x^9\) in
\([Z,S]_{e+58}+2x^5r_{19}T_{e+19}\) to vanish. Exact expansion gives
\[
\begin{aligned}
0&=\mu s(23+4s)c_{9-s}
 +\frac{\mu s}{5}(4s^2-25s-269)a_3b_0,\\
0&=\mu s(19+4s)c_{10-s}
 +\frac{\mu s}{5}(4s^2-29s-228)a_3b_1
 +\frac{\mu s}{5}(4s^2-25s-209)a_4b_0.
\end{aligned}
\]
Both coefficients of the \(c_i\) are nonzero in characteristic zero, so these
rows solve for \(c_{9-s}\) and \(c_{10-s}\) in every branch, including
\(s=9\). In the order \(s=6,7,8,9\), the resulting pairs are
\[
\begin{array}{c|cc}
s & c_{9-s} & c_{10-s}\\\hline
6&\frac{55}{47}a_3b_0&\frac65a_3b_1+a_4b_0\\
7&\frac{248}{255}a_3b_0&a_3b_1+\frac45a_4b_0\\
8&\frac{213}{275}a_3b_0&\frac45a_3b_1+\frac35a_4b_0\\
9&\frac{34}{59}a_3b_0&\frac35a_3b_1+\frac25a_4b_0.
\end{array}
\]
The formulas are replayed by
[magnus_75_125_tail_remainder_next_slice.py](../../src/magnus_75_125_tail_remainder_next_slice.py).
They are exact necessary conditions on the Q-tail slices, not exclusions.

## The next upper-y quotient row gives two more exact selectors

Use the same notation, but write
\[
r_{19}=x^3(a_3+a_4x),\qquad
r_{18}=\rho_0+\rho_1x+\rho_2x^2+\cdots,\qquad
S_{60}=x^{11-s}(u_0+u_1x+u_2x^2+\cdots).
\]
For \(s=7,8,9\), the earlier commutator row gives \(\rho_0=0\).  The
coefficient of \(y^{e+57}\) in \(R^2T=1-[Z,S]\), after eliminating the
two preceding quotient slices, has an unknown new quotient term multiplied
by \(x^{10}\).  The new \(S_{58}\) slice first contributes at \(x^{s-1}\).
Consequently, the coefficients below \(x^{s-1}\) are fixed obstructions.
The resulting exact equations are
\[
\begin{array}{c|l}
s & \text{fixed coefficient equations}\\\hline
6,7 & \text{none below the }S_{58}\text{ image}\\
8 & u_0(2529a_3^2-15125\rho_1)=0\\
9 & u_0(59a_3^2-240\rho_1)=0,\\
  & 55a_3^2u_1+114a_3a_4u_0-235\rho_1u_1-190\rho_2u_0=0.
\end{array}
\]
No division by \(u_0\), \(a_3\), or a coefficient of \(r_{18}\) is used, so
the zero branches remain included.  For example, the first equation at
\(s=9\) only fixes \(\rho_1\) when \(u_0\ne0\); if \(u_0=0\), the second
equation remains as written and must be followed through its own zero cases.
Resolving those zero cases gives the following exact leaves for these rows:
\[
\begin{array}{c|l}
s=8 & u_0=0\quad\text{or}\quad
u_0\ne0,\ \rho_1=\frac{2529}{15125}a_3^2,\\[2pt]
s=9 & u_0\ne0:\ \rho_1=\frac{59}{240}a_3^2,\quad
\rho_2=\frac{57}{95}a_3a_4-\frac{133}{9120}\frac{a_3^2u_1}{u_0};\\
& u_0=0,\ u_1=0:\ \text{no further condition from these rows};\\
& u_0=0,\ u_1\ne0:\ \rho_1=\frac{11}{47}a_3^2.
\end{array}
\]
The rational expression on the \(u_0\ne0\) leaf is used only after that
nonvanishing branch is explicit.  The checker verifies these substitutions
symbolically as well as the original product equations.
The direct bivariate bracket slices and \(R^2\) coefficients, followed by
the exact low-degree coefficient extraction, are replayed in
[magnus_75_125_tail_remainder_lowrows.py](../../src/magnus_75_125_tail_remainder_lowrows.py).
These selectors are additional necessary conditions; they do not exclude
any of the four high branches.

On the `s=9`, `e=36` leaf, the saved exact `z_6` commutator calculation
has one remaining row, `(612/5) mu sigma_0^2=0`. Since `mu!=0`, every
solution lies on `sigma_0=0`; the nonzero-`sigma_0` leaf is empty. The
`z_7` residual also forces `r_(13,8)=r_(13,9)=r_(13,10)=0`; its `x^3` row
reduces to `(594/125)mu rho_2^2(2a_4 rho_2-5sigma_1)=0`. Thus either
`rho_2=0` or `5sigma_1=2a_4 rho_2`. The depth-eight `x^1` row gives
`(72/25)mu c_2 rho_2^3=0`; on `rho_2!=0`, this forces `c_2=0`, and the
depth-seven `x^3` row then forces `u_2=0`. The `rho_2=0` leaf and higher
coefficients remain open. These are exact local splits on this Magnus
support vertex, still short of excluding either full F2 factor branch.

On the remaining `rho_2=0` leaf, the replay
`magnus_75_125_depth8_s9_rho2zero_reduction_20261021.py` specializes the
saved depth-seven/eight Laurent rows further. Its first new row is
`-(11493/25)mu*c_2*sigma_1^2=0`, splitting into `c_2=0` or `sigma_1=0`.
The depth-seven `x^4` row gives `sigma_1^2*u_2=0` on `c_2=0`, and
`a_4*c_2*tau_0=0` on `sigma_1=0`. If `c_2=u_2=0`, the depth-eight `x^3`
row also gives `c_3*sigma_1^2=0`. These necessary equations still leave
parameter strata and do not exclude the focused `2+1+1` branch.

## The first three upper-y commutator rows prune the four high branches

In the four branches just considered, write \(e=4s\), so \(s=6,7,8,9\),
and after translation write
\[
 R=x^5y^{20}+r_{19}(x)y^{19}+r_{18}(x)y^{18}+r_{17}(x)y^{17}+\cdots,
 \qquad [y^e]Z=\mu x^s,
\]
where \(\mu\ne0\) and each of \(r_{19},r_{18},r_{17}\) has x-degree at
most 4. The bound \(\deg_y[Z,R]\le39\) forces its coefficients at
\(y^{e+18},y^{e+17},y^{e+16}\) to vanish. The slice operator at depth
\(\ell\) is
\[
 z\longmapsto5x^4\bigl(4xz'-(e-\ell)z\bigr).
\]
For \(\ell=1,2,3\) it has no polynomial kernel, since \(e-\ell\) is not
divisible by 4. The unique slices are therefore the first coefficients of
the formal fractional power \(\mu R^{s/5}\):
\[
\begin{aligned}
[y^{e-1}]Z&=\mu\frac{s}{5}x^{s-5}r_{19},\\
[y^{e-2}]Z&=\mu\left(\frac{s}{5}x^{s-5}r_{18}
 +\binom{s/5}{2}x^{s-10}r_{19}^2\right),\\
[y^{e-3}]Z&=\mu\left(\frac{s}{5}x^{s-5}r_{17}
 +\frac{s(s-5)}{25}x^{s-10}r_{19}r_{18}
 +\binom{s/5}{3}x^{s-15}r_{19}^3\right).
\end{aligned}
\]
Every displayed slice must be a polynomial in \(x\). Let
\(\nu=\operatorname{ord}_{x=0}r_{19}\), with \(\nu=\infty\) if
\(r_{19}=0\). The second slice first gives
\(s-10+2\nu\ge0\). In the third slice, for each excluded smaller
valuation, the cubic term has strictly lower x-order than both the mixed
term and the \(r_{17}\) term, so it cannot cancel. The exact conclusions are:

| Magnus vertex \((x_0,e)\) | Required divisibility | Further exact split |
| --- | --- | --- |
| \((6,24)\) | \(x^3\mid r_{19}\) | \([x^3]r_{19}\,r_{18}(0)=0\) |
| \((7,28)\) | \(x^3\mid r_{19}\), \(x\mid r_{18}\) | \((5b_1-a_3^2)(5b_1-2a_3^2)=0\) |
| \((8,32)\) | \(x^3\mid r_{19}\), \(x\mid r_{18}\) | â€” |
| \((9,36)\) | \(x^3\mid r_{19}\), \(x\mid r_{18}\) | â€” |

Here \(a_3=[x^3]r_{19}\) and \(b_1=[x]r_{18}\) in the \((7,28)\) row. The
depth-four row is available for \(e=28,32,36\); its operator has the
one-dimensional kernel spanned by \(x^{s-1}\), which contributes no negative
x-powers. For \((7,28)\), the coefficient of \(x^{-1}\) in the formal
fractional-power slice is exactly
\(\frac7{625}(5b_1-a_3^2)(5b_1-2a_3^2)\). For \((8,32)\) and \((9,36)\),
the uniquely lowest negative terms force the additional divisibility
\(x\mid r_{18}\). For \((9,36)\), the same row strengthens
\(x^2\mid r_{19}\) to \(x^3\mid r_{19}\).

After the \(S_{60}\) and \(S_{59}\) conditions, the \((6,24)\) product
equation splits into two 907-variable children, with their 906-variable
intersection. The \((7,28)\) row splits into two 944-variable leaves, with
their 943-variable intersection. The other two high branches have 986 and
1030 variables after the divisibility conditions. These are still necessary
local prunings, not exclusions of any branch.

The exact bracket-slice identities, resonant depth-four coefficient, no-kernel
arithmetic, x-adic valuation splits, and variable counts are replayed by
[magnus_75_125_commutator_front.py](../../src/magnus_75_125_commutator_front.py).

The replay checker [magnus_75_125_q_tail.py](../../src/magnus_75_125_q_tail.py)
verifies the bracket identities, the exact support counts, and the
multiplicity-partition consequences. Before using the final source
translation, the variable counts are 793, 795, 809, 812, 837, 865, 896, 930,
967, 1007, and 1050 in the listed vertex order. After normalizing the root in
the last four branches, they become 793, 795, 809, 812, 837, 865, 896, 919,
955, 994, and 1036. These count the branch support variables before imposing
the remaining Jacobian equations. The Euclidean degree argument above is the
written proof of the symbolic reduction; the checker is not a unit-ideal
certificate and no branch is excluded here.

## One further \(R^2\)-remainder row

Write \(\sigma_j=[x^j]r_{17}\), and let \(c_j=[x^j]S_{59}\).  In the
coefficient of \(y^{e+56}\), the new \(S_{57}\) slice first contributes at
\(x^{s-1}\).  I use the preceding quotient rows and solve the coefficients
\(d_j=[x^j]S_{58}\), \(0\le j\le10-s\), needed to make the
\(y^{e+57}\) quotient slice polynomial.  The remaining coefficients below
\(x^{s-1}\) give these exact restrictions:

- **\(s=6\).** On the \(a_3=0\) part of \(a_3\rho_0=0\), this row gives no
  new fixed coefficient.  On the \(\rho_0=0\) part it gives
  \(a_3u_0(a_3^2-5\rho_1)=0\).
- **\(s=7\).** On either branch
  \(\rho_1=a_3^2/5\) or \(\rho_1=2a_3^2/5\), the first fixed equation is
  \(a_3^3u_0=0\).  If \(a_3=0\), the two root branches coincide at
  \(\rho_1=0\), and the next equation is \(\sigma_0u_0=0\).  If
  \(a_3\ne0\), then \(u_0=0\).
- **\(s=8\).** On \(u_0\ne0\), the row forces
  \(a_3=\rho_1=\sigma_0=0\), then
  \(\sigma_1=\frac25a_4\rho_2\).  On \(u_0=0\), it gives
  \[
  a_3u_1(a_3^2-5\rho_1)=0
  \]
  and the additional exact equation
  \[
  \begin{aligned}
  0={}&15869a_3^3u_2+26552a_3^2a_4u_1-14629a_3^2c_3
  -91520a_3\rho_1u_2-24748a_3\rho_2u_1\\
  &-90288a_4\rho_1u_1+87480c_3\rho_1+101520\sigma_0u_1.
  \end{aligned}
  \]
- **\(s=9\).** On \(u_0\ne0\), it forces
  \(a_3=\rho_1=\rho_2=\sigma_0=\sigma_1=0\), followed by
  \(\sigma_2=\frac35a_4\rho_3-\frac4{25}a_4^3\).  On
  \(u_0=0,u_1\ne0\), it forces \(a_3=\rho_1=\sigma_0=0\), followed by
  \[
  -289050a_4\rho_2u_1+263200c_2\rho_2+327825\sigma_1u_1=0.
  \]
  On \(u_0=u_1=0\), the two fixed equations left by this row are
  \[
  \begin{aligned}
  0={}&1828a_3^3u_2-1740a_3^2c_2-7620a_3\rho_1u_2+7075c_2\rho_1,\\
  0={}&1700a_3^3u_3+5244a_3^2a_4u_2-1620a_3^2c_3-3360a_3a_4c_2\\
  &-7325a_3\rho_1u_3-6580a_3\rho_2u_2-7200a_4\rho_1u_2
  +5600c_2\rho_2+6925c_3\rho_1+8300\sigma_0u_2.
  \end{aligned}
  \]

Here \(u_j=[x^j](x^{s-11}S_{60})\).  The direct bivariate bracket and
quotient-row expansion, including the coefficient substitutions and the
stated simplifications, is checked by
[magnus_75_125_tail_remainder_depth4.py](../../src/magnus_75_125_tail_remainder_depth4.py).
These are still local necessary conditions: the vanishing-
\(S_{60}\) leaves remain, and no \((75,125)\) branch has been excluded.

## The following \(R^2\)-remainder row

The next quotient numerator is the coefficient at \(y^{e+55}\). I retained
the fourth \(Z\)-slice, including its one-dimensional resonant kernel for
\(s=7,8,9\), and kept all the preceding fixed rows as branch equations.
The exact checker
[magnus_75_125_tail_remainder_depth5.py](../../src/magnus_75_125_tail_remainder_depth5.py)
replays the Laurent coefficient extraction and the direct bivariate bracket
identity. Its branch-selective mode was also used to audit the terminal
\(s=9, u_0=u_1=0\) leaf. These are further necessary equations, not a
certificate of emptiness.

The new fixed coefficients give the following pruning:

- **\(s=6\), branch \(a_3=0\).** The first new equation is
  \(\rho_0^2u_0=0\). If \(u_0\ne0\), it gives \(\rho_0=0\),
  and the next fixed coefficient then gives \(\rho_1=0\). The
  \(u_0=0\) subbranch is retained.
- **\(s=6\), branch \(\rho_0=0\).** The next selector is
  \[
  u_0(112167a_3^4-616510a_3^2\rho_1+339575\rho_1^2)=0.
  \]
  Together with the preceding \(a_3u_0(a_3^2-5\rho_1)=0\), the
  \(u_0\ne0\) branch forces \(a_3=\rho_1=0\); \(u_0=0\)
  remains a separate leaf.
- **\(s=7\).** On each root branch
  \(\rho_1=a_3^2/5\) and \(\rho_1=2a_3^2/5\), respectively, the
  new selector is a nonzero rational multiple of \(a_3^4u_0=0\). Thus
  each branch first splits into \(a_3=0\) and \(u_0=0\). On the
  \(u_0=0,\ a_3\ne0\) child, the next fixed coefficient forces \(u_1=0\). On
  \(a_3=0\), the two root labels coincide at \(\rho_1=0\), and the next
  coefficient gives \(\sigma_0u_0=0\). The \(x^4\) row is then
  \(\frac{28\mu u_0}{135}(2835a_4\sigma_0+32\rho_2^2)=0\).
  Thus \(u_0=0\) remains a child; if \(u_0\ne0\), then
  \(\sigma_0=0\) and \(\rho_2=0\). The remaining children are \(a_3=0,u_0=0\);
  \(a_3=0,u_0\ne0,\rho_2=0\); and \(a_3\ne0,u_0=u_1=0\).
- **\(s=8, u_0=0\).** The new equation is
  \[
  u_1(203a_3^4-990a_3^2\rho_1+25\rho_1^2)=0.
  \]
  With the earlier \(a_3u_1(a_3^2-5\rho_1)=0\), it forces either
  \(u_1=0\) or \(a_3=\rho_1=0\). On the latter child with \(u_1\ne0\),
  the next quotient row gives \(\sigma_0=0\), then \(\rho_2=0\), and the
  following row requires
  \[
  22180500\tau_0u_1+19139625c_3\sigma_1
  -19747800a_4\sigma_1u_1=0.
  \]
  This narrows that leaf but does not close it.
- **\(s=8, u_0\ne0\).** The new rows force
  \(a_3=\rho_1=\sigma_0=\rho_2=\sigma_1=0\). Deeper slices still
  contain free coefficients; the next fixed row additionally forces
  \(\tau_0=0\).
- **\(s=9, u_0\ne0\).** The new selector is a nonzero rational
  multiple of \(a_3^4u_0=0\), hence \(a_3=0\). The earlier rows then
  give \(\rho_1=\rho_2=\sigma_0=\sigma_1=0\), with later
  coefficients still to be treated.
- **\(s=9, u_0=0, u_1\ne0\).** The coefficient at \(x^3\) is a
  nonzero rational multiple of \(a_3^4u_1\), so \(a_3=0\). The
  previous row gives \(\rho_1=\sigma_0=0\) and leaves
  \(-289050a_4\rho_2u_1+263200c_2\rho_2+327825\sigma_1u_1=0\).
- **\(s=9, u_0=u_1=0\).** Let \(u_2=[x^2](x^{-2}S_{60})\) and
  \(c_2=[x^2]S_{59}\). The depth-4 row and the new depth-5 \(x^4\)
  coefficient form a homogeneous linear system in \(u_2,c_2\). The
  complete zero/nonzero split is
  \[
  (u_2,c_2)=(0,0)\quad\text{or}\quad
  (u_2,c_2)\ne(0,0)\ \text{and}\
  -172a_3^6+2925a_3^4\rho_1-26100a_3^2\rho_1^2
  +70750\rho_1^3=0.
  \]
  This determinant condition is **not** imposed on the zero-vector leaf.
  On the nonzero-vector leaf, \(a_3=0\) implies \(\rho_1=0\); if
  \(a_3\ne0\), then \(t=\rho_1/a_3^2\) satisfies
  \(70750t^3-26100t^2+2925t-172=0\). The cubic has three distinct
  complex roots. If \((u_2,c_2)=(0,0)\), the next depth-4 and depth-5
  rows give a second homogeneous system in \((u_3,c_3)\), hence the
  additional split
  \[
  (u_3,c_3)=(0,0),\quad\text{or}\quad
  (u_3,c_3)\ne(0,0),\ a_3=0,\ \rho_1=0,\quad\text{or}\quad
  (u_3,c_3)\ne(0,0),\ a_3\ne0,\
  69250t^3-24300t^2+2115t-52=0.
  \]
  This second cubic also has three distinct complex roots and has no root
  in common with the first cubic. These are finite algebraic branches, not
  an exclusion.

The exact coefficient-matrix determinant, the retained zero-vector branch,
the successive \((u_3,c_3)\) split, and both cubic discriminants are replayed by
[magnus_75_125_depth5_s9_selector_audit.py](../../src/magnus_75_125_depth5_s9_selector_audit.py).
The \(s=9\) leaves still require the higher fixed coefficients; the other
branches likewise remain open. No degree-pair exclusion follows yet.

## Sources and scope

- Magnus, Applegateâ€“Onishi, and Nagata degree criteria: the degree gcd
  exceptions do not include 25; the individual factor-count condition does
  not include degrees 75 or 125.
- Lee and Li, [*On the two-dimensional Jacobian conjecture: Magnus' formula
  revisited, IV*](https://arxiv.org/abs/2408.01279), especially Lemma 4.12,
  Theorem 4.10, Theorem 5.8, and Corollary 5.10. Conjectures C/D/E are not
  used here.
- The F2 support envelopes and their proof dependencies are recorded in
  [F2-75-125-band-model.md](../F2-75-125-band-model.md).
- The exact rational region enumeration is replayed by
  [magnus_75_125_region.py](../../src/magnus_75_125_region.py).
- The branchwise cube-support deductions are replayed by
  [magnus_75_125_cube_support.py](../../src/magnus_75_125_cube_support.py).
- The high-y Q-tail identities and exact support arithmetic are replayed by
  [magnus_75_125_q_tail.py](../../src/magnus_75_125_q_tail.py).
- The first three upper-y commutator rows are replayed by
  [magnus_75_125_commutator_front.py](../../src/magnus_75_125_commutator_front.py).
- The leading tail-remainder slice and its exact support reduction are
  replayed by
  [magnus_75_125_tail_remainder.py](../../src/magnus_75_125_tail_remainder.py).
- The next two (R^2)-remainder rows, including their explicit zero/nonzero
  branches, are replayed by
  [magnus_75_125_tail_remainder_lowrows.py](../../src/magnus_75_125_tail_remainder_lowrows.py)
  and
  [magnus_75_125_tail_remainder_depth4.py](../../src/magnus_75_125_tail_remainder_depth4.py).

### Universal s=9 ratio selector

The next \[Z,R\] slice was recomputed with only the previously forced
condition \(\rho_0=0\); no coefficients of \(S_{60}\) or \(S_{59}\)
were set to zero. Solving the slice coefficients of \(z_5\) leaves exactly
one residual Laurent coefficient, at \(x^3\):

\[
-\frac{63\mu}{625}a_3
  \bigl(11a_3^4-100a_3^2\rho_1+250\rho_1^2\bigr)=0.
\]

Thus, if \(a_3\ne0\), the ratio \(t=\rho_1/a_3^2\) must satisfy
\(250t^2-100t+11=0\). Its discriminant is \(-1000\), so this is a
nonempty quadratic over the algebraic closure; it cannot be discarded by a
real-root argument.

The earlier nonzero-vector \( (u_2,c_2)\) and \( (u_3,c_3)\) strata
require, respectively,
\[
70750t^3-26100t^2+2925t-172=0,
\qquad
69250t^3-24300t^2+2115t-52=0.
\]
The two exact resultants with the new quadratic are
\(295619000000\) and \(221699000000\), respectively, and both gcds
are 1. Hence \(a_3\ne0\) is impossible on either nonzero-vector stratum.
Setting \(a_3=0\) in either cubic then gives \(\rho_1=0\). The
nonzero-vector strata have therefore been narrowed to the face
\(a_3=\rho_1=0\), but that face is not excluded. The zero-vector leaf
also remains open.

The generic-slice identity and exact resultant checks are replayed by
`magnus_75_125_depth6_s9_universal_ratio_selector.py`; its receipt is
`magnus_75_125_depth6_s9_universal_ratio_selector_20261019.json`. A separate
continuation retains \(u_2,c_2,u_3,c_3\) on the remaining
\(a_3=\rho_1=0\) face.

### Depth-six continuation on the s=9 face

The continuation
`magnus_75_125_tail_remainder_depth6_s9_a3rho1zero.py` specializes only the
already forced face
\(a_3=\rho_1=u_0=u_1=c_0=c_1=0\), leaving
\(u_2,c_2,u_3,c_3\) free. It solves the next \([Z,R]\) slice with no
branch-parameter division and verifies the next coefficient of
\(R^2\mid 1-[Z,S]\). The exact output is preserved in
`magnus_75_125_tail_remainder_depth6_s9_a3rho1zero_20261020.log`.

The new rows imply two further local reductions:

1. On \((u_2,c_2)\ne(0,0)\), the prior fixed row
\(-36\mu(56c_2\rho_2+83\sigma_0u_2)/5=0\), together with
\(-612c_2\mu\sigma_0=0\) and the displayed depth-five \(x^6\) row,
   forces \(\sigma_0=\rho_2=0\). If \(c_2\ne0\), the first two rows
   give this directly. If \(c_2=0,u_2\ne0\), the prior row gives
   \(\sigma_0=0\), after which the depth-five \(x^6\) row gives
   \(\rho_2^2u_2=0\).
2. On \((u_2,c_2)=(0,0)\) and \(c_3\ne0\), the depth-five \(x^6\)
   row first gives \(\sigma_0=0\); the depth-six \(x^5\) row then
   reduces to a nonzero scalar multiple of \(c_3\rho_2^2=0\), so
   \(\rho_2=0\).

These are exact necessary-equation consequences. The exceptional leaf
\((u_2,c_2,c_3)=(0,0,0),u_3\ne0\), the zero-vector leaf, and the remaining
higher coefficients are still open. Even the reductions above do not make
the relevant Magnus strata empty. This continuation therefore does not
exclude the F2 two-factor branch or change either branch-coverage figure.

### Degree-25 leading carrier and its factor strata

For a hypothetical characteristic-zero Jacobian pair of degrees `(75,125)`,
the leading homogeneous forms have the common-power shape
`P_75 = a H^3`, `Q_125 = b H^5`, with `deg(H)=25` and nonzero scalars `a,b`.
Over the algebraic closure, factor `H` into distinct linear factors with
multiplicities. This gives an exhaustive split by the number of distinct
factors: one (a pure 25th power), two, or at least three. The specific choice
`H=x^25+y^25` is squarefree over characteristic zero and has 25 distinct
complex linear factors, so it lies only in the last stratum. In that case
the projective leading form has 25 points on the line at infinity, which is
incompatible with Abhyankar's two-point theorem for a planar Jacobian pair
(*Lectures on Expansion Techniques in Algebraic Geometry*, TIFR 57 (1977),
Theorem 18.13). In this theorem the number of points at infinity is the
number of distinct pairwise-coprime linear factors of the degree form. Since
the leading forms are `H^3` and `H^5`, they have exactly the same distinct
factor sets as `H`; therefore the theorem rules out the entire `r>=3`
leading-`H` stratum, not just this Fermat-type example.

The split is exhaustive. The `r>=3` case is excluded by Abhyankar's
two-point theorem. Shastri's Proposition 2 prints a pointwise equivalence for
`r=1`, but it cites Abhyankar's Theorem 19.4, whose condition is universally
quantified over every Jacobian pair. That global equivalence does not imply
that one particular pair with one point at infinity is automorphic. This
scope issue is decisive: a high-degree source shear `x -> x+y^k` can make the
ordinary degree form of any polynomial a pure power and, applied to both
components, preserves the Keller condition. Treating the printed pointwise
reading as a theorem would therefore settle the still-open Jacobian
conjecture. Accordingly, this application excludes only `r>=3`; both `r=1`
and `r=2` require separate arguments.

For that stratum, write `H=L^a M^b`, with `a,b>0` and `a+b=25`, after a
linear coordinate change. Moskowicz's Proposition 4.1, applied to the
monomial leading-form normal form (swapping the two components so the larger
degree is first), requires `gcd(a,b)>2`. Since `gcd(a,b)` divides 25 and
cannot equal 25 for positive `a,b`, this leaves `(a,b)=(5,20)` or `(10,15)`,
up to swapping. The GGV `F2,j=1` normalization used below has its own
standardized support and top form `H=x^5 y^20`; its polynomial coordinate
change need not preserve the original `H` factor count. Neither the
two-point theorem nor this exponent filter excludes the surviving normalized
two-factor system. A computer-assisted exclusion still needs to cover all
allowed lower coefficients and degenerations in that system.

The exact pair `(H^3,H^5)` has zero Jacobian, so these are statements about
the leading forms of a hypothetical pair, not a candidate pair by themselves.
The F2 quartic labels `2+2` and `2+1+1` describe a later normalized quartic
edge and should not be conflated with the factor count of `H`.

Sources: [Abhyankar, *Lectures on Expansion Techniques in Algebraic Geometry*,
Theorem 18.13](https://mathweb.tifr.res.in/Documents/Publications/Lectures/tifr57.pdf);
[Shastri, *A Discriminant Criterion for the Two-Dimensional Jacobian Problem*,
Proposition 2](https://ir.library.osaka-u.ac.jp/repo/ouka/all/6794/ojm27_01_11.pdf);
[Moskowicz, *A Variation on Magnus' Theorem and its Generalizations*,
Proposition 4.1](https://arxiv.org/abs/1810.08202).

### Moskowicz criteria checked on the standardized F2 monomial chart -- 2026-10-27

The exhaustive F2 support box and exact degree force the standardized ordinary top forms to be monomials `P_75 = a*x^15*y^60` and `Q_125 = b*x^25*y^100`. In Moskowicz's notation (Theorem 2.4), with `p=P`, `q=Q`, this gives `n=60`, `u=deg_x(a_n)=15`, `r=100`, and `v=deg_x(c_r)=25`; hence `A=gcd(n,u)=15` and `C=gcd(r,v)=25`. The reduced data are `(n/A,u/A)=(4,1)` and `(r/C,v/C)=(4,1)`, so Theorem 2.4's second mild condition fails. The equal-data alternative, Theorem 2.6, does not apply either: `gcd(A,C)=5` is neither 1 nor 2, and neither A nor C is 1, 4, or prime. Theorem 2.7 also does not apply because neither A nor C belongs to `{1,4} union primes`. Swapping the two variables or components leaves these gcd obstructions unchanged. Thus these Magnus-type criteria do not close the standardized F2 chart; this is an applicability audit, not an exclusion. Source: Moskowicz, *A Variation on Magnus' Theorem and its Generalizations*, Theorems 2.4, 2.6, 2.7: https://arxiv.org/abs/1810.08202.
