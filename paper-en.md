# Global crystallization of dilute alternating Gaussian chains

OpenAI Codex · Research manuscript · 26 November 2026

**Abstract.** A one-dimensional periodic chain has strictly alternating species, Gaussian attraction between like species and Gaussian repulsion between unlike species. We prove that equal spacing is the unique global minimizer for every even period when the dimensionless squared mean spacing satisfies $ah^2\ge4$. The admissible configurations have no lower bound on their individual gaps beyond strict ordering. At fixed cell length, deleting a closest neutral pair gives a uniform bound on its energy cost. A convexity estimate for the remaining gaps and a tangent inequality for the equidistant energy yield a quantitative lower bound uniform in the period. At fixed dilute density, any sequence whose energy per particle approaches the minimum has mean absolute gap deviation tending to zero.

## 1. Model and results

Let $a>0$, let $N\ge2$ be even, and let $x_{i+N}=x_i+Nh$ with $x_{i+1}>x_i$. The sign of site $i$ is $\sigma_i=(-1)^i$. Distinct particles $i,j$ interact through $-\sigma_i\sigma_j\exp(-a|x_i-x_j|^2)$. Dilation by $\sqrt a$ reduces the calculation to $a=1$ and changes the mean spacing to $\sqrt a\,h$. From here on $h$ and all gaps are in these rescaled units. Write $L=Nh$, $f(r)=e^{-r^2}$, $d_i=x_{i+1}-x_i$ and $S_{i,m}=\sum_{j=0}^{m-1}d_{i+j}$. The doubly counted energy per cell, its periodized-kernel representation and the uniform energy per particle are

$$
H_N(X)=2\sum_{i=1}^N\sum_{m\ge1}(-1)^{m+1}f(S_{i,m})
=N-\sum_{i,j=1}^N\sigma_i\sigma_j K_L(x_i-x_j),
$$

$$
K_L(z)=\sum_{k\in\mathbb Z} f(z+kL),\qquad
e(h)=2\sum_{m\ge1}(-1)^{m+1}f(mh).
$$

The subtracted diagonal in the interaction is exactly the zero-distance self term. Images of the same site in other periods are included. All series are absolutely convergent for each strictly ordered periodic cell.

**Theorem 1.** If $h\ge2$, every such configuration satisfies $H_N(X)\ge N e(h)$. Equality holds precisely when $d_i=h$ for every $i$, hence the minimizer is unique up to translation and indexing.

**Theorem 2.** Fix $h\ge2$. For any sequence of even periods $N\to\infty$ at this mean spacing,
$$
H_N(X)/N-e(h)\longrightarrow0
\quad\Longrightarrow\quad
\frac1N\sum_{i=1}^N |d_i-h|\longrightarrow0 .
$$
In particular the fraction of gaps outside any fixed neighborhood of $h$ vanishes. The statement also controls the total length carried by exceptional gaps.



## 2. Convexity on separated cells

Set $D=\sqrt{\log 3}$, so $f(D)=1/3$, and set $B=D^2$. For every even cell of $M\ge2$ gaps in $[D,\infty)$ and any gap variation $v\in\mathbb R^M$ we have
$$
\frac12 H_M''[v,v]\ \ge\ (1-R)\sum_j f''(d_j)v_j^2,\qquad
R=\sum_{\substack{m\ge2\\m\ {\rm even}}}m^2\frac{f''(mD)}{f''(D)}<\frac{49}{50}. \tag{1}
$$
Here the period can vary with the gaps. Differentiation of the gap series is justified locally by Gaussian decay, including after multiplication by any fixed power of $m$. For $r>1/\sqrt2$, $f''(r)=(4r^2-2)e^{-r^2}>0$, and its logarithm is strictly concave by differentiating $\log(4r^2-2)-r^2$. Also $f''$ decreases beyond $\sqrt{3/2}$. Keep the positive $m=1$ term and discard all other odd terms of the Hessian. If a gap $d_j$ occurs in a window of length $m\ge2$, then
$$
f''(S_{i,m})\le f''(d_j+(m-1)D)
\le f''(d_j)\,\frac{f''(mD)}{f''(D)} .
$$
The first step uses $d_j+(m-1)D\ge2D>\sqrt{3/2}$. For the second, the ratio $f''(t+c)/f''(t)$ decreases with $t\ge D$ for $c>0$. Bound the square of the sum of $m$ variations by $m$ times their sum of squares. After summing over all starting indices, each gap occurs $m$ times, including multiplicities for windows longer than the period. This proves the first inequality in (1).

The summand of $R$ is $m^2(2Bm^2-1)3^{1-m^2}/(2B-1)$. The series for $2\operatorname{arctanh}(1/2)$ gives $B>1+1/12=13/12$. The factor $(2Bm^2-1)/(2B-1)$ decreases with $B$ for $m\ge2$, so the $m=2$ term is below $184/189$. For $m\ge4$ the summand is below $2m^4\,3^{1-m^2}$. Successive terms of this majorant for all integer $m\ge4$ have ratio at most $(5/4)^4 3^{-9}<1/2$. Thus $R<184/189+1024/3^{15}<49/50$.

Write $F(t;s)=f(t)-f(s)-f'(s)(t-s)$, the error above the tangent at $s$. When $t,s\ge D$, $F(t;s)\ge0$ with equality exactly at $t=s$. By (1), $H_M-2(1-R)\sum_j f(d_j)$ is convex and cyclically invariant. Average all cyclic shifts of a gap vector with mean $s$ and apply Jensen's inequality:
$$
H_M(Y)\ge M e(s)+2(1-R)\sum_j F(d_j;s). \tag{2}
$$
Taking uniform gap vectors and variations $v_j=1$ in (1) gives $e''(s)\ge2(1-R)f''(s)>0$ for every $s\ge D$.

## 3. Removing a closest neutral pair

**Lemma.** In a cell of length $L\ge4$, choose a shortest gap $d\le D$ and remove its two endpoints, including all their periodic images. Denote the remaining alternating cell at the same $L$ by $Y$, with energy zero if empty. Then $H(X)-H(Y)\ge2/3$.

Translate the endpoints to $0,d$ with respective signs $+,-$, reversing all signs if necessary. The part of the energy difference from the pair and its own periodic copies is
$$
2-2[K_L(0)-K_L(d)]
=2 f(d)+2\sum_{k\ge1}\{ f(kL-d)+f(kL+d)-2f(kL)\}\ge2 f(d).
$$
Indeed $d\le L/2$, hence $kL-d\ge2$, in the convex range of $f$.

For the interaction with the remaining sites of the infinite periodic extension, measure the distances $s_1<s_2<\cdots$ to the nearer endpoint, separately on each side of $[0,d]$. The first remaining site on either side has opposite sign to its nearer endpoint, and the remaining signs alternate. Deleting copies of an adjacent pair preserves this parity. Thus the cross energy is $2B_\ell+2B_r$, where each side gives a sum
$$
B_*=\sum_{j\ge1}(-1)^{j+1}g_d(s_j),\qquad
g_d(s)=f(s)-f(s+d),\qquad s_1\ge d .
$$
For an empty residual cell both sums vanish. If $1/2\le d\le D$, $g_d$ decreases on $[d,\infty)$ because $(1+d/s)e^{-2sd-d^2}\le2e^{-3d^2}\le2e^{-3/4}<1$. Grouping consecutive terms gives $B_*\ge0$. In this case the cost per deleted particle is at least $f(d)\ge1/3$.

For $0<d<1/2$, pair successive terms in the alternating sum. The intervals $[s_{2j-1},s_{2j}]$ are disjoint, so
$$
B_*\ge-\int_d^\infty (g_d'(s))_+\,ds
\ge-d\int_d^\infty (-f''(x))_+\,dx
=-d\big(\sqrt{2/e}-2d e^{-d^2}\big).
$$
The second inequality uses $g'_d(s)=\int_s^{s+d}-f''(x)\,dx$ and the bound $d$ on the measure of contributing $s$ for fixed $x$. The last equality uses the inflection point $1/\sqrt2$. Including both sides, the cost per deleted particle is at least
$$
(1+4d^2)e^{-d^2}-2\sqrt{2/e}\,d
\ge 1+2d^2-2d\ge \tfrac12.
$$
For the first bound use $e^{-d^2}\ge1-d^2$, $4d^2\le1$ and $\sqrt{2/e}<1$. This proves the lemma in both ranges, including a two-site cell.

## 4. Return to the prescribed density

Start with $h\ge2$, so $L=Nh\ge4$. While the shortest gap is below $D$, apply the lemma at fixed $L$. Let $M$ be the number of remaining sites. Then $H_N(X)\ge H_M(Y)+(N-M)/3$. For $M>0$, the survivors have gaps $d'_j\ge D$ and mean $s=L/M\ge h$. Define $\mu(h)=e(h)-h e'(h)$. Since $e''>0$ on $[D,\infty)$, $\mu'(h)=-h e''(h)<0$ there. Moreover
$$
\mu(2)=2\sum_{m\ge1}(-1)^{m+1}(1+8m^2)e^{-4m^2}<18e^{-4}<1/3 .
$$
The alternating terms decrease in magnitude, and the Taylor sum for $e^4$ through degree nine exceeds 54. The tangent bound $e(s)\ge e(h)+e'(h)(s-h)$ yields $M e(s)\ge N e(h)-(N-M)\mu(h)$. Combining with (2),
$$
H_N(X)-N e(h)\ge (N-M)[1/3-\mu(h)]
+2(1-R)\sum_{j=1}^M F(d'_j;s). \tag{3}
$$
For $M=0$ the same bound holds with the last term omitted: $H_N\ge N/3$ and $e'(h)<0$, since $m^2e^{-m^2h^2}$ strictly decreases with $m$ when $h\ge2$. The first coefficient in (3) is bounded below by $\gamma=1/3-18e^{-4}>0$, independently of period and mean spacing in this range. Any deletion therefore gives a strict energy excess. If there is no deletion, equality in (2) requires every gap to equal $h$. This proves Theorem 1.

## 5. Rigidity of near minimizers

Put $\delta_N=H_N/N-e(h)$ and $k=N-M$. If $\delta_N\to0$ at fixed $h\ge2$, (3) gives $k/N\le\delta_N/\gamma\to0$. Thus $M>0$ eventually, $s=Nh/M\to h$ and $\sum_j F(d'_j;s)=o(N)$. For any fixed $\varepsilon>0$, the quotient $F(t;s)/|t-s|$ has a positive lower bound for $t\ge D$, $|t-s|\ge\varepsilon$, uniformly for $s$ in a small compact neighborhood of $h$. Indeed strict convexity gives positivity on compact sets, and as $t\to\infty$ the quotient tends uniformly to $-f'(s)>0$. Splitting at $\varepsilon$ and then letting $\varepsilon\downarrow0$ gives $\sum_j|d'_j-s|=o(N)$, hence $\sum_j|d'_j-h|=o(N)$.

The surviving gaps partition the original cyclic list of gaps into consecutive blocks. At most $k$ blocks contain deleted sites, and together they contain at most $2k$ original gaps: a block with $r$ deleted sites contains $r+1$ original gaps. Their total length is at most $kh+\sum_j|d'_j-h|=o(N)$. The contribution to the original absolute deviation on these blocks is bounded by this length plus $2kh$. All remaining gaps are unchanged. Theorem 2 follows.

## References and scope

Bétermin, Knüpfer and Nolte proposed all-scale crystallization for alternating species with a kernel completely monotone in squared distance [1, Conjecture 2.7]. A four-site example in [2] rules it out for a Gaussian. The follow-up [3] gives lower-energy finite periodic competitors for every $0<b<b_d$, where $b=ah^2$ in original units and $2.47416700<b_d<2.47416702$. Even the strict local minimum at all fixed periods for $b_s\le b<b_d$, $b_s\simeq1.081826$, is beaten by a finite rearrangement there. Theorem 1 proves eventual global crystallization and places its onset, defined as the infimum of $B$ such that all $b\ge B$ have the uniform minimizer, in $[b_d,4]$. We do not determine the minimizers throughout the intervening interval. The short-cell infimum obtained in [4] concerns a different regime, specified by the full period length.

1. L. Bétermin, H. Knüpfer and F. Nolte, *Note on Crystallization for Alternating Particle Chains*, Journal of Statistical Physics 181 (2020), [doi:10.1007/s10955-020-02603-2](https://doi.org/10.1007/s10955-020-02603-2).
2. Codex and G. Jiang, *A Gaussian counterexample and a Riesz bound for alternating particle chains*, [upstream manuscript](https://github.com/xiujiang4/alternating-particle-chains/tree/a48fe4c2623a29dd45edd227410a07baba400250).
3. OpenAI Codex, *An exact crystallization boundary and a long-wave instability in alternating chains*, [companion research draft](https://github.com/gamehhh7/alternating-chain-stability/tree/3a7012e285a7482a9caadc5130cc21bbe8057172).
4. *The energy infimum and the exact cost of a minimum spacing for alternating Gaussian chains*, [upstream contribution 1](https://github.com/xiujiang4/alternating-particle-chains/pull/1).

This AI-generated research manuscript has not undergone independent peer review. Its claims are the stated mathematical results. Literature priority and the sharp global transition remain open to independent checking.
