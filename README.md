# Maxwell's conjecture for three positive point charges

**Claim.** Three positive point charges in ℝ³ with the Coulomb potential (1/r) have at most
**4** isolated equilibrium points. This is Maxwell's (n−1)² bound for n = 3.

**Status of the argument.** Rigorous mathematics, plus a rigorous computer-assisted
verification of four explicit strict inequalities on a compact 4-parameter domain:
interval arithmetic with outward rounding and mean-value forms obtained by automatic
differentiation. Independent review is welcome. See *Caveats* below.

## Background

Maxwell (*Treatise*, 1873) conjectured that n point charges create at most (n−1)²
equilibrium points. Known results:

- For general n the bound is false: Arathoon, Ball and Kvalheim (2026) found 5 charges with 24.
- For n = 3:
  - Gabrielov, Novikov and Shapiro (2007) proved ≤ 12;
  - Gabrielov, Novikov, Novikov and Shapiro (arXiv:2607.28785, 2026) proved ≤ 6;
  - Tsai (2015) proved Maxwell's bound when the three charges have equal magnitudes.
- For n = 3, the bound **4 vs 6** for unequal positive charges was open.

## The idea in one paragraph

All equilibria lie inside the triangle. Their number is 2 + 2·(number of local minima of V
in the plane), so the bound is equivalent to "**at most one local minimum**". By Lami's
theorem, a point x is an equilibrium only for the charges q ∝ (λᵢ rᵢ^{s+2}), where λ are the
barycentric coordinates of x. The Hessian there is ∝ (s+2)T − I, with T = Σλᵢuᵢuᵢᵀ, so x is
a minimum iff |Z(x)| < s/(s+2), where Z = Σ λᵢ uᵢ².

Geometrically, V_q = ⟨q, Φ⟩ is a height function on the surface S = Φ(triangle), with
Φ = (r₁^{−s}, r₂^{−s}, r₃^{−s}). The charge map F is the Gauss map of S, the set
M = {|Z| < κ} is its elliptic region, and the number of minima equals the winding number of
the **caustic** C = F(∂M) around q.

The caustic is always a deltoid-like front. If it has no inflections, exactly three
ordinary cusps, and total tangent turning π, then a Sturm-type argument shows that at most
3 of its tangent lines pass through any point (*hedgehog lemma*). Hence the winding number
is ≤ 1, i.e. at most one minimum. These hypotheses reduce to four strict inequalities in an
explicit local chart, which are certified by interval arithmetic.

Details: [`notes/theory.md`](notes/theory.md).

## Repository map

| file | purpose |
|---|---|
| `notes/theory.md` | the reduction, the lemmas and the hedgehog lemma |
| `core.py`, `fast.py` | direct numerics: critical points, the region M, tracing the caustic |
| `localchart.py`, `closedform.py` | the local chart (θ, λ) and the closed-form invariants |
| `interval.py`, `ad.py` | outward-rounded interval arithmetic and interval forward-mode AD |
| `verify.py`, `certify.py` | branch-and-bound certification of conditions (J), (I), (T), (W) |
| `tests/test_enclosures.py` | enclosures must contain independent float evaluations |
| `scan_*.py`, `stress.py`, `dual_curve.py`, ... | numerical exploration (see below) |

## Numerical findings beyond the proof

- For the Riesz potentials r^{−s}, the hedgehog structure holds for every triangle when
  s ≲ 1.5. It fails above a shape-dependent threshold s*, which is 2.13 for the equilateral
  triangle and about 1.53 near degenerate shapes: the caustic acquires inflection pairs.
  Even then, all winding numbers stayed in {0, 1}, so we found no counterexample to the
  bound 4 for any s.
- Every point of M is the *global* minimizer over M for its own charges: the elliptic part
  of S is a convex cap.

## Reproduce

```
pip install numpy scipy matplotlib
python tests/test_enclosures.py
python certify.py 1.0 16 60     # s = 1 (Coulomb); writes cert_s1.0.json
```
