# Maxwell's conjecture for three positive point charges

> **Theorem (computer-assisted).** Three positive point charges in ℝ³ with the Coulomb
> potential have at most **4** isolated equilibrium points. Generically there are exactly
> 2 or 4.

This is Maxwell's (n−1)² bound for n = 3, for positive charges. Gabrielov, Novikov,
Novikov and Shapiro (arXiv:2607.28785, July 2026) brought the best known bound down from
12 to **6** and noted that closing the gap "would be equivalent to Maxwell's conjectured
bound for three charges". This repository closes the gap for the Coulomb potential.

**Status.** The argument combines rigorous mathematics ([`notes/theory.md`](notes/theory.md))
with a rigorous interval-arithmetic certificate. The certificate checks four explicit
strict inequalities on a compact 4-parameter domain:
- `cert_s1.0.json`: 9.97·10⁶ boxes, **0 unresolved**;
- `cert_s1.0_n10.json`: an independent rerun from a different subdivision, also **0 unresolved**.

It has **not** been peer-reviewed. Please read the caveats below.

## Background

- Maxwell (*Treatise*, 1873) conjectured at most (n−1)² equilibria for n point charges.
- For general n this is false: Arathoon, Ball and Kvalheim (2026) found 5 charges with ≥ 24.
- For n = 3:
  - Gabrielov, Novikov and Shapiro (2007) proved ≤ 12 for positive charges;
  - the 2026 preprint above proved ≤ 6;
  - Tsai (2015) proved ≤ 4 for charges of equal magnitude.

## The argument in brief

1. All equilibria lie inside the triangle of charges, and their number is
   2 + 2·(number of local minima of V in the plane). So the claim is
   **at most one local minimum**.
2. **Lami's theorem.** x is an equilibrium only for the charges q ∝ (λᵢ rᵢ^{s+2}), with λ
   the barycentric coordinates of x. The Hessian there is ∝ (s+2)T − I, with T = Σλᵢuᵢuᵢᵀ.
   So x is a minimum iff |Z(x)| < κ = s/(s+2) (κ = 1/3 for Coulomb), where
   **Z(x) = Σ λᵢ uᵢ²** in complex notation.
3. **Gauss-map picture.** V_q = ⟨q, Φ⟩ is a height function on the surface S = Φ(triangle),
   with Φ = (r₁^{−s}, r₂^{−s}, r₃^{−s}). So x ↦ [q] is the Gauss map of S, and M = {|Z| < κ}
   is its elliptic region. The number of minima equals the winding number of the
   **caustic** C = F(∂M) around q.
4. **Hedgehog lemma.** Suppose a closed front has no inflections, total tangent turning π,
   and exactly 3 ordinary cusps. Then it is simple, and its winding numbers are 0 or ±1.
   The proof uses a Sturm comparison for g + g″ = −ρ: at most 3 tangent lines pass through
   any point.
5. These hypotheses reduce to four strict inequalities in a local chart (θᵢ, λᵢ) of
   (triangle, point), where everything is in closed form:
   - (J) Z is orientation-reversing on the closure of M;
   - (I) no inflections of C;
   - (T) the cusps are simple;
   - (W) the total turning of C is π.

   They are certified by branch and bound with outward-rounded interval arithmetic and
   mean-value forms (interval forward-mode AD).

Along the way:
- Z is (numerically) a diffeomorphism from the triangle onto the unit disk;
- C\* (the dual of the caustic) is automatically one-sided;
- the flat-triangle limit is universal.

## Caveats

- **Floating-point rigor.** The interval arithmetic uses IEEE double with one-ulp outward
  rounding per operation. sin and cos are enclosed rigorously by Taylor series with explicit
  remainder bounds, so the certificate relies only on correctly rounded +, −, ×, ÷ and √.
  It does not depend on the accuracy of the platform libm. A port to MPFR or Arb would
  still be a welcome independent check.
- **Closed forms.** They were derived by hand and validated in three ways:
  - against an independent implementation (agreement to 10⁻¹⁴);
  - against finite differences along traced boundaries of concrete triangles;
  - by containment tests of the enclosures (`tests/test_enclosures.py`).
- **Scope.** Only s = 1 (Coulomb in ℝ³) is certified. For Riesz exponents s ≲ 1.5 the same
  structure holds numerically, so the same certification should go through.
- **Mixed signs are not covered.** One lemma is proved: in-plane maxima are impossible.
  Numerically the bound 4 holds there too. A certificate needs a compactified chart, which
  is not done yet (see §8 of the notes).

## Numerical findings

- **Riesz exponents.** For r^{−s}, the hedgehog structure holds for all triangles when
  s ≲ 1.5. Above a shape-dependent threshold (2.13 for equilateral, ≈ 1.53 near degenerate
  shapes) the caustic acquires inflection pairs. We still never found a winding number of
  2, i.e. no counterexample to the bound 4 for any s we tried.
- **Convex cap.** Every point of M is the global minimizer over M for its own charges.

## Repository map

| file | purpose |
|---|---|
| `notes/theory.md` | reduction, lemmas, hedgehog lemma, certificate design, mixed signs |
| `interval.py`, `ad.py` | outward-rounded interval arithmetic; interval forward-mode AD |
| `verify.py` | interval evaluation of ζ, the Jacobian, Ǐ, the turning rate, the cusp test and (W) |
| `certify.py` | parallel branch and bound; writes `cert_*.json` |
| `closedform.py`, `localchart.py` | float implementations of the closed forms (two independent derivations) |
| `core.py`, `fast.py` | direct numerics: critical points, Z, tracing ∂M |
| `tests/test_enclosures.py` | containment tests of the enclosures |
| `scan_*.py`, `stress.py`, `dual_curve.py`, `caustic_geometry.py`, `test_convexcap.py` | exploration |
| `mixed_caustic.py`, `collinear_mixed.py`, `mixed_explore.py` | mixed-sign exploration |
| `figs/` | figures (minimum regions and caustics) |

## Reproduce

```
pip install numpy scipy matplotlib
python tests/test_enclosures.py
python certify.py 1.0 16 60        # Coulomb, 16 workers  -> cert_s1.0.json   (~7 min)
python certify.py 1.0 16 80 10     # independent subdivision -> cert_s1.0_n10.json
```
