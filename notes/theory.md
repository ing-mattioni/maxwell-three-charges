# Maxwell's problem for three point charges — the reduction

Notation. Charges q₁,q₂,q₃ > 0 at the vertices p₁,p₂,p₃ of a non-degenerate triangle Δ ⊂ ℝ² ⊂ ℝ³.
Potential V(x) = Σ qᵢ |x − pᵢ|^{−s}, with s = 1 the Coulomb case
(s = 2α in the notation of Gabrielov–Novikov–Shapiro). For x in the plane:
rᵢ = |x − pᵢ|, uᵢ = (x − pᵢ)/rᵢ, λᵢ(x) = barycentric coordinates, θᵢ = ∠pⱼ x pₖ
(the angle subtended at x by the side opposite pᵢ), σᵢ = sin θᵢ.

## 1. Where equilibria can be

**Lemma 1.** All equilibria lie in the open triangle Δ° (in the plane of the charges).
For any signs of the charges, none lie off the plane.

*Proof.* Off the plane: if E(x) = 0 at height z ≠ 0, put wᵢ = qᵢ rᵢ^{−s−2}. Then the vertical
component gives Σwᵢ = 0, and the horizontal one gives Σwᵢ(x∥ − pᵢ) = 0, hence Σ wᵢ pᵢ = 0.
Since the pᵢ are affinely independent, w = 0, which is impossible. In the plane but outside
Δ° (positive charges): take a line separating x from Δ. Every term of the field has a
positive component along its normal. On an open edge, the third charge pushes off the edge. ∎

(Collinear charges are trivial: by convexity on each segment, exactly 2 equilibria on the
axis and none off it.)

## 2. Charges as a function of the equilibrium point

**Lemma 2 (Lami).** x ∈ Δ° is a critical point of V_q iff
  q ∝ F(x) := (λᵢ(x) rᵢ(x)^{s+2})ᵢ ∝ (σᵢ rᵢ^{s+1})ᵢ.

*Proof.* ∇V = −s Σ qᵢ rᵢ^{−s−2}(x − pᵢ), and Σ wᵢ(x − pᵢ) = 0 with w > 0 iff w ∝ λ(x).
The second form is Lami's theorem: λᵢ rᵢ ∝ σᵢ. ∎

**Lemma 3 (Hessian).** At a critical point,
  Hess V = s W [ (s+2) T(x) − I ],  W = Σ qᵢ rᵢ^{−s−2} > 0,  T(x) = Σ λᵢ uᵢ uᵢᵀ.

trace T = 1, and the eigenvalues of T are (1 ± |Z|)/2, where (complex notation)
  **Z(x) = Σ λᵢ(x) uᵢ(x)².**

Hence the critical point is
- a non-degenerate local minimum iff |Z(x)| < κ := s/(s+2) (κ = 1/3 for Coulomb),
- a saddle iff |Z| > κ,
- degenerate iff |Z| = κ.

It is never a maximum.

So M := {x ∈ Δ° : |Z(x)| < κ} is a purely geometric "minimum region", and
**#{minima of V_q} = #(F⁻¹([q]) ∩ M).**

**Lemma 4 (Morse count).** For Morse V_q: #saddles = #minima + 2. Hence #equilibria = 2 + 2·#minima.
(Use Poincaré–Hopf on S² = ℝ² ∪ {∞}: each charge and ∞ has index +1.)

So Maxwell's bound of 4 for three positive charges ⇔ **at most one local minimum**, i.e.
F|_M is injective.

## 3. Gauss-map picture

Put Φ(x) = (r₁^{−s}, r₂^{−s}, r₃^{−s}) ∈ ℝ³. Then V_q = ⟨q, Φ⟩ is a height function on the
surface S = Φ(Δ°). Its critical points are the points where q is normal to S, so [F] is the
Gauss map of S. M is the elliptic part of S, ∂M its parabolic curve, and
C := F(∂M) (the *caustic*, i.e. the bifurcation set in charge space) is the Gauss image of the
parabolic curve.

**Lemma 5.** F|_M is a local diffeomorphism with constant orientation.
For q ∉ C, #(F⁻¹(q) ∩ M) = |wind(C, q)|.

*Proof.* Implicit differentiation of Σ qᵢ ∇φᵢ(x) = 0 gives DF⁻¹ = −H⁻¹B. Here
B(dq) = Σ dqᵢ ∇φᵢ has det B ≠ 0 on Δ°, and H is positive definite on M. The count then
follows from degree theory on the closure of M, which is a closed disk by Lemma 7. ∎

## 4. Tangent lines of the caustic and the dual curve

For x₀ ∈ ∂M, let k = k(x₀) be the kernel of H. This is the eigenvector of T for the smaller
eigenvalue; if Z = κe^{iω} then k = i·e^{iω/2}.

**Lemma 6.** The tangent line of C at F(x₀) is { q : ∂_k V_q(x₀) = 0 } = { q : ⟨q, m(x₀)⟩ = 0 },
where mᵢ = rᵢ^{−s−1}(uᵢ·k). So C is a front, and its dual curve is C* = {[m(x₀)]}.
Since k(ω+2π) = −k(ω), we get m(ω+2π) = −m(ω): **C\* is automatically one-sided.**

Consequently, for q ∈ ℙ², n(q) := #{tangent lines of C through q} = #{x₀ ∈ ∂M : ∇V_q(x₀) ⊥ k(x₀)}.

## 5. The hedgehog lemma

**Lemma (hedgehog).** Let C be a closed front in an affine plane. Assume:
- C has no inflection points;
- its tangent line turns monotonically by total angle π;
- it has exactly 3 cusps, each ordinary.

Then wind(C, q) ∈ {0, ±1} for every q ∉ C, and C is a simple closed curve.

*Proof.* Parametrize C by the tangent direction: x(θ) = h n + h′ n′, with h(θ+π) = −h(θ).
Cusps ↔ sign changes of ρ = h + h″.

For q ∉ C, let g(θ) = ⟨q, n(θ)⟩ − h(θ). Then g + g″ = −ρ, and the zeros of g are the
tangent lines through q.

(i) Phase-plane count. Along C, the vector x(θ) − q in the frame (n, n′) is (−g, −g′).
It crosses the n′-axis only clockwise, once per zero of g. This gives
wind = (1 − n)/2, where n = #zeros of g in [0, π).

(ii) Sturm bound. Put W(θ) = g′ sin(θ−c) − g cos(θ−c). Then W′ = (g + g″) sin(θ − c).
On (c, c+π), Rolle applied to g/sin(θ−c) gives #zeros(g) ≤ #sign changes(ρ) + 2.
Choose c at a cusp: then ρ has only 2 sign changes in (c, c+π), so n ≤ 4. Since n is odd,
n ≤ 3, and hence |wind| ≤ 1.

(iii) Simplicity. A double point P = x(θ₁) = x(θ₂) makes θ₁ and θ₂ double zeros of g_P.
That is total multiplicity ≥ 4, and ≥ 5 by parity, contradicting the multiplicity-counting
version of (ii). ∎

With C simple, a non-Morse q (q ∈ C) has exactly one degenerate critical point: a fold or a
cusp of F. Perturbing q to the side where it splits shows that the number of *isolated*
equilibria is ≤ 4 for every q. Non-isolated equilibria would form a curve in ∂M mapped to a
single point, which injectivity of F|∂M excludes.

## 6. The local chart and the conditions to verify

Describe a pair (triangle, point x₀) by the angles θᵢ at x₀ and by λ. Lami gives rᵢ = σᵢ/λᵢ.
Then
  Z = Σ λᵢ vᵢ²   (vᵢ = direction from x₀ to pᵢ),
  ∂λᵢ/∂z = (i/2) ρᵢ (ρᵢ₊₂ v̄ᵢ₊₁ − ρᵢ₊₁ v̄ᵢ₊₂),  ρᵢ = 1/rᵢ = λᵢ/σᵢ,
and, with P = v₁v₂v₃, φⱼ = θⱼ₊₁ − θⱼ₊₂, bⱼ = σⱼ(ρⱼ₊₁ρⱼ₊₂ + ρⱼ² e^{iφⱼ}):
  Z_z̄ = P Σ bⱼ,   Z_z = −P Σ bⱼ v̄ⱼ².

Along ∂M, parametrized by ω = arg Z:
- t = dx₀/dω = (conj(Z_z)·iZ − Z_z̄·conj(iZ)) / (|Z_z|² − |Z_z̄|²);
- k′ = (i/2)k;
- in the rescaled frame q̃ = q / q(x₀), m̃ᵢ = σᵢ(uᵢ·k) and
  m̃′ᵢ = σᵢ(uᵢ·k′) − λᵢ[(s+2)(uᵢ·t)(uᵢ·k) − t·k].

Inflections of C ⇔ m̃′ ∥ m̃. After dividing by σ₃ (WLOG θ₃ = max θᵢ):

  Ǐ = ½σ₁σ₂ + (s+2)Q(v₁,v₃)[λ₁σ₃(v₁·ť) − σ₁λ₃(v₃·ť)] − λ₁σ₃Q(v₃,ť) + σ₁λ₃Q(v₁,ť),

where ť = t/σ₃, Q(a,b) = (a·k)(b·k) = ½Re(a b̄) + ½Re(a b k̄²), and k̄² = −ζ̄/|ζ|.

**Conditions** (domain: θ₃ = max θᵢ ≤ π, λ in the simplex, ζ = Σλᵢvᵢ²):
- **(J)** |Z_z|² − |Z_z̄|² < 0 wherever |ζ| ≤ κ.
- **(I)** Ǐ > 0 wherever |ζ| = κ.
- **(T)** d arg t/dω < ½ at every cusp point (t ∥ k) with |ζ| = κ.
- **(W)** |d arg u₁/dω| < ½ wherever |ζ| = κ.

**Lemma 7.** (J) ⇒ the closure of M is a closed disk, and Z maps it diffeomorphically onto
{|ζ| ≤ κ}. Proof: Z extends to ∂Δ with degree −1 (each edge maps to a point, the angle at
vertex pᵢ maps to an arc of length 2αᵢ), and every preimage of a point with |ζ| ≤ κ has
local degree −1.

Then:
- (T) ⇒ the relative angle arg t − arg k decreases strictly through every cusp. It changes
  by −3π over one turn, since ∂M is a Jordan curve traversed clockwise. So there are
  exactly 3 cusps, all ordinary.
- (I) ⇒ there are no inflections.
- (W) ⇒ u₁·k has exactly one zero, so n(e₁) = 1, and since wind(C, e₁) = 0 the total
  turning is Θ = π.

All four hypotheses of the hedgehog lemma then hold for every triangle. This gives
**at most 4 equilibria** for three positive Coulomb charges.

## 7. The computer-assisted part

`certify.py` covers the whole closed domain
  θ₃ ∈ [2π/3, π],  τ ∈ [−1, 1],  λ ∈ simplex,  with θ₁,₂ = (2π − θ₃)/2 ± τ(3θ₃ − 2π)/2.
This is exactly the set of local data with θ₃ = max θᵢ. Every other configuration is
equivalent to one of these by cyclic relabelling and, if needed, a reflection; all four
conditions are invariant under both.

The domain is covered by boxes, processed by branch and bound. For each box we compute
rigorous enclosures of ζ², the Jacobian, Ǐ, the turning rate and the cusp indicator
conj(t)²ζ. These are mean-value forms f(c) + ∇f(B)·(B − c), with ∇f(B) from interval
forward-mode AD, intersected with the naive enclosure. Arithmetic is IEEE double with
one-ulp outward rounding per operation. Endpoint values of sin and cos are enclosed by Taylor
series, after range reduction by an interval enclosure of π/2, with an explicit Lagrange
remainder; interior extrema are handled explicitly. No libm accuracy is assumed.

A box is closed only when one of the following holds, each a positively established
strict inequality (so NaN can never certify):
- |ζ|² > κ² on the whole box;
- |ζ|² < κ² and (J) holds;
- (J), (I), (T) and (W) all hold, where (T) may be replaced by
  "no cusp is possible in the box": Im(conj(t)²ζ) ≠ 0 or Re(conj(t)²ζ) > 0.

**Result for s = 1.** 9.97·10⁶ boxes, 0 unresolved, 402 s on 16 cores (`cert_s1.0.json`).
An independent rerun from a different initial subdivision also left 0 unresolved boxes
(`cert_s1.0_n10.json`).

Validation:
- the closed forms were checked against an independent linear-solve implementation, to
  10⁻¹⁴, and against finite differences along numerically traced ∂M of concrete triangles;
- the interval enclosures contain the float values at ~3·10⁵ random points
  (`tests/test_enclosures.py`).

## 8. Mixed signs (partial)

Let q = (+,+,−) (the other mixed case follows by global sign flip). Then:
- Equilibria lie in R_A = {λ₁,λ₂ > 0 > λ₃} (beyond edge p₁p₂) or in
  R_B = {λ₁,λ₂ < 0 < λ₃} (the vertical angle at p₃).
- **Maxima are impossible.** In R_B, p₃ lies inside the triangle (x, p₁, p₂). For e ⊥ u₃,
  eᵀTe = λ₁(u₁·e)² + λ₂(u₂·e)² ≤ 0, so λ_min(T) ≤ 0 < 1/(s+2).
- Poincaré–Hopf gives #equilibria = 2 + 2·#minima if Σq ≠ 0, and 1 + 2·#minima if Σq = 0.
  So Maxwell's bound ⇔ at most one minimum in M_A = R_A ∩ {|Z| < κ}.
- In the local chart the same closed forms hold, with rᵢ = c·σᵢ/λᵢ, c = −1, and
  Ǐ = ½σ₁σ₂ + c·{…}.
- Numerically (`mixed_caustic.py`): M_A is a single star-shaped disk, and the caustic
  winding is ≤ 1 for every triangle tried.
- The degenerate limit is near-collinear (+,−,+). There each ring of equilibria of the
  collinear configuration leaves exactly 2 in-plane equilibria under perturbation; for
  random samples, axis points + 2·rings never exceeded 4 (`collinear_mixed.py`).
- A certificate for the mixed case needs a compactified chart (λ is unbounded here) and
  is **not done yet**.

### 8.1 Mixed signs: numerical evidence that the same proof architecture applies

Measured with `generic_invariants.py` (chart-free formulas, validated against the closed forms
on positive configurations to 10⁻¹⁵) along ∂M_A for 38 random triangles, including very flat
ones and both kinds of near-collinear limits (`mixed_margins.py`):
- jac < 0 everywhere: Z is orientation-reversing on M_A;
- the inflection function never changes sign: no inflections;
- the turning rate is ≤ −0.0 everywhere, far below the ½ threshold;
- exactly 3 cusps on every caustic.

So the hedgehog lemma would apply verbatim, which would give ≤ 1 minimum and hence ≤ 4
equilibria for mixed signs as well.

What a certificate still needs:
1. **Non-compactness of the chart.** λ is unbounded on R_A. Three degenerate regimes must be
   blown up, with rescaled normal forms:
   - flat (+,−,+), with the negative charge between the positive ones;
   - flat (+,+,−), with the negative charge at an end;
   - the far-field regime Σq ≈ 0, where the minimum can sit at L = 1 − 2λ₃ ~ 10³.

   The inflection margin degenerates linearly in the flatness, as in the positive case.
2. **Total turning π.** The shortcut (W) fails here: |d arg u_j/dω| reaches ½ in flat limits
   and exceeds it in the far-field regime. Θ = π must instead come from continuity over the
   connected space of mixed configurations, since Θ is an integer that cannot jump without an
   inflection, together with an evaluation at one base configuration.
