# A second implementation

These files repeat the computer-assisted part of the proof with a program that shares no
code with the main verifier (`interval.py`, `ad.py`, `verify.py`, `certify.py`,
`base_triangle.py`). They were written on 9 November 2026 by Claude Fable 5.2 (Anthropic),
an AI system different from the one that wrote the main verifier, using the manuscript
only. The review they belong to is reported in
[`../paper/REVIEW-2026-11-09.md`](../paper/REVIEW-2026-11-09.md) (in Italian).

| file | purpose |
|---|---|
| `evaluator.py` | interval type, enclosures of sin and cos, and the rescaled quantities of Section 6 by plain interval evaluation |
| `certify_domain.py` | branch and bound for items 1–2 of Theorem 5.1 on the whole chart domain |
| `certify_base.py` | item 3 of Theorem 5.1 on the base triangle |
| `first_principles.py` | the closed forms of Section 6 against the definitions, in 50-digit arithmetic |
| `validate.py` | the enclosures of `evaluator.py` against 40-digit values |

## What is and is not shared with the main verifier

Not shared: the interval library, the enclosures of sin and cos, automatic differentiation
and mean-value forms (this program has neither), the initial subdivision, the bisection
rule, and the thresholds (|Z|² is compared with an enclosure of 1/9, without the 10⁻¹²
margins).

Shared: IEEE 754 binary64 arithmetic in NumPy with one-ulp outward rounding, and the
formulas of the manuscript. Each program checks that arithmetic on its own when it is
imported (`environment_check()` here, `check_arithmetic()` in `interval.py`, with different
test operations) and refuses to run under another rounding mode or with flush-to-zero. `first_principles.py` checks those formulas at sample points; it
does not prove them. Both programs were written by AI systems. This is a cross-check, not a
formal proof and not a review by a human expert.

## Running

`certify_domain.py` and `certify_base.py` need only NumPy; the other two also need mpmath.
The certifiers write their JSON summary in the current directory and exit with status 1 if
any box is unresolved.

```powershell
python independent/certify_domain.py 16
python independent/certify_base.py
python independent/validate.py
python independent/first_principles.py
```

## Recorded results

The runs are in [`../certificates/rerun-2026-11-10/independent/`](../certificates/rerun-2026-11-10/independent/). The first runs, made before the check of the arithmetic environment was added, are in `../certificates/review-2026-11-10/independent/` and have the same counts.

| run | boxes processed | unresolved |
|---|---:|---:|
| `certify_domain.py 16` (8⁴ initial boxes) | 84,329,540 | 0 |
| `certify_base.py` (one initial box) | 1,279 | 0 |

`certify_domain.py` closes a box when one of the following holds on the whole box, as a
strict inequality between enclosure endpoints:

- λ₃ < 0 (the box lies outside the domain);
- |Z|² > 1/9;
- |Z|² < 1/9 and σ₃²J < 0;
- σ₃²J < 0, Ǐ > 0, |ρ₁ Im t| < 1/2, and either ν < 1/2 or no cusp is possible in the box.

The count 84,329,540 does not depend on the number of workers: every box is processed
independently of the others.
