# Runs recorded on 10 November 2026, after the check of the arithmetic environment

These runs repeat those of [`../review-2026-11-10/`](../review-2026-11-10/) after the last
change suggested by the second review (rilievo R4 in
[`paper/REVIEW-2026-11-09.md`](../../paper/REVIEW-2026-11-09.md), in Italian). `interval.py` and
`independent/evaluator.py` now check, when they are imported, that NumPy's binary64 arithmetic
rounds to nearest with ties to even and keeps subnormal numbers. Each module evaluates
operations whose correctly rounded results are known exactly (18 in `interval.py`, 15 in
`evaluator.py`, different ones), compares the results bit for bit on scalars and on arrays, and
refuses to load if any differs. `tests/test_rigorous_arithmetic.py` checks the expected results
with exact rationals and, on Windows, that the import fails under the three other rounding
modes and under flush-to-zero or denormals-are-zero. Nothing else in the computation changed.
The hashes of the sources are in [`SOURCE_HASHES.md`](SOURCE_HASHES.md).

The runs were made by Claude (Opus 5.5) in the author's repository, at the author's request, on
the machine used for the earlier records: Windows, CPython 3.11.9, NumPy 2.4.6, mpmath 1.3.0,
16 worker processes.

## Main verifier

```powershell
python reproduce.py --workers 16 --output certificates/rerun-2026-11-10
```

| directory | run | boxes processed | unresolved | time |
|---|---|---:|---:|---:|
| `arithmetic/` | `tests/test_rigorous_arithmetic.py`, with the test of the environment check | – | – | 2.0 s |
| `enclosures/` | `tests/test_enclosures.py` | – | – | 6.0 s |
| `n12/` | `certify.py 1.0 16 60` | 9,965,924 | 0 | 510.3 s |
| `n10/` | `certify.py 1.0 16 80 10` | 9,966,320 | 0 | 490.7 s |
| `base64/` | `base_triangle.py 64` | 2,356 | 0 | < 1 s |
| `base48/` | `base_triangle.py 48` | 1,618 | 0 | < 1 s |

`reproduction.json` is the runner's own record, with the validation of each step. The box
statistics and the certified bounds of `n12/` and `n10/` are identical, number for number,
to those of all earlier records.

## Second implementation

```powershell
python independent/certify_domain.py 16      # independent_domain.json, domain.log
python independent/certify_base.py           # independent_base.json, base.log
python independent/validate.py               # validate.log
python independent/first_principles.py       # first_principles.log
```

| file | run | result |
|---|---|---|
| `independent_domain.json`, `domain.log` | items 1–2 of Theorem 5.1 on the whole chart domain | 84,329,540 boxes, 0 unresolved, 841 s |
| `independent_base.json`, `base.log` | item 3 of Theorem 5.1 on the base triangle | 1,279 boxes, 0 unresolved |
| `validate.log` | enclosures of the second evaluator against 40-digit values | no violation |
| `first_principles.log` | closed forms of Section 6 against the definitions, 50 digits | agreement to about 30 digits |

The counts of the closed boxes, by kind, are identical to those of the first runs in
`../review-2026-11-10/independent/`.
