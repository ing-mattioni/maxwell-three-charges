# Runs recorded on 10 November 2026

> **Historical record.** These runs were made with the sources of commit `3c9a49c`. Later the
> same day a check of the arithmetic environment was added to `interval.py` and
> `independent/evaluator.py`, and every run was repeated with identical counts; the current
> runs are in [`../rerun-2026-11-10/`](../rerun-2026-11-10/).

These runs were made after the second review ([`paper/REVIEW-2026-11-09.md`](../../paper/REVIEW-2026-11-09.md),
in Italian), with the sources as corrected on 10 November 2026: `sin_i` and `cos_i` reject
arguments outside the range on which they are enclosures, the prefilter `zeta2_naive` falls
back to [−1, 1] for intervals that are too wide, and an unused legacy code path was removed.
The hashes of the sources are in [`SOURCE_HASHES.md`](SOURCE_HASHES.md).

They were made by Claude Fable 5.2 in the author's repository, at the author's request, on
the machine used for the earlier records: Windows, CPython 3.11.9, NumPy 2.4.6,
mpmath 1.3.0, 16 worker processes.

## Main verifier

```powershell
python reproduce.py --workers 16 --output certificates/review-2026-11-10
```

| directory | run | boxes processed | unresolved | time |
|---|---|---:|---:|---:|
| `arithmetic/` | `tests/test_rigorous_arithmetic.py`, with the two new regression tests | – | – | 0.5 s |
| `enclosures/` | `tests/test_enclosures.py` | – | – | 5.8 s |
| `n12/` | `certify.py 1.0 16 60` | 9,965,924 | 0 | 464.8 s |
| `n10/` | `certify.py 1.0 16 80 10` | 9,966,320 | 0 | 502.9 s |
| `base64/` | `base_triangle.py 64` | 2,356 | 0 | < 1 s |
| `base48/` | `base_triangle.py 48` | 1,618 | 0 | < 1 s |

`reproduction.json` is the runner's own record, with the validation of each step. The box
statistics and the certified bounds of `n12/` and `n10/` are identical, number for number,
to those in `../review-2026-10-30/`. Only the timings and the source hashes differ.

## Second implementation

The directory `independent/` holds the runs of the program in the top-level directory of
the same name, which shares no code with the main verifier.

```powershell
python independent/certify_domain.py 16      # independent_domain.json, domain.log
python independent/certify_base.py           # independent_base.json, base.log
python independent/validate.py               # validate.log
python independent/first_principles.py       # first_principles.log
```

| file | run | result |
|---|---|---|
| `independent_domain.json`, `domain.log` | items 1–2 of Theorem 5.1 on the whole chart domain | 84,329,540 boxes, 0 unresolved |
| `independent_base.json`, `base.log` | item 3 of Theorem 5.1 on the base triangle | 1,279 boxes, 0 unresolved |
| `validate.log` | enclosures of the second evaluator against 40-digit values | no violation |
| `first_principles.log` | closed forms of Section 6 against the definitions, 50 digits | agreement to about 30 digits |

The first three were first run on 9 November 2026, during the review, from working copies
of the same code; the box counts are the same.
