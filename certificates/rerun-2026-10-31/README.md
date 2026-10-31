# Repeat of run A (n0 = 12), 31 October 2026

A repeat of the main certificate run, made by Claude in a session separate from the
review, to check that the recorded result reproduces.

- Command: `python certify.py 1.0 16 60`, run from a copy of the working tree in a scratch
  directory, so that the files in the repository root were not overwritten.
- Runtime: CPython 3.11.9, NumPy 2.4.6, Windows, 16 worker processes; 493 s.
- Sources: the SHA-256 hashes in `n12/cert_s1.0.json` are identical to those recorded in
  `../review-2026-10-30/reproduction.json` (see `../review-2026-10-30/SOURCE_HASHES.md`).
- Result: 9,965,924 boxes, 0 unresolved, exit status 0. The box statistics are identical,
  count for count, to the recorded run in `../review-2026-10-30/n12/`.

This is the same implementation on the same machine. It confirms reproducibility, not the
correctness of the implementation.
