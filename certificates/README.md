# Recorded runs

Each certificate JSON records SHA-256 hashes of the source files that produced it. The
sources changed twice on 10 November 2026, so the directories below refer to three sets of
sources. First, a range check was added for sin and cos and an unused legacy code path was
removed. Then a check of the arithmetic environment was added, which runs when the interval
modules are imported.

| directory | sources | contents |
|---|---|---|
| [`rerun-2026-11-10/`](rerun-2026-11-10/) | current | full `reproduce.py` run: both test suites, the n0=12 and n0=10 certificates, the base triangle at n0=64 and n0=48; [`SOURCE_HASHES.md`](rerun-2026-11-10/SOURCE_HASHES.md); runs of the second implementation in `independent/` |
| [`review-2026-11-10/`](review-2026-11-10/) | commit `3c9a49c` | the same set of runs, before the check of the arithmetic environment |
| [`review-2026-10-30/`](review-2026-10-30/) | commit `1d7c1f8` | the runs of the main verifier, recorded during the first review |
| [`rerun-2026-10-31/`](rerun-2026-10-31/) | commit `1d7c1f8` | a repeat of the n0=12 run in a separate session |

The box counts and the certified bounds of the main runs are identical in all four
directories, and so are the counts of the second implementation in the two directories that
have it. Only the timings and the source hashes differ.

The three JSON files in the repository root (`cert_s1.0.json`, `cert_s1.0_n10.json`,
`cert_base_triangle.json`) are copies of the n0=12, n0=10 and n0=64 results in
`rerun-2026-11-10/`.

[`HISTORY.json`](HISTORY.json) names the directories recorded with earlier sources and the
commit whose files their hashes refer to. The author's release check, kept with the
editorial materials outside this repository, uses it to compare those records with the
files of that commit instead of the current ones.
