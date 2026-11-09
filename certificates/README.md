# Recorded runs

Each certificate JSON records SHA-256 hashes of the source files that produced it. The
sources were changed on 10 November 2026 (range check for sin and cos, removal of an unused
legacy code path), so the directories below refer to two different sets of sources.

| directory | sources | contents |
|---|---|---|
| [`review-2026-11-10/`](review-2026-11-10/) | current | full `reproduce.py` run: both test suites, the n0=12 and n0=10 certificates, the base triangle at n0=64 and n0=48; [`SOURCE_HASHES.md`](review-2026-11-10/SOURCE_HASHES.md); runs of the second implementation in `independent/` |
| [`review-2026-10-30/`](review-2026-10-30/) | commit `1d7c1f8` | the same set of runs, recorded during the first review |
| [`rerun-2026-10-31/`](rerun-2026-10-31/) | commit `1d7c1f8` | a repeat of the n0=12 run in a separate session |

The box counts and the certified bounds of the main runs are identical in all three
directories; only the timings and the source hashes differ.

The three JSON files in the repository root (`cert_s1.0.json`, `cert_s1.0_n10.json`,
`cert_base_triangle.json`) are copies of the n0=12, n0=10 and n0=64 results in
`review-2026-11-10/`.

[`HISTORY.json`](HISTORY.json) names the directories recorded with earlier sources and the
commit whose files their hashes refer to. The release check in
`publication/prepare_bundle.py` uses it to compare those records with the files of that
commit instead of the current ones.
