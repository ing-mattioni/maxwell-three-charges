# KaTeX assets

Unmodified KaTeX 0.16.11 JavaScript, CSS and the WOFF2 font subset used by this
manuscript. Source distribution: <https://github.com/KaTeX/KaTeX/releases/tag/v0.16.11>.
The files were recovered from the browser cache of the paper's original build,
which requested `https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/`.

The upstream MIT license is retained in `LICENSE`. `SHA256SUMS.json` records
the exact assets used in the reviewed build; it is not an upstream signature.

The CSS references other fallback fonts, which are not included. All resources
required by the current manuscript are present. `paper/build.cjs` rejects failed
resource requests and blocks network access. A future change that needs an
additional font must vendor that font from the same release before building.
