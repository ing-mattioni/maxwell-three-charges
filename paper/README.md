# Building and reviewing the paper

`paper.html` is the source and `paper.pdf` the rendered manuscript. The source uses local
KaTeX 0.16.11 assets (`vendor/katex/`), so rendering does not need a network connection.

```powershell
& .\paper\build.ps1
```

If Node.js can load `playwright-core`, `build.ps1` runs `build.cjs`, which rejects
missing fonts or images, JavaScript or KaTeX errors, unresolved references, external
resource requests and overflowing display equations, and writes `paper.layout.json` next to
the PDF. Otherwise it prints the page with headless Microsoft Edge, without those checks.
In both cases the existing `paper.pdf` is replaced only by a complete PDF; on a failure or
timeout the script stops and leaves it unchanged.

Node does not have to be on `PATH`:

```powershell
$env:PLAYWRIGHT_MODULE = 'C:\path\to\node_modules\playwright-core'
& .\paper\build.ps1 -NodeExecutable 'C:\path\to\node.exe'
```

```powershell
npm install --no-save playwright-core@1.57.0
node paper/build.cjs              # or: node paper/build.cjs other-output.pdf
```

`build.cjs` was used with Node.js 24.14.0 and Playwright Core 1.57.0. On Windows it
defaults to Edge's standard installation path. Set `EDGE_EXECUTABLE` to use another
Chromium executable, or `PLAYWRIGHT_MODULE` to point to an existing `playwright-core`.

The automatic checks do not replace visual inspection, and they do not verify the
mathematics. After substantive edits, render the pages with Poppler and look at them:

```powershell
pdftoppm -scale-to 1500 -png paper/paper.pdf path/to/existing-directory/page
```

The certificate runs are recorded in `../certificates/review-2026-10-30/`.
`REVIEW.md` is the report of the AI-assisted review of 30–31 October 2026 (in Italian).
