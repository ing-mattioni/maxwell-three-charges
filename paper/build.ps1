# Build paper.pdf from paper.html. All assets (KaTeX, fonts, figures) are local.
# If Node.js can load playwright-core, build.cjs is used: it also checks fonts, KaTeX
# errors, unresolved references and equation overflow. Otherwise headless Microsoft
# Edge prints the page directly, without those checks. Either way paper.pdf is
# replaced only by a complete PDF.
# Usage (PowerShell, from any folder):  & .\paper\build.ps1
#   Node need not be on PATH:  & .\paper\build.ps1 -NodeExecutable C:\path\to\node.exe
#   (set PLAYWRIGHT_MODULE to the playwright-core folder if it is not installed globally)
param([string]$NodeExecutable = 'node',
      [string]$ProfileDir = "$env:TEMP\maxwell-paper-edge-$(Get-Random)",
      [int]$TimeoutSeconds = 180)
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$out = Join-Path $here 'paper.pdf'

$usePlaywright = $false
if (Get-Command $NodeExecutable -ErrorAction SilentlyContinue) {
    $module = if ($env:PLAYWRIGHT_MODULE) { $env:PLAYWRIGHT_MODULE } else { 'playwright-core' }
    # a failed require writes to stderr; do not let that become a terminating error
    $ErrorActionPreference = 'Continue'
    & $NodeExecutable -e "require.resolve(process.argv[1])" $module 2>$null | Out-Null
    $usePlaywright = ($LASTEXITCODE -eq 0)
    $ErrorActionPreference = 'Stop'
    if (-not $usePlaywright) { Write-Warning "playwright-core not found; using Edge without layout checks" }
}

if ($usePlaywright) {
    # build.cjs validates the page first and writes the PDF only if every check passes
    & $NodeExecutable (Join-Path $here 'build.cjs')
    if ($LASTEXITCODE -ne 0) { throw "PDF build failed with exit code $LASTEXITCODE; paper.pdf left unchanged" }
} else {
    $edge = "C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
    $url = "file:///" + ($here -replace '\\', '/') + "/paper.html"
    $tmp = Join-Path $here "paper_build.pdf"
    if (Test-Path $tmp) { Remove-Item -Force $tmp }
    $edgeArgs = @("--headless=new", "--disable-gpu", "--no-pdf-header-footer", "--user-data-dir=`"$ProfileDir`"",
                  "--virtual-time-budget=30000", "--print-to-pdf=`"$tmp`"", $url)
    $p = Start-Process -FilePath $edge -ArgumentList $edgeArgs -PassThru -NoNewWindow
    if (-not $p.WaitForExit($TimeoutSeconds * 1000)) {
        Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue
        if (Test-Path $tmp) { Remove-Item -Force $tmp }
        throw "Edge did not finish within $TimeoutSeconds s; paper.pdf left unchanged"
    }
    # the launcher can return before the browser has written the file: wait for a stable size
    $stable = $false
    $last = -1
    for ($i = 0; $i -lt $TimeoutSeconds; $i++) {
        Start-Sleep -Seconds 1
        if (Test-Path $tmp) {
            $len = (Get-Item $tmp).Length
            if ($len -gt 0 -and $len -eq $last) { $stable = $true; break }
            $last = $len
        }
    }
    if (-not $stable) {
        if (Test-Path $tmp) { Remove-Item -Force $tmp }
        throw "No complete PDF was written within $TimeoutSeconds s; paper.pdf left unchanged"
    }
    # accept only a complete PDF: %PDF- header and an %%EOF trailer
    $bytes = [IO.File]::ReadAllBytes($tmp)
    $head = [Text.Encoding]::ASCII.GetString($bytes, 0, [Math]::Min(5, $bytes.Length))
    $n = [Math]::Min(1024, $bytes.Length)
    $tail = [Text.Encoding]::ASCII.GetString($bytes, $bytes.Length - $n, $n)
    if ($head -ne '%PDF-' -or -not $tail.Contains('%%EOF')) {
        Remove-Item -Force $tmp
        throw "Edge produced an incomplete PDF; paper.pdf left unchanged"
    }
    Move-Item -Force $tmp $out
}
Get-Item -LiteralPath $out | Select-Object Name, Length, LastWriteTime
