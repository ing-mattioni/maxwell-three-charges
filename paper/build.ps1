# Build paper.pdf from paper.html. All assets (KaTeX, fonts, figures) are local.
# With Node.js and playwright-core available, build.cjs is used: it also checks fonts,
# KaTeX errors, unresolved references and equation overflow. Otherwise headless
# Microsoft Edge prints the page directly, without those checks.
# Usage (PowerShell, from any folder):  & .\paper\build.ps1
param([string]$NodeExecutable = 'node',
      [string]$ProfileDir = "$env:TEMP\maxwell-paper-edge-$(Get-Random)")
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$out = Join-Path $here 'paper.pdf'
if (Get-Command $NodeExecutable -ErrorAction SilentlyContinue) {
    & $NodeExecutable (Join-Path $here 'build.cjs')
    if ($LASTEXITCODE -ne 0) { throw "PDF build failed with exit code $LASTEXITCODE" }
} else {
    $edge = "C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
    $url = "file:///" + ($here -replace '\\', '/') + "/paper.html"
    $tmp = Join-Path $here "paper_build.pdf"
    if (Test-Path $tmp) { Remove-Item -Force $tmp }
    $edgeArgs = @("--headless=new", "--disable-gpu", "--no-pdf-header-footer", "--user-data-dir=`"$ProfileDir`"",
                  "--virtual-time-budget=30000", "--print-to-pdf=`"$tmp`"", $url)
    $p = Start-Process -FilePath $edge -ArgumentList $edgeArgs -PassThru -NoNewWindow
    $null = $p.WaitForExit(180000)
    # the launcher can return before the browser has written the file: wait for a stable size
    $last = -1
    for ($i = 0; $i -lt 120; $i++) {
        Start-Sleep -Seconds 1
        if (Test-Path $tmp) {
            $len = (Get-Item $tmp).Length
            if ($len -gt 0 -and $len -eq $last) { break }
            $last = $len
        }
    }
    if (-not (Test-Path $tmp)) { throw "PDF was not produced" }
    Move-Item -Force $tmp $out
}
Get-Item -LiteralPath $out | Select-Object Name, Length, LastWriteTime
