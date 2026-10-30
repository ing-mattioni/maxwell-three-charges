# Build paper.pdf from paper.html with headless Microsoft Edge (Chromium).
# KaTeX is loaded from cdn.jsdelivr.net, so a network connection is needed.
# Usage (PowerShell, from any folder):  & .\paper\build.ps1
param([string]$ProfileDir = "$env:TEMP\maxwell-paper-edge")
$edge = "C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
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
if (Test-Path $tmp) {
    Move-Item -Force $tmp (Join-Path $here "paper.pdf")
    Get-Item (Join-Path $here "paper.pdf") | Select-Object Name, Length, LastWriteTime
} else {
    Write-Error "PDF was not produced"
}
