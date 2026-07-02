<#
Install helper: copy Prodigal and HMMER from Downloads into C:\tools if present.
If not found, prints instructions where to download them.
Usage: run from PowerShell as admin if writing to C:\tools
#>
param(
    [string]$Downloads = "$env:USERPROFILE\Downloads",
    [string]$ToolsDir = 'C:\tools'
)

New-Item -ItemType Directory -Path $ToolsDir -Force | Out-Null

# Prodigal
$prodCandidate = Join-Path $Downloads 'prodigal.exe'
if (Test-Path $prodCandidate) {
    Copy-Item $prodCandidate -Destination $ToolsDir -Force
    Write-Output "Copied prodigal.exe -> $ToolsDir"
} else {
    Write-Output "prodigal.exe not found in $Downloads. Download from https://github.com/hyattpd/Prodigal/releases and place prodigal.exe into $Downloads then re-run this script."
}

# HMMER
$hmmerZip = Get-ChildItem -Path $Downloads -Filter 'hmmer*.zip' -File -ErrorAction SilentlyContinue | Select-Object -First 1
if ($hmmerZip) {
    $extractDir = Join-Path $ToolsDir 'hmmer_tmp'
    Expand-Archive $hmmerZip.FullName -DestinationPath $extractDir -Force
    # try to find bin\hmmsearch.exe
    $hmmSearchPath = Get-ChildItem -Path $extractDir -Recurse -Filter 'hmmsearch.exe' -File -ErrorAction SilentlyContinue | Select-Object -First 1
    if ($hmmSearchPath) {
        Copy-Item $hmmSearchPath.FullName -Destination $ToolsDir -Force
        Copy-Item (Join-Path $hmmSearchPath.DirectoryName 'hmmscan.exe') -Destination $ToolsDir -Force -ErrorAction SilentlyContinue
        Write-Output "Copied hmmsearch.exe (and hmmscan.exe if present) -> $ToolsDir"
    } else { Write-Output "Could not find hmmsearch.exe inside $($hmmerZip.FullName)." }
    Remove-Item -Recurse -Force $extractDir
} else {
    Write-Output "No hmmer zip found in $Downloads. Download hmmer-3.3.2-win64.zip from http://hmmer.org/download.html and put it into $Downloads, then re-run this script."
}

Write-Output "Done. Verify with: Get-ChildItem $ToolsDir`\n & 'C:\tools\prodigal.exe' -h`\n & 'C:\tools\hmmsearch.exe' -h"