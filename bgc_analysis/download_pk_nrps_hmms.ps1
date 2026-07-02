# Download PKS/NRPS HMMs into bgc_analysis\hmms
$base = Split-Path -Path $MyInvocation.MyCommand.Definition -Parent
Set-Location $base
New-Item -ItemType Directory -Path hmms -Force | Out-Null
$uri = "https://raw.githubusercontent.com/eggnogdb/eggnog-mapper/master/data/hmms/PKS_NRPS.hmm"
$dest = Join-Path $base "hmms\PKS_NRPS.hmm"
Write-Output "Downloading $uri -> $dest"
try {
    Invoke-WebRequest -Uri $uri -OutFile $dest -UseBasicParsing -ErrorAction Stop
    Write-Output "Downloaded to: $dest"
} catch {
    Write-Error "Download failed: $_. Please download manually from: $uri"
    exit 1
}
