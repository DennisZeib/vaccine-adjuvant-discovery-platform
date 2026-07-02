<#
PowerShell wrapper to download a genome, run Prodigal, and scan proteins with HMMER.
Requires: prodigal.exe and hmmsearch.exe available in PATH or pass full paths.
Usage examples:
  .\run_prodigal_and_hmmer.ps1
  .\run_prodigal_and_hmmer.ps1 -Prodigal C:\tools\prodigal.exe -HMMSEARCH C:\tools\hmmsearch.exe
#>
param(
    [string]$GenomeUrl = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/585/GCF_000009585.1_ASM958v1/GCF_000009585.1_ASM958v1_genomic.fna.gz",
    [string]$OutDir = ".\genome",
    [string]$Prodigal = "prodigal.exe",
    [string]$HMMSEARCH = "hmmsearch.exe",
    [string]$HMM = ".\hmms\PKS_NRPS.hmm"
)

New-Item -ItemType Directory -Path $OutDir -Force | Out-Null
Set-Location $OutDir

# Download genome gz
$gz = Split-Path -Leaf $GenomeUrl
$gzPath = Join-Path (Get-Location) $gz
if (-not (Test-Path $gzPath)) {
    Write-Output "Downloading genome: $GenomeUrl"
    try { Invoke-WebRequest -Uri $GenomeUrl -OutFile $gzPath -UseBasicParsing -ErrorAction Stop }
    catch { Write-Error "Genome download failed: $_"; exit 1 }
} else { Write-Output "Genome archive exists: $gzPath" }

# Decompress .gz -> .fna
$fna = [System.IO.Path]::ChangeExtension($gzPath, $null)  # remove .gz
if ($gzPath -like "*.gz" -and -not (Test-Path $fna)) {
    Write-Output "Decompressing $gzPath -> $fna"
    try {
        $inF = [System.IO.File]::OpenRead($gzPath)
        $outF = [System.IO.File]::Create($fna)
        $gzip = New-Object System.IO.Compression.GzipStream($inF,[System.IO.Compression.CompressionMode]::Decompress)
        $buffer = New-Object byte[] 4096
        while (($read = $gzip.Read($buffer,0,$buffer.Length)) -gt 0) { $outF.Write($buffer,0,$read) }
        $gzip.Close(); $inF.Close(); $outF.Close()
    } catch { Write-Error "Decompression failed: $_"; exit 1 }
} elseif (Test-Path $fna) { Write-Output "Genome FASTA exists: $fna" } else { Write-Error "No .gz input found and expected $gzPath"; exit 1 }

# Run Prodigal
$prots = Join-Path (Get-Location) "proteins.faa"
$gbk = Join-Path (Get-Location) "genes.gbk"
function Resolve-Executable($name) {
    # If path exists, return it
    if (Test-Path $name) { return (Resolve-Path $name).Path }
    # Try Get-Command for executables on PATH
    try { $cmd = Get-Command $name -ErrorAction Stop; return $cmd.Source } catch { return $null }
}

$prodPath = Resolve-Executable $Prodigal
if (-not $prodPath) { Write-Error "Prodigal executable not found: $Prodigal. Place prodigal.exe in C:\tools or pass full path via -Prodigal."; exit 1 }
Write-Output "Running Prodigal: $prodPath -i $fna -a $prots -o $gbk"
try {
    & $prodPath -i $fna -a $prots -o $gbk -p single
} catch {
    Write-Error "Prodigal execution failed: $_"; exit 1
}
if ($LASTEXITCODE -ne 0) { Write-Error "Prodigal failed with exit code $LASTEXITCODE"; exit 1 }

# Run HMMER
if (-not (Test-Path $HMM)) { Write-Error "HMM file not found: $HMM"; exit 1 }
$hmmPath = Resolve-Executable $HMMSEARCH
if (-not $hmmPath) { Write-Error "hmmsearch executable not found: $HMMSEARCH. Place hmmsearch.exe in C:\tools or pass full path via -HMMSEARCH."; exit 1 }
$hitsTbl = Join-Path (Get-Location) "hits.tbl"
$hitsOut = Join-Path (Get-Location) "hits.out"
Write-Output "Running HMMER: $hmmPath --tblout $hitsTbl $HMM $prots > $hitsOut"
try {
    & $hmmPath --tblout $hitsTbl $HMM $prots > $hitsOut
} catch {
    Write-Error "hmmsearch execution failed: $_"; exit 1
}
if ($LASTEXITCODE -ne 0) { Write-Error "hmmsearch failed with exit code $LASTEXITCODE"; exit 1 }

Write-Output "Done. Outputs: $prots, $gbk, $hitsTbl, $hitsOut"
