<#
Runs the full lightweight BGC discovery pipeline:
- download PKS/NRPS HMMs
- download genome (optional)
- run Prodigal -> proteins.faa
- run HMMER (hmmsearch)
- parse HMMER hits -> candidates.csv/faa
- score candidates with micromol_scorer
- optional BLAST vs MIBiG (requires local MIBiG proteins fasta and blastp)

Usage:
  .\run_full_bgc_discovery.ps1 -Prodigal C:\tools\prodigal.exe -HMMSEARCH C:\tools\hmmsearch.exe
  .\run_full_bgc_discovery.ps1 -RunBlast -Blastp C:\blast\bin\blastp.exe -MibigFasta C:\data\mibig_proteins.faa

Notes: run from project root or the script will `Set-Location` to its folder.
#>
param(
    [string]$Prodigal = "prodigal.exe",
    [string]$HMMSEARCH = "hmmsearch.exe",
    [string]$GenomeUrl = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/585/GCF_000009585.1_ASM958v1/GCF_000009585.1_ASM958v1_genomic.fna.gz",
    [string]$OutDir = "genome",
    [switch]$ForceDownloadGenome = $false,
    [switch]$RunBlast = $false,
    [string]$Blastp = "blastp.exe",
    [string]$MibigFasta = ""
)

$scriptDir = Split-Path -Path $MyInvocation.MyCommand.Definition -Parent
Set-Location $scriptDir

# 1. Download HMMs
Write-Output "Step 1: Downloading PKS/NRPS HMMs"
if ((Test-Path .\hmms\PKS_NRPS.hmm -PathType Leaf -ErrorAction SilentlyContinue) -and (-not $ForceDownloadGenome)) {
    Write-Output "HMM already exists: hmms\PKS_NRPS.hmm"
} else {
    powershell -ExecutionPolicy Bypass -File .\download_pk_nrps_hmms.ps1
    if ($LASTEXITCODE -ne 0) { Write-Error "Failed to download HMMs"; exit 1 }
}

# 2. Run Prodigal + HMMER
Write-Output "Step 2: Running Prodigal + HMMER"
$runArgs = @()
$runArgs += "-Prodigal"; $runArgs += $Prodigal
$runArgs += "-HMMSEARCH"; $runArgs += $HMMSEARCH
$runArgs += "-OutDir"; $runArgs += $OutDir
$runArgs += "-GenomeUrl"; $runArgs += $GenomeUrl

# Call the existing runner
powershell -ExecutionPolicy Bypass -File .\run_prodigal_and_hmmer.ps1 @runArgs
if ($LASTEXITCODE -ne 0) { Write-Error "Prodigal+HMMER step failed"; exit 1 }

# 3. Parse HMMER hits
Write-Output "Step 3: Parsing HMMER hits to candidates"
$hits = Join-Path $OutDir 'hits.tbl'
$prots = Join-Path $OutDir 'proteins.faa'
if (-not (Test-Path $hits) -or -not (Test-Path $prots)) {
    Write-Error "Expected outputs not found: $hits or $prots"
    exit 1
}
python .\parse_hmmer_hits.py --hits $hits --proteins $prots --out "$OutDir/candidates"
if ($LASTEXITCODE -ne 0) { Write-Error "parse_hmmer_hits.py failed"; exit 1 }

# 4. Score candidates
Write-Output "Step 4: Scoring candidate proteins"
Set-Location $scriptDir
python .\score_candidates.py
if ($LASTEXITCODE -ne 0) { Write-Error "score_candidates.py failed"; exit 1 }

# 5. Optional BLAST vs MIBiG
if ($RunBlast) {
    if ((-not $MibigFasta) -or (-not (Test-Path $MibigFasta))) {
        Write-Error "RunBlast requested but MIBiG fasta not provided or not found. Provide -MibigFasta path to local MIBiG proteins fasta."; exit 1
    }
    if (-not (Get-Command $Blastp -ErrorAction SilentlyContinue)) {
        Write-Error "blastp not found at $Blastp. Install BLAST+ and provide path via -Blastp."; exit 1
    }
    $candFaa = Join-Path $OutDir 'candidates.faa'
    if (-not (Test-Path $candFaa)) { Write-Error "Candidate proteins not found: $candFaa"; exit 1 }
    $blastOut = Join-Path $OutDir 'candidates_vs_mibig.tsv'
    Write-Output "Running blastp of candidates vs MIBiG (this may take a while)"
    & $Blastp -query $candFaa -db $MibigFasta -outfmt '6 qseqid sseqid pident length evalue bitscore' -max_target_seqs 10 -out $blastOut
    if ($LASTEXITCODE -ne 0) { Write-Error "blastp failed"; exit 1 }
    Write-Output "BLAST done: $blastOut"
}

Write-Output "All done. Outputs are in: $scriptDir\$OutDir (candidates.faa, candidates.csv, candidate_scores.csv)"
