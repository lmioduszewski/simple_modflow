param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$PytestArgs
)

$repoRoot = Split-Path -Parent $PSScriptRoot
$tempRoot = Join-Path $repoRoot ".pytest-work\custom_tmp_runs"
$runTemp = Join-Path $tempRoot ([guid]::NewGuid().ToString())

New-Item -ItemType Directory -Force -Path $tempRoot | Out-Null

Get-ChildItem $tempRoot -Directory -ErrorAction SilentlyContinue |
    ForEach-Object {
        try {
            Remove-Item -LiteralPath $_.FullName -Recurse -Force -ErrorAction Stop
        } catch {
            # Ignore locked temp folders; a fresh folder is used for this run.
        }
    }

Write-Host "Using pytest temp dir: $runTemp"

Set-Location $repoRoot
$env:SIMPLE_MODFLOW_PYTEST_TMP_ROOT = $runTemp
mf-env
python -m pytest @PytestArgs
exit $LASTEXITCODE
