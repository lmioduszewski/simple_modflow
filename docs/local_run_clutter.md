# Local Run Clutter

This repo keeps generated runs and scratch work out of the importable package tree, but local runs still
create files you may want to clean up from time to time.

## Main Locations

- `examples/artifacts/`
  Older tracked example-run location from earlier workflows. This should be treated as disposable run output.
- `examples/mf6/artifacts/`
  Example model runs, MP3DU outputs, PEST workspaces, logs, plots, and timestamped run folders.
- `examples/mf6/*_latest.txt`
  Small pointer files that record the latest artifact folder for a workflow.
- `.pytest-work/custom_tmp_runs/`
  Temporary pytest workspace used by `scripts/pytest_local.cmd`.
- `.pytest_tmp/` and `.pytest-tmp/`
  Older pytest temp roots from earlier local runs.
- `manual_mp3du_check*/`
  Ad hoc MP3DU scratch folders from manual validation runs.
- `mp3du_test_tmp/`
  Older MP3DU temp workspace.
- `PstFrom.log` and `examples/mf6/notebooks/PstFrom.log`
  Generated PEST logs from local experiments.
- `.build-tmp/` and `pytest-cache-files*/`
  Local scratch space from packaging/tests.

## Safe Cleanup

When no test run or model run is active, these locations are safe to delete:

- `examples/artifacts/`
- `examples/mf6/artifacts/` subfolders you no longer need
- `examples/mf6/*_latest.txt`
- `.pytest-work/`
- `.pytest_tmp/`
- `.pytest-tmp/`
- `manual_mp3du_check*/`
- `mp3du_test_tmp/`
- `PstFrom.log`
- `examples/mf6/notebooks/PstFrom.log`
- `.build-tmp/`
- `pytest-cache-files*/`

## Useful Commands

Delete all example run artifacts:

```powershell
if (Test-Path .\examples\artifacts) { Get-ChildItem .\examples\artifacts | Remove-Item -Recurse -Force }
Get-ChildItem .\examples\mf6\artifacts | Remove-Item -Recurse -Force
```

Delete old pytest temp workspaces:

```powershell
Remove-Item -LiteralPath .\.pytest-work -Recurse -Force
Remove-Item -LiteralPath .\.pytest_tmp -Recurse -Force
Remove-Item -LiteralPath .\.pytest-tmp -Recurse -Force
```

Delete old MP3DU scratch folders:

```powershell
Get-ChildItem .\manual_mp3du_check* | Remove-Item -Recurse -Force
Remove-Item -LiteralPath .\mp3du_test_tmp -Recurse -Force
```

## Practical Rule

If a folder contains timestamped run names, diagnostics, shapefiles, binary model outputs, or test temp files,
assume it is disposable local run clutter unless you want to keep those results for comparison.
