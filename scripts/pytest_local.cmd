@echo off
powershell -ExecutionPolicy Bypass -File "%~dp0pytest_local.ps1" %*
exit /b %ERRORLEVEL%
