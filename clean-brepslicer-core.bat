@echo off
setlocal
rem Re-vendor code-only core from sibling L2-v6 (does not wipe sibling apps/build).
set SRC=E:\learning\slicer\BrepSlicer-core-L2-v6
call "%~dp0integrate-l2-v6.bat" "%SRC%"
exit /b %ERRORLEVEL%
