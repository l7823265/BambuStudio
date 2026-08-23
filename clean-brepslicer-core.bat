@echo off
setlocal
set SRC=E:\learning\slicer\BrepSlicer-core-L2-v3
set CORE=E:\learning\slicer\BambuStudio\src\brepslicer\core
set BREP=E:\learning\slicer\BambuStudio\src\brepslicer

echo ==== clean non-code brepslicer core ====

rem purge non-code from sibling L2-v3
for %%D in (third_party _zip_extract dist example out_uvmatch build _bin scripts apps) do (
  if exist "%SRC%\%%D" rmdir /s /q "%SRC%\%%D"
)
del /q "%SRC%\*.md" "%SRC%\PACKAGE.txt" "%SRC%\*.pdf" "%SRC%\build.bat" 2>nul

rem purge non-code from vendored core
if exist "%CORE%" rmdir /s /q "%CORE%"
mkdir "%CORE%\include"
mkdir "%CORE%\src"
if exist "%SRC%\cmake" (
  mkdir "%CORE%\cmake"
  robocopy "%SRC%\cmake" "%CORE%\cmake" /E /NFL /NDL /NJH /NJS /nc /ns /np
)
robocopy "%SRC%\include" "%CORE%\include" /E /NFL /NDL /NJH /NJS /nc /ns /np
robocopy "%SRC%\src" "%CORE%\src" /E /XF *.tmp /NFL /NDL /NJH /NJS /nc /ns /np
if exist "%SRC%\CMakeLists.txt" copy /Y "%SRC%\CMakeLists.txt" "%CORE%\CMakeLists.txt"
if exist "%SRC%\CMakePresets.json" copy /Y "%SRC%\CMakePresets.json" "%CORE%\CMakePresets.json"

for %%D in (core.bak_before_l2v3 core.bak_before_restore) do (
  if exist "%BREP%\%%D" rmdir /s /q "%BREP%\%%D"
)

if not exist "%CORE%\include\brepslicer\Slicer.h" (
  echo CLEAN_FAIL
  exit /b 1
)
echo CLEAN_OK
exit /b 0
