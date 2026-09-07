@echo off
setlocal
set ZIP=%~1
if "%ZIP%"=="" set ZIP=E:\learning\slicer\BrepSlicer-core-L2-v3\dist\BrepSlicer-core-L2-src-20260907.zip
set SRC=E:\learning\slicer\BrepSlicer-core-L2-v3
set EXTRACT=%SRC%\_zip_extract
set CORE=E:\learning\slicer\BambuStudio\src\brepslicer\core
set LOG=E:\learning\slicer\BambuStudio\integrate-l2-v3.log

echo ==== integrate L2-v3 %DATE% %TIME% ==== > "%LOG%"
echo ZIP=%ZIP% >> "%LOG%"

if not exist "%ZIP%" (
  echo MISSING_ZIP >> "%LOG%"
  echo MISSING_ZIP: %ZIP%
  exit /b 1
)

if exist "%EXTRACT%" rmdir /s /q "%EXTRACT%"
mkdir "%EXTRACT%"
powershell -NoProfile -Command "Expand-Archive -Path '%ZIP%' -DestinationPath '%EXTRACT%' -Force" >> "%LOG%" 2>&1
if not exist "%EXTRACT%\include\brepslicer\Slicer.h" (
  echo EXTRACT_FAIL >> "%LOG%"
  exit /b 2
)

rem --- sync code-only into sibling L2-v3 (no third_party, examples, outputs, docs) ---
if not exist "%SRC%\include" mkdir "%SRC%\include"
if not exist "%SRC%\src" mkdir "%SRC%\src"
if not exist "%SRC%\cmake" mkdir "%SRC%\cmake"
robocopy "%EXTRACT%\include" "%SRC%\include" /E /NFL /NDL /NJH /NJS /nc /ns /np >> "%LOG%" 2>&1
robocopy "%EXTRACT%\src" "%SRC%\src" /E /XF *.tmp /NFL /NDL /NJH /NJS /nc /ns /np >> "%LOG%" 2>&1
if exist "%EXTRACT%\cmake" robocopy "%EXTRACT%\cmake" "%SRC%\cmake" /E /NFL /NDL /NJH /NJS /nc /ns /np >> "%LOG%" 2>&1
if exist "%EXTRACT%\CMakeLists.txt" copy /Y "%EXTRACT%\CMakeLists.txt" "%SRC%\CMakeLists.txt" >> "%LOG%" 2>&1
if exist "%EXTRACT%\CMakePresets.json" copy /Y "%EXTRACT%\CMakePresets.json" "%SRC%\CMakePresets.json" >> "%LOG%" 2>&1
echo SYNC_CODE_TO_L2V3=OK >> "%LOG%"

if not exist "%SRC%\include\brepslicer\Slicer.h" (
  echo L2V3_SYNC_FAIL >> "%LOG%"
  exit /b 3
)
if not exist "%SRC%\src\intersect\UvMatch.cpp" (
  echo MISSING_UvMatch >> "%LOG%"
  exit /b 4
)

rem --- purge non-code from sibling L2-v3 (keep dist/ for release zips) ---
for %%D in (third_party _zip_extract example out_uvmatch build _bin scripts apps) do (
  if exist "%SRC%\%%D" rmdir /s /q "%SRC%\%%D"
)
del /q "%SRC%\*.md" "%SRC%\PACKAGE.txt" "%SRC%\*.pdf" "%SRC%\build.bat" 2>nul

rem --- fresh code-only vendored core for Studio ---
if exist "%CORE%" rmdir /s /q "%CORE%"
mkdir "%CORE%\include"
mkdir "%CORE%\src"
if exist "%SRC%\cmake" (
  mkdir "%CORE%\cmake"
  robocopy "%SRC%\cmake" "%CORE%\cmake" /E /NFL /NDL /NJH /NJS /nc /ns /np >> "%LOG%" 2>&1
)
robocopy "%SRC%\include" "%CORE%\include" /E /NFL /NDL /NJH /NJS /nc /ns /np >> "%LOG%" 2>&1
robocopy "%SRC%\src" "%CORE%\src" /E /XF *.tmp /NFL /NDL /NJH /NJS /nc /ns /np >> "%LOG%" 2>&1
if exist "%SRC%\CMakeLists.txt" copy /Y "%SRC%\CMakeLists.txt" "%CORE%\CMakeLists.txt" >> "%LOG%" 2>&1
if exist "%SRC%\CMakePresets.json" copy /Y "%SRC%\CMakePresets.json" "%CORE%\CMakePresets.json" >> "%LOG%" 2>&1
echo ROBOCOPY_TO_STUDIO=OK >> "%LOG%"

rem --- remove old backups and stray non-code under brepslicer ---
set BREP=E:\learning\slicer\BambuStudio\src\brepslicer
for %%D in (core.bak_before_l2v3 core.bak_before_restore) do (
  if exist "%BREP%\%%D" rmdir /s /q "%BREP%\%%D"
)

if not exist "%CORE%\include\brepslicer\Slicer.h" (
  echo STUDIO_SYNC_FAIL >> "%LOG%"
  exit /b 5
)
if not exist "%CORE%\src\intersect\UvMatch.cpp" (
  echo STUDIO_MISSING_UvMatch >> "%LOG%"
  exit /b 6
)

echo INTEGRATE_OK >> "%LOG%"
echo INTEGRATE_OK
type "%LOG%"
exit /b 0
