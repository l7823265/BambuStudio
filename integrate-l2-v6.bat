@echo off
setlocal
rem Vendor BrepSlicer L2-v6 code-only into Studio.
rem Does NOT purge apps/build/_bin/dist from the sibling tree.
set SRC=%~1
if "%SRC%"=="" set SRC=E:\learning\slicer\BrepSlicer-core-L2-v6
set CORE=E:\learning\slicer\BambuStudio\src\brepslicer\core
set BREP=E:\learning\slicer\BambuStudio\src\brepslicer
set LOG=E:\learning\slicer\BambuStudio\integrate-l2-v6.log

echo ==== integrate L2-v6 %DATE% %TIME% ==== > "%LOG%"
echo SRC=%SRC% >> "%LOG%"

if not exist "%SRC%\include\brepslicer\Slicer.h" (
  echo MISSING_SRC >> "%LOG%"
  echo MISSING_SRC: %SRC%\include\brepslicer\Slicer.h
  exit /b 1
)
if not exist "%SRC%\src\intersect\UvMatch.cpp" (
  echo MISSING_UvMatch >> "%LOG%"
  exit /b 2
)
if not exist "%SRC%\src\topo\SolidAdjacency.cpp" (
  echo MISSING_SolidAdjacency >> "%LOG%"
  exit /b 3
)
if not exist "%SRC%\src\occ\OccTopology.cpp" (
  echo MISSING_OccTopology >> "%LOG%"
  exit /b 4
)

rem --- fresh code-only vendored core for Studio ---
if exist "%CORE%" rmdir /s /q "%CORE%"
mkdir "%CORE%\include"
mkdir "%CORE%\src"
if exist "%SRC%\cmake" (
  mkdir "%CORE%\cmake"
  robocopy "%SRC%\cmake" "%CORE%\cmake" /E /NFL /NDL /NJH /NJS /nc /ns /np >> "%LOG%" 2>&1
)
robocopy "%SRC%\include" "%CORE%\include" /E /NFL /NDL /NJH /NJS /nc /ns /np >> "%LOG%" 2>&1
robocopy "%SRC%\src" "%CORE%\src" /E /XF *.tmp *.bak /NFL /NDL /NJH /NJS /nc /ns /np >> "%LOG%" 2>&1
if exist "%SRC%\CMakeLists.txt" copy /Y "%SRC%\CMakeLists.txt" "%CORE%\CMakeLists.txt" >> "%LOG%" 2>&1
if exist "%SRC%\CMakePresets.json" copy /Y "%SRC%\CMakePresets.json" "%CORE%\CMakePresets.json" >> "%LOG%" 2>&1
echo ROBOCOPY_TO_STUDIO=OK >> "%LOG%"

rem --- remove old backups under brepslicer ---
for %%D in (core.bak_before_l2v3 core.bak_before_l2v6 core.bak_before_restore) do (
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
if not exist "%CORE%\src\topo\SolidAdjacency.cpp" (
  echo STUDIO_MISSING_SolidAdjacency >> "%LOG%"
  exit /b 7
)
if not exist "%CORE%\src\occ\OccTopology.cpp" (
  echo STUDIO_MISSING_OccTopology >> "%LOG%"
  exit /b 8
)

echo INTEGRATE_OK >> "%LOG%"
echo INTEGRATE_OK
type "%LOG%"
exit /b 0
