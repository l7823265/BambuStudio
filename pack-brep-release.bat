@echo off
setlocal EnableDelayedExpansion
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat" || exit /b 1

set ROOT=E:\learning\slicer\BambuStudio
set BIN=%ROOT%\build\src\RelWithDebInfo
set STAGE=%ROOT%\dist\bambu-studio-brep-release
set OUT=%ROOT%\dist\BambuStudio-brep-L2v3-%DATE:~0,4%%DATE:~5,2%%DATE:~8,2%.zip
set OUT=%OUT:/=%
set OUT=%OUT: =%

echo ==== pack release %DATE% %TIME% ====
if not exist "%BIN%\bambu-studio.exe" (
  echo MISSING_BUILD: %BIN%\bambu-studio.exe
  exit /b 1
)
if not exist "%BIN%\BambuStudio.dll" (
  echo MISSING_BUILD: %BIN%\BambuStudio.dll
  exit /b 2
)
if not exist "%BIN%\resources" (
  echo MISSING_RESOURCES: %BIN%\resources
  exit /b 3
)

if exist "%STAGE%" rmdir /s /q "%STAGE%"
mkdir "%STAGE%"

copy /Y "%BIN%\bambu-studio.exe" "%STAGE%\" >nul
copy /Y "%BIN%\BambuStudio.dll" "%STAGE%\" >nul
for %%F in ("%BIN%\*.dll") do (
  set "name=%%~nxF"
  if /I not "!name!"=="BambuStudio.dll" copy /Y "%%F" "%STAGE%\" >nul
)
robocopy "%BIN%\resources" "%STAGE%\resources" /E /NFL /NDL /NJH /NJS /nc /ns /np >nul

if exist "%OUT%" del /F /Q "%OUT%"
powershell -NoProfile -Command "Compress-Archive -Path '%STAGE%\*' -DestinationPath '%OUT%' -Force"
if not exist "%OUT%" (
  echo ZIP_FAIL
  exit /b 4
)

for %%A in ("%OUT%") do echo PACK_OK %%~fA size=%%~zA bytes
exit /b 0
