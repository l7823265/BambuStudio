@echo off
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat" || exit /b 1
cd /d E:\learning\slicer\BambuStudio\build
cmake --build . --config RelWithDebInfo --target brepslicer libslic3r BambuStudio_app_gui -j 8
echo BUILD_EXIT=%ERRORLEVEL%
exit /b %ERRORLEVEL%
