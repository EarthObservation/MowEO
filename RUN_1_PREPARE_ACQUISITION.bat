@echo off
cd /d "%~dp0"
python 00_check_environment.py
if errorlevel 1 goto :error
python 01_prepare_acquisition.py
if errorlevel 1 goto :error
pause
exit /b 0
:error
echo.
echo PREPARATION STOPPED WITH AN ERROR.
pause
exit /b 1
