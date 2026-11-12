@echo off
cd /d "%~dp0"
set /p WORKERS=Workers [4]: 
if "%WORKERS%"=="" set WORKERS=4
python 14_run_processing.py --workers %WORKERS%
if errorlevel 1 goto :error
pause
exit /b 0
:error
echo.
echo PROCESSING STOPPED WITH AN ERROR.
pause
exit /b 1
