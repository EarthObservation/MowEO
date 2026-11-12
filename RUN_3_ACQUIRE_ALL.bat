@echo off
cd /d "%~dp0"
python 03_acquire.py --run --parallel 2
pause
