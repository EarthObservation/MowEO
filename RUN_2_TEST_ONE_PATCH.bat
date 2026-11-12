@echo off
cd /d "%~dp0"
python 03_acquire.py --run --limit 1 --parallel 1
pause
