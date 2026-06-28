@echo off
cd /d "%~dp0"
del error.log 2>nul
python gui.py
if exist error.log (
    type error.log
    pause
)
