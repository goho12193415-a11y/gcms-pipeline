@echo off
echo ============================================
echo   GC-MS Auto-Processing v2.1 - Setup
echo   Copyright (c) 2026 go ho
echo ============================================
echo.
echo Installing Python dependencies...
pip install numpy scipy pandas openpyxl pymzml pyms-nist-search
echo   olefile: read Shimadzu .qgd ; pythonnet: read Thermo .RAW natively
pip install olefile pythonnet
echo.
echo ============================================
echo   Setup complete!
echo   Double-click 软件\启动.bat to launch.
echo ============================================
pause
