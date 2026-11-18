@echo off
rem ifacekit: drop a structure file (.cif or .pdb) on this file, or double-click it and type the path.
rem Makes an HTML report of every chain next to the file and opens it. A file with several models (an NMR ensemble) is analysed model by model.
setlocal EnableExtensions
set "HOME_DIR=%LOCALAPPDATA%\ifacekit"
set "VENV=%HOME_DIR%\venv"
if not exist "%VENV%\Scripts\ifacekit.exe" (
  echo ifacekit is not installed yet. Double-click launchers\windows\ifacekit-portal.bat once, close it, then try again.
  goto end
)
set "F=%~1"
if "%F%"=="" set /p F=Path of the structure file ^(or drag it into this window and press Enter^): 
set "F=%F:"=%"
if not exist "%F%" ( echo File not found: %F% & goto end )
for %%i in ("%F%") do set "OUT=%%~dpni_ifacekit_report.html" & set "RAMA=%%~dpni_ramachandran.html"
for /f %%n in ('"%VENV%\Scripts\python.exe" -c "import gemmi,sys;print(len(gemmi.read_structure(sys.argv[1])))" "%F%"') do set "NM=%%n"
if "%NM%"=="" set "NM=1"
if %NM% GTR 1 (
  echo %NM% models found: analysing every model.
  "%VENV%\Scripts\ifacekit.exe" "%F%" --chains all --models all --html "%OUT%" --rama-html "%RAMA%" --workers 4 --quiet
) else (
  "%VENV%\Scripts\ifacekit.exe" "%F%" --chains all --html "%OUT%" --rama-html "%RAMA%" --quiet
)
if errorlevel 1 ( echo The analysis failed, see the message above. & goto end )
echo Report: %OUT%
start "" "%OUT%"
:end
pause
