@echo off
rem ifacekit: start the web portal on this computer and open it in the browser (double-click this file).
rem First start: makes a private Python environment in %LOCALAPPDATA%\ifacekit and installs ifacekit into it (internet needed once).
setlocal EnableExtensions
cd /d "%~dp0..\.."
set "HOME_DIR=%LOCALAPPDATA%\ifacekit"
set "VENV=%HOME_DIR%\venv"
set "PORT=8000"
if defined IFACEKIT_PORT set "PORT=%IFACEKIT_PORT%"

set "PY="
if defined IFACEKIT_PYTHON set "PY=%IFACEKIT_PYTHON%"
if not defined PY ( where py >nul 2>nul && set "PY=py -3" )
if not defined PY ( where python >nul 2>nul && set "PY=python" )
if not defined PY goto nopython
%PY% -c "import sys; sys.exit(0 if sys.version_info>=(3,9) else 1)" 2>nul
if errorlevel 1 goto nopython

for /f "usebackq delims=" %%v in (`%PY% -c "import re;print(re.search(r'__version__ = .([0-9.]+).',open(r'ifacekit\__init__.py').read()).group(1))"`) do set "WANT=%%v"
if not exist "%VENV%\Scripts\python.exe" (
  echo Making the Python environment, this happens once ...
  %PY% -m venv "%VENV%"
  if errorlevel 1 goto fail
)
set "HAVE="
if exist "%HOME_DIR%\installed_version.txt" set /p HAVE=<"%HOME_DIR%\installed_version.txt"
if not "%HAVE%"=="%WANT%" (
  echo Installing ifacekit %WANT% ^(internet needed, a minute or two^) ...
  "%VENV%\Scripts\python.exe" -m pip install --quiet --upgrade pip
  "%VENV%\Scripts\python.exe" -m pip install --quiet ".[server]"
  if errorlevel 1 goto fail
  >"%HOME_DIR%\installed_version.txt" echo %WANT%
)
if defined IFACEKIT_LAUNCHER_TEST ( echo ready: ifacekit %WANT% & exit /b 0 )

set "IFACEKIT_DATA=%HOME_DIR%\jobs"
echo.
echo ifacekit %WANT% is running at http://127.0.0.1:%PORT%
echo Keep this window open while you work; close it to stop ifacekit.
start "" cmd /c "timeout /t 4 /nobreak >nul & start http://127.0.0.1:%PORT%"
"%VENV%\Scripts\python.exe" server\app.py --port %PORT%
goto end

:nopython
echo Python 3.9 or newer was not found.
echo Install it from https://www.python.org/downloads/ ^(tick "Add python.exe to PATH"^), then double-click this file again.
goto end
:fail
echo The installation did not work. Read the messages above, or send them to the author.
:end
pause
