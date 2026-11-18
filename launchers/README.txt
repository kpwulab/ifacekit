Double-click launchers for ifacekit (no command line needed)

These files start ifacekit without typing commands. The first start makes a private Python environment (Windows: %LOCALAPPDATA%\ifacekit,
macOS: ~/.ifacekit) and installs ifacekit into it, which needs Python 3.9 or newer and an internet connection once. Later starts take a few seconds.
Nothing is installed anywhere else, and the web portal listens only on this computer (http://127.0.0.1:8000).

WINDOWS  (folder launchers\windows)
  ifacekit-portal.bat      double-click: starts the web portal and opens it in the browser; close the window to stop it
  ifacekit-check-file.bat  drop a .cif or .pdb file on it (or double-click and type the path): HTML report of every chain next to the file, opened for you
  Python: if it is missing, install it from https://www.python.org/downloads/ and tick "Add python.exe to PATH".
  Windows may say "Windows protected your PC" for a downloaded .bat file: choose More info, then Run anyway.

MACOS  (folder launchers/mac)
  ifacekit-portal.command      double-click: starts the web portal and opens it in the browser; close the window to stop it
  ifacekit-check-file.command  double-click, drag a .cif or .pdb file into the window, press Enter: HTML report of every chain next to the file
  macOS blocks downloaded files the first time: right-click the file, choose Open, then Open again. If it says that the file cannot be
  executed, open Terminal once and run:   chmod +x launchers/mac/*.command   (and, for a quarantined download:   xattr -dr com.apple.quarantine <the ifacekit folder>)
  Python: macOS asks to install the developer tools when python3 is first used; accept, or install Python from https://www.python.org/downloads/.

A file with several models (an NMR ensemble) is analysed model by model. Everything else (options, the command line) is in README.md and docs/USAGE.md.
These launchers were written for Windows 10/11 and macOS 12 or newer; the author's test machine is Linux, so please report problems.
