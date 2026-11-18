#!/bin/bash
# ifacekit: start the web portal on this computer and open it in the browser (double-click this file in Finder).
# First start: makes a private Python environment in ~/.ifacekit and installs ifacekit into it (internet needed once).
cd "$(dirname "$0")/../.." || exit 1
HOME_DIR="$HOME/.ifacekit"; VENV="$HOME_DIR/venv"; PORT="${IFACEKIT_PORT:-8000}"
pause() { echo; read -n 1 -s -r -p "Press any key to close this window ..."; echo; }

PY="${IFACEKIT_PYTHON:-}"
if [ -z "$PY" ]; then for c in python3 python; do if command -v "$c" >/dev/null 2>&1; then PY="$c"; break; fi; done; fi
if [ -z "$PY" ] || ! "$PY" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)' 2>/dev/null; then
  echo "Python 3.9 or newer was not found."
  echo "Install it from https://www.python.org/downloads/ (or: xcode-select --install), then double-click this file again."
  pause; exit 1
fi
WANT=$(sed -n 's/^__version__ = "\(.*\)"/\1/p' ifacekit/__init__.py)
if [ ! -x "$VENV/bin/python" ]; then
  echo "Making the Python environment, this happens once ..."
  mkdir -p "$HOME_DIR"
  "$PY" -m venv "$VENV" || { echo "The Python environment could not be made."; pause; exit 1; }
fi
HAVE=$(cat "$HOME_DIR/installed_version.txt" 2>/dev/null)
if [ "$HAVE" != "$WANT" ]; then
  echo "Installing ifacekit $WANT (internet needed, a minute or two) ..."
  "$VENV/bin/python" -m pip install --quiet --upgrade pip
  "$VENV/bin/python" -m pip install --quiet ".[server]" || { echo "The installation did not work. Read the messages above."; pause; exit 1; }
  echo "$WANT" > "$HOME_DIR/installed_version.txt"
fi
if [ -n "$IFACEKIT_LAUNCHER_TEST" ]; then echo "ready: ifacekit $WANT"; exit 0; fi

export IFACEKIT_DATA="$HOME_DIR/jobs"
echo
echo "ifacekit $WANT is running at http://127.0.0.1:$PORT"
echo "Keep this window open while you work; close it (or press Control-C) to stop ifacekit."
"$VENV/bin/python" server/app.py --port "$PORT" &
SERVER=$!
for _ in $(seq 1 30); do curl -s -o /dev/null "http://127.0.0.1:$PORT/healthz" && break; sleep 1; done
if command -v open >/dev/null 2>&1; then open "http://127.0.0.1:$PORT"; elif command -v xdg-open >/dev/null 2>&1; then xdg-open "http://127.0.0.1:$PORT" >/dev/null 2>&1; fi
wait "$SERVER"
