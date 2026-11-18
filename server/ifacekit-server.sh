#!/bin/bash
# Start the ifacekit portal.  Usage: ifacekit-server.sh [--host 127.0.0.1] [--port 8000] [--data DIR]
# With gunicorn installed it runs as a production server (one worker process with many threads, so that the limit on
# simultaneous analyses is global); otherwise it falls back to Flask's built-in server (development only).
# Settings by environment variable: IFACEKIT_DATA, IFACEKIT_KEEP_DAYS (7), IFACEKIT_PURGE_MINUTES (60),
# IFACEKIT_MAX_JOBS (4), IFACEKIT_QUEUE_WAIT (600 s), IFACEKIT_MEM_MB (8192, 0 = no limit).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="$ROOT/.venv/bin/python"; [ -x "$PY" ] || PY="$(command -v python3)"
HOST=127.0.0.1; PORT=8000
while [ $# -gt 0 ]; do
  case "$1" in
    --host) HOST="$2"; shift 2;;
    --port) PORT="$2"; shift 2;;
    --data) export IFACEKIT_DATA="$2"; shift 2;;
    *) echo "unknown option $1" >&2; exit 2;;
  esac
done
export IFACEKIT_DATA="${IFACEKIT_DATA:-$ROOT/jobs}"
mkdir -p "$IFACEKIT_DATA"
if "$PY" -c "import gunicorn" 2>/dev/null; then
  exec "$PY" -m gunicorn --chdir "$HERE" -w 1 --threads "${IFACEKIT_THREADS:-16}" --timeout 1500 \
       --access-logfile - -b "$HOST:$PORT" app:app
else
  echo "gunicorn not installed: using Flask's development server" >&2
  exec "$PY" "$HERE/app.py" --host "$HOST" --port "$PORT" --data "$IFACEKIT_DATA"
fi
