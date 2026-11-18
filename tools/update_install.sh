#!/bin/bash
# Update an installed ifacekit from its git checkout.
#   tools/update_install.sh            go to the newest release tag (vX.Y.Z) and reinstall
#   tools/update_install.sh v0.6.0     go to a given tag
#   tools/update_install.sh main       follow a branch (development version)
# It refuses to run when you have uncommitted changes. To return to the development branch afterwards:  git checkout main
set -euo pipefail
cd "$(dirname "$0")/.."
[ -d .git ] || { echo "this folder is not a git checkout" >&2; exit 1; }
if [ -n "$(git status --porcelain)" ]; then echo "uncommitted changes: commit or stash them first" >&2; exit 1; fi
git fetch --tags --prune origin
REF="${1:-$(git tag --list 'v[0-9]*' --sort=-v:refname | head -1)}"
[ -n "$REF" ] || { echo "no release tag found: give a tag or a branch name" >&2; exit 1; }
echo "now at $(git describe --tags --always); going to $REF"
git checkout "$REF"
if git show-ref --verify --quiet "refs/remotes/origin/$REF"; then git merge --ff-only "origin/$REF"; fi
PY="${PYTHON:-.venv/bin/python}"; [ -x "$PY" ] || PY="$(command -v python3)"
"$PY" -m pip install -q --upgrade ".[ligand,pdf,server]"
rm -rf build
"$PY" -m ifacekit --version
echo "If the portal runs as a service, restart it (sudo systemctl restart ifacekit) or reload gunicorn (kill -HUP <master process id>)."
