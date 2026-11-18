#!/bin/bash
# ifacekit: double-click, drag a structure file (.cif or .pdb) into the window and press Enter (or give the file as the first argument).
# Makes an HTML report of every chain next to the file and opens it. A file with several models (an NMR ensemble) is analysed model by model.
HOME_DIR="$HOME/.ifacekit"; VENV="$HOME_DIR/venv"
pause() { echo; read -n 1 -s -r -p "Press any key to close this window ..."; echo; }
if [ ! -x "$VENV/bin/ifacekit" ]; then
  echo "ifacekit is not installed yet. Double-click launchers/mac/ifacekit-portal.command once, close it, then try again."
  pause; exit 1
fi
F="$1"
if [ -z "$F" ]; then read -r -p "Drag the structure file into this window and press Enter: " F; fi
F=$(printf '%s' "$F" | sed -e "s/^[[:space:]]*//" -e "s/[[:space:]]*$//" -e "s/^'//" -e "s/'$//" -e 's/^"//' -e 's/"$//' -e 's/\\ / /g')
[ -f "$F" ] || { echo "File not found: $F"; pause; exit 1; }
BASE="${F%.*}"; OUT="${BASE}_ifacekit_report.html"; RAMA="${BASE}_ramachandran.html"
NM=$("$VENV/bin/python" -c 'import gemmi,sys;print(len(gemmi.read_structure(sys.argv[1])))' "$F" 2>/dev/null || echo 1)
if [ "$NM" -gt 1 ]; then
  echo "$NM models found: analysing every model."
  "$VENV/bin/ifacekit" "$F" --chains all --models all --html "$OUT" --rama-html "$RAMA" --workers 4 --quiet || { pause; exit 1; }
else
  "$VENV/bin/ifacekit" "$F" --chains all --html "$OUT" --rama-html "$RAMA" --quiet || { pause; exit 1; }
fi
echo "Report: $OUT"
if command -v open >/dev/null 2>&1; then open "$OUT"; elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$OUT" >/dev/null 2>&1; fi
pause
