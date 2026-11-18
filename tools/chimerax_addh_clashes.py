"""Run inside ChimeraX (started by tools/compare_chimerax_hydrogens.py, not by hand): add hydrogens to every prepared model with
`addh hbond true`, find the clashes with the clash finder of ChimeraX, and save each model with its hydrogens.

Reads the folder from the environment variable CX_WORK (ids.txt, <id>_clean.cif); writes <id>_cx.cif and clashes_cx.json there.
"""
import json
import os

from chimerax.clashes import clashes as cc
from chimerax.core.commands import run

work = os.environ["CX_WORK"]
ids = open(os.path.join(work, "ids.txt")).read().split()
out = {}


def label(x):
    r = x.residue
    return [r.chain_id, str(r.number) + (r.insertion_code or ""), r.name, x.name, x.element.name]


for cid in ids:
    src = os.path.join(work, cid + "_clean.cif")
    try:
        m = run(session, f"open {src}")[0]                      # noqa: F821  (session is provided by ChimeraX)
        run(session, "addh hbond true")                          # noqa: F821
        res = cc.find_clashes(session, m.atoms, attr_name="cxclash", clash_threshold=0.4, hbond_allowance=0.6,   # noqa: F821
                              bond_separation=4, intra_res=True, intra_mol=True, restrict="both")
        rows, seen = [], set()
        for a, d in res.items():
            for b, ov in d.items():
                key = tuple(sorted([id(a), id(b)]))
                if key not in seen:
                    seen.add(key)
                    rows.append([label(a), label(b), float(ov)])
        out[cid] = rows
        run(session, f"save {os.path.join(work, cid + '_cx.cif')} format mmcif models #1")   # noqa: F821
    except Exception as e:                                       # one odd entry must not stop the batch
        out[cid] = {"error": str(e)}
    try:
        run(session, "close")                                    # noqa: F821
    except Exception:
        pass
json.dump(out, open(os.path.join(work, "clashes_cx.json"), "w"))
