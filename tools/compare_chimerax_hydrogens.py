#!/usr/bin/env python3
"""Compare hydrogens added by ChimeraX with ifacekit's own, by the clashes they give against MolProbity (a validation experiment; ChimeraX is
not used by ifacekit).

Input: the folder of a MolProbity comparison (tools/compare_molprobity_batch.py; one folder per case with m.cif and molprobity_coot.py), the list of
cases, and a work folder. For every case the model is cleaned (no hydrogens, first conformer only), ChimeraX adds hydrogens (`addh hbond true`) and finds
clashes with its own finder (tools/chimerax_addh_clashes.py), and four ways of finding clashes are compared with MolProbity's pairs (first conformer,
residue-pair level, as in tools/compare_hydrogen_clashes.py):
  ifacekit            hydrogens and clash rules of ifacekit
  chimerax_as_placed  ChimeraX hydrogens, ifacekit clash rules (ChimeraX puts C-H at 1.09 A, the nuclear distance)
  chimerax_shortened  ChimeraX hydrogens moved to the shortened X-ray X-H distances, ifacekit clash rules
  chimerax_clashes    the `clashes` command of ChimeraX (overlap 0.4 A, hydrogen-bond allowance 0.6 A, 4 bonds)
Usage: python tools/compare_chimerax_hydrogens.py CASES.txt MP_DIR WORK_DIR OUT.json [--chimerax /usr/bin/chimerax] [--jobs 4] [--skip-chimerax]
Writes OUT.json = {"cases": [...], "summary": {...}}; the work folder keeps the cleaned models, the ChimeraX models and the ChimeraX clash lists.
"""
import argparse
import collections
import json
import multiprocessing as mp
import os
import subprocess
import sys

import gemmi
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_hydrogen_clashes as chc  # noqa: E402

from ifacekit import hydrogens as hy  # noqa: E402
from ifacekit import model  # noqa: E402

VARIANTS = (("ifacekit", "ifacekit's hydrogens and clash rules"),
            ("chimerax_as_placed", "ChimeraX hydrogens as placed, ifacekit clash rules"),
            ("chimerax_shortened", "ChimeraX hydrogens moved to X-ray X-H distances, ifacekit clash rules"),
            ("chimerax_clashes", "ChimeraX hydrogens and the clash finder of ChimeraX"))


def prepare(ids, work):
    os.makedirs(work, exist_ok=True)
    for cid in ids:
        st = gemmi.read_structure(os.path.join(ARGS.mp_dir, cid, "m.cif"))
        st.remove_hydrogens()
        st.remove_alternative_conformations()
        st.setup_entities()
        st.make_mmcif_document().write_file(os.path.join(work, cid + "_clean.cif"))
    open(os.path.join(work, "ids.txt"), "w").write("\n".join(ids))


def env_from_h(m, h_path, shorten):
    """Hydrogens read from a ChimeraX model, attached to the nearest heavy atom of the same residue."""
    st = gemmi.read_structure(h_path)
    amap = hy._atom_map(m)
    key2res = {(r["chain"], r["resnum"], r["icode"]): k for k, r in enumerate(m.residues)}
    H = hy.Hydrogens()
    for ch in st[0]:
        for r in ch:
            ri = key2res.get((ch.name, r.seqid.num, r.seqid.icode.strip()))
            if ri is None:
                continue
            for a in r:
                if not a.element.is_hydrogen:
                    continue
                p = np.array(a.pos.tolist())
                dmin, par = min((np.linalg.norm(m.xyz[k] - p), k) for k in amap[ri].values())
                if dmin > 1.3:
                    continue
                pel, nm, rn = m.element[par], m.atom_name[par], m.residues[ri]["resname"]
                polar = pel in ("N", "O", "S")
                if shorten and pel in hy.XH:
                    v = p - m.xyz[par]
                    p = m.xyz[par] + v / np.linalg.norm(v) * hy.XH[pel]
                rad = hy.H_POLAR if (polar or nm in hy.AROMATIC_C.get(rn, ())) else hy.H_ALI
                H.add(p, par, a.name.strip(), ri, polar, rad, bool(m.charge[par] > 0 and polar))
    H.finish()
    H.active = np.ones(H.n, bool)
    return H, hy.Env(m, H)


def one(cid):
    try:
        work, cx = ARGS.work, CX
        mpall, _, _ = chc.mp_pairs(os.path.join(ARGS.mp_dir, cid))
        m = model.load(os.path.join(work, cid + "_clean.cif"))

        def keys(H, env, cl):
            def k(x):
                r = m.residues[m.res_idx[x]] if x < env.n else m.residues[H.res[x - env.n]]
                return (r["chain"], str(r["resnum"]) + r["icode"])
            return {frozenset([k(x["i"]), k(x["j"])]) for x in cl}

        res = {"id": cid, "mp": len(mpall)}
        H, env, _ = hy.add_hydrogens(m)
        got = {"ifacekit": keys(H, env, hy.clashes(m, H, env, heavy=True))}
        res["n_h_ifacekit"] = int(H.active.sum())
        for name, shorten in (("chimerax_as_placed", False), ("chimerax_shortened", True)):
            H, env = env_from_h(m, os.path.join(work, cid + "_cx.cif"), shorten)
            got[name] = keys(H, env, hy.clashes(m, H, env, heavy=True))
            if not shorten:
                res["n_h_chimerax"] = int(H.n)
        got["chimerax_clashes"] = {frozenset([(a[0], a[1]), (b[0], b[1])]) for a, b, _ in cx[cid]}
        for name, pairs in got.items():
            res[name] = [len(pairs), len(pairs & mpall)]
        return res
    except Exception as e:
        return {"id": cid, "error": str(e)}


def init(args, cx):
    global ARGS, CX
    ARGS, CX = args, cx


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("cases", "mp_dir", "work", "out"):
        ap.add_argument(a)
    ap.add_argument("--chimerax", default="chimerax")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--skip-chimerax", action="store_true", help="use the ChimeraX files already in the work folder")
    args = ap.parse_args()
    global ARGS
    ARGS = args
    ids = [c for c in open(args.cases).read().split() if os.path.exists(os.path.join(args.mp_dir, c, "molprobity_coot.py"))]
    version = ""
    if not args.skip_chimerax:
        prepare(ids, args.work)
        env = dict(os.environ, CX_WORK=os.path.abspath(args.work))
        subprocess.run([args.chimerax, "--nogui", "--exit", "--script", os.path.join(HERE, "chimerax_addh_clashes.py")], env=env, check=True)
    try:
        version = subprocess.run([args.chimerax, "--version"], capture_output=True, text=True).stdout.split("version:")[-1].split("\n")[0].strip()
    except OSError:
        pass
    cx = json.load(open(os.path.join(args.work, "clashes_cx.json")))
    ids = [c for c in ids if isinstance(cx.get(c), list) and os.path.exists(os.path.join(args.work, c + "_cx.cif"))]
    with mp.Pool(args.jobs, initializer=init, initargs=(args, cx)) as pool:
        rows = pool.map(one, ids, chunksize=2)
    ok = [r for r in rows if "error" not in r]
    total_mp = sum(r["mp"] for r in ok)
    summary = {"cases": len(ok), "mp_pairs": total_mp, "chimerax_version": version,
               "hydrogens_chimerax": sum(r["n_h_chimerax"] for r in ok), "hydrogens_ifacekit": sum(r["n_h_ifacekit"] for r in ok), "variants": {}}
    for name, label in VARIANTS:
        n, sh = sum(r[name][0] for r in ok), sum(r[name][1] for r in ok)
        summary["variants"][name] = {"label": label, "pairs": n, "shared": sh, "found_pct": round(100.0 * sh / max(total_mp, 1), 1),
                                     "confirmed_pct": round(100.0 * sh / max(n, 1), 1)}
    json.dump({"cases": rows, "summary": summary}, open(args.out, "w"), indent=1)
    print(json.dumps(summary, indent=1))
    if len(ok) != len(rows):
        print("errors:", [r for r in rows if "error" in r][:3])


if __name__ == "__main__":
    main()
