#!/usr/bin/env python3
"""Compare the clashes ifacekit finds with hydrogens added against the clash list of MolProbity (probe), case by case.

Input: the folder of a MolProbity comparison made with tools/compare_molprobity_batch.py (it holds one folder per PDB entry with the
model `m.cif` and `molprobity_coot.py`, the lists of phenix.molprobity), and the list of cases.
For every case ifacekit places hydrogens (ifacekit.hydrogens), takes all clashes (hydrogen and heavy-atom pairs) and compares them with
the pairs of the probe list at the level of residue pairs (two atoms of the same two residues are one pair).
Alternate conformers are left out: the hydrogen check works on the first conformer of every atom, so the pairs of the probe list that involve
an atom of any other conformer are not counted (their number is reported as `mp_alt_excluded`).
Usage: python tools/compare_hydrogen_clashes.py CASES.txt MP_DIR OUT.json [--jobs 4] [--dev-cases N]
  --dev-cases N: the first N cases are the ones the settings were tuned on; the summary reports them apart from the others.
Writes OUT.json = {"cases": [...], "summary": {...}} and prints the summary.
"""
import argparse
import json
import multiprocessing as mp
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_molprobity_batch as cm  # noqa: E402

from ifacekit import hydrogens as hy  # noqa: E402
from ifacekit import model  # noqa: E402


def first_conformers(path):
    """(chain, residue number + insertion code) -> letter of the first alternate conformer of that residue."""
    import gemmi
    st = gemmi.read_structure(path)
    fa = {}
    for ch in st[0]:
        for r in ch:
            a = next((x.altloc for x in r if x.altloc != "\0"), None)
            if a:
                fa[(ch.name, str(r.seqid.num) + r.seqid.icode.strip())] = a
    return fa


def mp_pairs(case_dir):
    """The clash pairs of MolProbity (residue pairs) for one case, first conformer only.

    Returns (all pairs, pairs with a hydrogen, pairs left out because an atom belongs to a conformer other than the first)."""
    coot = cm.read_coot(os.path.join(case_dir, "molprobity_coot.py"))
    fa = first_conformers(os.path.join(case_dir, "m.cif"))
    mpall, mph, excluded = set(), set(), set()
    for a, b, o, _ in coot.get("probe", []):
        ra, rb = cm.probe_atom(a), cm.probe_atom(b)
        key = frozenset([(ra[0][0], ra[0][1]), (rb[0][0], rb[0][1])])
        if any(lab[7] != " " and fa.get(r[0]) != lab[7] for lab, r in ((a, ra), (b, rb))):
            excluded.add(key)
            continue
        mpall.add(key)
        if cm.is_h(ra[2]) or cm.is_h(rb[2]):
            mph.add(key)
    return mpall, mph, excluded


def one(args):
    cid, mp_dir = args
    d = os.path.join(mp_dir, cid)
    if not os.path.exists(os.path.join(d, "molprobity_coot.py")):
        return {"id": cid, "failed": True}
    try:
        m = model.load(os.path.join(d, "m.cif"))
        t = time.time()
        H, env, groups = hy.add_hydrogens(m)
        cl = hy.clashes(m, H, env, heavy=True)
        sec = time.time() - t
    except Exception as e:                       # one odd entry must not stop the batch
        return {"id": cid, "failed": True, "error": str(e)}
    n = env.n

    def rkey(k):
        if k < n:
            r = m.residues[m.res_idx[k]]
        else:
            r = m.residues[H.res[k - n]]
        return (r["chain"], str(r["resnum"]) + r["icode"])

    mine, mine_h = set(), set()
    for x in cl:
        a, b = rkey(x["i"]), rkey(x["j"])
        key = frozenset([a, b])
        mine.add(key)
        if x["i"] >= n or x["j"] >= n:
            mine_h.add(key)
    mpall, mph, excluded = mp_pairs(d)
    n_flip = sum(1 for g in groups if (g.kind == "flip" and g.state == 1) or (g.kind == "his" and g.state[0]))
    return {"id": cid, "n_atoms": int(n), "n_h": int(H.active.sum()), "mp": len(mpall), "mp_h": len(mph), "mp_heavy": len(mpall - mph),
            "ours": len(mine), "ours_h": len(mine_h), "shared": len(mine & mpall), "shared_h": len(mine & mph),
            "shared_heavy": len(mine & (mpall - mph)), "flips": n_flip, "mp_alt_excluded": len(excluded - mpall), "seconds": round(sec, 1)}


def summarize(rows):
    ok = [r for r in rows if not r.get("failed")]
    if not ok:
        return {}
    s = {k: sum(r[k] for r in ok) for k in ("mp", "mp_h", "mp_heavy", "ours", "ours_h", "shared", "shared_h", "shared_heavy", "flips", "n_atoms", "mp_alt_excluded")}
    s["cases"] = len(ok)
    s["found_pct"] = round(100.0 * s["shared"] / max(s["mp"], 1), 1)
    s["found_h_pct"] = round(100.0 * s["shared_h"] / max(s["mp_h"], 1), 1)
    s["found_heavy_pct"] = round(100.0 * s["shared_heavy"] / max(s["mp_heavy"], 1), 1)
    s["confirmed_pct"] = round(100.0 * s["shared"] / max(s["ours"], 1), 1)
    xs, ys = [r["mp"] for r in ok], [r["ours"] for r in ok]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    var = (sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys)) ** 0.5
    s["correlation"] = round(cov / var, 2) if var else None
    s["median_seconds"] = sorted(r["seconds"] for r in ok)[len(ok) // 2]
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cases")
    ap.add_argument("mp_dir")
    ap.add_argument("out")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--dev-cases", type=int, default=0)
    a = ap.parse_args()
    ids = open(a.cases).read().split()
    with mp.Pool(a.jobs) as pool:
        rows = pool.map(one, [(c, a.mp_dir) for c in ids], chunksize=2)
    dev, test = rows[:a.dev_cases], rows[a.dev_cases:]
    out = {"cases": rows, "summary": {"all": summarize(rows), "tuning": summarize(dev), "independent": summarize(test) if a.dev_cases else {}}}
    json.dump(out, open(a.out, "w"), indent=1)
    for k, v in out["summary"].items():
        print(k, v)
    print("failed (no MolProbity result):", [r["id"] for r in rows if r.get("failed")])


if __name__ == "__main__":
    main()
