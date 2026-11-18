#!/usr/bin/env python3
"""Compare ifacekit's interface areas with PISA (CCP4) on many structures.

Input: the folder of a MolProbity comparison (one folder per case with m.cif; tools/compare_molprobity_batch.py) and the list of cases.
For every case PISA analyses the model (`pisa NAME -analyse FILE`, `pisa NAME -xml interfaces`); every interface between two chains of the model
itself (both molecules the identity copy, no symmetry mate) is compared with `ifacekit --pair A B`: interface area, number of interface residues,
hydrogen bonds and salt bridges. Interface area is the mean of the area buried on the two sides, in both programs.
Usage: python tools/compare_pisa_batch.py CASES.txt MP_DIR WORK_DIR OUT.json [--jobs 4] [--min-area 100] [--max-pairs 4]
       (needs the CCP4 environment: set CCP4_SETUP to the path of ccp4.setup-sh of your CCP4 installation)
Writes OUT.json = {"pairs": [...], "summary": {...}}. A pair counts only when PISA gives an area of at least --min-area; at most --max-pairs
(the largest) per structure are used, so that one large complex does not dominate.
"""
import argparse
import json
import multiprocessing as mp
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_pisa as cp  # noqa: E402

SETUP = os.environ.get("CCP4_SETUP", "ccp4.setup-sh")


def run_pisa(cid, src, work):
    xml = os.path.join(work, cid + "_pisa_interfaces.xml")
    if not (os.path.exists(xml) and os.path.getsize(xml) > 100):
        shell = f'. "{SETUP}" >/dev/null 2>&1; cd "{work}"; pisa {cid} -analyse "{src}" >/dev/null 2>&1; pisa {cid} -xml interfaces > "{xml}" 2>/dev/null'
        subprocess.run(["bash", "-c", shell], check=False, timeout=900)
    return xml


def one(args):
    cid, mp_dir, work, min_area, max_pairs = args
    from ifacekit.core import analyze
    src = os.path.join(mp_dir, cid, "m.cif")
    out = []
    try:
        xml = run_pisa(cid, src, work)
        pairs = [p for p in cp.pisa_pairs(xml) if p["area"] >= min_area and "[" not in p["a"] and "[" not in p["b"]]    # PISA lists ligands as molecules too: chain pairs only
    except Exception as e:
        return [{"id": cid, "error": "PISA: " + str(e)[:120]}]
    pairs.sort(key=lambda p: -p["area"])
    for p in pairs[:max_pairs]:
        try:
            r = analyze(src, pair=(p["a"], p["b"]), n_points=960)
            if not r["interfaces"]:
                out.append({"id": cid, "a": p["a"], "b": p["b"], "pisa_area": p["area"], "ifk_area": 0.0})
                continue
            it = r["interfaces"][0]
            hb = it.get("hbonds") or {}
            out.append({"id": cid, "a": p["a"], "b": p["b"], "pisa_area": round(p["area"], 1), "ifk_area": it["area"],
                        "pisa_nres": p["nres"], "ifk_nres": it["n_residues_a"] + it["n_residues_b"],
                        "pisa_hb35": p["hb35"], "ifk_hb": hb.get("direct"), "pisa_sb": p["sb"],
                        "ifk_sb": it["contacts"].get("salt_bridge", 0)})
        except Exception as e:
            out.append({"id": cid, "a": p["a"], "b": p["b"], "error": "ifacekit: " + str(e)[:120]})
    return out


def stats(rows, key_p, key_i, rel=True):
    import statistics
    x = [r[key_p] for r in rows]
    y = [r[key_i] for r in rows]
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y))
    r = cov / ((sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)) ** 0.5)
    slope = cov / sum((a - mx) ** 2 for a in x)
    d = [100.0 * (b - a) / a for a, b in zip(x, y) if a]
    ad = sorted(abs(v) for v in d)
    return {"n": n, "pearson_r": round(r, 4), "slope": round(slope, 3), "median_diff_pct": round(statistics.median(d), 2),
            "median_abs_diff_pct": round(statistics.median(ad), 2), "within_5pct": round(100.0 * sum(v <= 5 for v in ad) / len(ad), 1),
            "within_10pct": round(100.0 * sum(v <= 10 for v in ad) / len(ad), 1), "p90_abs_diff_pct": round(ad[int(0.9 * (len(ad) - 1))], 1)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("cases", "mp_dir", "work", "out"):
        ap.add_argument(a)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--min-area", type=float, default=100.0)
    ap.add_argument("--max-pairs", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    ids = [c for c in open(a.cases).read().split() if os.path.exists(os.path.join(a.mp_dir, c, "m.cif"))]
    with mp.Pool(a.jobs) as pool:
        res = pool.map(one, [(c, a.mp_dir, a.work, a.min_area, a.max_pairs) for c in ids], chunksize=1)
    rows = [r for lst in res for r in lst]
    ok = [r for r in rows if "error" not in r and r.get("pisa_area") and r.get("ifk_area") is not None]
    summary = {"structures_with_a_compared_pair": len({r["id"] for r in ok}), "pairs": len(ok), "errors": len(rows) - len(ok),
               "area": stats(ok, "pisa_area", "ifk_area"),
               "area_ge_300": stats([r for r in ok if r["pisa_area"] >= 300], "pisa_area", "ifk_area") if any(r["pisa_area"] >= 300 for r in ok) else None,
               "area_ge_1000": stats([r for r in ok if r["pisa_area"] >= 1000], "pisa_area", "ifk_area") if any(r["pisa_area"] >= 1000 for r in ok) else None,
               "residues": stats([r for r in ok if r.get("pisa_nres")], "pisa_nres", "ifk_nres"),
               "min_area": a.min_area, "max_pairs": a.max_pairs}
    json.dump({"pairs": rows, "summary": summary}, open(a.out, "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
