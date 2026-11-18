"""Compare ifacekit with MolProbity (phenix.molprobity) on many structures (development only).

Step 1, outside this script: run `phenix.molprobity m.cif` in one folder per case (the folder holds m.cif, molprobity_coot.py and
molprobity.out; the large molprobity_probe.txt may be deleted).
Step 2: python tools/compare_molprobity_batch.py CASES_DIR OUT_DIR [--cases cases.txt] [--dates release_dates.tsv --cutoff 2024-01-01
        --rama-pre grid.npz --rota-pre grid.npz] [--procs 16]
Per case: ifacekit.analyze (clashes), ifacekit.rama and ifacekit.rotamer against MolProbity's severe clashes ('probe' list, overlap of
0.4 A or more, with hydrogens added by Reduce), Ramachandran outliers and rotamer outliers. With --rama-pre / --rota-pre the
Ramachandran and rotamer references built only from entries released before the cutoff are used as well, so that entries released
after the cutoff give a test that is independent of the reference.
Writes OUT_DIR/cases.json (one record per case) and OUT_DIR/summary.txt.
"""
import argparse
import ast
import json
import multiprocessing as mp
import os
import re
import sys
import time

import numpy as np
from scipy.stats import spearmanr

import ifacekit
from ifacekit import rama, rotamer

CLASH_REASONS = ("steric_overlap", "close_contact")


def read_coot(path):
    d = {}
    for line in open(path):
        m = re.match(r"data\['(\w+)'\] = (.*)$", line.strip())
        if m:
            d[m.group(1)] = ast.literal_eval(m.group(2))
    return d


def read_summary(path):
    t = open(path).read()
    g = lambda pat: float(re.search(pat, t).group(1)) if re.search(pat, t) else None
    return {"rama_out": g(r"Ramachandran outliers\s*=\s*([\d.]+)"), "rama_fav": g(r"favored\s*=\s*([\d.]+)"),
            "rota_out": g(r"Rotamer outliers\s*=\s*([\d.]+)"), "clashscore": g(r"Clashscore\s*=\s*([\d.]+)")}


def probe_atom(s):
    """Fixed-width Probe label: chain(2) resnum(4) icode(1) altloc(1) resname(3) space atom."""
    return (s[0:2].strip(), s[2:6].strip() + s[6].strip()), s[8:11], s[12:].strip()


def is_h(name):
    return name[:1] in ("H", "D")


def res_key(chain, num):
    return (chain.strip(), str(num).strip())


def split_atom(a):
    ch, num, rn, an = a.split(":")[:4]
    return res_key(ch, num), rn, an


def clash_part(cif, coot):
    r = ifacekit.analyze(cif, n_points=240)
    flagged, unf, flagged_clash, unf_clash = set(), set(), set(), set()
    res_any = set()
    for c in r["contacts"]:
        if c["verdict"] not in ("borderline", "unfavorable"):
            continue
        ka, _, _ = split_atom(c["atom1"])
        kb, _, _ = split_atom(c["atom2"])
        key = frozenset([ka, kb])
        flagged.add(key)
        res_any.update([ka, kb])
        clash = c["reason"] in CLASH_REASONS
        if clash:
            flagged_clash.add(key)
        if c["verdict"] == "unfavorable":
            unf.add(key)
            if clash:
                unf_clash.add(key)
    mp_all = []
    for a, b, ov, _ in coot.get("probe", []):
        (ra, na, aa), (rb, nb, ab) = probe_atom(a), probe_atom(b)
        alt = (a[7:8] not in (" ", "")) or (b[7:8] not in (" ", ""))
        wat = a[8:11] == "HOH" or b[8:11] == "HOH"
        try:
            adj = ra[0] == rb[0] and abs(int(re.sub(r"\D", "", ra[1])) - int(re.sub(r"\D", "", rb[1]))) == 1
        except ValueError:
            adj = False
        cat = "same residue" if ra == rb else "water" if wat else "adjacent residue" if adj else "alternate conformer" if alt else "other residue pair"
        el = "-".join(sorted([aa[:1], ab[:1]]))
        mp_all.append({"pair": frozenset([ra, rb]), "ov": abs(ov), "heavy": not (is_h(aa) or is_h(ab)), "same": ra == rb, "res": (ra, rb), "cat": cat, "el": el})
    mp_keys = {x["pair"] for x in mp_all}
    mp_heavy_keys = {x["pair"] for x in mp_all if x["heavy"]}
    heavy = [x for x in mp_all if x["heavy"]]
    hyd = [x for x in mp_all if not x["heavy"]]
    out = {"n_atoms": r["summary"]["n_non_water_atoms"], "ifk_unf": len(unf), "ifk_border": len(flagged) - len(unf),
           "ifk_unf_clash": len(unf_clash), "mp_total": len(mp_all), "mp_heavy": len(heavy), "mp_h": len(hyd),
           "mp_heavy_inter": sum(1 for x in heavy if not x["same"])}
    for tag, th in (("04", 0.4), ("06", 0.6), ("08", 0.8)):
        h = [x for x in heavy if x["ov"] >= th]
        out[f"heavy_n_{tag}"] = len(h)
        out[f"heavy_hit_{tag}"] = sum(1 for x in h if x["pair"] in flagged)
        out[f"heavy_hit_clash_{tag}"] = sum(1 for x in h if x["pair"] in flagged_clash)
    inter = [x for x in heavy if not x["same"]]
    out["heavy_inter_hit"] = sum(1 for x in inter if x["pair"] in flagged)
    cats = {}
    for x in heavy:
        c = cats.setdefault(x["cat"], [0, 0])
        c[0] += 1
        c[1] += 1 if x["pair"] in flagged else 0
    out["heavy_cat"] = cats
    els = {}
    for x in heavy:
        c = els.setdefault(x["el"], [0, 0])
        c[0] += 1
        c[1] += 1 if x["pair"] in flagged else 0
    out["heavy_el"] = els
    out["h_n"] = len(hyd)
    out["h_pair_hit"] = sum(1 for x in hyd if x["pair"] in flagged)
    out["h_res_hit"] = sum(1 for x in hyd if x["res"][0] in res_any or x["res"][1] in res_any)
    out["unf_in_mp"] = sum(1 for k in unf if k in mp_keys)
    out["unf_in_mp_heavy"] = sum(1 for k in unf if k in mp_heavy_keys)
    out["unf_clash_in_mp"] = sum(1 for k in unf_clash if k in mp_keys)
    out["border_in_mp"] = sum(1 for k in (flagged - unf) if k in mp_keys)
    return out


def outlier_part(kind, cif, coot, ref_path):
    mod = rama if kind == "rama" else rotamer
    mod.use_reference(ref_path)
    summ, rows = mod.analyze(cif)
    ifk_names = {res_key(r["residue"].split(":")[0], r["residue"].split(":")[1]): r["residue"].split(":")[2] for r in rows if r["category"] == "outlier"}
    ifk = set(ifk_names)
    ent = coot.get("rama" if kind == "rama" else "rota", [])
    mp_names = {res_key(e[0], e[1]): e[2] for e in ent}
    mp_set = set(mp_names)
    names = {}
    for k in ifk | mp_set:
        nm = ifk_names.get(k) or mp_names.get(k)
        c = names.setdefault(nm, [0, 0, 0])
        c[0] += 1 if (k in ifk and k in mp_set) else 0
        c[1] += 1 if (k in ifk and k not in mp_set) else 0
        c[2] += 1 if (k in mp_set and k not in ifk) else 0
    return {"n": summ["n"], "ifk_out": len(ifk), "mp_out": len(mp_set), "both": len(ifk & mp_set), "by_res": names,
            "ifk_fav_pct": summ["percent"]["favored"], "ifk_out_pct": summ["percent"]["outlier"]}


def one(args):
    d, name, pre = args
    t = time.time()
    cif = os.path.join(d, "m.cif")
    try:
        coot = read_coot(os.path.join(d, "molprobity_coot.py"))
        rec = {"id": name, "mp": read_summary(os.path.join(d, "molprobity.out"))}
        rec["clash"] = clash_part(cif, coot)
        for kind in ("rama", "rota"):
            rec[kind] = {"packaged": outlier_part(kind, cif, coot, None)}
            ref = pre.get(kind)
            if ref:
                rec[kind]["pre"] = outlier_part(kind, cif, coot, ref)
        rec["seconds"] = round(time.time() - t, 1)
        return rec
    except Exception as e:  # keep going, report the failure
        return {"id": name, "error": f"{type(e).__name__}: {e}"}


def pct(a, b):
    return f"{100.0 * a / b:.1f}%" if b else "n/a"


def summarize(recs, dates, cutoff):
    ok = [r for r in recs if "error" not in r]
    post = [r for r in ok if dates and dates.get(r["id"], "0") >= cutoff]
    pre_ = [r for r in ok if dates and dates.get(r["id"], "0") < cutoff]
    L = [f"cases {len(recs)}, analysed {len(ok)}, failed {len(recs) - len(ok)}"
         + (f"; released {cutoff} or later: {len(post)}, earlier: {len(pre_)}" if dates else "")]
    L.append("")
    L.append("CLASHES  (MolProbity severe clashes: overlap 0.4 A or more after Reduce adds hydrogens; ifacekit: heavy atoms only)")
    for tag, rs in (("all cases", ok),) + ((("released after the cutoff", post), ("released before the cutoff", pre_)) if dates else ()):
        c = lambda k: sum(r["clash"][k] for r in rs)
        L.append(f"  [{tag}] structures {len(rs)}")
        L.append(f"    MolProbity severe clashes: {c('mp_total')} (heavy-atom pairs {c('mp_heavy')}, involving hydrogen {c('mp_h')})")
        for tag2, th in (("04", "0.4"), ("06", "0.6"), ("08", "0.8")):
            L.append(f"    heavy-atom clashes with overlap >= {th} A: {c('heavy_n_' + tag2)}; residue pair also flagged by ifacekit (borderline or unfavorable, any reason): "
                     f"{c('heavy_hit_' + tag2)} ({pct(c('heavy_hit_' + tag2), c('heavy_n_' + tag2))}); flagged as a steric clash: {pct(c('heavy_hit_clash_' + tag2), c('heavy_n_' + tag2))}")
        L.append(f"    heavy-atom clashes between different residues: {c('mp_heavy_inter')}; flagged {pct(c('heavy_inter_hit'), c('mp_heavy_inter'))}")
        L.append(f"    clashes involving hydrogen: {c('h_n')}; residue pair flagged {pct(c('h_pair_hit'), c('h_n'))}; one of the two residues flagged anywhere {pct(c('h_res_hit'), c('h_n'))}")
        L.append(f"    ifacekit unfavorable residue pairs: {c('ifk_unf')}; also in MolProbity's severe list {c('unf_in_mp')} ({pct(c('unf_in_mp'), c('ifk_unf'))}), "
                 f"in its heavy-atom list {pct(c('unf_in_mp_heavy'), c('ifk_unf'))}; steric-clash type only: {pct(c('unf_clash_in_mp'), c('ifk_unf_clash'))} of {c('ifk_unf_clash')}")
        L.append(f"    ifacekit borderline residue pairs: {c('ifk_border')}; in MolProbity's list {pct(c('border_in_mp'), c('ifk_border'))}")
    L.append("  heavy-atom MolProbity clashes by kind (flagged by ifacekit / all):")
    cat, el = {}, {}
    for r in ok:
        for k, (n, h) in r["clash"]["heavy_cat"].items():
            c = cat.setdefault(k, [0, 0]); c[0] += n; c[1] += h
        for k, (n, h) in r["clash"]["heavy_el"].items():
            c = el.setdefault(k, [0, 0]); c[0] += n; c[1] += h
    for k, (n, h) in sorted(cat.items(), key=lambda kv: -kv[1][0]):
        L.append(f"    {k:22s} {h:5d} / {n:5d}  ({pct(h, n)})")
    L.append("  by element pair: " + ", ".join(f"{k} {pct(h, n)} of {n}" for k, (n, h) in sorted(el.items(), key=lambda kv: -kv[1][0])[:8]))
    cs = [r["mp"]["clashscore"] for r in ok]
    pk = lambda k: [1000.0 * r["clash"][k] / max(r["clash"]["n_atoms"], 1) for r in ok]
    mh = [r["clash"]["mp_heavy"] for r in ok]
    L.append(f"  per structure (Spearman): ifacekit unfavorable per 1000 atoms vs MolProbity clashscore {spearmanr(pk('ifk_unf'), cs)[0]:.2f}; "
             f"(unfavorable + borderline) per 1000 atoms vs clashscore "
             f"{spearmanr([a + b for a, b in zip(pk('ifk_unf'), pk('ifk_border'))], cs)[0]:.2f}; ifacekit unfavorable count vs MolProbity heavy-atom clash count "
             f"{spearmanr([r['clash']['ifk_unf'] for r in ok], mh)[0]:.2f}")
    L.append(f"  structures with no MolProbity severe clash: {sum(1 for r in ok if r['clash']['mp_total'] == 0)}; "
             f"with no ifacekit unfavorable pair: {sum(1 for r in ok if r['clash']['ifk_unf'] == 0)}; "
             f"median clashscore {np.median(cs):.1f}, median ifacekit unfavorable per 1000 atoms {np.median(pk('ifk_unf')):.2f}")
    for kind, title in (("rama", "RAMACHANDRAN (ramalyze outliers)"), ("rota", "ROTAMERS (rotalyze outliers)")):
        L.append("")
        L.append(title)
        sets = [("packaged reference, all cases", ok, "packaged")]
        if dates:
            sets = [(f"independent: released {cutoff} or later, reference built from earlier entries only", post, "pre"),
                    (f"released {cutoff} or later, packaged reference (the reference contains these entries)", post, "packaged"),
                    (f"released before {cutoff}, packaged reference (in the reference set)", pre_, "packaged")]
        for tag, rs, which in sets:
            rs = [r for r in rs if which in r[kind]]
            if not rs:
                continue
            g = lambda k: sum(r[kind][which][k] for r in rs)
            mpo, ifo, both = g("mp_out"), g("ifk_out"), g("both")
            L.append(f"  [{tag}] structures {len(rs)}, residues scored by ifacekit {g('n')}")
            L.append(f"    MolProbity outliers {mpo}, ifacekit outliers {ifo}, both {both}: found by ifacekit {pct(both, mpo)}, also outliers for MolProbity {pct(both, ifo)}")
            if kind == "rama":
                a = [r["mp"]["rama_fav"] for r in rs]
                b = [r[kind][which]["ifk_fav_pct"] for r in rs]
                L.append(f"    favored % per structure: MolProbity mean {np.mean(a):.1f}, ifacekit mean {np.mean(b):.1f}; mean absolute difference {np.mean(np.abs(np.array(a) - np.array(b))):.2f}; Spearman {spearmanr(a, b)[0]:.2f}")
            if kind == "rota" and which in ("pre", "packaged") and tag.startswith("independent"):
                bt = {}
                for r in rs:
                    for nm, (b_, i_, m_) in r[kind][which]["by_res"].items():
                        c = bt.setdefault(nm, [0, 0, 0]); c[0] += b_; c[1] += i_; c[2] += m_
                L.append("    by residue type (both / only ifacekit / only MolProbity):")
                for nm, (b_, i_, m_) in sorted(bt.items(), key=lambda kv: -(kv[1][0] + kv[1][1] + kv[1][2])):
                    L.append(f"      {nm}: {b_} / {i_} / {m_}   (found {pct(b_, b_ + m_)}, confirmed {pct(b_, b_ + i_)})")
            a = [r["mp"][f"{kind}_out"] for r in rs]
            b = [r[kind][which]["ifk_out_pct"] for r in rs]
            L.append(f"    outlier % per structure: MolProbity mean {np.mean(a):.2f}, ifacekit mean {np.mean(b):.2f}; Spearman {spearmanr(a, b)[0]:.2f}")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cases_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--cases")
    ap.add_argument("--dates")
    ap.add_argument("--cutoff", default="2024-01-01")
    ap.add_argument("--rama-pre")
    ap.add_argument("--rota-pre")
    ap.add_argument("--procs", type=int, default=16)
    a = ap.parse_args()
    names = [l.strip() for l in open(a.cases)] if a.cases else sorted(os.listdir(a.cases_dir))
    names = [n for n in names if n and os.path.exists(os.path.join(a.cases_dir, n, "molprobity_coot.py"))]
    dates = dict(l.split() for l in open(a.dates)) if a.dates else None
    pre = {"rama": a.rama_pre, "rota": a.rota_pre}
    os.makedirs(a.out_dir, exist_ok=True)
    with mp.Pool(a.procs) as pool:
        recs = pool.map(one, [(os.path.join(a.cases_dir, n), n, pre) for n in names], chunksize=2)
    json.dump(recs, open(os.path.join(a.out_dir, "cases.json"), "w"))
    text = summarize(recs, dates, a.cutoff)
    open(os.path.join(a.out_dir, "summary.txt"), "w").write(text + "\n")
    for r in recs:
        if "error" in r:
            print("failed:", r["id"], r["error"], file=sys.stderr)
    print(text)


if __name__ == "__main__":
    main()
