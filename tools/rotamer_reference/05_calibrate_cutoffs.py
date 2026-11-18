"""Step 5: calibrate the rotamer cut-offs (and the smoothing) against MolProbity's rotalyze, with a clean split.

Input: residues.pkl, rows (case, table index, our chi angles, rotalyze score, rotalyze category 0/1/2, rotalyze chi angles, rotamer name) made by
extracting the chi angles of ifacekit and reading the output of `phenix.rotalyze m.cif outliers_only=False` for many structures.
Split: CALIBRATION cases (released before the cut-off date, not in the reference used for scoring them) and TEST cases (released on or after the
cut-off date). For calibration the reference grids are built from the entries released before the cut-off WITHOUT the calibration entries; for
the test the grids are built from all entries before the cut-off. So the test entries are never part of anything that was fitted.
For every residue type and every smoothing factor the two cut-offs (the density levels enclosing a fraction FAV and ALW of the reference
residues) are searched: ALW for the best F1 of the outlier class against rotalyze, then FAV for the best three-class agreement. Older and more recent
calibration entries (before / on or after the cut-off date) get the same weight: the score of a setting is the mean of its score on the two groups, because
older entries have more rotamer outliers (for Asn, Gln and His also unflipped groups) and a calibration on a mixture would favour the larger group.
A type with fewer than MIN_OUT rotalyze outliers in the calibration set keeps the cut-offs given with the baseline file.
A second test (TEST B) on older entries can be added: those entries are left out of the grids that score them.
Usage: python 05_calibrate_cutoffs.py residues.pkl chi.npz release_dates.tsv cal_cases.txt test_cases.txt baseline_enclosure.json OUT.json [CUTOFF=2024-01-01] [test_b_cases.txt]
Writes OUT.json = {"scale": s, "enclosure": {"arg": [fav, alw], ...}, "report": ...} and prints the tables.
"""
import json
import math
import pickle
import sys

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

from ifacekit import rotamer

MIN_OUT = 40
SCALES = (1.0, 1.5, 2.0, 3.0)
ALW_GRID = (0.90, 0.95, 0.96, 0.97, 0.975, 0.98, 0.985, 0.99, 0.992, 0.994, 0.996, 0.998, 0.999, 0.9995, 0.9999)
FAV_GRID = (0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99)


def build(z, mask, ti, scale):
    t = rotamer.TABLES[ti]
    lay = rotamer.LAYOUT[t]
    m = (z["table"] == ti) & mask
    data = z["chi"][m][:, :len(lay)].astype(np.float64)
    h, _ = np.histogramdd(data, bins=[nb for _, nb in lay], range=[(0, rng) for rng, _ in lay])
    sig = [scale * rotamer.SIGMA_BY_DIM[len(lay)] / (rng / nb) for rng, nb in lay]
    d = gaussian_filter(h, sig, mode="wrap")
    return d / d.max()


def thresholds(d, fav, alw):
    order = np.sort(d.ravel())[::-1]
    cum = np.cumsum(order) / order.sum()
    pick = lambda e: order[min(np.searchsorted(cum, e), len(order) - 1)]
    return pick(fav), pick(alw)


def density(d, ti, chi):
    lay = rotamer.LAYOUT[rotamer.TABLES[ti]]
    coords = [chi[:, k] / (rng / nb) - 0.5 for k, (rng, nb) in enumerate(lay)]
    return map_coordinates(d, coords, order=1, mode="grid-wrap")


def classify(dens, tf, ta):
    return np.where(dens >= tf, 0, np.where(dens >= ta, 1, 2))


def f1(pred_out, true_out):
    tp = int((pred_out & true_out).sum())
    p = tp / max(int(pred_out.sum()), 1)
    r = tp / max(int(true_out.sum()), 1)
    return 0.0 if tp == 0 else 2 * p * r / (p + r), p, r


def residues_of(rows, names):
    sel = [r for r in rows if r[0] in names]
    by = {}
    for ti in range(len(rotamer.TABLES)):
        rs = [r for r in sel if r[1] == ti]
        if not rs:
            continue
        chi = np.array([r[2] for r in rs], dtype=np.float64)
        by[ti] = (chi, np.array([r[4] for r in rs]))
    return by


def main():
    rows = pickle.load(open(sys.argv[1], "rb"))
    z = np.load(sys.argv[2])
    dates = dict(l.split() for l in open(sys.argv[3]))
    cal = set(open(sys.argv[4]).read().split())
    test = set(open(sys.argv[5]).read().split())
    base = json.load(open(sys.argv[6]))
    cutoff = sys.argv[8] if len(sys.argv) > 8 else "2024-01-01"
    ids = z["entry_ids"]
    entry_date = np.array([dates.get(i, "0") for i in ids])
    before = (entry_date < cutoff)[z["entry"]]
    in_cal = np.array([i in cal for i in ids])[z["entry"]]
    cal_old = residues_of(rows, {n for n in cal if dates.get(n, "0") < cutoff})
    cal_new = residues_of(rows, {n for n in cal if dates.get(n, "0") >= cutoff})
    cal_res, test_res = residues_of(rows, cal), residues_of(rows, test)
    print(f"calibration residues {sum(len(v[1]) for v in cal_res.values())} (rotalyze outliers {sum(int((v[1] == 2).sum()) for v in cal_res.values())}), "
          f"test residues {sum(len(v[1]) for v in test_res.values())} (outliers {sum(int((v[1] == 2).sum()) for v in test_res.values())})")
    best_scale, best_total, results = None, -1, {}
    for scale in SCALES:
        total, per = 0.0, {}
        for ti, (chi, lab) in cal_res.items():
            t = rotamer.TABLES[ti]
            d = build(z, before & ~in_cal, ti, scale)
            groups = []
            for g in (cal_old, cal_new):
                if ti in g and (g[ti][1] == 2).sum() >= 10:
                    groups.append((density(d, ti, g[ti][0]), g[ti][1]))
            if len(groups) < 2:                                   # too few outliers in one group: use all residues
                groups = [(density(d, ti, chi), lab)]
            n_out = int((lab == 2).sum())
            def score_out(tf, ta):
                return float(np.mean([f1(classify(dn, tf, ta) == 2, lb == 2)[0] for dn, lb in groups]))
            def score_acc(tf, ta):
                return float(np.mean([(classify(dn, tf, ta) == lb).mean() for dn, lb in groups]))
            tf0, ta0 = thresholds(d, *base[t])
            b = score_out(tf0, ta0)
            if n_out < MIN_OUT:
                per[t] = (base[t], b, b, n_out, False)
                total += n_out * b
                continue
            bestA, bestF = base[t][1], -1
            for alw in ALW_GRID:
                tf, ta = thresholds(d, 0.9, alw)
                sc = score_out(tf, ta)
                if sc > bestF + 1e-9 or (abs(sc - bestF) <= 1e-9 and abs(alw - base[t][1]) < abs(bestA - base[t][1])):
                    bestF, bestA = sc, alw
            bestV, bestAcc = base[t][0], -1
            for fav in FAV_GRID:
                if fav >= bestA:
                    continue
                tf, ta = thresholds(d, fav, bestA)
                acc = score_acc(tf, ta)
                if acc > bestAcc:
                    bestAcc, bestV = acc, fav
            per[t] = ([bestV, bestA], b, bestF, n_out, True)
            total += n_out * bestF
        results[scale] = per
        print(f"smoothing factor {scale}: sum of (rotalyze outliers x F1) = {total:.1f}")
        if total > best_total:
            best_total, best_scale = total, scale
    per = results[best_scale]
    print(f"\nchosen smoothing factor: {best_scale}")
    print(f"{'type':8s} {'rotalyze outliers':>17s} {'F1 baseline':>12s} {'F1 calibrated':>14s}  fav / alw enclosure")
    enc = {}
    for t in rotamer.TABLES:
        if t in per:
            e, b, c, n, used = per[t]
            enc[t] = [round(float(e[0]), 4), round(float(e[1]), 4)]
            print(f"{t:8s} {n:17d} {b:12.2f} {c:14.2f}  {enc[t][0]} / {enc[t][1]}" + ("" if used else "   (too few outliers, baseline kept)"))
        else:
            enc[t] = base[t]
    def evaluate(title, res, grid_mask):
        print("\n" + title)
        tot = {"base": [0, 0, 0, 0, 0], "cal": [0, 0, 0, 0, 0]}
        per_t = {}
        for ti, (chi, lab) in res.items():
            t = rotamer.TABLES[ti]
            for tag, scale, e in (("base", 1.5, base[t]), ("cal", best_scale, enc[t])):
                d = build(z, grid_mask, ti, scale)
                tf, ta = thresholds(d, *e)
                pred = classify(density(d, ti, chi), tf, ta)
                po, to = pred == 2, lab == 2
                c = tot[tag]
                c[0] += int((po & to).sum()); c[1] += int(po.sum()); c[2] += int(to.sum()); c[3] += int((pred == lab).sum()); c[4] += len(lab)
                per_t.setdefault(t, {})[tag] = (int((po & to).sum()), int(po.sum()), int(to.sum()))
        for tag, name in (("base", "before calibration (cut-offs from the Top8000 enclosure, smoothing 1.5)"), ("cal", f"after calibration (smoothing {best_scale})")):
            a_, b_, c_, acc, n = tot[tag]
            print(f"  {name}: rotalyze outliers {c_}, ifacekit outliers {b_}, shared {a_}: found {100 * a_ / c_:.1f}%, confirmed {100 * a_ / max(b_, 1):.1f}%, "
                  f"three-class agreement {100 * acc / n:.1f}% of {n} residues")
        print(f"  {'type':8s} {'rotalyze':>9s} {'before: ifk / shared':>22s} {'after: ifk / shared':>22s}")
        for t in sorted(per_t, key=lambda t: -per_t[t]['base'][2]):
            b_, c_ = per_t[t]["base"], per_t[t]["cal"]
            print(f"  {t:8s} {b_[2]:9d} {b_[1]:12d} / {b_[0]:<7d} {c_[1]:12d} / {c_[0]:<7d}")

    evaluate("TEST A (independent: entries released on or after the cut-off; grids built from entries before it)", test_res, before)
    if len(sys.argv) > 9:
        test_b = set(open(sys.argv[9]).read().split())
        in_b = np.array([i in test_b for i in ids])[z["entry"]]
        evaluate("TEST B (older entries held out: left out of the grids that score them)", residues_of(rows, test_b), before & ~in_b)
    json.dump({"scale": best_scale, "enclosure": enc, "cutoff": cutoff, "min_outliers": MIN_OUT}, open(sys.argv[7], "w"), indent=1)


if __name__ == "__main__":
    main()
