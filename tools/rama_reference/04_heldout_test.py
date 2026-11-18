"""Held-out test of the Ramachandran reference.

Grids are built (same method and sigma as the packaged reference) from entries released BEFORE the cutoff date only.
The residues of entries released ON OR AFTER the cutoff are then classified. If the regions generalise, the held-out
residues should fall into favored / allowed / outlier at about the rates used to define them (98% / ~2% / <0.05%).
Optional: a Top8000 grid file (dev only, from tools/build_rama_data.py) is scored on the same residues.
Usage: python 04_heldout_test.py phipsi.npz release_dates.tsv CUTOFF [SIGMA=4] [TOP8000.npz]
"""
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

CLASSES = ("general", "glycine", "cispro", "transpro", "prepro", "ileval")
ALLOWED8000 = {"general": 0.0005, "glycine": 0.001, "cispro": 0.002, "transpro": 0.001, "prepro": 0.001, "ileval": 0.001}
z = np.load(sys.argv[1])
dates = dict(l.split() for l in open(sys.argv[2]))
cutoff = sys.argv[3]
sigma = float(sys.argv[4]) if len(sys.argv) > 4 else 4.0
ids = z["entry_ids"]
held_entry = np.array([dates.get(i, "0") >= cutoff for i in ids])
held = held_entry[z["entry"]]
print(f"cutoff {cutoff}: {int(held_entry.sum())} held-out entries ({int(held.sum())} residues), "
      f"{int((~held_entry).sum())} reference entries ({int((~held).sum())} residues)")


def build(mask, ci):
    m = (z["cls"] == ci) & mask
    h, _, _ = np.histogram2d(z["phi"][m], z["psi"][m], bins=180, range=[[-180, 180]] * 2)
    d = gaussian_filter(h, sigma / 2.0, mode="wrap")
    d /= d.sum()
    order = np.sort(d.ravel())[::-1]
    cum = np.cumsum(order)
    peak = d.max()
    return d / peak, order[np.searchsorted(cum, 0.98)] / peak, order[np.searchsorted(cum, 0.9995)] / peak


def density(g, phi, psi):
    x, y = (phi + 180) / 2 - 0.5, (psi + 180) / 2 - 0.5
    i, j = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x - i, y - j
    i0, i1, j0, j1 = i % 180, (i + 1) % 180, j % 180, (j + 1) % 180
    return g[i0, j0] * (1 - fx) * (1 - fy) + g[i1, j0] * fx * (1 - fy) + g[i0, j1] * (1 - fx) * fy + g[i1, j1] * fx * fy


def score(d, fav, alw):
    return np.where(d >= fav, 0, np.where(d >= alw, 1, 2))     # 0 favored, 1 allowed, 2 outlier


top = np.load(sys.argv[5]) if len(sys.argv) > 5 else None
rows = []
tot = {"train": np.zeros(3), "held": np.zeros(3), "top8000 on held": np.zeros(3)}
for ci, name in enumerate(CLASSES):
    g, fav, alw = build(~held, ci)
    for tag, mask in (("train", ~held), ("held", held)):
        m = (z["cls"] == ci) & mask
        c = np.bincount(score(density(g, z["phi"][m], z["psi"][m]), fav, alw), minlength=3)
        tot[tag] += c
        rows.append((name, tag, int(m.sum()), *(100.0 * c / max(c.sum(), 1))))
    if top is not None:
        m = (z["cls"] == ci) & held
        tf, ta = 0.02, ALLOWED8000[name]
        c = np.bincount(score(density(top[name].astype(float), z["phi"][m], z["psi"][m]), tf, ta), minlength=3)
        tot["top8000 on held"] += c
        rows.append((name, "top8000 on held", int(m.sum()), *(100.0 * c / max(c.sum(), 1))))
print(f"{'class':9s} {'set':16s} {'n':>8s} {'favored%':>9s} {'allowed%':>9s} {'outlier%':>9s}")
for r in rows:
    print(f"{r[0]:9s} {r[1]:16s} {r[2]:8d} {r[3]:9.2f} {r[4]:9.2f} {r[5]:9.3f}")
for k, c in tot.items():
    if c.sum():
        p = 100.0 * c / c.sum()
        print(f"{'ALL':9s} {k:16s} {int(c.sum()):8d} {p[0]:9.2f} {p[1]:9.2f} {p[2]:9.3f}")
