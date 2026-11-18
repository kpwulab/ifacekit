"""Step 3: turn the extracted phi/psi values into reference grids (own Ramachandran distributions).

Per residue class: wrapped (periodic) Gaussian kernel density on a 2-degree grid. The favored and allowed
cut-offs are the density levels that enclose 98% and 99.95% of the reference residues (MolProbity's
convention). Output npz holds for every class: grid (180x180, peak = 1), favored threshold, allowed threshold,
residue count. Usage: python 03_build_grids.py phipsi.npz OUT.npz [SIGMA_DEG=6] [CUTOFF release_dates.tsv]
With the last two arguments only entries released before the cutoff date are used (grids for tests on newer entries).
"""
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

CLASSES = ("general", "glycine", "cispro", "transpro", "prepro", "ileval")
src = np.load(sys.argv[1])
sigma = float(sys.argv[3]) if len(sys.argv) > 3 else 6.0
keep = np.ones(len(src["cls"]), dtype=bool)
if len(sys.argv) > 5:
    dates = dict(l.split() for l in open(sys.argv[5]))
    held = np.array([dates.get(i, "0") >= sys.argv[4] for i in src["entry_ids"]])
    keep = ~held[src["entry"]]
out, info = {}, {}
for ci, name in enumerate(CLASSES):
    m = (src["cls"] == ci) & keep
    phi, psi = src["phi"][m], src["psi"][m]
    h, _, _ = np.histogram2d(phi, psi, bins=180, range=[[-180, 180], [-180, 180]])
    d = gaussian_filter(h, sigma / 2.0, mode="wrap")       # grid step is 2 degrees
    d /= d.sum()
    order = np.sort(d.ravel())[::-1]
    cum = np.cumsum(order)
    t_fav = order[np.searchsorted(cum, 0.98)]
    t_alw = order[np.searchsorted(cum, 0.9995)]
    peak = d.max()
    out[name] = (d / peak).astype(np.float32)
    out[name + "_fav"], out[name + "_alw"], out[name + "_n"] = t_fav / peak, t_alw / peak, int(m.sum())
    info[name] = (int(m.sum()), round(float(t_fav / peak), 4), round(float(t_alw / peak), 5))
np.savez_compressed(sys.argv[2], sigma=sigma, bmax=float(src["bmax"]), files=int(src["files"]), **out)
for k, v in info.items():
    print(f"{k:9s} n={v[0]:>8d}  favored>={v[1]}  allowed>={v[2]}")
