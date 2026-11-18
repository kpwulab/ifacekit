"""Step 2: smoothed density grids per rotamer table (periodic Gaussian kernel, peak = 1).
Layout and smoothing widths come from ifacekit.rotamer (LAYOUT, SIGMA_BY_DIM).
Cut-offs: per type, the density levels that enclose FAV_ENC and ALW_ENC (fractions, from enclosure.json) of the reference
residues. enclosure.json = {"arg": [0.865, 0.981], ...}; without the file, peak-relative 0.02 / 0.003 are stored instead.
Usage: python 02_build_grids.py chi.npz OUT.npz [enclosure.json|-] [HELDOUT_CUTOFF release_dates.tsv]
With the two extra arguments only entries released before the cutoff date are used (for the held-out test)."""
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

from ifacekit import rotamer

import json
import os
SCALE = float(os.environ.get("ROTA_SIGMA_SCALE", "1"))   # multiplies the smoothing widths of ifacekit.rotamer
z = np.load(sys.argv[1])
enc = json.load(open(sys.argv[3])) if len(sys.argv) > 3 and sys.argv[3] != "-" else None
mask = np.ones(len(z["table"]), dtype=bool)
if len(sys.argv) > 5:
    dates = dict(l.split() for l in open(sys.argv[5]))
    held = np.array([dates.get(i, "0") >= sys.argv[4] for i in z["entry_ids"]])
    mask = ~held[z["entry"]]
out, info = {}, {}
for ti, t in enumerate(rotamer.TABLES):
    lay = rotamer.LAYOUT[t]
    m = (z["table"] == ti) & mask
    data = z["chi"][m][:, :len(lay)].astype(np.float64)
    h, _ = np.histogramdd(data, bins=[nb for _, nb in lay], range=[(0, rng) for rng, _ in lay])
    sig = [SCALE * rotamer.SIGMA_BY_DIM[len(lay)] / (rng / nb) for rng, nb in lay]
    d = gaussian_filter(h, sig, mode="wrap")
    d = d / d.max()
    out[t] = d.astype(np.float16)      # half precision keeps the file small; cut-offs are far above its resolution
    if enc:
        order = np.sort(d.ravel())[::-1]
        cum = np.cumsum(order) / order.sum()
        out[t + "_fav"] = order[min(np.searchsorted(cum, enc[t][0]), len(order) - 1)]
        out[t + "_alw"] = order[min(np.searchsorted(cum, enc[t][1]), len(order) - 1)]
    else:
        out[t + "_fav"], out[t + "_alw"] = rotamer.FAVORED, rotamer.ALLOWED
    info[t] = int(m.sum())
np.savez_compressed(sys.argv[2], bmax=float(z["bmax"]), rfree=float(z["rfree"]),
                    n_entries=int(z["files"]), sigma_scale=SCALE, **out)
print(info)
