"""Held-out test of the rotamer reference.
Grids built only from entries released before the cutoff (02_build_grids.py with cutoff and dates); the residues of entries
released on or after the cutoff are scored. If the reference generalises, the held-out residues fall into favored /
allowed / outlier at about the rates of the reference residues themselves.
Usage: python 04_heldout_test.py chi.npz grids_before_cutoff.npz release_dates.tsv CUTOFF"""
import sys

import numpy as np
from scipy.ndimage import map_coordinates

from ifacekit import rotamer

z = np.load(sys.argv[1])
rotamer.use_reference(sys.argv[2])
dates = dict(l.split() for l in open(sys.argv[3]))
held_entry = np.array([dates.get(i, "0") >= sys.argv[4] for i in z["entry_ids"]])
held = held_entry[z["entry"]]
print(f"cutoff {sys.argv[4]}: {int(held_entry.sum())} held-out entries ({int(held.sum())} residues), {int((~held_entry).sum())} reference entries ({int((~held).sum())} residues)")
G = rotamer.grids()
tot = {"train": np.zeros(3), "held": np.zeros(3)}
print(f"{'type':8s} {'set':6s} {'n':>8s} {'favored%':>9s} {'allowed%':>9s} {'outlier%':>9s}")
for ti, t in enumerate(rotamer.TABLES):
    lay = rotamer.LAYOUT[t]
    fav, alw = rotamer._thr[t]
    for tag, mask in (("train", ~held), ("held", held)):
        m = (z["table"] == ti) & mask
        chi = z["chi"][m][:, :len(lay)].astype(np.float64)
        coords = [chi[:, k] / (rng / nb) - 0.5 for k, (rng, nb) in enumerate(lay)]
        d = map_coordinates(G[t], coords, order=1, mode="grid-wrap")
        c = np.array([(d >= fav).sum(), ((d >= alw) & (d < fav)).sum(), (d < alw).sum()], dtype=float)
        tot[tag] += c
        p = 100 * c / max(c.sum(), 1)
        print(f"{t:8s} {tag:6s} {int(m.sum()):8d} {p[0]:9.2f} {p[1]:9.2f} {p[2]:9.3f}")
for k, c in tot.items():
    p = 100 * c / c.sum()
    print(f"{'ALL':8s} {k:6s} {int(c.sum()):8d} {p[0]:9.2f} {p[1]:9.2f} {p[2]:9.3f}")
