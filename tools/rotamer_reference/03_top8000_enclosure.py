"""Development helper: the enclosure implied by MolProbity's cut-offs on the Top8000 rotamer tables.
For every table, the fraction of the table's density mass at or above 0.02 (favored) and 0.003 (allowed).
Usage: python 03_top8000_enclosure.py /path/to/rotarama_data OUT.json   (only the two fractions per type are kept)"""
import glob
import json
import os
import sys

import numpy as np

FILES = {"phetyr": "rota8000-phetyr.data"}
out = {}
for t in "arg asn asp cys gln glu his ile leu lys met phetyr pro ser thr trp val".split():
    v = np.array([float(l.split()[-1]) for l in open(f"{sys.argv[1]}/{FILES.get(t, f'rota8000-{t}.data')}") if not l.startswith("#") and l.strip()])
    out[t] = [round(float(v[v >= 0.02].sum() / v.sum()), 4), round(float(v[v >= 0.003].sum() / v.sum()), 4)]
json.dump(out, open(sys.argv[2], "w"), indent=1)
print(out)
