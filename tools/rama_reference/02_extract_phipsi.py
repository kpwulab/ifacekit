"""Step 2: extract phi/psi of high-quality residues from the downloaded mmCIF files.

A residue is kept when it and its two neighbours are linked (C-N < 2.0 A), the backbone atoms N, CA, C, O
of the residue have occupancy 1, no alternate location and B < BMAX, and it is an L-amino acid with all of N, CA, C.
Residue classes follow ifacekit.rama (general, glycine, cis/trans-proline, pre-proline, Ile/Val).
Usage: python 02_extract_phipsi.py SRC_DIR OUT.npz [BMAX=30] [NPROC=16] [RFREE_MAX=none]
With RFREE_MAX, entries without an R-free value or above it are skipped (read from _refine.ls_R_factor_R_free).
"""
import glob
import math
import multiprocessing as mp
import sys

import gemmi
import numpy as np

from ifacekit import rama

BMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 30.0
BB = ("N", "CA", "C", "O")
RFREE = float(sys.argv[5]) if len(sys.argv) > 5 else None


def rfree_ok(path):
    if RFREE is None:
        return True
    try:
        v = gemmi.cif.read(path).sole_block().find_value("_refine.ls_R_factor_R_free")
        return v not in (None, "?", ".") and float(v) <= RFREE
    except Exception:
        return False


def good(res):
    for n in BB:
        a = res.find_atom(n, "*")
        if a is None or a.altloc != "\0" or a.occ < 0.99 or a.b_iso >= BMAX:
            return False
    return True


def one(arg):
    idx, path = arg
    rows = []
    if not rfree_ok(path):
        return rows
    try:
        st = gemmi.read_structure(path)
        if st.resolution and st.resolution > 1.85:
            return rows
        ch_list = st[0]
    except Exception:
        return rows
    for ch in ch_list:
        pol = [r for r in ch if (gemmi.find_tabulated_residue(r.name) or None) is not None
               and gemmi.find_tabulated_residue(r.name).is_amino_acid() and r.het_flag != "H"]
        for k in range(1, len(pol) - 1):
            prev, res, nxt = pol[k - 1], pol[k], pol[k + 1]
            if not (rama._linked(prev, res) and rama._linked(res, nxt)) or not good(res):
                continue
            pts = [rama._atom(prev, "C"), rama._atom(res, "N"), rama._atom(res, "CA"), rama._atom(res, "C"), rama._atom(nxt, "N")]
            if any(p is None for p in pts):
                continue
            phi = math.degrees(gemmi.calculate_dihedral(*pts[0:4]))
            psi = math.degrees(gemmi.calculate_dihedral(*pts[1:5]))
            cis = False
            if res.name == "PRO":
                ca0 = rama._atom(prev, "CA")
                if ca0 is None:
                    continue
                cis = -30 < math.degrees(gemmi.calculate_dihedral(ca0, pts[0], pts[1], pts[2])) < 30
            rows.append((rama.CLASSES.index(rama.residue_class(res, nxt, cis)), phi, psi, idx))
    return rows


if __name__ == "__main__":
    files = sorted(glob.glob(sys.argv[1] + "/*.cif.gz"))
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    allrows, done = [], 0
    with mp.Pool(nproc) as pool:
        for r in pool.imap_unordered(one, list(enumerate(files)), chunksize=8):
            allrows += r
            done += 1
            if done % 1000 == 0:
                print(f"{done}/{len(files)} files, {len(allrows)} residues", file=sys.stderr, flush=True)
    a = np.array(allrows, dtype=np.float32)
    import os
    ids = np.array([os.path.basename(f)[:-7] for f in files])
    np.savez_compressed(sys.argv[2], rfree=(RFREE or 0.0), cls=a[:, 0].astype(np.int8), phi=a[:, 1], psi=a[:, 2], entry=a[:, 3].astype(np.int32), entry_ids=ids, bmax=BMAX, files=len(files))
    print("residues per class:", {rama.CLASSES[i]: int((a[:, 0] == i).sum()) for i in range(6)})
