"""Step 1: extract side-chain chi angles of high-quality residues from the downloaded mmCIF files.

Same entry filter as the Ramachandran reference (X-ray, R-free <= RFREE_MAX, resolution <= 1.85 A). A residue is kept
when every atom used in its chi definitions has occupancy 1, no alternate location and B < BMAX. Output npz: table index,
chi array (n x 4, NaN where absent, 0-360 and folded as in ifacekit.rotamer), entry index and entry ids.
Usage: python 01_extract_chi.py SRC_DIR OUT.npz [BMAX=30] [NPROC=16] [RFREE_MAX=0.23]
"""
import glob
import multiprocessing as mp
import os
import sys

import gemmi
import numpy as np

from ifacekit import rotamer

BMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 30.0
RFREE = float(sys.argv[5]) if len(sys.argv) > 5 else 0.23


def one(arg):
    idx, path = arg
    rows = []
    try:
        v = gemmi.cif.read(path).sole_block().find_value("_refine.ls_R_factor_R_free")
        if v in (None, "?", ".") or float(v) > RFREE:
            return rows
        st = gemmi.read_structure(path)
        if st.resolution and st.resolution > 1.85:
            return rows
        model = st[0]
    except Exception:
        return rows
    for ch in model:
        for res in ch:
            if res.name not in rotamer.CHI or res.het_flag == "H":
                continue
            c = rotamer.chi_values(res, strict=True, bmax=BMAX)
            if c is None:
                continue
            row = [rotamer.TABLES.index(rotamer.TABLE_OF[res.name])] + c + [float("nan")] * (4 - len(c)) + [idx]
            rows.append(row)
    return rows


if __name__ == "__main__":
    files = sorted(glob.glob(sys.argv[1] + "/*.cif.gz"))
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    allrows = []
    with mp.Pool(nproc) as pool:
        for k, r in enumerate(pool.imap_unordered(one, list(enumerate(files)), chunksize=8), 1):
            allrows += r
            if k % 2000 == 0:
                print(f"{k}/{len(files)} files, {len(allrows)} residues", file=sys.stderr, flush=True)
    a = np.array(allrows, dtype=np.float32)
    ids = np.array([os.path.basename(f)[:-7] for f in files])
    np.savez_compressed(sys.argv[2], table=a[:, 0].astype(np.int8), chi=a[:, 1:5], entry=a[:, 5].astype(np.int32),
                        entry_ids=ids, bmax=BMAX, rfree=RFREE, files=len(files))
    print("residues per table:", {t: int((a[:, 0] == i).sum()) for i, t in enumerate(rotamer.TABLES)})
