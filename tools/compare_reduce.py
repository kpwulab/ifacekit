"""Compare ifacekit's Asn/Gln/His flip check with Reduce (development only).

For each mmCIF file: write a PDB without hydrogens, run Reduce (-BUILD, which optimises flips), read the "USER  MOD" lines
(":FLIP" = Reduce flipped the group; K / C / X / F = keep / clashes / uncertain / flip, with the o and f scores), and run
ifacekit.flips on the same file. Output: one TSV row per residue.
Usage: python compare_reduce.py OUT.tsv FILE.cif[.gz] ...   (env REDUCE and REDUCE_DB: the reduce program and its het dictionary; set REDUCE_DB to the full path of the dictionary)
"""
import multiprocessing as mp
import os
import re
import subprocess
import sys
import tempfile

import gemmi

from ifacekit import flips

REDUCE = os.environ.get("REDUCE", "reduce")
DB = os.environ.get("REDUCE_DB", "reduce_wwPDB_het_dict.txt")
LINE = re.compile(r"USER  MOD (?:Set\s*\d+\.\d+|Single\s*)\s*:\s*(\S)\s+(-?\d+)([A-Za-z]?)\s+(ASN|GLN|HIS)\s+:(FLIP)?")
SCORE = re.compile(r"([KCXF])\(o=(-?[\d.]+)!?,f=(-?[\d.]+)!?\)")


def one(path):
    rows = []
    try:
        st = gemmi.read_structure(path)
        st.setup_entities()
        st.remove_hydrogens()
        if any(len(c.name) > 1 for c in st[0]):
            return rows
        with tempfile.TemporaryDirectory() as td:
            pdb = os.path.join(td, "m.pdb")
            st.write_pdb(pdb)
            r = subprocess.run([REDUCE, "-BUILD", "-DB", DB, pdb], capture_output=True, text=True, timeout=600)
        ref = {}
        for l in r.stdout.splitlines():
            m = LINE.match(l)
            if not m:
                continue
            s = SCORE.search(l)
            ref[(m.group(1).strip(), int(m.group(2)), m.group(4))] = (1 if m.group(5) else 0, s.group(1) if s else "?",
                                                                      float(s.group(2)) if s else 0.0, float(s.group(3)) if s else 0.0)
        _, mine = flips.analyze(path)
        for x in mine:
            ch, num, name = x["residue"].split(":")
            k = (ch, int("".join(c for c in num if c.isdigit() or c == "-")), name)
            if k in ref:
                f, cat, so, sf = ref[k]
                rows.append((os.path.basename(path).split(".")[0], x["residue"], x["type"], f, cat, so, sf,
                             x["score_modelled"], x["score_flipped"], x["delta"]))
    except Exception as e:
        print("failed", path, e, file=sys.stderr)
    return rows


if __name__ == "__main__":
    with mp.Pool(8) as pool, open(sys.argv[1], "w") as out:
        out.write("entry\tresidue\ttype\treduce_flip\treduce_cat\treduce_o\treduce_f\tours_modelled\tours_flipped\tours_delta\n")
        for rows in pool.imap_unordered(one, sys.argv[2:]):
            for r in rows:
                out.write("\t".join(str(x) for x in r) + "\n")
