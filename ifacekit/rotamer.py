"""Side-chain rotamer analysis with ifacekit's own reference distributions (no external program).

For each residue with chi angles the chi values are scored against a smoothed density built from high-resolution PDB
entries (tools/rotamer_reference/, summary in ifacekit/data/README.txt). Residue types and chi definitions follow the
MolProbity convention so results can be compared (Phe and Tyr share one table; chi2 of Asp and Phe/Tyr and chi3 of Glu are
folded to 0-180 degrees because the end groups are symmetric). A residue is favored when the density at its chi
values is at or above the favored level of its type, allowed at or above the allowed level, and an outlier below. The levels are
stored with the reference grids (density relative to the peak; FAVORED and ALLOWED are the defaults for files without levels). Ala and Gly have no chi and are not scored.
No hydrogens are used and Asn, Gln and His are not flipped, so a call on these three can differ from a program that
tests flips.
"""
import math
import os

import gemmi
import numpy as np
from scipy.ndimage import map_coordinates

FAVORED = 0.02
ALLOWED = 0.003

# residue name -> list of chi definitions (4 atom names each)
CHI = {
    "ARG": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "NE"), ("CG", "CD", "NE", "CZ")],
    "ASN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "ASP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "CYS": [("N", "CA", "CB", "SG")],
    "GLN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
    "GLU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
    "HIS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")],
    "ILE": [("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")],
    "LEU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "LYS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "CE"), ("CG", "CD", "CE", "NZ")],
    "MET": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "SD"), ("CB", "CG", "SD", "CE")],
    "PHE": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "PRO": [("N", "CA", "CB", "CG")],
    "SER": [("N", "CA", "CB", "OG")],
    "THR": [("N", "CA", "CB", "OG1")],
    "TRP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "TYR": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "VAL": [("N", "CA", "CB", "CG1")],
}
TABLE_OF = {k: k.lower() for k in CHI}
TABLE_OF["PHE"] = TABLE_OF["TYR"] = "phetyr"
TABLES = ("arg", "asn", "asp", "cys", "gln", "glu", "his", "ile", "leu", "lys", "met", "phetyr", "pro", "ser", "thr", "trp", "val")
# per table: (range in degrees, number of bins) for each chi; range 180 = folded (symmetric end group)
LAYOUT = {
    "arg": [(360, 36)] * 4, "lys": [(360, 36)] * 4,
    "gln": [(360, 45)] * 3, "met": [(360, 45)] * 3, "glu": [(360, 45), (360, 45), (180, 23)],
    "asn": [(360, 72)] * 2, "his": [(360, 72)] * 2, "ile": [(360, 72)] * 2, "leu": [(360, 72)] * 2, "trp": [(360, 72)] * 2,
    "asp": [(360, 72), (180, 36)], "phetyr": [(360, 72), (180, 36)],
    "cys": [(360, 360)], "pro": [(360, 360)], "ser": [(360, 360)], "thr": [(360, 360)], "val": [(360, 360)],
}
# smoothing width (degrees) by number of dimensions, used by the build script
SIGMA_BY_DIM = {1: 5.0, 2: 5.0, 3: 8.0, 4: 10.0}
LABEL = {t: t.capitalize() if t != "phetyr" else "Phe/Tyr" for t in TABLES}
_PACKAGED = os.path.join(os.path.dirname(__file__), "data", "rotamer_ifacekit.npz")
_grids, _name, _thr = None, "", {}


def use_reference(path=None, name=None):
    """Select the reference grids; None = the packaged ifacekit reference."""
    global _grids, _name, _thr
    z = np.load(path or _PACKAGED)
    _grids = {t: z[t].astype(np.float64) for t in TABLES}
    _thr = {t: ((float(z[t + "_fav"]), float(z[t + "_alw"])) if t + "_fav" in z else (FAVORED, ALLOWED)) for t in TABLES}
    if name:
        _name = name
    elif "n_entries" in z:
        _name = (f"ifacekit reference ({int(z['n_entries'])} PDB entries, R-free <= {float(z['rfree']):g}, "
                 f"side-chain B < {float(z['bmax']):g}; per-type cut-offs calibrated against MolProbity)")
    else:
        _name = os.path.basename(path)


def grids():
    if _grids is None:
        use_reference()
    return _grids


def reference_name():
    grids()
    return _name


def table_for(resname):
    return TABLE_OF.get(resname)


def conformers(res):
    """Alternate-conformer letters of a residue (sorted), or [None] when it has none."""
    return sorted({a.altloc for a in res if a.altloc != "\0"}) or [None]


def label(row):
    """Residue id of a result row, with the conformer letter when it is an alternate conformer."""
    return row["residue"] + (f" alt {row['altloc']}" if row.get("altloc") else "") + (f" model {row['model']}" if row.get("model") else "")


def chi_values(res, strict=False, bmax=30.0, alt=None):
    """Chi angles (degrees, 0-360, folded where the end group is symmetric) or None when an atom is missing.
    strict: every atom of the chain must have occupancy 1, no alternate location and B below bmax (reference building).
    alt: conformer letter; atoms of that conformer, and the atoms shared by all conformers (default: the first conformer of every atom)."""
    defs = CHI.get(res.name)
    if defs is None:
        return None
    cache, out = {}, []
    for d in defs:
        pts = []
        for n in d:
            if n not in cache:
                a = res.find_atom(n, alt or "*")
                if a is None or (strict and (a.altloc != "\0" or a.occ < 0.99 or a.b_iso >= bmax)):
                    return None
                cache[n] = a.pos
            pts.append(cache[n])
        out.append(math.degrees(gemmi.calculate_dihedral(*pts)) % 360.0)
    t = TABLE_OF[res.name]
    for k, (rng, _) in enumerate(LAYOUT[t][:len(out)]):
        if rng == 180:
            out[k] %= 180.0
    return out


def density(table, chis):
    """Multilinear interpolation on the periodic grid (bin centres at (i + 0.5) * width)."""
    g = grids()[table]
    coords = [[c / (rng / nb) - 0.5] for c, (rng, nb) in zip(chis, LAYOUT[table])]
    return float(map_coordinates(g, coords, order=1, mode="grid-wrap")[0])


def classify_density(table, d):
    grids()
    fav, alw = _thr[table]
    return "favored" if d >= fav else "allowed" if d >= alw else "outlier"


def _is_aa(res):
    t = gemmi.find_tabulated_residue(res.name)
    return t is not None and t.is_amino_acid()


FLIP_CHI = {"ASN": 1, "HIS": 1, "GLN": 2}      # index of the chi angle that changes by 180 degrees when the group is flipped


def analyze_structure(st, model_index=0, chains=None, flip_verdicts=None):
    """Per-residue rotamer rows. For Asn, Gln and His the rotamer of the flipped group is scored as well
    (category_flipped); flip_verdicts maps residue id to the verdict of ifacekit.flips."""
    out = []
    for ch in st[model_index]:
        if chains and ch.name not in chains:
            continue
        for res in ch:
            if res.name not in CHI or not _is_aa(res) or res.het_flag == "H":
                continue
            for conf in conformers(res):               # every alternate conformer is scored on its own, as MolProbity does
                chis = chi_values(res, alt=conf)
                if chis is None:
                    continue
                t = TABLE_OF[res.name]
                d = density(t, chis)
                rid = f"{ch.name}:{res.seqid.num}{res.seqid.icode.strip()}:{res.name}"
                row = {"residue": rid, "chain": ch.name, "table": t, "chi": [round(x, 1) for x in chis],
                       "density": round(d, 5), "category": classify_density(t, d), "altloc": conf or ""}
                if res.name in FLIP_CHI:
                    k = FLIP_CHI[res.name]
                    other = list(chis)
                    other[k] = (other[k] + 180.0) % 360.0
                    da = density(t, other)
                    row.update({"density_flipped": round(da, 5), "category_flipped": classify_density(t, da),
                                "flip": (flip_verdicts or {}).get(rid)})
                out.append(row)
    return out


def _after_flips(rows):
    """Categories after applying the flips that ifacekit.flips suggests with confidence (verdict 'flip')."""
    cat = [r["category_flipped"] if r.get("flip") == "flip" and "category_flipped" in r else r["category"] for r in rows]
    n = len(rows)
    cnt = {c: cat.count(c) for c in ("favored", "allowed", "outlier")}
    return {"counts": cnt, "percent": {k: round(100.0 * v / n, 2) if n else 0.0 for k, v in cnt.items()},
            "n_flips_applied": sum(1 for r in rows if r.get("flip") == "flip")}


def summarize(rows):
    n = len(rows)
    cnt = {c: sum(1 for r in rows if r["category"] == c) for c in ("favored", "allowed", "outlier")}
    pct = {k: (100.0 * v / n if n else 0.0) for k, v in cnt.items()}
    by = {}
    for t in TABLES:
        sub = [r for r in rows if r["table"] == t]
        if sub:
            by[t] = {"n": len(sub), **{k: sum(1 for r in sub if r["category"] == k) for k in ("favored", "allowed", "outlier")}}
    return {"n": n, "counts": cnt, "percent": {k: round(v, 2) for k, v in pct.items()}, "by_type": by,
            "outliers": [r for r in rows if r["category"] == "outlier"], "after_suggested_flips": _after_flips(rows),
            "reference": reference_name()}


def analyze(path, model_index=0, chains=None, models=None):
    """models: indices (from 0) of the models to analyse together, for an NMR ensemble; every row then carries "model" (counted from 1)."""
    st = gemmi.read_structure(path)
    st.setup_entities()
    from . import flips
    rows = []
    for i in ([model_index] if models is None else models):
        fv = {r["residue"]: r["verdict"] for r in flips.analyze_structure(st, i, chains)}
        part = analyze_structure(st, i, chains, fv)
        if models is not None:
            for r in part:
                r["model"] = i + 1
        rows += part
    summ = summarize(rows)
    if models is not None:
        summ["n_models"] = len(models)
    return summ, rows
