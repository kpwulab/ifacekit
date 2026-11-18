"""Ramachandran analysis with ifacekit's own reference distributions (no external program).

Residue classes follow the MolProbity convention so that results can be compared: general, glycine,
cis-proline, trans-proline, pre-proline, isoleucine/valine. The reference was built from high-resolution PDB
entries (tools/rama_reference/, summary in ifacekit/data/README.txt). A residue is favored when its interpolated
density is above the level enclosing 98% of the reference residues of its class, allowed above the 99.95% level,
and an outlier below. use_reference(path) selects another grid file (for example Top8000 tables packed by
tools/build_rama_data.py, which carry no thresholds and use the fixed 0.02 / class cut-offs instead).
"""
import math
import os

import gemmi
import numpy as np

CLASSES = ("general", "glycine", "cispro", "transpro", "prepro", "ileval")
LABEL = {"general": "General", "glycine": "Glycine", "cispro": "cis-Pro", "transpro": "trans-Pro",
         "prepro": "Pre-Pro", "ileval": "Ile/Val"}
FAVORED = 0.02                              # fixed cut-offs, used only for reference files without thresholds (Top8000)
ALLOWED = {"general": 0.0005, "glycine": 0.001, "cispro": 0.002, "transpro": 0.001, "prepro": 0.001, "ileval": 0.001}
FAV_T, ALW_T = {}, {}                       # active thresholds per class
_grids = None
_name = ""
_PACKAGED = os.path.join(os.path.dirname(__file__), "data", "rama_ifacekit.npz")


def use_reference(path=None, name=None):
    """Select the reference grids; None = the packaged ifacekit reference."""
    global _grids, _name
    z = np.load(path or _PACKAGED)
    _grids = {k: z[k].astype(float) for k in CLASSES}
    for k in CLASSES:
        if k + "_fav" in z:
            FAV_T[k], ALW_T[k] = float(z[k + "_fav"]), float(z[k + "_alw"])
        else:
            FAV_T[k], ALW_T[k] = FAVORED, ALLOWED[k]
    if name:
        _name = name
    elif "n_entries" in z:
        _name = (f"ifacekit reference ({int(z['n_entries'])} PDB entries, R-free <= {float(z['rfree']):g}, "
                 f"B < {float(z['bmax']):g}, sigma {float(z['sigma']):g} deg; 98% / 99.95% contours)")
    else:
        _name = os.path.basename(path)


def grids():
    if _grids is None:
        use_reference()
    return _grids


def density(cls, phi, psi):
    """Bilinear interpolation on the periodic 2-degree grid (bin centres at -179, -177, ...)."""
    g = grids()[cls]
    x, y = (phi + 180.0) / 2.0 - 0.5, (psi + 180.0) / 2.0 - 0.5
    i, j = int(math.floor(x)), int(math.floor(y))
    fx, fy = x - i, y - j
    i0, i1, j0, j1 = i % 180, (i + 1) % 180, j % 180, (j + 1) % 180
    return float(g[i0, j0] * (1 - fx) * (1 - fy) + g[i1, j0] * fx * (1 - fy)
                 + g[i0, j1] * (1 - fx) * fy + g[i1, j1] * fx * fy)


def classify_density(cls, d):
    grids()
    if d >= FAV_T[cls]:
        return "favored"
    return "allowed" if d >= ALW_T[cls] else "outlier"


def residue_class(res, nxt, cis):
    n = res.name
    if n == "GLY":
        return "glycine"
    if n == "PRO":
        return "cispro" if cis else "transpro"
    if nxt is not None and nxt.name == "PRO":
        return "prepro"
    if n in ("ILE", "VAL"):
        return "ileval"
    return "general"


def _atom(res, name):
    a = res.find_atom(name, "*")
    return a.pos if a else None


def _linked(a, b):
    c, n = _atom(a, "C"), _atom(b, "N")
    return c is not None and n is not None and c.dist(n) < 2.0


def analyze_structure(st, model_index=0, chains=None):
    """Return list of per-residue dicts for polymer residues with a defined phi and psi."""
    out = []
    model = st[model_index]
    for ch in model:
        if chains and ch.name not in chains:
            continue
        pol = [r for r in ch if gemmi.find_tabulated_residue(r.name) is not None
               and gemmi.find_tabulated_residue(r.name).is_amino_acid()]
        for k in range(1, len(pol) - 1):
            prev, res, nxt = pol[k - 1], pol[k], pol[k + 1]
            if not (_linked(prev, res) and _linked(res, nxt)):
                continue
            pts = [_atom(prev, "C"), _atom(res, "N"), _atom(res, "CA"), _atom(res, "C"), _atom(nxt, "N")]
            if any(p is None for p in pts):
                continue
            phi = math.degrees(gemmi.calculate_dihedral(*pts[0:4]))
            psi = math.degrees(gemmi.calculate_dihedral(*pts[1:5]))
            cis = False
            if res.name == "PRO":
                ca0 = _atom(prev, "CA")
                if ca0 is not None:
                    om = math.degrees(gemmi.calculate_dihedral(ca0, pts[0], pts[1], pts[2]))
                    cis = -30 < om < 30
            cls = residue_class(res, nxt, cis)
            d = density(cls, phi, psi)
            out.append({"residue": f"{ch.name}:{res.seqid.num}{res.seqid.icode.strip()}:{res.name}", "chain": ch.name,
                        "class": cls, "phi": round(phi, 1), "psi": round(psi, 1),
                        "density": round(d, 5), "category": classify_density(cls, d)})
    return out


def summarize(rows):
    n = len(rows)
    cnt = {c: sum(1 for r in rows if r["category"] == c) for c in ("favored", "allowed", "outlier")}
    pct = {k: (100.0 * v / n if n else 0.0) for k, v in cnt.items()}
    by_class = {}
    for c in CLASSES:
        sub = [r for r in rows if r["class"] == c]
        if sub:
            by_class[c] = {"n": len(sub), **{k: sum(1 for r in sub if r["category"] == k)
                                             for k in ("favored", "allowed", "outlier")}}
    return {"n": n, "counts": cnt, "percent": {k: round(v, 2) for k, v in pct.items()},
            "by_class": by_class, "outliers": [r for r in rows if r["category"] == "outlier"],
            "reference": _name}


def analyze(path, model_index=0, chains=None, models=None):
    """models: indices (from 0) of the models to analyse together, for an NMR ensemble; every row then carries "model" (counted from 1)."""
    st = gemmi.read_structure(path)
    st.setup_entities()
    if models is None:
        rows = analyze_structure(st, model_index, chains)
        return summarize(rows), rows
    rows = []
    for i in models:
        for r in analyze_structure(st, i, chains):
            r["model"] = i + 1
            rows.append(r)
    summ = summarize(rows)
    summ["n_models"] = len(models)
    return summ, rows
