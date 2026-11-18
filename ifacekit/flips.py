"""Flip check for Asn, Gln and His side chains (heavy atoms only, no hydrogens).

The two end atoms of an Asn or Gln amide (O and N) cannot be told apart in an electron density map, and the imidazole of His
can be turned by 180 degrees about chi2. For each such side chain the environment is scored twice, as modelled and with the
end atoms swapped, and a flip is suggested when the swapped assignment makes clearly better hydrogen bonds and fewer
polar mismatches or overlaps. Hydrogens of an amide nitrogen are placed from the geometry of the group. This is a simpler test than a program that adds
hydrogens and optimises the whole hydrogen-bond network (Reduce). It was checked against Reduce (tools/compare_reduce.py): on 2024+
entries a suggested flip (gain of 2 or more) was also a Reduce flip in 29% of cases, and 22% of Reduce's flips were found, so it
only suggests; look at the density. It uses the model only: no crystal contacts.
A suggested flip changes the chi2 (Asn, His) or chi3 (Gln) of the rotamer by 180 degrees, which ifacekit.rotamer also reports.
"""
import gemmi
import numpy as np
from scipy.spatial import cKDTree

from . import chem

# end atoms that swap: names whose positions exchange roles when the group is flipped
AMIDE = {"ASN": ("OD1", "ND2"), "GLN": ("OE1", "NE2")}
HIS_N = ("ND1", "NE2")                   # modelled ring nitrogens; after the flip the nitrogens sit at CD2 and CE1
HIS_C = ("CD2", "CE1")
HB_MAX = 3.5
MARGIN = {"ASN": (1.0, 2.0), "GLN": (1.0, 2.0), "HIS": (1.0, 2.0)}   # (uncertain from, flip from); set on 70 pre-2024 entries


def _environment(model):
    """Heavy atoms of the model (first conformer): positions, element, role, residue key."""
    pos, rows = [], []
    for ch in model:
        for res in ch:
            seen = set()
            water = res.name in chem.WATER_NAMES
            prot = (gemmi.find_tabulated_residue(res.name) is not None and gemmi.find_tabulated_residue(res.name).is_amino_acid())
            nuc = res.name in chem.STANDARD_NUC
            for a in res:
                if a.name in seen or a.element.name == "H":
                    continue
                seen.add(a.name)
                d, ac, q, amb = chem.atom_role(res.name, a.name, a.element.name, prot, nuc, water)
                pos.append([a.pos.x, a.pos.y, a.pos.z])
                rows.append((ch.name, res.seqid.num, res.seqid.icode, res.name, a.name, a.element.name, d, ac, amb))
    return np.array(pos), rows


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else v


def _angle(a, b, c):
    """Angle at b between a and c, degrees."""
    u, v = _unit(a - b), _unit(c - b)
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v), -1.0, 1.0))))


def _amide_h(c, n, o):
    """The two hydrogens of an amide nitrogen placed in the plane of the group (C-N-H 120 degrees, N-H 1.01 A)."""
    u = _unit(n - c)
    normal = _unit(np.cross(o - c, n - c))
    w = _unit(np.cross(normal, u))
    return [n + 1.01 * (0.5 * u + 0.866 * w), n + 1.01 * (0.5 * u - 0.866 * w)]


def _score(point, kind, own, tree, pos, rows, geom=None):
    """Score one end atom placed at point.
    kind: 'A' amide oxygen (acceptor), 'D' amide nitrogen (donor), 'B' His ring nitrogen (both), 'C' carbon.
    geom: for 'A' the carbon bonded to it; for 'D' the two hydrogen positions; for 'B' the outward direction."""
    s = 0.0
    for j in tree.query_ball_point(point, 3.7):
        r = rows[j]
        if (r[0], r[1], r[2]) == own:
            continue
        x = pos[j]
        d = float(np.linalg.norm(x - point))
        polar_n = r[5] in ("N", "O")
        nd, na, amb = r[6], r[7], r[8]
        if kind == "C":
            if (polar_n and d < 3.1) or (not polar_n and d < 3.2):
                s -= 0.5
            continue
        if not polar_n:
            if d < 3.1:
                s -= 1.0                      # non-polar atom pressing on a polar end atom
            continue
        if d < 2.5:
            s -= 1.0
            continue
        gain = 0.0
        if kind == "A":
            if nd:
                if _angle(geom, point, x) >= 100.0:
                    gain = 1.0 if d <= 3.2 else 0.5 if d <= HB_MAX else 0.0
            elif not amb and na and d < 3.3:
                s -= 1.0
        elif kind == "D":
            if na:
                for h in geom:
                    dh = float(np.linalg.norm(x - h))
                    if dh <= 2.5 and _angle(point, h, x) >= 120.0:
                        gain = max(gain, 1.0)
                    elif dh <= 2.8 and _angle(point, h, x) >= 110.0:
                        gain = max(gain, 0.5)
            elif not amb and nd and d < 3.3:
                s -= 1.0
        else:                                   # His ring nitrogen: as donor (H along geom) or as acceptor (lone pair along geom)
            if na or amb:
                h = point + 1.0 * geom
                dh = float(np.linalg.norm(x - h))
                if dh <= 2.5 and _angle(point, h, x) >= 120.0:
                    gain = max(gain, 1.0)
                elif dh <= 2.8 and _angle(point, h, x) >= 110.0:
                    gain = max(gain, 0.5)
            if nd or amb:
                if _angle(point + geom, point, x) <= 60.0 and d <= HB_MAX:
                    gain = max(gain, 1.0 if d <= 3.2 else 0.5)
        s += gain
    return s


def analyze_structure(st, model_index=0, chains=None):
    model = st[model_index]
    pos, rows = _environment(model)
    tree = cKDTree(pos)
    out = []
    for ch in model:
        if chains and ch.name not in chains:
            continue
        for res in ch:
            if res.name not in ("ASN", "GLN", "HIS"):
                continue
            own = (ch.name, res.seqid.num, res.seqid.icode)
            g = {}
            for n in (HIS_N + HIS_C if res.name == "HIS" else AMIDE[res.name]):
                a = res.find_atom(n, "*")
                if a is not None:
                    g[n] = np.array([a.pos.x, a.pos.y, a.pos.z])
            if res.name == "HIS":
                if len(g) < 4:
                    continue
                ctr = (g["ND1"] + g["NE2"] + g["CD2"] + g["CE1"]) / 4.0
                ring = [(g["ND1"], "B"), (g["NE2"], "B"), (g["CD2"], "C"), (g["CE1"], "C")]
                flip = [(g["CD2"], "B"), (g["CE1"], "B"), (g["ND1"], "C"), (g["NE2"], "C")]
                so = sum(_score(p, k, own, tree, pos, rows, _unit(p - ctr) if k == "B" else None) for p, k in ring)
                sf = sum(_score(p, k, own, tree, pos, rows, _unit(p - ctr) if k == "B" else None) for p, k in flip)
            else:
                o, n = AMIDE[res.name]
                cname = "CG" if res.name == "ASN" else "CD"
                ca = res.find_atom(cname, "*")
                if o not in g or n not in g or ca is None:
                    continue
                c = np.array([ca.pos.x, ca.pos.y, ca.pos.z])
                po, pn = g[o], g[n]
                so = _score(po, "A", own, tree, pos, rows, c) + _score(pn, "D", own, tree, pos, rows, _amide_h(c, pn, po))
                sf = _score(pn, "A", own, tree, pos, rows, c) + _score(po, "D", own, tree, pos, rows, _amide_h(c, po, pn))
            delta = sf - so
            u, f = MARGIN[res.name]
            verdict = "flip" if delta >= f else "uncertain" if delta >= u else "keep"
            out.append({"residue": f"{ch.name}:{res.seqid.num}{res.seqid.icode.strip()}:{res.name}", "chain": ch.name,
                        "type": res.name, "score_modelled": round(so, 2), "score_flipped": round(sf, 2),
                        "delta": round(delta, 2), "verdict": verdict})
    return out


def summarize(rows):
    cnt = {k: sum(1 for r in rows if r["verdict"] == k) for k in ("keep", "uncertain", "flip")}
    return {"n": len(rows), "counts": cnt, "suggested": [r for r in rows if r["verdict"] == "flip"],
            "uncertain": [r for r in rows if r["verdict"] == "uncertain"],
            "method": "heavy-atom hydrogen-bond and overlap score, modelled against flipped (no hydrogens, no crystal contacts)"}


def analyze(path, model_index=0, chains=None):
    st = gemmi.read_structure(path)
    st.setup_entities()
    rows = analyze_structure(st, model_index, chains)
    return summarize(rows), rows
