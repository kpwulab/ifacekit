"""Two-dimensional protein-ligand interaction diagram (a LigPlot-like picture) drawn as SVG.

Steps. (1) The ligand graph is rebuilt from the 3D coordinates (bonds by covalent radii, bond orders
from planarity and valence), because a model file seldom carries bonds. (2) RDKit draws the ligand in 2D.
(3) That drawing is rotated or mirrored to lie as close as possible to the ligand's own 3D projection, so
that residues can be placed on the side where they really are. (4) Contacting residues are placed around the
ligand and spread apart so that they do not overlap. (5) Hydrogen bonds and salt bridges are dashed lines with
their distance, other contacts are thin lines.
"""
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdDetermineBonds

from . import chem
from . import model as mdl

PX = 44.0                      # pixels per Angstrom in the ligand drawing
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}
POSITIVE = {"LYS", "ARG", "HIS"}
NEGATIVE = {"ASP", "GLU"}
STD_AA = HYDROPHOBIC | {"LYS", "ARG", "HIS", "ASP", "GLU", "SER", "THR", "ASN", "GLN", "TYR", "GLY"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "GLY"}
FILL = {"hydrophobic": "#cfe8c6", "polar": "#cfe3f1", "positive": "#c9d3f5", "negative": "#f6cfcf", "other": "#e6e6e6"}
ELEM_COL = {"O": "#c0392b", "N": "#2b6cb0", "S": "#b7791f", "P": "#d35400", "F": "#2f855a", "CL": "#2f855a", "BR": "#8e44ad"}


def residue_class(name):
    n = name.upper()
    return ("hydrophobic" if n in HYDROPHOBIC else "positive" if n in POSITIVE else "negative" if n in NEGATIVE
            else "polar" if n in POLAR else "other")


# ---------------------------------------------------------------- ligand graph
def _angle(a, b, c):
    v1, v2 = a - b, c - b
    return math.degrees(math.acos(max(-1, min(1, v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)))))


def dictionary_bonds(path, comp_id):
    """Bonds of a chemical component as written in the file itself ({frozenset(names): order}), or None.
    Files from the wwPDB carry this dictionary, so the true bond orders need no guessing."""
    try:
        import gemmi
        doc = gemmi.cif.read(str(path))
        order = {"sing": 1, "doub": 2, "trip": 3, "arom": 1.5}
        for blk in doc:
            t = blk.find("_chem_comp_bond.", ["comp_id", "atom_id_1", "atom_id_2", "value_order", "?pdbx_aromatic_flag"])
            if not t:
                continue
            out = {}
            for row in t:
                if row[0].strip("'\" ").upper() == comp_id.upper():
                    arom = len(row) > 4 and row[4].strip().upper() == "Y"
                    out[frozenset((row[1].strip("'\" "), row[2].strip("'\" ")))] = 1.5 if arom else order.get(row[3].lower(), 1)
            if out:
                return out
    except Exception:
        pass
    return None


def ligand_mol(elements, xyz, names=None, bonds=None):
    """RDKit molecule. With a bond dictionary (from the file) the chemistry is exact; otherwise bonds and bond
    orders are guessed from heavy-atom coordinates."""
    n = len(elements)
    block = f"{n}\n\n" + "\n".join(f"{e.capitalize()} {x:.4f} {y:.4f} {z:.4f}" for e, (x, y, z) in zip(elements, xyz))
    mol = Chem.MolFromXYZBlock(block)
    if bonds and names:
        pos = {nm: i for i, nm in enumerate(names)}
        rw = Chem.RWMol(mol)
        for key, o in bonds.items():
            if len(key) == 2 and all(k in pos for k in key):
                a, b = (pos[k] for k in key)
                rw.AddBond(a, b, Chem.BondType.AROMATIC if o == 1.5 else Chem.BondType.DOUBLE if o == 2 else Chem.BondType.TRIPLE if o == 3 else Chem.BondType.SINGLE)
                if o == 1.5:
                    rw.GetAtomWithIdx(a).SetIsAromatic(True)
                    rw.GetAtomWithIdx(b).SetIsAromatic(True)
                    rw.GetBondBetweenAtoms(a, b).SetIsAromatic(True)
        m = rw.GetMol()
        if m.GetNumBonds() >= n - 1:
            for a in m.GetAtoms():
                a.SetNoImplicit(False)
                a.SetNumRadicalElectrons(0)
            try:
                Chem.SanitizeMol(m)
            except Exception:
                m.UpdatePropertyCache(strict=False)
                Chem.FastFindRings(m)
            m.SetProp("bond_source", "file dictionary")
            return m
        mol = Chem.MolFromXYZBlock(block)
    rdDetermineBonds.DetermineConnectivity(mol)
    nb = [[x.GetIdx() for x in a.GetNeighbors()] for a in mol.GetAtoms()]
    need = [0] * n                       # number of pi bonds the atom should take
    for i, e in enumerate(elements):
        k = len(nb[i])
        e = e.upper()
        if e == "C":
            if k == 3:
                angs = [_angle(xyz[a], xyz[i], xyz[b]) for ai, a in enumerate(nb[i]) for b in nb[i][ai + 1:]]
                need[i] = 1 if sum(angs) > 350 else 0
            elif k == 2:
                need[i] = 1 if 115 < _angle(xyz[nb[i][0]], xyz[i], xyz[nb[i][1]]) < 165 else 0
            elif k == 1:
                need[i] = 0
        elif e == "N":
            if k == 2 and 110 < _angle(xyz[nb[i][0]], xyz[i], xyz[nb[i][1]]) < 165:
                need[i] = 1
        elif e in ("O", "S") and k == 1 and elements[nb[i][0]].upper() in ("C", "S", "P", "N"):
            need[i] = 1
    edges = [(a, b) for a in range(n) for b in nb[a] if a < b and need[a] and need[b]]
    # terminal O or S first, so a carbonyl is always a C=O
    taken, double = set(), []
    for a, b in edges:
        t, o = (a, b) if len(nb[a]) == 1 else (b, a) if len(nb[b]) == 1 else (None, None)
        if t is not None and t not in taken and o not in taken and elements[t].upper() in ("O", "S"):
            taken |= {t, o}
            double.append((a, b))
    rest = [(a, b) for a, b in edges if a not in taken and b not in taken]
    best = []

    def search(k, used, chosen):
        nonlocal best
        if len(chosen) > len(best):
            best = list(chosen)
        if k == len(rest) or len(best) * 2 >= sum(1 for i in range(n) if need[i] and i not in taken):
            return
        a, b = rest[k]
        if a not in used and b not in used:
            chosen.append((a, b))
            search(k + 1, used | {a, b}, chosen)
            chosen.pop()
        search(k + 1, used, chosen)
    if len(rest) < 60:
        search(0, set(), [])
    else:                                  # greedy for very large systems
        used = set()
        for a, b in rest:
            if a not in used and b not in used:
                best.append((a, b))
                used |= {a, b}
    rw = Chem.RWMol(mol)
    for a, b in double + best:
        rw.GetBondBetweenAtoms(a, b).SetBondType(Chem.BondType.DOUBLE)
    m = rw.GetMol()
    for a in m.GetAtoms():
        a.SetNoImplicit(False)
        a.SetNumRadicalElectrons(0)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        m.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(m)
    return m


def kekule(mol):
    """Copy of the molecule with aromatic bonds turned into alternating single and double bonds, for drawing."""
    km = Chem.Mol(mol)
    try:
        Chem.Kekulize(km, clearAromaticFlags=True)
        return km
    except Exception:
        return mol


def depict(mol, xyz):
    """2D coordinates (Angstrom-like units) rotated and, if needed, mirrored to match the 3D projection."""
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    D = np.array([[mol.GetConformer().GetAtomPosition(i).x, mol.GetConformer().GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())])
    c = xyz.mean(axis=0)
    u, s, vt = np.linalg.svd(xyz - c)
    frame = vt[:3]                          # rows: plane axes e1, e2 and normal e3
    P = (xyz - c) @ frame[:2].T
    D = D - D.mean(axis=0)
    best = None
    for flip in (1, -1):
        Df = D * np.array([flip, 1])
        U, S, Vt = np.linalg.svd(Df.T @ P)
        R = U @ Vt
        if np.linalg.det(R) < 0:
            continue
        err = np.sqrt(((Df @ R - P) ** 2).sum(axis=1).mean())
        if best is None or err < best[0]:
            best = (err, Df @ R)
    if best is None:
        best = (0, D)
    return best[1], frame, c


# ---------------------------------------------------------------- layout and drawing
def build(path, result, ligand_unit=None, protein_units=None, vdw_max=4.0, model_index=0, detail="full"):
    m = mdl.load(path, model_index)
    lig = [u for u in m.units if u["kind"] == "ligand"]
    if ligand_unit is None:
        ligand_unit = lig[0]["id"]
    uid = next(k for k, u in enumerate(m.units) if u["id"] == ligand_unit)
    li = np.where(m.unit_idx == uid)[0]
    elems = [m.element[k] for k in li]
    xyz = m.xyz[li]
    names = [m.atom_name[k] for k in li]
    rn0 = m.residues[m.res_idx[li[0]]]["resname"]
    mol = ligand_mol(elems, xyz, names, dictionary_bonds(path, rn0))
    D, frame, cen = depict(mol, xyz)
    atom_pos = {int(k): D[j] * PX for j, k in enumerate(li)}      # atom index in model -> 2D pixel position
    local = {int(k): j for j, k in enumerate(li)}
    # contacts of this ligand with everything else
    res = {}
    for c in result["contacts"]:
        if not c["inter_unit"] or ligand_unit not in (c["unit1"], c["unit2"]):
            continue
        a_l, a_p = (c["atom1"], c["atom2"]) if c["unit1"] == ligand_unit else (c["atom2"], c["atom1"])
        if c["type"] in ("vdw", "hydrophobic") and c["distance"] > vdw_max:
            continue
        ch, num, rn, an = a_p.split(":")[:4]
        lname = a_l.split(":")[-1]
        r = res.setdefault((ch, num, rn), {"contacts": [], "xyz": []})
        r["contacts"].append({"lig": lname, "atom": an, "d": c["distance"], "type": c["type"], "verdict": c["verdict"]})
    if detail == "simple":
        # a simpler picture: a water only when it makes a hydrogen bond to the ligand; a residue only with a polar contact, or with two
        # or more contacts of which the closest is within 3.8 A. The full set stays in the tables of the report.
        keep = {}
        for key, r in res.items():
            if key[2] in chem.WATER_NAMES:
                ok = any(c["type"] == "hbond" for c in r["contacts"])
            else:
                ok = any(c["type"] in ("hbond", "salt_bridge", "metal") for c in r["contacts"]) or \
                     (len(r["contacts"]) >= 2 and min(c["d"] for c in r["contacts"]) <= 3.8)
            if ok:
                keep[key] = r
        res = keep
    name_to_local = {names[j]: j for j in range(len(names))}
    # atom coordinates of the contact atoms of each residue
    index = {(m.residues[m.res_idx[k]]["chain"], str(m.residues[m.res_idx[k]]["resnum"]), m.atom_name[k]): k for k in range(m.n_atoms)
             if m.residues[m.res_idx[k]]["kind"] == "polymer"}
    nodes = []
    for (ch, num, rn), r in res.items():
        pts = [m.xyz[index[(ch, num, c["atom"])]] for c in r["contacts"] if (ch, num, c["atom"]) in index]
        lig_pts = np.array([D[name_to_local[c["lig"]]] for c in r["contacts"]])
        lig_mean = lig_pts.mean(axis=0)
        outward = lig_mean / (np.linalg.norm(lig_mean) + 1e-9) if np.linalg.norm(lig_mean) > 0.3 else np.array([1.0, 0.0])
        q = ((np.mean(pts, axis=0) - cen) @ frame[:2].T) if pts else lig_mean
        qd = q / (np.linalg.norm(q) + 1e-9)
        w = 0.35 if np.linalg.norm(q) > 3.0 else 0.0
        direction = (1 - w) * outward + w * qd
        direction = direction / (np.linalg.norm(direction) + 1e-9)
        ideal = (lig_mean + 3.0 * direction) * PX
        has_polar = any(c["type"] in ("hbond", "salt_bridge", "metal") for c in r["contacts"])
        nodes.append({"ch": ch, "num": num, "rn": rn, "contacts": r["contacts"], "ideal": ideal, "pos": ideal.copy(),
                      "cls": residue_class(rn), "polar": has_polar,
                      "plddt": float(np.mean([m.bfactor[k] for k in np.where((m.res_idx == m.res_idx[index[(ch, num, r["contacts"][0]["atom"])]]))[0]]))
                      if (ch, num, r["contacts"][0]["atom"]) in index else None})
    relax(nodes, list(atom_pos.values()))
    return {"mol": mol, "atom_pos": atom_pos, "local": local, "names": names, "elems": elems, "nodes": nodes,
            "li": li, "name_to_local": name_to_local, "D": D * PX, "ligand": ligand_unit,
            "resname": m.residues[m.res_idx[li[0]]]["resname"], "predicted": bool(result["flags"].get("predicted_model")),
            "path": path, "model_index": model_index,
            "bond_source": mol.GetProp("bond_source") if mol.HasProp("bond_source") else "geometry guess"}


def relax(nodes, ligand_pts, iters=600):
    if not nodes:
        return
    rng = np.random.default_rng(1)
    pos = np.array([n["pos"] for n in nodes]) + rng.normal(0, 1.0, (len(nodes), 2))   # break exact ties
    ideal = np.array([n["ideal"] for n in nodes])
    L = np.array(ligand_pts)
    for t in range(iters):
        f = 0.06 * (ideal - pos)
        d = pos[:, None, :] - pos[None, :, :]
        sd = np.linalg.norm(d / np.array([100.0, 38.0]), axis=2) + np.eye(len(pos)) * 1e6 + 1e-6   # cards are wide and flat
        push = np.clip(1.0 - sd, 0, None)[:, :, None] * d / sd[:, :, None] * 60
        f += 0.3 * push.sum(axis=1)
        dl = pos[:, None, :] - L[None, :, :]
        dd = np.linalg.norm(dl / np.array([58.0, 30.0]), axis=2) + 1e-6          # a card must clear every ligand atom
        f += 0.25 * (np.clip(1.0 - dd, 0, None)[:, :, None] * dl / dd[:, :, None] * 9).sum(axis=1)
        pos += np.clip(f, -6, 6) * (1.0 - 0.8 * t / iters)      # capped, cooling steps keep the layout from flying apart
    for n, p in zip(nodes, pos):
        n["pos"] = p


def svg(d, title=None):
    nodes, ap, mol = d["nodes"], d["atom_pos"], d["mol"]
    li = list(d["li"])
    pts = list(ap.values()) + [n["pos"] for n in nodes]
    P = np.array(pts)
    pad = 70
    x0, y0 = P.min(axis=0) - pad
    x1, y1 = P.max(axis=0) + pad
    W, H = x1 - x0, y1 - y0
    if W < 860:                      # keep the legend row inside the picture, centre the drawing
        x0 -= (860 - W) / 2
        W = 860
    sx = lambda p: p[0] - x0
    sy = lambda p: -(p[1]) - (-y1)           # SVG y axis points down, flip so the picture is not mirrored
    X = lambda p: (p[0] - x0, y1 - p[1])
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.0f} {H + 70:.0f}" width="100%" style="max-width:{max(W, 520):.0f}px" role="img" '
         f'aria-label="Interaction diagram of {d["resname"]}"><style>text{{font-family:-apple-system,Helvetica,Arial,sans-serif}}'
         f'.bond{{stroke:var(--fg,#222);stroke-width:2;stroke-linecap:round}}.lab{{font-size:15px;font-weight:600;text-anchor:middle;dominant-baseline:central}}'
         f'.rl{{font-size:12.5px;text-anchor:middle;fill:#1d2326}}.rs{{font-size:10.5px;text-anchor:middle;fill:#444}}.dst{{font-size:11px;text-anchor:middle;fill:#2f855a;font-weight:700}}.halo2{{fill:#fff;stroke:#fff;stroke-width:4px}}'
         f'.hb{{stroke:#2f855a;stroke-width:2.2;stroke-dasharray:6 4;fill:none}}.vd{{stroke:#c1824a;stroke-width:1.1;stroke-dasharray:1.5 3;fill:none;opacity:.8}}'
         f'.sb{{stroke:#8e44ad;stroke-width:2.2;stroke-dasharray:6 4;fill:none}}</style>']
    # contact lines first
    for n in nodes:
        px, py = X(n["pos"])
        shown = 0
        for c in sorted(n["contacts"], key=lambda c: (c["type"] not in ("hbond", "salt_bridge", "metal"), c["d"])):
            gl = d["name_to_local"].get(c["lig"])
            if gl is None:
                continue
            ax, ay = X(ap[int(d["li"][gl])])
            polar = c["type"] in ("hbond", "salt_bridge", "metal")
            if not polar and shown >= 1:
                continue
            shown += 0 if polar else 1
            v = np.array([ax - px, ay - py])
            ln = np.linalg.norm(v) + 1e-9
            sx0, sy0 = px + v[0] / ln * 26, py + v[1] / ln * 20
            ex, ey = ax - v[0] / ln * (13 if polar else 6), ay - v[1] / ln * (13 if polar else 6)
            cls = "sb" if c["type"] == "salt_bridge" else "hb" if polar else "vd"
            o.append(f'<line class="{cls}" x1="{sx0:.1f}" y1="{sy0:.1f}" x2="{ex:.1f}" y2="{ey:.1f}"><title>{c["lig"]} to {n["rn"]}{n["num"]} {c["atom"]} {c["d"]} A ({c["type"]})</title></line>')
            if polar:
                mx, my = (sx0 + ex) / 2, (sy0 + ey) / 2
                o.append(f'<text class="dst halo2" x="{mx:.1f}" y="{my - 4:.1f}">{c["d"]:.2f}</text>'
                         f'<text class="dst" x="{mx:.1f}" y="{my - 4:.1f}">{c["d"]:.2f}</text>')       # outline layer, then number
    # ligand bonds
    conf = {int(k): X(ap[int(k)]) for k in li}
    for b in kekule(mol).GetBonds():             # aromatic rings are drawn as alternating single and double bonds
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        (x1_, y1_), (x2_, y2_) = conf[int(li[i])], conf[int(li[j])]
        v = np.array([x2_ - x1_, y2_ - y1_])
        L = np.linalg.norm(v) + 1e-9
        nrm = np.array([-v[1], v[0]]) / L
        shr = lambda k: 11 if d["elems"][k].upper() != "C" else 0
        a = np.array([x1_, y1_]) + v / L * shr(i)
        c = np.array([x2_, y2_]) - v / L * shr(j)
        order = b.GetBondTypeAsDouble()
        if order == 2:
            off = 3.0
            o.append(f'<line class="bond" x1="{a[0] + nrm[0] * off:.1f}" y1="{a[1] + nrm[1] * off:.1f}" x2="{c[0] + nrm[0] * off:.1f}" y2="{c[1] + nrm[1] * off:.1f}"/>')
            o.append(f'<line class="bond" x1="{a[0] - nrm[0] * off:.1f}" y1="{a[1] - nrm[1] * off:.1f}" x2="{c[0] - nrm[0] * off:.1f}" y2="{c[1] - nrm[1] * off:.1f}"/>')
        else:
            o.append(f'<line class="bond" x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{c[0]:.1f}" y2="{c[1]:.1f}"/>')
    for k, (e, nm) in enumerate(zip(d["elems"], d["names"])):
        if e.upper() != "C":
            x, y = conf[int(li[k])]
            o.append(f'<text class="lab" x="{x:.1f}" y="{y:.1f}" fill="{ELEM_COL.get(e.upper(), "#333")}"><title>{nm}</title>{e.capitalize()}</text>')
    # residues
    for n in nodes:
        x, y = X(n["pos"])
        edge = "#2f855a" if n["polar"] else "#9aa5ad"
        o.append(f'<ellipse cx="{x:.1f}" cy="{y:.1f}" rx="31" ry="21" fill="{FILL[n["cls"]]}" stroke="{edge}" stroke-width="{2.4 if n["polar"] else 1.2}"/>'
                 f'<text class="rl" x="{x:.1f}" y="{y - 1:.1f}"><tspan font-weight="600">{n["rn"].title() if n["rn"] in STD_AA else n["rn"]}{n["num"]}</tspan></text>'
                 f'<text class="rs" x="{x:.1f}" y="{y + 12:.1f}">{n["ch"]}</text>')
    # legend
    ly = H + 12
    items = [("hydrophobic", "hydrophobic"), ("polar", "polar"), ("positive", "positive"), ("negative", "negative")]
    lx = 10
    for k, t in items:
        o.append(f'<ellipse cx="{lx + 10:.0f}" cy="{ly + 8:.0f}" rx="10" ry="7" fill="{FILL[k]}" stroke="#9aa5ad"/><text x="{lx + 25:.0f}" y="{ly + 12:.0f}" font-size="11" fill="currentColor">{t}</text>')
        lx += 108
    o.append(f'<line class="hb" x1="{lx:.0f}" y1="{ly + 8:.0f}" x2="{lx + 30:.0f}" y2="{ly + 8:.0f}"/><text x="{lx + 36:.0f}" y="{ly + 12:.0f}" font-size="11" fill="currentColor">H-bond (A)</text>')
    lx += 112
    o.append(f'<line class="vd" x1="{lx:.0f}" y1="{ly + 8:.0f}" x2="{lx + 30:.0f}" y2="{ly + 8:.0f}"/><text x="{lx + 36:.0f}" y="{ly + 12:.0f}" font-size="11" fill="currentColor">other contact (up to 4 A)</text>')
    o.append("</svg>")
    return "".join(o)


# ---------------------------------------------------------------- modern style
TINT = {"hydrophobic": ("#b7791f", "rgba(214,158,46,.20)"), "polar": ("#0f8a8a", "rgba(15,138,138,.16)"),
        "positive": ("#5b3fc4", "rgba(91,63,196,.16)"), "negative": ("#c53030", "rgba(197,48,48,.16)"),
        "other": ("#8a949b", "rgba(138,148,155,.18)")}
CLS_LETTER = {"hydrophobic": "●", "polar": "◆", "positive": "+", "negative": "−", "other": "○"}
CLS_WORD = {"hydrophobic": "hydrophobic", "polar": "polar", "positive": "positive +", "negative": "negative \u2212", "other": "other"}


def svg_modern(d):
    conf_label = "pLDDT" if d.get("predicted") else "B-factor"
    """Card-like diagram: rounded residue chips, curved edges with distance pills, hover focus, theme aware."""
    nodes, ap, mol = d["nodes"], d["atom_pos"], d["mol"]
    li = list(d["li"])
    P = np.array(list(ap.values()) + [n["pos"] for n in nodes])
    pad = 60
    x0, y0 = P.min(axis=0) - pad
    x1, y1 = P.max(axis=0) + pad
    W, H = x1 - x0, y1 - y0
    if W < 780:                      # centre the picture so the legend row fits without a blank side
        x0 -= (780 - W) / 2
        W = 780
    X = lambda p: (p[0] - x0, y1 - p[1])
    conf = {int(k): X(ap[int(k)]) for k in li}
    css = ("text{font-family:-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif}"
           ".bd{stroke:var(--fg,#222);stroke-width:2.6;stroke-linecap:round;fill:none}"
           ".al{font-size:14px;font-weight:700;text-anchor:middle;dominant-baseline:central}"
           ".chip{stroke-width:1.2}.rn{font-size:13px;font-weight:700;text-anchor:start;fill:var(--fg,#1d2326)}"
           ".mrs{font-size:10px;fill:var(--mut,#667078);text-anchor:start}.cl{font-size:12px;font-weight:700;text-anchor:middle;dominant-baseline:central}"
           ".mhb{stroke:#2f855a;stroke-width:2.6;fill:none;stroke-linecap:round}.msb{stroke:#8e44ad;stroke-width:2.6;fill:none;stroke-linecap:round}"
           ".ot{stroke:var(--mut,#8a949b);stroke-width:1.4;fill:none;opacity:.55;stroke-linecap:round}"
           ".pill{fill:var(--card,#fff);stroke:#2f855a;stroke-width:1.2}.pt{font-size:11px;font-weight:700;fill:#2f855a;text-anchor:middle;dominant-baseline:central}"
           ".res{transition:opacity .15s;cursor:default}svg:has(.res:hover) .res:not(:hover){opacity:.18}.res:hover .chip{stroke-width:2.4}"
           ".halo{fill:var(--acc,#2b6cb0);opacity:.07}.bur{fill:#c0392b;opacity:.0}")
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.0f} {H + 96:.0f}" width="100%" style="max-width:{W:.0f}px" role="img" '
         f'aria-label="Interaction diagram of {d["resname"]}"><style>{css}</style>']
    # soft pocket halo behind the ligand
    for k in li:
        x, y = conf[int(k)]
        o.append(f'<circle class="halo" cx="{x:.1f}" cy="{y:.1f}" r="40"/>')
    # per residue group: edges, then chip
    groups = []
    for n in nodes:
        px, py = X(n["pos"])
        col, fill = TINT[n["cls"]]
        g = [f'<g class="res"><title>{n["rn"]}{n["num"]} chain {n["ch"]}, {len(n["contacts"])} contacts, closest {min(c["d"] for c in n["contacts"]):.2f} A'+(f', {conf_label} {n["plddt"]:.0f}' if n.get("plddt") is not None else '')+'</title>']
        shown = 0
        for c in sorted(n["contacts"], key=lambda c: (c["type"] not in ("hbond", "salt_bridge", "metal"), c["d"])):
            gl = d["name_to_local"].get(c["lig"])
            if gl is None:
                continue
            polar = c["type"] in ("hbond", "salt_bridge", "metal")
            if not polar and shown >= 1:
                continue
            shown += 0 if polar else 1
            ax, ay = X(ap[int(d["li"][gl])])
            v = np.array([ax - px, ay - py])
            ln = np.linalg.norm(v) + 1e-9
            s = np.array([px, py]) + v / ln * min(34, ln * 0.4)
            e = np.array([ax, ay]) - v / ln * (14 if polar else 7)
            nrm = np.array([-v[1], v[0]]) / ln
            bend = 0.12 * ln
            mid = (s + e) / 2 + nrm * bend
            cls = "msb" if c["type"] == "salt_bridge" else "mhb" if polar else "ot"
            g.append(f'<path class="{cls}" d="M{s[0]:.1f},{s[1]:.1f} Q{mid[0]:.1f},{mid[1]:.1f} {e[0]:.1f},{e[1]:.1f}"/>')
            if polar:
                t = 0.5
                bx = (1 - t) ** 2 * s[0] + 2 * (1 - t) * t * mid[0] + t ** 2 * e[0]
                by = (1 - t) ** 2 * s[1] + 2 * (1 - t) * t * mid[1] + t ** 2 * e[1]
                g.append(f'<rect class="pill" x="{bx - 19:.1f}" y="{by - 9:.1f}" width="38" height="18" rx="9"/><text class="pt" x="{bx:.1f}" y="{by:.1f}">{c["d"]:.2f}</text>')
        if n["rn"] == "HOH":
            g.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="13" fill="rgba(80,150,220,.18)" stroke="#4a90d9" stroke-width="1.4"/>'
                     f'<text class="mrs" x="{px:.1f}" y="{py + 0.5:.1f}" style="text-anchor:middle;dominant-baseline:central;font-weight:700;fill:#2b6cb0">{n["num"]}</text>')
        else:
            nm = n["rn"].title() if n["rn"] in STD_AA else n["rn"]
            w = 86
            g.append(f'<rect class="chip" x="{px - w / 2:.1f}" y="{py - 14:.1f}" width="{w}" height="28" rx="9" fill="{fill}" stroke="{col}"/>')
            g.append(f'<circle cx="{px - w / 2 + 14:.1f}" cy="{py:.1f}" r="6" fill="{col}"/>')
            g.append(f'<text class="rn" x="{px - w / 2 + 26:.1f}" y="{py + 0.5:.1f}" style="dominant-baseline:central">{nm}{n["num"]}<tspan style="font-weight:400;fill:var(--mut,#667078)"> {n["ch"]}</tspan></text>')
        g.append("</g>")
        groups.append("".join(g))
    o.extend(groups)
    # ligand on top, ball-and-stick with a light from the upper left
    km = kekule(mol)
    ecol = lambda e: ELEM_COL.get(e.upper(), "#59636b")
    defs = ['<defs><filter id="sh" x="-20%" y="-20%" width="140%" height="140%"><feDropShadow dx="1.5" dy="2.5" stdDeviation="2.2" flood-color="#000" flood-opacity=".28"/></filter>']
    for e in sorted({x.upper() for x in d["elems"]}):
        c = "#6b7279" if e == "C" else ecol(e)
        defs.append(f'<radialGradient id="g{e}" cx="35%" cy="30%" r="75%"><stop offset="0" stop-color="#fff" stop-opacity=".95"/>'
                    f'<stop offset=".28" stop-color="{c}"/><stop offset="1" stop-color="{c}"/></radialGradient>'
                    f'<radialGradient id="s{e}" cx="35%" cy="30%" r="80%"><stop offset="0" stop-color="{c}" stop-opacity=".55"/><stop offset="1" stop-color="#111" stop-opacity=".75"/></radialGradient>')
    defs.append("</defs>")
    o.insert(1, "".join(defs))
    o.append('<g filter="url(#sh)">')

    gid = [0]

    def tube(p1, p2, w):
        """One stroke shaded across its width like a cylinder, lit from the upper left (no separate highlight line)."""
        x1_, y1_ = p1
        x2_, y2_ = p2
        v = np.array([x2_ - x1_, y2_ - y1_])
        n = np.array([-v[1], v[0]]) / (np.linalg.norm(v) + 1e-9)
        if n[0] + n[1] > 0:                       # point the normal toward the light (up and left)
            n = -n
        mx, my = (x1_ + x2_) / 2, (y1_ + y2_) / 2
        gid[0] += 1
        o.append(f'<linearGradient id="t{gid[0]}" gradientUnits="userSpaceOnUse" x1="{mx + n[0] * w / 2:.1f}" y1="{my + n[1] * w / 2:.1f}" '
                 f'x2="{mx - n[0] * w / 2:.1f}" y2="{my - n[1] * w / 2:.1f}"><stop offset="0" stop-color="#39414a"/><stop offset=".3" stop-color="#cfd5da"/>'
                 f'<stop offset=".55" stop-color="#8b949c"/><stop offset="1" stop-color="#2b3137"/></linearGradient>')
        o.append(f'<path d="M{x1_:.1f},{y1_:.1f} L{x2_:.1f},{y2_:.1f}" stroke="url(#t{gid[0]})" stroke-width="{w:.1f}" stroke-linecap="butt" fill="none"/>')
    for b in km.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        (xa, ya), (xb, yb) = conf[int(li[i])], conf[int(li[j])]
        v = np.array([xb - xa, yb - ya])
        Ln = np.linalg.norm(v) + 1e-9
        nrm = np.array([-v[1], v[0]]) / Ln
        pa, pb = np.array([xa, ya]), np.array([xb, yb])
        if b.GetBondType() == Chem.BondType.DOUBLE:
            for sgn in (1, -1):
                q = nrm * 4.4 * sgn
                tube(pa + q, pb + q, 5.0)
        else:
            tube(pa, pb, 7.4)
    for k, (e, nm) in enumerate(zip(d["elems"], d["names"])):
        x, y = conf[int(li[k])]
        E = e.upper()
        if E == "C":
            o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="6.2" fill="url(#gC)"><title>{nm}</title></circle>')
        else:
            r_ = 13 if E in ("S", "CL", "P", "BR") else 11.5
            o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r_}" fill="url(#g{E})"><title>{nm}</title></circle>'
                     f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r_}" fill="none" stroke="#000" stroke-opacity=".25" stroke-width="1"/>'
                     f'<text class="al" x="{x:.1f}" y="{y + 0.5:.1f}" fill="#fff" style="font-size:12.5px">{e.capitalize()}</text>')
    o.append("</g>")
    # legend, two rows: residue types (with water) and contact types
    ly = H + 14
    lx = 10
    o.append(f'<text x="{lx}" y="{ly + 13}" font-size="11" font-weight="700" fill="var(--mut,#667078)">Residues</text>')
    lx += 72
    for k in ("hydrophobic", "polar", "positive", "negative"):
        col, fill = TINT[k]
        o.append(f'<rect x="{lx}" y="{ly}" width="34" height="18" rx="6" fill="{fill}" stroke="{col}" stroke-width="1.4"/><circle cx="{lx + 11}" cy="{ly + 9}" r="5" fill="{col}"/><text x="{lx + 40}" y="{ly + 13}" font-size="11.5" fill="currentColor">{k}</text>')
        lx += 118
    o.append(f'<circle cx="{lx + 9}" cy="{ly + 9}" r="9" fill="rgba(80,150,220,.18)" stroke="#4a90d9"/><text x="{lx + 24}" y="{ly + 13}" font-size="11.5" fill="currentColor">water</text>')
    ly += 28
    lx = 10
    o.append(f'<text x="{lx}" y="{ly + 13}" font-size="11" font-weight="700" fill="var(--mut,#667078)">Contacts</text>')
    lx += 72
    o.append(f'<path class="mhb" d="M{lx},{ly + 9} L{lx + 30},{ly + 9}"/><text x="{lx + 38}" y="{ly + 13}" font-size="11.5" fill="currentColor">hydrogen bond, distance in A</text>')
    lx += 215
    o.append(f'<path class="ot" d="M{lx},{ly + 9} L{lx + 30},{ly + 9}"/><text x="{lx + 38}" y="{ly + 13}" font-size="11.5" fill="currentColor">closest other contact (under 4 A)</text>')
    o.append(f'<text x="10" y="{ly + 40}" font-size="10.5" fill="var(--mut,#667078)">Residue label = amino acid, number, chain. Hover a residue to isolate its contacts and see distance and {conf_label}.</text>')
    o.append("</svg>")
    return "".join(o)


# ---------------------------------------------------------------- figure files written by the program
def standalone_svg(svg_text):
    """The in-page SVG as a file that opens anywhere: theme variables resolved, white background, explicit size."""
    import re
    vb = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg_text)
    w, h = float(vb.group(1)), float(vb.group(2))
    t = re.sub(r"var\(--[a-z-]+,\s*([^)]+)\)", r"\1", svg_text).replace("currentColor", "#1d2326")
    t = re.sub(r'width="100%" style="[^"]*"', f'width="{w:.0f}" height="{h:.0f}"', t, count=1)
    t = re.sub(r"(<svg[^>]*>)", lambda m: m.group(1) + f'<rect x="0" y="0" width="{w:.0f}" height="{h:.0f}" fill="#ffffff"/>', t, count=1)
    return t


def write_figures(d, base, dpi=300, views=None):
    """Write BASE_modern.svg and BASE_classic.svg (or the views given as (name, function) pairs), and PNG and PDF of both when cairosvg is installed.
    Returns the list of files written and a note about anything that was skipped."""
    out, note = [], ""
    try:
        import cairosvg
    except Exception:
        cairosvg = None
        note = "PNG and PDF were not written because cairosvg is not installed (pip install cairosvg). The SVG files and the download buttons in the HTML page still work."
    for view, fn in (views or (("modern", svg_modern), ("classic", svg))):
        t = standalone_svg(fn(d))
        path = f"{base}_{view}.svg"
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
        out.append(path)
        if cairosvg:
            cairosvg.svg2png(bytestring=t.encode("utf-8"), write_to=f"{base}_{view}.png", dpi=dpi, scale=dpi / 96.0)
            cairosvg.svg2pdf(bytestring=t.encode("utf-8"), write_to=f"{base}_{view}.pdf")
            out += [f"{base}_{view}.png", f"{base}_{view}.pdf"]
    return out, note
