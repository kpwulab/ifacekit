"""Detailed ligand interaction diagram: residues that hydrogen-bond or salt-bridge the ligand are drawn as chemical structures.

The layout of the ordinary diagram is reused (ligand in 2D, residues on the side where they are in 3D). Here each such
residue is drawn as its side chain (plus the backbone atoms that take part in a contact), turned so that the contacting
atom faces the ligand. The residues are pushed apart by their real size. Residues with only van der Waals or hydrophobic
contacts, and waters, stay small chips so the picture does not fill up.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor

from .ligplot import ELEM_COL, PX, STD_AA, TINT, kekule

BOND_LEN = 30.0                    # pixels per bond in a residue drawing (ligand bonds are about 66 px)
_ORDER = {"-": Chem.BondType.SINGLE, "=": Chem.BondType.DOUBLE}

# side chains from CA outwards; "-" single and "=" double, aromatic rings written as alternating bonds
SIDE = {
    "GLY": "", "ALA": "CA-CB", "SER": "CA-CB CB-OG", "CYS": "CA-CB CB-SG", "THR": "CA-CB CB-OG1 CB-CG2",
    "VAL": "CA-CB CB-CG1 CB-CG2", "LEU": "CA-CB CB-CG CG-CD1 CG-CD2", "ILE": "CA-CB CB-CG1 CB-CG2 CG1-CD1",
    "MET": "CA-CB CB-CG CG-SD SD-CE", "PRO": "CA-CB CB-CG CG-CD CD-N N-CA",
    "ASP": "CA-CB CB-CG CG=OD1 CG-OD2", "GLU": "CA-CB CB-CG CG-CD CD=OE1 CD-OE2",
    "ASN": "CA-CB CB-CG CG=OD1 CG-ND2", "GLN": "CA-CB CB-CG CG-CD CD=OE1 CD-NE2",
    "LYS": "CA-CB CB-CG CG-CD CD-CE CE-NZ", "ARG": "CA-CB CB-CG CG-CD CD-NE NE-CZ CZ-NH1 CZ=NH2",
    "HIS": "CA-CB CB-CG CG-ND1 ND1-CE1 CE1=NE2 NE2-CD2 CD2=CG",
    "PHE": "CA-CB CB-CG CG=CD1 CD1-CE1 CE1=CZ CZ-CE2 CE2=CD2 CD2-CG",
    "TYR": "CA-CB CB-CG CG=CD1 CD1-CE1 CE1=CZ CZ-CE2 CE2=CD2 CD2-CG CZ-OH",
    "TRP": "CA-CB CB-CG CG=CD1 CD1-NE1 NE1-CE2 CE2=CD2 CD2-CG CD2-CE3 CE3=CZ3 CZ3-CH2 CH2=CZ2 CZ2-CE2",
}
BACKBONE_N = "N-CA"
BACKBONE_C = "CA-C C=O"
POLAR_TYPES = ("hbond", "salt_bridge", "metal")
_CACHE = {}


def _fragment(resname, contact_atoms):
    """2D drawing of CA + side chain (+ the backbone atoms in contact): names, bonds, coordinates centred on their mean."""
    extra = set()
    if "N" in contact_atoms:
        extra.add("N")
    if contact_atoms & {"C", "O", "OXT"}:
        extra.add("C")
    key = (resname, frozenset(extra))
    if key in _CACHE:
        return _CACHE[key]
    spec = SIDE.get(resname)
    if spec is None:
        return None
    parts = spec.split()
    if "N" in extra:
        parts.append(BACKBONE_N)
    if "C" in extra:
        parts += BACKBONE_C.split()
    names, bonds = [], []
    for p in parts:
        sep = "=" if "=" in p else "-"
        a, b = p.split(sep)
        for x in (a, b):
            if x not in names:
                names.append(x)
        if not any({a, b} == {u, v} for u, v, _ in bonds):
            bonds.append((a, b, sep))
    if not names:
        names = ["CA"]
    rw = Chem.RWMol()
    for n in names:
        rw.AddAtom(Chem.Atom({"C": 6, "N": 7, "O": 8, "S": 16}[n[0]]))
    for a, b, s in bonds:
        rw.AddBond(names.index(a), names.index(b), _ORDER[s])
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    xy = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(len(names))])
    if len(names) > 1:
        ln = np.mean([np.linalg.norm(xy[names.index(a)] - xy[names.index(b)]) for a, b, _ in bonds])
        xy = xy / ln * BOND_LEN
    xy = xy - xy.mean(axis=0)
    out = {"names": names, "bonds": bonds, "xy": xy}
    _CACHE[key] = out
    return out


def _turn(xy, ang, mirror):
    if mirror:
        xy = xy * np.array([-1.0, 1.0])
    R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    return xy @ R.T


def _placed(frag, contact_atoms, toward, pos, obstacles):
    """Fragment turned (and mirrored if that helps) so that the contacting atoms face the ligand and the atoms keep clear of
    the ligand and of what is already drawn. Tries turns of up to 100 degrees either side of the straight-on direction."""
    xy = frag["xy"]
    idx = [frag["names"].index(a) for a in contact_atoms if a in frag["names"]]
    if not idx or len(frag["names"]) < 2:
        return xy
    obs = np.array(obstacles)
    best = None
    for mirror in (False, True):
        base = _turn(xy, 0.0, mirror)
        c = base[idx].mean(axis=0)
        if np.linalg.norm(c) < 1.0:
            return xy
        a0 = np.arctan2(toward[1], toward[0]) - np.arctan2(c[1], c[0])
        for k in range(-7, 8):
            off = np.radians(15 * k)
            cand = _turn(xy, a0 + off, mirror)
            pts = pos + cand
            dd = np.linalg.norm(pts[:, None, :] - obs[None, :, :], axis=2)
            pen = np.clip(30.0 - dd, 0, None).sum() + 0.35 * abs(15 * k) + (3.0 if mirror else 0.0)
            ca = pos + cand[idx].mean(axis=0)
            pen += 0.5 * np.linalg.norm(ca - (pos + toward * np.linalg.norm(c)))      # contact atom on the ligand side
            if best is None or pen < best[0]:
                best = (pen, cand)
    return best[1]


def _crosses(p, q, a, b):
    def ccw(A, B, C):
        return (C[1] - A[1]) * (B[0] - A[0]) - (B[1] - A[1]) * (C[0] - A[0])
    return ccw(p, q, a) * ccw(p, q, b) < 0 and ccw(a, b, p) * ccw(a, b, q) < 0


def _seg_dist(pt, a, b):
    ab = b - a
    t = np.clip(np.dot(pt - a, ab) / (np.dot(ab, ab) + 1e-9), 0.0, 1.0)
    return float(np.linalg.norm(pt - (a + t * ab)))


def _line_cost(pos, targets, atoms, bonds):
    """Price of the lines from `pos` to the ligand atoms it contacts: a line that crosses a bond of the ligand, or runs close
    to a ligand atom other than its own, is hard to read (the situation of a contact seen from behind the ligand)."""
    cost = 0.0
    for idx, tp in targets:
        for i, j, pa, pb in bonds:
            if idx in (i, j):
                continue
            if _crosses(pos, tp, pa, pb):
                cost += 60.0
        for k, pk in enumerate(atoms):
            if k != idx and _seg_dist(pk, pos, tp) < 20.0:
                cost += 40.0
    return cost


def _choose_side(items, atoms, bonds, obstacles):
    """Move each picture round the ligand (keeping its distance) to the side where its contact lines stay clear of the
    rest of the ligand, as near as possible to the side it really is on in 3D."""
    obs = np.array(obstacles)
    for it in items:
        if not it.get("targets"):
            continue
        anchor = it["anchor"]
        v = it["ideal"] - anchor
        R = np.linalg.norm(v)
        a0 = np.arctan2(v[1], v[0])
        best = None
        for k in range(36):
            th = a0 + np.radians(10 * ((k + 1) // 2) * (1 if k % 2 else -1))
            pos = anchor + R * np.array([np.cos(th), np.sin(th)])
            cost = _line_cost(pos, it["targets"], atoms, bonds) + 0.5 * abs(np.degrees(th - a0))
            near = np.linalg.norm((obs - pos) / np.array([it["rx"] + 14.0, it["ry"] + 14.0]), axis=1)
            cost += 25.0 * np.count_nonzero(near < 1.0)
            if best is None or cost < best[0] - 1e-9:
                best = (cost, pos)
        it["ideal"] = best[1]


def _relax(items, ligand_pts, iters=1500):
    """Push the pictures apart by their size, keep them near their ideal place and away from every ligand atom."""
    pos = np.array([it["pos"] for it in items], dtype=float)
    rng = np.random.default_rng(3)
    pos += rng.normal(0, 1.0, pos.shape)
    ideal = np.array([it["ideal"] for it in items], dtype=float)
    rad = np.array([[it["rx"], it["ry"]] for it in items], dtype=float)
    L = np.array(ligand_pts)
    for t in range(iters):
        f = 0.02 * (ideal - pos)
        d = pos[:, None, :] - pos[None, :, :]
        lim = rad[:, None, :] + rad[None, :, :]
        sd = np.linalg.norm(d / lim, axis=2) + np.eye(len(pos)) * 1e6 + 1e-6
        f += 1.2 * (np.clip(1.0 - sd, 0, None)[:, :, None] * d / sd[:, :, None] * 60).sum(axis=1)
        dl = pos[:, None, :] - L[None, :, :]
        lim_l = rad[:, None, :] + np.array([16.0, 16.0])
        dd = np.linalg.norm(dl / lim_l, axis=2) + 1e-6
        f += 0.4 * (np.clip(1.0 - dd, 0, None)[:, :, None] * dl / dd[:, :, None] * 30).sum(axis=1)
        pos += np.clip(f, -7, 7) * (1.0 - 0.85 * t / iters)
    for it, p in zip(items, pos):
        it["pos"] = p


def svg_detail(d, all_residues=False):
    """The detailed diagram as an SVG string (same data `d` as the other diagrams, from ligplot.build).
    By default only residues that hydrogen-bond or salt-bridge the ligand are drawn as structures; all_residues=True draws every
    standard residue of the diagram that way."""
    nodes, ap, mol = d["nodes"], d["atom_pos"], d["mol"]
    li = list(d["li"])
    lig_px = {int(k): np.array(ap[int(k)]) for k in li}
    lig_pts = list(lig_px.values())
    obstacles = list(lig_pts)
    for b in mol.GetBonds():                      # points along every ligand bond, so that a box does not sit on a bond
        pa, pb = lig_px[int(li[b.GetBeginAtomIdx()])], lig_px[int(li[b.GetEndAtomIdx()])]
        obstacles += [pa + (pb - pa) * t for t in (0.33, 0.5, 0.67)]
    centre = np.mean(lig_pts, axis=0)
    items = []
    for n in nodes:
        polar = [c for c in n["contacts"] if c["type"] in POLAR_TYPES and d["name_to_local"].get(c["lig"]) is not None]
        full = (bool(polar) or all_residues) and n["rn"] in STD_AA and n["rn"] != "HOH"
        anchors = [lig_px[int(li[d["name_to_local"][c["lig"]]])] for c in (polar or n["contacts"]) if d["name_to_local"].get(c["lig"]) is not None]
        anchor = np.mean(anchors, axis=0) if anchors else centre
        v = np.array(n["ideal"]) - anchor
        direction = v / (np.linalg.norm(v) + 1e-9)
        closest = min(c["d"] for c in n["contacts"])
        it = {"node": n, "full": False, "anchor": anchor, "dir": direction}
        tgt = polar or sorted(n["contacts"], key=lambda c: c["d"])[:1]
        it["targets"] = [(d["name_to_local"][c["lig"]], lig_px[int(li[d["name_to_local"][c["lig"]]])])
                         for c in tgt if d["name_to_local"].get(c["lig"]) is not None]
        if full:
            drawn = polar or sorted(n["contacts"], key=lambda c: c["d"])[:1]       # the contacts that get a line
            cat = {c["atom"] for c in drawn}
            it["ring"] = {c["atom"] for c in polar}
            frag = _fragment(n["rn"], cat)
            if frag is not None and all(a in frag["names"] for a in cat):
                it["full"], it["frag"], it["cat"] = True, frag, cat
        if it["full"]:
            r = np.max(np.linalg.norm(it["frag"]["xy"], axis=1)) + 22 if len(it["frag"]["names"]) > 1 else 30
            it["rx"] = it["ry"] = r
            it["ideal"] = anchor + direction * (r + 48 + 6 * max(0.0, closest - 2.8))
        elif n["rn"] == "HOH":
            it["rx"] = it["ry"] = 16
            it["ideal"] = np.array(n["ideal"])
        else:
            it["rx"], it["ry"] = 52, 19
            it["ideal"] = np.array(n["ideal"])
        it["pos"] = it["ideal"].copy()
        items.append(it)
    if items:
        atoms = [lig_px[int(k)] for k in li]
        bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), lig_px[int(li[b.GetBeginAtomIdx()])], lig_px[int(li[b.GetEndAtomIdx()])])
                 for b in mol.GetBonds()]
        _choose_side(items, atoms, bonds, obstacles)
        for it in items:
            it["pos"] = it["ideal"].copy()
        _relax(items, obstacles)
    P = np.array(lig_pts + [it["pos"] + s * np.array([it["rx"] + 30, it["ry"] + 30]) for it in items for s in (1, -1)])
    x0, y0 = P.min(axis=0) - 40
    x1, y1 = P.max(axis=0) + 40
    W, H = x1 - x0, y1 - y0
    if W < 780:
        x0 -= (780 - W) / 2
        W = 780
    X = lambda p: (p[0] - x0, y1 - p[1])
    ecol = lambda e: ELEM_COL.get(e.upper(), "#59636b")
    css = ("text{font-family:-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif}"
           ".dd-b{stroke:var(--mut,#667078);stroke-width:2;stroke-linecap:round;fill:none}"
           ".dd-lb{stroke:var(--fg,#222);stroke-width:3;stroke-linecap:round;fill:none}"
           ".dd-al{font-size:12px;font-weight:700;text-anchor:middle;dominant-baseline:central}"
           ".dd-rn{font-size:13px;font-weight:700;fill:var(--fg,#1d2326);text-anchor:middle;dominant-baseline:central}"
           ".dd-mrs{font-size:10px;fill:var(--mut,#667078)}"
           ".dd-hb{stroke:#2f855a;stroke-width:2.4;stroke-dasharray:7 5;fill:none;stroke-linecap:round}"
           ".dd-sb{stroke:#8e44ad;stroke-width:2.4;stroke-dasharray:7 5;fill:none;stroke-linecap:round}"
           ".dd-ot{stroke:var(--mut,#667078);stroke-width:2.4;stroke-dasharray:3 4;fill:none;opacity:1;stroke-linecap:round}"
           ".dd-pill{fill:var(--card,#fff);stroke:#2f855a;stroke-width:1.2}.dd-pt{font-size:11px;font-weight:700;fill:#2f855a;text-anchor:middle;dominant-baseline:central}"
           ".dd-bg{fill:var(--card,#fff)}")
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.0f} {H + 84:.0f}" width="100%" style="max-width:{W:.0f}px" role="img" '
         f'aria-label="Detailed interaction diagram of {d["resname"]}"><style>{css}</style>']
    # contact lines first, so the structures lie on top of them
    labels = []
    taken = list(obstacles)                      # ligand points, then every chip and structure placed so far
    for it in items:
        if not it["full"]:
            taken += [it["pos"] + np.array([dx, 0.0]) for dx in (-30, -15, 0, 15, 30)]
    for it in sorted(items, key=lambda i: -len(i["frag"]["names"]) if i["full"] else 0):
        n = it["node"]
        pos = it["pos"]
        toward = it["anchor"] - pos
        toward = toward / (np.linalg.norm(toward) + 1e-9)
        if it["full"]:
            xy = _placed(it["frag"], it["cat"], toward, pos, taken)
            taken += [pos + q for q in xy]
            it["xy"] = xy
            it["atom_xy"] = {nm: pos + xy[i] for i, nm in enumerate(it["frag"]["names"])}
        shown = 0
        for c in sorted(n["contacts"], key=lambda c: (c["type"] not in POLAR_TYPES, c["d"])):
            gl = d["name_to_local"].get(c["lig"])
            if gl is None:
                continue
            is_polar = c["type"] in POLAR_TYPES
            if not is_polar and shown >= 1:
                continue
            shown += 0 if is_polar else 1
            end = lig_px[int(li[gl])]
            if it["full"] and c["atom"] in it["atom_xy"]:
                start = it["atom_xy"][c["atom"]]
                gap0 = 9
            else:
                start = pos
                gap0 = it["rx"] * 0.7 if n["rn"] != "HOH" else 16
            v = end - start
            ln = np.linalg.norm(v) + 1e-9
            s = start + v / ln * gap0
            e = end - v / ln * (14 if is_polar else 8)
            cls = "dd-sb" if c["type"] == "salt_bridge" else "dd-hb" if is_polar else "dd-ot"
            sx, sy = X(s)
            ex, ey = X(e)
            o.append(f'<path class="{cls}" d="M{sx:.1f},{sy:.1f} L{ex:.1f},{ey:.1f}"/>')
            if is_polar:
                labels.append(((sx + ex) / 2, (sy + ey) / 2, c["d"]))
    # residues
    for it in items:
        n = it["node"]
        px, py = X(it["pos"])
        title = f'{n["rn"]}{n["num"]} chain {n["ch"]}, closest {min(c["d"] for c in n["contacts"]):.2f} A'
        if it["full"]:
            frag, xy = it["frag"], it["xy"]
            col = TINT[n["cls"]][0]
            g = [f'<g><title>{title}</title>']
            for a, b, sep in frag["bonds"]:
                pa, pb = xy[frag["names"].index(a)], xy[frag["names"].index(b)]
                (xa, ya), (xb, yb) = X(it["pos"] + pa), X(it["pos"] + pb)
                v = np.array([xb - xa, yb - ya])
                nrm = np.array([-v[1], v[0]]) / (np.linalg.norm(v) + 1e-9)
                if sep == "=":
                    for sg in (1, -1):
                        q = nrm * 2.6 * sg
                        g.append(f'<path class="dd-b" d="M{xa + q[0]:.1f},{ya + q[1]:.1f} L{xb + q[0]:.1f},{yb + q[1]:.1f}"/>')
                else:
                    g.append(f'<path class="dd-b" d="M{xa:.1f},{ya:.1f} L{xb:.1f},{yb:.1f}"/>')
            for i, nm in enumerate(frag["names"]):
                x, y = X(it["pos"] + xy[i])
                if nm in it["ring"]:
                    g.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="11" fill="none" stroke="#2f855a" stroke-width="2" opacity=".55"/>')
                if nm == "CA":
                    g.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.2" fill="var(--mut,#8a949b)"/>')
                elif nm[0] != "C":
                    g.append(f'<circle class="dd-bg" cx="{x:.1f}" cy="{y:.1f}" r="8"/>'
                             f'<text class="dd-al" x="{x:.1f}" y="{y + 0.5:.1f}" fill="{ecol(nm[0])}">{nm[0]}</text>')
            # name of the residue on the side away from the ligand
            allp = it["pos"] + xy                          # beyond the atom that is farthest from the ligand
            far = int(np.argmax(np.linalg.norm(allp - it["anchor"], axis=1)))
            out_dir = allp[far] - it["anchor"]
            out_dir = out_dir / (np.linalg.norm(out_dir) + 1e-9)
            lab = allp[far] + out_dir * 20
            lx, ly = X(lab)
            g.append(f'<text class="dd-rn" x="{lx:.1f}" y="{ly:.1f}" fill="{col}">{n["rn"].title()}{n["num"]}'
                     f'<tspan class="dd-mrs" dx="4">{n["ch"]}</tspan></text></g>')
            o.append("".join(g))
        elif n["rn"] == "HOH":
            o.append(f'<g><title>{title}</title><circle cx="{px:.1f}" cy="{py:.1f}" r="14" fill="rgba(80,150,220,.18)" stroke="#4a90d9" stroke-width="1.4"/>'
                     f'<text class="dd-mrs" x="{px:.1f}" y="{py + 0.5:.1f}" style="text-anchor:middle;dominant-baseline:central;font-weight:700;fill:#2b6cb0">{n["num"]}</text></g>')
        else:
            col, fill = TINT[n["cls"]]
            nm = n["rn"].title() if n["rn"] in STD_AA else n["rn"]
            o.append(f'<g><title>{title}</title><rect x="{px - 40:.1f}" y="{py - 13:.1f}" width="80" height="26" rx="9" fill="{fill}" stroke="{col}" stroke-width="1.2"/>'
                     f'<text class="dd-rn" x="{px:.1f}" y="{py + 0.5:.1f}" style="font-size:12px">{nm}{n["num"]}<tspan class="dd-mrs" dx="3">{n["ch"]}</tspan></text></g>')
    # ligand, flat style so that it matches the residue drawings
    km = kekule(mol)
    for b in km.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        (xa, ya), (xb, yb) = X(lig_px[int(li[i])]), X(lig_px[int(li[j])])
        v = np.array([xb - xa, yb - ya])
        nrm = np.array([-v[1], v[0]]) / (np.linalg.norm(v) + 1e-9)
        if b.GetBondType() == Chem.BondType.DOUBLE:
            for sg in (1, -1):
                q = nrm * 3.2 * sg
                o.append(f'<path class="dd-lb" d="M{xa + q[0]:.1f},{ya + q[1]:.1f} L{xb + q[0]:.1f},{yb + q[1]:.1f}"/>')
        else:
            o.append(f'<path class="dd-lb" d="M{xa:.1f},{ya:.1f} L{xb:.1f},{yb:.1f}"/>')
    for k, (e, nm) in enumerate(zip(d["elems"], d["names"])):
        x, y = X(lig_px[int(li[k])])
        if e.upper() != "C":
            o.append(f'<circle class="dd-bg" cx="{x:.1f}" cy="{y:.1f}" r="11"/><text class="dd-al" x="{x:.1f}" y="{y + 0.5:.1f}" fill="{ecol(e)}" style="font-size:15px"><title>{nm}</title>{e.capitalize()}</text>')
    for x, y, dist in labels:
        o.append(f'<rect class="dd-pill" x="{x - 19:.1f}" y="{y - 9:.1f}" width="38" height="18" rx="9"/><text class="dd-pt" x="{x:.1f}" y="{y:.1f}">{dist:.2f}</text>')
    ly = H + 12
    o.append(f'<path class="dd-hb" d="M10,{ly + 8} L40,{ly + 8}"/><text x="48" y="{ly + 12}" font-size="11.5" fill="currentColor">hydrogen bond, distance in A</text>'
             f'<path class="dd-sb" d="M250,{ly + 8} L280,{ly + 8}"/><text x="288" y="{ly + 12}" font-size="11.5" fill="currentColor">salt bridge</text>'
             f'<path class="dd-ot" d="M390,{ly + 8} L420,{ly + 8}"/><text x="428" y="{ly + 12}" font-size="11.5" fill="currentColor">closest other contact</text>')
    first = ("Residues are drawn as structures (side chain; grey dot = C-alpha; ring = atom of a hydrogen bond or salt bridge)." if all_residues else
             "Residues that hydrogen-bond or salt-bridge the ligand are drawn as structures (side chain; grey dot = C-alpha; ring = contacting atom).")
    second = ("Waters and modified residues stay small. Residue label = amino acid, number, chain." if all_residues else
              "Other residues are small boxes. Residue label = amino acid, number, chain.")
    o.append(f'<text x="10" y="{ly + 34}" font-size="10.5" fill="var(--mut,#667078)">{first}</text>'
             f'<text x="10" y="{ly + 49}" font-size="10.5" fill="var(--mut,#667078)">{second}</text>')
    o.append("</svg>")
    return "".join(o)
