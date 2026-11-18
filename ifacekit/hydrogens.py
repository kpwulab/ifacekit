"""Place hydrogens on a protein model and check all-atom contacts.

Hydrogens are not read from the file (the loader drops them): they are built from the heavy atoms with ideal geometry, so that
every model is judged the same way, whatever the deposited file contains. Placement covers the 20 standard amino acids, the
N terminus and the backbone amide. Waters, ligands, nucleic acids and modified residues get no hydrogens.

Rotatable groups (methyl, NH3+, Ser/Thr/Cys/Tyr hydroxyl and thiol) are turned to the position with the best contacts, and the
amide groups of Asn and Gln and the imidazole of His are flipped and (His) protonated to the state with the best hydrogen
bonds and fewest overlaps, in a few passes over the whole model. The flip and protonation choices are made on the contacts of
the model as it is; no heavy atom is ever moved.

Radii follow the electron-cloud values of the Probe program used by MolProbity (aliphatic H 1.22, aromatic and polar H 1.05,
C 1.70, aromatic C 1.75, N 1.55, O 1.40, S 1.80); X-H distances are the shortened X-ray values (C-H 0.97, N-H 0.86,
O-H 0.84, S-H 1.20), the same Reduce uses for X-ray models. A clash is an overlap of 0.4 A or more between atoms that are more than four bonds apart; a polar hydrogen
next to an acceptor is a hydrogen bond, not a clash, up to an overlap of 1.0 A (1.2 A if one of the two atoms is charged); the limits were
calibrated against the clash lists of MolProbity.
"""
import numpy as np
from scipy.spatial import cKDTree

from . import chem

XH = {"C": 0.97, "N": 0.86, "O": 0.84, "S": 1.20}
TET = np.radians(109.5)
OVERLAP_CLASH = 0.40
HB_CAP = 1.00                          # calibrated against MolProbity (see the design document)
HB_CAP_CHARGED = 1.20
HEAVY_POLAR_CAP = 0.55                 # N/O pairs closer than this overlap are short hydrogen bonds (no hydrogens are involved, so no angle test)
FLIP_COST = 4.0                        # energy a flip must gain before it is accepted (calibrated against Reduce)
HH_HOPS = 3                            # two hydrogens are not compared when their parent atoms are this many bonds apart or fewer
H_ALI, H_POLAR = 1.22, 1.05            # aliphatic hydrogen; aromatic and polar hydrogen
HEAVY_RADIUS = {"C": 1.70, "N": 1.55, "O": 1.40, "S": 1.80, "SE": 1.90, "P": 1.80}
AROMATIC_C = {"PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}, "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
              "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"}, "HIS": {"CG", "CD2", "CE1"}}

# kind, parent, references, hydrogen names
#   t1 one H on a tetrahedral C (three neighbours)   t2 two H on a CH2 (two neighbours)   p1 one H on an sp2 atom (two neighbours)
#   m3 three H on a methyl or NH3 (neighbour, reference)   o1 one H on a hydroxyl or thiol (neighbour, reference)
#   n2 two H on an sp2 NH2 (neighbour, in-plane reference)
_CH2 = lambda x, a, b, n: ("t2", x, (a, b), n)
SPEC = {
    "ALA": [("m3", "CB", ("CA", "N"), ("HB1", "HB2", "HB3"))],
    "ARG": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), _CH2("CG", "CB", "CD", ("HG2", "HG3")), _CH2("CD", "CG", "NE", ("HD2", "HD3")),
            ("p1", "NE", ("CD", "CZ"), ("HE",)), ("n2", "NH1", ("CZ", "NE"), ("HH11", "HH12")), ("n2", "NH2", ("CZ", "NE"), ("HH21", "HH22"))],
    "ASN": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), ("n2", "ND2", ("CG", "CB"), ("HD21", "HD22"))],
    "ASP": [_CH2("CB", "CA", "CG", ("HB2", "HB3"))],
    "CYS": [_CH2("CB", "CA", "SG", ("HB2", "HB3")), ("o1", "SG", ("CB", "CA"), ("HG",))],
    "GLN": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), _CH2("CG", "CB", "CD", ("HG2", "HG3")), ("n2", "NE2", ("CD", "CG"), ("HE21", "HE22"))],
    "GLU": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), _CH2("CG", "CB", "CD", ("HG2", "HG3"))],
    "GLY": [],
    "HIS": [_CH2("CB", "CA", "CG", ("HB2", "HB3"))],
    "ILE": [("t1", "CB", ("CA", "CG1", "CG2"), ("HB",)), _CH2("CG1", "CB", "CD1", ("HG12", "HG13")),
            ("m3", "CG2", ("CB", "CA"), ("HG21", "HG22", "HG23")), ("m3", "CD1", ("CG1", "CB"), ("HD11", "HD12", "HD13"))],
    "LEU": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), ("t1", "CG", ("CB", "CD1", "CD2"), ("HG",)),
            ("m3", "CD1", ("CG", "CB"), ("HD11", "HD12", "HD13")), ("m3", "CD2", ("CG", "CB"), ("HD21", "HD22", "HD23"))],
    "LYS": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), _CH2("CG", "CB", "CD", ("HG2", "HG3")), _CH2("CD", "CG", "CE", ("HD2", "HD3")),
            _CH2("CE", "CD", "NZ", ("HE2", "HE3")), ("m3", "NZ", ("CE", "CD"), ("HZ1", "HZ2", "HZ3"))],
    "MET": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), _CH2("CG", "CB", "SD", ("HG2", "HG3")), ("m3", "CE", ("SD", "CG"), ("HE1", "HE2", "HE3"))],
    "PHE": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), ("p1", "CD1", ("CG", "CE1"), ("HD1",)), ("p1", "CD2", ("CG", "CE2"), ("HD2",)),
            ("p1", "CE1", ("CD1", "CZ"), ("HE1",)), ("p1", "CE2", ("CD2", "CZ"), ("HE2",)), ("p1", "CZ", ("CE1", "CE2"), ("HZ",))],
    "PRO": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), _CH2("CG", "CB", "CD", ("HG2", "HG3")), _CH2("CD", "CG", "N", ("HD2", "HD3"))],
    "SER": [_CH2("CB", "CA", "OG", ("HB2", "HB3")), ("o1", "OG", ("CB", "CA"), ("HG",))],
    "THR": [("t1", "CB", ("CA", "OG1", "CG2"), ("HB",)), ("o1", "OG1", ("CB", "CA"), ("HG1",)), ("m3", "CG2", ("CB", "CA"), ("HG21", "HG22", "HG23"))],
    "TRP": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), ("p1", "CD1", ("CG", "NE1"), ("HD1",)), ("p1", "NE1", ("CD1", "CE2"), ("HE1",)),
            ("p1", "CE3", ("CD2", "CZ3"), ("HE3",)), ("p1", "CZ2", ("CE2", "CH2"), ("HZ2",)), ("p1", "CZ3", ("CE3", "CH2"), ("HZ3",)),
            ("p1", "CH2", ("CZ2", "CZ3"), ("HH2",))],
    "TYR": [_CH2("CB", "CA", "CG", ("HB2", "HB3")), ("p1", "CD1", ("CG", "CE1"), ("HD1",)), ("p1", "CD2", ("CG", "CE2"), ("HD2",)),
            ("p1", "CE1", ("CD1", "CZ"), ("HE1",)), ("p1", "CE2", ("CD2", "CZ"), ("HE2",)), ("o1", "OH", ("CZ", "CE1"), ("HH",))],
    "VAL": [("t1", "CB", ("CA", "CG1", "CG2"), ("HB",)), ("m3", "CG1", ("CB", "CA"), ("HG11", "HG12", "HG13")),
            ("m3", "CG2", ("CB", "CA"), ("HG21", "HG22", "HG23"))],
}
_HA = ("t1", "CA", ("N", "C", "CB"), ("HA",))
_HA2 = ("t2", "CA", ("N", "C"), ("HA2", "HA3"))


def _unit(v):
    return v / np.linalg.norm(v)


def _p1(x, a, b, length):
    return x + length * _unit(-(_unit(a - x) + _unit(b - x)))


def _t1(x, a, b, c, length):
    return x + length * _unit(-(_unit(a - x) + _unit(b - x) + _unit(c - x)))


def _t2(x, a, b, length):
    ua, ub = _unit(a - x), _unit(b - x)
    d = _unit(-(ua + ub))
    n = _unit(np.cross(ua, ub))
    al = np.radians(54.75)
    return x + length * (d * np.cos(al) + n * np.sin(al)), x + length * (d * np.cos(al) - n * np.sin(al))


def _rotor(x, a, b, length, azimuth):
    """H on X bonded to a, at the tetrahedral angle, azimuth (radians) measured from the direction of b around the a-X axis."""
    u = _unit(x - a)
    w = (b - a) - np.dot(b - a, u) * u
    w = _unit(w) if np.linalg.norm(w) > 1e-6 else _unit(np.cross(u, np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])))
    v = np.cross(u, w)
    pol = np.pi - TET
    return x + length * (np.cos(pol) * u + np.sin(pol) * (np.cos(azimuth) * w + np.sin(azimuth) * v))


def _n2(x, a, b, length):
    u = _unit(a - x)
    w = (b - a) - np.dot(b - a, u) * u
    w = _unit(w)
    c, s = np.cos(np.radians(120)), np.sin(np.radians(120))
    return x + length * (c * u + s * w), x + length * (c * u - s * w)


class Hydrogens:
    """The hydrogens of one model, as flat arrays (index i of these arrays is atom n_heavy + i of the combined atom list)."""

    def __init__(self):
        self.xyz, self.parent, self.name, self.res = [], [], [], []
        self.polar, self.rad, self.charged = [], [], []
        self.group = []                  # index into groups, -1 for hydrogens that never move
        self.groups = []
        self.n = 0

    def add(self, xyz, parent, name, res, polar, rad, charged=False, group=-1):
        self.xyz.append(xyz)
        self.parent.append(parent)
        self.name.append(name)
        self.res.append(res)
        self.polar.append(polar)
        self.rad.append(rad)
        self.charged.append(charged)
        self.group.append(group)
        self.n += 1
        return self.n - 1

    def finish(self):
        self.xyz = np.array(self.xyz, dtype=float).reshape(-1, 3)
        for k in ("parent", "res", "group"):
            setattr(self, k, np.array(getattr(self, k), dtype=int))
        for k in ("polar", "charged"):
            setattr(self, k, np.array(getattr(self, k), dtype=bool))
        self.rad = np.array(self.rad, dtype=float)
        self.name = np.array(self.name, dtype=object)
        return self


def heavy_radii(m):
    """Probe-style radii of the heavy atoms."""
    rad = np.array([HEAVY_RADIUS.get(e, chem.vdw_radius(e)) for e in m.element], dtype=float)
    for k in range(m.n_atoms):
        r = m.residues[m.res_idx[k]]
        if m.element[k] == "C" and m.atom_name[k] in AROMATIC_C.get(r["resname"], ()):
            rad[k] = 1.75
    return rad


def _atom_map(m):
    amap = [dict() for _ in m.residues]
    for k in range(m.n_atoms):
        amap[m.res_idx[k]][m.atom_name[k]] = k
    return amap


def _hop_table(m, depth=3):
    """hops[i] = {j: bonds between heavy atoms i and j} for up to `depth` bonds."""
    adj = [[] for _ in range(m.n_atoms)]
    for a, b in m.bonds:
        adj[a].append(b)
        adj[b].append(a)
    hops = []
    for s in range(m.n_atoms):
        seen = {s: 0}
        front = [s]
        for d in range(1, depth + 1):
            nxt = []
            for u in front:
                for v in adj[u]:
                    if v not in seen:
                        seen[v] = d
                        nxt.append(v)
            front = nxt
        hops.append(seen)
    return hops


class _Rotor:
    """Methyl, NH3+, hydroxyl or thiol: the hydrogens turn around the bond X-a."""
    def __init__(self, res, x, a, b, hidx, length, step_deg, base, one):
        self.res, self.x, self.a, self.b, self.hidx, self.length, self.one = res, x, a, b, hidx, length, one
        self.cands = [np.radians(base + k * step_deg) for k in range(int(round((120 if not one else 360) / step_deg)))] if step_deg else [np.radians(base)]
        self.phi = self.cands[0]
        self.kind = "rot"


class _Flip:
    """Asn or Gln amide: the NH2 sits on ND2/NE2 (state 0) or on OD1/OE1 (state 1)."""
    def __init__(self, res, n_atom, o_atom, a, b, hidx, length):
        self.res, self.atoms, self.a, self.b, self.hidx, self.length = res, (n_atom, o_atom), a, b, hidx, length
        self.state = 0
        self.kind = "flip"


class _His:
    """Histidine ring: 2 orientations x 3 protonation states, four hydrogen slots (ND1, CE1, NE2, CD2 positions)."""
    def __init__(self, res, ring, hidx):
        self.res, self.ring, self.hidx = res, ring, hidx    # ring: atom indices CG, ND1, CD2, CE1, NE2
        self.state = (0, 1)                                 # (flipped, protonation: 0 HID, 1 HIE, 2 HIP)
        self.kind = "his"


# roles around the ring in order, with their neighbours: slot 0 = ND1, 1 = CE1, 2 = NE2, 3 = CD2
_HIS_NB = {0: ("CG", 1), 1: (0, 2), 2: (1, 3), 3: (2, "CG")}


def _his_roles(ring_xyz, flipped):
    """Position of each role (CG, 0 ND1, 1 CE1, 2 NE2, 3 CD2) and the index of the atom sitting there."""
    pos = {"CG": ring_xyz["CG"], 0: ring_xyz["ND1"], 1: ring_xyz["CE1"], 2: ring_xyz["NE2"], 3: ring_xyz["CD2"]}
    at = {"CG": "CG", 0: "ND1", 1: "CE1", 2: "NE2", 3: "CD2"}
    if flipped:
        pos[0], pos[3] = pos[3], pos[0]
        pos[1], pos[2] = pos[2], pos[1]
        at[0], at[3] = at[3], at[0]
        at[1], at[2] = at[2], at[1]
    return pos, at


def build(m):
    """Place all hydrogens with default rotor and flip states. Returns (H, groups)."""
    H = Hydrogens()
    amap = _atom_map(m)
    xyz = m.xyz
    n_res = len(m.residues)
    groups = []
    last = {}
    for ri, r in enumerate(m.residues):
        if not r["protein"] or r["resname"] not in SPEC:
            continue
        a = amap[ri]
        name = r["resname"]
        prev_c = None
        if (r["chain"],) in last:
            pr = last[(r["chain"],)]
            pa = amap[pr]
            if "C" in pa and "N" in a and np.linalg.norm(xyz[pa["C"]] - xyz[a["N"]]) < 2.0:
                prev_c = pa["C"]
        last[(r["chain"],)] = ri

        def has(*names):
            return all(nm in a for nm in names)

        def arom_rad(parent_name):
            return H_POLAR if parent_name in AROMATIC_C.get(name, ()) else H_ALI

        # backbone amide / N terminus
        if "N" in a and "CA" in a:
            if name == "PRO":
                if prev_c is None and "CD" in a:
                    h1, h2 = _t2(xyz[a["N"]], xyz[a["CA"]], xyz[a["CD"]], XH["N"])
                    for nm, p in (("H2", h1), ("H3", h2)):
                        H.add(p, a["N"], nm, ri, True, H_POLAR, True)
            elif prev_c is not None:
                H.add(_p1(xyz[a["N"]], xyz[prev_c], xyz[a["CA"]], XH["N"]), a["N"], "H", ri, True, H_POLAR)
            elif "C" in a:
                hid = [H.add(xyz[a["N"]], a["N"], f"H{k + 1}", ri, True, H_POLAR, True) for k in range(3)]
                g = _Rotor(ri, a["N"], a["CA"], a["C"], hid, XH["N"], 10, 60, False)
                groups.append(g)
                for h in hid:
                    H.group[h] = len(groups) - 1
        specs = list(SPEC[name])
        if "CA" in a:
            specs.append(_HA2 if name == "GLY" else _HA)
        for kind, parent, refs, hn in specs:
            if parent not in a or not all(rf in a for rf in refs):
                continue
            px = xyz[a[parent]]
            pel = m.element[a[parent]]
            L = XH.get(pel, 1.09)
            rf = [xyz[a[x]] for x in refs]
            polar = pel in ("N", "O", "S")
            rad = H_POLAR if polar else arom_rad(parent)
            chg = bool(m.charge[a[parent]] > 0) and polar
            if kind == "p1":
                H.add(_p1(px, rf[0], rf[1], L), a[parent], hn[0], ri, polar, rad, chg)
            elif kind == "t1":
                H.add(_t1(px, rf[0], rf[1], rf[2], L), a[parent], hn[0], ri, polar, rad, chg)
            elif kind == "t2":
                p1, p2 = _t2(px, rf[0], rf[1], L)
                H.add(p1, a[parent], hn[0], ri, polar, rad, chg)
                H.add(p2, a[parent], hn[1], ri, polar, rad, chg)
            elif kind == "n2":
                p1, p2 = _n2(px, rf[0], rf[1], L)
                hid = [H.add(p1, a[parent], hn[0], ri, True, H_POLAR, chg), H.add(p2, a[parent], hn[1], ri, True, H_POLAR, chg)]
                if name in ("ASN", "GLN"):
                    o_name = "OD1" if name == "ASN" else "OE1"
                    if o_name in a:
                        groups.append(_Flip(ri, a[parent], a[o_name], a[refs[0]], a[refs[1]], hid, L))
                        for h in hid:
                            H.group[h] = len(groups) - 1
            elif kind == "m3":
                hid = [H.add(px, a[parent], hn[k], ri, polar, rad, chg) for k in range(3)]
                g = _Rotor(ri, a[parent], a[refs[0]], a[refs[1]], hid, L, 10, 60, False)
                if pel == "C":                                  # methyl: staggered, never turned (Reduce does the same)
                    for h, pos in zip(hid, _rot_positions(m, g, np.radians(60))):
                        H.xyz[h] = pos
                    continue
                groups.append(g)
                for h in hid:
                    H.group[h] = len(groups) - 1
            elif kind == "o1":
                if name == "CYS" and any(m.element[k] == "S" and k != a["SG"] for k in _bonded(m, a["SG"])):
                    continue                                # disulfide: no thiol hydrogen
                hid = [H.add(px, a[parent], hn[0], ri, True, H_POLAR, False)]
                if name == "TYR":
                    g = _Rotor(ri, a[parent], a[refs[0]], a[refs[1]], hid, L, 180, 0, True)
                else:
                    g = _Rotor(ri, a[parent], a[refs[0]], a[refs[1]], hid, L, 15, 0, True)
                groups.append(g)
                H.group[hid[0]] = len(groups) - 1
        if name == "HIS" and has("CG", "ND1", "CD2", "CE1", "NE2"):
            slots = [H.add(xyz[a["ND1"]], a["ND1"], "HD1", ri, True, H_POLAR), H.add(xyz[a["CE1"]], a["CE1"], "HE1", ri, False, H_POLAR),
                     H.add(xyz[a["NE2"]], a["NE2"], "HE2", ri, True, H_POLAR), H.add(xyz[a["CD2"]], a["CD2"], "HD2", ri, False, H_POLAR)]
            groups.append(_His(ri, {k: a[k] for k in ("CG", "ND1", "CD2", "CE1", "NE2")}, slots))
            for h in slots:
                H.group[h] = len(groups) - 1
    H.finish()
    H.active = np.ones(H.n, bool)
    return H, groups


_bond_cache = {}


def _bonded(m, k):
    key = id(m)
    if key not in _bond_cache or _bond_cache[key][0] is not m:
        adj = {}
        for x, y in m.bonds:
            adj.setdefault(int(x), []).append(int(y))
            adj.setdefault(int(y), []).append(int(x))
        _bond_cache.clear()
        _bond_cache[key] = (m, adj)
    return _bond_cache[key][1].get(int(k), [])


class Env:
    """All atoms (heavy atoms first, then the hydrogens) with the attributes the contact scoring needs. Rotor and flip states
    change entries of these arrays, never the heavy-atom coordinates."""

    def __init__(self, m, H):
        self.m, self.H, self.n = m, H, m.n_atoms
        N = self.n + H.n
        self.P = np.vstack([m.xyz, H.xyz])
        self.rad = np.concatenate([heavy_radii(m), H.rad])
        self.is_h = np.concatenate([np.zeros(self.n, bool), np.ones(H.n, bool)])
        self.polar_h = np.concatenate([np.zeros(self.n, bool), H.polar])
        self.acc = np.concatenate([m.acceptor.copy(), np.zeros(H.n, bool)])
        self.chg = np.concatenate([m.charge.astype(int).copy(), H.charged.astype(int)])
        self.parent = np.concatenate([np.arange(self.n), H.parent])
        self.active = np.concatenate([np.ones(self.n, bool), H.active])
        self.hops = _hop_table(m, 4)

    def sync(self, H):
        self.P[self.n:] = H.xyz
        self.rad[self.n:] = H.rad
        self.polar_h[self.n:] = H.polar
        self.chg[self.n:] = H.charged.astype(int)
        self.parent[self.n:] = H.parent
        self.active[self.n:] = H.active


def _pair_energy(P1, rad1, polar1, acc1, chg1, isH1, parent1, ppos1, P2, rad2, polar2, acc2, chg2, isH2, ppos2, excl):
    """Contact energy between item sets 1 (own) and 2 (environment). excl[k, l] True: the pair is not scored."""
    d = np.linalg.norm(P1[:, None, :] - P2[None, :, :], axis=2)
    ov = rad1[:, None] + rad2[None, :] - d
    hb1 = polar1[:, None] & acc2[None, :] & isH1[:, None] & ~isH2[None, :]          # own polar H, env acceptor
    hb2 = (polar2[None, :] & isH2[None, :]) & acc1[:, None] & ~isH1[:, None]        # own acceptor, env polar H
    hb = hb1 | hb2
    # angle at the hydrogen: X-H...A must be 110 degrees or more
    u1 = ppos1 - P1
    u1 = u1 / np.maximum(np.linalg.norm(u1, axis=1, keepdims=True), 1e-9)
    u2 = ppos2 - P2
    u2 = u2 / np.maximum(np.linalg.norm(u2, axis=1, keepdims=True), 1e-9)
    v12 = (P2[None, :, :] - P1[:, None, :])
    v12 = v12 / np.maximum(np.linalg.norm(v12, axis=2, keepdims=True), 1e-9)
    cos1 = np.einsum("kc,klc->kl", u1, v12)                  # at own H toward env atom
    cos2 = np.einsum("lc,klc->kl", u2, -v12)                 # at env H toward own atom
    ang_ok = np.where(hb1, cos1 <= -0.34, np.where(hb2, cos2 <= -0.34, False))
    charged = (chg1[:, None] != 0) | (chg2[None, :] != 0)
    cap = np.where(charged, HB_CAP_CHARGED, HB_CAP)
    use_hb = hb & ang_ok
    pen_plain = np.where(ov > 0, 5.0 * ov ** 2 + 50.0 * np.maximum(ov - OVERLAP_CLASH, 0) ** 2, 0.0)
    pen_hb = np.where(ov > cap, 5.0 * ov ** 2 + 50.0 * np.maximum(ov - cap, 0) ** 2, 0.0)
    reward = 2.0 * np.clip((2.7 - d) / 0.7, 0, 1)
    e = np.where(use_hb, pen_hb - reward, pen_plain)
    e = np.where(excl, 0.0, e)
    return float(e.sum())


def _exclusion(env, own_parent, own_is_h, env_idx):
    """Pairs that are bonded or too close in the bond graph (four bonds or fewer, hydrogens included)."""
    k = len(own_parent)
    excl = np.zeros((k, len(env_idx)), bool)
    for r in range(k):
        hp = env.hops[own_parent[r]]
        for c, j in enumerate(env_idx):
            if j < env.n:
                if own_is_h[r]:
                    excl[r, c] = hp.get(j, 99) <= 3
                else:
                    excl[r, c] = True                   # heavy-heavy pairs are not scored here
            else:
                z = env.parent[j]
                excl[r, c] = hp.get(z, 99) <= (2 if own_is_h[r] else 3)
    return excl


def _rot_positions(m, g, phi):
    x, a, b = m.xyz[g.x], m.xyz[g.a], m.xyz[g.b]
    ph = [phi + k * 2 * np.pi / 3 for k in range(len(g.hidx))]
    return np.array([_rotor(x, a, b, g.length, p) for p in ph])


def _flip_state(m, H, env, g, s):
    n_atom, o_atom = (g.atoms[0], g.atoms[1]) if s == 0 else (g.atoms[1], g.atoms[0])
    p1, p2 = _n2(m.xyz[n_atom], m.xyz[g.a], m.xyz[g.b], g.length)
    return n_atom, o_atom, np.array([p1, p2])


def _his_state(m, g, flipped, prot):
    ring = g.ring
    xyz = {k: m.xyz[v] for k, v in ring.items()}
    pos, at = _his_roles(xyz, flipped)
    H_pos, H_par, H_act, H_pol = [], [], [], []
    for k in range(4):
        nb1, nb2 = _HIS_NB[k]
        n1 = pos[nb1] if not isinstance(nb1, int) else pos[nb1]
        n2 = pos[nb2]
        is_n = k in (0, 2)
        L = XH["N"] if is_n else XH["C"]
        H_pos.append(_p1(pos[k], n1, n2, L))
        H_par.append(ring[at[k]])
        act = True if not is_n else (prot in (0, 2) if k == 0 else prot in (1, 2))
        H_act.append(act)
    # heavy atom attributes per ring atom index
    attrs = {ring["CG"]: (1.75, False, 0)}
    for k in range(4):
        idx = ring[at[k]]
        if k in (0, 2):
            has_h = H_act[k]
            attrs[idx] = (1.55, not has_h, 1 if prot == 2 else 0)
        else:
            attrs[idx] = (1.75, False, 0)
    return np.array(H_pos), H_par, H_act, attrs


def _apply(m, H, env, g, state):
    """Write the chosen state of a group into H and env."""
    if g.kind == "rot":
        g.phi = state
        H.xyz[g.hidx] = _rot_positions(m, g, state)
    elif g.kind == "flip":
        g.state = state
        n_atom, o_atom, pos = _flip_state(m, H, env, g, state)
        H.xyz[g.hidx] = pos
        H.parent[g.hidx] = n_atom
        env.rad[n_atom], env.acc[n_atom] = 1.55, False
        env.rad[o_atom], env.acc[o_atom] = 1.40, True
    else:
        flipped, prot = state
        g.state = state
        pos, par, act, attrs = _his_state(m, g, flipped, prot)
        H.xyz[g.hidx] = pos
        H.parent[g.hidx] = par
        H.active[g.hidx] = act
        H.charged[g.hidx] = prot == 2
        for idx, (rad, acc, chg) in attrs.items():
            env.rad[idx], env.acc[idx], env.chg[idx] = rad, acc, chg
    env.sync(H)


def _candidates(m, H, env, g):
    """Yield (state, own arrays) for every state of the group."""
    if g.kind == "rot":
        for phi in g.cands:
            pos = _rot_positions(m, g, phi)
            yield phi, 0.0, pos, [g.x] * len(g.hidx), np.ones(len(g.hidx), bool), None
    elif g.kind == "flip":
        for s in (0, 1):
            n_atom, o_atom, pos = _flip_state(m, H, env, g, s)
            heavy = {n_atom: (1.55, False, 0), o_atom: (1.40, True, 0)}
            yield s, FLIP_COST * s, pos, [n_atom, n_atom], np.ones(2, bool), heavy
    else:
        for flipped in (0, 1):
            for prot in (0, 1, 2):
                pos, par, act, attrs = _his_state(m, g, flipped, prot)
                extra = FLIP_COST * flipped + (1.0 if prot == 2 else 0.0)
                yield (flipped, prot), extra, pos, par, np.array(act), attrs


def optimise(m, H, groups, env, passes=3):
    """Choose the state of every rotor and flip group, a few passes, each group against the current state of all others."""
    n = env.n
    for g in groups:
        _apply(m, H, env, g, g.phi if g.kind == "rot" else g.state)
    for _ in range(passes):
        changed = 0
        tree_h = cKDTree(env.P[:n]) if _ == 0 else tree_h
        tree_hh = cKDTree(env.P[n:]) if H.n else None
        for g in groups:
            hid = np.array(g.hidx)
            center = env.P[n + hid].mean(axis=0)
            nb_h = tree_h.query_ball_point(center, 6.0)
            nb_hh = [n + k for k in tree_hh.query_ball_point(center, 6.0)] if tree_hh is not None else []
            own_heavy = set()
            if g.kind == "flip":
                own_heavy = set(g.atoms)
            elif g.kind == "his":
                own_heavy = set(g.ring.values())
            own_hset = set(n + hid)
            env_idx = np.array([j for j in nb_h if j not in own_heavy] + [j for j in nb_hh if j not in own_hset and env.active[j]], dtype=int)
            if len(env_idx) == 0:
                continue
            best, best_state = None, None
            e_cache = {}
            for state, extra, pos, par, act, heavy in _candidates(m, H, env, g):
                items_P, items_rad, items_polar, items_acc, items_chg, items_isH, items_par, items_ppos = [], [], [], [], [], [], [], []
                for k, hi in enumerate(hid):
                    if not act[k]:
                        continue
                    items_P.append(pos[k])
                    items_rad.append(H.rad[hi])
                    items_polar.append(H.polar[hi])
                    items_acc.append(False)
                    chg_h = (state[1] == 2) if g.kind == "his" else H.charged[hi]
                    items_chg.append(int(chg_h))
                    items_isH.append(True)
                    items_par.append(par[k] if g.kind != "rot" else g.x)
                    items_ppos.append(m.xyz[par[k] if g.kind != "rot" else g.x])
                if heavy:
                    for idx, (rad, acc, chg) in heavy.items():
                        items_P.append(m.xyz[idx])
                        items_rad.append(rad)
                        items_polar.append(False)
                        items_acc.append(acc)
                        items_chg.append(chg)
                        items_isH.append(False)
                        items_par.append(idx)
                        items_ppos.append(m.xyz[idx])
                P1 = np.array(items_P)
                excl = _exclusion(env, np.array(items_par), np.array(items_isH), env_idx)
                e = _pair_energy(P1, np.array(items_rad), np.array(items_polar), np.array(items_acc), np.array(items_chg),
                                 np.array(items_isH), np.array(items_par), np.array(items_ppos),
                                 env.P[env_idx], env.rad[env_idx], env.polar_h[env_idx], env.acc[env_idx], env.chg[env_idx],
                                 env.is_h[env_idx], env.P[env.parent[env_idx]], excl) + extra
                if best is None or e < best - 1e-9:
                    best, best_state = e, state
            cur = g.phi if g.kind == "rot" else g.state
            if best_state != cur:
                changed += 1
            _apply(m, H, env, g, best_state)
        if not changed:
            break
    return H


def add_hydrogens(m, passes=3):
    """Hydrogens of a loaded Model, optimised. Returns (H, env, groups)."""
    H, groups = build(m)
    env = Env(m, H)
    if groups:
        optimise(m, H, groups, env, passes)
    return H, env, groups


def clashes(m, H, env, min_overlap=OVERLAP_CLASH, hb_cap=HB_CAP, hb_cap_charged=HB_CAP_CHARGED, hh_hops=HH_HOPS, heavy_hops=3, raw=False, heavy=False):
    """All-atom clashes that involve at least one hydrogen. Returns a list of dicts (indices in the combined atom list).

    raw=True keeps every pair with an overlap of 0.2 A or more and adds the fields `hb` (polar hydrogen next to an acceptor at a
    hydrogen-bond angle) and `charged`, so that the thresholds can be calibrated. heavy=True also returns the clashes between
    two heavy atoms, with the radii used here."""
    n = env.n
    keep = np.where(env.active)[0]
    P = env.P[keep]
    tree = cKDTree(P)
    pairs = tree.query_pairs(r=3.4, output_type="ndarray")
    if len(pairs) == 0:
        return []
    a, b = keep[pairs[:, 0]], keep[pairs[:, 1]]
    sel = (a >= n) | (b >= n) if not heavy else np.ones(len(a), bool)
    a, b = a[sel], b[sel]
    swap = (a < n) & (b >= n)                                   # the hydrogen goes first
    a, b = np.where(swap, b, a), np.where(swap, a, b)
    d = np.linalg.norm(env.P[a] - env.P[b], axis=1)
    ov = env.rad[a] + env.rad[b] - d
    ok = ov >= (0.2 if raw else min_overlap - 1e-9)
    a, b, d, ov = a[ok], b[ok], d[ok], ov[ok]
    out = []
    for x, y, dist, o in zip(a, b, d, ov):
        px = env.parent[x]
        hop = env.hops[px].get(y if y < n else env.parent[y], 99)
        if x < n:                                               # two heavy atoms (only with heavy=True)
            if hop <= 4 or m.is_metal[x] or m.is_metal[y] or (m.is_water[x] and m.is_water[y]):
                continue
            if not raw and m.element[x] in ("N", "O") and m.element[y] in ("N", "O") and o < HEAVY_POLAR_CAP:
                continue
            out.append({"i": int(x), "j": int(y), "d": float(dist), "overlap": float(o), "heavy": True, "hb": False, "charged": False, "hh": False, "polar": False, "hop": 99})
            continue
        if y < n:
            if hop <= (1 if raw else heavy_hops):
                continue
        elif hop <= (0 if raw else hh_hops):
            continue
        hb = False
        if env.polar_h[x] and y < n and env.acc[y]:                        # hydrogen bond
            u = m.xyz[px] - env.P[x]
            v = env.P[y] - env.P[x]
            cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-9)
            hb = bool(cos <= -0.34)
        charged = bool(env.chg[x] != 0 or (y < n and env.chg[y] != 0))
        if not raw:
            if hb and o < (hb_cap_charged if charged else hb_cap):
                continue
        rec = {"i": int(x), "j": int(y), "d": float(dist), "overlap": float(o)}
        if raw:
            rec.update(hb=hb, charged=charged, hh=bool(y >= n), polar=bool(env.polar_h[x]), hop=int(hop))
        out.append(rec)
    return out
