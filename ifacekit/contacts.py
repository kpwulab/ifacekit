"""Atom-pair table: one row per pair of non-bonded atoms that are close.

Contacts and clashes live in the same table. Every row gets a *type*
(what the pair looks like) and a *verdict* (is it chemically reasonable):

    verdict   meaning
    -------   -------------------------------------------------------------
    favorable good hydrogen bond or salt bridge geometry
    ok        ordinary van der Waals / hydrophobic / metal contact
    borderline  close contact, long-range odd polar pair, like charges
    unfavorable  steric overlap, impossibly short polar contact

The engine uses heavy atoms only. Without hydrogens hydrogen-bond angles
cannot be tested, so hydrogen bonds are judged on donor-acceptor distance.
Thresholds are provisional and meant to be calibrated on high resolution
structures.
"""
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from .model import within_n_bonds

# Distance thresholds (angstrom)
CANDIDATE_CUTOFF = 4.5
# Pairs up to this many covalent bonds apart are never compared. Four (1-5)
# is needed because O(i)...C(i+1) across a peptide bond is about 2.8-3.0 A in
# helices, which would otherwise look like a steric overlap.
EXCLUDE_BOND_DEPTH = 4
HB_MAX = 3.5
# Hydrogen bond with geometry (no hydrogens are used): donor X with heavy neighbours, acceptor A with heavy neighbours.
# Calibrated against PISA, PyMOL and ChimeraX on 280 interfaces (tools and manuscript/hbond_variants_heavy.py).
HBGEO_MAX = 3.3                  # X...A distance
HBGEO_DONOR_DIR = 50.0           # donor with two or more neighbours: angle between the expected N-H direction and X->A, at most
HBGEO_DONOR_NEIGHBOUR = 100.0    # donor with one neighbour: angle N-X...A, at least
HBGEO_ACCEPTOR = 100.0           # every angle Y-A...X at the acceptor, at least
HB_OK_MAX = 3.35
HB_SHORT_BAD = 2.30
HB_SHORT_BAD_WATER = 2.10
HB_SHORT_WARN = 2.55
SALT_MAX = 4.0
HYDROPHOBIC_MAX = 4.0
METAL_MAX = 3.2
MISMATCH_MAX = 3.1
MISMATCH_BAD = 2.5
LIKE_CHARGE_MAX = 4.0
LIKE_CHARGE_BAD = 3.0
OVERLAP_BAD = 0.60
OVERLAP_WARN = 0.25
VDW_MARGIN = 0.5

TYPES = ["metal", "salt_bridge", "hbond", "polar_mismatch", "like_charge",
         "hydrophobic", "vdw"]


@dataclass
class ContactSet:
    i: np.ndarray
    j: np.ndarray
    d: np.ndarray
    overlap: np.ndarray
    type: np.ndarray
    verdict: np.ndarray
    reason: np.ndarray
    hbond: np.ndarray
    inter_unit: np.ndarray
    satisfied: np.ndarray      # per atom: has an h-bond / salt bridge / metal partner
    hbond_geometry: np.ndarray = None   # hydrogen bond with geometry (see _hbond_geometry); the hydrogen bonds that are counted
    counts: dict = None        # tallies over ALL typed pairs, recorded or not

    def __len__(self):
        return len(self.i)


def _angle(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return 0.0 if nu < 1e-9 or nv < 1e-9 else float(np.degrees(np.arccos(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))))


def _neighbours(m):
    """Bonded heavy neighbours of every atom. Atoms on the same position (duplicates, for example two chains placed on top of each other) are no neighbours."""
    nb = [[] for _ in range(m.n_atoms)]
    for a, b in m.bonds:
        if np.linalg.norm(m.xyz[a] - m.xyz[b]) < 0.8:
            continue
        nb[a].append(int(b))
        nb[b].append(int(a))
    return nb


def _geometry_ok(m, nb, don, acc):
    """Geometry of a donor atom and an acceptor atom without hydrogens: the donor points to the acceptor and the acceptor accepts from the donor's side."""
    x, a = m.xyz[don], m.xyz[acc]
    xa = a - x
    dn = nb[don]
    if len(dn) >= 2:
        h = -sum((m.xyz[k] - x) / max(np.linalg.norm(m.xyz[k] - x), 1e-9) for k in dn)       # expected direction of the N-H bond
        if np.linalg.norm(h) > 1e-6 and _angle(h, xa) > HBGEO_DONOR_DIR:
            return False
    elif len(dn) == 1 and _angle(m.xyz[dn[0]] - x, xa) < HBGEO_DONOR_NEIGHBOUR:
        return False
    return all(_angle(m.xyz[y] - a, x - a) >= HBGEO_ACCEPTOR for y in nb[acc])


def hbond_geometry(m, i, j, d, hb_cand):
    """Boolean array over the candidate pairs: a hydrogen bond with geometry. A pair needs complementary roles, X...A <= HBGEO_MAX and, when both roles are certain,
    the geometry of _geometry_ok in either direction (a water has no neighbours and passes; atoms with an uncertain role are judged by the distance only)."""
    out = np.zeros(len(i), bool)
    idx = np.where(hb_cand & (d <= HBGEO_MAX))[0]
    if len(idx) == 0:
        return out
    nb = _neighbours(m)
    for k in idx:
        a, b = int(i[k]), int(j[k])
        if m.ambiguous[a] or m.ambiguous[b]:
            out[k] = True
            continue
        ok = (m.donor[a] and m.acceptor[b] and _geometry_ok(m, nb, a, b)) or (m.donor[b] and m.acceptor[a] and _geometry_ok(m, nb, b, a))
        out[k] = bool(ok)
    return out


def find_contacts(m, cutoff=CANDIDATE_CUTOFF, record_intra="problems", restrict=None):
    """Build the atom-pair table for a loaded Model.

    restrict: optional boolean array over the atoms; only pairs with at least one of these atoms are kept and tallied
    (used to add the contacts of alternate conformers).
    record_intra: "problems" keeps intra-unit pairs only when they are
    borderline/unfavorable; "all" keeps every typed intra-unit pair.
    """
    n = m.n_atoms
    empty = ContactSet(*(np.zeros(0, dtype=t) for t in
                         (int, int, float, float, object, object, object, bool, bool)),
                       satisfied=np.zeros(n, bool), hbond_geometry=np.zeros(0, bool),
                       counts={"verdicts": {}, "verdicts_between_units": {},
                               "verdicts_within_units": {}, "types": {}})
    if n < 2:
        return empty
    tree = cKDTree(m.xyz)
    pairs = tree.query_pairs(r=cutoff, output_type="ndarray")
    if len(pairs) == 0:
        return empty
    i, j = pairs[:, 0].astype(np.int64), pairs[:, 1].astype(np.int64)

    # drop water-water and 1-2/1-3/1-4 bonded pairs
    keep = ~(m.is_water[i] & m.is_water[j])
    excl = within_n_bonds(m.bonds, n, depth=EXCLUDE_BOND_DEPTH)
    if len(excl):
        keep &= ~np.isin(i * n + j, excl)
    i, j = i[keep], j[keep]
    if len(i) == 0:
        return empty

    d = np.linalg.norm(m.xyz[i] - m.xyz[j], axis=1)
    overlap = m.radius[i] + m.radius[j] - d
    ui, uj = m.unit_idx[i], m.unit_idx[j]
    inter = ~((ui == uj) & (ui >= 0))

    polar_i = np.isin(m.element[i], ["N", "O"])
    polar_j = np.isin(m.element[j], ["N", "O"])
    both_polar = polar_i & polar_j
    metal = m.is_metal[i] | m.is_metal[j]
    both_metal = m.is_metal[i] & m.is_metal[j]
    chg = m.charge[i] * m.charge[j]

    hb_cand = both_polar & ((m.donor[i] & m.acceptor[j]) | (m.acceptor[i] & m.donor[j]))
    definite = ~m.ambiguous[i] & ~m.ambiguous[j]
    don_only_i = m.donor[i] & ~m.acceptor[i]
    don_only_j = m.donor[j] & ~m.acceptor[j]
    acc_only_i = m.acceptor[i] & ~m.donor[i]
    acc_only_j = m.acceptor[j] & ~m.donor[j]
    mismatch = both_polar & definite & ((don_only_i & don_only_j) | (acc_only_i & acc_only_j))

    metal_partner = np.isin(m.element[i], ["N", "O", "S", "CL"]) | np.isin(m.element[j], ["N", "O", "S", "CL"])
    t_metal = metal & ~both_metal & metal_partner & (d <= METAL_MAX)
    nm = ~metal
    t_salt = nm & both_polar & (chg == -1) & (d <= SALT_MAX)
    t_hb = nm & hb_cand & (d <= HB_MAX)
    t_mis = nm & mismatch & (d < MISMATCH_MAX)
    t_like = nm & both_polar & (chg == 1) & (d < LIKE_CHARGE_MAX)
    t_hyd = nm & m.apolar[i] & m.apolar[j] & (d <= HYDROPHOBIC_MAX)
    t_vdw = nm & (d <= m.radius[i] + m.radius[j] + VDW_MARGIN)

    type_idx = np.full(len(i), -1)
    for k, mask in reversed(list(enumerate([t_metal, t_salt, t_hb, t_mis, t_like, t_hyd, t_vdw]))):
        type_idx[mask] = k
    typed = type_idx >= 0
    if restrict is not None:
        typed = typed & (restrict[i] | restrict[j])

    # satisfied partners (for buried unsatisfied polar atoms), all typed pairs
    satisfied = np.zeros(n, bool)
    sat_pair = (t_metal | t_salt | (t_hb & (d >= HB_SHORT_WARN)))
    satisfied[i[sat_pair]] = True
    satisfied[j[sat_pair]] = True

    verdict = np.full(len(i), "ok", dtype=object)
    reason = np.full(len(i), "", dtype=object)

    def setv(mask, v, r=""):
        verdict[mask] = v
        reason[mask] = r

    # hydrogen bond
    hbm = type_idx == 2
    setv(hbm & (d >= HB_SHORT_WARN) & (d <= HB_OK_MAX), "favorable")
    setv(hbm & (d > HB_OK_MAX), "ok", "long_hbond")
    setv(hbm & (d < HB_SHORT_WARN), "borderline", "short_hbond")
    wat = m.is_water[i] | m.is_water[j]
    setv(hbm & (d < np.where(wat, HB_SHORT_BAD_WATER, HB_SHORT_BAD)), "unfavorable", "short_hbond")
    # salt bridge
    sm = type_idx == 1
    setv(sm & (d >= HB_SHORT_WARN), "favorable")
    setv(sm & (d < HB_SHORT_WARN), "borderline", "short_polar")
    setv(sm & (d < HB_SHORT_BAD), "unfavorable", "short_polar")
    # polar mismatch (donor-donor or acceptor-acceptor)
    mm = type_idx == 3
    setv(mm, "borderline", "polar_mismatch")
    setv(mm & (d < MISMATCH_BAD), "unfavorable", "polar_mismatch")
    # a carboxylate or phosphate oxygen near another acceptor may simply be protonated
    # (hydrogens are not modelled here), so such pairs are never called unfavorable
    acid = (m.charge[i] == -1) | (m.charge[j] == -1)
    setv(mm & acid, "borderline", "possible_protonation")
    # like charges
    lm = type_idx == 4
    setv(lm, "borderline", "like_charge")
    setv(lm & (d < LIKE_CHARGE_BAD), "unfavorable", "like_charge")
    # hydrophobic and generic vdW: steric overlap rules
    gm = (type_idx == 5) | (type_idx == 6)
    setv(gm & (overlap >= OVERLAP_WARN), "borderline", "close_contact")
    setv(gm & (overlap >= OVERLAP_BAD), "unfavorable", "steric_overlap")
    # hydrogen-bond capable pairs that ended up in the generic class because
    # they are beyond HB_MAX never overlap, so nothing more to do.

    type_names = np.array(TYPES + [""], dtype=object)[type_idx]
    hb_flag = t_hb | (t_salt & hb_cand & (d <= HB_MAX))
    hb_geo = hbond_geometry(m, i, j, d, hb_cand & nm & both_polar) if (typed.any()) else np.zeros(len(i), bool)

    problem = np.isin(verdict, ["borderline", "unfavorable"])
    if record_intra == "all":
        rec = typed
    else:
        rec = typed & (inter | problem)

    def tally(mask):
        vals, cnt = np.unique(verdict[mask].astype(str), return_counts=True)
        return {str(v): int(c) for v, c in zip(vals, cnt)}

    tvals, tcnt = np.unique(type_names[typed].astype(str), return_counts=True)
    counts = {"verdicts": tally(typed),
              "verdicts_between_units": tally(typed & inter),
              "verdicts_within_units": tally(typed & ~inter),
              "types": {str(v): int(c) for v, c in zip(tvals, tcnt)}}

    return ContactSet(i=i[rec], j=j[rec], d=d[rec], overlap=overlap[rec],
                      type=type_names[rec], verdict=verdict[rec],
                      reason=reason[rec], hbond=hb_flag[rec], hbond_geometry=hb_geo[rec],
                      inter_unit=inter[rec], satisfied=satisfied, counts=counts)
