"""Aromatic ring interactions (stacking, T-shaped, cation-pi) for standard residues."""
import numpy as np
from scipy.spatial import cKDTree

from . import chem

STACK_MAX = 5.5
TSHAPE_MAX = 5.5
CATION_PI_MAX = 6.0


def _plane(points):
    c = points.mean(axis=0)
    u, s, vt = np.linalg.svd(points - c)
    n = vt[2]
    return c, n / np.linalg.norm(n)


def _angle(a, b):
    cosv = abs(float(np.dot(a, b)))
    return float(np.degrees(np.arccos(min(1.0, cosv))))


def find_pi(m):
    """Return a list of dict records of ring based interactions."""
    # atom lookup per residue
    by_res = {}
    for k in range(m.n_atoms):
        by_res.setdefault(int(m.res_idx[k]), {})[m.atom_name[k]] = k

    rings = []       # (res_idx, centroid, normal, ring_label)
    cations = []     # (res_idx, atom_idx)
    for ridx, atoms in by_res.items():
        r = m.residues[ridx]
        if r["kind"] != "polymer" or not r["standard"]:
            continue
        for rl, names in enumerate(chem.RINGS.get(r["resname"], [])):
            if all(nm in atoms for nm in names):
                pts = m.xyz[[atoms[nm] for nm in names]]
                c, nrm = _plane(pts)
                rings.append((ridx, c, nrm, rl))
        for (rn, an) in chem.CATIONS:
            if r["resname"] == rn and an in atoms:
                cations.append((ridx, atoms[an]))

    out = []
    if len(rings) >= 2:
        cents = np.array([x[1] for x in rings])
        tree = cKDTree(cents)
        for a, b in tree.query_pairs(r=max(STACK_MAX, TSHAPE_MAX), output_type="ndarray"):
            ra, ca, na, _ = rings[a]
            rb, cb, nb, _ = rings[b]
            if ra == rb:
                continue
            v = cb - ca
            dist = float(np.linalg.norm(v))
            ang = _angle(na, nb)
            perp = abs(float(np.dot(v, na)))
            offset = float(np.sqrt(max(dist * dist - perp * perp, 0.0)))
            if ang <= 30.0 and dist <= STACK_MAX and offset <= 3.5 and perp >= 3.0:
                kind = "pi_stacking"
            elif ang >= 60.0 and dist <= TSHAPE_MAX:
                kind = "pi_tshaped"
            else:
                continue
            out.append({"type": kind, "res1": m.residues[ra]["id"],
                        "res2": m.residues[rb]["id"], "distance": round(dist, 2),
                        "angle": round(ang, 1), "offset": round(offset, 2),
                        "inter_unit": m.residues[ra]["unit"] != m.residues[rb]["unit"]})
    if rings and cations:
        cent_tree = cKDTree(np.array([x[1] for x in rings]))
        for cres, ai in cations:
            for rk in cent_tree.query_ball_point(m.xyz[ai], CATION_PI_MAX):
                rres, cc, nrm, _ = rings[rk]
                if rres == cres or m.residues[rres]["resname"] == "HIS":
                    continue
                v = m.xyz[ai] - cc
                dist = float(np.linalg.norm(v))
                ang = _angle(nrm, v / dist)
                if ang <= 40.0:
                    out.append({"type": "cation_pi", "res1": m.residues[cres]["id"],
                                "res2": m.residues[rres]["id"],
                                "distance": round(dist, 2), "angle": round(ang, 1),
                                "inter_unit": m.residues[cres]["unit"] != m.residues[rres]["unit"]})
    return out
