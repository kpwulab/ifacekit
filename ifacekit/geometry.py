"""Backbone bond-length and bond-angle check against approximate Engh and Huber values.

Only protein backbone and C-beta are covered, using fixed ideal values (no restraint
library is read). Values are approximate and are meant to find gross distortions,
not to replace a refinement program's restraint report.
"""
import numpy as np

# (atom1, atom2): (ideal, sigma); '+' marks the next residue
BONDS = {
    ("N", "CA"): (1.459, 0.020), ("CA", "C"): (1.525, 0.026), ("C", "O"): (1.229, 0.019),
    ("CA", "CB"): (1.530, 0.020), ("C", "+N"): (1.336, 0.023),
}
# (a, center, c): (ideal deg, sigma)
ANGLES = {
    ("N", "CA", "C"): (111.2, 2.8), ("N", "CA", "CB"): (110.5, 1.7), ("CB", "CA", "C"): (110.1, 1.9),
    ("CA", "C", "O"): (120.1, 2.1), ("CA", "C", "+N"): (117.2, 2.2), ("O", "C", "+N"): (122.7, 1.6),
    ("C", "+N", "+CA"): (121.7, 1.8),
}
GLY_N_CA_C = (113.1, 2.5)
PRO_N_CA_C = (112.1, 2.6)
Z_OUTLIER = 4.0
BREAK_DIST = 2.0


def backbone_geometry(m, residue_ids=None):
    prot = [k for k, r in enumerate(m.residues) if r["protein"]]
    atoms = {}
    for a in range(m.n_atoms):
        atoms.setdefault(int(m.res_idx[a]), {})[m.atom_name[a]] = a
    pos = lambda a: m.xyz[a]
    chain_res = {}
    for k in prot:
        chain_res.setdefault(m.residues[k]["chain"], []).append(k)
    bond_z, angle_z, outliers, breaks = [], [], [], []
    for ch, ks in chain_res.items():
        for n, k in enumerate(ks):
            r = m.residues[k]
            if residue_ids is not None and r["id"] not in residue_ids:
                continue
            nxt = None
            if n + 1 < len(ks):
                cand = ks[n + 1]
                A, B = atoms.get(k, {}), atoms.get(cand, {})
                if "C" in A and "N" in B:
                    dist = float(np.linalg.norm(pos(A["C"]) - pos(B["N"])))
                    if dist <= BREAK_DIST:
                        nxt = cand
                    else:
                        gap = m.residues[cand]["resnum"] - r["resnum"]
                        breaks.append({"after": r["id"], "before": m.residues[cand]["id"], "C_N_distance": round(dist, 2),
                                       "kind": "stretched_bond" if gap == 1 else "gap_in_numbering"})
            here = atoms.get(k, {})
            there = atoms.get(nxt, {}) if nxt is not None else {}

            def get(name):
                if name.startswith("+"):
                    return there.get(name[1:])
                return here.get(name)
            for (a, b), (ideal, sig) in BONDS.items():
                ia, ib = get(a), get(b)
                if ia is None or ib is None:
                    continue
                v = float(np.linalg.norm(pos(ia) - pos(ib)))
                z = (v - ideal) / sig
                bond_z.append(z)
                if abs(z) > Z_OUTLIER:
                    outliers.append({"residue": r["id"], "kind": "bond", "atoms": f"{a}-{b}", "value": round(v, 3),
                                     "ideal": ideal, "z": round(z, 1)})
            for (a, c, e), (ideal, sig) in ANGLES.items():
                ia, ic, ie = get(a), get(c), get(e)
                if ia is None or ic is None or ie is None:
                    continue
                if (a, c, e) == ("N", "CA", "C"):
                    ideal, sig = GLY_N_CA_C if r["resname"] == "GLY" else PRO_N_CA_C if r["resname"] == "PRO" else (ideal, sig)
                v1, v2 = pos(ia) - pos(ic), pos(ie) - pos(ic)
                ang = float(np.degrees(np.arccos(np.clip(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1))))
                z = (ang - ideal) / sig
                angle_z.append(z)
                if abs(z) > Z_OUTLIER:
                    outliers.append({"residue": r["id"], "kind": "angle", "atoms": f"{a}-{c}-{e}", "value": round(ang, 1),
                                     "ideal": ideal, "z": round(z, 1)})
    outliers.sort(key=lambda x: -abs(x["z"]))
    rms = lambda v: round(float(np.sqrt(np.mean(np.square(v)))), 2) if v else None
    return {"n_bonds": len(bond_z), "n_angles": len(angle_z), "rmsz_bonds": rms(bond_z), "rmsz_angles": rms(angle_z),
            "n_outliers": len(outliers), "outliers": outliers[:200], "chain_breaks": breaks,
            "note": "approximate Engh-Huber backbone values, |Z| > 4 listed"}
