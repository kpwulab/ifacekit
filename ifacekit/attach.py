"""Decide which polymer chain each ligand, glycan and water belongs to.

Chain labels are not reliable for this: refinement programs often put all
waters in their own chain ("S", "W") and ligands in a chain of their own.
Ownership is therefore taken from the structure itself:

  1. A ligand or glycan that is covalently bonded to a polymer residue belongs
     to that polymer's chain. Residues linked to each other by covalent bonds
     (a glycan tree, a multi-residue ligand) are one group and move together.
  2. Otherwise a group belongs to the polymer chain(s) it touches, that is,
     the chains with a heavy atom within CONTACT_CUTOFF of any atom of the
     group, ordered by number of such atom pairs.
  3. A water belongs to the polymer chain(s) with an atom within
     CONTACT_CUTOFF, nearest first.

One group may touch several chains. The first chain in its list is its owner
when a single owner is needed (merged ligand units in --pair mode).
"""
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import cKDTree

from .model import load

CONTACT_CUTOFF = 4.0


def plan_attachments(path, model_index=0):
    """Return {residue_id: [chain, ...]} for every ligand, glycan and water."""
    m = load(path, model_index)
    res_kind = np.array([r["kind"] for r in m.residues])
    res_chain = np.array([r["chain"] for r in m.residues])
    atom_kind = res_kind[m.res_idx]
    atom_chain = res_chain[m.res_idx]
    poly = np.where(atom_kind == "polymer")[0]
    owners = {}
    if len(poly) == 0:
        return owners
    tree = cKDTree(m.xyz[poly])

    # ---- groups of covalently linked non-polymer residues
    parent = list(range(len(m.residues)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    nonpoly = {"ligand", "branched"}
    covalent = defaultdict(Counter)       # residue index -> Counter(chain)
    for a, b in m.bonds:
        ra, rb = int(m.res_idx[a]), int(m.res_idx[b])
        if ra == rb:
            continue
        ka, kb = res_kind[ra], res_kind[rb]
        if ka in nonpoly and kb in nonpoly:
            parent[find(ra)] = find(rb)
        elif ka in nonpoly and kb == "polymer":
            covalent[ra][res_chain[rb]] += 1
        elif kb in nonpoly and ka == "polymer":
            covalent[rb][res_chain[ra]] += 1

    groups = defaultdict(list)
    for ri, k in enumerate(res_kind):
        if k in nonpoly:
            groups[find(ri)].append(ri)
    atoms_of = defaultdict(list)
    for a in range(m.n_atoms):
        atoms_of[int(m.res_idx[a])].append(a)

    for members in groups.values():
        cov = Counter()
        for ri in members:
            cov.update(covalent.get(ri, {}))
        contact = Counter()
        for ri in members:
            for a in atoms_of[ri]:
                for p in tree.query_ball_point(m.xyz[a], CONTACT_CUTOFF):
                    contact[atom_chain[poly[p]]] += 1
        ordered = [c for c, _ in sorted(cov.items(), key=lambda kv: (-kv[1], kv[0]))]
        ordered += [c for c, _ in sorted(contact.items(), key=lambda kv: (-kv[1], kv[0]))
                    if c not in ordered]
        for ri in members:
            owners[m.residues[ri]["id"]] = ordered

    # ---- waters: nearest chain first
    for ri, k in enumerate(res_kind):
        if k != "water":
            continue
        best = {}
        for a in atoms_of[ri]:
            for p in tree.query_ball_point(m.xyz[a], CONTACT_CUTOFF):
                d = float(np.linalg.norm(m.xyz[a] - m.xyz[poly[p]]))
                c = atom_chain[poly[p]]
                if c not in best or d < best[c]:
                    best[c] = d
        owners[m.residues[ri]["id"]] = [c for c, _ in sorted(best.items(), key=lambda kv: kv[1])]
    return owners
