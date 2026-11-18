"""AlphaFold 3 confidence files and the interface scores built from them.

Summary file (`*_summary_confidences_*.json`): ipTM, pTM, ranking score, chain pair ipTM and
minimum PAE. Full-data file (`*_full_data_*.json`): the token PAE matrix, needed for ipSAE and
interface PAE. Without the full-data file only the summary numbers are available.

ipSAE follows Dunbrack's definition (PAE cutoff 10, d0 from the number of partner residues
below the cutoff, maximum over the aligned residues). This is a re-implementation, so compare
against the author's script before quoting numbers in a paper. pDockQ follows Bryant et al.
(2022) with CB contacts within 8 A.
"""
import json
import math

import numpy as np

from .model import load


def _unique(seq):
    out = []
    for x in seq:
        if x not in out:
            out.append(x)
    return out


def load_summary(path):
    with open(path) as fh:
        d = json.load(fh)
    chains = _unique(d.get("chain_ids", []))
    return {
        "chains": chains,
        "iptm": d.get("iptm"), "ptm": d.get("ptm"), "ranking_score": d.get("ranking_score"),
        "has_clash": d.get("has_clash"), "fraction_disordered": d.get("fraction_disordered"),
        "chain_pair_iptm": d.get("chain_pair_iptm"), "chain_pair_pae_min": d.get("chain_pair_pae_min"),
        "chain_iptm": d.get("chain_iptm"), "chain_ptm": d.get("chain_ptm"),
    }


def load_full(path):
    with open(path) as fh:
        d = json.load(fh)
    return {"pae": np.array(d["pae"], dtype=float),
            "chain": np.array(d["token_chain_ids"]),
            "res": np.array(d["token_res_ids"], dtype=int)}


def _d0(n):
    n = max(float(n), 27.0)
    return max(1.0, 1.24 * (n - 15.0) ** (1.0 / 3.0) - 1.8)


def ipsae(full, a, b, cutoff=10.0):
    """ipSAE for the direction a -> b (rows are aligned on chain a)."""
    pae, ch = full["pae"], full["chain"]
    ia, ib = np.where(ch == a)[0], np.where(ch == b)[0]
    if len(ia) == 0 or len(ib) == 0:
        return None
    sub = pae[np.ix_(ia, ib)]
    best = 0.0
    for row in sub:
        ok = row < cutoff
        if not ok.any():
            continue
        d0 = _d0(ok.sum())
        best = max(best, float(np.mean(1.0 / (1.0 + (row[ok] / d0) ** 2))))
    return best


def interface_pae(full, a, b, res_a, res_b):
    """Mean PAE between the interface residues of a and b, both directions, plus the best value."""
    pae, ch, rs = full["pae"], full["chain"], full["res"]
    ia = [k for k in range(len(ch)) if ch[k] == a and rs[k] in res_a]
    ib = [k for k in range(len(ch)) if ch[k] == b and rs[k] in res_b]
    if not ia or not ib:
        return None
    ab = pae[np.ix_(ia, ib)]
    ba = pae[np.ix_(ib, ia)]
    return {"a_to_b": round(float(ab.mean()), 2), "b_to_a": round(float(ba.mean()), 2),
            "mean": round(float((ab.mean() + ba.mean()) / 2), 2),
            "fraction_below_5": round(float(((ab < 5).sum() + (ba < 5).sum()) / (ab.size + ba.size)), 3)}


def pdockq(path, a, b, model_index=0):
    """pDockQ from CB contacts (CA for Gly) within 8 A and the mean pLDDT (B-factor) of those residues."""
    m = load(path, model_index, select=lambda ch, kind, rid: kind == "polymer" and ch in (a, b))
    from scipy.spatial import cKDTree
    cb = {}
    for k in range(m.n_atoms):
        r = m.residues[m.res_idx[k]]
        if not r["protein"]:
            continue
        key = int(m.res_idx[k])
        if m.atom_name[k] == "CB" or (m.atom_name[k] == "CA" and r["resname"] == "GLY"):
            cb[key] = k
    ra = [k for k in cb if m.residues[k]["chain"] == a]
    rb = [k for k in cb if m.residues[k]["chain"] == b]
    if not ra or not rb:
        return None
    ta = cKDTree(m.xyz[[cb[k] for k in ra]])
    tb = cKDTree(m.xyz[[cb[k] for k in rb]])
    pairs = ta.query_ball_tree(tb, 8.0)
    n = sum(len(p) for p in pairs)
    if n == 0:
        return {"pdockq": 0.0, "n_contacts": 0, "mean_plddt": None}
    inter = {ra[i] for i, p in enumerate(pairs) if p} | {rb[j] for p in pairs for j in p}
    bs = []
    for rr in inter:
        bs.append(float(np.mean([m.bfactor[k] for k in range(m.n_atoms) if m.res_idx[k] == rr])))
    mean_b = float(np.mean(bs))
    x = mean_b * math.log10(n)
    val = 0.724 / (1.0 + math.exp(-0.052 * (x - 152.611))) + 0.018
    return {"pdockq": round(val, 3), "n_contacts": int(n), "mean_plddt": round(mean_b, 1)}


def residue_plddt(path, model_index=0):
    """Per-residue mean B-factor (pLDDT in an AF model) for every polymer residue, in file order."""
    m = load(path, model_index)
    acc = {}
    for k in range(m.n_atoms):
        acc.setdefault(int(m.res_idx[k]), []).append(m.bfactor[k])
    rows = []
    for rr, r in enumerate(m.residues):
        if r["kind"] == "polymer" and rr in acc:
            rows.append({"chain": r["chain"], "resnum": r["resnum"], "resname": r["resname"],
                         "b": round(float(np.mean(acc[rr])), 1)})
    return rows


def plddt_bands(values):
    v = np.array(values, dtype=float)
    if len(v) == 0:
        return None
    return {"mean": round(float(v.mean()), 1), "n": int(len(v)),
            "very_high_ge90": round(float((v >= 90).mean()), 3),
            "confident_70_90": round(float(((v >= 70) & (v < 90)).mean()), 3),
            "low_50_70": round(float(((v >= 50) & (v < 70)).mean()), 3),
            "very_low_lt50": round(float((v < 50).mean()), 3)}
