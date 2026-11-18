#!/usr/bin/env python3
"""Sort the clash pairs on which ifacekit (with hydrogens added) and MolProbity disagree by the most likely cause.

Same input and scope as tools/compare_hydrogen_clashes.py (first conformer, residue pairs). For every pair that only MolProbity has, and every pair
that only ifacekit has, one cause is assigned, in this order (the rules are listed in the code; they are attributions, not proofs):
  MolProbity only:
    hydrogen of a residue ifacekit gives no hydrogens (water, ligand, nucleic acid, modified residue)
    left out by the bond-count rule (hydrogens of one residue or neighbours, parent atoms three bonds apart or fewer)
    called a hydrogen bond by ifacekit (polar hydrogen next to an acceptor at a hydrogen-bond angle)
    two N/O atoms below the short hydrogen-bond overlap
    overlap just below the limit in ifacekit (0.30 to 0.40 A): the position of a hydrogen differs a little
    no or small overlap in ifacekit: the hydrogen is somewhere else (turned group, flip, protonation, or another conformer choice)
    two heavy atoms, one of them in a water or ligand / both in the protein (no overlap of 0.3 A or more in ifacekit)
    other (a rule of ifacekit not named above)
  ifacekit only:
    two heavy atoms (with a water / with a ligand or ion / protein only)
    overlap of 0.40 to 0.45 A (just above the limit)
    a hydrogen whose position is a choice (hydroxyl, thiol, NH3, amide, His ring)
    a hydrogen with fixed geometry (C-H, backbone N-H)
Usage: python tools/explain_hydrogen_disagreement.py CASES.txt MP_DIR OUT.json [--jobs 4]
"""
import argparse
import json
import multiprocessing as mp
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_hydrogen_clashes as chc  # noqa: E402
import compare_molprobity_batch as cm  # noqa: E402

from ifacekit import chem  # noqa: E402
from ifacekit import hydrogens as hy  # noqa: E402
from ifacekit import model  # noqa: E402

CHOICE_H = {("SER", "HG"), ("THR", "HG1"), ("CYS", "HG"), ("TYR", "HH"), ("LYS", "HZ1"), ("LYS", "HZ2"), ("LYS", "HZ3"),
            ("ASN", "HD21"), ("ASN", "HD22"), ("GLN", "HE21"), ("GLN", "HE22"), ("HIS", "HD1"), ("HIS", "HE2")}
MP_ONLY = ["ligand_water_hydrogen", "bond_rule", "hydrogen_bond_rule", "short_polar_rule", "near_threshold", "placement", "heavy_water_ligand", "heavy_protein", "other"]
OURS_ONLY = ["heavy_water", "heavy_ligand", "heavy_protein", "near_threshold", "choice_hydrogen", "fixed_hydrogen"]


def one(args):
    cid, mp_dir = args
    d = os.path.join(mp_dir, cid)
    if not os.path.exists(os.path.join(d, "molprobity_coot.py")):
        return None
    m = model.load(os.path.join(d, "m.cif"))
    H, env, _ = hy.add_hydrogens(m)
    n = env.n

    def rkey(k):
        r = m.residues[m.res_idx[k]] if k < n else m.residues[H.res[k - n]]
        return (r["chain"], str(r["resnum"]) + r["icode"])

    def rname(k):
        return m.residues[m.res_idx[k]]["resname"] if k < n else m.residues[H.res[k - n]]["resname"]

    final = {}
    for x in hy.clashes(m, H, env, heavy=True):
        k = frozenset([rkey(x["i"]), rkey(x["j"])])
        if k not in final or x["overlap"] > final[k]["overlap"]:
            final[k] = x
    raw = {}
    for x in hy.clashes(m, H, env, heavy=True, raw=True):
        raw.setdefault(frozenset([rkey(x["i"]), rkey(x["j"])]), []).append(x)
    # MolProbity atoms per residue pair (first conformer only), strongest first
    coot = cm.read_coot(os.path.join(d, "molprobity_coot.py"))
    fa = chc.first_conformers(os.path.join(d, "m.cif"))
    mpa = {}
    for a, b, o, _ in coot.get("probe", []):
        ra, rb = cm.probe_atom(a), cm.probe_atom(b)
        if any(lab[7] != " " and fa.get(r[0]) != lab[7] for lab, r in ((a, ra), (b, rb))):
            continue
        mpa.setdefault(frozenset([(ra[0][0], ra[0][1]), (rb[0][0], rb[0][1])]), []).append((ra, rb, -o))
    out = {"mp": {c: 0 for c in MP_ONLY}, "ours": {c: 0 for c in OURS_ONLY}, "shared": 0, "mp_total": len(mpa), "ours_total": len(final)}
    for k, atoms in mpa.items():
        if k in final:
            out["shared"] += 1
            continue
        atoms.sort(key=lambda t: -t[2])
        ra, rb, _ = atoms[0]
        if any(cm.is_h(r[2]) and r[1] not in chem.STANDARD_AA for r in (ra, rb)):
            cat = "ligand_water_hydrogen"
        else:
            R = [x for x in raw.get(k, []) if x["overlap"] >= hy.OVERLAP_CLASH]
            if R:                                    # every one of these atom pairs was removed by a rule; name the rule of the largest overlap
                x = max(R, key=lambda t: t["overlap"])
                if x.get("heavy"):
                    cat = "short_polar_rule"
                elif x.get("hb") and x["overlap"] < (hy.HB_CAP_CHARGED if x.get("charged") else hy.HB_CAP):
                    cat = "hydrogen_bond_rule"
                elif x["hop"] <= (hy.HH_HOPS if x["hh"] else 3):
                    cat = "bond_rule"
                else:
                    cat = "other"
            else:
                best = max([x["overlap"] for x in raw.get(k, [])] + [0.0])
                if best >= 0.30:
                    cat = "near_threshold"
                elif any(cm.is_h(r[2]) for r in (ra, rb)):
                    cat = "placement"
                elif any(r[1] not in chem.STANDARD_AA for r in (ra, rb)):
                    cat = "heavy_water_ligand"
                else:
                    cat = "heavy_protein"
        out["mp"][cat] += 1
    for k, x in final.items():
        if k in mpa:
            continue
        i, j = x["i"], x["j"]
        if i < n and j < n:
            water = bool(m.is_water[i] or m.is_water[j])
            lig = any(m.residues[m.res_idx[t]]["kind"] != "polymer" for t in (i, j))
            cat = "heavy_water" if water else "heavy_ligand" if lig else "heavy_protein"
        elif x["overlap"] < hy.OVERLAP_CLASH + 0.05:
            cat = "near_threshold"
        else:
            hs = [t for t in (i, j) if t >= n]
            cat = "choice_hydrogen" if any((rname(t), H.name[t - n]) in CHOICE_H or H.name[t - n] in ("H1", "H2", "H3") for t in hs) else "fixed_hydrogen"
        out["ours"][cat] += 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cases")
    ap.add_argument("mp_dir")
    ap.add_argument("out")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    ids = open(a.cases).read().split()
    with mp.Pool(a.jobs) as pool:
        rows = [r for r in pool.map(one, [(c, a.mp_dir) for c in ids], chunksize=2) if r]
    tot = {"cases": len(rows), "shared": sum(r["shared"] for r in rows), "mp_total": sum(r["mp_total"] for r in rows), "ours_total": sum(r["ours_total"] for r in rows),
           "mp_only": {c: sum(r["mp"][c] for r in rows) for c in MP_ONLY}, "ours_only": {c: sum(r["ours"][c] for r in rows) for c in OURS_ONLY}}
    json.dump(tot, open(a.out, "w"), indent=1)
    print(json.dumps(tot, indent=1))


if __name__ == "__main__":
    main()
