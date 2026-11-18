"""Compare ifacekit contact verdicts with the clashes MolProbity (phenix.molprobity) reports.

usage: python tools/compare_molprobity.py model.cif molprobity_coot.py
The Coot script written by phenix.molprobity carries the outlier lists as Python data.
Only the 'probe' list (severe clashes, overlap in A) and the rotamer, Ramachandran
and C-beta lists are read.
"""
import ast
import re
import sys

import ifacekit

ATOM = re.compile(r"\s*(\w)\s*(-?\d+[A-Za-z]?)\s+(\w{3})\s+(\S+)")


def read_data(path):
    data = {}
    for line in open(path):
        m = re.match(r"data\['(\w+)'\] = (.*)$", line.strip())
        if m:
            data[m.group(1)] = ast.literal_eval(m.group(2))
    return data


def parse_atom(s):
    """Fixed-width Probe atom label: chain(2) resnum(4) icode(1) altloc(1) resname(3) space atom."""
    return s[0:2].strip(), s[2:6].strip() + s[6].strip(), s[8:11], s[12:].strip()


def is_h(name):
    return name[0] == "H" or name[0] == "D"


def main():
    cif, coot = sys.argv[1], sys.argv[2]
    d = read_data(coot)
    probe = d.get("probe", [])
    mp = []
    for a, b, ov, _ in probe:
        pa, pb = parse_atom(a), parse_atom(b)
        mp.append(dict(a=pa, b=pb, overlap=-ov if ov < 0 else ov, h=is_h(pa[3]) or is_h(pb[3])))
    res = ifacekit.analyze(cif)
    ours = {}
    ours_res = {}
    for c in res["contacts"]:
        if c["verdict"] not in ("unfavorable", "borderline"):
            continue
        def key(s):
            ch, num, rn, an = s.split(":")
            return (ch, num, rn, an)
        ka, kb = key(c["atom1"]), key(c["atom2"])
        ours[frozenset([ka, kb])] = c
        ours_res.setdefault(frozenset([ka[:2], kb[:2]]), []).append(c)
    heavy = [x for x in mp if not x["h"]]
    hyd = [x for x in mp if x["h"]]
    print(f"MolProbity severe clashes: {len(mp)}  (heavy-atom pairs {len(heavy)}, involving hydrogen {len(hyd)})")
    print(f"ifacekit unfavorable+borderline pairs: {len(ours)}  "
          f"(unfavorable {sum(c['verdict'] == 'unfavorable' for c in ours.values())})")
    # Match on the residue pair, not on atom names: ring and carboxylate atoms are symmetric
    # (Phe CD1/CD2, Asp OD1/OD2) and MolProbity may name them the other way round.
    rp = lambda x: frozenset([x["a"][:2], x["b"][:2]])
    hit = [x for x in heavy if rp(x) in ours_res]
    print(f"heavy-atom MolProbity clashes whose residue pair ifacekit also flags: {len(hit)} of {len(heavy)}")
    hit_at = [x for x in heavy if frozenset([x["a"], x["b"]]) in ours]
    print(f"  (same atom pair, exact names: {len(hit_at)})")
    miss = [x for x in heavy if rp(x) not in ours_res]
    for x in miss[:20]:
        print("  missed by ifacekit:", ":".join(x["a"]), "-", ":".join(x["b"]), f"overlap {x['overlap']:.2f}")
    # residue level for hydrogen clashes
    ores = {(k[0], k[1]) for p in ours for k in p}
    hres = [x for x in hyd if (x["a"][0], x["a"][1]) in ores or (x["b"][0], x["b"][1]) in ores]
    print(f"hydrogen-involving MolProbity clashes where a residue is flagged by ifacekit anyway: {len(hres)} of {len(hyd)}")
    mp_res = {rp(x) for x in mp}
    only = [p for p, c in ours.items() if c["verdict"] == "unfavorable" and
            frozenset(k[:2] for k in p) not in mp_res]
    print(f"ifacekit unfavorable pairs whose residue pair is not in MolProbity's severe list: {len(only)}")
    from collections import Counter
    print("  by reason:", dict(Counter(ours[p]["reason"] for p in only)))
    for kind in ("rota", "rama", "cbeta"):
        print(f"MolProbity {kind} outliers: {len(d.get(kind, []))}")


main()
