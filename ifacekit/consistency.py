"""Compare metal sites between chains that have the same sequence.

For a homo-oligomer every copy of a site should look alike. Chains are
grouped by identical polymer sequence. A metal site is matched between chains
by its residue name and number. Its signature is the set of protein ligand
atoms plus the number of waters in the first shell. Chains whose signature
differs from the most common one are reported.
"""
from collections import Counter, defaultdict

from .model import load
from .waters import analyze_waters


def _strip(atom_id):
    return atom_id.split(":", 1)[1]


def metal_consistency(path, model_index=0, full=None):
    m = load(path, model_index)
    if full is None:
        full = analyze_waters(m)
    seqs = defaultdict(list)
    chain_res = defaultdict(list)
    for r in m.residues:
        if r["kind"] == "polymer":
            chain_res[r["chain"]].append(r["resname"])
    for ch, rn in chain_res.items():
        seqs[tuple(rn)].append(ch)
    groups = []
    for chains in seqs.values():
        if len(chains) < 2:
            continue
        sites = defaultdict(dict)
        for x in full["metals"]:
            ch, num, name = x["metal"].split(":")[:3]
            if ch not in chains:
                continue
            prot = tuple(sorted(_strip(l["atom"]) for l in x["ligands"] if not l["water"]))
            sites[(name, num)][ch] = dict(
                prot=prot, n_water=sum(l["water"] for l in x["ligands"]),
                n=x["n_ligands"], dmax=max((l["d"] for l in x["ligands"]), default=0),
                dmin=min((l["d"] for l in x["ligands"]), default=0),
                oct=x["octahedral_rms_deg"], water_ids=sorted(_strip(l["atom"]) for l in x["ligands"] if l["water"]))
        out_sites = []
        for (name, num), per in sorted(sites.items()):
            sig = Counter((v["prot"], v["n_water"]) for v in per.values())
            (cprot, cwat), cn = sig.most_common(1)[0]
            dev = []
            for ch, v in sorted(per.items()):
                if (v["prot"], v["n_water"]) != (cprot, cwat):
                    dev.append({"chain": ch, "n_ligands": v["n"], "n_water": v["n_water"],
                                "missing_protein_ligands": sorted(set(cprot) - set(v["prot"])),
                                "extra_protein_ligands": sorted(set(v["prot"]) - set(cprot))})
            octs = [v["oct"] for v in per.values() if v["oct"] is not None]
            out_sites.append({
                "site": f"{name}{num}", "chains_with_site": sorted(per),
                "missing_in": sorted(set(chains) - set(per)),
                "consensus": {"protein_ligands": list(cprot), "n_water": cwat, "n_chains": cn},
                "deviants": dev,
                "oct_dev_range": [round(min(octs), 1), round(max(octs), 1)] if octs else None,
                "dmin_dmax": [min(v["dmin"] for v in per.values()), max(v["dmax"] for v in per.values())]})
        groups.append({"chains": sorted(chains), "sites": out_sites})
    return {"groups": groups}


def consistency_text(c):
    if not c["groups"]:
        return ""
    lines = ["METAL SITE CONSISTENCY between chains with identical sequence"]
    for g in c["groups"]:
        lines.append(f"  chains {', '.join(g['chains'])}")
        for s in g["sites"]:
            con = s["consensus"]
            same = f"{con['n_chains']} of {len(s['chains_with_site'])} chains agree"
            lines.append(f"    {s['site']}: {same}; protein ligands {len(con['protein_ligands'])}, waters {con['n_water']}"
                         + (f"; octahedral deviation {s['oct_dev_range'][0]}-{s['oct_dev_range'][1]} deg" if s["oct_dev_range"] else ""))
            if s["missing_in"]:
                lines.append(f"      site missing in chains {', '.join(s['missing_in'])}")
            for d in s["deviants"]:
                bits = []
                if d["missing_protein_ligands"]:
                    bits.append("missing " + ", ".join(d["missing_protein_ligands"]))
                if d["extra_protein_ligands"]:
                    bits.append("extra " + ", ".join(d["extra_protein_ligands"]))
                bits.append(f"{d['n_water']} waters (consensus {con['n_water']})")
                lines.append(f"      chain {d['chain']}: {d['n_ligands']} ligands; " + "; ".join(bits))
    return "\n".join(lines)
