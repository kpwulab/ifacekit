"""Water and metal-site analysis.

Every water oxygen is described by the polar atoms around it and given one
category, plus any warning flags. Metal ions get a coordination table.

Categories (by number of polar partners within PARTNER_MAX, waters included):
  metal_bound    a metal ion within METAL_BOUND_MAX
  well_supported 3 or more partners
  supported      2 partners
  weak           1 partner
  unsupported    no polar partner at all (floating water, possibly noise)

Flags: too_close (a non-metal polar partner nearer than 2.4 A),
ww_close (another water nearer than 2.5 A), close_carbon (carbon nearer than
3.0 A), low_occupancy (occupancy below 1), high_b (B-factor above twice the
mean of its partners), partners_in_one_unit / bridging (partners from more than
one residue or unit). Thresholds are provisional.
"""
import numpy as np
from scipy.spatial import cKDTree

PARTNER_MAX = 3.5
METAL_BOUND_MAX = 2.6
TOO_CLOSE = 2.4
WW_CLOSE = 2.5
CARBON_CLOSE = 3.0
HIGH_B_RATIO = 2.0

# expected metal-oxygen/nitrogen distance windows (A) and coordination numbers
METAL_RULES = {
    "MG": dict(cn=(6,), dist=(1.9, 2.3)),
    "ZN": dict(cn=(4, 5, 6), dist=(1.9, 2.3)),
    "CA": dict(cn=(6, 7, 8), dist=(2.2, 2.7)),
    "MN": dict(cn=(6,), dist=(2.0, 2.4)),
}
METAL_LIG_MAX = {"MG": 2.6, "ZN": 2.6, "CA": 3.0, "MN": 2.7}


def _octahedral_deviation(center, lig):
    """RMS deviation (degrees) of the ligand angles from an ideal octahedron."""
    v = lig - center
    v /= np.linalg.norm(v, axis=1)[:, None]
    n = len(v)
    devs = []
    for i in range(n):
        ang = np.degrees(np.arccos(np.clip(v[i] @ v.T, -1, 1)))
        ang[i] = np.nan
        # each ligand should have one partner near 180 and four near 90
        a = np.sort(np.abs(ang[~np.isnan(ang)] ))
        ideal = np.array([90, 90, 90, 90, 180])
        if len(a) == 5:
            devs.extend(a - ideal)
    return float(np.sqrt(np.mean(np.square(devs)))) if devs else None


def analyze_waters(m):
    res_ids = [r["id"] for r in m.residues]
    aid = lambda k: f"{res_ids[m.res_idx[k]]}:{m.atom_name[k]}"
    wat = np.where(m.is_water & (m.element == "O"))[0]
    nonw = np.where(~m.is_water)[0]
    out = {"records": [], "metals": []}
    if len(wat) == 0 and not m.is_metal.any():
        return out
    el = m.element
    polar = (el == "N") | (el == "O")
    carbon = el == "C"
    metal = m.is_metal & ~m.is_water
    tree_all = cKDTree(m.xyz)

    rows = []
    for k in wat:
        near = [j for j in tree_all.query_ball_point(m.xyz[k], 3.6) if j != k and m.res_idx[j] != m.res_idx[k]]
        d = {j: float(np.linalg.norm(m.xyz[j] - m.xyz[k])) for j in near}
        mt = [j for j in near if metal[j] and d[j] <= METAL_BOUND_MAX]
        partners = [j for j in near if polar[j] and not metal[j] and d[j] <= PARTNER_MAX]
        nonwater_partners = [j for j in partners if not m.is_water[j]]
        water_partners = [j for j in partners if m.is_water[j]]
        if mt:
            cat = "metal_bound"
        else:
            n = len(partners)
            cat = "well_supported" if n >= 3 else "supported" if n == 2 else "weak" if n == 1 else "unsupported"
        flags = []
        if any(d[j] < TOO_CLOSE for j in nonwater_partners):
            flags.append("too_close")
        if any(d[j] < WW_CLOSE for j in water_partners):
            flags.append("ww_close")
        if any(carbon[j] and d[j] < CARBON_CLOSE for j in near):
            flags.append("close_carbon")
        if m.occupancy is not None and m.occupancy[k] < 0.999:
            flags.append("low_occupancy")
        if nonwater_partners:
            mb = float(np.mean([m.bfactor[j] for j in nonwater_partners]))
            if mb > 0 and m.bfactor[k] > HIGH_B_RATIO * mb:
                flags.append("high_b")
        res_set = {int(m.res_idx[j]) for j in nonwater_partners}
        unit_set = {int(m.unit_idx[j]) for j in nonwater_partners if m.unit_idx[j] >= 0}
        bridging = len(res_set) >= 2
        between_units = len(unit_set) >= 2
        rows.append({
                "water": res_ids[m.res_idx[k]], "category": cat, "flags": flags,
                "bridges_units": bool(between_units),
                "n_partners": len(partners),
                "partners": [{"atom": aid(j), "d": round(d[j], 2)} for j in sorted(partners, key=lambda x: d[x])][:6],
                "metals": [{"atom": aid(j), "d": round(d[j], 2)} for j in mt],
                "b": round(float(m.bfactor[k]), 1),
                "occupancy": round(float(m.occupancy[k]), 2) if m.occupancy is not None else None,
                "severity": 2 if ("too_close" in flags or "ww_close" in flags) else 1,
            })
    out["records"] = rows

    # ---- metal sites
    for k in np.where(metal)[0]:
        e = el[k]
        lim = METAL_LIG_MAX.get(e, 2.8)
        lig = [j for j in tree_all.query_ball_point(m.xyz[k], lim)
               if j != k and el[j] in ("O", "N", "S") and not metal[j]]
        dists = [float(np.linalg.norm(m.xyz[j] - m.xyz[k])) for j in lig]
        rule = METAL_RULES.get(e)
        verdict, notes = "info", []
        if rule:
            verdict = "ok"
            if len(lig) not in rule["cn"]:
                verdict = "borderline"
                notes.append(f"{len(lig)} ligands, expected {'/'.join(map(str, rule['cn']))}")
            out_of = [x for x in dists if not (rule["dist"][0] <= x <= rule["dist"][1])]
            if out_of:
                verdict = "borderline"
                notes.append(f"{len(out_of)} distance(s) outside {rule['dist'][0]}-{rule['dist'][1]} A")
        geom = None
        if e == "MG" and len(lig) == 6:
            geom = _octahedral_deviation(m.xyz[k].copy(), m.xyz[lig])
            if geom is not None and geom > 15:
                verdict = "borderline"
                notes.append(f"octahedral deviation {geom:.0f} deg")
        out["metals"].append({
            "metal": aid(k), "element": e, "n_ligands": len(lig), "verdict": verdict, "notes": notes,
            "octahedral_rms_deg": round(geom, 1) if geom is not None else None,
            "ligands": sorted(({"atom": aid(j), "water": bool(m.is_water[j]), "d": round(x, 2)}
                               for j, x in zip(lig, dists)), key=lambda r: r["d"]),
        })
    return out


def subset(full, water_ids, metal_ids, max_listed=200):
    """Summary of the waters and metals that belong to one analysis scope.

    `full` comes from analyze_waters on the whole model, so a water whose
    partners sit on a chain that is not part of the current scope is still
    judged with all its neighbours."""
    wset, mset = set(water_ids), set(metal_ids)
    recs = [r for r in full["records"] if r["water"] in wset]
    counts, flag_counts = {}, {}
    for r in recs:
        counts[r["category"]] = counts.get(r["category"], 0) + 1
        for f in r["flags"]:
            flag_counts[f] = flag_counts.get(f, 0) + 1
        if r["bridges_units"]:
            flag_counts["bridges_units"] = flag_counts.get("bridges_units", 0) + 1
    listed = [r for r in recs if r["flags"] or r["category"] == "unsupported"]
    listed.sort(key=lambda r: (-r["severity"], r["category"] != "unsupported", r["water"]))
    return {"n": len(recs), "counts": counts, "flag_counts": flag_counts,
            "listed": listed[:max_listed], "n_listed_total": len(listed),
            "metals": [x for x in full["metals"] if x["metal"].rsplit(":", 1)[0] in mset or x["metal"] in mset]}
