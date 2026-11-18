"""Top level analysis: load a model, build contacts, compute interfaces."""
import os
import time
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import cKDTree

from . import __version__
from .attach import plan_attachments
from .contacts import find_contacts
from .model import load
from .rings import find_pi
from .geometry import backbone_geometry
from . import rama as _rama
from . import flips as _flips
from . import hbonds as _hbonds
from . import rotamer as _rotamer
from .waters import analyze_waters, subset as water_subset
from .sasa import PROBE, atom_sasa

SCHEMA_VERSION = "0.12"
NEAR_CUTOFF = 6.5     # an atom farther than this from the partner cannot lose area
UNIT_PAIR_CUTOFF = 6.0
BURIED_DELTA_MIN = 1.0
BURIED_REMAIN_MAX = 0.5
_ORDER = {"ok": 0, "favorable": 0, "borderline": 1, "unfavorable": 2}
# Side-chain atoms whose identity is often swapped (amide or ring flip)
_FLIP_ATOMS = {("ASN", "OD1"), ("ASN", "ND2"), ("GLN", "OE1"), ("GLN", "NE2"),
               ("HIS", "ND1"), ("HIS", "NE2")}


def _atom_id(m, k):
    return f"{m.residues[m.res_idx[k]]['id']}:{m.atom_name[k]}"


def _unit_label(m, u):
    return "water" if u < 0 else m.units[u]["id"]


_WATER_CACHE = {}


def _full_waters(path, model_index, m):
    """Water and metal analysis on the whole model, computed once per file."""
    key = (str(path), model_index, os.path.getmtime(path))
    if key not in _WATER_CACHE:
        full = m if m.flags.get("scope") == "full" else load(path, model_index)
        _WATER_CACHE.clear()
        _WATER_CACHE[key] = analyze_waters(full)
    return _WATER_CACHE[key]


def list_chains(path, model_index=0):
    """One row per chain: what it holds. Used to help the user pick chains."""
    m = load(path, model_index, merge_by_chain=True)
    rows = {}
    for r in m.residues:
        row = rows.setdefault(r["chain"], {"chain": r["chain"], "polymer_type": None,
                                           "n_polymer_residues": 0, "n_atoms": 0,
                                           "ligands": [], "n_waters": 0})
        if r["kind"] == "polymer":
            row["n_polymer_residues"] += 1
        elif r["kind"] == "water":
            row["n_waters"] += 1
        elif r["resname"] not in row["ligands"]:
            row["ligands"].append(r["resname"])
    for u in m.units:
        if u["chain"] in rows:
            rows[u["chain"]]["polymer_type"] = u["polymer_type"]
    for row in rows.values():
        row["has_polymer"] = row["n_polymer_residues"] > 0
    for k in range(m.n_atoms):
        if not m.is_water[k]:
            rows[m.residues[m.res_idx[k]]["chain"]]["n_atoms"] += 1
    return list(rows.values())


def _rama_summary(path, model_index, chains):
    try:
        return _rama.analyze(path, model_index, chains)[0]
    except Exception as e:  # never let the extra check stop the main analysis
        return {"error": str(e)}


def _rotamer_summary(path, model_index, chains):
    try:
        return _rotamer.analyze(path, model_index, chains)[0]
    except Exception as e:  # never let the extra check stop the main analysis
        return {"error": str(e)}


def _flips_summary(path, model_index, chains):
    try:
        return _flips.analyze(path, model_index, chains)[0]
    except Exception as e:
        return {"error": str(e)}


def _contact_records(m, cs, altloc=None):
    out = []
    for r in range(len(cs)):
        i, j = int(cs.i[r]), int(cs.j[r])
        c = {
            "atom1": _atom_id(m, i), "atom2": _atom_id(m, j),
            "distance": round(float(cs.d[r]), 2),
            "overlap": round(float(cs.overlap[r]), 2),
            "type": cs.type[r], "verdict": cs.verdict[r], "reason": cs.reason[r],
            "hbond": bool(cs.hbond[r]),
            "hbond_geometry": bool(cs.hbond_geometry[r]) if cs.hbond_geometry is not None and len(cs.hbond_geometry) else False,
            "unit1": _unit_label(m, int(m.unit_idx[i])),
            "unit2": _unit_label(m, int(m.unit_idx[j])),
            "inter_unit": bool(cs.inter_unit[r]),
            "b1": round(float(m.bfactor[i]), 1), "b2": round(float(m.bfactor[j]), 1),
        }
        if altloc:
            c["altloc"] = altloc
        out.append(c)
    return out


def _issue_items(m, cs, mode, altloc=None):
    """Residue id -> list of borderline / unfavorable items, in the order of the residues in the file."""
    by = defaultdict(list)
    for r in range(len(cs)):
        if cs.verdict[r] not in ("borderline", "unfavorable"):
            continue
        if mode == "pair" and not cs.inter_unit[r]:
            continue          # binary mode reports interactions between the chains only
        i, j = int(cs.i[r]), int(cs.j[r])
        flip = cs.reason[r] == "polar_mismatch" and any(
            (m.residues[m.res_idx[k]]["resname"], m.atom_name[k]) in _FLIP_ATOMS for k in (i, j))
        for me, other in ((i, j), (j, i)):
            item = {"atom": m.atom_name[me], "with": _atom_id(m, other),
                    "type": cs.type[r], "reason": cs.reason[r],
                    "verdict": cs.verdict[r], "distance": round(float(cs.d[r]), 2),
                    "inter_unit": bool(cs.inter_unit[r]), "possible_flip": bool(flip)}
            if altloc:
                item["altloc"] = altloc
            by[int(m.res_idx[me])].append(item)
    return {m.residues[rr]["id"]: items for rr, items in sorted(by.items())}


def _hydrogen_block(m, mode):
    """All-atom clash check on a copy of the model with hydrogens placed (the loaded model itself is not changed)."""
    from . import hydrogens as hy
    t0 = time.time()
    H, env, groups = hy.add_hydrogens(m)
    n = env.n

    def label(k):
        if k < n:
            return _atom_id(m, k), int(m.unit_idx[k]), float(m.bfactor[k])
        h = k - n
        r = int(H.res[h])
        return f"{m.residues[r]['id']}:{H.name[h]}", int(m.residues[r]["unit"]), float(m.bfactor[int(H.parent[h])])

    rows, seen = [], set()
    for x in hy.clashes(m, H, env, heavy=True):
        a1, u1, b1 = label(x["i"])
        a2, u2, b2 = label(x["j"])
        key = frozenset([a1, a2])
        if key in seen:
            continue
        seen.add(key)
        inter = not (u1 == u2 and u1 >= 0)
        if mode == "pair" and not inter:
            continue
        rows.append({"atom1": a1, "atom2": a2, "distance": round(x["d"], 2), "overlap": round(x["overlap"], 2),
                     "unit1": _unit_label(m, u1), "unit2": _unit_label(m, u2), "inter_unit": bool(inter),
                     "hydrogen": bool(x["i"] >= n or x["j"] >= n), "b1": round(b1, 1), "b2": round(b2, 1)})
    rows.sort(key=lambda r: -r["overlap"])
    flips = []
    for g in groups:
        if g.kind == "flip" and g.state == 1:
            flips.append({"residue": m.residues[g.res]["id"], "what": "amide flipped"})
        elif g.kind == "his":
            fl, prot = g.state
            flips.append({"residue": m.residues[g.res]["id"],
                          "what": ("ring flipped, " if fl else "") + ("HID: H on ND1", "HIE: H on NE2", "HIP: both protonated")[prot]})
    n_all = int((~m.is_water).sum()) + int(H.active.sum()) + int(m.is_water.sum())
    return {"n_hydrogens": int(H.active.sum()), "n_clashes": len(rows),
            "n_clashes_with_hydrogen": sum(1 for r in rows if r["hydrogen"]),
            "clashscore": round(1000.0 * len(rows) / max(n_all, 1), 2),
            "flips": flips, "clashes": rows, "seconds": round(time.time() - t0, 1),
            "limits": {"min_overlap": hy.OVERLAP_CLASH, "h_bond_overlap_allowed": hy.HB_CAP, "h_bond_overlap_allowed_charged": hy.HB_CAP_CHARGED,
                       "polar_heavy_overlap_allowed": hy.HEAVY_POLAR_CAP, "flip_cost": hy.FLIP_COST}}


def analyze(path, model_index=0, min_area=10.0, n_points=960,
            record_intra="problems", progress=None, chain=None, pair=None,
            attach="contact", plan=None, alt_conformers=True, hydrogens=False):
    """Analyse a model.

    chain   one chain name: quality inspection of that chain alone (its own
            ligands and waters included, nothing else in the file)
    pair    (chain_a, chain_b): interface analysis between exactly two chains
    neither: every unit in the file (ligands as separate units)

    attach  "contact" (default): ligands, glycans and waters join the polymer
            chain they are bonded to or touch, whatever their chain label says;
            "label": only residues carrying the chain's own label join it
    plan    result of attach.plan_attachments, to avoid recomputing it
    alt_conformers  also check the contacts of every alternate conformer (default); the first conformer is always checked
    hydrogens  also check the model for clashes with hydrogens added (see ifacekit.hydrogens). Only this check uses the
            hydrogens; every other result comes from the model as given
    """
    t0 = time.time()
    log = progress or (lambda msg: None)
    mode, select, merge, chains, unit_chain, owners = "all", None, False, None, None, {}
    if chain is not None and pair is not None:
        raise ValueError("give either chain or pair, not both")
    if chain is not None or pair is not None:
        chains = [chain] if chain is not None else list(pair)
        if pair is not None and (len(chains) != 2 or chains[0] == chains[1]):
            raise ValueError("pair needs two different chains")
        rows = list_chains(path, model_index)
        have = [r["chain"] for r in rows]
        missing = [c for c in chains if c not in have]
        if missing:
            raise ValueError(f"chain {', '.join(missing)} not found; chains in file: {', '.join(have)}")
        no_poly = [r["chain"] for r in rows if r["chain"] in chains and not r["has_polymer"]]
        if no_poly:
            raise ValueError(
                f"chain {', '.join(no_poly)} holds only ligands or waters. They are attached to the "
                "polymer chains they touch, so inspect those chains instead. Polymer chains: "
                + ", ".join(r["chain"] for r in rows if r["has_polymer"]))
        merge = True
        if attach == "contact":
            owners = plan if plan is not None else plan_attachments(path, model_index)
        else:
            owners = {}
        if chain is not None:
            mode = "quality"
            def select(ch, kind, rid):
                if kind in ("polymer",):
                    return ch == chain
                return ch == chain or chain in owners.get(rid, ())
            unit_chain = lambda ch, kind, rid: chain
        else:
            mode = "pair"
            def _owner(rid, ch):
                cand = [c for c in owners.get(rid, ()) if c in chains]
                return cand[0] if cand else (ch if ch in chains else None)
            def select(ch, kind, rid):
                if kind == "polymer":
                    return ch in chains
                if kind == "water":
                    return ch in chains or any(c in chains for c in owners.get(rid, ()))
                return _owner(rid, ch) is not None
            unit_chain = lambda ch, kind, rid: _owner(rid, ch) or ch
    m = load(path, model_index, select=select, merge_by_chain=merge, unit_chain=unit_chain)
    if mode == "pair":
        shared = []
        for r in m.residues:
            if r["kind"] in ("ligand", "branched"):
                touch = [c for c in owners.get(r["id"], ()) if c in chains]
                if len(touch) > 1:
                    shared.append({"residue": r["id"], "touches": touch, "assigned_to": m.units[r["unit"]]["id"]})
        m.flags["ligands_touching_both_chains"] = shared
    log(f"loaded {m.n_atoms} atoms, {len(m.units)} units")

    cs = find_contacts(m, record_intra=record_intra)
    log(f"{len(cs)} contact records")
    pi = find_pi(m)

    # ---- candidate unit pairs
    nw = np.where(m.unit_idx >= 0)[0]
    pair_set = set()
    if len(nw) > 1:
        tree = cKDTree(m.xyz[nw])
        pr = tree.query_pairs(r=UNIT_PAIR_CUTOFF, output_type="ndarray")
        if len(pr):
            ua = m.unit_idx[nw[pr[:, 0]]]
            ub = m.unit_idx[nw[pr[:, 1]]]
            diff = ua != ub
            lo = np.minimum(ua[diff], ub[diff])
            hi = np.maximum(ua[diff], ub[diff])
            pair_set = set(zip(lo.tolist(), hi.tolist()))

    # ---- SASA of each unit on its own (only units that touch another unit)
    alone = np.zeros(m.n_atoms)
    unit_atoms = {u: np.where(m.unit_idx == u)[0] for u in range(len(m.units))}
    for u in sorted({x for pr_ in pair_set for x in pr_}):
        idx = unit_atoms[u]
        alone[idx] = atom_sasa(m.xyz[idx], m.radius[idx], n_points=n_points)
    log("unit SASA done")

    # ---- index contact records by unit pair
    ui_arr = m.unit_idx[cs.i] if len(cs) else np.zeros(0, int)
    uj_arr = m.unit_idx[cs.j] if len(cs) else np.zeros(0, int)
    by_pair = defaultdict(list)
    water_units = defaultdict(set)
    for r in range(len(cs)):
        a, b = int(ui_arr[r]), int(uj_arr[r])
        if a >= 0 and b >= 0 and a != b:
            by_pair[(min(a, b), max(a, b))].append(r)
        elif (a < 0) != (b < 0) and cs.hbond[r]:
            widx = cs.i[r] if a < 0 else cs.j[r]
            water_units[int(m.res_idx[widx])].add(b if a < 0 else a)

    polar_el = np.isin(m.element, ["N", "O"]) | m.is_metal
    needs_partner = (m.donor | m.acceptor) & ~m.is_metal

    interfaces = []
    for (a, b) in sorted(pair_set):
        ia, ib = unit_atoms[a], unit_atoms[b]
        tb = cKDTree(m.xyz[ib])
        ta = cKDTree(m.xyz[ia])
        near_a = ia[tb.query(m.xyz[ia])[0] < NEAR_CUTOFF]
        near_b = ib[ta.query(m.xyz[ib])[0] < NEAR_CUTOFF]
        env = np.concatenate([ia, ib])
        pos = {int(k): p for p, k in enumerate(env)}
        near = np.concatenate([near_a, near_b])
        targets = np.array([pos[int(k)] for k in near])
        pair_sasa = atom_sasa(m.xyz[env], m.radius[env], targets=targets, n_points=n_points)
        delta = np.clip(alone[near] - pair_sasa, 0.0, None)
        remain = pair_sasa
        area_a = float(delta[: len(near_a)].sum())
        area_b = float(delta[len(near_a):].sum())
        area = 0.5 * (area_a + area_b)
        if area < min_area:
            continue
        pol = float(delta[polar_el[near]].sum())
        tot = float(delta.sum())

        res_bsa = defaultdict(float)
        for k, dv in zip(near, delta):
            res_bsa[int(m.res_idx[k])] += float(dv)

        rec_idx = by_pair.get((a, b), [])
        type_counts = Counter(cs.type[r] for r in rec_idx)
        verdict_counts = Counter(cs.verdict[r] for r in rec_idx)

        res_contacts = defaultdict(lambda: [0, 0])
        for r in rec_idx:
            for k in (cs.i[r], cs.j[r]):
                rr = int(m.res_idx[k])
                res_contacts[rr][0] += 1
                if cs.verdict[r] in ("borderline", "unfavorable"):
                    res_contacts[rr][1] += 1

        res_b = {}
        for k in range(m.n_atoms):
            res_b.setdefault(int(m.res_idx[k]), []).append(m.bfactor[k])

        def side(unit):
            rows = []
            for rr, v in sorted(res_bsa.items()):
                if m.residues[rr]["unit"] != unit or v < 0.5:
                    continue
                nc, npb = res_contacts.get(rr, (0, 0))
                rows.append({"residue": m.residues[rr]["id"], "bsa": round(v, 1),
                             "contacts": nc, "problem_contacts": npb,
                             "mean_b": round(float(np.mean(res_b[rr])), 1)})
            return rows

        buns = []
        for k, dv, rem in zip(near, delta, remain):
            if needs_partner[k] and dv >= BURIED_DELTA_MIN and rem <= BURIED_REMAIN_MAX \
                    and not cs.satisfied[k]:
                buns.append({"atom": _atom_id(m, k), "buried": round(float(dv), 1),
                             "role_uncertain": bool(m.ambiguous[k])})

        bridging = sorted(
            m.residues[w]["id"] for w, us in water_units.items() if a in us and b in us)

        interfaces.append({
            "unit_a": m.units[a]["id"], "unit_b": m.units[b]["id"],
            "kind_a": m.units[a]["kind"], "kind_b": m.units[b]["kind"],
            "area": round(area, 1), "area_a": round(area_a, 1), "area_b": round(area_b, 1),
            "polar_fraction": round(pol / tot, 3) if tot > 0 else None,
            "n_residues_a": len(side(a)), "n_residues_b": len(side(b)),
            "contacts": dict(type_counts), "verdicts": dict(verdict_counts),
            "buried_unsatisfied_polar": buns,
            "bridging_waters": bridging,
            "residues_a": side(a), "residues_b": side(b),
        })
        for key in ("a", "b"):
            bs = [x["mean_b"] for x in interfaces[-1]["residues_" + key]]
            interfaces[-1]["mean_b_" + key] = round(float(np.mean(bs)), 1) if bs else None
            interfaces[-1]["n_low_b_" + key] = int(sum(1 for x in bs if x < 70))
    interfaces.sort(key=lambda x: -x["area"])
    log(f"{len(interfaces)} interfaces")

    # ---- contact records for output
    contacts = _contact_records(m, cs)

    # ---- residue level issues
    issues_by_id = _issue_items(m, cs, mode)

    # ---- contacts of alternate conformers: scan every other conformer and add the contacts that involve its own atoms
    counts = {k: dict(v) for k, v in cs.counts.items()}
    letters = m.flags.get("altloc_letters", []) if alt_conformers else []
    scanned = []
    for L in letters:
        m2 = load(path, model_index, select=select, merge_by_chain=merge, unit_chain=unit_chain, conformer=L)
        if m2.alt_extra is None or not m2.alt_extra.any():
            continue
        cs2 = find_contacts(m2, record_intra=record_intra, restrict=m2.alt_extra)
        contacts += _contact_records(m2, cs2, altloc=L)
        for rid, items in _issue_items(m2, cs2, mode, altloc=L).items():
            issues_by_id.setdefault(rid, []).extend(items)
        for k, v in cs2.counts.items():
            for kk, n in v.items():
                counts[k][kk] = counts[k].get(kk, 0) + n
        scanned.append(L)
    m.flags["alt_conformers_scanned"] = scanned
    issues = [{"residue": rid, "worst": max((x["verdict"] for x in items), key=lambda v: _ORDER[v]), "items": items}
              for rid, items in issues_by_id.items()]

    xyz_of = _hbonds.xyz_lookup(m, _atom_id)
    for itf in interfaces:                       # H-bonds without and with water (see ifacekit.hbonds)
        ua, ub = itf["unit_a"], itf["unit_b"]
        itf["hbonds"] = _hbonds.split(contacts, lambda u, ua=ua: u == ua, lambda u, ub=ub: u == ub, xyz=xyz_of)

    hyd = _hydrogen_block(m, mode) if hydrogens else None
    if hyd:
        log(f"{hyd['n_hydrogens']} hydrogens placed, {hyd['n_clashes']} all-atom clashes")

    all_verdicts = counts["verdicts"]
    n_heavy = int((~m.is_water).sum())
    unfav = all_verdicts.get("unfavorable", 0)
    summary = {
        "n_atoms": m.n_atoms, "n_non_water_atoms": n_heavy,
        "n_units": len(m.units), "n_interfaces": len(interfaces),
        "contact_records": len(contacts), "pairs_evaluated": sum(all_verdicts.values()),
        "verdicts": dict(all_verdicts),
        "verdicts_between_units": counts["verdicts_between_units"],
        "verdicts_within_units": counts["verdicts_within_units"],
        "types": counts["types"],
        "alt_conformer_contacts": sum(1 for c in contacts if c.get("altloc")),
        "unfavorable_per_1000_atoms": round(1000.0 * unfav / max(n_heavy, 1), 2),
        "note": "heavy-atom engine; not comparable to MolProbity clashscore",
    }
    if hyd:
        summary["all_atom_clashes"] = hyd["n_clashes"]
        summary["all_atom_clashscore"] = hyd["clashscore"]
    return {
        "schema_version": SCHEMA_VERSION,
        "tool": {"name": "ifacekit", "version": __version__},
        "input": {"path": m.path, "format": m.fmt, "model_index": m.model_index,
                  "mode": mode, "chains": chains, "attach": attach if mode != "all" else None,
                  "probe_radius": PROBE, "sphere_points": n_points},
        "flags": m.flags,
        "units": m.units,
        "interfaces": interfaces,
        "pi_interactions": pi,
        "ramachandran": _rama_summary(path, model_index, chains if mode != "all" else None),
        "rotamer": _rotamer_summary(path, model_index, chains if mode != "all" else None),
        "flips": _flips_summary(path, model_index, chains if mode != "all" else None),
        "geometry": backbone_geometry(m, {r["id"] for r in m.residues if r["protein"]}),
        "waters": water_subset(_full_waters(path, model_index, m),
                               [r["id"] for r in m.residues if r["kind"] == "water"],
                               [r["id"] for r in m.residues if r["kind"] != "water"]),
        "hydrogens": hyd,
        "contacts": contacts,
        "issues": issues,
        "summary": summary,
        "runtime_seconds": round(time.time() - t0, 2),
    }
