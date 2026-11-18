"""Plain text reports. One layout per mode: chain quality, chain pair, whole file."""
from collections import Counter

_VERDICT_ORDER = ["favorable", "ok", "borderline", "unfavorable"]


def _notes(res):
    fl = res["flags"]
    out = []
    if fl.get("hydrogens_removed"):
        out.append(f"{fl['hydrogens_removed']} hydrogens ignored (heavy-atom analysis)")
    if fl.get("altloc_atoms_dropped"):
        out.append(f"{fl['altloc_atoms_dropped']} alternate-conformer atoms dropped (first conformer kept)")
    if fl.get("nonstandard_residues"):
        out.append("non-standard residues: " + ", ".join(fl["nonstandard_residues"]))
    if fl.get("unknown_ligand_roles"):
        out.append("some ligand or modified atoms have guessed donor/acceptor roles")
    if fl.get("n_models_in_file", 1) > 1:
        if fl.get("ensemble"):
            out.append(f"{fl['n_models_in_file']} models in file; this is the first model, the numbers of every model are in the ensemble table")
        else:
            out.append(f"{fl['n_models_in_file']} models in file, only model {res['input']['model_index']} analysed (use --ensemble for every model)")
    return out


def _head(res, title):
    inp = res["input"]
    return [f"ifacekit {res['tool']['version']}   {inp['path']}  ({inp['format']})",
            title,
            f"runtime {res['runtime_seconds']} s"] + ["note: " + n for n in _notes(res)]


def _tally_line(label, v):
    return "{}  favorable {}  ok {}  borderline {}  unfavorable {}".format(
        label, v.get("favorable", 0), v.get("ok", 0), v.get("borderline", 0), v.get("unfavorable", 0))


def _table(rows, header):
    widths = [max(len(str(x)) for x in col) for col in zip(header, *rows)] if rows else [len(h) for h in header]
    fmt = "  ".join("{:<%d}" % w for w in widths)
    out = [fmt.format(*header)]
    out += [fmt.format(*[str(x) for x in r]) for r in rows]
    return out


def _geometry_lines(res, max_rows=8):
    g = res.get("geometry") or {}
    if not g or not g.get("n_bonds"):
        return []
    lines = ["", "BACKBONE GEOMETRY (approximate Engh-Huber values, |Z| > 4 listed)",
             f"  bonds {g['n_bonds']} rmsZ {g['rmsz_bonds']};  angles {g['n_angles']} rmsZ {g['rmsz_angles']};  outliers {g['n_outliers']}"]
    rows = [[o["residue"], o["kind"], o["atoms"], o["value"], o["ideal"], o["z"]] for o in g["outliers"][:max_rows]]
    if rows:
        lines += ["  " + x for x in _table(rows, ["residue", "kind", "atoms", "value", "ideal", "Z"])]
    st = [b for b in g["chain_breaks"] if b["kind"] == "stretched_bond"]
    gaps = [b for b in g["chain_breaks"] if b["kind"] == "gap_in_numbering"]
    if st:
        lines.append("  stretched peptide bonds: " + ", ".join(f"{b['after']}-{b['before'].split(':')[1]} ({b['C_N_distance']} A)" for b in st[:8]))
    if gaps:
        lines.append(f"  gaps in numbering (missing residues): {len(gaps)}")
    return lines


def _rama_lines(res, max_rows=12):
    g = res.get("ramachandran") or {}
    if not g.get("n"):
        return []
    from .rama import LABEL
    p, c = g["percent"], g["counts"]
    lines = ["", f"RAMACHANDRAN (ifacekit reference, {g['n']} residues with phi and psi)",
             f"  favored {c['favored']} ({p['favored']}%)  allowed {c['allowed']} ({p['allowed']}%)  "
             f"outliers {c['outlier']} ({p['outlier']}%)   goals: favored > 98%, outliers < 0.2%"]
    rows = [[o["residue"], LABEL[o["class"]], o["phi"], o["psi"]] for o in g["outliers"][:max_rows]]
    if rows:
        lines += ["  " + x for x in _table(rows, ["residue", "class", "phi", "psi"])]
    if len(g["outliers"]) > max_rows:
        lines.append(f"  ... {len(g['outliers']) - max_rows} more (see JSON)")
    return lines


def _rotamer_lines(res, max_rows=12):
    g = res.get("rotamer") or {}
    if not g.get("n"):
        return []
    from .rotamer import LABEL
    p, c = g["percent"], g["counts"]
    lines = ["", f"ROTAMERS (ifacekit reference, {g['n']} side chains with chi angles)",
             f"  favored {c['favored']} ({p['favored']}%)  allowed {c['allowed']} ({p['allowed']}%)  "
             f"outliers {c['outlier']} ({p['outlier']}%)   goal: outliers < 0.3%"]
    from .rotamer import label as _rlabel
    rows = [[_rlabel(o), LABEL[o["table"]], " ".join(str(x) for x in o["chi"])] for o in g["outliers"][:max_rows]]
    if rows:
        lines += ["  " + x for x in _table(rows, ["residue", "type", "chi angles"])]
    if len(g["outliers"]) > max_rows:
        lines.append(f"  ... {len(g['outliers']) - max_rows} more (see JSON)")
    af = g.get("after_suggested_flips") or {}
    if af.get("n_flips_applied"):
        lines.append(f"  after the {af['n_flips_applied']} flip(s) suggested below: favored {af['percent']['favored']}%, "
                     f"outliers {af['counts']['outlier']} ({af['percent']['outlier']}%)")
    return lines


def _flip_lines(res, max_rows=12):
    f = res.get("flips") or {}
    if not f.get("n"):
        return []
    c = f["counts"]
    lines = ["", f"ASN / GLN / HIS FLIPS ({f['n']} groups, heavy-atom hydrogen-bond and overlap score, no hydrogens)",
             f"  flip suggested {c['flip']}, uncertain {c['uncertain']}, keep {c['keep']}"]
    rows = [[x["residue"], x["score_modelled"], x["score_flipped"], x["delta"]] for x in f["suggested"][:max_rows]]
    if rows:
        lines += ["  " + x for x in _table(rows, ["residue", "score as modelled", "score flipped", "gain"])]
    if len(f["suggested"]) > max_rows:
        lines.append(f"  ... {len(f['suggested']) - max_rows} more (see JSON)")
    lines.append("  a suggestion, not a decision: the check is weaker than Reduce (see the design document); look at the density")
    return lines


def _water_lines(res, max_rows=12):
    w = res.get("waters") or {}
    lines = []
    if not w or (not w.get("n") and not w.get("metals")):
        return lines
    lines += ["", "WATERS AND METAL SITES"]
    if w["n"]:
        c = w["counts"]
        lines.append(f"  waters {w['n']}: " + ", ".join(
            f"{k} {c[k]}" for k in ("metal_bound", "well_supported", "supported", "weak", "unsupported") if c.get(k)))
        fc = w["flag_counts"]
        parts = [f"{k} {fc[k]}" for k in ("too_close", "ww_close", "close_carbon", "high_b", "low_occupancy") if fc.get(k)]
        lines.append("  flags: " + (", ".join(parts) if parts else "none")
                     + (f"; waters touching two or more units: {fc['bridges_units']}" if fc.get("bridges_units") else ""))
        lines.append("  partner = N or O within 3.5 A (waters included); 'weak' = one partner, 'unsupported' = none")
        rows = []
        for r in w["listed"][:max_rows]:
            near = ", ".join(f"{p['atom'].split(':', 1)[1]} {p['d']}" for p in r["partners"][:2]) or "-"
            rows.append([r["water"], r["category"], ",".join(r["flags"]) or "-", near, r["b"]])
        if rows:
            lines.append(f"  waters to check ({w.get('n_listed_total', len(w['listed']))} flagged, shown {len(rows)}):")
            lines += ["    " + x for x in _table(rows, ["water", "category", "flags", "nearest partners", "B"])]
    ms = w.get("metals") or []
    if ms:
        bad = [x for x in ms if x["verdict"] == "borderline"]
        lines.append(f"  metal sites {len(ms)}: ok or informational {len(ms) - len(bad)}, to check {len(bad)}")
        rows = []
        for x in ms:
            if x["verdict"] == "borderline" or len(ms) <= 4:
                ds = [l["d"] for l in x["ligands"]]
                rows.append([x["metal"], x["n_ligands"], f"{min(ds)}-{max(ds)}" if ds else "-",
                             sum(l["water"] for l in x["ligands"]),
                             x["octahedral_rms_deg"] if x["octahedral_rms_deg"] is not None else "-",
                             "; ".join(x["notes"]) or x["verdict"]])
        if rows:
            lines += ["    " + x for x in _table(rows, ["metal", "ligands", "d range", "waters", "oct dev", "note"])]
    return lines


# ---------------------------------------------------------------- quality
def _hydrogen_lines(res, units=None, max_rows=10):
    """All-atom clashes with hydrogens added (only when the analysis was run with hydrogens)."""
    h = res.get("hydrogens")
    if not h:
        return []
    rows = h["clashes"]
    if units is not None:
        rows = [c for c in rows if {c["unit1"], c["unit2"]} == set(units)]
    lines = ["", f"All-atom clashes with hydrogens added (ifacekit places {h['n_hydrogens']} hydrogens; only this check uses them): {len(rows)}"
             + ("" if units is not None else f", clashscore {h['clashscore']} per 1000 atoms"),
             "  clash = overlap of 0.4 A or more; hydrogens of ligands, waters and nucleic acids are not placed"]
    for c in rows[:max_rows]:
        lines.append(f"  {c['atom1']:<22} {c['atom2']:<22} d={c['distance']:<5} overlap={c['overlap']}")
    if len(rows) > max_rows:
        lines.append(f"  ... {len(rows) - max_rows} more (see JSON)")
    return lines


def report_quality(res, max_rows=25):
    u = res["units"][0] if res["units"] else {}
    chain = res["input"]["chains"][0]
    s = res["summary"]
    lines = _head(res, f"QUALITY of chain {chain}: {u.get('polymer_type') or 'no polymer'}, "
                       f"{u.get('n_residues', 0)} residues, {u.get('n_atoms', 0)} atoms"
                       + (f", ligands {', '.join(x.split(':')[-1] for x in u.get('ligands', []))}" if u.get("ligands") else ""))
    lines.append(f"  waters included: {res['flags'].get('n_water_residues', 0)}; ligands and glycans are "
                 + ("assigned by bond or contact, not by chain label" if res["input"].get("attach") == "contact"
                    else "taken by chain label"))
    lines.append("")
    w, allv = s["verdicts_within_units"], s["verdicts"]
    n_in, n_all = sum(w.values()), s["pairs_evaluated"]
    lines.append(f"atom pairs evaluated {n_all} (inside the chain {n_in}, with waters {n_all - n_in})")
    lines.append(_tally_line("  ", allv))
    lines.append(f"  unfavorable per 1000 non-water atoms: {s['unfavorable_per_1000_atoms']}")
    lines.append("")
    reasons = Counter()
    flips = 0
    seen = set()
    for x in res["issues"]:
        for it in x["items"]:
            key = tuple(sorted((f"{x['residue']}:{it['atom']}", it["with"])))
            if key in seen:
                continue
            seen.add(key)
            reasons[(it["verdict"], it["reason"])] += 1
            flips += it["possible_flip"]
    if reasons:
        lines.append("Problem contacts by reason (each pair counted once)")
        for (verdict, reason), c in sorted(reasons.items(), key=lambda kv: (kv[0][0] != "unfavorable", -kv[1])):
            lines.append(f"  {verdict:<12}{reason:<18}{c}")
        if flips:
            lines.append(f"  {flips} polar_mismatch contacts involve Asn/Gln/His side-chain atoms: check for a flipped amide or ring")
        lines.append("")
    bad = [x for x in res["issues"] if x["worst"] == "unfavorable"]
    lines.append(f"Residues with an unfavorable contact: {len(bad)}")
    rows = []
    for x in bad[:max_rows]:
        it = next(t for t in x["items"] if t["verdict"] == "unfavorable")
        rows.append([x["residue"], it["atom"] + (f" (alt {it['altloc']})" if it.get("altloc") else ""), it["reason"], it["distance"], it["with"],
                     "flip?" if it["possible_flip"] else ""])
    if rows:
        lines += ["  " + r for r in _table(rows, ["residue", "atom", "reason", "d", "with", ""])]
    if len(bad) > max_rows:
        lines.append(f"  ... {len(bad) - max_rows} more (see JSON)")
    nb = sum(1 for x in res["issues"] if x["worst"] == "borderline")
    lines.append(f"Residues with only borderline contacts: {nb} (listed in JSON)")
    lines += _hydrogen_lines(res)
    lines += _geometry_lines(res)
    lines += _rama_lines(res)
    lines += _rotamer_lines(res)
    lines += _flip_lines(res)
    lines += _water_lines(res)
    return "\n".join(lines)


# ---------------------------------------------------------------- pair
def report_pair(res, max_rows=30):
    a, b = res["input"]["chains"]
    lines = _head(res, f"INTERACTION between chain {a} and chain {b}")
    for u in res["units"]:
        lines.append(f"  chain {u['id']}: {u['polymer_type'] or 'ligand'}, {u['n_residues']} residues, {u['n_atoms']} atoms"
                     + (f", ligands {', '.join(x.split(':')[-1] for x in u['ligands'])}" if u["ligands"] else ""))
    if res["flags"].get("predicted_model") and res["interfaces"]:
        it0 = res["interfaces"][0]
        lines.append(f"  predicted model: the B-factor column is read as pLDDT. Interface residues, mean pLDDT "
                     f"{a} {it0['mean_b_a']}, {b} {it0['mean_b_b']}; residues below 70: {a} {it0['n_low_b_a']}, {b} {it0['n_low_b_b']}")
    sh = res["flags"].get("ligands_touching_both_chains") or []
    if sh:
        lines.append("  ligands touching both chains (counted with the first one listed): "
                     + ", ".join(f"{x['residue']} -> {x['assigned_to']}" for x in sh))
    lines.append("")
    if not res["interfaces"]:
        lines.append(f"No interface: no buried surface of at least the minimum area between {a} and {b}.")
        return "\n".join(lines)
    it = res["interfaces"][0]
    ua, ub = it["unit_a"], it["unit_b"]
    lines.append(f"Interface area {it['area']:.0f} A^2  (side {ua} {it['area_a']:.0f}, side {ub} {it['area_b']:.0f}),"
                 f" polar fraction {100 * (it['polar_fraction'] or 0):.0f}%")
    lines.append(f"Interface residues: {ua} {it['n_residues_a']}, {ub} {it['n_residues_b']}")
    c = it["contacts"]
    pis = [p for p in res["pi_interactions"] if p["inter_unit"]]
    pc = Counter(p["type"] for p in pis)
    hb = it["hbonds"]
    lines.append("Hydrogen bonds: direct {} (no water), via water {}{} (waters bonded to both sides), connections in total {}{}".format(
        hb["direct"], hb["water_mediated"],
        f" [{hb['water_mediated_clear']} clear, {hb['water_mediated_doubtful']} doubtful]" if hb["water_mediated"] else "",
        hb["total_connections"], f" ({hb['total_connections_clear']} counting clear bridges only)" if hb["water_mediated_doubtful"] else ""))
    lines.append("Contacts: hydrogen bonds {}, salt bridges {}, hydrophobic {}, other vdW {}, metal {}, "
                 "pi stacking {}, T-shaped {}, cation-pi {}".format(
                     c.get("hbond", 0), c.get("salt_bridge", 0), c.get("hydrophobic", 0), c.get("vdw", 0),
                     c.get("metal", 0), pc.get("pi_stacking", 0), pc.get("pi_tshaped", 0), pc.get("cation_pi", 0)))
    lines.append(_tally_line("Verdicts:", it["verdicts"]))
    lines.append(f"Buried polar atoms with no partner: {len(it['buried_unsatisfied_polar'])}"
                 f"   waters bridging the chains: {len(it['bridging_waters'])}")
    lines.append("")

    pair_contacts = [c for c in res["contacts"] if {c["unit1"], c["unit2"]} == {ua, ub}]

    def ordered(c):
        # always print chain a first
        if c["unit1"] == ua:
            return c["atom1"], c["atom2"]
        return c["atom2"], c["atom1"]

    polar = sorted((c for c in pair_contacts if c["type"] in ("hbond", "salt_bridge", "metal")),
                   key=lambda c: c["distance"])
    lines.append(f"Hydrogen bonds, salt bridges and metal contacts ({len(polar)}), shortest first")
    rows = []
    for c in polar[:max_rows]:
        x, y = ordered(c)
        rows.append([x, y, c["distance"], c["type"], c["verdict"]])
    if rows:
        lines += ["  " + r for r in _table(rows, [f"chain {ua}", f"chain {ub}", "d", "type", "verdict"])]
    if len(polar) > max_rows:
        lines.append(f"  ... {len(polar) - max_rows} more (see JSON)")
    lines.append("")

    if hb["mediated_waters"]:
        lines.append(f"Hydrogen bonds via water: bridging waters only, i.e. bonded to both chains ({hb['water_mediated']}: {hb['water_mediated_clear']} clear, "
                     f"{hb['water_mediated_doubtful']} doubtful): water, atoms bonded on chain {ua}, atoms bonded on chain {ub}, check")
        from .hbonds import FLAG_TEXT
        rows = [[x["water"], " ".join(f"{y['atom'].split(':', 2)[-1]} {y['distance']}" for y in x["side_a"]),
                 " ".join(f"{y['atom'].split(':', 2)[-1]} {y['distance']}" for y in x["side_b"]),
                 "clear" if not x["flags"] else "doubtful: " + ", ".join(x["flags"])] for x in hb["mediated_waters"][:max_rows]]
        lines += ["  " + r for r in _table(rows, ["water", f"on {ua}", f"on {ub}", "check"])]
        if hb["water_mediated_doubtful"]:
            lines.append("  doubtful = " + "; ".join(f"{k}: {v}" for k, v in FLAG_TEXT.items()) + " (simple limits, not calibrated)")
        if len(hb["mediated_waters"]) > max_rows:
            lines.append(f"  ... {len(hb['mediated_waters']) - max_rows} more (see JSON)")
        lines.append("")
    elif not res["flags"].get("n_water_residues"):
        lines.append("Hydrogen bonds via water: none can be counted, the model has no waters")
        lines.append("")
    prob = sorted((c for c in pair_contacts if c["verdict"] in ("borderline", "unfavorable")),
                  key=lambda c: (c["verdict"] != "unfavorable", c["distance"]))
    lines.append(f"Contacts that look wrong ({len(prob)}): unfavorable first")
    rows = []
    for c in prob[:max_rows]:
        x, y = ordered(c)
        rows.append([x, y, c["distance"], c["reason"], c["verdict"]])
    if rows:
        lines += ["  " + r for r in _table(rows, [f"chain {ua}", f"chain {ub}", "d", "reason", "verdict"])]
    elif not prob:
        lines.append("  none")
    lines.append("")

    if it["buried_unsatisfied_polar"]:
        lines.append("Buried polar atoms without a hydrogen-bond, salt-bridge or metal partner")
        for x in it["buried_unsatisfied_polar"][:max_rows]:
            lines.append(f"  {x['atom']:<18} buried {x['buried']} A^2" + ("  (role guessed)" if x["role_uncertain"] else ""))
        lines.append("")
    if it["bridging_waters"]:
        lines.append("Waters in contact with both chains: " + ", ".join(it["bridging_waters"][:20]))
        lines.append("")

    for side, key in ((ua, "residues_a"), (ub, "residues_b")):
        top = sorted(it[key], key=lambda r: -r["bsa"])[:12]
        lines.append(f"Largest buried residues, chain {side}")
        rows = [[r["residue"], r["bsa"], r["contacts"], r["problem_contacts"]] for r in top]
        lines += ["  " + r for r in _table(rows, ["residue", "buried A^2", "contacts", "problems"])]
        lines.append("")
    lines += _hydrogen_lines(res, units=(a, b))
    lines += _water_lines(res, max_rows=8)
    return "\n".join(lines).rstrip()


# ---------------------------------------------------------------- whole file
def report_all(res, max_issues=15):
    s = res["summary"]
    lines = _head(res, f"ALL units: atoms {s['n_atoms']}  units {s['n_units']}  interfaces {s['n_interfaces']}")
    lines.append("")
    v, vb = s["verdicts"], s["verdicts_between_units"]
    lines.append(f"atom pairs evaluated {s['pairs_evaluated']}")
    lines.append(_tally_line("  all:", v))
    lines.append(_tally_line("  between units / water:", vb))
    lines.append(f"  unfavorable per 1000 non-water atoms: {s['unfavorable_per_1000_atoms']}")
    lines.append("")
    if res["interfaces"]:
        lines.append("Interfaces (area = mean buried area of the two sides, A^2)")
        rows = []
        for it in res["interfaces"]:
            cc, vd = it["contacts"], it["verdicts"]
            pf = "-" if it["polar_fraction"] is None else f"{100 * it['polar_fraction']:.0f}"
            rows.append([it["unit_a"], it["unit_b"], f"{it['area']:.0f}", pf, cc.get("hbond", 0),
                         cc.get("salt_bridge", 0), cc.get("hydrophobic", 0), vd.get("borderline", 0),
                         vd.get("unfavorable", 0), len(it["buried_unsatisfied_polar"]), len(it["bridging_waters"])])
        lines += ["  " + r for r in _table(rows, ["A", "B", "area", "polar%", "HB", "salt", "hphob",
                                                 "bord", "unfav", "BUNS", "wat"])]
    else:
        lines.append("No interfaces found.")
    lines.append("")
    bad = [x for x in res["issues"] if x["worst"] == "unfavorable"]
    lines.append(f"Residues with unfavorable contacts: {len(bad)}")
    for x in bad[:max_issues]:
        it = next(t for t in x["items"] if t["verdict"] == "unfavorable")
        lines.append(f"  {x['residue']:<14} {it['atom']:<5} {it['reason']:<16} d={it['distance']}  with {it['with']}")
    if len(bad) > max_issues:
        lines.append(f"  ... {len(bad) - max_issues} more (see JSON)")
    lines += _hydrogen_lines(res)
    lines += _geometry_lines(res)
    lines += _water_lines(res)
    return "\n".join(lines)


def text_report(res):
    mode = res["input"].get("mode", "all")
    return {"quality": report_quality, "pair": report_pair}.get(mode, report_all)(res)


def quality_overview(runs):
    """One row per chain, for a multi-chain quality run."""
    head = ["chain", "type", "residues", "atoms", "pairs", "unfav", "bord", "unfav/1000", "flip?", "waters", "wat chk", "metal chk"]
    rows = []
    for r in runs:
        u = r["units"][0] if r["units"] else {}
        s = r["summary"]
        v = s["verdicts"]
        flips = len({tuple(sorted((f"{x['residue']}:{it['atom']}", it["with"])))
                     for x in r["issues"] for it in x["items"] if it["possible_flip"]})
        rows.append([r["input"]["chains"][0], u.get("polymer_type") or "-", u.get("n_residues", 0),
                     u.get("n_atoms", 0), s["pairs_evaluated"], v.get("unfavorable", 0),
                     v.get("borderline", 0), s["unfavorable_per_1000_atoms"], flips,
                     (r.get("waters") or {}).get("n", 0), (r.get("waters") or {}).get("n_listed_total", 0),
                     sum(1 for x in (r.get("waters") or {}).get("metals", []) if x["verdict"] == "borderline")])
    lines = ["Overview: each chain checked on its own (contacts between chains are not included here,",
             "use --pair for those). unfav/bord count atom pairs; flip? counts Asn/Gln/His polar mismatches.", ""]
    lines += ["  " + x for x in _table(rows, head)]
    return "\n".join(lines)


def chain_list_text(path, rows):
    lines = [f"{path}", "", "Chains in this file:"]
    table = []
    for r in rows:
        table.append([r["chain"], r["polymer_type"] or "-", r["n_polymer_residues"], r["n_atoms"],
                      ",".join(r["ligands"]) or "-", r["n_waters"]])
    lines += ["  " + x for x in _table(table, ["chain", "type", "residues", "atoms", "ligands", "waters"])]
    if any(not r["has_polymer"] for r in rows):
        lines += ["", "Chains without polymer hold only ligands or waters; they are attached to the polymer "
                      "chains they touch (--attach contact)."]
    lines += ["",
              "Choose what to inspect:",
              "  ifacekit FILE --chains A [B ...]   quality of each chain on its own",
              "  ifacekit FILE --chains all         quality of every chain, with a comparison table",
              "  ifacekit FILE --pair A B           interactions between two chains",
              "  ifacekit FILE --all                every chain and ligand in the file"]
    return "\n".join(lines)
