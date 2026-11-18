"""Self-contained HTML report for the quality of single chains (same look as the pair and ligand reports).

One overview (comparison table, B-factor or pLDDT plot, attention list) and one section per chain with contact
problems, backbone geometry, Ramachandran plots, waters and metal sites.
"""
import html

from . import af3, rama
from .about import comparison_html
from .banner import banner_html
from .attach import plan_attachments
from .core import analyze
from . import htmlreport as _hr
from .htmlreport import (BODY_END, BODY_START, CSS, NAV_JS, card, chem_table, clash_table, esc, hclash_html, plddt_plot, split_problems, st_hi, st_lo,
                         table)
from .ramaplot import fragment as rama_fragment
from .rotamerhtml import fragment as rot_fragment


def collect(path, chains, model_index=0, predicted=None, pdf_name=None, hydrogens=False, ensemble=False):
    kw = dict(model_index=model_index, n_points=480, record_intra="problems", attach="contact", hydrogens=hydrogens)
    if len(chains) > 1:
        kw["plan"] = plan_attachments(path, model_index)
    runs = [analyze(path, chain=c, **kw) for c in chains]
    pred = bool(runs[0]["flags"].get("predicted_model")) if predicted is None else predicted
    if ensemble:
        from .ensemble import residue_rmsd
        pred = False
    info = []
    for c, r in zip(chains, runs):
        rm = r.get("ramachandran") or {}
        clashes, chem = split_problems(r["contacts"])
        info.append({"chain": c, "run": r, "rama": rm, "clashes": clashes, "chem": chem,
                     "unit": r["units"][0] if r["units"] else {}})
    return {"pdf": pdf_name, "path": path, "chains": chains, "model_index": model_index, "info": info, "predicted": pred,
            "ensemble": bool(ensemble), "kind": "predicted" if pred else "experimental", "iface": None, "n_models": runs[0]["flags"].get("n_models_in_file", 1), "model_index": model_index,
            "residues": [r for r in residue_rmsd(path) if r["chain"] in chains] if ensemble else af3.residue_plddt(path, model_index)}


def _status(x):
    return {"bad": "bad", "warn": "warn", "good": "good"}.get(x, "")


def _chain_state(i):
    r = i["run"]
    v = r["summary"]["verdicts"]
    rm = i["rama"]
    g = r.get("geometry") or {}
    rot = i["run"].get("rotamer") or {}
    out = {"rot_n": rot.get("n", 0), "rot_out_pct": (rot.get("percent") or {}).get("outlier"), "rot_out": (rot.get("counts") or {}).get("outlier", 0),
           "unfav": v.get("unfavorable", 0), "bord": v.get("borderline", 0),
           "fav_pct": (rm.get("percent") or {}).get("favored"), "out_pct": (rm.get("percent") or {}).get("outlier"),
           "out_n": (rm.get("counts") or {}).get("outlier", 0), "geo": g.get("n_outliers", 0)}
    w = r.get("waters") or {}
    out["wat"] = w.get("n_listed_total", len(w.get("listed", [])))
    out["met"] = sum(1 for x in w.get("metals", []) if x["verdict"] == "borderline")
    return out


def _water_html(r, limit=15):
    w = r.get("waters") or {}
    if not w or (not w.get("n") and not w.get("metals")):
        return "<p class='note'>no waters or metal ions attached to this chain</p>"
    out = []
    if w["n"]:
        c, fc = w["counts"], w["flag_counts"]
        out.append("<div class='cards'>" + "".join(card(k.replace("_", " "), c[k], "", "bad" if k == "unsupported" else "warn" if k == "weak" else "")
                                                  for k in ("metal_bound", "well_supported", "supported", "weak", "unsupported") if c.get(k)) + "</div>")
        parts = [f"{k} {fc[k]}" for k in ("too_close", "ww_close", "close_carbon", "high_b", "low_occupancy") if fc.get(k)]
        out.append(f"<p class='note'>{w['n']} waters. Flags: {esc(', '.join(parts) or 'none')}"
                   + (f"; waters touching two or more units: {fc['bridges_units']}" if fc.get("bridges_units") else "")
                   + ". Partner = N or O within 3.5 &Aring; (waters included); weak = one partner, unsupported = none.</p>")
        rows = [[x["water"], x["category"], ", ".join(x["flags"]) or "-",
                 ", ".join(f"{p['atom'].split(':', 1)[1]} {p['d']}" for p in x["partners"][:2]) or "-", x["b"]]
                for x in w["listed"][:limit]]
        if rows:
            out.append(f"<p class='note'>Waters to check ({w.get('n_listed_total', len(w['listed']))} flagged, {len(rows)} shown):</p>"
                       + table(["water", "category", "flags", "nearest partners", "B"], rows))
    ms = w.get("metals") or []
    if ms:
        bad = [x for x in ms if x["verdict"] == "borderline"]
        out.append(f"<p class='note'>Metal sites {len(ms)}: {len(ms) - len(bad)} ok or informational, {len(bad)} to check.</p>")
        rows = []
        for x in ms:
            ds = [l["d"] for l in x["ligands"]]
            rows.append([x["metal"], x["n_ligands"], f"{min(ds)}-{max(ds)}" if ds else "-", sum(l["water"] for l in x["ligands"]),
                         x["octahedral_rms_deg"] if x["octahedral_rms_deg"] is not None else "-", "; ".join(x["notes"]) or x["verdict"]])
        out.append(table(["metal", "ligands", "d range (&Aring;)", "waters", "octahedral dev (deg)", "note"], rows))
    return "".join(out)


def _geometry_html(g, limit=12):
    if not g or not g.get("n_bonds"):
        return ""
    out = [f"<h3>Backbone geometry</h3><p class='note'>bond rmsZ {g['rmsz_bonds']}, angle rmsZ {g['rmsz_angles']}, "
           f"outliers (|Z| over 4): {g['n_outliers']}. Approximate Engh and Huber values.</p>"]
    if g["outliers"]:
        out.append(table(["residue", "kind", "atoms", "value", "ideal", "Z"],
                         [[o["residue"], o["kind"], o["atoms"], o["value"], o["ideal"], o["z"]] for o in g["outliers"][:limit]]))
    st = [b for b in g["chain_breaks"] if b["kind"] == "stretched_bond"]
    gaps = [b for b in g["chain_breaks"] if b["kind"] == "gap_in_numbering"]
    if st:
        out.append("<p class='note'>Stretched peptide bonds: " + esc(", ".join(f"{b['after']}-{b['before'].split(':')[1]} ({b['C_N_distance']} A)" for b in st[:8])) + "</p>")
    if gaps:
        out.append(f"<p class='note'>Gaps in numbering (missing residues): {len(gaps)}</p>")
    return "".join(out)


def render(d, molstar=False):
    _hr._MS[0] = bool(molstar)        # table rows are clickable (data-a) only for the page with the 3D panel
    try:
        return _render(d, molstar)
    finally:
        _hr._MS[0] = False


def _render(d, molstar):
    name = d["path"].replace("\\", "/").split("/")[-1]
    infos, pred = d["info"], d["predicted"]
    r0 = infos[0]["run"]
    chains = d["chains"]
    states = {i["chain"]: _chain_state(i) for i in infos}
    bname = "pLDDT" if pred else "B-factor"
    out = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
           f"<title>{esc(name)} chain quality</title><style>{CSS}</style></head><body><main>",
           banner_html() + f"<h1>{esc(name)}</h1>",
           f"<p class='sub'>Quality of chain{'s' if len(chains) > 1 else ''} {esc(', '.join(chains))}, each checked on its own with its ligands, glycans and waters. "
           f"Report made by ifacekit {esc(r0['tool']['version'])} on {esc(r0['input']['format'])}. "
           f"Model type: {'prediction (confidence column read as pLDDT)' if pred else ('experimental ensemble of ' + str(r0['flags'].get('n_models_in_file', '')) + ' models: the tables and the 3D view switch with the model chosen above; the plot and the summary show all models' if d.get('ensemble_page') else 'experimental ensemble: the tables are for the first model, the plot shows all ' + str(r0['flags'].get('n_models_in_file', '')) + ' models') if d.get('ensemble') else 'experimental (confidence column read as B-factor)'}. "
           "Heavy-atom analysis, no hydrogens; how it compares with MolProbity is stated at the end of the report. Contacts between chains are not included: use the pair report for those.</p>"]
    out.append(BODY_START)
    nav = ["<nav class='top' aria-label='Report sections'><a href='#s1' data-s='s1'>Overview</a>"]
    for c in chains:
        nav.append(f"<a href='#c-{esc(c)}' data-s='c-{esc(c)}'>Chain {esc(c)}</a>")
    if molstar:
        nav.append("<a href='#' id='ms-open' class='active'>3D view</a>")
    nav.append("<span class='sp'></span><button class='theme' onclick=\"var r=document.documentElement;r.dataset.theme=(r.dataset.theme==='dark'?'light':'dark')\">light / dark</button></nav>")
    out.append("".join(nav))

    # chips and attention list
    chips, attn = [], []
    for c in chains:
        s = states[c]
        st_c = "bad" if s["unfav"] else "warn" if s["bord"] else "good"
        chips.append(f"<a href='#c-{esc(c)}' class='st-{st_c}'><b>{s['unfav']}</b><span>chain {esc(c)} unfavorable pairs</span></a>")
        if s["fav_pct"] is not None:
            st_r = st_hi(s["fav_pct"], 98, 95) if not s["out_n"] else "bad" if (s["out_pct"] or 0) > 0.2 else "warn"
            chips.append(f"<a href='#c-{esc(c)}-rama' class='st-{st_r}'><b>{s['fav_pct']:.1f}%</b><span>chain {esc(c)} Ramachandran favored</span></a>")
            if st_r != "good":
                attn.append((st_r, f"#c-{esc(c)}-rama", f"chain {esc(c)}: {s['fav_pct']:.1f}% favored, {s['out_n']} Ramachandran outlier(s) (goals: over 98% and under 0.2%)"))
        if s["rot_n"]:
            st_o = st_lo(s["rot_out_pct"], 0.3, 1.0)
            chips.append(f"<a href='#c-{esc(c)}-rot' class='st-{st_o}'><b>{s['rot_out_pct']:.1f}%</b><span>chain {esc(c)} rotamer outliers</span></a>")
            if st_o != "good":
                attn.append((st_o, f"#c-{esc(c)}-rot", f"chain {esc(c)}: {s['rot_out']} rotamer outlier(s), {s['rot_out_pct']:.1f}% (goal under 0.3%)"))
        if s["unfav"]:
            attn.append(("bad", f"#c-{esc(c)}", f"chain {esc(c)}: {s['unfav']} unfavorable atom pair(s)"))
        if s["geo"]:
            attn.append(("warn", f"#c-{esc(c)}", f"chain {esc(c)}: {s['geo']} backbone geometry outlier(s)"))
        if s["met"]:
            attn.append(("warn", f"#c-{esc(c)}", f"chain {esc(c)}: {s['met']} metal site(s) to check"))
    out.append("<div class='chips'>" + "".join(chips) + "</div>")
    if attn:
        worst = "bad" if any(a[0] == "bad" for a in attn) else "warn"
        out.append(f"<div class='attn {worst}'><b>Worth a look</b><ul>" + "".join(f"<li><a href='{a[1]}'>{a[2]}</a></li>" for a in attn) + "</ul></div>")
    else:
        out.append("<div class='attn'><b>Nothing flagged</b> by the checks in this report.</div>")

    # overview
    out.append('<section id="s1" class="sec s1"><h2>1. Overview</h2>')
    rows = []
    for i in infos:
        s, u = states[i["chain"]], i["unit"]
        rows.append([i["chain"], u.get("polymer_type") or "-", u.get("n_residues", 0), u.get("n_atoms", 0), s["unfav"], s["bord"],
                     i["run"]["summary"]["unfavorable_per_1000_atoms"],
                     "-" if s["fav_pct"] is None else f"{s['fav_pct']:.1f}", "-" if s["out_pct"] is None else s["out_n"], "-" if not s["rot_n"] else s["rot_out"], s["geo"], s["wat"], s["met"]])
    out.append(table(["chain", "type", "residues", "atoms", "unfavorable", "borderline", "unfav / 1000 atoms", "Rama favored %", "Rama outliers", "rotamer outliers",
                      "geometry outliers", "waters to check", "metals to check"], rows,
                     cls={4: lambda v: "unfavorable" if v else "", 5: lambda v: "borderline" if v else ""}))
    out.append("<p class='note'>Unfavorable and borderline count atom pairs inside the chain and with its waters and ligands. "
               "Unfavorable per 1000 atoms is not a MolProbity clashscore.</p>")
    iface = set()
    if d.get("ensemble"):
        out.append("<h3 id='plddt'>Backbone deviation along the sequence</h3>" + plddt_plot(d["residues"], chains, iface, exp=True, label="deviation (A)")
                   + "<p class='legend'>Line: mean over the models of the distance of each residue from the mean structure (root of the mean squared distance of the backbone atoms N, CA, C, O to their mean positions, "
                   "each chain superposed on its least mobile 70% of residues). Dark band: mean &plusmn; 1 standard deviation; light band: lowest to highest model; the RMSD is in the tooltip. "
                   "Scroll over the plot to zoom, drag to pan, and hover to read a residue.</p>")
    else:
        out.append(f"<h3 id='plddt'>{bname} along the sequence</h3>" + plddt_plot(d["residues"], chains, iface, exp=not pred)
                   + "<p class='legend'>" + ("Dashed lines at 90, 70 and 50. " if pred else "") + "Scroll over the plot to zoom, drag to pan, and hover to read a residue.</p>")
    out.append("</section>")

    # one section per chain
    for k, i in enumerate(infos):
        c, r, s = i["chain"], i["run"], states[i["chain"]]
        sm = r["summary"]
        v = sm["verdicts"]
        u = i["unit"]
        lig = ", ".join(x.split(":")[-1] for x in u.get("ligands", []))
        out.append(f'<section id="c-{esc(c)}" class="sec s2"><h2>{k + 2}. Chain {esc(c)}</h2>')
        out.append(f"<p class='note'>{esc(u.get('polymer_type') or 'no polymer')}, {u.get('n_residues', 0)} residues, {u.get('n_atoms', 0)} atoms"
                   + (f", with {esc(lig)}" if lig else "") + f", {r['flags'].get('n_water_residues', 0)} waters.</p>")
        out.append("<div class='cards'>" + "".join([
            card("unfavorable pairs", v.get("unfavorable", 0), "steric overlap, short H-bond, mismatch", "bad" if v.get("unfavorable", 0) else "good"),
            card("borderline pairs", v.get("borderline", 0), "", "warn" if v.get("borderline", 0) else "good"),
            card("favorable pairs", v.get("favorable", 0), ""),
            card("unfavorable per 1000 atoms", sm["unfavorable_per_1000_atoms"], "not a MolProbity clashscore"),
            card("Ramachandran favored", "-" if s["fav_pct"] is None else f"{s['fav_pct']:.1f}%", "goal over 98%",
                 "" if s["fav_pct"] is None else st_hi(s["fav_pct"], 98, 95)),
            card("Ramachandran outliers", s["out_n"], "goal under 0.2%", "" if s["out_pct"] is None else st_lo(s["out_pct"], 0.2, 0.5)),
            card("rotamer outliers", "-" if not s["rot_n"] else f"{s['rot_out']} ({s['rot_out_pct']:.1f}%)", "goal under 0.3%",
                 "" if not s["rot_n"] else st_lo(s["rot_out_pct"], 0.3, 1.0)),
            card("geometry outliers", s["geo"], "|Z| over 4", "bad" if s["geo"] else "good")]) + "</div>")
        out.append(f"<h3>Clashes ({len(i['clashes'])}, of which {sum(1 for x in i['clashes'] if x['verdict'] == 'unfavorable')} unfavorable)</h3>" + clash_table(i["clashes"], 20) + hclash_html(r.get("hydrogens"), limit=20))
        out.append(f"<h3>Other contact issues ({len(i['chem'])})</h3>" + chem_table(i["chem"], 20))
        out.append(_geometry_html(r.get("geometry")))
        try:
            out.append(rama_fragment(d["path"], [c], d["model_index"], anchor=f"c-{c}-rama", pdf_name=d.get("pdf")))
        except Exception as e:  # never let the extra section break the report
            out.append(f"<p class='note'>Ramachandran analysis failed: {esc(e)}</p>")
        try:
            out.append(rot_fragment(d["path"], [c], d["model_index"], anchor=f"c-{c}-rot"))
        except Exception as e:
            out.append(f"<p class='note'>Rotamer analysis failed: {esc(e)}</p>")
        out.append("<h3>Waters and metal sites</h3>" + _water_html(r))
        out.append("</section>")
    out.append("<h3>What this report does not cover</h3><p class='note'>Hydrogen-based clashes and Asn, Gln and His flips are not evaluated. "
               "Use MolProbity or the wwPDB validation report for those. " + comparison_html() + "</p>")
    out.append(NAV_JS)
    out.append(BODY_END)
    from .about import footer_html
    out.append(footer_html())
    if molstar and not d.get("ensemble_page"):
        out.append(_hr.molstar_block({**d, "residues": d["residues"]}))
    out.append("</main></body></html>")
    return "".join(out)


def build(path, chains, out_path, model_index=0, predicted=None, pdf_name=None, hydrogens=False, ensemble=False):
    d = collect(path, chains, model_index, predicted, pdf_name, hydrogens, ensemble)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(render(d))
    return d
