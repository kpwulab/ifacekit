"""HTML for the rotamer analysis: summary, plots per residue type, table by type, outliers, and the Asn/Gln/His flip check.
Used as a block in the quality and pair reports and as a stand-alone page (--rotamer-html)."""
import base64
import html

import numpy as np

from . import flips, rotamer
from .banner import banner_html
from .ramaplot import _COL, _png

FRAG_CSS = ("<style>.rgrid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:12px;margin:8px 0}"
            ".rgrid2 .rc{border:1px solid var(--line);border-radius:10px;background:var(--card);padding:8px 10px}"
            ".rgrid2 .rc b{font-size:14px}.rgrid2 .rc p{margin:0;color:var(--mut);font-size:12.5px}</style>")


def _footer():
    from .about import footer_html
    return footer_html()


def _page_css():
    from .htmlreport import CSS
    return CSS


def _thr(t):
    rotamer.grids()
    return rotamer._thr[t]


def _background(t):
    g = rotamer.grids()[t]
    proj = g.max(axis=tuple(range(2, g.ndim))) if g.ndim > 2 else g       # highest density over the remaining chi angles
    fav, alw = _thr(t)
    img = np.zeros(proj.shape[::-1] + (4,), dtype=np.uint8)               # rows = chi2 (top = largest), cols = chi1
    for mask, c in ((proj >= alw, (43, 108, 176, 38)), (proj >= fav, (43, 108, 176, 105))):
        img[mask.T[::-1]] = c
    return "data:image/png;base64," + base64.b64encode(_png(img)).decode()


def _panel(t, rows, size=250):
    pad = 34
    tot = size + 2 * pad
    lay = rotamer.LAYOUT[t]
    sx = lambda a: pad + a / 360.0 * size
    out = [f'<svg viewBox="0 0 {tot} {tot}" role="img" aria-label="{rotamer.LABEL[t]} rotamer plot">']
    if len(lay) >= 2:
        r2 = lay[1][0]
        sy = lambda a: pad + (1 - a / r2) * size
        out.append(f'<image href="{_background(t)}" x="{pad}" y="{pad}" width="{size}" height="{size}" preserveAspectRatio="none" style="image-rendering:pixelated"/>')
        yt = (0, 90, 180, 270, 360) if r2 == 360 else (0, 90, 180)
        ylab = "chi2"
        yval = lambda r: r["chi"][1]
    else:
        fav, alw = _thr(t)
        g = rotamer.grids()[t]
        ymax = 1.0
        sy = lambda v: pad + (1 - min(v, ymax)) * size
        pts = " ".join(f"{sx((i + 0.5) * 360.0 / len(g)):.1f},{sy(float(v)):.1f}" for i, v in enumerate(g))
        out.append(f'<polyline points="{pts}" fill="none" stroke="#2b6cb0" stroke-width="1.4"/>')
        for lv, lab in ((fav, "favored"), (alw, "allowed")):
            out.append(f'<line x1="{pad}" y1="{sy(lv):.1f}" x2="{pad + size}" y2="{sy(lv):.1f}" stroke="var(--mut)" stroke-dasharray="3 3" stroke-width=".7"/>'
                       f'<text x="{pad + size - 2}" y="{sy(lv) - 2:.1f}" text-anchor="end" font-size="9" fill="var(--mut)">{lab}</text>')
        yt = (0, 0.5, 1)
        ylab = "density"
        yval = lambda r: r["density"]
    out.append(f'<rect x="{pad}" y="{pad}" width="{size}" height="{size}" fill="none" stroke="var(--mut)"/>')
    for v in (0, 90, 180, 270, 360):
        out.append(f'<line x1="{sx(v)}" y1="{pad}" x2="{sx(v)}" y2="{pad + size}" stroke="var(--line)" stroke-width=".5"/>'
                   f'<text x="{sx(v)}" y="{pad + size + 13}" text-anchor="middle" font-size="10" fill="var(--mut)">{v}</text>')
    for v in yt:
        out.append(f'<line x1="{pad}" y1="{sy(v):.1f}" x2="{pad + size}" y2="{sy(v):.1f}" stroke="var(--line)" stroke-width=".5"/>'
                   f'<text x="{pad - 5}" y="{sy(v) + 3:.1f}" text-anchor="end" font-size="10" fill="var(--mut)">{v:g}</text>')
    out.append(f'<text x="{pad + size / 2}" y="{tot - 3}" text-anchor="middle" font-size="11" fill="var(--fg)">chi1</text>'
               f'<text x="9" y="{pad + size / 2}" font-size="11" fill="var(--fg)" transform="rotate(-90 9 {pad + size / 2})" text-anchor="middle">{ylab}</text>')
    for cat in ("favored", "allowed", "outlier"):
        for r in rows:
            if r["category"] == cat:
                out.append(f'<circle cx="{sx(r["chi"][0]):.1f}" cy="{sy(yval(r)):.1f}" r="{(2.8 if "model" in r else 3.4) if cat == "outlier" else (1.6 if "model" in r else 2.0)}" fill="{_COL[cat]}" '
                           f'fill-opacity="{(.9 if "model" in r else 1) if cat == "outlier" else (.45 if "model" in r else .7)}"><title>{html.escape(rotamer.label(r))} chi {" ".join(str(x) for x in r["chi"])} ({cat})</title></circle>')
    out.append("</svg>")
    return "".join(out)


def _plots(rows):
    cards = []
    for t in rotamer.TABLES:
        sub = [r for r in rows if r["table"] == t]
        if not sub:
            continue
        k = {x: sum(1 for r in sub if r["category"] == x) for x in ("favored", "allowed", "outlier")}
        cards.append(f'<div class="rc"><b>{rotamer.LABEL[t]}</b><p>{len(sub)} residues: {k["favored"]} favored, {k["allowed"]} allowed, {k["outlier"]} outlier</p>{_panel(t, sub)}</div>')
    return ("<p class='note'>Plots by residue type: chi1 against chi2 (types with two or more chi angles) or the reference density along chi1 "
            "(one chi angle). Shading is the highest reference density over the remaining chi angles: favored dark, allowed light "
            "(for three and four chi angles the true region in the full space is smaller). Point colour is the category in the full space. Hover for the residue.</p>"
            f"<div class='rgrid2'>{''.join(cards)}</div>")


def flips_block(path_in, chains, model_index=0):
    summ, rows = flips.analyze(path_in, model_index, chains)
    if not summ["n"]:
        return ""
    from .htmlreport import table
    c = summ["counts"]
    lst = sorted([r for r in rows if r["verdict"] != "keep"], key=lambda r: -r["delta"])
    tbl = table(["residue", "score as modelled", "score flipped", "gain", "verdict"],
                [[r["residue"], r["score_modelled"], r["score_flipped"], r["delta"], r["verdict"]] for r in lst[:20]]) if lst else ""
    return (f"<h3>Asn, Gln and His flips (chains {html.escape(', '.join(chains))})</h3>"
            f"<p class='note'>{summ['n']} groups: flip suggested for <b>{c['flip']}</b>, uncertain for <b>{c['uncertain']}</b>, keep for {c['keep']}. "
            "The two end atoms of an Asn or Gln amide, and the two sides of the His ring, cannot be told apart in the map. Each group is scored with the "
            "hydrogen bonds and overlaps of the model as built and with the group turned by 180 degrees (hydrogens placed from the geometry). "
            "A suggestion is not a decision: against Reduce on 2024 and later entries under a third of the suggested flips were also made by Reduce, "
            f"and about a fifth of Reduce's flips were found (see the design document). Check the density. {summ['method']}.</p>" + tbl
            + (f"<p class='note'>{len(lst) - 20} more in the JSON</p>" if len(lst) > 20 else ""))


def fragment(path_in, chains, model_index=0, anchor="rotamer", models=None, with_flips=True):
    summ, rows = rotamer.analyze(path_in, model_index, chains, models)
    if not summ["n"]:
        return ""
    from .htmlreport import table
    c, p = summ["counts"], summ["percent"]
    by = [[rotamer.LABEL[t], v["n"], v["favored"], v["allowed"], v["outlier"]] for t, v in summ["by_type"].items()]
    outl = [[rotamer.label(o), rotamer.LABEL[o["table"]], " ".join(str(x) for x in o["chi"]), o["density"],
             (o.get("category_flipped") if "category_flipped" in o else "-")] for o in summ["outliers"][:20]]
    af = summ.get("after_suggested_flips") or {}
    aft = (f" After the {af['n_flips_applied']} flip(s) suggested below: {af['percent']['favored']}% favored, "
           f"{af['counts']['outlier']} outliers ({af['percent']['outlier']}%)." if af.get("n_flips_applied") else "")
    return (FRAG_CSS + f"<h3 id='{anchor}'>Rotamers (chains {html.escape(', '.join(chains))})</h3>"
            f"<p class='note'>{summ['n']} side chains with chi angles: <b>{c['favored']}</b> favored ({p['favored']}%), "
            f"<b>{c['allowed']}</b> allowed ({p['allowed']}%), <b>{c['outlier']}</b> outliers ({p['outlier']}%). Goal: outliers under 0.3%.{aft} "
            "Scored against the ifacekit reference built from high-resolution PDB entries (density at the chi angles, compared with the "
            "levels of that residue type). No hydrogens are used. Asn, Gln and His are also scored with the group flipped (last column of the outlier table). "
            "Leu has a higher outlier rate than other types by the MolProbity convention that the cut-offs follow.</p>"
            + _plots(rows)
            + table(["type", "residues", "favored", "allowed", "outlier"], by, cls={4: lambda v: "unfavorable" if v else ""})
            + ("<p class='note'>Rotamer outliers (density, chi angles in degrees, category if flipped):</p>"
               + table(["residue", "type", "chi angles", "density (relative to peak)", "if flipped"], outl) if outl else "")
            + (f"<p class='note'>{len(summ['outliers']) - 20} more in the JSON</p>" if len(summ["outliers"]) > 20 else "")
            + (f"<p class='note'><b>All {len(models)} models of the ensemble in one plot</b> (hover for the model); the numbers count side chains in models.</p>" if models is not None else "")
            + (flips_block(path_in, chains, model_index) if with_flips else ""))


def build(path_in, out_path, chains=None, model_index=0, title=None, models=None):
    chains = chains or []
    body = fragment(path_in, chains, model_index, models=models) if chains else ""
    doc = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
           f'<title>Rotamer analysis</title><style>{_page_css()}</style></head><body><main>'
           + banner_html() + f'<h1>Rotamer analysis{": " + html.escape(title) if title else ""}</h1>{body}'
           '<p style="color:var(--mut);font-size:12.5px;margin-top:14px">Calculated by ifacekit with its own reference, no external program.</p>' + _footer() + '</main></body></html>')
    with open(out_path, "w") as fh:
        fh.write(doc)
