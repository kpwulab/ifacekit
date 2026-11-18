"""Self-contained HTML page with Ramachandran plots (one panel per residue class)."""
import base64
import html
import struct
import zlib

import numpy as np

from . import rama
from .banner import banner_html

_COL = {"favored": "#2f855a", "allowed": "#b7791f", "outlier": "#c0392b"}


def _png(rgba):
    h, w, _ = rgba.shape
    raw = b"".join(b"\x00" + rgba[y].tobytes() for y in range(h))

    def chunk(t, d):
        c = struct.pack(">I", len(d)) + t + d
        return c + struct.pack(">I", zlib.crc32(t + d) & 0xFFFFFFFF)
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


def _background(cls):
    g = rama.grids()[cls]
    img = np.zeros((180, 180, 4), dtype=np.uint8)           # rows = psi (top is +180), cols = phi
    fav, alw = g >= rama.FAV_T[cls], g >= rama.ALW_T[cls]
    rgba_f, rgba_a = (43, 108, 176, 105), (43, 108, 176, 38)
    for (mask, c) in ((alw, rgba_a), (fav, rgba_f)):
        m = mask.T[::-1]
        img[m] = c
    return "data:image/png;base64," + base64.b64encode(_png(img)).decode()


def _panel(cls, rows, size=300):
    pad = 34
    tot = size + 2 * pad
    sx = lambda phi: pad + (phi + 180) / 360 * size
    sy = lambda psi: pad + (180 - psi) / 360 * size
    out = [f'<svg viewBox="0 0 {tot} {tot}" role="img" aria-label="{rama.LABEL[cls]} Ramachandran plot">',
           f'<image href="{_background(cls)}" x="{pad}" y="{pad}" width="{size}" height="{size}" '
           'preserveAspectRatio="none" style="image-rendering:pixelated"/>',
           f'<rect x="{pad}" y="{pad}" width="{size}" height="{size}" fill="none" stroke="var(--mut)"/>']
    for v in (-180, -90, 0, 90, 180):
        out.append(f'<line x1="{sx(v)}" y1="{pad}" x2="{sx(v)}" y2="{pad + size}" stroke="var(--line)" stroke-width=".5"/>')
        out.append(f'<line x1="{pad}" y1="{sy(v)}" x2="{pad + size}" y2="{sy(v)}" stroke="var(--line)" stroke-width=".5"/>')
        out.append(f'<text x="{sx(v)}" y="{pad + size + 14}" text-anchor="middle" font-size="10" fill="var(--mut)">{v}</text>')
        out.append(f'<text x="{pad - 6}" y="{sy(v) + 3}" text-anchor="end" font-size="10" fill="var(--mut)">{v}</text>')
    out.append(f'<text x="{pad + size / 2}" y="{tot - 4}" text-anchor="middle" font-size="11" fill="var(--fg)">phi</text>')
    out.append(f'<text x="10" y="{pad + size / 2}" font-size="11" fill="var(--fg)" transform="rotate(-90 10 {pad + size / 2})" text-anchor="middle">psi</text>')
    for cat in ("favored", "allowed", "outlier"):           # outliers drawn last, on top
        for r in rows:
            if r["category"] == cat:
                ens = "model" in r                    # several models in one plot: smaller, see-through points
                rad = (3.0 if ens else 3.6) if cat == "outlier" else (1.7 if ens else 2.2)
                out.append(f'<circle cx="{sx(r["phi"]):.1f}" cy="{sy(r["psi"]):.1f}" r="{rad}" fill="{_COL[cat]}" '
                           f'fill-opacity="{(.9 if ens else 1) if cat == "outlier" else (.45 if ens else .75)}"><title>{html.escape(r["residue"])}'
                           f'{" model " + str(r["model"]) if ens else ""} phi {r["phi"]} psi {r["psi"]} ({cat})</title></circle>')
    out.append("</svg>")
    return "".join(out)


CSS = """:root{--bg:#fbfbfa;--fg:#1d2326;--mut:#667078;--card:#fff;--line:#dde1e3}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#14181a;--fg:#e6e9ea;--mut:#98a2a8;--card:#1c2124;--line:#2e3539}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif}
main{max-width:1060px;margin:0 auto;padding:16px 16px 48px}h1{font-size:22px;margin:.2em 0}h2{font-size:15px;margin:0 0 4px}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:14px;margin-top:14px}
.card{border:1px solid var(--line);border-radius:10px;background:var(--card);padding:10px 12px}.card p{margin:0;color:var(--mut);font-size:13px}
table{border-collapse:collapse;font-size:13px}td,th{padding:3px 12px 3px 0;text-align:left}th{color:var(--mut);font-weight:600}
.lg span{display:inline-block;width:10px;height:10px;border-radius:50%;margin:0 4px 0 12px}"""


def _footer():
    from .about import footer_html
    return footer_html()


FRAG_CSS = ("<style>.rgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:12px;margin:8px 0}"
            ".rgrid .rc{border:1px solid var(--line);border-radius:10px;background:var(--card);padding:8px 10px}"
            ".rgrid .rc b{font-size:14px}.rgrid .rc p{margin:0;color:var(--mut);font-size:12.5px}"
            "p.pdfdl{margin:6px 0 10px}a.pdfdl,a.pdfdl:visited{display:inline-block;font-size:16px;font-weight:600;color:var(--acc);text-decoration:none;"
            "padding:7px 14px;border:1.5px solid var(--acc);border-radius:8px;background:var(--card)}a.pdfdl:hover{background:var(--chip)}</style>")


def fragment(path_in, chains, model_index=0, anchor="rama", pdf_name=None, models=None):
    """HTML block (heading, summary, plots, outlier table) to embed in the pair report. models: indices of all models of an ensemble to plot together."""
    summ, rows = rama.analyze(path_in, model_index, chains, models)
    if not summ["n"]:
        return ""
    c, p = summ["counts"], summ["percent"]
    cards = []
    for cls in rama.CLASSES:
        sub = [r for r in rows if r["class"] == cls]
        if sub:
            k = {x: sum(1 for r in sub if r["category"] == x) for x in ("favored", "allowed", "outlier")}
            cards.append(f'<div class="rc"><b>{rama.LABEL[cls]}</b><p>{len(sub)} residues: {k["favored"]} favored, '
                         f'{k["allowed"]} allowed, {k["outlier"]} outlier</p>{_panel(cls, sub, 240)}</div>')
    mc = models is not None
    outl = "".join(f"<tr><td>{html.escape(o['residue'])}</td>{'<td>' + str(o['model']) + '</td>' if mc else ''}<td>{rama.LABEL[o['class']]}</td><td>{o['phi']}</td><td>{o['psi']}</td></tr>"
                   for o in summ["outliers"][:20])
    tbl = (f"<table><tr><th>residue</th>{'<th>model</th>' if mc else ''}<th>class</th><th>phi</th><th>psi</th></tr>{outl}</table>" if outl else "")
    dl = (f"<p class='pdfdl'><a class='pdfdl' href='{html.escape(pdf_name)}' download>Download all Ramachandran and rotamer plots as one PDF</a></p>"
          if pdf_name else "")
    return (FRAG_CSS + f"<h3 id='{anchor}'>Ramachandran (chains {html.escape(', '.join(chains))})</h3>" + dl
            + (f"<p class='note'><b>All {len(models)} models of the ensemble in one plot</b>: every point is one residue in one model (hover for the model); the numbers count residue-model pairs.</p>" if mc else "") +
            f"<p class='note'>{summ['n']} residues: <b>{c['favored']}</b> favored ({p['favored']}%), <b>{c['allowed']}</b> allowed "
            f"({p['allowed']}%), <b>{c['outlier']}</b> outliers ({p['outlier']}%). Goals: favored &gt; 98%, outliers &lt; 0.2%. "
            f"Shaded areas are the favored (dark) and allowed (light) regions of the ifacekit reference (98% and 99.95% contours); hover a point for the residue.</p>"
            f"<div class='rgrid'>{''.join(cards)}</div>{tbl}")


def build(path_in, out_path, chains=None, model_index=0, title=None, models=None):
    summ, rows = rama.analyze(path_in, model_index, chains, models)
    mc = models is not None
    c, p = summ["counts"], summ["percent"]
    cards = []
    for cls in rama.CLASSES:
        sub = [r for r in rows if r["class"] == cls]
        if not sub:
            continue
        k = {x: sum(1 for r in sub if r["category"] == x) for x in ("favored", "allowed", "outlier")}
        cards.append(f'<div class="card"><h2>{rama.LABEL[cls]}</h2><p>{len(sub)} residues: {k["favored"]} favored, '
                     f'{k["allowed"]} allowed, {k["outlier"]} outlier</p>{_panel(cls, sub)}</div>')
    outl = "".join(f"<tr><td>{html.escape(o['residue'])}</td>{'<td>' + str(o['model']) + '</td>' if mc else ''}<td>{rama.LABEL[o['class']]}</td><td>{o['phi']}</td><td>{o['psi']}</td></tr>"
                   for o in summ["outliers"]) or '<tr><td colspan="5">none</td></tr>'
    doc = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
           f'<title>Ramachandran plot</title><style>{CSS}</style></head><body><main>'
           + banner_html() + f'<h1>Ramachandran analysis{": " + html.escape(title) if title else ""}</h1>'
           f'<p>{summ["n"]} residues with phi and psi{" in chain(s) " + ", ".join(chains) if chains else ""}{" in " + str(len(models)) + " models (every point is one residue in one model; hover for the model)" if mc else ""}. '
           f'<b>{c["favored"]}</b> favored ({p["favored"]}%), <b>{c["allowed"]}</b> allowed ({p["allowed"]}%), '
           f'<b>{c["outlier"]}</b> outliers ({p["outlier"]}%). Goals: favored &gt; 98%, outliers &lt; 0.2%.</p>'
           f'<p class="lg" style="color:var(--mut);font-size:13px">Shaded: favored (dark) and allowed (light) regions of the ifacekit reference (98% and 99.95% contours). '
           f'<span style="background:{_COL["favored"]}"></span>favored <span style="background:{_COL["allowed"]}"></span>allowed '
           f'<span style="background:{_COL["outlier"]}"></span>outlier. Hover a point for the residue.</p>'
           f'<div class="grid">{"".join(cards)}</div>'
           f'<div class="card" style="margin-top:14px"><h2>Outliers</h2><table><tr><th>residue</th>{"<th>model</th>" if mc else ""}<th>class</th><th>phi</th><th>psi</th></tr>{outl}</table></div>'
           f'<p style="color:var(--mut);font-size:12.5px;margin-top:14px">Reference: {summ["reference"]}. Calculated by ifacekit, no external program.</p>'
           + _footer() + '</main></body></html>')
    with open(out_path, "w") as fh:
        fh.write(doc)
    return summ
