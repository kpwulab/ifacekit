"""One PDF with all Ramachandran and rotamer plots of the chosen chains (A4 pages).

The plots are the same SVG panels as in the HTML reports, with the page colours written out (an SVG in a PDF has no CSS
variables), laid out as 2 x 3 panels per page: a summary page, then for every chain its Ramachandran pages, its rotamer pages
and a list of outliers. Needs cairosvg (with the system library libcairo) and pypdf: pip install "ifacekit[pdf]".
"""
import html
import io
import textwrap

from . import flips, rama, rotamer
from . import ramaplot, rotamerhtml

W, H = 595, 842                       # A4 in points
COL = {"var(--fg)": "#1d2326", "var(--mut)": "#667078", "var(--line)": "#dde1e3", "var(--card)": "#ffffff"}
PANEL = 150                           # plot size; the panel including axes is PANEL + 68
CELL_W, CELL_H = PANEL + 68, PANEL + 68 + 28


def _need():
    try:
        import cairosvg
        import pypdf
    except Exception as e:             # ImportError, or OSError when libcairo is missing
        raise RuntimeError("the PDF needs cairosvg (with libcairo) and pypdf: pip install 'ifacekit[pdf]' (%s)" % e)
    return cairosvg, pypdf


def _resolve(svg):
    for k, v in COL.items():
        svg = svg.replace(k, v)
    return svg


def _panel_at(svg, x, y, size):
    return svg.replace("<svg viewBox", f'<svg x="{x}" y="{y}" width="{size}" height="{size}" viewBox', 1)


def _text(x, y, t, size=11, weight="normal", fill="#1d2326", anchor="start"):
    return (f'<text x="{x}" y="{y}" font-family="DejaVu Sans, Helvetica, Arial, sans-serif" font-size="{size}" '
            f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{html.escape(str(t))}</text>')


def _page(parts):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{W}pt" height="{H}pt" '
            f'viewBox="0 0 {W} {H}"><rect width="{W}" height="{H}" fill="#ffffff"/>' + "".join(parts) + "</svg>")


def _grid_pages(title, subtitle, panels):
    """panels: list of (caption, svg). Six per page, two columns."""
    pages = []
    for k in range(0, len(panels), 6):
        parts = [_text(30, 34, title, 15, "bold"), _text(30, 52, subtitle, 10, fill="#667078")]
        for n, (cap, svg) in enumerate(panels[k:k + 6]):
            cx, cy = 30 + (n % 2) * (CELL_W + 20), 68 + (n // 2) * (CELL_H + 6)
            parts.append(_text(cx + 6, cy + 12, cap, 11, "bold"))
            parts.append(_panel_at(_resolve(svg), cx, cy + 18, CELL_W))
        pages.append(_page(parts))
    return pages


def _wrap(lines, width=98):
    out = []
    for l in lines:
        lead = len(l) - len(l.lstrip())
        out += textwrap.wrap(l, width, subsequent_indent=" " * (lead + 6), replace_whitespace=False, drop_whitespace=False) or [""]
    return out


def _lines_pages(title, lines, per_page=52):
    lines = _wrap(lines)
    pages = []
    for k in range(0, max(len(lines), 1), per_page):
        parts = [_text(30, 34, title, 15, "bold")]
        for n, l in enumerate(lines[k:k + per_page]):
            parts.append(_text(30, 62 + n * 13.5, l, 9.5))
        pages.append(_page(parts))
    return pages


def build(path, chains, out_pdf, model_index=0, title=None, models=None):
    cairosvg, pypdf = _need()
    name = title or str(path).replace("\\", "/").split("/")[-1]
    pages, summ_lines = [], []
    per_chain = []
    for c in chains:
        rs, rrows = rama.analyze(path, model_index, [c], models)
        ts, trows = rotamer.analyze(path, model_index, [c], models)
        fs, _ = flips.analyze(path, model_index, [c])
        per_chain.append((c, rs, rrows, ts, trows, fs))
    from . import __version__
    summ_lines += [f"ifacekit {__version__}. Chains: {', '.join(chains)}.", "",
                   f"Ramachandran reference: {per_chain[0][1].get('reference', '')}", f"Rotamer reference: {per_chain[0][3].get('reference', '')}", ""]
    for c, rs, _, ts, _, fs in per_chain:
        rc, rp = rs.get("counts", {}), rs.get("percent", {})
        tc, tp = ts.get("counts", {}), ts.get("percent", {})
        summ_lines.append(f"Chain {c}")
        if rs.get("n"):
            summ_lines.append(f"   Ramachandran: {rs['n']} residues; favored {rc['favored']} ({rp['favored']}%), allowed {rc['allowed']} ({rp['allowed']}%), "
                              f"outliers {rc['outlier']} ({rp['outlier']}%)   goals: favored over 98%, outliers under 0.2%")
        if ts.get("n"):
            summ_lines.append(f"   Rotamers: {ts['n']} side chains; favored {tc['favored']} ({tp['favored']}%), allowed {tc['allowed']} ({tp['allowed']}%), "
                              f"outliers {tc['outlier']} ({tp['outlier']}%)   goal: outliers under 0.3%")
        if fs.get("n"):
            fc = fs["counts"]
            summ_lines.append(f"   Asn/Gln/His flip check: {fs['n']} groups; flip suggested {fc['flip']}, uncertain {fc['uncertain']}, keep {fc['keep']} (a weak test, see the design document)")
        summ_lines.append("")
    summ_lines += ["How to read the plots", "   Ramachandran: shaded = favored (dark) and allowed (light) regions of the ifacekit reference; one panel per residue class.",
                   "   Rotamers: chi1 against chi2 per residue type, shaded with the highest reference density over the other chi angles;",
                   "   one chi angle: the reference density along chi1. Points are coloured by category in the full space:",
                   "   green favored, orange allowed, red outlier."]
    pages += _lines_pages(f"{name}: Ramachandran and rotamer plots", summ_lines)
    for c, rs, rrows, ts, trows, fs in per_chain:
        pan = []
        for cls in rama.CLASSES:
            sub = [r for r in rrows if r["class"] == cls]
            if sub:
                k = {x: sum(1 for r in sub if r["category"] == x) for x in ("favored", "allowed", "outlier")}
                pan.append((f"{rama.LABEL[cls]}: {len(sub)} res, {k['outlier']} outlier", ramaplot._panel(cls, sub, PANEL)))
        pages += _grid_pages(f"Chain {c}: Ramachandran", name, pan)
        pan = []
        for t in rotamer.TABLES:
            sub = [r for r in trows if r["table"] == t]
            if sub:
                k = {x: sum(1 for r in sub if r["category"] == x) for x in ("favored", "allowed", "outlier")}
                pan.append((f"{rotamer.LABEL[t]}: {len(sub)} res, {k['outlier']} outlier", rotamerhtml._panel(t, sub, PANEL)))
        pages += _grid_pages(f"Chain {c}: rotamers", name, pan)
        lines = []
        if rs.get("outliers"):
            lines += ["Ramachandran outliers (residue, class, phi, psi)"] + [f"   {o['residue']:<14s} {rama.LABEL[o['class']]:<10s} {o['phi']:>7} {o['psi']:>7}" for o in rs["outliers"]] + [""]
        if ts.get("outliers"):
            lines += ["Rotamer outliers (residue, type, chi angles, category if flipped)"] + [
                f"   {rotamer.label(o):<18s} {rotamer.LABEL[o['table']]:<8s} {' '.join(str(x) for x in o['chi']):<28s} {o.get('category_flipped', '-')}" for o in ts["outliers"]] + [""]
        if fs.get("suggested"):
            lines += ["Asn/Gln/His flips suggested (residue, score as modelled, score flipped)"] + [
                f"   {o['residue']:<14s} {o['score_modelled']:>6} {o['score_flipped']:>6}" for o in fs["suggested"]]
        if lines:
            pages += _lines_pages(f"Chain {c}: outliers and suggested flips", lines)
    writer = pypdf.PdfWriter()
    for svg in pages:
        writer.append(pypdf.PdfReader(io.BytesIO(cairosvg.svg2pdf(bytestring=svg.encode("utf-8")))))
    writer.add_metadata({"/Title": f"{name}: Ramachandran and rotamer plots", "/Creator": f"ifacekit {__version__}"})
    with open(out_pdf, "wb") as fh:
        writer.write(fh)
    return len(pages)
