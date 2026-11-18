"""HTML reports of the models of a multi-model file (NMR ensemble).

One model:      the normal report of that model (chain quality or interface report) with the plots of the ensemble.
Several models: one page. It keeps the frame of the normal report (banner, title, notes) and has a bar with a button for the ensemble summary and one for each model; a button
                swaps the tables and sections below it for those of that model. The summary holds the table of all models with mean and standard deviation, the band plot of the
                deviation of each residue, and the Ramachandran and rotamer plots with every model in one plot. With molstar=True the page has the Mol* panel (one structure file with
                all models): it follows the model chosen, and a button shows all models together. separate=True writes one file per model instead.
"""
import json
import os
from concurrent.futures import ProcessPoolExecutor

from . import htmlreport as hr
from . import qualityreport as qr
from .about import footer_html
from .banner import banner_html  # noqa: F401
from .ensemble import COLUMNS, parse_models  # noqa: F401  (parse_models is re-exported for the command line)

EXTRA_CSS = """
.etabs{display:flex;flex-wrap:wrap;gap:6px;margin:10px 0 14px;position:sticky;top:0;z-index:5;background:var(--bg);padding:6px 0}
.etabs button{font:inherit;padding:4px 10px;border:1px solid var(--line);border-radius:6px;background:var(--card);color:var(--fg);cursor:pointer}
.etabs button.on{background:var(--acc,#2c6fbb);color:#fff;border-color:transparent}
.etabs .sp{flex:1}#ens-body nav.top .theme,#ens-body nav.top .sp{display:none}
"""


def _collect(path, chains, pair, index, n_models, hydrogens, molstar, molstar_url):
    if pair is not None:
        d = hr.collect(path, pair[0], pair[1], None, None, index, hydrogens, True)
        d["kind"] = "experimental"
        d["pdf"] = None
    else:
        d = qr.collect(path, chains, index, False, None, hydrogens, True)
    d["molstar_url"] = molstar_url
    d["n_models"] = n_models
    return d


def _render(d, pair, molstar):
    return hr.render(d, molstar=molstar) if pair is not None else qr.render(d, molstar=molstar)


def _page(args):
    """One model: the complete report (for a single model) and the part of it that is specific to the model (the tables and sections), as strings. Runs in a worker process."""
    path, chains, pair, index, n_models, hydrogens, molstar, molstar_url, single = args
    d = _collect(path, chains, pair, index, n_models, hydrogens, molstar, molstar_url)
    d["ensemble_page"] = not single
    html = _render(d, pair, molstar)
    out = {"index": index}
    if single:
        note = (f"<p class='note'><b>Model {index + 1} of {n_models}</b> of the ensemble. The plot of the deviations and, in the Ramachandran and rotamer sections, "
                "the numbers of this model are shown; the summary of all models is in the ensemble page.</p>")
        i = html.find("</h1>")
        out["html"] = (html[:i + 5] + note + html[i + 5:] if i >= 0 else note + html).replace("<title>", f"<title>model {index + 1} - ", 1)
        return out
    head, _, rest = html.partition(hr.BODY_START)
    body, _, _ = rest.partition(hr.BODY_END)
    out["head"], out["body"] = head, body
    it = d.get("iface")
    out["iface"] = [[r["residue"].split(":")[0], int(r["residue"].split(":")[1]), r["residue"].split(":")[2]] for key in ("residues_a", "residues_b") for r in (it[key] if it else [])]
    out["residues"] = d["residues"]
    out["data"] = {k: d[k] for k in ("path", "residues", "iface", "kind", "molstar_url", "n_models", "model_index") if k in d}
    out["data"].update({k: d[k] for k in ("a", "b") if k in d})
    return out


def _summary_html(ens_list, chains, pair, path=None, n_models=1):
    out = []
    for e in ens_list:
        cols = [(k, lab) for k, lab in COLUMNS if any(k in r for r in e["models"])]
        av = e["average"]
        rows = [[r["model"]] + [("-" if r.get(k) is None else r.get(k)) for k, _ in cols] for r in e["models"]]
        rows.append(["mean"] + [av.get(k, {}).get("mean", "-") for k, _ in cols])
        rows.append(["sd"] + [av.get(k, {}).get("sd", "-") for k, _ in cols])
        out.append(f"<h3>{hr.esc(e['target'])}: {len(e['models'])} of {e['n_models_in_file']} models</h3>" + hr.table(["model"] + [lab for _, lab in cols], rows))
    res = [r for e in ens_list for r in e["residue_rmsd"]]
    names = list(dict.fromkeys(r["chain"] for r in res))
    if res:
        out.append("<h3 id='plddt'>Backbone deviation along the sequence</h3>" + hr.plddt_plot(res, names, set(), True, label="deviation (A)")
                   + "<p class='legend'>Line: mean over the models of the distance of each residue from the mean structure (root of the mean squared distance of the backbone atoms N, CA, C, O to their mean positions, "
                     "each chain superposed on its least mobile 70% of residues). Dark band: mean &plusmn; 1 standard deviation; light band: lowest to highest model; the RMSD is in the tooltip. "
                     "Scroll over the plot to zoom, drag to pan, and hover to read a residue.</p>")
    if path:        # the Ramachandran and rotamer plots with every model of the ensemble in one plot
        from .ramaplot import fragment as rama_fragment
        from .rotamerhtml import fragment as rot_fragment
        allm = list(range(n_models))
        for grp in ([[c] for c in chains] if chains else [list(pair)]):
            out.append(rama_fragment(path, grp, 0, anchor=f"ens-rama-{'-'.join(grp)}", models=allm))
            out.append(rot_fragment(path, grp, 0, anchor=f"ens-rot-{'-'.join(grp)}", models=allm, with_flips=False))
    return "".join(out)


JS = ("(function(){var P=JSON.parse(document.getElementById('ens-data').textContent),B=document.querySelectorAll('.etabs button[data-v]'),C=document.getElementById('ens-body');"
      "function put(h){C.innerHTML=h;[].forEach.call(C.querySelectorAll('script'),function(o){var n=document.createElement('script');n.text=o.textContent;o.parentNode.replaceChild(n,o)})}"
      "function show(v){B.forEach(function(b){b.classList.toggle('on',b.dataset.v===v)});"
      "if(v==='sum'){put(P.sum)}else{put(P.m[v].body);if(window.IFK&&window.IFK.setModel)window.IFK.setModel(parseInt(v,10)-1,P.m[v].iface,null)}"
      "try{history.replaceState(null,'','#'+(v==='sum'?'summary':'model'+v))}catch(e){}}"
      "B.forEach(function(b){b.onclick=function(){show(b.dataset.v)}});"
      "var h=location.hash.match(/^#model(\\d+)$/);show(h&&P.m[h[1]]?h[1]:'sum')})();")


def combined(path, ens_list, parts, n_models, chains=None, pair=None, molstar=False, molstar_url=None):
    """The ensemble page: the frame of the report of the first chosen model, the bar of buttons, the model-specific part, and (optionally) the Mol* panel."""
    first = parts[0]
    head = first["head"].replace("<title>", "<title>ensemble - ", 1)
    summary = _summary_html(ens_list, chains, pair, path, n_models)
    data = {"sum": summary, "m": {str(p["index"] + 1): {"body": p["body"], "iface": p["iface"]} for p in parts}}
    blob = json.dumps(data).replace("</", "<\\/")
    buttons = ["<button class='on' data-v='sum'>Ensemble summary</button>"] + [f"<button data-v='{p['index'] + 1}'>Model {p['index'] + 1}</button>" for p in parts]
    bar = ("<p class='note'>NMR ensemble of " + str(n_models) + " models. Choose the summary (every model in one table, mean and standard deviation, the deviation of every residue, "
           "the Ramachandran and rotamer plots with all models) or one model: the tables and sections below change to those of that model"
           + (", and the 3D view follows it (the button <b>all models</b> in the 3D panel shows them together)" if molstar else "") + ".</p>"
           f"<style>{EXTRA_CSS}</style><div class='etabs'>" + "".join(buttons)
           + "<span class='sp'></span><button type='button' class='theme' onclick=\"var r=document.documentElement;r.dataset.theme=(r.dataset.theme==='dark'?'light':'dark')\">light / dark</button>"
           + "</div><div id='ens-body'></div>")
    tail = f"<script type='application/json' id='ens-data'>{blob}</script><script>{JS}</script>" + footer_html()
    if molstar:
        dm = dict(first["data"])
        dm["molstar_url"] = molstar_url
        dm["n_models"] = n_models
        dm["model_index"] = first["index"]
        tail += hr.molstar_block(dm)
    return head + bar + tail + "</main></body></html>"


def build(path, ens_list, out_path, chains=None, pair=None, indices=None, separate=False, hydrogens=False, workers=1, molstar=False, molstar_url=None):
    """Write the HTML of the chosen models (indices counted from 0). Returns the list of files written.
    With molstar the page is also written with the Mol* panel, as FILE_3d.html (the plain page without Mol* is always written)."""
    n = ens_list[0]["n_models_in_file"]
    indices = [0] if indices is None else list(indices)
    base = out_path[:-5] if out_path.lower().endswith(".html") else out_path
    files = []
    for mol in ((False, True) if molstar else (False,)):
        single = len(indices) == 1 or separate
        jobs = [(path, chains, pair, i, n, hydrogens, mol, molstar_url, single) for i in indices]
        if workers and workers > 1 and len(jobs) > 1:
            with ProcessPoolExecutor(max_workers=workers) as ex:
                parts = list(ex.map(_page, jobs))
        else:
            parts = [_page(j) for j in jobs]
        suffix = "_3d" if mol else ""
        if len(indices) == 1:
            f = out_path if not mol else f"{base}_3d.html"
            open(f, "w", encoding="utf-8").write(parts[0]["html"])
            files.append(f)
        elif separate:
            for p in parts:
                f = f"{base}_model{p['index'] + 1:02d}{suffix}.html"
                open(f, "w", encoding="utf-8").write(p["html"])
                files.append(f)
        else:
            f = out_path if not mol else f"{base}_3d.html"
            open(f, "w", encoding="utf-8").write(combined(path, ens_list, parts, n, chains, pair, mol, molstar_url))
            files.append(f)
    return files
