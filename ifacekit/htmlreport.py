"""One self-contained HTML page for a designed or predicted complex.

Three sections: AI confidence metrics, interface interactions, prediction quality.
No external files or libraries are used, so the page can be opened offline or published as is.
"""
import base64
import html
import json

import numpy as np

from . import af3
from .about import comparison_html
from .banner import banner_html
from .core import analyze

CSS = """
:root{--bg:#fbfbfa;--fg:#1d2326;--mut:#667078;--card:#ffffff;--line:#dde1e3;--acc:#2b6cb0;--good:#2f855a;--warn:#b7791f;--bad:#c0392b;--chip:#eef1f3;--badbg:#fdf0ee;--t1:#eef4fb;--t2:#eef7f1;--t3:#fbf5e9}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#14181a;--fg:#e6e9ea;--mut:#98a2a8;--card:#1c2124;--line:#2e3539;--acc:#63a4e8;--good:#5cc08a;--warn:#e3b04b;--bad:#ef6f61;--chip:#262d31;--badbg:#341c19;--t1:#182430;--t2:#18281f;--t3:#2a2417}}
:root[data-theme=dark]{--bg:#14181a;--fg:#e6e9ea;--mut:#98a2a8;--card:#1c2124;--line:#2e3539;--acc:#63a4e8;--good:#5cc08a;--warn:#e3b04b;--bad:#ef6f61;--chip:#262d31;--badbg:#341c19;--t1:#182430;--t2:#18281f;--t3:#2a2417}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
main{max-width:1080px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:18px;margin:36px 0 12px;padding-top:8px;border-top:1px solid var(--line)}h3{font-size:14px;margin:20px 0 8px;color:var(--mut);font-weight:600}
.sub{color:var(--mut);margin:0 0 16px}
.cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:10px}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px 12px}
.card .k{font-size:12px;color:var(--mut)}.card .v{font-size:22px;font-weight:600;font-variant-numeric:tabular-nums}.card .n{font-size:12px;color:var(--mut)}
table{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums}th,td{text-align:left;padding:5px 8px;border-bottom:1px solid var(--line)}th{color:var(--mut);font-weight:500}
.wrap{overflow-x:auto}
.tag{display:inline-block;padding:1px 8px;border-radius:10px;font-size:12px;background:var(--chip)}
.card.st-good{border-left:5px solid var(--good)}.card.st-warn{border-left:5px solid var(--warn)}.card.st-bad{border-left:5px solid var(--bad);background:var(--badbg)}
.chips a.st-good{border-left:5px solid var(--good)}.chips a.st-warn{border-left:5px solid var(--warn)}.chips a.st-bad{border-left:5px solid var(--bad);background:var(--badbg)}
.attn{border:1px solid var(--line);border-left:6px solid var(--good);border-radius:8px;background:var(--card);padding:8px 14px;margin:10px 0;font-size:14px}.attn.bad{border-left-color:var(--bad)}.attn.warn{border-left-color:var(--warn)}.attn ul{margin:4px 0 0 18px;padding:0}.attn a{color:var(--acc)}
.favorable,.ok{color:var(--good)}.borderline{color:var(--warn)}.unfavorable{color:var(--bad)}
.note{color:var(--mut);font-size:13px;margin:8px 0}
.row{display:grid;grid-template-columns:1fr 1fr;gap:16px}@media(max-width:760px){.row{grid-template-columns:1fr}}
canvas,svg{max-width:100%;height:auto}
.bar{display:flex;height:14px;border-radius:7px;overflow:hidden;background:var(--chip)}
.legend{font-size:12px;color:var(--mut)}
html{scroll-behavior:smooth}
nav.top{position:sticky;top:0;z-index:10;background:var(--bg);border-bottom:1px solid var(--line);margin:0 -16px 16px;padding:8px 16px;display:flex;gap:6px;flex-wrap:wrap;align-items:center}
nav.top a{color:var(--fg);text-decoration:none;padding:5px 12px;border-radius:16px;font-size:14px;border:1px solid var(--line);background:var(--card)}
nav.top a.sub{font-size:12px;color:var(--mut);padding:3px 9px}
nav.top a.active{background:var(--acc);color:#fff;border-color:var(--acc)}
nav.top .sp{flex:1}
.plotbar{display:flex;gap:6px;align-items:center;margin-bottom:6px}.plotbar button{font-size:12px;padding:3px 10px;border:1px solid var(--line);background:var(--card);color:var(--fg);border-radius:6px;cursor:pointer}
.chips{display:flex;gap:8px;flex-wrap:wrap;margin:6px 0 12px}
.chips a{text-decoration:none;color:var(--fg);background:var(--card);border:1px solid var(--line);border-radius:8px;padding:6px 10px;font-size:13px}
.chips a b{font-size:16px;font-variant-numeric:tabular-nums}.chips a span{color:var(--mut);margin-left:4px}
.chips a:hover{border-color:var(--acc)}
section.sec{border-radius:12px;padding:4px 18px 18px;margin:20px 0;border:1px solid var(--line);border-left-width:6px}
section.sec h2{border-top:none;margin-top:14px}
section.s1{background:var(--t1);border-left-color:#2b6cb0}section.s2{background:var(--t2);border-left-color:#2f855a}section.s3{background:var(--t3);border-left-color:#b7791f}
section.sec .card{background:var(--card)}
:target{scroll-margin-top:64px}h2,h3{scroll-margin-top:64px}
button.theme{float:right;background:var(--chip);color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:4px 10px;cursor:pointer}
"""
BAND_COL = ["#0053d6", "#65cbf3", "#ffdb13", "#ff7d45"]
CHAIN_COL = ["#2b6cb0", "#d9822b", "#2f855a", "#805ad5", "#c0392b", "#2c7a7b", "#b7791f", "#6b46c1"]   # one per chain in the sequence plot


def esc(x):
    return html.escape(str(x))


def st_hi(x, good, warn):
    """status for a value where high is good"""
    if x is None:
        return ""
    return "good" if x >= good else "warn" if x >= warn else "bad"


def st_lo(x, good, warn):
    if x is None:
        return ""
    return "good" if x <= good else "warn" if x <= warn else "bad"


def st_probs(items):
    """status from a list of contacts: any unfavorable -> bad, any borderline -> warn"""
    if any(c["verdict"] == "unfavorable" for c in items):
        return "bad"
    return "warn" if items else "good"


def card(k, v, n="", st=""):
    return f'<div class="card {("st-" + st) if st else ""}"><div class="k">{esc(k)}</div><div class="v">{esc(v)}</div><div class="n">{esc(n)}</div></div>'


def fmt(x, nd=2):
    return "-" if x is None else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


_MS = [False]   # set while rendering the Mol* page
BODY_START = "<!--IFK-BODY-->"        # markers around the part of a report that belongs to one model (the ensemble page swaps it)
BODY_END = "<!--IFK-BODY-END-->"


def table(head, rows, cls=None, atoms=None):
    """atoms: optional list parallel to rows; each item is a list of one or two atom or residue ids to show in Mol*."""
    out = ['<div class="wrap"><table' + (' class="ms"' if _MS[0] and atoms else "") + '><thead><tr>' + "".join(f"<th>{esc(h)}</th>" for h in head) + "</tr></thead><tbody>"]
    for k, r in enumerate(rows):
        at = f' data-a="{esc("|".join(atoms[k]))}"' if _MS[0] and atoms and atoms[k] else ""
        out.append(f"<tr{at}>" + "".join(
            (f'<td class="{cls[i](c)}">{esc(c)}</td>' if cls and cls.get(i) else f"<td>{esc(c)}</td>") for i, c in enumerate(r)) + "</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def plddt_plot(rows, chains, iface, exp=False, label=None):
    """Interactive per-residue pLDDT plot: wheel to zoom, drag to pan, hover for the residue.
    label: name of the plotted number when it is not a B-factor or pLDDT (the RMSD of an NMR ensemble, in A)."""
    if not rows:
        return ""
    rr = [r for r in rows if r["chain"] in chains]
    order = []
    for r in rr:
        if r["chain"] not in order:
            order.append(r["chain"])
    ymax = 100
    band = bool(rr) and all("dev_sd" in r for r in rr)
    if label:
        top = max([r.get("dev_max", r["b"]) for r in rr] + [1.0])
        ymax = max(4, int(-(-top // 4) * 4))          # a multiple of 4 A, so that the quarter lines carry whole numbers
    elif exp:
        top = max([r["b"] for r in rr] + [10])
        ymax = int(-(-top // 10) * 10)
    data = {"y": ymax, "x": 1 if exp else 0, "c": [r["chain"] for r in rr], "n": [r["resnum"] for r in rr], "t": [r["resname"] for r in rr],
            "b": [r["b"] for r in rr], "i": [1 if (r["chain"], r["resnum"]) in iface else 0 for r in rr],
            "k": [order.index(r["chain"]) for r in rr], "l": label or ""}
    if band:        # mean line with a band: +-1 SD (dark) and the range of the models (light)
        data.update({"lo": [r["dev_min"] for r in rr], "hi": [r["dev_max"] for r in rr], "s1": [max(r["dev_mean"] - r["dev_sd"], 0) for r in rr],
                     "s2": [r["dev_mean"] + r["dev_sd"] for r in rr], "sd": [r["dev_sd"] for r in rr], "r": [r["rmsd_backbone"] or 0 for r in rr]})
    return ("<div class='plotbox'><div class='plotbar'><button type='button' id='pl-in'>zoom in</button><button type='button' id='pl-out'>zoom out</button>"
            "<button type='button' id='pl-reset'>reset</button><span id='pl-info' class='note'>Scroll to zoom, drag to pan, hover for the residue.</span></div>"
            "<canvas id='plcv' height='268' style='width:100%;height:268px;cursor:grab;border:1px solid var(--line);border-radius:6px;background:var(--card)'></canvas></div>"
            "<script>(function(){var D=" + json.dumps(data) + ",C=" + json.dumps(BAND_COL) + ",P=" + json.dumps(CHAIN_COL) + ";var cv=document.getElementById('plcv'),x=cv.getContext('2d'),N=D.b.length;"
            "var v0=0,v1=N-1,W,H=268,L=36,B=52,T=8,hv=-1,dr=null;"
            "function css(k){return getComputedStyle(document.documentElement).getPropertyValue(k).trim()}"
            "function px(i){return L+(W-L-8)*(i-v0)/Math.max(v1-v0,1)}function py(v){return T+(H-B-T)*(1-v/D.y)}"
            "function draw(){var r=window.devicePixelRatio||1;W=cv.clientWidth;cv.width=W*r;cv.height=H*r;x.setTransform(r,0,0,r,0,0);x.clearRect(0,0,W,H);var fg=css('--fg'),mut=css('--mut');"
            "x.font='10px sans-serif';x.fillStyle=mut;x.textAlign='right';(D.x?[[D.y,mut],[D.y*.75,mut],[D.y*.5,mut],[D.y*.25,mut]]:[[90,C[0]],[70,C[1]],[50,C[2]]]).forEach(function(a){x.strokeStyle=a[1];x.setLineDash([3,4]);x.beginPath();x.moveTo(L,py(a[0]));x.lineTo(W-8,py(a[0]));x.stroke();x.fillText(Math.round(a[0]),L-4,py(a[0])+3)});x.setLineDash([]);"
            "x.textAlign='left';var a0=Math.max(0,Math.floor(v0)),a1=Math.min(N-1,Math.ceil(v1));"
            "for(var i=a0;i<=a1;i++){if(D.i[i]){x.fillStyle='rgba(192,57,43,.13)';var w=Math.max((W-L-8)/Math.max(v1-v0,1),1);x.fillRect(px(i)-w/2,T,w,H-B-T)}}"
            "if(D.lo){[[D.lo,D.hi,'rgba(43,108,176,.14)'],[D.s1,D.s2,'rgba(43,108,176,.32)']].forEach(function(b){x.fillStyle=b[2];for(var j=a0;j<a1;j++){if(D.c[j]!==D.c[j+1])continue;x.beginPath();x.moveTo(px(j),py(b[0][j]));x.lineTo(px(j+1),py(b[0][j+1]));x.lineTo(px(j+1),py(b[1][j+1]));x.lineTo(px(j),py(b[1][j]));x.closePath();x.fill()}})}"
            "x.lineWidth=1.8;var st=true;x.beginPath();for(i=a0;i<=a1;i++){var X=px(i),Y=py(D.b[i]);if(i>a0&&D.c[i]!==D.c[i-1]){x.strokeStyle=P[D.k[i-1]%P.length];x.stroke();x.beginPath();st=true}if(st){x.moveTo(X,Y);st=false}else x.lineTo(X,Y)}x.strokeStyle=P[D.k[a1]%P.length];x.stroke();"
            "var k=Math.max(1,Math.round((v1-v0)/12));x.fillStyle=mut;x.textAlign='center';for(i=a0;i<=a1;i++){if(i===a0||D.c[i]!==D.c[i-1]){x.strokeStyle=mut;x.globalAlpha=.6;x.lineWidth=2;x.beginPath();x.moveTo(px(i),T);x.lineTo(px(i),H-B+30);x.stroke();x.globalAlpha=1;x.lineWidth=1}"
            "if((i-a0)%k===0||D.i[i]&&(v1-v0)<80){x.fillText(D.n[i],px(i),H-B+13)}}"
            "var seg=[];for(i=a0;i<=a1;i++){if(i===a0||D.c[i]!==D.c[i-1])seg.push([i,i]);else seg[seg.length-1][1]=i}"
            "seg.forEach(function(s){var wi=Math.max((W-L-8)/Math.max(v1-v0,1),1),xa=Math.max(L,px(s[0])-wi/2),xb=Math.min(W-8,px(s[1])+wi/2),col=P[D.k[s[0]]%P.length];x.fillStyle=col;x.fillRect(xa,H-B+20,xb-xa,10);x.font='bold 15px sans-serif';x.textAlign='center';x.fillStyle=col;x.fillText((xb-xa)>80?'Chain '+D.c[s[0]]:D.c[s[0]],(xa+xb)/2,H-B+46);x.font='10px sans-serif'});x.textAlign='center';x.fillStyle=mut;"
            "for(i=a0;i<=a1;i++){if(D.i[i]){x.fillStyle=css('--bad');x.fillRect(px(i)-1.5,H-B-5,3,5)}}"
            "if(hv>=0&&hv>=a0&&hv<=a1){x.strokeStyle=fg;x.lineWidth=1;x.beginPath();x.moveTo(px(hv),T);x.lineTo(px(hv),H-B);x.stroke();x.fillStyle=css('--bad');x.beginPath();x.arc(px(hv),py(D.b[hv]),3.5,0,7);x.fill()}x.lineWidth=1}"
            "function info(i){document.getElementById('pl-info').textContent=i<0?'Scroll to zoom, drag to pan, hover for the residue.':D.c[i]+':'+D.t[i]+D.n[i]+(D.l?'   '+D.l+' ':D.x?'   B-factor ':'   pLDDT ')+D.b[i].toFixed(D.l?2:1)+(D.lo?' \u00b1 '+D.sd[i].toFixed(2)+' (range '+D.lo[i].toFixed(2)+' to '+D.hi[i].toFixed(2)+'), RMSD '+D.r[i].toFixed(2):'')+(D.i[i]?'   interface residue':'')}"
            "function at(e){var r=cv.getBoundingClientRect();return v0+(e.clientX-r.left-L)/(W-L-8)*(v1-v0)}"
            "function zoom(f,c){var n0=c-(c-v0)*f,n1=c+(v1-c)*f;if(n1-n0<8){var m=(n0+n1)/2;n0=m-4;n1=m+4}if(n1-n0>N-1){n0=0;n1=N-1}if(n0<0){n1-=n0;n0=0}if(n1>N-1){n0-=n1-(N-1);n1=N-1}v0=Math.max(0,n0);v1=n1;draw()}"
            "cv.addEventListener('wheel',function(e){e.preventDefault();zoom(e.deltaY<0?.8:1.25,at(e))},{passive:false});"
            "var dx=0;cv.addEventListener('mousedown',function(e){dr={x:e.clientX,a:v0,b:v1};dx=e.clientX;cv.style.cursor='grabbing'});window.addEventListener('mouseup',function(e){if(dr&&Math.abs(e.clientX-dx)<4&&hv>=0&&window.IFK)window.IFK.show({a:[[D.c[hv]+':'+D.n[hv]+':'+D.t[hv]]]});dr=null;cv.style.cursor='grab'});"
            "window.addEventListener('mousemove',function(e){if(dr){var d=(e.clientX-dr.x)/(W-L-8)*(dr.b-dr.a),a=dr.a-d,b=dr.b-d;if(a<0){b-=a;a=0}if(b>N-1){a-=b-(N-1);b=N-1}v0=Math.max(0,a);v1=b;draw()}});"
            "cv.addEventListener('mousemove',function(e){var i=Math.round(at(e));hv=(i>=0&&i<N)?i:-1;info(hv);draw()});cv.addEventListener('mouseleave',function(){hv=-1;info(-1);draw()});"
            "document.getElementById('pl-in').onclick=function(){zoom(.6,(v0+v1)/2)};document.getElementById('pl-out').onclick=function(){zoom(1.67,(v0+v1)/2)};document.getElementById('pl-reset').onclick=function(){v0=0;v1=N-1;draw()};"
            "window.addEventListener('resize',draw);new MutationObserver(draw).observe(document.documentElement,{attributes:true});matchMedia('(prefers-color-scheme:dark)').addEventListener('change',draw);draw();})();</script>")


def pae_canvas(full, chains):
    pae = full["pae"]
    n = pae.shape[0]
    q = np.clip(np.round(pae * 8), 0, 255).astype(np.uint8)
    b64 = base64.b64encode(q.tobytes()).decode()
    ch = full["chain"]
    bounds = [int(i) for i in range(1, n) if ch[i] != ch[i - 1]]
    return (f'<canvas id="paecv" width="{n}" height="{n}" style="width:100%;max-width:460px;image-rendering:pixelated;border:1px solid var(--line)"></canvas>'
            f'<script>(function(){{var n={n},b="{b64}",bd={json.dumps(bounds)};var s=atob(b),c=document.getElementById("paecv"),x=c.getContext("2d"),im=x.createImageData(n,n);'
            'for(var i=0;i<n*n;i++){var t=s.charCodeAt(i)/8/31.75;t=Math.min(1,t);var r=Math.round(20+235*t),g=Math.round(75+170*t),bb=Math.round(140+95*t);'
            'im.data[4*i]=r;im.data[4*i+1]=g;im.data[4*i+2]=bb;im.data[4*i+3]=255;}x.putImageData(im,0,0);'
            'x.strokeStyle="#e03c31";x.lineWidth=1;bd.forEach(function(p){x.beginPath();x.moveTo(p,0);x.lineTo(p,n);x.moveTo(0,p);x.lineTo(n,p);x.stroke();});})();</script>')


def band_bar(b):
    parts = [("very_high_ge90", BAND_COL[0], ">=90"), ("confident_70_90", BAND_COL[1], "70-90"),
             ("low_50_70", BAND_COL[2], "50-70"), ("very_low_lt50", BAND_COL[3], "<50")]
    segs = "".join(f'<div title="{lab} {b[k] * 100:.0f}%" style="width:{b[k] * 100:.1f}%;background:{col}"></div>' for k, col, lab in parts)
    return f'<div class="bar">{segs}</div>'


CLASH_REASONS = ("steric_overlap", "close_contact")


def _bridge_check_text(x):
    from .hbonds import FLAG_TEXT
    return "clear" if not x["flags"] else "doubtful: " + "; ".join(FLAG_TEXT[f] for f in x["flags"])


MEDIATED_NOTE = ("Each bridge is checked on its closest partner atom of each side: an H-bond shorter or longer than the limits, partner atoms that touch each "
                 "other directly (the water only sits beside that contact), or a narrow angle at the water make it doubtful. The limits are simple and not calibrated "
                 "(see ifacekit/hbonds.py), so read a doubtful bridge as a reason to look at the density, not as an error. A bridge through two waters is not counted. "
                 "Water positions in X-ray models are less certain than atoms of the protein.")


def _mediated_html(it, a, b, n_waters):
    """Table of the waters that bridge two sides by hydrogen bonds (ifacekit.hbonds)."""
    h = it["hbonds"]
    head = (f"<h3 id='mediated'>Hydrogen bonds via water: water bridges ({h['water_mediated']}"
            + (f": {h['water_mediated_clear']} clear, {h['water_mediated_doubtful']} doubtful" if h["water_mediated"] else "") + ")</h3>")
    if not n_waters:
        return head + "<p class='note'>This model has no waters (as in an AI prediction), so no bonds via water can be counted.</p>"
    if not h["mediated_waters"]:
        return head + "<p class='note'>No water makes a hydrogen bond to both sides.</p>"
    rows = [[x["water"], ", ".join(f"{y['atom']} {y['distance']}" for y in x["side_a"]), ", ".join(f"{y['atom']} {y['distance']}" for y in x["side_b"]),
             _bridge_check_text(x)] for x in h["mediated_waters"]]
    atoms = [[x["water"]] + [y["atom"] for y in x["side_a"] + x["side_b"]] for x in h["mediated_waters"]]
    return (head + f"<p class='note'><b>This table lists the bridging waters only</b>: waters that make a hydrogen bond to an atom of {esc(a)} <b>and</b> to an atom of {esc(b)} "
            f"(distance in &Aring; after the atom). A water bonded to only one of the two chains is not an interface bond and is not listed or counted here; "
            "it belongs to the hydration of that chain (see the model quality part). "
            + MEDIATED_NOTE + "</p>"
            + table(["water", f"bonded to {a}", f"bonded to {b}", "check"], rows, cls={3: lambda v: "borderline" if v.startswith("doubtful") else "ok"}, atoms=atoms))


def split_problems(contacts, only_between=False, units=None):
    """Atom pairs with a borderline or unfavorable verdict, each pair once, split into
    clashes (space) and chemistry problems (hydrogen-bond distance, polarity, charge)."""
    clashes, chem, seen = [], [], set()
    for c in contacts:
        if c["verdict"] not in ("borderline", "unfavorable"):
            continue
        if units is not None:                         # only contacts between exactly these two units (a water is never one of them)
            if {c["unit1"], c["unit2"]} != set(units):
                continue
        elif only_between and not c["inter_unit"]:      # between any two units, a water counts as a unit here
            continue
        key = frozenset([c["atom1"], c["atom2"]])
        if key in seen:
            continue
        seen.add(key)
        (clashes if c["reason"] in CLASH_REASONS else chem).append(c)
    clashes.sort(key=lambda c: -c["overlap"])
    chem.sort(key=lambda c: (c["verdict"] != "unfavorable", c["distance"]))
    return clashes, chem


def clash_table(rows, limit=30):
    if not rows:
        return "<p class='note'>none</p>"
    alt = any(c.get("altloc") for c in rows[:limit])
    return table(["atom", "atom", "distance (Å)", "overlap (Å)", "verdict"] + (["conformer"] if alt else []),
                 [[c["atom1"], c["atom2"], c["distance"], f"{c['overlap']:.2f}", c["verdict"]] + ([c.get("altloc", "")] if alt else []) for c in rows[:limit]],
                 cls={4: lambda v: v}, atoms=[[c["atom1"], c["atom2"]] for c in rows[:limit]]) + (f"<p class='note'>{len(rows) - limit} more in the JSON</p>" if len(rows) > limit else "")


HYDROGEN_NOTE = ("Clashes found after ifacekit placed hydrogens on the protein; every other table of this report uses the model as given. "
                 "A clash is an overlap of 0.4 &Aring; or more between atoms that are more than four bonds apart (radii of the Probe program used by MolProbity); "
                 "a polar hydrogen next to an acceptor is a hydrogen bond, not a clash. Hydrogens of ligands, waters, nucleic acids and modified residues are not placed, "
                 "so their clashes with hydrogens are not seen, and only the first conformer of every atom is checked. The settings were calibrated against MolProbity (see the fact sheet).")


def hclash_html(h, units=None, limit=30, anchor=None, what="All-atom clashes with hydrogens added"):
    """Table of the clashes found with hydrogens added (ifacekit.hydrogens); units: only clashes between these two units."""
    if not h:
        return ""
    rows = h["clashes"]
    if units is not None:
        rows = [c for c in rows if {c["unit1"], c["unit2"]} == set(units)]
    idattr = f" id='{anchor}'" if anchor else ""
    head = f"<h3{idattr}>{what} ({len(rows)})</h3><p class='note'>{HYDROGEN_NOTE}</p>"
    if not rows:
        return head + "<p class='note'>none</p>"
    cut = rows[:limit]
    return head + table(["atom", "atom", "distance (Å)", "overlap (Å)", "verdict"],
                        [[c["atom1"], c["atom2"], c["distance"], f"{c['overlap']:.2f}", "unfavorable" if c["overlap"] >= 0.6 else "borderline"] for c in cut],
                        cls={4: lambda v: v}, atoms=[[c["atom1"], c["atom2"]] for c in cut]) + (
        f"<p class='note'>{len(rows) - limit} more in the JSON</p>" if len(rows) > limit else "")


def chem_table(rows, limit=30):
    if not rows:
        return "<p class='note'>none</p>"
    alt = any(c.get("altloc") for c in rows[:limit])
    return table(["atom", "atom", "distance (Å)", "problem", "verdict"] + (["conformer"] if alt else []),
                 [[c["atom1"], c["atom2"], c["distance"], c["reason"].replace("possible_protonation", "possible protonation"), c["verdict"]] + ([c.get("altloc", "")] if alt else [])
                  for c in rows[:limit]], cls={4: lambda v: v}, atoms=[[c["atom1"], c["atom2"]] for c in rows[:limit]]) + (f"<p class='note'>{len(rows) - limit} more in the JSON</p>" if len(rows) > limit else "")


def collect(path, a, b, summary_path=None, full_path=None, model_index=0, hydrogens=False, ensemble=False):
    pair = analyze(path, model_index=model_index, pair=(a, b), n_points=480, hydrogens=hydrogens)
    whole = analyze(path, model_index=model_index, n_points=480, hydrogens=hydrogens)
    if ensemble:
        from .ensemble import residue_rmsd
        rows = [r for r in residue_rmsd(path) if r["chain"] in (a, b)]
    else:
        rows = af3.residue_plddt(path, model_index)
    data = {"ensemble": bool(ensemble), "path": path, "model_index": model_index, "a": a, "b": b, "pair": pair, "whole": whole, "residues": rows,
            "summary": af3.load_summary(summary_path) if summary_path else None,
            "full": af3.load_full(full_path) if full_path else None,
            "pdockq": af3.pdockq(path, a, b, model_index)}
    it = pair["interfaces"][0] if pair["interfaces"] else None
    data["iface"] = it
    if it and data["full"] is not None:
        ra = {int(x["residue"].split(":")[1]) for x in it["residues_a"]}
        rb = {int(x["residue"].split(":")[1]) for x in it["residues_b"]}
        data["ipsae"] = {"a_to_b": af3.ipsae(data["full"], a, b), "b_to_a": af3.ipsae(data["full"], b, a)}
        data["ipsae"]["max"] = max(data["ipsae"]["a_to_b"] or 0, data["ipsae"]["b_to_a"] or 0)
        data["ipae"] = af3.interface_pae(data["full"], a, b, ra, rb)
    return data


def nav_html(molstar=False, exp=False):
    three = "<a href='#' id='ms-open' class='active'>3D view</a>" if molstar else ""
    if exp:
        return ("<nav class='top' aria-label='Report sections'>"
                "<a href='#s2' data-s='s2'>1 Interactions</a><a class='sub' href='#polar'>H-bonds</a>"
                "<a class='sub' href='#iclash'>clashes</a><a class='sub' href='#iissue'>issues</a>"
                "<a href='#s3' data-s='s3'>2 Model quality</a><a class='sub' href='#plddt'>B-factors</a>"
                "<a class='sub' href='#wclash'>clashes</a>" + three + "<span class='sp'></span>"
                "<button class='theme' onclick=\"var r=document.documentElement;r.dataset.theme=(r.dataset.theme==='dark'?'light':'dark')\">light / dark</button></nav>")
    return ("<nav class='top' aria-label='Report sections'>"
            "<a href='#s1' data-s='s1'>1 AI metrics</a><a class='sub' href='#pae'>PAE</a>"
            "<a href='#s2' data-s='s2'>2 Interactions</a><a class='sub' href='#polar'>H-bonds</a>"
            "<a class='sub' href='#iclash'>clashes</a><a class='sub' href='#iissue'>issues</a>"
            "<a href='#s3' data-s='s3'>3 Prediction quality</a><a class='sub' href='#plddt'>pLDDT</a>"
            "<a class='sub' href='#wclash'>clashes</a>" + three + "<span class='sp'></span>"
            "<button class='theme' onclick=\"var r=document.documentElement;r.dataset.theme=(r.dataset.theme==='dark'?'light':'dark')\">light / dark</button></nav>")


def chips_html(d):
    it, s = d["iface"], d["summary"]
    chips, attn = [], []
    if d.get("ipsae"):
        v = d["ipsae"]["max"]
        chips.append(("#s1", fmt(v), "ipSAE", st_hi(v, .6, .3)))
        if st_hi(v, .6, .3) != "good":
            attn.append((st_hi(v, .6, .3), "#s1", f"ipSAE {fmt(v)} (good at 0.6 or more)"))
    elif s and s.get("iptm") is not None:
        chips.append(("#s1", fmt(s["iptm"]), "ipTM", st_hi(s["iptm"], .8, .6)))
    if d.get("ipae"):
        v = d["ipae"]["mean"]
        chips.append(("#pae", fmt(v), "interface PAE", st_lo(v, 5, 10)))
        if st_lo(v, 5, 10) != "good":
            attn.append((st_lo(v, 5, 10), "#pae", f"interface PAE {fmt(v)} A (good at 5 or less)"))
    if it:
        chips.append(("#s2", f"{it['area']:.0f}", "A\u00b2 interface", ""))
        pol = it["contacts"].get("hbond", 0) + it["contacts"].get("salt_bridge", 0)
        chips.append(("#polar", str(pol), "H-bonds + salt bridges (direct)", ""))
        hbw = it["hbonds"]
        chips.append(("#mediated", str(hbw["water_mediated"]), "H-bonds via water" + (f" ({hbw['water_mediated_doubtful']} doubtful)" if hbw["water_mediated_doubtful"] else ""),
                      "warn" if hbw["water_mediated_doubtful"] else ""))
        cl, ch = split_problems(d["pair"]["contacts"], units=(d["a"], d["b"]))
        for h, lst, t, what in (("#iclash", cl, "interface clashes", "interface clash"), ("#iissue", ch, "interface issues", "interface issue")):
            st = st_probs(lst)
            chips.append((h, str(len(lst)), t, st))
            if st != "good":
                nb = sum(1 for c in lst if c["verdict"] == "unfavorable")
                attn.append((st, h, f"{len(lst)} {what}{'s' if len(lst) != 1 else ''}" + (f", {nb} unfavorable" if nb else ", all borderline")))
    wc, wi = split_problems(d["whole"]["contacts"])
    st = st_probs(wc)
    chips.append(("#wclash", str(len(wc)), "clashes, whole model", st))
    if st != "good":
        nb = sum(1 for c in wc if c["verdict"] == "unfavorable")
        attn.append((st, "#wclash", f"{len(wc)} clashes in the whole model" + (f", {nb} unfavorable" if nb else ", all borderline")))
    html = "<div class='chips'>" + "".join(f"<a class='st-{st}' href='{h}'><b>{esc(v)}</b><span>{esc(t)}</span></a>" for h, v, t, st in chips) + "</div>"
    worst = "bad" if any(x[0] == "bad" for x in attn) else "warn" if attn else "good"
    if attn:
        box = "<b>Needs attention</b><ul>" + "".join(f"<li><a href='{h}'>{esc(t)}</a></li>" for st, h, t in sorted(attn, key=lambda x: x[0] != "bad")) + "</ul>"
    else:
        box = "<b>Nothing to flag</b> in the checks below."
    legend = ("<div class='note'>Colour key: <span class='good'>green</span> acceptable, <span class='borderline'>amber</span> look at it, "
              "<span class='unfavorable'>red</span> likely a problem. The cut-offs are rules of thumb, not validated limits.</div>")
    return html + f"<div class='attn {worst}'>{box}</div>" + legend


NAV_JS = ("<script>(function(){document.querySelectorAll('nav.top a.sub, .chips a').forEach(function(a){if(!document.getElementById(a.getAttribute('href').slice(1)))a.remove()});var links=document.querySelectorAll('nav.top a[data-s]');var secs=[].map.call(links,function(a){return document.getElementById(a.dataset.s)});"
          "function upd(){var y=window.scrollY+90,cur=0;secs.forEach(function(e,i){if(e&&e.offsetTop<=y)cur=i});links.forEach(function(a,i){a.classList.toggle('active',i===cur)})}"
          "window.addEventListener('scroll',upd,{passive:true});upd();})();</script>")


MOLSTAR_VER = "5.0.0"
MOLSTAR_CDN = "https://cdn.jsdelivr.net/npm/molstar@" + MOLSTAR_VER + "/build/viewer/"

MS_CSS = """
.mvpanel{position:fixed;right:16px;bottom:16px;width:520px;height:420px;max-width:calc(100vw - 32px);z-index:30;background:var(--card);border:1px solid var(--line);border-radius:10px;box-shadow:0 8px 28px rgba(0,0,0,.28);display:none;flex-direction:column;overflow:hidden}
.mvpanel.on{display:flex}.mvpanel.big{width:min(980px,calc(100vw - 32px));height:72vh}
.mvbar{display:flex;gap:6px;align-items:center;padding:6px 8px;border-bottom:1px solid var(--line);flex-wrap:wrap}
.mvbar b{font-size:13px}.mvbar button{font-size:12px;padding:2px 9px;border:1px solid var(--line);background:var(--card);color:var(--fg);border-radius:6px;cursor:pointer}
.mvbar button.on{background:var(--acc);color:#fff;border-color:var(--acc)}.mvbar .sp{flex:1}
#mvc .msp-viewport-top-left-controls,#mvc .msp-snapshot-description-wrapper{display:none!important}#mvt{font-size:12px;color:var(--mut);padding:2px 10px;min-height:16px}#mvc{position:relative;flex:1;min-height:0}
table.ms tbody tr{cursor:pointer}table.ms tbody tr:hover{background:var(--chip)}table.ms tbody tr.sel{outline:2px solid var(--acc);outline-offset:-2px}
"""

MS_JS = r"""
(function(){
var CIF=__CIF__,FMT=__FMT__,MI=__MI__,CH=__CH__,RES=__RES__,IF=__IF__,CDN=__CDN__,NM=__NM__;
var BAND=['#0053d6','#65cbf3','#ffdb13','#ff7d45'];
var COL={};var PAL=[['#a9c4de','#e07b24'],['#b8d9b5','#8e44ad'],['#d9d3a8','#2a9d8f'],['#d8b4b4','#c0392b']];
CH.forEach(function(c,i){COL[c]=PAL[Math.min(i,PAL.length-1)]});
var V=null,lib=null,url=null,last=null,plddt=false,panel=document.getElementById('mvp'),msg=document.getElementById('mvt');
function say(t){msg.textContent=t||''}
function load(){if(lib)return lib;lib=new Promise(function(ok,no){var l=document.createElement('link');l.rel='stylesheet';l.href=CDN+'molstar.css';document.head.appendChild(l);
 var s=document.createElement('script');s.src=CDN+'molstar.js';s.onload=ok;s.onerror=function(){no(new Error('Mol* could not be downloaded'))};document.head.appendChild(s)})
 .then(function(){return molstar.Viewer.create('mvc',{layoutIsExpanded:false,layoutShowControls:false,layoutShowRemoteState:false,layoutShowSequence:false,layoutShowLog:false,layoutShowLeftPanel:false,viewportShowExpand:false,viewportShowSelectionMode:false,viewportShowAnimation:false})})
 .then(function(v){V=v;url=URL.createObjectURL(new Blob([CIF]));return v});return lib}
function P(a){var q=a.split(':');return {c:q[0],n:parseInt(q[1],10),t:q[2]||'',at:q[3]||''}}
function sel(r){return {auth_asym_id:r.c,auth_seq_id:r.n}}
function open_(){panel.classList.add('on')}
async function draw(o){
 open_();say('Loading Mol* ...');
 try{await load()}catch(e){say('Mol* could not be loaded. The 3D view needs an internet connection.');lib=null;return}
 var M=molstar.PluginExtensions.mvs,b=M.MVSData.createBuilder();
 var pa=b.download({url:url}).parse({format:FMT});
 if(o.all){for(var k=0;k<NM;k++){var s2=pa.modelStructure({model_index:k});CH.forEach(function(c){s2.component({selector:{auth_asym_id:c}}).representation({type:'cartoon'}).color({color:COL[c][0]}).opacity({opacity:.45})})}
  await M.loadMVSData(V.plugin,b.getState(),'mvsj',{keepCamera:false});say('All '+NM+' models');return}
 var st=pa.modelStructure({model_index:MI});
 if(plddt){var g=[[],[],[],[]];RES.forEach(function(r){var k=r[2]>=90?0:r[2]>=70?1:r[2]>=50?2:3;g[k].push({auth_asym_id:r[0],auth_seq_id:r[1]})});
  g.forEach(function(l,k){if(l.length)st.component({selector:l}).representation({type:'cartoon'}).color({color:BAND[k]})})}
 else CH.forEach(function(c){st.component({selector:{auth_asym_id:c}}).representation({type:'cartoon'}).color({color:COL[c][0]})});
 st.component({selector:'ligand'}).representation({type:'ball_and_stick'});
 var rs=[],seen={};(o.res||[]).forEach(function(r){var k=r.c+':'+r.n;if(!seen[k]){seen[k]=1;rs.push(r)}});
 if(rs.length){
  rs.forEach(function(r){var c=st.component({selector:sel(r)});c.representation({type:'ball_and_stick'}).color({color:(COL[r.c]||PAL[3])[1]});
   if(rs.length<=14)c.label({text:r.t+r.n+' ('+r.c+')'})});
  st.component({selector:rs.map(sel)}).focus({radius_factor:rs.length>2?1.5:1,radius_extent:rs.length>2?2:6})}
 if(o.pair){var A=o.pair[0],B=o.pair[1];
  st.primitives({}).distance({start:{auth_asym_id:A.c,auth_seq_id:A.n,auth_atom_id:A.at},end:{auth_asym_id:B.c,auth_seq_id:B.n,auth_atom_id:B.at},radius:0.07,color:'#e03c31',dash_length:0.12,label_size:0.9})}
 await M.loadMVSData(V.plugin,b.getState(),'mvsj',{keepCamera:false});say(o.title||'')}
function show(o){var res=[],pair=null,title='';if(o.all){last={all:1};draw(last);return}
 if(o.a){var p=o.a.map(P);res=p;title=p.map(function(x){return x.c+':'+x.t+x.n+(x.at?' '+x.at:'')}).join('  to  ');
  if(p.length==2&&p[0].at&&p[1].at)pair=[p[0],p[1]]}
 if(o.iface){res=IF.map(function(r){return {c:r[0],n:r[1],t:r[2]}});title='Interface residues ('+res.length+')'}
 last={res:res,pair:pair,title:title};draw(last)}
window.IFK={show:show,setModel:function(i,f,r){MI=i;if(f)IF=f;if(r)RES=r;if(panel.classList.contains('on'))draw({res:[],title:'Model '+(i+1)})}};
document.addEventListener('click',function(e){var tr=e.target.closest&&e.target.closest('table.ms tr[data-a]');if(!tr)return;
 document.querySelectorAll('table.ms tr.sel').forEach(function(x){x.classList.remove('sel')});tr.classList.add('sel');show({a:tr.dataset.a.split('|')})});
var mi=document.getElementById('mv-if');if(mi)mi.onclick=function(){show({iface:1})};var ma=document.getElementById('mv-all');if(ma)ma.onclick=function(){show({all:1})};
var pl=document.getElementById('mv-pl');if(pl)pl.onclick=function(){plddt=!plddt;this.classList.toggle('on',plddt);draw(last||{res:[]})};
document.getElementById('mv-big').onclick=function(){panel.classList.toggle('big');if(V)setTimeout(function(){V.plugin.canvas3d&&V.plugin.canvas3d.handleResize()},50)};
document.getElementById('mv-x').onclick=function(){panel.classList.remove('on')};
draw({res:[]});var op=document.getElementById('ms-open');if(op)op.onclick=function(e){e.preventDefault();if(panel.classList.contains('on'))panel.classList.remove('on');else draw(last||{res:[]})};
})();
"""


def molstar_block(d):
    """Floating Mol* panel plus the code that drives it. The structure is embedded; Mol* itself comes from a CDN."""
    path = d["path"]
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    fmt = "pdb" if path.lower().endswith((".pdb", ".ent")) else "mmcif"
    chains = []
    for r in d["residues"]:
        if r["chain"] not in chains:
            chains.append(r["chain"])
    # the two chains of the pair first so they get the first two colours
    first = [c for c in (d.get("a"), d.get("b")) if c in chains]
    chains = first + [c for c in chains if c not in first]
    res = [[r["chain"], r["resnum"], r["b"]] for r in d["residues"]]
    it = d.get("iface")
    iface = []
    if it:
        for key in ("residues_a", "residues_b"):
            for r in it[key]:
                ch, num, name = r["residue"].split(":")[:3]
                iface.append([ch, int(num), name])

    def js(x):
        return json.dumps(x).replace("</", "<\\/")
    code = (MS_JS.replace("__CIF__", js(text)).replace("__FMT__", js(fmt)).replace("__MI__", str(int(d.get("model_index", 0))))
            .replace("__NM__", str(int(d.get("n_models", 1)))).replace("__CH__", js(chains)).replace("__RES__", js(res)).replace("__IF__", js(iface)).replace("__CDN__", js(d.get("molstar_url") or MOLSTAR_CDN)))
    return (f"<style>{MS_CSS}</style><div class='mvpanel' id='mvp'><div class='mvbar'><b>3D view (Mol*)</b>"
            + ("<button type='button' id='mv-if'>interface</button>" if it else "") + ("<button type='button' id='mv-all'>all models</button>" if d.get("n_models", 1) > 1 else "") + ("" if d.get("kind") == "experimental" else "<button type='button' id='mv-pl'>pLDDT colours</button>") + "<span class='sp'></span>"
            "<button type='button' id='mv-big'>larger</button><button type='button' id='mv-x'>close</button></div>"
            "<div id='mvt'></div><div id='mvc'></div></div><script>" + code + "</script>")


def render(d, molstar=False):
    _MS[0] = bool(molstar)
    try:
        return _render(d, molstar)
    finally:
        _MS[0] = False


def _render(d, molstar):
    exp = d.get("kind") == "experimental"
    a, b = d["a"], d["b"]
    it, s, pair, whole = d["iface"], d["summary"], d["pair"], d["whole"]
    name = d["path"].replace("\\", "/").split("/")[-1]
    out = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
           f"<title>{esc(name)} interface report</title><style>{CSS}</style></head><body><main>",
           banner_html() + f"<h1>{esc(name)}" + (" <span class='tag'>3D view</span>" if molstar else "") + "</h1>",
           f"<p class='sub'>Chain {esc(a)} with chain {esc(b)}. Report made by ifacekit {esc(pair['tool']['version'])} "
           f"on {esc(pair['input']['format'])}. Heavy-atom analysis, no hydrogens; how it compares with MolProbity is stated at the end of the report.</p>"]
    units = {u["id"]: u for u in pair["units"]}
    chain_line = ", ".join(f"chain {u}: {units[u]['n_residues']} residues" for u in (a, b) if u in units)
    out.append(f"<p class='note'>{esc(chain_line)}</p>")
    out.append(BODY_START)
    out.append(nav_html(molstar, exp).replace(">B-factors</a>", ">RMSD</a>") if d.get("ensemble") else nav_html(molstar, exp))
    out.append(chips_html(d))
    if molstar:
        out.append("<p class='note'>3D view: click a row in any table, or a residue on the plot, to show it in Mol*.</p>")

    # ---------- 1 AI metrics
    if not exp:
        out.append('<section id="s1" class="sec s1"><h2>1. AI confidence metrics</h2>')
    cards = []
    if s and not exp:
        cp = s.get("chain_pair_iptm")
        ia, ib = (s["chains"].index(a), s["chains"].index(b)) if a in s["chains"] and b in s["chains"] else (None, None)
        pair_iptm = cp[ia][ib] if cp and ia is not None else None
        pae_min = (s.get("chain_pair_pae_min") or [[None]])[ia][ib] if ia is not None and s.get("chain_pair_pae_min") else None
        cards += [card("ipTM", fmt(s["iptm"]), "whole model", st_hi(s["iptm"], .8, .6)), card("pTM", fmt(s["ptm"]), "whole model", st_hi(s["ptm"], .8, .5)),
                  card("ranking score", fmt(s["ranking_score"]), "AF3"), card(f"ipTM {a}-{b}", fmt(pair_iptm), "chain pair", st_hi(pair_iptm, .8, .6)),
                  card(f"min PAE {a}-{b}", fmt(pae_min), "Angstrom, chain pair")]
        if s.get("has_clash") is not None:
            cards.append(card("AF3 clash flag", "yes" if s["has_clash"] else "no", f"disordered {fmt(s['fraction_disordered'])}", "bad" if s["has_clash"] else "good"))
    if d.get("ipsae") and not exp:
        ip = d["ipsae"]
        cards += [card("ipSAE", fmt(ip["max"]), f"{a}>{b} {fmt(ip['a_to_b'])}, {b}>{a} {fmt(ip['b_to_a'])}", st_hi(ip["max"], .6, .3)),
                  card("interface PAE", fmt(d["ipae"]["mean"]) if d["ipae"] else "-", "mean over interface residues", st_lo(d["ipae"]["mean"] if d["ipae"] else None, 5, 10)),
                  card("PAE < 5 at interface", f"{d['ipae']['fraction_below_5'] * 100:.0f}%" if d["ipae"] else "-", "fraction of residue pairs")]
    if d.get("pdockq") and not exp:
        cards.append(card("pDockQ", fmt(d["pdockq"]["pdockq"], 3), f"{d['pdockq']['n_contacts']} CB contacts"))
    if it and not exp:
        cards.append(card("interface pLDDT", fmt(round((it["mean_b_a"] + it["mean_b_b"]) / 2, 1), 1) if it.get("mean_b_a") is not None else "-",
                          f"{a} {fmt(it['mean_b_a'], 1)}, {b} {fmt(it['mean_b_b'], 1)}",
                          st_hi((it["mean_b_a"] + it["mean_b_b"]) / 2 if it.get("mean_b_a") is not None else None, 90, 70)))
    if not exp:
        out.append('<div class="cards">' + "".join(cards) + "</div>" if cards else "<p class='note'>No AF3 confidence files given, only the pLDDT stored in the model is shown.</p>")
    if not d.get("full") and s and not exp:
        out.append("<p class='note'>No full-data file was given, so ipSAE and interface PAE are not computed.</p>")
    if d.get("full") is not None and not exp:
        out.append("<h3 id='pae'>Predicted aligned error</h3><div class='row'><div>" + pae_canvas(d["full"], [a, b]) + "</div>"
                   "<div class='note'>Rows are the aligned residues, columns the scored residues. Dark blue is low error, white is high. "
                   "Red lines mark chain boundaries. The off-diagonal blocks show how well the two chains are placed relative to each other. "
                   "ipSAE here uses a PAE cutoff of 10 and is a re-implementation of the published definition. pDockQ is from CB contacts within 8 A and pLDDT.</div></div>")

    # ---------- 2 interactions
    out.append(("" if exp else "</section>") + f'<section id="s2" class="sec s2"><h2>{1 if exp else 2}. Interface interactions</h2>')
    if not it:
        out.append(f"<p class='note'>No interface of at least the minimum area between {esc(a)} and {esc(b)}.</p>")
    else:
        c = it["contacts"]
        out.append('<div class="cards">' + "".join([
            card("interface area", f"{it['area']:.0f} A²", f"{a} {it['area_a']:.0f}, {b} {it['area_b']:.0f}"),
            card("polar fraction", f"{(it['polar_fraction'] or 0) * 100:.0f}%", ""),
            card("interface residues", f"{it['n_residues_a']} + {it['n_residues_b']}", f"{a} + {b}"),
            card("hydrogen bonds", it["hbonds"]["direct"], "direct, no water"),
            card("H-bonds via water", it["hbonds"]["water_mediated"],
                 (f"waters bonded to both chains: {it['hbonds']['water_mediated_clear']} clear, {it['hbonds']['water_mediated_doubtful']} doubtful"
                  if it["hbonds"]["water_mediated"] else "waters bonded to both chains"),
                 "warn" if it["hbonds"]["water_mediated_doubtful"] else ""),
            card("H-bond connections", it["hbonds"]["total_connections"], f"direct + via water; clear only: {it['hbonds']['total_connections_clear']}"),
            card("salt bridges", c.get("salt_bridge", 0), ""),
            card("hydrophobic", c.get("hydrophobic", 0), "carbon pairs under 4 A"),
            card("other vdW", c.get("vdw", 0), ""),
            card("ring / cation-pi", sum(1 for x in pair["pi_interactions"] if x["inter_unit"]), "stacking, T-shaped, cation-pi"),
            card("buried polar, no partner", len(it["buried_unsatisfied_polar"]), "atoms", st_lo(len(it["buried_unsatisfied_polar"]), 2, 6))]) + "</div>")
        pol = [x for x in pair["contacts"] if {x["unit1"], x["unit2"]} == {a, b} and x["type"] in ("hbond", "salt_bridge", "metal")]
        pol.sort(key=lambda x: x["distance"])
        out.append(f"<h3 id='polar'>Hydrogen bonds, salt bridges and metal contacts between {esc(a)} and {esc(b)} ({len(pol)})</h3>"
                   "<p class='note'>Contacts between the two chains themselves. Contacts with waters are not listed here: waters that bridge the two "
                   "chains are in the table of bonds via water below, and the other protein-water contacts belong to the model quality part.</p>")
        out.append(table(["atom A", "atom B", "distance (Å)", "type", "verdict"],
                         [[x["atom1"], x["atom2"], x["distance"], x["type"], x["verdict"]] for x in pol],
                         cls={4: lambda v: v}, atoms=[[x["atom1"], x["atom2"]] for x in pol]))
        out.append(_mediated_html(it, a, b, pair["flags"].get("n_water_residues", 0)))
        out.append("<h3 id='buried'>Most buried residues</h3><div class='row'>")
        for side, key in ((a, "residues_a"), (b, "residues_b")):
            top = sorted(it[key], key=lambda r: -r["bsa"])[:10]
            out.append("<div><div class='note'>chain " + esc(side) + "</div>" + table(
                ["residue", "buried A²", "contacts", "B-factor" if exp else "pLDDT"], [[r["residue"], r["bsa"], r["contacts"], r["mean_b"]] for r in top],
                atoms=[[r["residue"]] for r in top]) + "</div>")
        out.append("</div>")
        if it["buried_unsatisfied_polar"]:
            out.append("<h3>Buried polar atoms without a partner</h3>" + table(
                ["atom", "buried A²"], [[x["atom"], x["buried"]] for x in it["buried_unsatisfied_polar"]],
                atoms=[[x["atom"]] for x in it["buried_unsatisfied_polar"]]))
        clashes, chem = split_problems(pair["contacts"], units=(a, b))
        n_unf = sum(1 for c in clashes if c["verdict"] == "unfavorable")
        out.append(f"<h3 id='iclash'>Interface clashes ({len(clashes)}, of which {n_unf} unfavorable)</h3>"
                   "<p class='note'>Atoms of the two chains that overlap in space. Unfavorable means an overlap of 0.6 A or more, "
                   "borderline 0.25 to 0.6 A. Overlap is the sum of the van der Waals radii minus the distance, heavy atoms only.</p>"
                   + clash_table(clashes) + hclash_html(pair.get("hydrogens"), units=(a, b), anchor="ihclash"))
        out.append(f"<h3 id='iissue'>Interface issues that are not clashes ({len(chem)})</h3>"
                   "<p class='note'>Hydrogen bonds or salt bridges that are too short, polar atoms of the same kind facing each other, "
                   "and like charges close together. A carboxylate near another acceptor can be protonated, which hydrogen-free "
                   "analysis cannot see, so these are hints, not errors.</p>" + chem_table(chem))

    # ---------- 3 quality
    out.append('</section><section id="s3" class="sec s3"><h2>' + ("2. Quality of the model" if exp else "3. Quality of the prediction") + '</h2>')
    res = d["residues"]
    iface = set()
    if it:
        for key in ("residues_a", "residues_b"):
            for r in it[key]:
                ch, num, _ = r["residue"].split(":")
                iface.add((ch, int(num)))
    if d.get("ensemble"):
        out.append("<h3 id='plddt'>Backbone deviation along the sequence</h3>" + plddt_plot(res, [a, b], iface, True, label="deviation (A)")
                   + "<p class='legend'>Line: mean over the models of the distance of each residue from the mean structure (root of the mean squared distance of the backbone atoms N, CA, C, O to their mean positions, "
                   "each chain superposed on its least mobile 70% of residues). Dark band: mean \u00b1 1 standard deviation; light band: lowest to highest model; the RMSD is in the tooltip. "
                   "Red ticks and shading mark interface residues of the first model. Scroll over the plot to zoom, drag to pan, and hover to read a residue.</p>")
    elif exp:
        out.append("<h3 id='plddt'>B-factor along the sequence</h3>" + plddt_plot(res, [a, b], iface, True)
                   + "<p class='legend'>Mean B-factor of each residue. Red ticks and shading mark interface residues. Scroll over the plot to zoom, drag to pan, and hover to read a residue.</p>")
    else:
        out.append("<h3 id='plddt'>pLDDT along the sequence</h3>" + plddt_plot(res, [a, b], iface)
                   + "<p class='legend'>Dashed lines at 90, 70 and 50. Red ticks and shading mark interface residues. Scroll over the plot to zoom, drag to pan, and hover to read a residue. In a predicted model the B-factor column holds pLDDT.</p>")
    out.append("<div class='row'>")
    for ch in (a, b):
        bands = None if exp else af3.plddt_bands([r["b"] for r in res if r["chain"] == ch])
        if bands:
            out.append(f"<div><div class='note'>chain {esc(ch)}: mean {bands['mean']}, residues {bands['n']}</div>{band_bar(bands)}"
                       f"<div class='legend'>&ge;90 {bands['very_high_ge90'] * 100:.0f}%, 70-90 {bands['confident_70_90'] * 100:.0f}%, "
                       f"50-70 {bands['low_50_70'] * 100:.0f}%, &lt;50 {bands['very_low_lt50'] * 100:.0f}%</div></div>")
    out.append("</div>")
    sm = whole["summary"]
    v = sm["verdicts"]
    out.append("<h3>Atom contacts in the whole model</h3><div class='cards'>" + "".join([
        card("unfavorable pairs", v.get("unfavorable", 0), "steric overlap, short H-bond, mismatch", "bad" if v.get("unfavorable", 0) else "good"),
        card("borderline pairs", v.get("borderline", 0), "", "warn" if v.get("borderline", 0) else "good"), card("favorable pairs", v.get("favorable", 0), ""),
        card("unfavorable per 1000 atoms", sm["unfavorable_per_1000_atoms"], "not a MolProbity clashscore")]) + "</div>")
    clashes, chem = split_problems(whole["contacts"])
    w_unf = sum(1 for c in clashes if c["verdict"] == "unfavorable")
    out.append(f"<h3 id='wclash'>Clashes in the whole model ({len(clashes)}, of which {w_unf} unfavorable)</h3>" + clash_table(clashes, 20) + hclash_html(whole.get("hydrogens"), limit=20, anchor="whclash"))
    out.append(f"<h3 id='wissue'>Other contact issues in the whole model ({len(chem)})</h3>" + chem_table(chem, 20))
    g = whole.get("geometry") or {}
    if g.get("n_bonds"):
        out.append(f"<h3>Backbone geometry</h3><p class='note'>bond rmsZ {g['rmsz_bonds']}, angle rmsZ {g['rmsz_angles']}, "
                   f"outliers (|Z| over 4): {g['n_outliers']}. Approximate Engh and Huber values.</p>")
        if g["outliers"]:
            out.append(table(["residue", "kind", "atoms", "value", "ideal", "Z"],
                             [[o["residue"], o["kind"], o["atoms"], o["value"], o["ideal"], o["z"]] for o in g["outliers"][:12]],
                             atoms=[[o["residue"]] for o in g["outliers"][:12]]))
    try:
        from .ramaplot import fragment as _rama_fragment
        out.append(_rama_fragment(d["path"], [a, b], d["model_index"], pdf_name=d.get("pdf")))
    except Exception as e:  # the extra section must never break the report
        out.append(f"<p class='note'>Ramachandran analysis failed: {esc(e)}</p>")
    try:
        from .rotamerhtml import fragment as _rot_fragment
        out.append(_rot_fragment(d["path"], [a, b], d["model_index"]))
    except Exception as e:
        out.append(f"<p class='note'>Rotamer analysis failed: {esc(e)}</p>")
    out.append("<h3>What this report does not cover</h3><p class='note'>Hydrogen-based clashes and Asn, Gln and His flips are not evaluated. "
               "Use MolProbity or the wwPDB validation report for those. " + comparison_html() + "</p>")
    out.append("</section>")
    out.append(NAV_JS)
    out.append(BODY_END)
    from .about import footer_html
    out.append(footer_html())
    if molstar:
        out.append(molstar_block(d))
    out.append("</main></body></html>")
    return "".join(out)


def build(path, a, b, out_path, summary_path=None, full_path=None, model_index=0, molstar_path=None, kind="predicted", pdf_name=None, molstar_url=None, hydrogens=False, ensemble=False):
    d = collect(path, a, b, summary_path, full_path, model_index, hydrogens, ensemble)
    d["kind"] = kind
    d["molstar_url"] = molstar_url
    d["pdf"] = pdf_name
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(render(d))
    if molstar_path:
        with open(molstar_path, "w", encoding="utf-8") as fh:
            fh.write(render(d, molstar=True))
    return d


def find_confidences(cif_path):
    """Guess the AF3 summary and full-data files next to a model by matching its name."""
    import glob
    import os
    d = os.path.dirname(os.path.abspath(cif_path))
    stem = os.path.splitext(os.path.basename(cif_path))[0].lower().replace("-", "_")
    out = {}
    for kind, pat in (("summary", "*summary_confidences*.json"), ("full", "*full_data*.json")):
        hits = [f for f in glob.glob(os.path.join(d, pat)) if stem in os.path.basename(f).lower().replace("-", "_")]
        if len(hits) == 1:
            out[kind] = hits[0]
    return out
