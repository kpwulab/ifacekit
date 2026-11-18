"""HTML report for one ligand bound to a protein: interaction diagram, tables and contact quality."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from . import hbonds
from . import ligdetail
from .banner import banner_html
from . import htmlreport as H
from . import ligplot
from . import model as mdl
from .about import comparison_short_html
from .chem import WATER_NAMES
from .core import analyze
from .sasa import atom_sasa


DG_CSS = (".dgbar{display:flex;flex-wrap:wrap;gap:12px;align-items:center;justify-content:space-between;margin:6px 0}"
          ".seg{display:inline-flex;border:1px solid var(--line,#c9d1d6);border-radius:8px;overflow:hidden}"
          ".seg button,.dl button{font:inherit;font-size:13px;padding:5px 12px;border:0;background:var(--card,#fff);color:inherit;cursor:pointer}"
          ".seg button.on{background:var(--acc,#2b6cb0);color:#fff}"
          ".dl{font-size:13px;color:var(--mut,#667078)}.dl button{border:1px solid var(--line,#c9d1d6);border-radius:7px;margin-left:6px}"
          ".dl button:hover,.seg button:hover:not(.on){background:rgba(120,140,160,.15)}[hidden]{display:none!important}")

DG_JS = """<script>
(function(){
var cur='modern',det='';function idOf(){return 'dg-'+cur+(det==='all'?'-all':'')}
function show(){['modern','classic','detailed'].forEach(function(k){['','-all'].forEach(function(x){var e=document.getElementById('dg-'+k+x);if(e)e.hidden=!(k===cur&&((x==='-all')===(det==='all')))});var b=document.getElementById('v-'+k);if(b)b.classList.toggle('on',k===cur)});document.getElementById('d-simple').classList.toggle('on',det!=='all');document.getElementById('d-all').classList.toggle('on',det==='all')}
window.dgView=function(v){cur=v;show()};window.dgDetail=function(x){det=(x==='all')?'all':'';show()};
function clean(){
  var svg=document.querySelector('#'+idOf()+' svg'), vb=svg.getAttribute('viewBox').split(' ').map(Number), w=vb[2], h=vb[3];
  var t=new XMLSerializer().serializeToString(svg);
  t=t.replace(/var\\(--[a-z-]+,\\s*([^)]+)\\)/g,'$1').replace(/currentColor/g,'#1d2326');
  t=t.replace(/<svg([^>]*?)width="100%"/,'<svg$1width="'+w+'"').replace(/<svg /,'<svg height="'+h+'" ');
  t=t.replace(/(<svg[^>]*>)/,'$1<rect x="0" y="0" width="'+w+'" height="'+h+'" fill="#ffffff"/>');
  return {t:t,w:w,h:h};
}
function save(blob,ext){var a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=document.getElementById('dgbox').dataset.name+'_'+cur+(det==='all'?'_all':'')+'.'+ext;document.body.appendChild(a);a.click();setTimeout(function(){URL.revokeObjectURL(a.href);a.remove()},500)}
function raster(cb,type,q){
  var c=clean(), k=300/96, img=new Image();
  img.onload=function(){var cv=document.createElement('canvas');cv.width=Math.round(c.w*k);cv.height=Math.round(c.h*k);var g=cv.getContext('2d');g.fillStyle='#fff';g.fillRect(0,0,cv.width,cv.height);g.drawImage(img,0,0,cv.width,cv.height);cv.toBlob(function(b){cb(b,cv.width,cv.height,c)},type,q)};
  img.src='data:image/svg+xml;charset=utf-8,'+encodeURIComponent(c.t);
}
var T=null;function crc(b){if(!T){T=[];for(var n=0;n<256;n++){var c=n;for(var k=0;k<8;k++)c=c&1?0xEDB88320^(c>>>1):c>>>1;T[n]=c>>>0}}var r=0xFFFFFFFF;for(var i=0;i<b.length;i++)r=T[(r^b[i])&255]^(r>>>8);return(r^0xFFFFFFFF)>>>0}
function dpi300(buf){var u=new Uint8Array(buf),d=new Uint8Array(21),v=new DataView(d.buffer),ppm=11811;v.setUint32(0,9);d.set([112,72,89,115],4);v.setUint32(8,ppm);v.setUint32(12,ppm);d[16]=1;v.setUint32(17,crc(d.subarray(4,17)));
  var o=new Uint8Array(u.length+21);o.set(u.subarray(0,33));o.set(d,33);o.set(u.subarray(33),54);return o}
function pdf(jpg,pw,ph,w,h){
  var enc=new TextEncoder(),parts=[],off=[],len=0;function add(x){var b=typeof x==='string'?enc.encode(x):x;parts.push(b);len+=b.length}
  function obj(n,s){off[n]=len;add(n+' 0 obj\\n'+s+'\\nendobj\\n')}
  var W=(w*0.75).toFixed(2),Hh=(h*0.75).toFixed(2);add('%PDF-1.4\\n');
  obj(1,'<< /Type /Catalog /Pages 2 0 R >>');obj(2,'<< /Type /Pages /Kids [3 0 R] /Count 1 >>');
  obj(3,'<< /Type /Page /Parent 2 0 R /MediaBox [0 0 '+W+' '+Hh+'] /Resources << /XObject << /Im0 5 0 R >> >> /Contents 4 0 R >>');
  var cs='q '+W+' 0 0 '+Hh+' 0 0 cm /Im0 Do Q';obj(4,'<< /Length '+cs.length+' >>\\nstream\\n'+cs+'\\nendstream');
  off[5]=len;add('5 0 obj\\n<< /Type /XObject /Subtype /Image /Width '+pw+' /Height '+ph+' /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length '+jpg.length+' >>\\nstream\\n');add(jpg);add('\\nendstream\\nendobj\\n');
  var x=len;add('xref\\n0 6\\n0000000000 65535 f \\n');for(var i=1;i<=5;i++)add(('0000000000'+off[i]).slice(-10)+' 00000 n \\n');
  add('trailer\\n<< /Size 6 /Root 1 0 R >>\\nstartxref\\n'+x+'\\n%%EOF');
  var out=new Uint8Array(len),p=0;parts.forEach(function(b){out.set(b,p);p+=b.length});return new Blob([out],{type:'application/pdf'})}
window.dgSave=function(kind){
  if(kind==='svg'){save(new Blob([clean().t],{type:'image/svg+xml'}),'svg');return}
  if(kind==='png'){raster(function(b){b.arrayBuffer().then(function(ab){save(new Blob([dpi300(ab)],{type:'image/png'}),'png')})},'image/png');return}
  raster(function(b,pw,ph,c){b.arrayBuffer().then(function(ab){save(pdf(new Uint8Array(ab),pw,ph,c.w,c.h),'pdf')})},'image/jpeg',0.95)
};
})();
</script>"""


def collect(path, ligand_unit=None, model_index=0, predicted=None, detailed=False, detailed_all=False):
    r = analyze(path, model_index=model_index, n_points=480)
    m = mdl.load(path, model_index)
    ligs = [u for u in m.units if u["kind"] == "ligand"]
    if not ligs:
        raise ValueError("no ligand in this file")
    lu = ligand_unit or ligs[0]["id"]
    if lu not in {u["id"] for u in ligs}:
        raise ValueError(f"ligand '{lu}' not found. Ligand units in this file (chain:name+residue number): "
                         + ", ".join(u["id"] for u in ligs))
    uid = next(k for k, u in enumerate(m.units) if u["id"] == lu)
    li = np.where(m.unit_idx == uid)[0]
    its = [i for i in r["interfaces"] if lu in (i["unit_a"], i["unit_b"]) and "polymer" in (i["kind_a"], i["kind_b"])]
    it = its[0] if its else None
    pu = None
    if it:
        pu = it["unit_a"] if it["unit_b"] == lu else it["unit_b"]
    free = atom_sasa(m.xyz[li], m.radius[li], n_points=480)
    comp = atom_sasa(m.xyz, m.radius, targets=li, n_points=480)
    contacts = [c for c in r["contacts"] if c["inter_unit"] and lu in (c["unit1"], c["unit2"])]
    if predicted is not None:                       # the user said what kind of model this is
        r["flags"]["predicted_model"] = bool(predicted)
        r["flags"]["model_source"] = "set by the user"
    d = ligplot.build(path, r, lu, model_index=model_index, detail="simple")
    d_all = ligplot.build(path, r, lu, model_index=model_index, detail="full")
    hb = hbonds.split(r["contacts"], lambda u: u == lu, lambda u: u not in (lu, "water"),
                      xyz=hbonds.xyz_lookup(m, lambda mm, k: f"{mm.residues[mm.res_idx[k]]['id']}:{mm.atom_name[k]}"))
    return {"path": path, "model_index": model_index, "result": r, "model": m, "lu": lu, "li": li, "iface": it, "pu": pu,
            "free": free, "comp": comp, "contacts": contacts, "diagram": d, "diagram_all": d_all, "hbonds": hb,
            "resname": m.residues[m.res_idx[li[0]]]["resname"], "a": pu or "A", "b": m.units[uid]["chain"],
            "residues": [], "kind": "ligand", "detailed": bool(detailed), "detailed_all": bool(detailed and detailed_all)}


def render(d, molstar=False):
    H._MS[0] = bool(molstar)
    try:
        return _render(d, molstar)
    finally:
        H._MS[0] = False


def _render(d, molstar):
    from . import af3
    esc, card, table, fmt = H.esc, H.card, H.table, H.fmt
    m, r, it, dg = d["model"], d["result"], d["iface"], d["diagram"]
    dga, hb = d["diagram_all"], d["hbonds"]
    name = d["path"].replace("\\", "/").split("/")[-1]
    li = d["li"]
    smi = Chem.MolToSmiles(dg["mol"])
    mol = dg["mol"]
    formula = rdMolDescriptors.CalcMolFormula(mol)
    free, comp = d["free"], d["comp"]
    buried = float(free.sum() - comp.sum())
    frac = buried / max(float(free.sum()), 1e-9)
    lp = [float(m.bfactor[k]) for k in li]
    pocket_res = {(n["ch"], n["num"]) for n in dga["nodes"]}
    pocket_b = [float(np.mean(m.bfactor[m.res_idx == m.res_idx[k]])) for k in range(m.n_atoms)
                if False]
    pb = []
    for n in dga["nodes"]:
        sel = [k for k in range(m.n_atoms) if m.residues[m.res_idx[k]]["chain"] == n["ch"] and str(m.residues[m.res_idx[k]]["resnum"]) == n["num"]]
        n["b"] = float(np.mean(m.bfactor[sel])) if sel else None
        if n["b"] is not None:
            pb.append(n["b"])
    pred = r["flags"].get("predicted_model")
    cs = d["contacts"]
    close = [c for c in cs if c["distance"] <= 4.0]
    nh = hb["direct"]
    ns = sum(1 for c in close if c["type"] == "salt_bridge")
    nm = sum(1 for c in close if c["type"] == "metal")
    nhy = sum(1 for c in close if c["type"] == "hydrophobic")
    prob = [c for c in cs if c["verdict"] in ("borderline", "unfavorable")]
    clashes, chem = H.split_problems(prob, only_between=True)
    unsat = it["buried_unsatisfied_polar"] if it else []
    out = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
           f"<title>{esc(name)} ligand report</title><style>{H.CSS}{DG_CSS}</style></head><body><main>",
           banner_html() + f"<h1>{esc(name)} · ligand {esc(d['resname'])}" + (" <span class='tag'>3D view</span>" if molstar else "") + "</h1>",
           f"<p class='sub'>Ligand {esc(d['resname'])} in chain {esc(d['b'])} against protein chain {esc(d['a'])}. Report made by ifacekit "
           f"{esc(r['tool']['version'])} on {esc(r['input']['format'])}. Heavy-atom analysis, no hydrogens. "
           f"Model type: {'prediction (confidence column read as pLDDT)' if r['flags'].get('predicted_model') else 'experimental (confidence column read as B-factor)'}"
           f", {esc(r['flags'].get('model_source', ''))}.</p>",
           f"<p class='note'>{len(li)} heavy atoms, formula without hydrogens: {esc(formula)}. "
           f"Structure used for the drawing: <code>{esc(smi)}</code>. Bonds and bond orders come from the " + ("chemical component dictionary written in the file, so they are exact." if dg.get("bond_source") == "file dictionary" else "3D geometry (the file has no chemistry), so they are a guess. Check the drawing against the real compound.") + "</p>"]
    nav = ("<nav class='top' aria-label='Report sections'><a href='#s1' data-s='s1'>1 Binding</a><a class='sub' href='#diagram'>diagram</a>"
           "<a href='#s2' data-s='s2'>2 Interactions</a><a class='sub' href='#polar'>H-bonds</a><a class='sub' href='#pocket'>pocket</a>"
           "<a href='#s3' data-s='s3'>3 Contact quality</a><a class='sub' href='#lclash'>clashes</a>"
           + ("<a href='#' id='ms-open' class='active'>3D view</a>" if molstar else "") + "<span class='sp'></span>"
           "<button class='theme' onclick=\"var r=document.documentElement;r.dataset.theme=(r.dataset.theme==='dark'?'light':'dark')\">light / dark</button></nav>")
    out.append(nav)
    # chips and attention
    attn = []
    if clashes:
        nb = sum(1 for c in clashes if c["verdict"] == "unfavorable")
        attn.append(("bad" if nb else "warn", "#lclash", f"{len(clashes)} ligand clash{'es' if len(clashes) != 1 else ''}" + (f", {nb} unfavorable" if nb else ", all borderline")))
    if chem:
        nb = sum(1 for c in chem if c["verdict"] == "unfavorable")
        attn.append(("bad" if nb else "warn", "#lissue", f"{len(chem)} other contact issue{'s' if len(chem) != 1 else ''}"))
    if unsat:
        attn.append(("warn", "#unsat", f"{len(unsat)} buried polar atom{'s' if len(unsat) != 1 else ''} without a partner"))
    lmean = float(np.mean(lp)) if lp else None
    if pred and lmean is not None and H.st_hi(lmean, 70, 50) != "good":
        attn.append((H.st_hi(lmean, 70, 50), "#s3", f"ligand pLDDT {lmean:.0f} (good at 70 or more)"))
    chips = [("#diagram", f"{buried:.0f}", "A² ligand buried", ""), ("#polar", str(nh + ns + nm), "H-bonds + salt bridges (direct)", ""), ("#mediated", str(hb["water_mediated"]), "H-bonds via water" + (f" ({hb['water_mediated_doubtful']} doubtful)" if hb["water_mediated_doubtful"] else ""), "warn" if hb["water_mediated_doubtful"] else ""),
             ("#pocket", str(len(dga["nodes"])), "pocket residues", ""),
             ("#lclash", str(len(clashes)), "ligand clashes", H.st_probs(clashes)),
             ("#lissue", str(len(chem)), "other issues", H.st_probs(chem))]
    out.append("<div class='chips'>" + "".join(f"<a class='st-{st}' href='{h}'><b>{esc(v)}</b><span>{esc(t)}</span></a>" for h, v, t, st in chips) + "</div>")
    worst = "bad" if any(a[0] == "bad" for a in attn) else "warn" if attn else "good"
    out.append(f"<div class='attn {worst}'>" + ("<b>Needs attention</b><ul>" + "".join(f"<li><a href='{h}'>{esc(t)}</a></li>" for s, h, t in sorted(attn, key=lambda x: x[0] != 'bad')) + "</ul>" if attn else "<b>Nothing to flag</b> in the checks below.") + "</div>")
    out.append("<div class='note'>Colour key: <span class='good'>green</span> acceptable, <span class='borderline'>amber</span> look at it, <span class='unfavorable'>red</span> likely a problem. Cut-offs are rules of thumb.</div>")
    if molstar:
        out.append("<p class='note'>3D view: click a row in any table to show it in Mol*.</p>")
    # 1 binding
    out.append('<section id="s1" class="sec s1"><h2>1. Ligand binding</h2><div class="cards">' + "".join([
        card("ligand buried", f"{frac * 100:.0f}%", f"{buried:.0f} of {free.sum():.0f} A² exposed when alone"),
        card("hydrogen bonds", nh, "ligand to protein, no water"), card("H-bonds to waters", hb["side_a_to_water"], "ligand to water"),
        card("H-bonds via water", hb["water_mediated"],
             (f"water bonded to ligand and protein: {hb['water_mediated_clear']} clear, {hb['water_mediated_doubtful']} doubtful" if hb["water_mediated"] else "water bonded to ligand and protein"),
             "warn" if hb["water_mediated_doubtful"] else ""), card("salt bridges", ns, ""), card("hydrophobic contacts", nhy, "carbon pairs under 4 A"),
        card("pocket residues", len(dga["nodes"]), "with a contact under 4 A"),
        card("buried polar, no partner", len(unsat), "ligand atoms", H.st_lo(len(unsat), 1, 3)),
        card("ligand pLDDT" if pred else "ligand B-factor", f"{lmean:.0f}", "mean over ligand atoms", H.st_hi(lmean, 70, 50) if pred else ""),
        card("pocket pLDDT" if pred else "pocket B-factor", f"{np.mean(pb):.0f}" if pb else "-", "mean over pocket residues", H.st_hi(float(np.mean(pb)), 90, 70) if pred and pb else "")]) + "</div>")
    out.append("<h3 id='diagram'>Interaction diagram</h3>")
    out.append("<div class='dgbar'><span class='seg' role='group' aria-label='Diagram style'>"
               "<button type='button' id='v-modern' class='on' onclick=\"dgView('modern')\">Modern</button>"
               "<button type='button' id='v-classic' onclick=\"dgView('classic')\">LigPlot-style</button>"
               + ("<button type='button' id='v-detailed' onclick=\"dgView('detailed')\">Detailed</button>" if d.get("detailed") else "") + "</span>"
               "<span class='seg' role='group' aria-label='Contacts drawn'>"
               "<button type='button' id='d-simple' class='on' onclick=\"dgDetail('simple')\">Simple</button>"
               "<button type='button' id='d-all' onclick=\"dgDetail('all')\">All contacts</button></span>"
               "<span class='dl'>Download <button type='button' onclick=\"dgSave('svg')\">SVG</button>"
               "<button type='button' onclick=\"dgSave('png')\">PNG 300 dpi</button>"
               "<button type='button' onclick=\"dgSave('pdf')\">PDF 300 dpi</button></span></div>")
    out.append("<div class='card' style='padding:10px' id='dgbox' data-name='" + esc(name.rsplit('.', 1)[0] + "_" + d["resname"]) + "'>"
               "<div id='dg-modern'>" + ligplot.svg_modern(dg) + "</div><div id='dg-classic' hidden>" + ligplot.svg(dg) + "</div>"
               "<div id='dg-modern-all' hidden>" + ligplot.svg_modern(dga) + "</div><div id='dg-classic-all' hidden>" + ligplot.svg(dga) + "</div>"
               + (("<div id='dg-detailed' hidden>" + ligdetail.svg_detail(dg, d.get("detailed_all")) + "</div><div id='dg-detailed-all' hidden>" + ligdetail.svg_detail(dga, d.get("detailed_all")) + "</div>")
                  if d.get("detailed") else "") + "</div>")
    keys = lambda g: {(n["ch"], n["num"], n["rn"]) for n in g["nodes"]}
    hid = keys(dga) - keys(dg)
    hw, hr = sum(1 for k in hid if k[2] in WATER_NAMES), sum(1 for k in hid if k[2] not in WATER_NAMES)
    out.append(f"<p class='legend'><b>Simple</b> view (default) draws a water only when it makes a hydrogen bond to the ligand, and a residue only when it has a polar contact "
               f"or two or more contacts of which the closest is within 3.8 A: {hr} residue{'s' if hr != 1 else ''} and {hw} water{'s' if hw != 1 else ''} are left out. "
               "<b>All contacts</b> draws everything within 4 A. The tables below always list every residue.</p>")
    out.append("<p class='legend'>The ligand is drawn flat and turned to match its orientation in the model, so residues sit on the side where they are in 3D. "
               "Modern view: coloured cards are residues (dot colour = residue type), green lines are hydrogen bonds with the heavy-atom distance in A, "
               "grey lines are the closest other contact of a residue, small blue circles are waters. Hover a residue to isolate its contacts. "
               "LigPlot-style view: dashed green lines are hydrogen bonds, dotted lines the closest other contact. "
               "Downloads are made in your browser from the view on screen, always on a white background. PNG and PDF are 300 dpi rasters.</p>")
    # 2 interactions
    out.append('</section><section id="s2" class="sec s2"><h2>2. Interactions</h2>')
    pol = sorted([c for c in close if c["type"] in ("hbond", "salt_bridge", "metal")], key=lambda c: c["distance"])
    out.append("<h3 id='polar'>Hydrogen bonds, salt bridges and metal contacts</h3>")
    if pol:
        rows = []
        atoms = []
        for c in pol:
            a_l, a_p = (c["atom1"], c["atom2"]) if c["unit1"] == d["lu"] else (c["atom2"], c["atom1"])
            rows.append([a_l.split(":")[-1], a_p, c["distance"], c["type"], c["verdict"], "water" if "water" in (c["unit1"], c["unit2"]) else "protein"])
            atoms.append([a_p, a_l])
        out.append(table(["ligand atom", "partner atom", "distance (Å)", "type", "verdict", "partner is"], rows, cls={4: lambda v: v}, atoms=atoms))
    else:
        out.append("<p class='note'>none</p>")
    out.append(f"<h3 id='mediated'>Hydrogen bonds via water: water bridges ({hb['water_mediated']}"
               + (f": {hb['water_mediated_clear']} clear, {hb['water_mediated_doubtful']} doubtful" if hb["water_mediated"] else "") + ")</h3>")
    if not r["flags"].get("n_water_residues", 0):
        out.append("<p class='note'>This model has no waters (as in an AI prediction), so no bonds via water can be counted.</p>")
    elif not hb["mediated_waters"]:
        out.append("<p class='note'>No water makes a hydrogen bond to both the ligand and the protein.</p>")
    else:
        out.append("<p class='note'><b>This table lists the bridging waters only</b>: waters that make a hydrogen bond to a ligand atom <b>and</b> to a protein atom "
                   "(distance in &Aring; after the atom). A water bonded to the ligand only is counted under H-bonds to waters (card above) and is drawn in the map when it "
                   "bonds the ligand, but it is not a bridge and is not listed here. " + H.MEDIATED_NOTE + "</p>"
                   + table(["water", "bonded to the ligand", "bonded to the protein", "check"],
                           [[x["water"], ", ".join(f"{y['atom'].split(':')[-1]} {y['distance']}" for y in x["side_a"]),
                             ", ".join(f"{y['atom']} {y['distance']}" for y in x["side_b"]), H._bridge_check_text(x)] for x in hb["mediated_waters"]],
                           cls={3: lambda v: "borderline" if v.startswith("doubtful") else "ok"},
                           atoms=[[x["water"]] + [y["atom"] for y in x["side_a"] + x["side_b"]] for x in hb["mediated_waters"]]))
    out.append("<h3 id='pocket'>Pocket residues</h3>")
    rows, atoms = [], []
    for n in sorted(dga["nodes"], key=lambda n: min(c["d"] for c in n["contacts"])):
        best = min(n["contacts"], key=lambda c: c["d"])
        kinds = sorted({c["type"] for c in n["contacts"]})
        rows.append([f"{n['ch']}:{n['num']}:{n['rn']}", len(n["contacts"]), best["d"], f"{best['atom']} to {best['lig']}", ", ".join(kinds),
                     round(n["b"], 1) if n.get("b") is not None else "-"])
        atoms.append([f"{n['ch']}:{n['num']}:{n['rn']}:{best['atom']}", f"{d['b']}:1:{d['resname']}:{best['lig']}"])
    out.append(table(["residue", "contacts", "closest (Å)", "atoms", "types", "pLDDT" if pred else "B-factor"], rows, atoms=atoms))
    out.append("<h3>Ligand atoms</h3>")
    cnt = {}
    for c in close:
        a_l = (c["atom1"] if c["unit1"] == d["lu"] else c["atom2"]).split(":")[-1]
        cnt.setdefault(a_l, [0, 0])
        cnt[a_l][0] += 1
        cnt[a_l][1] += 1 if c["type"] in ("hbond", "salt_bridge", "metal") else 0
    rows, atoms = [], []
    for j, k in enumerate(li):
        nm = m.atom_name[k]
        f0, f1 = float(free[j]), float(comp[j])
        rows.append([nm, m.element[k], round(f0, 1), round(f1, 1), f"{(1 - f1 / f0) * 100:.0f}%" if f0 > 0.01 else "-", cnt.get(nm, [0, 0])[0], cnt.get(nm, [0, 0])[1],
                     round(lp[j], 0) if pred else "-"])
        atoms.append([f"{d['b']}:1:{d['resname']}:{nm}"])
    out.append(table(["atom", "element", "exposed alone (Å²)", "exposed in complex (Å²)", "buried", "contacts", "polar contacts", "pLDDT" if pred else "B-factor"], rows, atoms=atoms))
    if unsat:
        out.append("<h3 id='unsat'>Buried polar ligand atoms without a partner</h3><p class='note'>A polar ligand atom that is mostly buried but has no hydrogen bond partner in the model. It may bind a water that is not modelled or it may be a weak point of the pose.</p>"
                   + table(["atom", "buried (Å²)"], [[x["atom"], x["buried"]] for x in unsat], atoms=[[x["atom"]] for x in unsat]))
    # 3 quality
    out.append('</section><section id="s3" class="sec s3"><h2>3. Contact quality</h2>')
    nclash = sum(1 for c in clashes if c["verdict"] == "unfavorable")
    out.append(f"<h3 id='lclash'>Ligand clashes ({len(clashes)}, of which {nclash} unfavorable)</h3><p class='note'>Ligand atoms that overlap protein atoms. Unfavorable means an overlap of 0.6 A or more, borderline 0.25 to 0.6 A.</p>" + H.clash_table(clashes))
    out.append(f"<h3 id='lissue'>Other contact issues ({len(chem)})</h3><p class='note'>Too short hydrogen bonds, polar atoms of the same kind facing each other, like charges close together.</p>" + H.chem_table(chem))
    out.append("<h3>What this report does not cover</h3><p class='note'>The geometry of the ligand itself (bond lengths, angles, planarity, stereochemistry, strain) is not checked, and the bond orders " + ("come from the file dictionary" if dg.get("bond_source") == "file dictionary" else "were guessed") + ". "
               "Hydrogens, protonation states, water molecules that are not in the model, and binding energy are not evaluated. " + comparison_short_html() + "</p></section>")
    out.append(H.NAV_JS)
    from .about import footer_html
    out.append(footer_html())
    out.append(DG_JS)
    if molstar:
        dd = {"path": d["path"], "model_index": d["model_index"], "a": d["a"], "b": d["b"], "kind": "ligand" if pred else "experimental",
              "residues": af3.residue_plddt(d["path"], d["model_index"]), "iface": it, "molstar_url": d.get("molstar_url")}
        out.append(H.molstar_block(dd))
    out.append("</main></body></html>")
    return "".join(out)


def build(path, out_path, ligand_unit=None, model_index=0, molstar_path=None, predicted=None, molstar_url=None, detailed=False, detailed_all=False):
    d = collect(path, ligand_unit, model_index, predicted, detailed, detailed_all)
    d["molstar_url"] = molstar_url
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(render(d))
    if molstar_path:
        with open(molstar_path, "w", encoding="utf-8") as fh:
            fh.write(render(d, True))
    return d
