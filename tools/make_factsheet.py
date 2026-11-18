"""Build the fact sheet page (ifacekit/data/factsheet.html): the comparison with MolProbity and the list of test cases (development only).

All numbers on the page are computed from the result files of tools/compare_molprobity_batch.py and of the rotamer calibration; nothing is typed by hand.
Usage: python tools/make_factsheet.py --after results/cases.json --before results_old/cases.json --cases cases.txt --dates release_dates.tsv \\
           --cases-dir mp_compare --cal-list cal_cases_all.txt --calibration calibration_summary.txt --out ifacekit/data/factsheet.html [--failed 9sqk,...] [--hydrogens hydrogen_clashes.json] [--chimerax chimerax_comparison.json] [--disagreement hydrogen_disagreement.json]
--disagreement is the result of tools/explain_hydrogen_disagreement.py (why the two clash lists differ).
--chimerax is the result of tools/compare_chimerax_hydrogens.py (section with hydrogens added by ChimeraX).
--hydrogens is the result of tools/compare_hydrogen_clashes.py (section and case-table column for the clash check with hydrogens added).
--before is the same comparison made with the previous version (optional). --cases-dir holds one folder per case with m.cif (for the resolution).
"""
import argparse
import html
import json
import os
import re

import gemmi

from ifacekit import __version__
from ifacekit.about import COPYRIGHT, LICENCE
from ifacekit.htmlreport import CSS

CUT = "2024-01-01"


def pct(a, b, nd=1):
    return f"{100.0 * a / b:.{nd}f}%" if b else "n/a"


def agg(recs, key, sub=None):
    out = {}
    for r in recs:
        d = r["clash"] if sub is None else r[key][sub]
        for k, v in d.items():
            if isinstance(v, (int, float)):
                out[k] = out.get(k, 0) + v
    return out


def cat_totals(recs):
    cat = {}
    for r in recs:
        for k, (n, h) in r["clash"]["heavy_cat"].items():
            c = cat.setdefault(k, [0, 0])
            c[0] += n
            c[1] += h
    return cat


def tr(cells, head=False):
    t = "th" if head else "td"
    return "<tr>" + "".join(f"<{t}>{c}</{t}>" for c in cells) + "</tr>"


def table(head, rows, cls=""):
    return f"<div class='wrap'><table class='{cls}'><thead>{tr(head, True)}</thead><tbody>{''.join(tr(r) for r in rows)}</tbody></table></div>"


def _lens(r1, r2, d):
    """Area shared by two circles of radii r1 and r2 whose centres are d apart."""
    import math
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    a = r1 * r1 * math.acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
    b = r2 * r2 * math.acos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2))
    c = 0.5 * math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    return a + b - c


def venn_svg(n_mp, n_ifk, shared, title, right="ifacekit"):
    """Two circles whose areas follow the numbers of clash pairs and whose overlap follows the shared pairs."""
    import math
    k = 95.0 / math.sqrt(max(n_mp, n_ifk) / math.pi)
    r1, r2 = k * math.sqrt(n_mp / math.pi), k * math.sqrt(n_ifk / math.pi)
    target = shared * k * k
    lo, hi = abs(r1 - r2), r1 + r2
    for _ in range(60):                                   # the overlap shrinks as the centres move apart
        mid = (lo + hi) / 2
        if _lens(r1, r2, mid) > target:
            lo = mid
        else:
            hi = mid
    d = (lo + hi) / 2
    w, h = 560, 300
    x1 = (w - (r1 + d + r2)) / 2 + r1
    x2 = x1 + d
    cy = 138
    left_edge, right_edge = x2 - r2, x1 + r1             # the lens lies between these two x values
    lens_x = (left_edge + right_edge) / 2
    num = lambda x, y, text: f'<text x="{x:.0f}" y="{y}" text-anchor="middle" fill="currentColor" font-size="18" font-weight="600">{text}</text>'
    return (f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{html.escape(title)}" style="max-width:{w}px;width:100%;height:auto;color:var(--fg)">'
            f'<circle cx="{x1:.1f}" cy="{cy}" r="{r1:.1f}" fill="#4c8bf5" fill-opacity="0.38" stroke="#4c8bf5" stroke-width="2"/>'
            f'<circle cx="{x2:.1f}" cy="{cy}" r="{r2:.1f}" fill="#f5a14c" fill-opacity="0.38" stroke="#f5a14c" stroke-width="2"/>'
            + num((x1 - r1 + left_edge) / 2, cy + 6, f"{n_mp - shared:,}") + num(lens_x, cy + 6, f"{shared:,}")
            + num((right_edge + x2 + r2) / 2, cy + 6, f"{n_ifk - shared:,}")
            + f'<text x="{lens_x:.0f}" y="{cy - 16}" text-anchor="middle" fill="currentColor" font-size="12">both</text>'
            + f'<text x="12" y="18" fill="#2f6fe0" font-size="13" font-weight="600">MolProbity {n_mp:,}</text>'
            + f'<text x="{w - 12}" y="18" text-anchor="end" fill="#d9801f" font-size="13" font-weight="600">{right} {n_ifk:,}</text>'
            + f'<text x="12" y="{h - 44}" fill="#2f6fe0" font-size="13" font-weight="600">{100.0 * shared / max(n_mp, 1):.1f}% of MolProbity&#8217;s pairs are also found by {right}</text>'
            + f'<text x="12" y="{h - 26}" fill="#d9801f" font-size="13" font-weight="600">{100.0 * shared / max(n_ifk, 1):.1f}% of the {right} pairs are also in MolProbity&#8217;s list</text>'
            + f'<text x="{w / 2:.0f}" y="{h - 6}" text-anchor="middle" fill="currentColor" font-size="11.5">left: MolProbity only &#183; middle: both &#183; right: {right} only</text></svg>')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--after", required=True)
    ap.add_argument("--before")
    ap.add_argument("--cases", required=True)
    ap.add_argument("--dates", required=True)
    ap.add_argument("--cases-dir", required=True)
    ap.add_argument("--cal-list", required=True)
    ap.add_argument("--calibration", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--failed", default="")
    ap.add_argument("--disagreement", help="result file of tools/explain_hydrogen_disagreement.py (why the two clash lists differ)")
    ap.add_argument("--chimerax", help="result file of tools/compare_chimerax_hydrogens.py (hydrogens added by ChimeraX instead of ifacekit's own)")
    ap.add_argument("--hydrogens", help="result file of tools/compare_hydrogen_clashes.py (clashes with hydrogens added)")
    a = ap.parse_args()
    dates = dict(l.split() for l in open(a.dates))
    cases = open(a.cases).read().split()
    failed = [x for x in a.failed.split(",") if x]
    after = {r["id"]: r for r in json.load(open(a.after)) if "error" not in r}
    before = {r["id"]: r for r in json.load(open(a.before)) if "error" not in r} if a.before else {}
    hy = json.load(open(a.hydrogens)) if a.hydrogens else None
    hy_rows = {r["id"]: r for r in hy["cases"] if not r.get("failed")} if hy else {}
    cxd = json.load(open(a.chimerax)) if a.chimerax else None
    cx_rows = {r["id"]: r for r in cxd["cases"] if "error" not in r} if cxd else {}
    recs = list(after.values())
    post = [r for r in recs if dates[r["id"]] >= CUT]
    older = [r for r in recs if dates[r["id"]] < CUT]
    cal = open(a.cal_list).read().split()
    calib_text = open(a.calibration).read()

    # ---------------- clashes
    c1, c0 = agg(recs, "clash"), agg(list(before.values()), "clash") if before else None
    cat1, cat0 = cat_totals(recs), cat_totals(list(before.values())) if before else {}
    kinds = [("between two different residues, not neighbours", "other residue pair"), ("with a water", "water"),
             ("involving an alternate conformer", "alternate conformer"), ("between neighbouring residues (i and i+1)", "adjacent residue"),
             ("inside one residue", "same residue")]
    rows = []
    for label, k in kinds:
        n, h = cat1.get(k, [0, 0])
        row = [label, f"{h:,} of {n:,}", pct(h, n)]
        if before:
            n0, h0 = cat0.get(k, [0, 0])
            row.insert(1, f"{pct(h0, n0)}")
        rows.append(row)
    tot_n, tot_h = c1["heavy_n_04"], c1["heavy_hit_04"]
    row = ["all heavy-atom clashes", f"{tot_h:,} of {tot_n:,}", pct(tot_h, tot_n)]
    if before:
        row.insert(1, pct(c0["heavy_hit_04"], c0["heavy_n_04"]))
    rows.append(row)
    clash_head = ["Heavy-atom clashes of MolProbity, by kind"] + (["flagged, previous version"] if before else []) + ["flagged by ifacekit", "share"]
    clash_tab = table(clash_head, rows)
    overlap_rows = [[f"{th} A or more", f"{c1['heavy_n_' + tag]:,}", f"{c1['heavy_hit_' + tag]:,}", pct(c1["heavy_hit_" + tag], c1["heavy_n_" + tag])]
                    for tag, th in (("04", "0.4"), ("06", "0.6"), ("08", "0.8"))]
    ov_tab = table(["Overlap of the MolProbity clash", "clashes", "flagged by ifacekit", "share"], overlap_rows)

    def spear(vals_a, vals_b):
        from scipy.stats import spearmanr
        return spearmanr(vals_a, vals_b)[0]
    cs = [r["mp"]["clashscore"] for r in recs]
    k1 = [1000.0 * r["clash"]["ifk_unf"] / max(r["clash"]["n_atoms"], 1) for r in recs]
    sp1 = spear(k1, cs)
    sp2 = spear([r["clash"]["ifk_unf"] for r in recs], [r["clash"]["mp_heavy"] for r in recs])
    no_ifk = sum(1 for r in recs if r["clash"]["ifk_unf"] == 0)
    no_mp = sum(1 for r in recs if r["clash"]["mp_total"] == 0)

    # ---------------- Ramachandran and rotamers (independent group: released after the cut-off, reference built from earlier entries only)
    def oc(rs, kind, which):
        g = lambda k: sum(r[kind][which][k] for r in rs)
        return g("mp_out"), g("ifk_out"), g("both"), g("n")
    rows = []
    for label, rs, which in (("recent entries, reference built from earlier entries only (independent)", post, "pre"),
                             ("older entries (they are part of the reference set)", older, "packaged")):
        for kind, name in (("rama", "Ramachandran"), ("rota", "Rotamers")):
            mpo, ifo, both, n = oc(rs, kind, which)
            rows.append([name, label, f"{len(rs)}", f"{n:,}", f"{mpo}", f"{ifo}", f"{both}", pct(both, mpo), pct(both, ifo)])
    out_tab = table(["", "set of structures", "structures", "residues scored", "MolProbity outliers", "ifacekit outliers", "shared", "found", "confirmed"], rows)
    types = {}
    for r in post:
        for nm, (b, i, m) in r["rota"]["pre"]["by_res"].items():
            c = types.setdefault(nm, [0, 0, 0])
            c[0] += b
            c[1] += i
            c[2] += m
    trows = [[nm, b, i, m, pct(b, b + m, 0), pct(b, b + i, 0)] for nm, (b, i, m) in sorted(types.items(), key=lambda kv: -(sum(kv[1])))]
    type_tab = table(["Residue type", "shared", "only ifacekit", "only MolProbity", "found", "confirmed"], trows)
    rb = None
    if before:
        pb = [before[r["id"]] for r in post if r["id"] in before]
        mpo, ifo, both, n = oc(pb, "rota", "pre")
        rb = (mpo, ifo, both)

    # ---------------- calibration tables parsed from the calibration output
    cal_rows = []
    for m in re.finditer(r"^\s+(before|after) calibration[^:]*: rotalyze outliers (\d+), ifacekit outliers (\d+), shared (\d+): found ([\d.]+)%, confirmed ([\d.]+)%, three-class agreement ([\d.]+)%",
                         calib_text, re.M):
        cal_rows.append(m.groups())
    tests = re.split(r"\n(?=TEST [AB])", calib_text)
    cal_html = []
    for blk in tests[1:]:
        title = blk.split("\n", 1)[0]
        lines = re.findall(r"^\s+(before|after) calibration[^:]*: rotalyze outliers (\d+), ifacekit outliers (\d+), shared (\d+): found ([\d.]+)%, confirmed ([\d.]+)%, three-class agreement ([\d.]+)%", blk, re.M)
        cal_html.append(f"<h4>{html.escape(title)}</h4>" + table(["", "MolProbity outliers", "ifacekit outliers", "shared", "found", "confirmed", "three-class agreement"],
                        [[("before calibration" if w == "before" else "after calibration"), mo, io, sh, f"{f}%", f"{c}%", f"{ag}%"] for w, mo, io, sh, f, c, ag in lines]))

    # ---------------- case list
    case_rows, csv = [], []
    head = ["PDB ID", "released", "resolution (A)", "set", "non-water atoms", "MolProbity clashscore", "ifacekit unfavorable per 1000 atoms",
            "heavy-atom clashes flagged", "Ramachandran outliers MolProbity / ifacekit (shared)", "rotamer outliers MolProbity / ifacekit (shared)"]
    if hy:
        head = head + ["all-atom clash pairs MolProbity / ifacekit (shared), hydrogens added"]
    if cxd:
        head = head + ["with ChimeraX hydrogens, pairs (shared): as placed / X-ray distances / ChimeraX finder"]
    csv.append(head)
    for i in sorted(cases):
        rel = dates.get(i, "")
        group = "recent, independent" if rel >= CUT else "older, in reference set"
        try:
            st = gemmi.read_structure(os.path.join(a.cases_dir, i, "m.cif"))
            res = f"{st.resolution:.2f}" if st.resolution else "-"
        except Exception:
            res = "-"
        r = after.get(i)
        link = f"<a href='https://www.rcsb.org/structure/{i.upper()}' target='_blank' rel='noopener'>{i.upper()}</a>"
        if r is None:
            case_rows.append([link, rel, res, group, "-", "MolProbity failed (Phenix error on TLS selection text)", "-", "-", "-", "-"] + (["-"] if hy else []) + (["-"] if cxd else []))
            csv.append([i.upper(), rel, res, group, "", "MolProbity failed", "", "", "", ""] + ([""] if hy else []) + ([""] if cxd else []))
            continue
        which = "pre" if rel >= CUT else "packaged"
        c = r["clash"]
        ra, ro = r["rama"][which], r["rota"][which]
        vals = [f"{c['n_atoms']:,}", f"{r['mp']['clashscore']:.1f}", f"{1000.0 * c['ifk_unf'] / max(c['n_atoms'], 1):.2f}", f"{c['heavy_hit_04']} of {c['heavy_n_04']}",
                f"{ra['mp_out']} / {ra['ifk_out']} ({ra['both']})", f"{ro['mp_out']} / {ro['ifk_out']} ({ro['both']})"]
        if hy:
            h = hy_rows.get(i)
            vals.append(f"{h['mp']} / {h['ours']} ({h['shared']})" if h else "-")
        if cxd:
            x = cx_rows.get(i)
            vals.append(" / ".join(f"{x[k][0]} ({x[k][1]})" for k in ("chimerax_as_placed", "chimerax_shortened", "chimerax_clashes")) if x else "-")
        case_rows.append([link, rel, res, group] + vals)
        csv.append([i.upper(), rel, res, group] + [v.replace(",", "") for v in vals])
    case_tab = table(head, case_rows, "cases")
    hy_html, hy_cards, hy_nav, cx_html, dis_html = "", "", "", "", ""
    if a.disagreement:
        D = json.load(open(a.disagreement))
        MPC = {"ligand_water_hydrogen": "MolProbity puts hydrogens on waters, ligands and other residues that are not protein; ifacekit gives them none.",
               "hydrogen_bond_rule": "ifacekit counts the pair as a hydrogen bond (a polar hydrogen next to an acceptor at a hydrogen-bond angle) and MolProbity does not, or reports it differently.",
               "near_threshold": "The overlap is 0.30 to 0.40 A in ifacekit and 0.40 A or more in MolProbity: a small difference in where a hydrogen sits, or in how overlap is measured (Probe measures it on dots on the atom surfaces).",
               "bond_rule": "Hydrogens of one residue or of neighbouring residues whose heavy atoms are three bonds apart or fewer: ifacekit leaves them out (its bond-count rule), MolProbity counts some of them.",
               "short_polar_rule": "Two N or O atoms that overlap by less than 0.55 A: ifacekit treats them as a short hydrogen bond, MolProbity reports some of them.",
               "placement": "The hydrogen sits elsewhere in ifacekit's model: a turned hydroxyl, thiol or NH3 group, a flipped amide or ring, or another protonation.",
               "heavy_water_ligand": "Two heavy atoms, one in a water or ligand, without an overlap of 0.3 A in ifacekit.",
               "heavy_protein": "Two heavy atoms of the protein, without an overlap of 0.3 A in ifacekit.",
               "other": "Not assigned."}
        OUC = {"near_threshold": "The overlap is 0.40 to 0.45 A in ifacekit; MolProbity's value for the same atoms is just under 0.40 A (the same kind of difference as in the first table).",
               "heavy_water": "Two heavy atoms, one of them a water. In MolProbity the water carries hydrogens, so such a pair is a hydrogen bond there; ifacekit has no water hydrogens and only exempts overlaps below 0.55 A.",
               "fixed_hydrogen": "A hydrogen whose position follows from the geometry (C-H, backbone N-H) overlaps another atom by 0.45 A or more in ifacekit but not in MolProbity: the neighbouring group is turned or flipped differently, or Probe gives the overlap to another atom.",
               "choice_hydrogen": "A hydroxyl, thiol, NH3, amide or His hydrogen whose state ifacekit chose differently from Reduce.",
               "heavy_ligand": "Two heavy atoms, one in a ligand or ion (different radii).",
               "heavy_protein": "Two heavy atoms of the protein."}
        def drows(counts, texts, total):
            return [[texts[k], f"{v:,}", f"{100.0 * v / max(total, 1):.0f}%"] for k, v in sorted(counts.items(), key=lambda kv: -kv[1]) if v]
        mo, uo = D["mp_total"] - D["shared"], D["ours_total"] - D["shared"]
        dis_html = (f"<h3>Why the two lists differ</h3><p class='note'>{D['shared']:,} pairs are in both lists. {mo:,} pairs ({100.0 * mo / D['mp_total']:.0f}% of MolProbity's) are only in MolProbity's list and "
                    f"{uo:,} ({100.0 * uo / D['ours_total']:.0f}% of ifacekit's) only in ifacekit's. Each pair was given the most likely cause by rules (tools/explain_hydrogen_disagreement.py); the causes are attributions, not proofs. "
                    "Both programs make choices: MolProbity places hydrogens with Reduce and finds overlaps with Probe's dots, ifacekit places the hydrogens itself and measures overlap between atom centres. "
                    "The two lists were never going to be identical.</p>"
                    f"<h4>{mo:,} pairs only MolProbity reports</h4>" + table(["Most likely cause", "pairs", "share"], drows(D["mp_only"], MPC, mo))
                    + f"<h4>{uo:,} pairs only ifacekit reports</h4>" + table(["Most likely cause", "pairs", "share"], drows(D["ours_only"], OUC, uo)))
    if cxd:
        cs = cxd["summary"]
        V = cs["variants"]
        crow = [[v["label"], f"{v['pairs']:,}", f"{v['shared']:,}", f"{v['found_pct']}%", f"{v['confirmed_pct']}%"] for v in V.values()]
        ctab = table(["Hydrogens and clash rules", "clash pairs found", "shared with MolProbity", "found", "confirmed"], crow)
        cven = "<div class='venns small'>" + "".join(
            f"<figure>{venn_svg(cs['mp_pairs'], v['pairs'], v['shared'], v['label'], right='this variant' if k != 'ifacekit' else 'ifacekit')}<figcaption>{html.escape(v['label'])}: <b>{v['shared']:,}</b> of {cs['mp_pairs']:,} = <b>{v['found_pct']}%</b></figcaption></figure>"
            for k, v in V.items()) + "</div>"
        cx_html = f"""<section id="chimerax" class="sec s3"><h2>Hydrogens added by ChimeraX instead</h2>
<p class="note">A check with another source of hydrogens, not part of ifacekit: ChimeraX {html.escape(cs['chimerax_version'])} added hydrogens (<code>addh hbond true</code>, which also protonates waters and ligands) to the same
{cs['cases']} models (first conformer, hydrogens and other conformers removed first), and the clashes were found in three ways and compared with the same {cs['mp_pairs']:,} MolProbity pairs. ChimeraX placed {cs['hydrogens_chimerax']:,} hydrogens,
ifacekit {cs['hydrogens_ifacekit']:,}. ChimeraX puts C-H at 1.09 A (the nuclear distance), while the radii of the clash rules belong to the shortened X-ray distances that MolProbity uses; the second line moves the ChimeraX
hydrogens to those distances. The last line is the clash finder of ChimeraX itself (overlap 0.4 A, hydrogen-bond allowance 0.6 A). ChimeraX does not flip Asn, Gln and His.</p>
{ctab}
{cven}
<p class="note">ChimeraX hydrogens do not give a better agreement than ifacekit's own: with the clash rules of ifacekit they find about as many of MolProbity's pairs only at the cost of many more pairs that MolProbity does not report. The
files (the models with hydrogens and the ChimeraX clash lists) are kept with the other comparison data. Method: <code>tools/compare_chimerax_hydrogens.py</code>.</p>
</section>
"""
    if hy:
        S = hy["summary"]
        al = S["all"]
        hrows = []
        for label, k in (("all cases", "all"), ("settings tuned on these (the first 40 cases)", "tuning"), ("independent: not used for tuning", "independent")):
            v = S.get(k)
            if v:
                hrows.append([label, v["cases"], f"{v['mp']:,}", f"{v['mp_h']:,}", f"{v['ours']:,}", f"{v['shared']:,}", f"{v['found_pct']}%", f"{v['found_h_pct']}%",
                              f"{v['found_heavy_pct']}%", f"{v['confirmed_pct']}%", v["correlation"]])
        hy_tab = table(["Structures", "cases", "MolProbity clash pairs", "of which with a hydrogen", "ifacekit clash pairs", "shared", "found", "found, with a hydrogen",
                        "found, heavy atoms only", "confirmed", "correlation of the counts per structure"], hrows)
        hy_cards = (f'<div class="card"><div class="k">clashes found, hydrogens added</div><div class="v">{al["found_pct"]:.0f}%</div><div class="n">{al["shared"]:,} of {al["mp"]:,} MolProbity clash pairs; '
                    f'{al["found_h_pct"]:.0f}% of those with a hydrogen</div></div>'
                    f'<div class="card"><div class="k">clashes confirmed, hydrogens added</div><div class="v">{al["confirmed_pct"]:.0f}%</div><div class="n">{al["shared"]:,} of {al["ours"]:,} ifacekit clash pairs are in MolProbity\'s list</div></div>')
        ind = S.get("independent") or {}
        def sentence(v, what):
            return (f"<p class='big'>{what}: ifacekit finds <b>{v['shared']:,}</b> of the <b>{v['mp']:,}</b> clash pairs of MolProbity = <b>{v['found_pct']}%</b>. "
                    f"Of the <b>{v['ours']:,}</b> pairs ifacekit reports, <b>{v['shared']:,}</b> ({v['confirmed_pct']}%) are also in MolProbity's list.</p>")
        venns = (f"<div class='venns'><figure>{venn_svg(al['mp'], al['ours'], al['shared'], 'all cases')}<figcaption>all {al['cases']} cases</figcaption></figure>"
                 + (f"<figure>{venn_svg(ind['mp'], ind['ours'], ind['shared'], 'independent cases')}<figcaption>the {ind['cases']} cases not used for tuning</figcaption></figure>" if ind else "")
                 + "</div>")
        hy_nav = '<a href="#hclashes">Hydrogens</a>'
        if cxd:
            hy_nav += '<a href="#chimerax">ChimeraX</a>'
        hy_html = f"""<section id="hclashes" class="sec s2"><h2>Clashes with hydrogens added</h2>
<p class="note">With the option <code>--hydrogens</code> (or the check box of the portal) ifacekit places hydrogens on the protein itself, does not use Reduce, and reports all-atom clashes: an overlap of
0.4 A or more between atoms more than four bonds apart (the radii of the Probe program), a polar hydrogen next to an acceptor counting as a hydrogen bond. The hydrogen positions of the non-rotating groups
agree with Reduce's (97% within 0.15 A on three structures, and 96% of the Asn and Gln flips on 40 cases, see the design document, section 5.13); methyl groups are staggered, hydroxyl, thiol and NH3 groups are turned, and Asn, Gln and His are flipped
and protonated to the state with the best contacts. Clashes are compared by residue pair (two atoms of the same two residues are one pair) with MolProbity's severe clashes, those between heavy atoms and those
that involve a hydrogen. The settings (hydrogen-bond overlap limits, the bond count that excludes close atoms, the gain a flip must show) were fitted on the first 40 cases; the other {hy['summary']['independent']['cases']
if hy['summary'].get('independent') else 0} cases were not used for any setting.</p>
{sentence(al, 'All ' + str(al['cases']) + ' cases')}
{sentence(ind, 'The ' + str(ind['cases']) + ' cases that were not used for tuning') if ind else ''}
{venns}
{hy_tab}
{dis_html}
<p class="note">Alternate conformers are left out of this comparison: the hydrogen check works on the first conformer of every atom, so the {al['mp_alt_excluded']:,} further clash pairs of MolProbity that involve an atom of another conformer are not counted. found = share of MolProbity's clash pairs that ifacekit also finds; confirmed = share of ifacekit's clash pairs that MolProbity also reports. The agreement is the same on the cases that were
not used for tuning. Not found: clashes of hydrogens that ifacekit does not place (ligands, waters, nucleic acids, modified residues), clashes between atoms of one residue and its neighbour that MolProbity counts
but the bond-count rule excludes, and groups where Reduce makes a different choice of flip, protonation or rotation. Clashes between two heavy atoms are found less well here (the heavy-atom contact tables of the
report cover them better, see the section above). This is not MolProbity's clashscore, only a comparison of the clash lists; the count per structure follows MolProbity's closely (correlation above).</p>
</section>
"""
    csv_js = json.dumps(csv)

    body = f"""
<h1>ifacekit fact sheet: comparison with MolProbity</h1>
<p class="sub">Valid for ifacekit {html.escape(__version__)}. Heavy-atom contacts, Ramachandran and rotamer outliers compared with MolProbity (Phenix 2.0, <code>phenix.molprobity</code>) on
{len(recs)} high-resolution X-ray structures. All numbers on this page are computed from the result files of the comparison. How the numbers were made, their limits and the
calibration are in the design document (sections 5.11 and 5.12).</p>
<nav class="top"><a href="#summary">Summary</a><a href="#clashes">Clashes</a>{hy_nav}<a href="#outliers">Ramachandran and rotamers</a><a href="#calibration">Calibration</a>
<a href="#limits">Limits</a><a href="#cases">Test cases</a><span class="sp"></span>
<button class="theme" onclick="var r=document.documentElement;r.dataset.theme=(r.dataset.theme==='dark'?'light':'dark')">light / dark</button></nav>

<section id="summary" class="sec s1"><h2>Summary</h2>
<div class="cards">
<div class="card"><div class="k">heavy-atom clashes flagged</div><div class="v">{pct(tot_h, tot_n, 0)}</div><div class="n">{tot_h:,} of {tot_n:,} MolProbity clashes between heavy atoms</div></div>
<div class="card"><div class="k">between two different residues</div><div class="v">{pct(c1['heavy_inter_hit'], c1['mp_heavy_inter'], 0)}</div><div class="n">{c1['heavy_inter_hit']:,} of {c1['mp_heavy_inter']:,} clashes, neighbours and waters included</div></div>
<div class="card"><div class="k">unfavorable pairs confirmed</div><div class="v">{pct(c1['unf_in_mp'], c1['ifk_unf'], 0)}</div><div class="n">{c1['unf_in_mp']:,} of {c1['ifk_unf']:,} are in MolProbity's severe list</div></div>
{hy_cards}
<div class="card"><div class="k">Ramachandran outliers found</div><div class="v">{pct(oc(post, 'rama', 'pre')[2], oc(post, 'rama', 'pre')[0], 0)}</div><div class="n">independent set; {pct(oc(post, 'rama', 'pre')[2], oc(post, 'rama', 'pre')[1], 0)} of ifacekit's confirmed</div></div>
<div class="card"><div class="k">rotamer outliers found</div><div class="v">{pct(oc(post, 'rota', 'pre')[2], oc(post, 'rota', 'pre')[0], 0)}</div><div class="n">independent set; {pct(oc(post, 'rota', 'pre')[2], oc(post, 'rota', 'pre')[1], 0)} of ifacekit's confirmed</div></div>
</div>
<p class="note">Method in short. Cases: random entries of a list of 13,115 PDB entries (X-ray, resolution 1.8 A or better, R-free 0.25 or lower, one per 30% sequence-identity cluster; files up to 450 kB, 1 to 8 chains), {len(post) + sum(1 for i in failed if dates.get(i, '') >= CUT)} released from 2024-01-01 and {len(older) + sum(1 for i in failed if dates.get(i, '') < CUT)} before; the structures used while the thresholds were set are excluded.
MolProbity ran on {len(recs)} of the {len(cases)} (it failed on {len(failed)}). Clashes are compared by residue pair with MolProbity's severe clashes (overlap 0.4 A or more after Reduce adds hydrogens), Ramachandran and rotamer outliers by residue.
For the recent entries the Ramachandran and rotamer references were built from entries released before 2024 only, so these two comparisons are independent of the references.</p>
</section>

<section id="clashes" class="sec s2"><h2>Clashes</h2>
<p class="note">MolProbity reported {c1['mp_total']:,} severe clashes in {len(recs)} structures: {c1['mp_heavy']:,} between heavy atoms and {c1['mp_h']:,} involving a hydrogen, which these heavy-atom tables do not evaluate (see the section on clashes with hydrogens added).
A clash counts as found when ifacekit flags (borderline or unfavorable) a contact between the same two residues.</p>
{clash_tab}
<h3>By overlap</h3>{ov_tab}
<p class="note">The share found falls with the overlap because the largest overlaps are dominated by clashes inside one residue and between neighbouring residues, which ifacekit does not evaluate (it skips atom pairs up to five bonds apart).</p>
<h3>The other way round, and per structure</h3>
<p class="note">ifacekit's {c1['ifk_unf']:,} unfavorable residue pairs: {c1['unf_in_mp']:,} ({pct(c1['unf_in_mp'], c1['ifk_unf'])}) are in MolProbity's severe list, {pct(c1['unf_in_mp_heavy'], c1['ifk_unf'])} in its heavy-atom list.
Its {c1['ifk_border']:,} borderline pairs agree in {pct(c1['border_in_mp'], c1['ifk_border'])}: borderline is a "look at it" category and not a clash list. Of the {c1['h_n']:,} MolProbity clashes that involve a hydrogen, the residue pair is
flagged in {pct(c1['h_pair_hit'], c1['h_n'])}. Per structure, ifacekit's unfavorable pairs per 1000 atoms follow MolProbity's clashscore with a Spearman correlation of {sp1:.2f}, and the count of unfavorable pairs follows MolProbity's heavy-atom clash count
with {sp2:.2f}; in {no_ifk} of {len(recs)} structures ifacekit has no unfavorable pair although MolProbity reports severe clashes in all but {no_mp}. ifacekit is not a clashscore.</p>
</section>

{hy_html}{cx_html}
<section id="outliers" class="sec s3"><h2>Ramachandran and rotamers</h2>
{out_tab}
<p class="note">found = share of MolProbity's outliers that ifacekit also reports; confirmed = share of ifacekit's outliers that MolProbity also reports. Rotamers: every alternate conformer is scored (MolProbity does the same);
Ramachandran uses the first conformer. There are few Ramachandran outliers in such structures (the numbers are small).</p>
{('<p class="note">Rotamers in the recent set with the previous version of ifacekit: MolProbity ' + str(rb[0]) + ' outliers, ifacekit ' + str(rb[1]) + ', shared ' + str(rb[2]) + ' (found ' + pct(rb[2], rb[0]) + ', confirmed ' + pct(rb[2], rb[1]) + ').</p>') if rb else ''}
<h3>Rotamer outliers by residue type (recent set)</h3>{type_tab}
</section>

<section id="calibration" class="sec s1"><h2>Calibration of the rotamer cut-offs</h2>
<p class="note">The rotamer cut-offs were calibrated against MolProbity's <code>rotalyze</code> on {len(cal)} calibration structures that are not among the comparison cases ({sum(1 for i in cal if dates.get(i, '') < CUT)} released before 2024, {sum(1 for i in cal if dates.get(i, '') >= CUT)} from 2024),
with reference grids built without the calibration entries, and then tested on {len(post) + sum(1 for i in failed if dates.get(i, '') >= CUT)} recent entries (test A) and {len(older) + sum(1 for i in failed if dates.get(i, '') < CUT)} older entries that were left out of the grids that scored them (test B). Residues without alternate conformers only.
The cut-offs are therefore fitted to MolProbity's calls.</p>
{''.join(cal_html)}
<details><summary>Full output of the calibration run (per type)</summary><pre>{html.escape(calib_text)}</pre></details>
</section>

<section id="limits" class="sec s2"><h2>What this comparison does not show</h2>
<ul>
<li>All structures are high-resolution X-ray entries of the kind that was used to build the references. Low-resolution models, cryo-EM models, AI predictions and nucleic acids were not tested, and the differences may be larger there.</li>
<li>MolProbity's clashes that involve hydrogens ({pct(c1['mp_h'], c1['mp_total'], 0)} of all of them) are not evaluated by the heavy-atom tables of the reports; the option <code>--hydrogens</code> evaluates them (section "Clashes with hydrogens added"), except for ligands, waters and nucleic acids, which get no hydrogens.</li>
<li>Clashes inside one residue and between neighbouring residues are only partly found, by design.</li>
<li>The rotamer cut-offs were calibrated against MolProbity; Asn is still over-called on recent entries; Cys, His, Pro and Trp keep uncalibrated cut-offs (too few outliers). Ramachandran, interface areas, bonds and the flip check use the first conformer.</li>
<li>MolProbity failed on {len(failed)} cases ({html.escape(', '.join(x.upper() for x in failed))}): Phenix stopped on TLS selections written in braces in the file header. They are listed below.</li>
</ul></section>

<section id="cases" class="sec s3"><h2>The test cases ({len(cases)})</h2>
<p class="note">Comparison cases, with the results per structure. "recent, independent": released from {CUT}, scored for Ramachandran and rotamers with references built from earlier entries only; "older, in reference set": released earlier and part of the reference set. The {len(cal)} calibration structures are listed at the end.
<button type="button" id="dl">Download the table as CSV</button> <input id="flt" placeholder="filter, for example 2024 or recent" size="26"></p>
{case_tab}
<details><summary>The {len(cal)} calibration structures</summary><p class="note">{', '.join(i.upper() for i in sorted(cal))}</p></details>
</section>
<style>.venns.small figure{{flex:1 1 220px}}p.big{{font-size:15px;margin:8px 0}}.venns{{display:flex;flex-wrap:wrap;gap:18px;margin:10px 0}}.venns figure{{margin:0;flex:1 1 300px;text-align:center}}.venns figcaption{{font-size:12.5px;color:var(--mut)}}table.cases{{font-size:12px}}table.cases td,table.cases th{{white-space:nowrap}}pre{{white-space:pre-wrap;font-size:12px;background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px}}
#dl{{font:inherit;font-size:13px;padding:4px 12px;border:1px solid var(--line);background:var(--card);color:var(--fg);border-radius:6px;cursor:pointer}}#flt{{font:inherit;padding:4px 8px;border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--fg)}}
a{{color:var(--acc)}}details{{margin:10px 0}}summary{{cursor:pointer}}</style>
<script>
var CSV={csv_js};
document.getElementById('dl').onclick=function(){{var t=CSV.map(function(r){{return r.map(function(x){{x=String(x);return /[",\\n]/.test(x)?'"'+x.replace(/"/g,'""')+'"':x}}).join(',')}}).join('\\n');var a=document.createElement('a');a.href=URL.createObjectURL(new Blob([t],{{type:'text/csv'}}));a.download='ifacekit_comparison_cases.csv';document.body.appendChild(a);a.click();a.remove()}};
document.getElementById('flt').oninput=function(){{var q=this.value.toLowerCase();document.querySelectorAll('table.cases tbody tr').forEach(function(r){{r.style.display=r.textContent.toLowerCase().indexOf(q)>=0?'':'none'}})}};
</script>
<footer style="margin-top:32px;padding-top:10px;border-top:1px solid var(--line);font-size:12.5px;color:var(--mut)">ifacekit {html.escape(__version__)} &middot; {LICENCE} &middot; {COPYRIGHT}.</footer>
"""
    doc = (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
           f"<title>ifacekit fact sheet</title><style>{CSS}</style></head><body><main>{body}</main></body></html>")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write(doc)
    print(f"fact sheet written to {a.out}: {len(cases)} cases ({len(recs)} with MolProbity results), {len(cal)} calibration structures, {len(doc) // 1024} kB")


if __name__ == "__main__":
    main()
