#!/usr/bin/env python3
"""Figures for the benchmarks of ifacekit against PISA and MolProbity (development only; needs matplotlib).

All numbers are computed from the result files of the comparisons; nothing is typed by hand.
  fig1_pisa            interface area and interface residues against PISA (tools/compare_pisa_batch.py)
  fig2_clashes         all-atom clashes with hydrogens added against MolProbity: overlap of the two lists, count per structure and the causes
                       of the differences (tools/compare_hydrogen_clashes.py, tools/explain_hydrogen_disagreement.py)
  fig3_summary         what share of MolProbity's calls ifacekit finds and what share of its own calls MolProbity confirms, for clashes,
                       Ramachandran outliers and rotamer outliers (tools/compare_molprobity_batch.py)
Colours: blue = MolProbity (and the share of its calls that ifacekit finds), orange = ifacekit (and the share of its calls that MolProbity confirms);
text and axes are in ink, not in the series colours. Output: PNG (300 dpi), SVG and PDF for every figure.
Usage: python tools/make_benchmark_figures.py --pisa pisa_comparison.json --hydrogens hydrogen_clashes.json --disagreement hydrogen_disagreement.json \\
           --after results/cases.json --dates release_dates.tsv --out figures/
"""
import argparse
import json
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

BLUE, ORANGE = "#2a78d6", "#eb6834"          # categorical slots 1 and 2 of the validated palette
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e0", "#ffffff"
CUT = "2024-01-01"

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 8.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK, "axes.linewidth": 0.8, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})


def save(fig, out, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def panel_label(ax, letter, title):
    ax.text(-0.02, 1.07, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="right", color=INK)
    ax.text(0.0, 1.07, title, transform=ax.transAxes, fontsize=9, va="bottom", ha="left", color=INK)


def recessive_grid(ax, axis="both"):
    ax.grid(True, axis=axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def lens(r1, r2, d):
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    a = r1 * r1 * math.acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
    b = r2 * r2 * math.acos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2))
    c = 0.5 * math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    return a + b - c


def circles(ax, n_mp, n_ifk, shared):
    """Two circles with areas that follow the counts and an overlap that follows the shared pairs."""
    k = 1.0 / math.sqrt(max(n_mp, n_ifk) / math.pi)
    r1, r2 = k * math.sqrt(n_mp / math.pi), k * math.sqrt(n_ifk / math.pi)
    lo, hi = abs(r1 - r2), r1 + r2
    for _ in range(60):
        mid = (lo + hi) / 2
        if lens(r1, r2, mid) > shared * k * k:
            lo = mid
        else:
            hi = mid
    d = (lo + hi) / 2
    x1 = -(r1 + d + r2) / 2 + r1
    x2 = x1 + d
    for x, r, c in ((x1, r1, BLUE), (x2, r2, ORANGE)):
        ax.add_patch(Circle((x, 0), r, facecolor=c, alpha=0.30, edgecolor="none", zorder=1))
        ax.add_patch(Circle((x, 0), r, facecolor="none", edgecolor=c, linewidth=1.4, zorder=2))
    left_edge, right_edge = x2 - r2, x1 + r1
    ax.text((left_edge + right_edge) / 2, 0, f"{shared:,}", ha="center", va="center", fontsize=10, fontweight="bold", color=INK, zorder=3)
    ax.text((left_edge + right_edge) / 2, 0.22 * r1, "both", ha="center", va="bottom", fontsize=7.5, color=MUTED, zorder=3)
    gap = 0.06 * max(r1, r2)
    ax.text(x1 - r1 - gap, 0.10 * r1, f"{n_mp - shared:,}", ha="right", va="bottom", fontsize=10, fontweight="bold", color=INK)
    ax.text(x1 - r1 - gap, -0.04 * r1, "only\nMolProbity", ha="right", va="top", fontsize=7.5, color=MUTED, linespacing=1.15)
    ax.text(x2 + r2 + gap, 0.10 * r2, f"{n_ifk - shared:,}", ha="left", va="bottom", fontsize=10, fontweight="bold", color=INK)
    ax.text(x2 + r2 + gap, -0.04 * r2, "only\nifacekit", ha="left", va="top", fontsize=7.5, color=MUTED, linespacing=1.15)
    mid = (x1 + x2) / 2
    top = max(r1, r2) * 1.08
    ax.text(mid - 0.04 * r1, top, f"MolProbity: {n_mp:,}", ha="right", va="bottom", fontsize=8.2, color=INK)
    ax.text(mid + 0.04 * r1, top, f"ifacekit: {n_ifk:,}", ha="left", va="bottom", fontsize=8.2, color=INK)
    ax.set_xlim(x1 - r1 * 2.05, x2 + r2 * 2.05)
    ax.set_ylim(-max(r1, r2) * 1.12, max(r1, r2) * 1.42)
    ax.set_aspect("equal")
    ax.axis("off")
    for k, (c, text) in enumerate(((BLUE, f"{100 * shared / n_mp:.1f}% of MolProbity's pairs are also found by ifacekit"),
                                   (ORANGE, f"{100 * shared / n_ifk:.1f}% of ifacekit's pairs are also in MolProbity's list"))):
        y = -0.03 - 0.075 * k
        ax.text(0.0, y, "\u25a0", transform=ax.transAxes, ha="left", va="top", fontsize=9, color=c)
        ax.text(0.045, y, text, transform=ax.transAxes, ha="left", va="top", fontsize=7.8, color=INK)


def scatter_log(ax, x, y, xlabel, ylabel, lim):
    x, y = np.array(x, float), np.array(y, float)
    ax.plot(lim, lim, color=MUTED, linewidth=0.9, zorder=1)
    ax.scatter(x, y, s=14, color=ORANGE, alpha=0.55, edgecolor=SURFACE, linewidth=0.4, zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    recessive_grid(ax)


def fig_pisa(rows, out):
    ok = [r for r in rows if "error" not in r and r.get("pisa_area") and r.get("ifk_area") is not None]
    pa, ia = np.array([r["pisa_area"] for r in ok]), np.array([r["ifk_area"] for r in ok])
    diff = 100.0 * (ia - pa) / pa
    r_area = np.corrcoef(pa, ia)[0, 1]
    pr = [(r["pisa_nres"], r["ifk_nres"]) for r in ok if r.get("pisa_nres")]
    r_res = np.corrcoef([a for a, _ in pr], [b for _, b in pr])[0, 1]
    n_struct = len({r["id"] for r in ok})
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.75), gridspec_kw={"wspace": 0.42})
    ax = axs[0]
    scatter_log(ax, pa, ia, "PISA interface area (Å$^2$)", "ifacekit interface area (Å$^2$)", (80, 12000))
    panel_label(ax, "a", "Interface area")
    ax.text(0.97, 0.05, f"n = {len(ok)} interfaces\nr = {r_area:.3f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=7.8, color=INK)
    ax = axs[1]
    edges = np.arange(-40, 40.01, 2.5)
    ax.axvspan(-5, 5, color=GRID, zorder=0)
    ax.hist(np.clip(diff, -39.99, 39.99), bins=edges, color=ORANGE, edgecolor=SURFACE, linewidth=0.8, zorder=2)
    ax.axvline(0, color=MUTED, linewidth=0.9, zorder=3)
    ax.set_xlabel("difference to PISA (%)")
    ax.set_ylabel("interfaces")
    ax.set_xlim(-40, 40)
    panel_label(ax, "b", "Difference in area")
    w5, w10 = 100 * np.mean(np.abs(diff) <= 5), 100 * np.mean(np.abs(diff) <= 10)
    ax.text(0.60, 0.95, f"median |\u0394| {np.median(np.abs(diff)):.1f}%\n{w5:.0f}% within 5%\n{w10:.0f}% within 10%", transform=ax.transAxes, ha="left", va="top", fontsize=7.6, color=INK)
    ax = axs[2]
    scatter_log(ax, [a for a, _ in pr], [b for _, b in pr], "PISA interface residues", "ifacekit interface residues", (8, 450))
    panel_label(ax, "c", "Interface residues")
    ax.text(0.97, 0.05, f"n = {len(pr)}\nr = {r_res:.3f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=7.8, color=INK)
    fig.text(0.0, -0.10, f"{n_struct} X-ray structures; interfaces between two chains of the model itself, PISA area of 100 Å$^2$ or more, at most the four largest per structure.\n"
             "Interface area is the mean of the area buried on the two sides. In b the grey band is within 5%; differences beyond 40% are drawn at the edge.",
             fontsize=6.8, color=MUTED, ha="left", va="top", linespacing=1.5)
    save(fig, out, "fig1_pisa")
    return {"n": len(ok), "structures": n_struct, "r_area": float(r_area), "median_abs_diff_pct": float(np.median(np.abs(diff))), "within5": float(w5), "within10": float(w10), "r_res": float(r_res)}


CAUSE_MP = {"ligand_water_hydrogen": "hydrogens of waters and ligands", "hydrogen_bond_rule": "hydrogen-bond rule", "bond_rule": "bond-count rule",
            "near_threshold": "overlap just under 0.40 Å", "short_polar_rule": "short N/O pair", "placement": "hydrogen placed elsewhere",
            "heavy_protein": "two heavy atoms (protein)", "heavy_water_ligand": "two heavy atoms (water, ligand)", "other": "other"}
CAUSE_IFK = {"near_threshold": "overlap 0.40–0.45 Å", "heavy_water": "heavy atoms, one a water", "fixed_hydrogen": "hydrogen of fixed geometry",
             "choice_hydrogen": "hydrogen of a turned or flipped group", "heavy_protein": "two heavy atoms (protein)", "heavy_ligand": "two heavy atoms (ligand, ion)"}


def hbars(ax, counts, labels, total, color, title, letter):
    items = [(labels[k], v) for k, v in sorted(counts.items(), key=lambda kv: -kv[1]) if v][::-1]
    ys = np.arange(len(items))
    ax.barh(ys, [v for _, v in items], height=0.62, color=color, edgecolor="none", zorder=2)
    ax.set_yticks(ys)
    ax.set_yticklabels([l for l, _ in items], fontsize=7.6)
    for y, (_, v) in zip(ys, items):
        ax.text(v + total * 0.012, y, f"{v:,} ({100 * v / total:.0f}%)", va="center", ha="left", fontsize=7.4, color=INK)
    ax.set_xlim(0, max(v for _, v in items) * 1.42)
    ax.set_xlabel("pairs")
    recessive_grid(ax, "x")
    ax.tick_params(axis="y", length=0)
    panel_label(ax, letter, title)


def fig_clashes(hy, dis, out):
    ok = [r for r in hy["cases"] if not r.get("failed")]
    al = hy["summary"]["all"]
    fig = plt.figure(figsize=(7.4, 6.3))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.62, height_ratios=[1.0, 1.0])
    ax = fig.add_subplot(gs[0, 0])
    circles(ax, al["mp"], al["ours"], al["shared"])
    panel_label(ax, "a", "Clash pairs found by both programs")
    ax = fig.add_subplot(gs[0, 1])
    x, y = np.array([r["mp"] for r in ok], float), np.array([r["ours"] for r in ok], float)
    lim = 400
    ax.plot([0, lim], [0, lim], color=MUTED, linewidth=0.9, zorder=1)
    ax.scatter(np.clip(x, 0, lim), np.clip(y, 0, lim), s=13, color=ORANGE, alpha=0.55, edgecolor=SURFACE, linewidth=0.4, zorder=2)
    ax.set_xscale("symlog", linthresh=10)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_xlim(-0.6, lim)
    ax.set_ylim(-0.6, lim)
    ticks = [0, 1, 3, 10, 30, 100, 300]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_yticklabels([str(t) for t in ticks])
    ax.set_xlabel("MolProbity clash pairs per structure")
    ax.set_ylabel("ifacekit clash pairs per structure")
    recessive_grid(ax)
    r = np.corrcoef(x, y)[0, 1]
    ax.text(0.04, 0.95, f"{len(ok)} structures\nr = {r:.2f}", transform=ax.transAxes, ha="left", va="top", fontsize=7.8, color=INK)
    panel_label(ax, "b", "Count per structure")
    mo, uo = dis["mp_total"] - dis["shared"], dis["ours_total"] - dis["shared"]
    ax = fig.add_subplot(gs[1, 0])
    hbars(ax, dis["mp_only"], CAUSE_MP, mo, BLUE, f"Only MolProbity reports ({mo:,} pairs)", "c")
    ax = fig.add_subplot(gs[1, 1])
    hbars(ax, dis["ours_only"], CAUSE_IFK, uo, ORANGE, f"Only ifacekit reports ({uo:,} pairs)", "d")
    fig.text(0.0, -0.005, "All-atom clashes (overlap of 0.4 Å or more) with hydrogens added by ifacekit, compared by residue pair with MolProbity's clash list on "
             f"{al['cases']} X-ray structures,\nfirst conformer only. Causes in c and d are assigned by fixed rules (tools/explain_hydrogen_disagreement.py): attributions, not proofs.",
             fontsize=6.8, color=MUTED, ha="left", va="top", linespacing=1.5)
    save(fig, out, "fig2_clashes")
    return {"pairs_mp": al["mp"], "pairs_ifk": al["ours"], "shared": al["shared"], "r_per_structure": float(r)}


def agg(rs, kind, which):
    g = lambda k: sum(r[kind][which][k] for r in rs)
    return g("mp_out"), g("ifk_out"), g("both")


def fig_summary(after, dates, hy, out):
    recs = [r for r in json.load(open(after)) if "error" not in r]
    post = [r for r in recs if dates[r["id"]] >= CUT]
    c = lambda k: sum(r["clash"][k] for r in recs)
    rows = [("All-atom clashes\n(hydrogens added)", 100 * hy["summary"]["all"]["found_pct"] / 100, hy["summary"]["all"]["confirmed_pct"]),
            ("Heavy-atom clashes\n(contact tables)", 100.0 * c("heavy_hit_04") / c("heavy_n_04"), 100.0 * c("unf_in_mp") / c("ifk_unf"))]
    for kind, name in (("rama", "Ramachandran outliers"), ("rota", "Rotamer outliers")):
        mpo, ifo, both = agg(post, kind, "pre")
        rows.append((f"{name}\n(independent set)", 100.0 * both / mpo, 100.0 * both / ifo))
    rows = rows[::-1]
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    ys = np.arange(len(rows))
    h = 0.34
    ax.barh(ys + h / 2 + 0.01, [r[1] for r in rows], height=h, color=BLUE, edgecolor="none", zorder=2, label="share of MolProbity's calls that ifacekit also makes (found)")
    ax.barh(ys - h / 2 - 0.01, [r[2] for r in rows], height=h, color=ORANGE, edgecolor="none", zorder=2, label="share of ifacekit's calls that MolProbity also makes (confirmed)")
    for y, r in zip(ys, rows):
        ax.text(r[1] + 1.2, y + h / 2 + 0.01, f"{r[1]:.0f}%", va="center", fontsize=8, color=INK)
        ax.text(r[2] + 1.2, y - h / 2 - 0.01, f"{r[2]:.0f}%", va="center", fontsize=8, color=INK)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xlim(0, 112)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("percent")
    recessive_grid(ax, "x")
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.17), frameon=False, fontsize=7.6, handlelength=1.2)
    fig.text(0.0, -0.17, f"{len(recs)} high-resolution X-ray structures; the independent set is the {len(post)} entries released from {CUT[:4]}, scored against references built from earlier entries only.\n"
             "Heavy-atom clashes: found = MolProbity's clashes between heavy atoms flagged by ifacekit; confirmed = ifacekit's unfavorable residue pairs that are in MolProbity's list.",
             fontsize=6.8, color=MUTED, ha="left", va="top", linespacing=1.5)
    save(fig, out, "fig3_summary")
    return {name.split("\n")[0]: (round(f, 1), round(cf, 1)) for (name, f, cf) in rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pisa", required=True)
    ap.add_argument("--hydrogens", required=True)
    ap.add_argument("--disagreement", required=True)
    ap.add_argument("--after", required=True, help="cases.json of the MolProbity comparison (tools/compare_molprobity_batch.py)")
    ap.add_argument("--dates", required=True, help="release_dates.tsv: PDB ID and release date")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    dates = dict(l.split() for l in open(a.dates))
    hy = json.load(open(a.hydrogens))
    res = {"fig1": fig_pisa(json.load(open(a.pisa))["pairs"], a.out),
           "fig2": fig_clashes(hy, json.load(open(a.disagreement)), a.out),
           "fig3": fig_summary(a.after, dates, hy, a.out)}
    json.dump(res, open(os.path.join(a.out, "figure_numbers.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
