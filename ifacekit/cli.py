import argparse
import json
import sys

from . import __version__
from .attach import plan_attachments
from .consistency import consistency_text, metal_consistency
from .core import analyze, list_chains
from .report import chain_list_text, quality_overview, report_quality, text_report


def _fig_base(html):
    b = html[:-5] if html.lower().endswith(".html") else html
    return (b[:-7] if b.endswith("_report") else b) + "_map"


def _plots_pdf(structure, chains, out, model, models=None):
    from .plotspdf import build as pdf_build
    try:
        n = pdf_build(structure, chains, out, model, models=models)
        print(f"Plots PDF written to {out} ({n} pages)", file=sys.stderr)
        return True
    except RuntimeError as e:
        print(f"ifacekit: PDF not written: {e}", file=sys.stderr)
        return False


def _pdf_link(pdf_ok, plots_pdf, html_path):
    """Name of the PDF for a download link inside the HTML report: only when it was written next to the report."""
    import os
    if pdf_ok and plots_pdf and html_path and os.path.dirname(os.path.abspath(plots_pdf)) == os.path.dirname(os.path.abspath(html_path)):
        return os.path.basename(plots_pdf)
    return None


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="ifacekit",
        description="Contact, clash and interface analysis for PDB/mmCIF models. "
                    "Say what to inspect: --chains (quality of chains) or --pair (interaction of two chains).")
    p.add_argument("structure", nargs="?", help="PDB or mmCIF file")
    what = p.add_mutually_exclusive_group()
    what.add_argument("--chains", nargs="+", metavar="ID",
                      help="inspect the quality of each listed chain on its own; "
                           "use the word all for every chain")
    what.add_argument("--pair", nargs=2, metavar=("A", "B"),
                      help="inspect the interaction between two chains")
    what.add_argument("--all", action="store_true",
                      help="analyse every chain and ligand in the file")
    p.add_argument("--attach", choices=["contact", "label"], default="contact",
                   help="how ligands, glycans and waters are assigned to chains: "
                        "contact (default; by covalent bond or proximity, ignoring their chain label) "
                        "or label (only those carrying the chain's own label)")
    p.add_argument("--html", metavar="FILE", help="with --pair: write a self-contained HTML report "
                   "(AI metrics, interactions, prediction quality); with --chains: HTML quality report of the chains")
    p.add_argument("--experimental", action="store_true", help="with --html: report for an experimental structure "
                   "(interactions and model quality, no AI metrics, B-factors instead of pLDDT). The model type is "
                   "detected from the file (pair and ligand reports), and this option overrides it")
    p.add_argument("--predicted", action="store_true", help="with --pair or --ligand: treat the file as a predicted model "
                   "(read the confidence column as pLDDT), overriding the automatic detection")
    p.add_argument("--detailed", action="store_true", help="with --ligand: also draw a detailed diagram in which the residues that "
                   "hydrogen-bond or salt-bridge the ligand are shown as chemical structures (a third view in the report, and figure files)")
    p.add_argument("--detailed-all", action="store_true", help="with --detailed: draw every standard residue of the diagram as a structure, "
                   "also those with only van der Waals or hydrophobic contacts (default: only hydrogen-bonding and salt-bridging residues)")
    p.add_argument("--ligand", nargs="?", const="auto", metavar="ID", help="with --html: report for one ligand bound to the protein "
                   "(interaction diagram, tables, contact quality). ID is the ligand unit as shown in the chain list; "
                   "the first ligand is used when it is left out")
    p.add_argument("--molstar", action="store_true", help="with --html: also write FILE_3d.html, the same report with a "
                   "Mol* 3D panel (click a table row to see it in 3D; needs internet for Mol* itself). The plain report is unchanged")
    p.add_argument("--hydrogens", action="store_true", help="also check clashes with hydrogens added: ifacekit places the hydrogens on the protein "
                   "(waters, ligands and nucleic acids get none) and reports the all-atom clashes. Only this check uses the hydrogens; "
                   "every other result comes from the file as given. Calibrated against MolProbity (see the fact sheet)")
    p.add_argument("--rama-html", metavar="FILE", help="with --chains: write an HTML page with Ramachandran plots "
                   "for the chosen chains (ifacekit reference regions)")
    p.add_argument("--rotamer-html", metavar="FILE", help="with --chains: write an HTML page with the rotamer plots, "
                   "outliers and the Asn/Gln/His flip check for the chosen chains")
    p.add_argument("--plots-pdf", metavar="FILE", help="with --chains or --pair: write one PDF with all Ramachandran and rotamer "
                   "plots of the chosen chains (needs cairosvg and pypdf: pip install 'ifacekit[pdf]')")
    p.add_argument("--molstar-url", metavar="URL", help="with --molstar: base address of a local copy of Mol* (molstar.js and molstar.css) "
                   "instead of the jsDelivr CDN, e.g. /static/molstar/ on the portal (the page then works without internet)")
    p.add_argument("--factsheet-url", metavar="URL", help="address of the fact sheet page (the comparison with MolProbity and the list of test cases); "
                   "the notes of the HTML reports then link to it. The portal serves it at /factsheet")
    p.add_argument("--factsheet", metavar="FILE", help="write the fact sheet page to FILE (open it in a browser) and stop; no structure needed")
    p.add_argument("--af3-summary", metavar="JSON", help="AlphaFold 3 summary_confidences file for the HTML report")
    p.add_argument("--af3-full", metavar="JSON", help="AlphaFold 3 full_data file (PAE) for ipSAE and interface PAE")
    p.add_argument("-o", "--output", help="write full JSON result here")
    p.add_argument("--model", type=int, default=0, help="model index (default 0)")
    p.add_argument("--ensemble", action="store_true", help="for a file with several models (NMR ensemble): analyse every model and give "
                   "the numbers of each model, their mean and standard deviation, and the complete analysis of the first model. "
                   "Works with --chains, --pair and --all; the HTML, Ramachandran, rotamer and PDF pages are for the first model")
    p.add_argument("--models", metavar="LIST", help="which models get a report: first (default), all, or model numbers as the PDB file counts them (MODEL 1, 2, ...), "
                   "for example 3, 1,4,7 or 2-5. One model: the report of that model. Several models: with --html one page with a summary of the ensemble, a button for each model and one "
                   "for all models below each other. Implies --ensemble")
    p.add_argument("--separate", action="store_true", help="with --models and --html: write one HTML file per model (FILE_model03.html) instead of one page")
    p.add_argument("--ensemble-csv", metavar="FILE", help="with --ensemble: write one row per model plus mean, sd, min and max to a CSV file "
                   "(with several chains, FILE gets the chain name before the extension)")
    p.add_argument("--workers", type=int, default=1, help="with --ensemble: models analysed in parallel (default 1)")
    p.add_argument("--min-area", type=float, default=10.0,
                   help="smallest interface area to report, A^2 (default 10)")
    p.add_argument("--points", type=int, default=960,
                   help="sphere points for SASA (default 960)")
    p.add_argument("--all-intra", action="store_true",
                   help="also record ordinary intra-chain contacts (large output)")
    p.add_argument("-q", "--quiet", action="store_true", help="no text summary")
    p.add_argument("--version", action="version", version=__version__)
    a = p.parse_args(argv)
    if a.factsheet:
        import shutil
        from .about import factsheet_path
        shutil.copyfile(factsheet_path(), a.factsheet)
        print(f"fact sheet written to {a.factsheet}", file=sys.stderr)
        return 0
    if not a.structure:
        p.error("the structure file is required (or use --factsheet FILE)")
    if a.factsheet_url:
        from . import about
        about.FACTSHEET_URL = a.factsheet_url

    if a.ligand and a.html:
        from .ligreport import build as lig_build
        ms = (a.html[:-5] if a.html.lower().endswith(".html") else a.html) + "_3d.html" if a.molstar else None
        try:
            ld = lig_build(a.structure, a.html, None if a.ligand == "auto" else a.ligand, a.model, ms, molstar_url=a.molstar_url,
                           predicted=True if a.predicted else False if a.experimental else None, detailed=a.detailed,
                           detailed_all=a.detailed_all)
            from .ligplot import write_figures
            figs, note = write_figures(ld["diagram"], _fig_base(a.html))
            figs_all, _ = write_figures(ld["diagram_all"], _fig_base(a.html) + "_all")
            if a.detailed:
                from .ligdetail import svg_detail
                dv = (("detailed", lambda g: svg_detail(g, a.detailed_all)),)
                figs += write_figures(ld["diagram"], _fig_base(a.html), views=dv)[0]
                figs_all += write_figures(ld["diagram_all"], _fig_base(a.html) + "_all", views=dv)[0]
            print("Figures written: " + ", ".join(figs + figs_all), file=sys.stderr)
            if note:
                print("ifacekit: " + note, file=sys.stderr)
        except ValueError as e:
            print(f"ifacekit: {e}", file=sys.stderr)
            return 2
        print(f"Ligand report written to {a.html}" + (f" and {ms}" if ms else ""), file=sys.stderr)
        return 0
    if not (a.chains or a.pair or a.all):
        print(chain_list_text(a.structure, list_chains(a.structure, a.model)))
        return 2

    if a.html and not (a.pair or a.chains):
        print("ifacekit: --html needs --pair A B or --chains", file=sys.stderr)
        return 2
    kw = dict(model_index=a.model, min_area=a.min_area, n_points=a.points,
              record_intra="all" if a.all_intra else "problems", attach=a.attach, hydrogens=a.hydrogens)
    ens_list = None
    all_models = None            # in ensemble mode the Ramachandran and rotamer pages and the PDF show every model of the file
    if a.models:
        a.ensemble = True
    if a.ensemble:
        from .ensemble import analyze_ensemble, count_models, parse_models
        try:
            sel = parse_models(a.models, count_models(a.structure))
        except ValueError as e:
            print(f"ifacekit: {e}", file=sys.stderr)
            return 2
        kw["keep"] = True
        all_models = list(range(count_models(a.structure)))
        try:
            if a.chains:
                tg = ([r["chain"] for r in list_chains(a.structure, 0) if r["has_polymer"]] if [c.lower() for c in a.chains] == ["all"]
                      else list(dict.fromkeys(a.chains)))
                ens_list = []
                for c in tg:
                    if len(tg) > 1 and not a.quiet:
                        print(f"checking chain {c} in every model ...", file=sys.stderr, flush=True)
                    ens_list.append(analyze_ensemble(a.structure, chain=c, workers=a.workers, **{k: v for k, v in kw.items() if k != "model_index"}))
            elif a.pair:
                ens_list = [analyze_ensemble(a.structure, pair=tuple(a.pair), workers=a.workers, **{k: v for k, v in kw.items() if k != "model_index"})]
            else:
                ens_list = [analyze_ensemble(a.structure, workers=a.workers, **{k: v for k, v in kw.items() if k != "model_index"})]
        except ValueError as e:
            print(f"ifacekit: {e}", file=sys.stderr)
            return 2
        if ens_list[0]["n_models_in_file"] == 1:
            print("ifacekit: the file has one model; the ensemble numbers are those of that model", file=sys.stderr)
        if a.ensemble_csv:
            from .ensemble import write_csv
            base, dot, ext = a.ensemble_csv.rpartition(".")
            for e in ens_list:
                write_csv(e, a.ensemble_csv if len(ens_list) == 1 else f"{base}_{e['target']}{dot}{ext}" if dot else f"{a.ensemble_csv}_{e['target']}")
            from .ensemble import write_residue_csv
            rb, rd, rx = a.ensemble_csv.rpartition(".")
            write_residue_csv(ens_list[0]["residue_rmsd"] if len(ens_list) == 1 else [r for e in ens_list for r in e["residue_rmsd"]], f"{rb}_residues.{rx}" if rd else f"{a.ensemble_csv}_residues")
            print(f"Ensemble table written to {a.ensemble_csv} (and the RMSD of each residue to the file with _residues)", file=sys.stderr)
        kw.pop("keep", None)
        for e in ens_list:
            e["_results"] = e.pop("model_results")
            e["selected_model"] = sel[0] + 1
        if a.html and (a.chains or a.pair):
            from .ensemblereport import build as ens_build
            try:
                files = ens_build(a.structure, ens_list, a.html, chains=tg if a.chains else None, pair=tuple(a.pair) if a.pair else None, indices=sel,
                                  separate=a.separate, hydrogens=a.hydrogens, workers=a.workers, molstar=a.molstar, molstar_url=a.molstar_url)
            except ValueError as e:
                print(f"ifacekit: {e}", file=sys.stderr)
                return 2
            print("Report written to " + ", ".join(files), file=sys.stderr)
            a.html = None
        files_asked = a.html or a.rama_html or a.rotamer_html or a.plots_pdf
        first_sel = sel[0]
        if not files_asked:
            a.model = 0
            for e in ens_list:
                e["first_model"] = e["_results"][e["models_analysed"].index(first_sel + 1)] if first_sel + 1 in e["models_analysed"] else e["first_model"]
            out = {"schema_version": ens_list[0]["schema_version"], "mode": "ensemble", "ensembles": [{k: v for k, v in e.items() if k != "_results"} for e in ens_list]}
            if a.output:
                with open(a.output, "w") as fh:
                    json.dump(out, fh, indent=1)
            if not a.quiet:
                from .ensemble import text_table
                sep = "\n\n" + "=" * 78 + "\n\n"
                parts = []
                for e in ens_list:
                    parts.append(text_table(e))
                    first = e["first_model"]
                    parts.append(f"Model {first_sel + 1}\n\n" + (report_quality(first, max_rows=10) if first["input"]["mode"] == "quality" else text_report(first)))
                print(sep.join(parts))
            return 0
        a.model = sel[0] if len(sel) == 1 else 0   # the pages below are for one model: the chosen one, or the first
        for e in ens_list:
            e.pop("_results", None)
    try:
        if a.chains:
            if [c.lower() for c in a.chains] == ["all"]:
                chains = [r["chain"] for r in list_chains(a.structure, a.model) if r["has_polymer"]]
            else:
                chains = list(dict.fromkeys(a.chains))
            runs = []
            if a.attach == "contact" and len(chains) > 1:
                kw["plan"] = plan_attachments(a.structure, a.model)
            for c in chains:
                if len(chains) > 1 and not a.quiet:
                    print(f"checking chain {c} ...", file=sys.stderr, flush=True)
                runs.append(analyze(a.structure, chain=c, **kw))
            out = {"schema_version": runs[0]["schema_version"], "mode": "quality", "runs": runs}
            pdf_ok = bool(a.plots_pdf) and _plots_pdf(a.structure, chains, a.plots_pdf, a.model, all_models)
            if a.html:
                from .qualityreport import build as quality_build
                quality_build(a.structure, chains, a.html, a.model, predicted=True if a.predicted else False if a.experimental else None,
                              pdf_name=_pdf_link(pdf_ok, a.plots_pdf, a.html), hydrogens=a.hydrogens, ensemble=a.ensemble)
                print(f"Quality report written to {a.html}", file=sys.stderr)
            if a.rotamer_html:
                from .rotamerhtml import build as rot_build
                rot_build(a.structure, a.rotamer_html, chains, a.model, title=a.structure.split("/")[-1], models=all_models)
                print(f"Rotamer page written to {a.rotamer_html}", file=sys.stderr)
            if a.rama_html:
                from .ramaplot import build as rama_build
                rama_build(a.structure, a.rama_html, chains, a.model, title=a.structure.split("/")[-1], models=all_models)
                print(f"Ramachandran page written to {a.rama_html}", file=sys.stderr)
            if len(chains) > 1:
                out["consistency"] = metal_consistency(a.structure, a.model)
        elif a.pair:
            runs = [analyze(a.structure, pair=tuple(a.pair), **kw)]
            out = runs[0]
            pdf_ok = bool(a.plots_pdf) and _plots_pdf(a.structure, list(a.pair), a.plots_pdf, a.model, all_models)
            if a.html:
                from .htmlreport import build, find_confidences
                # model type: the user's option wins, otherwise it is detected from the file (as for --ligand)
                exp = a.experimental or (not a.predicted and not runs[0]["flags"].get("predicted_model"))
                if not a.experimental and not a.predicted:
                    print(f"Model type: {'experimental' if exp else 'prediction'} ({runs[0]['flags'].get('model_source', '')}); "
                          "use --experimental or --predicted to override", file=sys.stderr)
                found = {} if exp else find_confidences(a.structure)
                summ = None if exp else (a.af3_summary or found.get("summary"))
                full = None if exp else (a.af3_full or found.get("full"))
                ms_path = (a.html[:-5] if a.html.lower().endswith(".html") else a.html) + "_3d.html" if a.molstar else None
                build(a.structure, a.pair[0], a.pair[1], a.html, summ, full, a.model, ms_path,
                      kind="experimental" if exp else "predicted", pdf_name=_pdf_link(pdf_ok, a.plots_pdf, a.html),
                      molstar_url=a.molstar_url, hydrogens=a.hydrogens, ensemble=a.ensemble)
                if ms_path:
                    print(f"3D report written to {ms_path}", file=sys.stderr)
                print(f"HTML report written to {a.html}"
                      + (f" (confidences: {summ or 'none'}, {full or 'no PAE'})"), file=sys.stderr)
        else:
            runs = [analyze(a.structure, **kw)]
            out = runs[0]
    except ValueError as e:
        print(f"ifacekit: {e}", file=sys.stderr)
        return 2

    if ens_list:
        out = {"schema_version": out["schema_version"], "mode": "ensemble", "ensembles": ens_list, "first_model_report_runs": out}
    if a.output:
        with open(a.output, "w") as fh:
            json.dump(out, fh, indent=1)
    if ens_list and not a.quiet:
        from .ensemble import text_table
        print("\n\n".join(text_table(e) for e in ens_list), end="\n\n" + "=" * 78 + "\n\n")
    if not a.quiet:
        sep = "\n\n" + "=" * 78 + "\n\n"
        if a.chains and len(runs) > 1:
            parts = [quality_overview(runs)]
            cons = consistency_text(out["consistency"])
            if cons:
                parts.append(cons)
            parts += [report_quality(r, max_rows=10) for r in runs]
            print(sep.join(parts))
        else:
            print(sep.join(text_report(r) for r in runs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
