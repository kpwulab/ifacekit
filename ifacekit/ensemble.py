"""Analysis of every model of a multi-model file (NMR ensembles, MD snapshots).

For one target (a chain, a pair of chains or the whole file) the analysis of ifacekit.core is run on each model. The result holds
  models       one row of numbers per model (all models)
  average      mean, standard deviation, minimum, maximum and the number of models for each number
  first_model  the complete result of the first model, as ifacekit gives it for a single model
The numbers are taken from the result of each model; nothing is re-computed here.
"""
import csv
import math

import numpy as np
from concurrent.futures import ProcessPoolExecutor

import gemmi

from . import core
from .attach import plan_attachments

# column name, label for the text table
COLUMNS = [
    ("rama_favored_pct", "Rama fav %"),
    ("rama_outlier_pct", "Rama out %"),
    ("rotamer_outlier_pct", "Rot out %"),
    ("rmsz_bonds", "RMSZ bond"),
    ("rmsz_angles", "RMSZ angle"),
    ("heavy_unfavorable", "unfav"),
    ("heavy_borderline", "border"),
    ("all_atom_clashes", "H clash"),
    ("all_atom_clashscore", "H score"),
    ("interface_area", "area A2"),
    ("hbonds_direct", "H bonds"),
]


def parse_models(spec, n):
    """Model numbers as the PDB file counts them (MODEL 1, 2, ...) -> indices counted from 0.
    spec: "first", "all", or numbers and ranges such as "3", "1,4,7", "2-5"."""
    spec = (spec or "first").strip().lower()
    if spec == "first":
        return [0]
    if spec == "all":
        return list(range(n))
    out = []
    for part in spec.split(","):
        part = part.strip()
        try:
            if "-" in part:
                a, b = part.split("-")
                nums = range(int(a), int(b) + 1)
            else:
                nums = [int(part)]
        except ValueError:
            raise ValueError(f"cannot read the model list '{spec}' (use first, all, or numbers like 3, 1,4,7 or 2-5)")
        for k in nums:
            if not 1 <= k <= n:
                raise ValueError(f"model {k} does not exist: the file has {n} models (numbered 1 to {n})")
            if k - 1 not in out:
                out.append(k - 1)
    if not out:
        raise ValueError("no model chosen")
    return out


def count_models(path):
    return len(gemmi.read_structure(path))


def metrics(res):
    """The numbers of one analysis result (one model) as a flat dict; entries that do not apply are left out."""
    m = {}
    s = res.get("summary", {})
    m["n_atoms"] = s.get("n_atoms")
    ra, ro, g = res.get("ramachandran") or {}, res.get("rotamer") or {}, res.get("geometry") or {}
    if "percent" in ra:
        m["rama_favored_pct"], m["rama_allowed_pct"], m["rama_outlier_pct"] = (ra["percent"][k] for k in ("favored", "allowed", "outlier"))
    if "percent" in ro:
        m["rotamer_favored_pct"], m["rotamer_outlier_pct"] = ro["percent"]["favored"], ro["percent"]["outlier"]
    for k in ("rmsz_bonds", "rmsz_angles", "n_outliers"):
        if g.get(k) is not None:
            m["geometry_outliers" if k == "n_outliers" else k] = g[k]
    v = s.get("verdicts", {})
    m["heavy_unfavorable"], m["heavy_borderline"] = v.get("unfavorable", 0), v.get("borderline", 0)
    if "all_atom_clashes" in s:
        m["all_atom_clashes"], m["all_atom_clashscore"] = s["all_atom_clashes"], s["all_atom_clashscore"]
    fl = (res.get("flips") or {}).get("counts")
    if fl:
        m["flips_suggested"] = fl.get("flip", 0)
    ifs = res.get("interfaces") or []
    if ifs:
        m["interface_area"] = round(sum(i["area"] for i in ifs), 1)
        m["hbonds_direct"] = sum((i.get("hbonds") or {}).get("direct", 0) for i in ifs)
        m["n_interfaces"] = len(ifs)
    return m


def _one(args):
    path, index, chain, pair, kw = args
    kw = dict(kw)
    if kw.get("attach", "contact") == "contact" and (chain is not None or pair is not None):
        kw["plan"] = plan_attachments(path, index)
    return core.analyze(path, model_index=index, chain=chain, pair=pair, **kw)


def summarize(rows):
    """Mean, standard deviation (n - 1), minimum, maximum and n of every number that occurs in the rows."""
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys and k != "model"]
    out = {}
    for k in keys:
        vals = [r[k] for r in rows if isinstance(r.get(k), (int, float))]
        if not vals:
            continue
        mean = sum(vals) / len(vals)
        sd = math.sqrt(sum((x - mean) ** 2 for x in vals) / (len(vals) - 1)) if len(vals) > 1 else 0.0
        out[k] = {"mean": round(mean, 3), "sd": round(sd, 3), "min": min(vals), "max": max(vals), "n": len(vals)}
    return out


def analyze_ensemble(path, chain=None, pair=None, models=None, workers=1, keep=False, **kw):
    """Run the analysis on every model of the file (or on the model indices in `models`, counted from 0).

    Returns {"mode": "ensemble", "target", "n_models_in_file", "models": [...], "average": {...}, "first_model": <full result>}.
    kw is passed to core.analyze (min_area, n_points, record_intra, attach, hydrogens, ...).
    keep: also return the complete result of every model as "model_results" (not for the JSON file).
    """
    n = count_models(path)
    idx = list(range(n)) if models is None else [i for i in models if 0 <= i < n]
    if not idx:
        raise ValueError(f"no model to analyse (the file has {n} models)")
    jobs = [(path, i, chain, pair, kw) for i in idx]
    if workers and workers > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(_one, jobs))
    else:
        results = [_one(j) for j in jobs]
    results[0].setdefault("flags", {})["ensemble"] = True
    rows = [{"model": i + 1, **metrics(r)} for i, r in zip(idx, results)]
    target = chain if chain is not None else "+".join(pair) if pair is not None else "all"
    rr = [r for r in residue_rmsd(path, models=idx) if chain is None and pair is None or r["chain"] in ([chain] if chain is not None else list(pair or []))]
    extra = {"model_results": results} if keep else {}
    return {**extra, "schema_version": results[0]["schema_version"], "mode": "ensemble", "target": target, "residue_rmsd": rr, "n_models_in_file": n, "models_analysed": [i + 1 for i in idx],
            "models": rows, "average": summarize(rows), "first_model": results[0]}


def write_csv(ens, path):
    """One row per model plus the mean and the standard deviation, for a spreadsheet."""
    keys = []
    for r in ens["models"]:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["target"] + keys)
        for r in ens["models"]:
            w.writerow([ens["target"]] + [r.get(k, "") for k in keys])
        for lab in ("mean", "sd", "min", "max"):
            w.writerow([ens["target"], lab] + [ens["average"].get(k, {}).get(lab, "") for k in keys[1:]])


def _fmt(x):
    if x is None or x == "":
        return "-"
    return f"{x:.2f}" if isinstance(x, float) and abs(x) < 100 else f"{x:.1f}" if isinstance(x, float) else str(x)


def text_table(ens):
    """Per-model table with the mean and standard deviation underneath."""
    cols = [(k, lab) for k, lab in COLUMNS if any(k in r for r in ens["models"])]
    head = ["model"] + [lab for _, lab in cols]
    lines = [f"Ensemble analysis of {ens['target']}: {len(ens['models'])} of {ens['n_models_in_file']} models", ""]
    rows = [[str(r["model"])] + [_fmt(r.get(k)) for k, _ in cols] for r in ens["models"]]
    av = ens["average"]
    rows.append(["mean"] + [_fmt(av.get(k, {}).get("mean")) for k, _ in cols])
    rows.append(["sd"] + [_fmt(av.get(k, {}).get("sd")) for k, _ in cols])
    width = [max(len(head[i]), *(len(r[i]) for r in rows)) for i in range(len(head))]
    fmt = lambda r: "  ".join(c.rjust(w) for c, w in zip(r, width))
    lines.append(fmt(head))
    lines += [fmt(r) for r in rows[:-2]]
    lines.append("-" * len(fmt(head)))
    lines += [fmt(r) for r in rows[-2:]]
    return "\n".join(lines)


BACKBONE = ("N", "CA", "C", "O")
CORE_FRACTION = 0.7        # share of the residues (lowest CA spread) that define the superposition frame
CORE_ROUNDS = 3


def _core(ca_by_model, n_res):
    """Indices of the residues used for the superposition: iterated choice of the least mobile 70% (at least 10 residues)."""
    core = np.arange(n_res)
    k = max(10, int(round(CORE_FRACTION * n_res)))
    for _ in range(CORE_ROUNDS if n_res > 10 else 0):
        ref = np.mean([_fit(m[core], ca_by_model[0][core], m) for m in ca_by_model], axis=0)
        dev = np.mean([np.sum((_fit(m[core], ref[core], m) - ref) ** 2, axis=1) for m in ca_by_model], axis=0)
        core = np.argsort(dev)[:min(k, n_res)]
    return np.sort(core)


def _fit(mobile_core, ref_core, mobile_all):
    """mobile_all moved onto the reference by the best rotation and translation of the core atoms (Kabsch)."""
    mc, rc = mobile_core.mean(axis=0), ref_core.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile_core - mc).T @ (ref_core - rc))
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    return (mobile_all - mc) @ rot + rc


def residue_rmsd(path, models=None):
    """Spread of every polymer residue over the models of a file (the precision of an NMR ensemble).

    Every chain is superposed on its own, on the CA atoms of its least mobile 70% of the residues (iterated three times), model 1 as the first reference.
    For each residue and each heavy atom the mean position over the models is taken; the RMSD of a residue is the root of the mean squared distance of its atoms to their mean positions,
    for the backbone atoms (N, CA, C, O), for the CA atom alone and for all heavy atoms. A residue has to be present with the same atoms in every model analysed.
    Besides the RMSD, every model has its own distance from the mean structure (root of the mean squared distance of the backbone atoms of the residue); dev_mean, dev_sd (n - 1), dev_min and dev_max
    are the mean, standard deviation and range of that distance over the models, and b (the value the plots draw as a line) is dev_mean.
    Returns rows in file order {chain, resnum, resname, b, dev_mean, dev_sd, dev_min, dev_max, rmsd_backbone, rmsd_ca, rmsd_all, n_models}.
    """
    import gemmi
    st = gemmi.read_structure(path)
    st.setup_entities()
    idx = list(range(len(st))) if models is None else [i for i in models if 0 <= i < len(st)]
    rows = []
    for ch in st[idx[0]]:
        recs = []        # per residue of the first model: key, name, {atom: index}
        for r in ch:
            if gemmi.find_tabulated_residue(r.name) is None or not gemmi.find_tabulated_residue(r.name).is_amino_acid() and not gemmi.find_tabulated_residue(r.name).is_nucleic_acid():
                continue
            recs.append(r)
        if not recs:
            continue
        per_model = []
        ok = True
        for i in idx:
            try:
                chm = st[i][ch.name]
            except Exception:
                ok = False
                break
            d = {}
            for r in chm:
                key = (r.seqid.num, r.seqid.icode.strip(), r.name)
                d[key] = {a.name: np.array(a.pos.tolist()) for a in r if a.element.name != "H"}
            per_model.append(d)
        if not ok:
            continue
        keys = [(r.seqid.num, r.seqid.icode.strip(), r.name) for r in recs if all((r.seqid.num, r.seqid.icode.strip(), r.name) in d for d in per_model)]
        keys = [k for k in keys if all("CA" in d[k] or "C1'" in d[k] for d in per_model)]
        if len(keys) < 3:
            continue
        cakey = lambda d, k: d[k]["CA"] if "CA" in d[k] else d[k]["C1'"]
        ca = [np.array([cakey(d, k) for k in keys]) for d in per_model]
        core = _core(ca, len(keys)) if len(keys) > 3 else np.arange(len(keys))
        ref = ca[0]
        for _ in range(2):
            ref = np.mean([_fit(m[core], ref[core], m) for m in ca], axis=0)
        # rotation of every model onto the mean frame, applied to all atoms
        moved = []
        for d, m in zip(per_model, ca):
            mc, rc = m[core].mean(axis=0), ref[core].mean(axis=0)
            u, _, vt = np.linalg.svd((m[core] - mc).T @ (ref[core] - rc))
            rot = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ vt))]) @ vt
            moved.append({k: {a: (x - mc) @ rot + rc for a, x in d[k].items()} for k in keys})
        for k in keys:
            names = [a for a in moved[0][k] if all(a in mv[k] for mv in moved)]

            def dev(sel):
                """Distance of each model from the mean structure for the atoms in sel: root of the mean squared distance of those atoms (one number per model)."""
                sel = [a for a in names if a in sel] if sel is not None else names
                if not sel:
                    return None
                pts = np.array([[mv[k][a] for a in sel] for mv in moved])          # models x atoms x 3
                return np.sqrt(np.mean(np.sum((pts - pts.mean(axis=0)) ** 2, axis=2), axis=1))

            def rms(d):
                return None if d is None else round(float(math.sqrt(np.mean(d ** 2))), 3)
            db, dc, da = dev(BACKBONE), dev(("CA",) if "CA" in names else ("C1'",)), dev(None)
            use = db if db is not None else da
            sd = float(np.std(use, ddof=1)) if len(use) > 1 else 0.0
            rows.append({"chain": ch.name, "resnum": k[0], "resname": k[2], "b": round(float(np.mean(use)), 3), "dev_mean": round(float(np.mean(use)), 3), "dev_sd": round(sd, 3),
                         "dev_min": round(float(np.min(use)), 3), "dev_max": round(float(np.max(use)), 3),
                         "rmsd_backbone": rms(db), "rmsd_ca": rms(dc), "rmsd_all": rms(da), "n_models": len(idx)})
    return rows


def write_residue_csv(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chain", "resnum", "resname", "rmsd_backbone", "rmsd_ca", "rmsd_all", "dev_mean", "dev_sd", "dev_min", "dev_max", "n_models"])
        for r in rows:
            w.writerow([r["chain"], r["resnum"], r["resname"], r["rmsd_backbone"], r["rmsd_ca"], r["rmsd_all"], r["dev_mean"], r["dev_sd"], r["dev_min"], r["dev_max"], r["n_models"]])
