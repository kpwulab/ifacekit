"""ifacekit web portal (internal use): serves server/index.html and the API behind it.

Run:  python server/app.py [--host 127.0.0.1] [--port 8000] [--data DIR]

No archive is ever accepted. The user uploads a structure (.cif/.pdb) and, for AI predictions, the
two AlphaFold JSON files as separate files. Every upload is checked (size, plain text, expected
content, JSON parsed). The structure is stored under its own name reduced to letters, digits, dot, dash and
underscore (a PDB download as e.g. 8AFB.cif), with the extension set from the checked content, so that the reports show
what was analysed; the JSON files keep fixed names.

  POST /api/inspect   multipart: structure, [summary], [full]
        -> {job, chains:[{chain,type,residues,atoms,waters}], ligands:[{id,name,chain,atoms}],
            predicted, has_ai_files}
  POST /api/fetch     form: pdb_id (4 characters), source = rcsb | pdbe | pdbj | auto   -> same as inspect
  POST /api/analyze   multipart/form: job, type = predicted | interface | ligand | quality,
        chain_a, chain_b (optional for predicted/interface; quality: chains = "A,B" or empty for all),
        ligand (unit id from inspect), molstar = 0|1, hydrogens = 0|1 (interface and chain quality: also check clashes with hydrogens added), detailed = 0|1 (ligand only: also the detailed diagram), detailed_all = 0|1 (with detailed: every residue as a structure)
        -> {"report": url|null, "report_3d": url|null, "files":[{label,url}], "warnings":[...]}  or {"error":...}
  GET  /r/<job>/<run>/<file>   results
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
import threading
import urllib.error
import urllib.parse
import urllib.request
import uuid
import zlib

from flask import Flask, Response, abort, jsonify, request, send_from_directory

HERE = os.path.dirname(os.path.abspath(__file__))
MAX_STRUCT = 100 * 1024 * 1024
MAX_SUMMARY = 5 * 1024 * 1024
MAX_FULL = 100 * 1024 * 1024
TIMEOUT = 600
# Settings from the environment (a service file sets them); the command-line options of this file override DATA.
KEEP_DAYS = float(os.environ.get("IFACEKIT_KEEP_DAYS", "7"))              # job folders older than this are deleted
PURGE_MINUTES = float(os.environ.get("IFACEKIT_PURGE_MINUTES", "60"))     # how often the server looks for old jobs
MAX_JOBS = int(os.environ.get("IFACEKIT_MAX_JOBS", "4"))                  # analyses running at the same time (per server process)
QUEUE_WAIT = int(os.environ.get("IFACEKIT_QUEUE_WAIT", "600"))            # seconds a request waits for a free slot
MEM_MB = int(os.environ.get("IFACEKIT_MEM_MB", "8192"))                   # address-space limit per analysis process, 0 = none
CHAIN_RE = re.compile(r"^[A-Za-z0-9]{1,4}$")
JOB_RE = re.compile(r"^[0-9a-f]{32}$")
TYPES = ("predicted", "interface", "ligand", "quality")

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = MAX_FULL + MAX_STRUCT + MAX_SUMMARY + 1024 * 1024
DATA = os.environ.get("IFACEKIT_DATA", os.path.join(HERE, "jobs"))
_slots = threading.BoundedSemaphore(MAX_JOBS)
_state = {"running": 0, "waiting": 0, "purger": False}
_lock = threading.Lock()
MOLSTAR_LOCAL = os.path.join(HERE, "static", "molstar", "molstar.js")


class UserError(Exception):
    pass


def read_upload(field, limit, required=False):
    f = request.files.get(field)
    if f is None or not f.filename:
        if required:
            raise UserError("The %s file is missing." % field)
        return None, None
    raw = f.read(limit + 1)
    if len(raw) > limit:
        raise UserError("%s file is larger than %d MB." % (field, limit // 2**20))
    if not raw.strip():
        raise UserError("%s file is empty." % field)
    if b"\x00" in raw[:65536] or raw[:2] in (b"PK", b"\x1f\x8b", b"MZ", b"#!") or raw[:4] == b"\x7fELF":
        raise UserError("%s: only plain-text structure or JSON files are accepted (no archives or binaries)." % field)
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        raise UserError("%s is not a text file." % field)
    return f.filename, text


def check_structure(name, text):
    if re.search(r"^data_", text, re.M) and re.search(r"^(loop_|_atom_site\.)", text, re.M):
        if "_atom_site." not in text:
            raise UserError("No atom coordinates (_atom_site) found in the mmCIF file.")
        return "cif"
    if re.search(r"^(ATOM  |HETATM)", text, re.M):
        return "pdb"
    raise UserError("Not recognised as a PDB or mmCIF file.")


def check_json(field, text):
    try:
        d = json.loads(text)
    except ValueError:
        raise UserError("%s is not valid JSON." % field)
    if not isinstance(d, dict):
        raise UserError("%s is not an AlphaFold JSON object." % field)
    if field == "summary" and "iptm" not in d and "ranking_score" not in d:
        raise UserError("summary does not look like an AlphaFold summary_confidences file.")
    if field == "full" and "pae" not in d:
        raise UserError("full does not look like an AlphaFold full_data file (no 'pae').")


def _limited(cmd):
    """Start the command under a memory limit (ulimit -v, in kB), so that one huge structure cannot take the whole machine."""
    if MEM_MB > 0 and os.path.exists("/bin/sh"):
        return ["/bin/sh", "-c", 'ulimit -v "$0" 2>/dev/null; exec "$@"', str(MEM_MB * 1024)] + cmd
    return cmd


def run_limited(cmd, cwd, timeout):
    """Run a command in a free analysis slot (at most MAX_JOBS at the same time), with time and memory limits."""
    with _lock:
        _state["waiting"] += 1
    got = _slots.acquire(timeout=QUEUE_WAIT)
    with _lock:
        _state["waiting"] -= 1
        if got:
            _state["running"] += 1
    if not got:
        raise UserError("The server is busy with other analyses. Please try again in a few minutes.")
    try:
        return subprocess.run(_limited(cmd), cwd=cwd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise UserError("The analysis took longer than %d s and was stopped." % timeout)
    finally:
        with _lock:
            _state["running"] -= 1
        _slots.release()


OOM_WORDS = ("MemoryError", "Memory allocation", "Cannot allocate memory", "bad_alloc", "Killed")


def out_of_memory(r):
    return bool(r.returncode) and any(w in (r.stderr or "") for w in OOM_WORDS)


def run(args, cwd):
    r = run_limited([sys.executable, "-I", "-m", "ifacekit", "--factsheet-url", "/factsheet"] + args, cwd, TIMEOUT)
    if out_of_memory(r):
        raise UserError("The structure is too large for the memory limit of %d MB per analysis." % MEM_MB)
    return r


def last_error(r, default):
    msg = [l for l in r.stderr.splitlines() if l.startswith("ifacekit:")]
    return msg[-1][9:].strip() if msg else default


def job_dir(jid):
    if not JOB_RE.match(jid or "") or not os.path.isdir(os.path.join(DATA, jid)):
        raise UserError("Unknown job. Upload the file again.")
    return os.path.join(DATA, jid)


def safe_structure_name(name, fmt):
    """Own file name for the stored structure: letters, digits, dot, dash, underscore only; extension from the content."""
    base = os.path.basename(name or "")
    stem = "" if base.startswith(".") and base.count(".") == 1 else os.path.splitext(base)[0]    # ".cif" alone has no name
    stem = re.sub(r"[^A-Za-z0-9._-]", "_", stem).strip("._-")[:60] or "structure"
    return stem + "." + fmt


def structure_path(job):
    try:
        return os.path.join(job, json.load(open(os.path.join(job, "info.json")))["file"])
    except (OSError, ValueError, KeyError):
        pass
    for e in ("cif", "pdb"):                              # jobs made before the name was kept
        p = os.path.join(job, "input." + e)
        if os.path.exists(p):
            return p


@app.before_request
def _start_purger():
    """Start (once per server process) the thread that deletes old job folders, so that the 7 days of the footer hold also on a server that runs for months."""
    if not _state["purger"]:
        with _lock:
            if not _state["purger"]:
                _state["purger"] = True
                threading.Thread(target=_purge_loop, daemon=True).start()


def _purge_loop():
    while True:
        try:
            purge_old()
        except Exception:
            app.logger.exception("purge failed")
        time.sleep(max(PURGE_MINUTES, 1) * 60)


@app.route("/")
def home():
    from ifacekit import __version__
    from ifacekit.about import comparison_short_html
    with open(os.path.join(HERE, "index.html"), encoding="utf-8") as fh:
        page = fh.read().replace("__VERSION__", __version__).replace("__COMPARISON_SHORT__", comparison_short_html("/factsheet"))
    return Response(page, mimetype="text/html")


@app.route("/factsheet")
def factsheet():
    from ifacekit.about import factsheet_path
    return send_from_directory(os.path.dirname(factsheet_path()), "factsheet.html")


@app.route("/healthz")
def healthz():
    from ifacekit import __version__
    try:
        free = shutil.disk_usage(DATA if os.path.isdir(DATA) else HERE).free / 2**30
        jobs = sum(1 for d in os.listdir(DATA) if JOB_RE.match(d)) if os.path.isdir(DATA) else 0
    except OSError:
        free, jobs = -1, -1
    return jsonify({"status": "ok", "version": __version__, "running": _state["running"], "waiting": _state["waiting"],
                    "max_jobs": MAX_JOBS, "jobs_stored": jobs, "free_disk_gb": round(free, 1),
                    "local_molstar": os.path.exists(MOLSTAR_LOCAL)})


@app.route("/r/<jid>/<run_id>/<name>")
def result(jid, run_id, name):
    if not (JOB_RE.match(jid) and re.fullmatch(r"[0-9a-f]{8}", run_id)
            and re.fullmatch(r"[A-Za-z0-9._-]+\.(html|svg|png|pdf|json|txt|csv)", name)):
        abort(404)
    resp = send_from_directory(os.path.join(DATA, jid, "out", run_id), name)
    resp.headers["X-Content-Type-Options"] = "nosniff"
    return resp


def create_job(text, stext=None, ftext=None, name=None):
    """Check a structure (and optional AF JSON texts), store them, read chains/ligands. name: the structure's own name."""
    fmt = check_structure(None, text)
    if stext:
        check_json("summary", stext)
    if ftext:
        check_json("full", ftext)
    jid = uuid.uuid4().hex
    job = os.path.join(DATA, jid)
    os.makedirs(os.path.join(job, "out"))
    fname = safe_structure_name(name, fmt)
    for fn, t in ((fname, text), ("summary.json", stext), ("full.json", ftext)):
        if t:
            with open(os.path.join(job, fn), "w") as fh:
                fh.write(t)
    try:
        r = run_limited([sys.executable, "-I", os.path.join(HERE, "inspect_file.py"), os.path.join(job, fname)], job, 120)
    except UserError:
        shutil.rmtree(job, ignore_errors=True)
        raise
    if r.returncode:
        shutil.rmtree(job, ignore_errors=True)
        if out_of_memory(r):
            raise UserError("The structure is too large for the memory limit of %d MB per analysis." % MEM_MB)
        raise UserError("The structure could not be read: " + (r.stderr.strip().splitlines() or ["unknown error"])[-1][:200])
    info = json.loads(r.stdout)
    if not info["chains"]:
        shutil.rmtree(job, ignore_errors=True)
        raise UserError("No polymer chains found in the file.")
    info["has_ai_files"] = bool(stext and ftext)
    info["file"] = fname
    info["job"] = jid
    with open(os.path.join(job, "info.json"), "w") as fh:
        json.dump(info, fh)
    return info


@app.route("/api/inspect", methods=["POST"])
def inspect_api():
    try:
        name, text = read_upload("structure", MAX_STRUCT, required=True)
        sname, stext = read_upload("summary", MAX_SUMMARY)
        fname, ftext = read_upload("full", MAX_FULL)
        return jsonify(create_job(text, stext, ftext, name))
    except UserError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        app.logger.exception("inspect failed")
        return jsonify({"error": "Internal error."}), 500


# Fixed download addresses; only the validated 4-character ID is inserted. Nothing else from the
# client reaches the URL, and redirects are followed only inside the same list of hosts.
SOURCES = {
    "rcsb": ("RCSB PDB", "https://files.rcsb.org/download/{U}.cif", False),
    "pdbe": ("PDBe", "https://www.ebi.ac.uk/pdbe/entry-files/download/{L}.cif", False),
    "pdbj": ("PDBj", "https://pdbj.org/rest/newweb/fetch/file?format=mmcif&id={L}", True),   # gzip
}
HOSTS = {"files.rcsb.org", "www.ebi.ac.uk", "pdbj.org", "data.pdbj.org"}
PDB_ID_RE = re.compile(r"[0-9][A-Za-z0-9]{3}")


class _SameHosts(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        p = urllib.parse.urlparse(newurl)
        if p.scheme != "https" or p.hostname not in HOSTS:
            raise UserError("The server redirected to an address that is not allowed.")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def download(url, gz):
    opener = urllib.request.build_opener(_SameHosts)
    req = urllib.request.Request(url, headers={"User-Agent": "ifacekit", "Accept-Encoding": "identity"})
    with opener.open(req, timeout=30) as r:
        raw = r.read(MAX_STRUCT + 1)
    if len(raw) > MAX_STRUCT:
        raise UserError("The downloaded file is larger than %d MB." % (MAX_STRUCT // 2**20))
    if raw[:2] == b"\x1f\x8b":
        if not gz:
            raise UserError("Unexpected compressed download.")
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        raw = d.decompress(raw, MAX_STRUCT + 1)
        if len(raw) > MAX_STRUCT or d.unconsumed_tail:
            raise UserError("The downloaded file is larger than %d MB when unpacked." % (MAX_STRUCT // 2**20))
    elif raw[:2] == b"PK" or b"\x00" in raw[:65536]:
        raise UserError("The download is not a text file.")
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        raise UserError("The download is not a text file.")


@app.route("/api/fetch", methods=["POST"])
def fetch_api():
    try:
        pid = request.form.get("pdb_id", "").strip()
        if not PDB_ID_RE.fullmatch(pid):
            raise UserError("A PDB ID has 4 characters, for example 8AFB.")
        src = request.form.get("source", "rcsb")
        order = list(SOURCES) if src == "auto" else [src]
        if any(o not in SOURCES for o in order):
            raise UserError("Unknown source.")
        tried = []
        for o in order:
            name, tmpl, gz = SOURCES[o]
            url = tmpl.replace("{U}", pid.upper()).replace("{L}", pid.lower())
            try:
                text = download(url, gz)
                check_structure(None, text)
            except urllib.error.HTTPError as e:
                tried.append("%s: %s" % (name, "entry not found" if e.code == 404 else "HTTP %d" % e.code))
                continue
            except (urllib.error.URLError, OSError, zlib.error) as e:
                tried.append("%s: not reachable" % name)
                continue
            except UserError as e:
                tried.append("%s: %s" % (name, e))
                continue
            info = create_job(text, name=pid.upper())
            info["source"], info["pdb_id"] = name, pid.upper()
            return jsonify(info)
        raise UserError("Could not get %s. %s" % (pid.upper(), "; ".join(tried)))
    except UserError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        app.logger.exception("fetch failed")
        return jsonify({"error": "Internal error."}), 500


def pick_pair(cif, job, chains):
    out = os.path.join(job, "auto.json")
    r = run([cif, "--all", "-q", "-o", out], job)
    if r.returncode or not os.path.exists(out):
        raise UserError(last_error(r, "Could not analyse the structure."))
    ifs = [i for i in json.load(open(out)).get("interfaces", [])
           if i.get("kind_a") == "polymer" and i.get("kind_b") == "polymer"]
    if not ifs:
        raise UserError("No chain-chain interface found. The file may have one chain: use ligand or chain quality.")
    b = max(ifs, key=lambda i: i["area"])
    return b["unit_a"], b["unit_b"]


def ensemble_files(out, base):
    """Download links of the CSV tables of an ensemble analysis (one table per chain when several chains were chosen)."""
    res = []
    for n in sorted(os.listdir(out)):
        if n == "ensemble_residues.csv":
            res.append({"label": "RMSD of each residue over the models (CSV)", "url": base + n})
        elif re.fullmatch(r"ensemble(_[A-Za-z0-9]+)?\.csv", n):
            res.append({"label": "Ensemble table, one row per model (CSV)" + ("" if n == "ensemble.csv" else " " + n[9:-4]), "url": base + n})
    return res


@app.route("/api/analyze", methods=["POST"])
def analyze_api():
    try:
        job = job_dir(request.form.get("job"))
        info = json.load(open(os.path.join(job, "info.json")))
        cif = structure_path(job)
        typ = request.form.get("type", "")
        if typ not in TYPES:
            raise UserError("Unknown analysis type.")
        names = [c["chain"] for c in info["chains"]]
        ms = request.form.get("molstar") == "1"
        rid = uuid.uuid4().hex[:8]
        out = os.path.join(job, "out", rid)
        os.makedirs(out)
        base = "/r/%s/%s/" % (os.path.basename(job), rid)
        warnings, files = [], []
        resp = {"report": None, "report_3d": None}

        if typ in ("predicted", "interface"):
            ca, cb = request.form.get("chain_a", "").strip(), request.form.get("chain_b", "").strip()
            if bool(ca) != bool(cb):
                raise UserError("Give both chains, or neither for the largest interface.")
            if ca:
                if ca == cb or ca not in names or cb not in names:
                    raise UserError("Chains must be two different chains of the file: " + ", ".join(names))
            else:
                if len(names) < 2:
                    raise UserError("The file has one chain, so there is no interface. Use ligand or chain quality.")
                ca, cb = pick_pair(cif, job, names)
                warnings.append("Chain pair chosen automatically (largest interface): %s and %s." % (ca, cb))
            args = [cif, "--pair", ca, cb, "--html", os.path.join(out, "report.html"), "--plots-pdf", os.path.join(out, "plots.pdf")]
            if typ == "interface":
                args.append("--experimental")
            else:
                if os.path.exists(os.path.join(job, "summary.json")):
                    args += ["--af3-summary", os.path.join(job, "summary.json")]
                if os.path.exists(os.path.join(job, "full.json")):
                    args += ["--af3-full", os.path.join(job, "full.json")]
                if not info["has_ai_files"]:
                    warnings.append("Both AlphaFold JSON files were not given: no AI metrics or PAE plot.")
            if request.form.get("hydrogens") == "1":
                args.append("--hydrogens")
            ens = typ == "interface" and request.form.get("ensemble") == "1" and info.get("n_models", 1) > 1
            if ens:
                args += ["--ensemble", "--models", "all", "--workers", "2", "--ensemble-csv", os.path.join(out, "ensemble.csv")]
            if ms:
                args.append("--molstar")
                if os.path.exists(MOLSTAR_LOCAL):
                    args += ["--molstar-url", "/static/molstar/"]      # served by this portal, so the 3D page works without internet
            r = run(args, job)
            if r.returncode or not os.path.exists(os.path.join(out, "report.html")):
                raise UserError(last_error(r, "The analysis failed."))
            resp["report"] = base + "report.html"
            if os.path.exists(os.path.join(out, "report_3d.html")):
                resp["report_3d"] = base + "report_3d.html"
            if os.path.exists(os.path.join(out, "plots.pdf")):
                files.append({"label": "Ramachandran and rotamer plots (PDF)", "url": base + "plots.pdf"})
            if ens:
                files += ensemble_files(out, base)

        elif typ == "ligand":
            lig = request.form.get("ligand", "")
            if lig not in [l["id"] for l in info["ligands"]]:
                raise UserError("Choose a ligand from the list.")
            args = [cif, "--ligand", lig, "--html", os.path.join(out, "report.html"),
                    "--predicted" if info["predicted"] else "--experimental"]
            if request.form.get("detailed") == "1":
                args.append("--detailed")
                if request.form.get("detailed_all") == "1":
                    args.append("--detailed-all")
            if ms:
                args.append("--molstar")
                if os.path.exists(MOLSTAR_LOCAL):
                    args += ["--molstar-url", "/static/molstar/"]      # served by this portal, so the 3D page works without internet
            r = run(args, job)
            if r.returncode or not os.path.exists(os.path.join(out, "report.html")):
                raise UserError(last_error(r, "The ligand analysis failed."))
            resp["report"] = base + "report.html"
            if os.path.exists(os.path.join(out, "report_3d.html")):
                resp["report_3d"] = base + "report_3d.html"
            for n in sorted(os.listdir(out)):
                m = re.fullmatch(r"report_map(_all)?_(modern|classic|detailed)\.(png|svg|pdf)", n)
                if m:
                    files.append({"label": "Interaction map, %s%s (%s)" % ("all contacts, " if m.group(1) else "", m.group(2), m.group(3)), "url": base + n})

        else:  # quality
            sel = [c.strip() for c in request.form.get("chains", "").split(",") if c.strip()] or ["all"]
            if sel != ["all"] and not all(c in names for c in sel):
                raise UserError("Chains must be in the file: " + ", ".join(names))
            ens = request.form.get("ensemble") == "1" and info.get("n_models", 1) > 1
            ens_args = (["--ensemble", "--models", "all", "--workers", "2", "--ensemble-csv", os.path.join(out, "ensemble.csv"), "--molstar"]
                        + (["--molstar-url", "/static/molstar/"] if os.path.exists(MOLSTAR_LOCAL) else [])) if ens else []
            r = run([cif, "--chains"] + sel + ens_args + ["-o", os.path.join(out, "quality.json"),
                      "--html", os.path.join(out, "report.html"),
                      "--rama-html", os.path.join(out, "ramachandran.html"),
                      "--rotamer-html", os.path.join(out, "rotamer.html"), "--plots-pdf", os.path.join(out, "plots.pdf")]
                      + (["--hydrogens"] if request.form.get("hydrogens") == "1" else []), job)
            if r.returncode:
                raise UserError(last_error(r, "The analysis failed."))
            with open(os.path.join(out, "quality.txt"), "w") as fh:
                fh.write(r.stdout)
            resp["report"] = base + "report.html" if os.path.exists(os.path.join(out, "report.html")) else base + "quality.txt"
            if os.path.exists(os.path.join(out, "report_3d.html")):
                resp["report_3d"] = base + "report_3d.html"
            files.append({"label": "Text report", "url": base + "quality.txt"})
            if os.path.exists(os.path.join(out, "ramachandran.html")):
                files.append({"label": "Ramachandran plots", "url": base + "ramachandran.html"})
            if os.path.exists(os.path.join(out, "plots.pdf")):
                files.append({"label": "Ramachandran and rotamer plots (PDF)", "url": base + "plots.pdf"})
            if os.path.exists(os.path.join(out, "rotamer.html")):
                files.append({"label": "Rotamer plots", "url": base + "rotamer.html"})
            if ens:
                files += ensemble_files(out, base)
            files.append({"label": "Full result (JSON)", "url": base + "quality.json"})

        resp["files"], resp["warnings"] = files, warnings
        return jsonify(resp)
    except UserError as e:
        return jsonify({"error": str(e)}), 400
    except Exception:
        app.logger.exception("analysis failed")
        return jsonify({"error": "Internal error."}), 500


@app.errorhandler(413)
def too_big(_):
    return jsonify({"error": "Upload is too large."}), 413


def purge_old():
    now = time.time()
    if not os.path.isdir(DATA):
        return
    for d in os.listdir(DATA):
        p = os.path.join(DATA, d)
        if JOB_RE.match(d) and now - os.path.getmtime(p) > KEEP_DAYS * 86400:
            shutil.rmtree(p, ignore_errors=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--data", default=DATA)
    a = ap.parse_args()
    DATA = os.path.abspath(a.data)
    os.makedirs(DATA, exist_ok=True)
    purge_old()
    app.run(host=a.host, port=a.port, threaded=True)
