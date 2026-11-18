"""Public PDB entries for the tests: downloaded from RCSB on first use and kept in a cache folder, so the repository carries no structure files.

need("4QZS.cif") or need("1LCD.pdb") gives the path of the entry; a test is skipped when the entry cannot be fetched (no internet).
The cache is tests/.pdb_cache, or the folder in the environment variable IFACEKIT_TEST_CACHE (point it to a folder that holds the files to run offline).
"""
import os
import time
import urllib.request
from pathlib import Path

import pytest

CACHE = Path(os.environ.get("IFACEKIT_TEST_CACHE", Path(__file__).resolve().parent / ".pdb_cache"))
URL = "https://files.rcsb.org/download/{name}"


def _valid(path, ext):
    head = path.read_bytes()[:20000].decode("latin-1")
    return head.startswith("data_") if ext == "cif" else ("ATOM" in head or "HEADER" in head)


def need(name):
    """Path of a public PDB entry, for example "4QZS.cif" or "1LCD.pdb"."""
    code, ext = name.rsplit(".", 1)
    fname = f"{code.upper()}.{ext.lower()}"
    path = CACHE / fname
    if path.exists() and path.stat().st_size > 1000 and _valid(path, ext.lower()):
        return str(path)
    CACHE.mkdir(parents=True, exist_ok=True)
    err = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(URL.format(name=fname), headers={"User-Agent": "ifacekit-tests"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            tmp = path.with_suffix(path.suffix + ".part")
            tmp.write_bytes(data)
            if len(data) > 1000 and _valid(tmp, ext.lower()):
                tmp.replace(path)
                return str(path)
            tmp.unlink()
            err = "the file is not a structure"
        except Exception as e:                      # offline, or the server is busy
            err = str(e)
        time.sleep(2 * (attempt + 1))
    pytest.skip(f"{fname} could not be downloaded from RCSB ({err})")


def _without(src_name, out_name, cats, marker=False):
    import gemmi
    path = CACHE / out_name
    if path.exists():
        return str(path)
    doc = gemmi.cif.read_file(need(src_name))
    block = doc.sole_block()
    for cat in cats:
        table = block.find_mmcif_category(cat)
        if table:
            table.erase()
    if marker:
        block.set_pair("_ma_qa_metric.id", "1")
    doc.write_file(str(path))
    return str(path)


def no_dictionary():
    """8JUX without the bond table of its ligand V00, so the bond orders have to be rebuilt from the geometry."""
    return _without("8JUX.cif", "8JUX_no_ligand_dictionary.cif", ("_chem_comp_bond.",))


def prediction_like():
    """4QZS changed to look like an AI prediction: no bond table for the ligand, no experimental records, and the marker of the model archive format."""
    return _without("4QZS.cif", "4QZS_prediction_like.cif", ("_chem_comp_bond.", "_exptl.", "_refine.", "_reflns.", "_refine_hist.", "_diffrn.", "_diffrn_source."), marker=True)
