"""Run in a separate process: print JSON with chains, ligands and model type of one structure file."""
import json
import re
import sys

from ifacekit import model as M
from ifacekit.core import list_chains

path = sys.argv[1]
m = M.load(path)
chains = [{"chain": r["chain"], "type": r["polymer_type"], "residues": r["n_polymer_residues"],
           "atoms": r["n_atoms"], "waters": r["n_waters"]} for r in list_chains(path) if r["has_polymer"]]
ligs = [{"id": u["id"], "name": re.sub(r"\d+$", "", u["id"].split(":")[-1]),
         "chain": u.get("chain"), "atoms": u["n_atoms"]}
        for u in m.units if u["kind"] == "ligand"]
print(json.dumps({"chains": chains, "ligands": ligs, "predicted": bool(m.flags.get("predicted_model")), "n_models": int(m.flags.get("n_models_in_file", 1))}))
