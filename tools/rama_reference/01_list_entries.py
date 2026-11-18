"""Step 1: list candidate PDB entries for the Ramachandran reference (RCSB Search API, no structure download).

Criteria: X-ray, resolution <= RES, R-free <= RFREE, protein polymer entities, one representative per
30% sequence-identity group. Output: entries.txt (one PDB ID per line) and entries_summary.json.
Usage: python 01_list_entries.py OUTDIR [RES=1.8] [RFREE=0.25]
"""
import json
import sys
import urllib.request

out = sys.argv[1]
RES = float(sys.argv[2]) if len(sys.argv) > 2 else 1.8
RFREE = float(sys.argv[3]) if len(sys.argv) > 3 else 0.25
URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def query(start, rows):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal", "value": RES}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "refine.ls_R_factor_R_free", "operator": "less_or_equal", "value": RFREE}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "entity_poly.rcsb_entity_polymer_type", "operator": "exact_match", "value": "Protein"}}]},
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": start, "rows": rows},
                            "group_by": {"aggregation_method": "sequence_identity", "similarity_cutoff": 30, "ranking_criteria_type": {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}},
                            "group_by_return_type": "representatives"}}
    req = urllib.request.Request(URL, data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


ids, start, total = [], 0, None
while total is None or start < total:
    d = query(start, 10000)
    total = d["group_by_count"] + d.get("ungrouped_count", 0)   # representatives: one per group plus singletons
    ids += [x["identifier"].split("_")[0] for x in d["result_set"]]
    start += 10000
    print(f"{len(ids)} / {total}", file=sys.stderr)
ids = sorted(set(ids))
open(f"{out}/entries.txt", "w").write("\n".join(ids) + "\n")
json.dump({"resolution_max": RES, "rfree_max": RFREE, "groups": total, "entries": len(ids)}, open(f"{out}/entries_summary.json", "w"), indent=1)
print(f"{len(ids)} entries written to {out}/entries.txt")
