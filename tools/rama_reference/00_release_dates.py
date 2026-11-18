"""Fetch the initial release date of every listed entry (RCSB data API, GraphQL, batches of 500).
Usage: python 00_release_dates.py entries.txt OUT.tsv"""
import json
import sys
import urllib.request

ids = [l.strip() for l in open(sys.argv[1]) if l.strip()]
out = open(sys.argv[2], "w")
for k in range(0, len(ids), 500):
    q = {"query": "query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id rcsb_accession_info{initial_release_date}}}",
         "variables": {"ids": ids[k:k + 500]}}
    req = urllib.request.Request("https://data.rcsb.org/graphql", data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    for e in json.load(urllib.request.urlopen(req, timeout=120))["data"]["entries"]:
        if e:
            out.write(f"{e['rcsb_id'].lower()}\t{e['rcsb_accession_info']['initial_release_date'][:10]}\n")
out.close()
