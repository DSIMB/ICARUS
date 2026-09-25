#!/usr/bin/env python3
"""
SCOP lineage (class, fold, superfamily, family sunids) of the SCOP40 domain set
used by the Foldseek benchmark (scop40pdb.tar.gz), from RCSB annotations.

Usage: scop40_labels.py SCOP40_PDB_DIR OUT.tsv
"""

import concurrent.futures as cf
import json
import os
import sys
import time
import urllib.request

GRAPHQL = "https://data.rcsb.org/graphql"


def fetch(query, tries=6):
    for i in range(tries):
        try:
            req = urllib.request.Request(GRAPHQL, data=json.dumps({"query": query}).encode(),
                                         headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read())
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 ** i)


def batch(entries):
    q = """{ entries(entry_ids: [%s]) { rcsb_id polymer_entities { polymer_entity_instances {
             rcsb_polymer_instance_annotation { type annotation_id annotation_lineage { id depth } } } } } }""" % (
        ",".join('"%s"' % e.upper() for e in entries))
    out = {}
    d = fetch(q)
    for ent in d["data"]["entries"] or []:
        for pe in ent["polymer_entities"] or []:
            for inst in pe["polymer_entity_instances"] or []:
                for ann in inst["rcsb_polymer_instance_annotation"] or []:
                    if ann["type"] != "SCOP" or not ann["annotation_lineage"]:
                        continue
                    lin = {x["depth"]: x["id"] for x in ann["annotation_lineage"]}
                    out[ann["annotation_id"]] = [lin.get(k, "") for k in (1, 2, 3, 4)]
    return out


def main():
    sids = sorted(os.listdir(sys.argv[1]))
    entries = sorted({s[1:5] for s in sids})
    chunks = [entries[i:i + 100] for i in range(0, len(entries), 100)]
    labels = {}
    with cf.ThreadPoolExecutor(4) as ex:
        for res in ex.map(batch, chunks):
            labels.update(res)
    missing = 0
    with open(sys.argv[2], "w") as f:
        f.write("sid\tclass\tfold\tsuperfamily\tfamily\n")
        for s in sids:
            if s in labels:
                f.write("%s\t%s\n" % (s, "\t".join(labels[s])))
            else:
                missing += 1
    print("labelled %d / %d domains (%d missing)" % (len(sids) - missing, len(sids), missing))


if __name__ == "__main__":
    main()
