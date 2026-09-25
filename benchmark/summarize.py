#!/usr/bin/env python3
"""
Print markdown tables summarising run_benchmark.py outputs.

Usage: summarize.py OUT_DIR label=file.tsv [label=file.tsv ...]
"""

import csv
import os
import statistics as st
import sys


def main():
    out_dir = sys.argv[1]
    print("| method | pairs | mean TM-score | median | ref. agreement | mean s/pair |")
    print("|---|---|---|---|---|---|")
    for spec in sys.argv[2:]:
        label, fname = spec.split("=", 1)
        rows = list(csv.DictReader(open(os.path.join(out_dir, fname)), delimiter="\t"))
        tms = [float(r["tm"]) for r in rows if r["tm"] != "NA"]
        ag = [float(r["ref_agreement"]) for r in rows if r["ref_agreement"] not in ("NA", "")]
        secs = [float(r["seconds"]) for r in rows if r["seconds"] not in ("NA", "")]
        print("| %s | %d | %.3f | %.3f | %s | %s |" % (
            label, len(tms), st.mean(tms), st.median(tms),
            "%.1f%%" % (100 * st.mean(ag)) if ag else "-",
            "%.3f" % st.mean(secs) if secs else "-"))


if __name__ == "__main__":
    main()
