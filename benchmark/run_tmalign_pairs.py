#!/usr/bin/env python3
"""
Run TM-align (or US-align) on a list of pairs in parallel and write a TSV with
the TM-scores normalised by each chain.

Usage: run_tmalign_pairs.py EXE PAIRS.tsv STRUCT_DIR OUT.tsv [threads]
"""

import concurrent.futures as cf
import os
import re
import subprocess
import sys


def run(args):
    exe, d, q, t = args
    out = subprocess.run([exe, os.path.join(d, q), os.path.join(d, t)], capture_output=True, text=True).stdout
    tms = re.findall(r"TM-score= ([0-9.]+) \(if normalized by length of Chain_([12])", out)
    s = {c: float(v) for v, c in tms}
    return q, t, s.get("1", float("nan")), s.get("2", float("nan"))


def main():
    exe, pairs, d, out = sys.argv[1:5]
    threads = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    todo = [(exe, d, *l.split()[:2]) for l in open(pairs) if l.strip()]
    with cf.ThreadPoolExecutor(threads) as ex, open(out, "w") as f:
        f.write("query\ttarget\ttm_q\ttm_t\ttm_max\n")
        for q, t, a, b in ex.map(run, todo, chunksize=64):
            f.write("%s\t%s\t%.4f\t%.4f\t%.4f\n" % (q, t, a, b, min(a, b)))


if __name__ == "__main__":
    main()
