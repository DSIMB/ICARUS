#!/usr/bin/env python3
"""
Summarise an ICARUS 2 proteome run (output of `icarus search`).

Usage: analyze_proteome.py flexible_alignments.tsv [top_n]
Reports throughput, score distributions and the pairs whose alignment gains
most from connected flexibility (tm_conn - tm_rigid_max), i.e. candidate
hinge motions / domain rearrangements between paralogues.
"""

import csv
import statistics as st
import sys


def main():
    rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
    top = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    rows = [r for r in rows if r["query"] != r["target"]]
    ms = [float(r["ms"]) for r in rows]
    print("pairs aligned: %d" % len(rows))
    print("per-pair time (1 thread): mean %.1f ms, median %.1f ms" % (st.mean(ms), st.median(ms)))
    for col in ("tm_rigid_max", "tm_conn", "tm_flex"):
        v = sorted(float(r[col]) for r in rows)
        print("%-13s median %.3f  >=0.5: %d pairs" % (col, v[len(v) // 2], sum(x >= 0.5 for x in v)))
    gain = []
    for r in rows:
        g = float(r["tm_conn"]) - float(r["tm_rigid_max"])
        if int(r["n_conn"]) >= 2:
            gain.append((g, r))
    gain.sort(key=lambda x: -x[0])
    print("pairs with connected flexible gain >= 0.1: %d" % sum(g >= 0.1 for g, _ in gain))
    print("\nlargest connected flexible gains:")
    print("%-24s %-24s %5s %5s %6s %6s %6s %s" % ("query", "target", "len_q", "len_t", "rigid", "conn", "gain", "bodies"))
    for g, r in gain[:top]:
        print("%-24s %-24s %5s %5s %6.3f %6.3f %6.3f %s" % (
            r["query"], r["target"], r["len_q"], r["len_t"], float(r["tm_rigid_max"]), float(r["tm_conn"]), g,
            r["bodies"]))


if __name__ == "__main__":
    main()
