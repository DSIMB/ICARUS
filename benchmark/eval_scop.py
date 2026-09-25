#!/usr/bin/env python3
"""
Sensitivity up to the first false positive on the sampled SCOP40 benchmark.

For each query, targets are ranked by a score; the fraction of family,
superfamily and fold level true positives ranked before the first false
positive (different fold) is averaged over queries (Foldseek/Reseek protocol).

Usage: eval_scop.py LABELS.tsv PAIRS.tsv RESULTS.tsv:COLUMN[:NAME] [...]
RESULTS.tsv needs 'query' and 'target' columns (or a header-less m8 file where
COLUMN is a 0-based index). Missing pairs score -inf.
"""

import csv
import sys
from collections import defaultdict


def level(lq, lt):
    if lq[1] != lt[1]:
        return "fp"
    if lq[2] != lt[2]:
        return "fold"
    if lq[3] != lt[3]:
        return "superfamily"
    return "family"


def strip_ext(name):
    for ext in (".gz", ".pdb", ".cif", ".ent"):
        if name.endswith(ext):
            name = name[: -len(ext)]
    return name


def load_scores(spec):
    parts = spec.split(":")
    path, col = parts[0], parts[1]
    name = parts[2] if len(parts) > 2 else path.split("/")[-1] + ":" + col
    scores = {}
    with open(path) as f:
        first = f.readline()
        f.seek(0)
        if first.startswith("query\t"):
            for r in csv.DictReader(f, delimiter="\t"):
                try:
                    scores[(r["query"], r["target"])] = float(r[col])
                except (ValueError, KeyError):
                    pass
        else:
            ci = int(col)
            for line in f:
                p = line.rstrip("\n").split("\t")
                q, t = strip_ext(p[0]), strip_ext(p[1])
                try:
                    v = float(p[ci])
                except ValueError:
                    continue
                if v > scores.get((q, t), float("-inf")):
                    scores[(q, t)] = v
    # symmetric scores (e.g. all-vs-all within one database) may be stored in
    # either orientation
    for (q, t), v in list(scores.items()):
        scores.setdefault((t, q), v)
    return name, scores


def main():
    labels = {}
    for line in open(sys.argv[1]):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sid":
            labels[f[0]] = tuple(f[1:5])
    targets = defaultdict(list)
    for line in open(sys.argv[2]):
        q, t = line.split()
        targets[q].append(t)
    print("%-28s %8s %12s %8s %8s" % ("method", "family", "superfamily", "fold", "pairs"))
    for spec in sys.argv[3:]:
        name, scores = load_scores(spec)
        sens = defaultdict(list)
        for q, ts in targets.items():
            ranked = sorted(ts, key=lambda t: -scores.get((q, t), float("-inf")))
            counts = defaultdict(int)
            totals = defaultdict(int)
            for t in ts:
                totals[level(labels[q], labels[t])] += 1
            for t in ranked:
                lv = level(labels[q], labels[t])
                if lv == "fp":
                    break
                counts[lv] += 1
            for lv in ("family", "superfamily", "fold"):
                if totals[lv]:
                    sens[lv].append(counts[lv] / totals[lv])
        found = sum(1 for q, ts in targets.items() for t in ts if (q, t) in scores)
        print("%-28s %8.3f %12.3f %8.3f %8d" % (
            name, *(sum(sens[l]) / len(sens[l]) if sens[l] else float("nan") for l in ("family", "superfamily", "fold")),
            found))


if __name__ == "__main__":
    main()
