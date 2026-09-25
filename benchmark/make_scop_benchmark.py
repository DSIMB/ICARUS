#!/usr/bin/env python3
"""
Build a sampled SCOP40 homology-discrimination benchmark.

Queries have at least one same-superfamily/different-family member. For each
query the targets are its same-fold domains (capped) and random decoys from
other folds. Levels (Foldseek convention): family = same family;
superfamily = same superfamily, other family; fold = same fold, other
superfamily; false positive = other fold.

Usage: make_scop_benchmark.py LABELS.tsv OUT_PAIRS.tsv [n_queries] [n_decoys] [seed]
"""

import random
import sys
from collections import defaultdict


def main():
    labels = {}
    for line in open(sys.argv[1]):
        f = line.rstrip("\n").split("\t")
        if f[0] == "sid" or len(f) < 5 or not all(f[1:5]):
            continue
        labels[f[0]] = tuple(f[1:5])
    nq = int(sys.argv[3]) if len(sys.argv) > 3 else 250
    nd = int(sys.argv[4]) if len(sys.argv) > 4 else 250
    rng = random.Random(int(sys.argv[5]) if len(sys.argv) > 5 else 1)
    by_fold, by_sf = defaultdict(list), defaultdict(list)
    for s, (c, fo, sf, fa) in labels.items():
        by_fold[fo].append(s)
        by_sf[sf].append(s)
    # classes a-d only (all-alpha, all-beta, alpha/beta, alpha+beta)
    main_classes = {"46456", "48724", "51349", "53931"}
    cands = [s for s, lab in labels.items() if lab[0] in main_classes
             and any(labels[o][3] != lab[3] for o in by_sf[lab[2]] if o != s)]
    queries = rng.sample(sorted(cands), min(nq, len(cands)))
    everyone = sorted(labels)
    with open(sys.argv[2], "w") as out:
        for q in queries:
            fo = labels[q][1]
            same = [s for s in by_fold[fo] if s != q]
            rng.shuffle(same)
            same = same[:60]
            decoys = []
            while len(decoys) < nd:
                s = rng.choice(everyone)
                if labels[s][1] != fo and s not in decoys:
                    decoys.append(s)
            for t in same + decoys:
                out.write("%s\t%s\n" % (q, t))
    print("queries:", len(queries))


if __name__ == "__main__":
    main()
