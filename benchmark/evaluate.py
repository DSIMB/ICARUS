#!/usr/bin/env python3
"""
Uniform evaluation of structural alignments.

Every method is scored the same way from its superposed coordinates:
the moved structure (for flexible methods: the chimera of transformed rigid
bodies, ordered along the target) and the fixed structure are aligned by the
gdt2.pl dynamic programme (`icarus gdt`), and the TM-score is normalised by the
length of the shorter of the two *original* structures. Agreement with the
reference alignments of Mayr et al. (2007) is the fraction of reference residue
pairs reproduced by that alignment.

Usage (as a module): see run_benchmark.py
"""

import os
import re
import subprocess
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ICARUS_BIN = os.environ.get("ICARUS_BIN", os.path.join(os.path.dirname(HERE), "icarus-rs", "target", "release", "icarus"))


def n_residues(path):
    """Number of residues with a CA atom (first chain, first model)."""
    seen = set()
    chain = None
    for line in open(path):
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")) and line[12:16] == " CA ":
            if chain is None:
                chain = line[21]
            if line[21] != chain:
                continue
            seen.add(line[22:27])
    return len(seen)


def gdt(moved, fixed, lnorm):
    """Returns (tm, pairs) where pairs are (resid_moved, resid_fixed, dist)."""
    out = subprocess.run([ICARUS_BIN, "gdt", moved, fixed, "--len", str(lnorm), "--pairs"],
                         capture_output=True, text=True, check=True).stdout.split("\n")
    tm = float(out[0].split("\t")[0])
    pairs = []
    for line in out[1:]:
        f = line.split("\t")
        if len(f) == 3:
            pairs.append((f[0], f[1], float(f[2])))
    return tm, pairs


def split_models(path, out1, out2):
    """Split a two-model PDB (blocks separated by END/ENDMDL/TER+HEADER) in two files."""
    blocks = [[]]
    for line in open(path):
        if line.startswith(("END", "ENDMDL")) and not line.startswith("ENDROOT"):
            if blocks[-1]:
                blocks.append([])
            continue
        if line.startswith(("ATOM", "HETATM")):
            blocks[-1].append(line)
    blocks = [b for b in blocks if b]
    with open(out1, "w") as f:
        f.writelines(blocks[0])
    with open(out2, "w") as f:
        f.writelines(blocks[1] if len(blocks) > 1 else [])


def read_reference(path):
    """Mayr et al. reference alignments: {(dom1, dom2): [((num, icode), (num, icode)), ...]}"""
    refs = defaultdict(list)
    key = None
    for line in open(path):
        line = line.strip()
        if line.startswith("#d"):
            a, b = line[1:].split("-")
            key = (a, b)
        elif line and not line.startswith("#") and key:
            f = line.split("\t")
            if len(f) == 2:
                def parse(x):
                    p = x.split(".")
                    return (int(p[1]), "" if p[2] == "_" else p[2])
                refs[key].append((parse(f[0]), parse(f[1])))
    return refs


def resid_key(s):
    m = re.match(r"(-?\d+)([A-Za-z]?)$", s)
    return (int(m.group(1)), m.group(2))


def residue_names(path):
    """{(num, icode): resname} of the first chain."""
    names = {}
    chain = None
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) and line[12:16] == " CA ":
            if chain is None:
                chain = line[21]
            if line[21] != chain:
                continue
            names[(int(line[22:26]), line[26].strip())] = line[17:20]
    return names


def reference_offset(ref_pairs, names1, names2):
    """Residue-number offsets (o1, o2) maximising residue-name agreement with the
    reference (handles entries renumbered since 2007, e.g. 1JWY)."""
    best = []
    for side, names in ((0, names1), (1, names2)):
        keys = [p[side] for p in ref_pairs]
        scores = {}
        for off in range(-1000, 1001):
            ok = sum(1 for (n, ic) in keys if (n + off, ic) in names)
            if ok:
                scores[off] = ok
        best.append(max(scores, key=lambda o: (scores[o], -abs(o))) if scores else 0)
    return best


def agreement(ref_pairs, aln_pairs, swap, off=(0, 0)):
    """Fraction of reference pairs present in the alignment. aln_pairs are
    (resid_moved, resid_fixed); swap=True if the moved structure is the second
    domain of the reference pair."""
    got = set()
    for a, b, _ in aln_pairs:
        ka, kb = resid_key(a), resid_key(b)
        got.add((kb, ka) if swap else (ka, kb))
    hit = 0
    for (n1, i1), (n2, i2) in ref_pairs:
        if ((n1 + off[0], i1), (n2 + off[1], i2)) in got:
            hit += 1
    return hit / len(ref_pairs) if ref_pairs else float("nan")
