#!/usr/bin/env python3
"""
Run the original (Python/KPAX) ICARUS v1 on a list of structure pairs and record
its score, wall-clock runtime and the path of the best solution PDB.

Usage: run_icarus_v1.py PAIRS.tsv STRUCT_DIR OUT_DIR [--level L] [--cpu N]
PAIRS.tsv: tab-separated, two columns with structure ids (file STRUCT_DIR/<id>.pdb)
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import time

ICARUS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "icarus.py")


def first_chain(path):
    for line in open(path):
        if line.startswith("ATOM"):
            return line[21] if line[21] != " " else "A"
    return "A"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pairs")
    ap.add_argument("struct_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--level", type=int, default=2)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--cols", default="0,1", help="columns of the two ids in PAIRS")
    a = ap.parse_args()
    c1, c2 = (int(x) for x in a.cols.split(","))
    os.makedirs(a.out_dir, exist_ok=True)
    res_path = os.path.join(a.out_dir, "results_level%d.tsv" % a.level)
    done = set()
    if os.path.exists(res_path):
        for line in open(res_path):
            f = line.split("\t")
            if len(f) > 2 and f[2] != "NA":
                done.add(tuple(f[:2]))
    pairs = []
    for line in open(a.pairs):
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(c1, c2) or p[c1] in ("Domain1", "Mol1"):
            continue
        pairs.append((p[c1], p[c2]))
    for q, t in pairs:
        if (q, t) in done:
            continue
        wd = os.path.join(a.out_dir, "L%d_%s_%s" % (a.level, q, t))
        shutil.rmtree(wd, ignore_errors=True)
        os.makedirs(wd)
        for s in (q, t):
            shutil.copy(os.path.join(a.struct_dir, s + ".pdb"), os.path.join(wd, s + ".pdb"))
        t0 = time.time()
        proc = subprocess.run(["python3", ICARUS, "-p1", q + ".pdb", "-p2", t + ".pdb", "-l", str(a.level),
                               "-c1", first_chain(os.path.join(wd, q + ".pdb")),
                               "-c2", first_chain(os.path.join(wd, t + ".pdb")),
                               "-c", str(a.cpu), "-f"], cwd=wd, capture_output=True, text=True)
        dt = time.time() - t0
        with open(os.path.join(wd, "stdout.txt"), "w") as f:
            f.write(proc.stdout + "\n-----STDERR-----\n" + proc.stderr)
        score = re.findall(r"Score: (?:\x1b\[92m)?([0-9.]+)", proc.stdout)
        sols = sorted(glob.glob(os.path.join(wd, "icarus_output", "results", "*", "solution_*-on-*.pdb")))
        sols = [s for s in sols if not s.endswith("_renum.pdb")]
        with open(res_path, "a") as f:
            f.write("%s\t%s\t%s\t%.2f\t%s\n" % (q, t, score[0] if score else "NA", dt, sols[0] if sols else "NA"))
        print(q, t, score[0] if score else "NA", "%.1fs" % dt, flush=True)


if __name__ == "__main__":
    main()
