#!/usr/bin/env python3
"""
Run and score structural aligners on a pair set (RIPC or SISY).

  run_benchmark.py METHOD --pairs PAIRS.tsv --cols 1,2 --structs DIR --out OUT_DIR [--refs REF.txt]

METHOD: icarus2 | icarus_v1 | tmalign | usalign
Each method produces a superposed "moved" structure per pair; all are scored with
evaluate.gdt (gdt2.pl-equivalent), normalised by the shorter original structure.
Writes OUT_DIR/<METHOD>[<tag>].tsv with one line per pair.
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time

import evaluate as ev

HERE = os.path.dirname(os.path.abspath(__file__))


def read_pairs(path, cols):
    c1, c2 = cols
    out = []
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(c1, c2) or p[c1] in ("Domain1", "Mol1"):
            continue
        out.append((p[c1], p[c2]))
    return out


def apply_matrix_pdb(src, dst, rot, trans):
    """Write src transformed by x' = t + U x (TM-align -m convention)."""
    with open(src) as f, open(dst, "w") as g:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                nx = trans[0] + rot[0][0] * x + rot[0][1] * y + rot[0][2] * z
                ny = trans[1] + rot[1][0] * x + rot[1][1] * y + rot[1][2] * z
                nz = trans[2] + rot[2][0] * x + rot[2][1] * y + rot[2][2] * z
                g.write("%s%8.3f%8.3f%8.3f%s" % (line[:30], nx, ny, nz, line[54:]))


def parse_tm_matrix(path):
    rot, trans = [], []
    lines = open(path).read().split("\n")
    for i, line in enumerate(lines):
        if line.strip().startswith("m ") and "t[m]" in line:
            for k in range(3):
                f = lines[i + 1 + k].split()
                trans.append(float(f[1]))
                rot.append([float(f[2]), float(f[3]), float(f[4])])
            break
    return rot, trans


def run_tmalign_like(exe, q, t, qpath, tpath, work):
    mat = os.path.join(work, "%s_%s.mat" % (q, t))
    t0 = time.time()
    subprocess.run([exe, qpath, tpath, "-m", mat], capture_output=True, text=True)
    dt = time.time() - t0
    rot, trans = parse_tm_matrix(mat)
    moved = os.path.join(work, "%s_%s_moved.pdb" % (q, t))
    apply_matrix_pdb(qpath, moved, rot, trans)
    return [(moved, tpath, False, dt)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("method")
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--cols", default="0,1")
    ap.add_argument("--structs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--refs")
    ap.add_argument("--tag", default="")
    ap.add_argument("--v1-results", help="results_levelN.tsv of run_icarus_v1.py")
    ap.add_argument("--exe", help="executable for tmalign/usalign")
    ap.add_argument("--icarus-args", default="", help="extra arguments for icarus pairs")
    ap.add_argument("--threads", type=int, default=0)
    a = ap.parse_args()
    cols = tuple(int(x) for x in a.cols.split(","))
    pairs = read_pairs(a.pairs, cols)
    os.makedirs(a.out, exist_ok=True)
    work = os.path.join(a.out, "work_" + a.method + a.tag)
    os.makedirs(work, exist_ok=True)
    refs = ev.read_reference(a.refs) if a.refs else {}
    spath = lambda s: os.path.join(a.structs, s + ".pdb")

    # method -> {(q, t): (moved_pdb, fixed_pdb, moved_is_second, seconds, extra)}
    results = {}
    if a.method == "icarus2":
        plist = os.path.join(work, "pairs.txt")
        with open(plist, "w") as f:
            for q, t in pairs:
                f.write("%s\t%s\n" % (q, t))
        tsv = os.path.join(work, "icarus2.tsv")
        models = os.path.join(work, "models")
        cmd = [ev.ICARUS_BIN, "pairs", plist, "--dir", a.structs, "--ext", ".pdb", "-o", tsv, "--models", models,
               "-t", str(a.threads)] + a.icarus_args.split()
        subprocess.run(cmd, check=True)
        for line in open(tsv):
            f = line.rstrip("\n").split("\t")
            if f[0] == "query":
                hdr = f
                continue
            rec = dict(zip(hdr, f))
            q, t = rec["query"], rec["target"]
            second = rec["peeled"] == "2"
            moved = os.path.join(models, "%s__%s.pdb" % (q, t))
            fixed = spath(q) if second else spath(t)
            results[(q, t)] = (moved, fixed, second, float(rec["ms"]) / 1000.0, rec["n_bodies"])
    elif a.method == "icarus_v1":
        for line in open(a.v1_results):
            f = line.rstrip("\n").split("\t")
            q, t, score, dt, sol = f[0], f[1], f[2], float(f[3]), f[4]
            if sol == "NA" or not os.path.exists(sol):
                continue
            base = os.path.basename(sol)
            moved_name = base.split("_", 2)[2].split("-level_")[0]
            second = moved_name != q
            mv = os.path.join(work, "%s_%s_moved.pdb" % (q, t))
            fx = os.path.join(work, "%s_%s_fixed.pdb" % (q, t))
            ev.split_models(sol, mv, fx)
            results[(q, t)] = (mv, spath(q) if second else spath(t), second, dt, score)
    elif a.method in ("tmalign", "usalign"):
        for q, t in pairs:
            (mv, fx, second, dt) = run_tmalign_like(a.exe, q, t, spath(q), spath(t), work)[0]
            results[(q, t)] = (mv, fx, second, dt, "")
    else:
        sys.exit("unknown method " + a.method)

    out_path = os.path.join(a.out, a.method + a.tag + ".tsv")
    with open(out_path, "w") as out:
        out.write("query\ttarget\ttm\tref_agreement\tseconds\textra\n")
        for q, t in pairs:
            if (q, t) not in results:
                out.write("%s\t%s\tNA\tNA\tNA\t\n" % (q, t))
                continue
            moved, fixed, second, dt, extra = results[(q, t)]
            lnorm = min(ev.n_residues(spath(q)), ev.n_residues(spath(t)))
            tm, aln = ev.gdt(moved, fixed, lnorm)
            agr = "NA"
            if (q, t) in refs:
                names_q, names_t = ev.residue_names(spath(q)), ev.residue_names(spath(t))
                off = ev.reference_offset(refs[(q, t)], names_q, names_t)
                agr = "%.3f" % ev.agreement(refs[(q, t)], aln, second, off)
            out.write("%s\t%s\t%.4f\t%s\t%.3f\t%s\n" % (q, t, tm, agr, dt, extra))
    print("wrote", out_path)


if __name__ == "__main__":
    main()
