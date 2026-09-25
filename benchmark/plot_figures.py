#!/usr/bin/env python3
"""
Figures for the ICARUS 2 benchmarks (requires matplotlib).

Usage: plot_figures.py --ripc RIPC_OUT --sisy SISY_OUT --scop-icarus icarus2.tsv
                       --scop-foldseek foldseek.m8 --ecoli flexible_alignments.tsv
                       --out figures/
RIPC_OUT / SISY_OUT are run_benchmark.py output directories; the file names
used for each method are listed in METHODS below. Missing inputs skip a figure.
"""

import argparse
import csv
import os
import statistics as st
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

import eval_scop  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# (label, run_benchmark.py output file, colour, marker)
METHODS = [
    ("ICARUS 2", "icarus2_final_mb6.tsv", "#c0392b", "*"),
    ("ICARUS 2 --fast", "icarus2_final_fast.tsv", "#e67e22", "*"),
    ("ICARUS 2 database mode", "icarus2_final_fast_hp02.tsv", "#f1c40f", "*"),
    ("ICARUS v1 (level 2)", "icarus_v1_L2.tsv", "#8e44ad", "D"),
    ("KPAX flexible", "kpax_flex.tsv", "#2980b9", "s"),
    ("FATCAT flexible", "fatcat_flex.tsv", "#16a085", "s"),
    ("KPAX", "kpax.tsv", "#7fb3d5", "o"),
    ("FATCAT rigid", "fatcat_rigid.tsv", "#76d7c4", "o"),
    ("TM-align", "tmalign.tsv", "#7f8c8d", "o"),
    ("Foldseek TM-align", "foldseek_tm.tsv", "#34495e", "^"),
    ("Foldseek LoL-align", "foldseek_lol.tsv", "#5d6d7e", "v"),
    ("Foldseek 3Di+AA", "foldseek_3di.tsv", "#aab7b8", "<"),
]


def read_bench(path):
    rows = {}
    if not os.path.exists(path):
        return rows
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["tm"] != "NA":
            rows[(r["query"], r["target"])] = (float(r["tm"]), float(r["seconds"]) if r["seconds"] else float("nan"))
    return rows


def fig_overview(out):
    fig, ax = plt.subplots(figsize=(14, 3.2))
    ax.axis("off")
    steps = [
        ("Preprocessing", "once per structure:\nDSSP + Protein Peeling\n→ tree of Protein Units"),
        ("Seeds", "8-residue fragment\npairs with similar\ndistance matrices"),
        ("Scoring", "distance grid +\nprefix sums: every\nseed on every PU"),
        ("Placement", "banded DP\n⇄ IRLS superposition\n(monotone TM ascent)"),
        ("Assembly", "exact DP over\nPU-tree cuts with\ndisjoint target spans"),
        ("Finalisation", "chimera + gdt2 DP,\nrefit on free target,\nhinge refinement"),
    ]
    w, gap = 0.145, 0.024
    for i, (title, body) in enumerate(steps):
        x = 0.005 + i * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 0.2), w, 0.66, boxstyle="round,pad=0.005,rounding_size=0.015",
                                    fc="#fdebd0" if i == 0 else "#eaf2f8", ec="#34495e", lw=1.2,
                                    transform=ax.transAxes))
        ax.text(x + w / 2, 0.75, title, ha="center", va="center", fontsize=11, fontweight="bold",
                transform=ax.transAxes)
        ax.text(x + w / 2, 0.45, body, ha="center", va="center", fontsize=9, transform=ax.transAxes)
        if i < len(steps) - 1:
            ax.annotate("", xy=(x + w + gap - 0.002, 0.53), xytext=(x + w + 0.002, 0.53),
                        xycoords="axes fraction", arrowprops=dict(arrowstyle="->", lw=1.4, color="#34495e"))
    ax.text(0.5, 0.06, "each protein is peeled in turn (one direction with --fast); "
            "the best solution of either direction is reported",
            ha="center", fontsize=9.5, style="italic", transform=ax.transAxes)
    fig.savefig(os.path.join(out, "fig1_overview.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_quality_speed(ripc, sisy, out):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, (name, d) in zip(axes, (("RIPC (40 domain pairs)", ripc), ("SISY (68 chain pairs)", sisy))):
        for label, fname, col, mk in METHODS:
            rows = read_bench(os.path.join(d, fname))
            if not rows:
                continue
            tm = st.mean(v[0] for v in rows.values())
            sec = st.mean(v[1] for v in rows.values())
            ax.scatter(sec, tm, s=160 if mk == "*" else 60, c=col, marker=mk, edgecolors="k", linewidths=0.5,
                       zorder=3, label=label)
        ax.set_xscale("log")
        ax.set_xlabel("mean time per pair (s, one core, log scale)")
        ax.set_ylabel("mean TM-score (gdt2, shorter chain)")
        ax.set_title(name)
        ax.grid(alpha=0.3, which="both")
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=4, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, -0.12))
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig2_quality_speed.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_v1(ripc, out):
    a = read_bench(os.path.join(ripc, "icarus2_final_mb5.tsv"))
    b = read_bench(os.path.join(ripc, "icarus_v1_L2.tsv"))
    keys = sorted(set(a) & set(b))
    if not keys:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
    ax = axes[0]
    ax.scatter([b[k][0] for k in keys], [a[k][0] for k in keys], c="#c0392b", edgecolors="k", linewidths=0.5)
    ax.plot([0.4, 1], [0.4, 1], "k--", lw=0.8)
    ax.set_xlim(0.4, 1)
    ax.set_ylim(0.4, 1)
    ax.set_xlabel("ICARUS v1, level 2 (TM-score)")
    ax.set_ylabel("ICARUS 2, ≤5 bodies (TM-score)")
    ax.set_title("A. RIPC, same rigid-body budget")
    ax.grid(alpha=0.3)
    ax = axes[1]
    data, labels = [], []
    for label, fname, _, _ in METHODS:
        rows = read_bench(os.path.join(ripc, fname))
        if rows and label in ("ICARUS 2", "ICARUS 2 --fast", "ICARUS v1 (level 2)", "KPAX flexible",
                              "FATCAT flexible", "TM-align"):
            data.append([v[1] for v in rows.values()])
            labels.append(label.replace(" (level 2)", "").replace("ICARUS 2 ", "ICARUS 2\n"))
    ax.boxplot(data, orientation="horizontal", widths=0.6)
    ax.set_yticks(range(1, len(labels) + 1))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xscale("log")
    ax.set_xlabel("time per pair (s, one core, log scale)")
    ax.set_title("B. RIPC run times")
    ax.grid(alpha=0.3, axis="x", which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig3_v1_and_times.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def scop_sensitivity(labels, targets, scores):
    sens = defaultdict(list)
    for q, ts in targets.items():
        ranked = sorted(ts, key=lambda t: -scores.get((q, t), float("-inf")))
        counts, totals = defaultdict(int), defaultdict(int)
        for t in ts:
            totals[eval_scop.level(labels[q], labels[t])] += 1
        for t in ranked:
            lv = eval_scop.level(labels[q], labels[t])
            if lv == "fp":
                break
            counts[lv] += 1
        for lv in ("family", "superfamily", "fold"):
            if totals[lv]:
                sens[lv].append(counts[lv] / totals[lv])
    return [st.mean(sens[lv]) for lv in ("family", "superfamily", "fold")]


def fig_scop(icarus, foldseek, out):
    labels = {}
    for line in open(os.path.join(HERE, "data", "scop40_labels.tsv")):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sid":
            labels[f[0]] = tuple(f[1:5])
    targets = defaultdict(list)
    for line in open(os.path.join(HERE, "data", "scop40_bench_pairs.tsv")):
        q, t = line.split()
        targets[q].append(t)
    specs = [("ICARUS 2 tm_rigid_max", icarus + ":tm_rigid_max", "#c0392b"),
             ("ICARUS 2 tm_conn", icarus + ":tm_conn", "#e67e22"),
             ("ICARUS 2 tm_flex", icarus + ":tm_flex", "#f5b041"),
             ("Foldseek bits", foldseek + ":3", "#34495e"),
             ("Foldseek E-value", None, "#85929e")]
    res = []
    for name, spec, col in specs:
        if spec is None:
            _, sc = eval_scop.load_scores(foldseek + ":2")
            sc = {k: -v for k, v in sc.items()}
        else:
            _, sc = eval_scop.load_scores(spec)
        res.append((name, scop_sensitivity(labels, targets, sc), col))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    ax = axes[0]
    w = 0.16
    for i, (name, vals, col) in enumerate(res):
        ax.bar([x + (i - 2) * w for x in range(3)], vals, w, label=name, color=col, edgecolor="k", lw=0.4)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["family", "superfamily", "fold"])
    ax.set_ylabel("sensitivity up to the 1st false positive")
    ax.set_title("A. SCOP40 homology discrimination")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    ax = axes[1]
    rows = list(csv.DictReader(open(icarus), delimiter="\t"))
    lab = {eval_scop.strip_ext(k): v for k, v in labels.items()}
    dec, fam = defaultdict(list), defaultdict(list)
    for r in rows:
        q, t = eval_scop.strip_ext(r["query"]), eval_scop.strip_ext(r["target"])
        if q not in lab or t not in lab or q == t:
            continue
        g = dec if lab[q][1] != lab[t][1] else fam if lab[q][3] == lab[t][3] else None
        if g is not None:
            for c in ("tm_flex", "tm_rigid_max"):
                g[c].append(float(r[c]))
    bins = [i / 50 for i in range(51)]
    ax.hist(dec["tm_flex"], bins, density=True, alpha=0.45, color="#f5b041", label="tm_flex, other-fold decoys")
    ax.hist(fam["tm_flex"], bins, density=True, histtype="step", lw=1.6, color="#b9770e",
            label="tm_flex, same family")
    ax.hist(dec["tm_rigid_max"], bins, density=True, alpha=0.45, color="#c0392b",
            label="tm_rigid_max, other-fold decoys")
    ax.hist(fam["tm_rigid_max"], bins, density=True, histtype="step", lw=1.6, color="#7b241c",
            label="tm_rigid_max, same family")
    ax.set_xlabel("score")
    ax.set_ylabel("density")
    ax.set_title("B. Flexible scores inflate on unrelated pairs")
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig4_scop40.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    return res


def fig_ecoli(path, out):
    rows = [r for r in csv.DictReader(open(path), delimiter="\t") if r["query"] != r["target"]]
    rig = [float(r["tm_rigid_max"]) for r in rows]
    con = [float(r["tm_conn"]) for r in rows]
    nb = [int(r["n_bodies"]) for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
    ax = axes[0]
    h = ax.hist2d(rig, con, bins=50, range=[[0, 1], [0, 1]], cmap="magma_r", cmin=1,
                  norm=matplotlib.colors.LogNorm())
    fig.colorbar(h[3], ax=ax, label="pairs")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.plot([0, 0.9], [0.1, 1], ":", c="#2980b9", lw=1)
    ax.set_xlabel("rigid TM-score (tm_rigid_max)")
    ax.set_ylabel("connected flexible TM-score (tm_conn)")
    ax.set_title("A. E. coli paralogue pairs (Foldseek E ≤ 1e-3)")
    ax = axes[1]
    ks = sorted(set(nb))
    ax.bar(ks, [nb.count(k) for k in ks], color="#c0392b", edgecolor="k", lw=0.4)
    ax.set_xlabel("rigid bodies in the flexible solution")
    ax.set_ylabel("pairs")
    ax.set_title("B. Body count (hinge penalty 0.02)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig5_ecoli.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ripc")
    ap.add_argument("--sisy")
    ap.add_argument("--scop-icarus")
    ap.add_argument("--scop-foldseek")
    ap.add_argument("--ecoli")
    ap.add_argument("--out", default=os.path.join(HERE, "figures"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fig_overview(a.out)
    if a.ripc and a.sisy:
        fig_quality_speed(a.ripc, a.sisy, a.out)
    if a.ripc:
        fig_v1(a.ripc, a.out)
    if a.scop_icarus and a.scop_foldseek:
        for name, vals, _ in fig_scop(a.scop_icarus, a.scop_foldseek, a.out):
            print("%-24s family %.3f superfamily %.3f fold %.3f" % (name, *vals))
    if a.ecoli:
        fig_ecoli(a.ecoli, a.out)


if __name__ == "__main__":
    main()
