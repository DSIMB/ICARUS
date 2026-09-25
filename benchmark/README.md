# ICARUS benchmarks

Scripts to reproduce the comparison of ICARUS 2 (`icarus-rs/`) with ICARUS v1
and other structural aligners.

## Datasets

* **RIPC** (Mayr et al. 2007, BMC Struct Biol 7:50): 40 pairs of SCOP domains
  related by Repetitions, Insertions/deletions, circular Permutations and
  Conformational variability; 23 pairs have reference alignments
  (`data/ripc_pairs.tsv`, `data/ripc_reference_alignments.txt`, from the
  paper's open-access supplementary files).
* **SISY** (same paper): 68 pairs of PDB chains (`data/sisy_pairs.tsv`; the pair
  involving the obsoleted entry 1VLF is excluded).
* **SCOP40 discrimination**: 250 queries from SCOP classes a-d that have a
  same-superfamily/other-family member; targets are their same-fold domains
  (up to 60) plus 250 random other-fold decoys each
  (`data/scop40_bench_pairs.tsv`, labels in `data/scop40_labels.tsv`). Domains
  come from the Foldseek benchmark set `scop40pdb.tar.gz`.

```bash
python3 build_datasets.py STRUCT_DIR              # RIPC + SISY structures from RCSB (SCOPe mapping)
python3 scop40_labels.py SCOP40_PDB_DIR data/scop40_labels.tsv
python3 make_scop_benchmark.py data/scop40_labels.tsv data/scop40_bench_pairs.tsv 250 250 1
```

## Scoring

All methods are scored identically from their superposed coordinates
(`evaluate.py`): the moved structure (for flexible methods, the chimera of
transformed rigid bodies ordered along the fixed structure) and the fixed
structure are aligned with the gdt2.pl dynamic programme (`icarus gdt`, an exact
re-implementation), and the TM-score is normalised by the shorter of the two
original structures, as in the ICARUS paper. Reference agreement is the
fraction of reference residue pairs reproduced by that alignment.

```bash
B="--pairs data/ripc_pairs.tsv --cols 1,2 --structs STRUCT_DIR/ripc --out OUT --refs data/ripc_reference_alignments.txt"
python3 run_benchmark.py icarus2 $B --icarus-args "--max-bodies 5"
python3 run_benchmark.py tmalign $B --exe TMalign
python3 run_benchmark.py kpax_flex $B
python3 run_benchmark.py fatcat_flex $B --exe "BIOJAVA_CLASSPATH"   # fatcat/RunFatcat.java
python3 run_benchmark.py foldseek_tm $B --exe foldseek
python3 run_icarus_v1.py data/ripc_pairs.tsv STRUCT_DIR/ripc V1_OUT --level 2 --cols 1,2
python3 run_benchmark.py icarus_v1 $B --v1-results V1_OUT/results_level2.tsv
```

Discrimination: sensitivity up to the first false positive (Foldseek protocol),
averaged over queries, at family / superfamily / fold level:

```bash
icarus createdb SCOP40_DIR scop40.icdb
icarus search scop40.icdb scop40.icdb --pairs data/scop40_bench_pairs.tsv -o icarus2.tsv
python3 eval_scop.py data/scop40_labels.tsv data/scop40_bench_pairs.tsv icarus2.tsv:tm_rigid_max icarus2.tsv:tm_conn
```

## Figures

`plot_figures.py` (needs matplotlib) draws the figures in `figures/` from the
outputs above: quality versus speed on RIPC and SISY, ICARUS 2 versus v1,
SCOP40 discrimination and score distributions, and a proteome run.

```bash
python3 plot_figures.py --ripc RIPC_OUT --sisy SISY_OUT --scop-icarus icarus2.tsv \
    --scop-foldseek foldseek.m8 --ecoli flexible_alignments.tsv --out figures/
```
