# ICARUS 2 — fast flexible structural alignment based on Protein Units

ICARUS 2 is a Rust re-implementation and algorithmic redesign of
[ICARUS](https://doi.org/10.1093/bioinformatics/btad459). It keeps the ICARUS
principle — the flexible protein is cut into compact **Protein Units** (PUs,
Protein Peeling) that are superposed independently onto the rigid partner, in
any order, so that hinge motions, circular permutations and domain swaps are
handled — but replaces the exhaustive KPAX/subprocess exploration with an
in-process search that runs in milliseconds, so that it can be used at
proteome scale.

* single static binary, no Perl/Python/KPAX/DSSP dependencies
* PDB and mmCIF input, optionally gzip-compressed (AlphaFold DB files work as is)
* `align` (one pair, report + superposed model), `pairs` (list of pairs),
  `createdb` + `search` (preprocessed database, all-vs-all or prefiltered
  candidate pairs, e.g. from Foldseek), `gdt` (gdt2.pl-compatible scoring)

## Build

```bash
cd icarus-rs
cargo build --release          # binary: target/release/icarus
cargo test --release
```

## Usage

```bash
# one pair: text report, flexibly superposed model in target order
icarus align query.pdb target.cif.gz --out-pdb moved.pdb

# a list of pairs (tab-separated ids), structures read from a directory
icarus pairs pairs.tsv --dir structures/ --ext .pdb -o results.tsv -t 16

# proteome scale: preprocess once, then align candidate pairs
icarus createdb AFDB_proteome_dir/ proteome.icdb --min-plddt 70
foldseek easy-search AFDB_proteome_dir/ AFDB_proteome_dir/ hits.m8 tmp --exhaustive-search 0
icarus search proteome.icdb proteome.icdb --pairs hits.m8 -o flexible.tsv --min-rigid 0.3 \
    --fast --hinge-penalty 0.02
```

`pipeline/proteome_flexdb.sh STRUCT_DIR OUT_DIR [threads] [min_plddt] [min_rigid_tm]`
runs the whole flow (preprocessing, Foldseek prefilter, flexible alignment).

Main options: `--max-bodies` (maximum number of rigid bodies, default 6),
`--hinge-penalty` (TM-score cost per extra body), `--min-pu-size` (default 15),
`--one-direction` (peel only the first structure; ~2× faster), `--fast`
(one direction and a smaller seed/candidate budget; ~2.5× faster),
`--min-plddt` (mask low-confidence residues of predicted models).

For a database, use `--hinge-penalty 0.02`: each extra rigid body must then
improve the TM-score by 0.02. The mean TM-score on RIPC drops only from 0.755
to 0.738, but solutions are simpler (4.4 instead of 5.6 bodies on average) and
easier to interpret. For example, NarL vs RcsB from *E. coli* comes out as two
bodies (receiver domain + HTH domain, TM 0.87 vs 0.24 rigid) instead of six
shuffled β-α units.

### Output columns (`pairs`, `search`)

| column | meaning |
|---|---|
| `tm_flex` | flexible TM-score normalised by the shorter chain (the ICARUS alignment score) |
| `tm_flex_q`, `tm_flex_t` | same alignment normalised by query / target length |
| `tm_rigid` | best rigid-body TM-score found (normalised by the shorter chain) |
| `tm_rigid_max` | rigid TM-score normalised by the longer chain |
| `tm_conn` | connectivity-aware flexible TM-score, normalised by the longer chain (see below) |
| `n_conn` | rigid bodies in the connected run scored by `tm_conn` |
| `n_bodies` | number of rigid bodies (PUs) in the flexible solution |
| `n_aligned`, `n_core`, `rmsd_core` | aligned pairs, pairs within 5 Å, their RMSD |
| `peeled` | which structure was cut into PUs (1 or 2) |
| `bodies` | `qstart-qend:tstart-tend` per body, in target order (author numbering) |
| `transforms` | (`search --transforms`) per-body rotation + translation |

**Which score for what.** `tm_flex` measures how well two structures
superpose when their PUs move independently: it is the score to compare
alignments of related proteins (the ICARUS benchmarks). It is *not* a homology
statistic: with up to six free bodies, small PUs of unrelated proteins always
find somewhere to fit (on SCOP40 decoy pairs from different folds the median
`tm_flex` is 0.67). For database searches and homology decisions use
`tm_rigid_max` (and the prefilter E-value): on the SCOP40 benchmark below it
ranks homologues as well as Foldseek. `tm_conn` — the flexible score restricted
to runs of consecutive bodies whose junctions remain chain-connected — and the
gap `tm_conn - tm_rigid_max` flag homologues related by hinge motions or
rearrangements; they do not improve homology detection.

## Algorithm

For each direction (each protein is peeled in turn, as in ICARUS):

1. **Preprocessing (once per structure).** DSSP secondary structure and
   Protein Peeling (both vendored from SWORD3's Rust port) give the PU
   hierarchy: every PU of every Peeling iteration is a node of a tree.
2. **Seeds.** Aligned fragment pairs are found by comparing intra-fragment
   C-alpha distance matrices (8-residue fragments; rigid-motion invariant, no
   superposition needed), merged along diagonals into gapless blocks; each
   block gives a candidate superposition.
3. **Scoring every seed on every PU at once.** A distance-transform grid of the
   rigid protein gives the distance to the nearest C-alpha with one memory
   access; prefix sums over the query turn the TM kernel of a superposition
   into an O(1) score for every PU.
4. **Placement refinement.** The best diverse candidates of each PU (more for
   small PUs, plus the placements of its parent) are refined by alternating a
   banded dynamic programme — the band is centred on a longest increasing
   subsequence of the nearest-residue map, then on the previous path — and an
   iteratively reweighted superposition that increases the TM objective
   monotonically (minorise-maximise).
5. **Assembly.** An exact dynamic programme selects disjoint PUs from any level
   of the hierarchy (bit masks over the finest PUs), with non-overlapping
   target regions, maximising the summed score minus an optional hinge
   penalty. This generalises the per-level permutation search of ICARUS.
6. **Finalisation.** The rigid bodies are ordered along the target into a
   chimera aligned with the gdt2.pl dynamic programme (the ICARUS scoring);
   each body is re-superposed on its pairs until convergence, then re-placed
   against the target residues left free by the others (the ICARUS
   "updated target" idea), and hinge positions are shifted by up to 12
   residues when that improves the score.

The TM-score reported is exactly the gdt2.pl score used by ICARUS v1
(`icarus gdt` reproduces gdt2.pl), normalised by the shorter chain.

## Provenance

`src/peeling/` and `src/dssp/` are vendored from SWORD3
(DSIMB/sword3 @ e47f837, CeCILL-2.1): Protein Peeling (Gelly et al. 2006) and
DSSP (Kabsch & Sander 1983) ports, trimmed to what ICARUS needs.
