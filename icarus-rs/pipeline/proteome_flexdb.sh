#!/usr/bin/env bash
# Build a precomputed flexible-alignment database for a set of structures
# (e.g. an AlphaFold DB proteome): Foldseek prefilter, then ICARUS 2 flexible
# alignment of every candidate pair.
#
# usage: proteome_flexdb.sh STRUCT_DIR OUT_DIR [threads] [min_plddt] [min_rigid_tm]
# env:   ICARUS (default: icarus), FOLDSEEK (default: foldseek),
#        FOLDSEEK_ARGS (extra search arguments, default: -e 10 --max-seqs 1000),
#        ICARUS_ARGS (extra alignment arguments, default: --fast --hinge-penalty 0.02:
#        one peeling direction, and each extra rigid body must gain 0.02 TM-score,
#        which keeps solutions interpretable, e.g. two domains around one hinge)
set -euo pipefail
IN=$1
OUT=$2
THREADS=${3:-0}
PLDDT=${4:-70}
MINTM=${5:-0.0}
ICARUS=${ICARUS:-icarus}
FOLDSEEK=${FOLDSEEK:-foldseek}
FOLDSEEK_ARGS=${FOLDSEEK_ARGS:--e 10 --max-seqs 1000}
ICARUS_ARGS=${ICARUS_ARGS:---fast --hinge-penalty 0.02}
FS_THREADS=$THREADS
[ "$FS_THREADS" = "0" ] && FS_THREADS=$(nproc)
mkdir -p "$OUT"
log() { echo "[$(date +%H:%M:%S)] (+${SECONDS}s) $*" >&2; }

log "1/3 preprocessing structures (DSSP, Protein Peeling; pLDDT >= $PLDDT)"
"$ICARUS" createdb "$IN" "$OUT/structures.icdb" --min-plddt "$PLDDT" -t "$THREADS" \
    2> >(tee "$OUT/createdb.log" >&2)

log "2/3 candidate pairs with the Foldseek prefilter"
# AlphaFold DB proteomes ship each model as .cif.gz and .pdb.gz: index one copy
INCLUDE='.*'
if ls "$IN" | grep -q 'cif.gz$' && ls "$IN" | grep -q 'pdb.gz$'; then INCLUDE='cif.gz$'; fi
"$FOLDSEEK" createdb "$IN" "$OUT/fsdb" --threads "$FS_THREADS" --mask-bfactor-threshold "$PLDDT" \
    --file-include "$INCLUDE" -v 1
"$FOLDSEEK" search "$OUT/fsdb" "$OUT/fsdb" "$OUT/fsaln" "$OUT/fstmp" --threads "$FS_THREADS" \
    $FOLDSEEK_ARGS -v 1 2> "$OUT/foldseek_search.log"
"$FOLDSEEK" convertalis "$OUT/fsdb" "$OUT/fsdb" "$OUT/fsaln" "$OUT/candidates.m8" \
    --format-output query,target,evalue,bits --threads "$FS_THREADS" -v 1
log "   $(wc -l < "$OUT/candidates.m8") candidate hits"

log "3/3 flexible alignment of candidate pairs (reporting rigid TM >= $MINTM)"
"$ICARUS" search "$OUT/structures.icdb" "$OUT/structures.icdb" --pairs "$OUT/candidates.m8" \
    --min-rigid "$MINTM" -t "$THREADS" --transforms $ICARUS_ARGS -o "$OUT/flexible_alignments.tsv" 2> >(tee "$OUT/search.log" >&2)
log "done: $OUT/flexible_alignments.tsv"
