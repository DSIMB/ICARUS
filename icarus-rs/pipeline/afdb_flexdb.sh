#!/usr/bin/env bash
# Flexible-alignment database for AlphaFold DB sets on a large machine:
# Swiss-Prot, the model-organism / global-health proteomes, or any directory of
# structures. Each step is resumable: re-running the same command skips what is
# already done, and alignment runs chunk by chunk.
#
#   afdb_flexdb.sh list                          bulk archives on the AFDB FTP site
#   afdb_flexdb.sh fetch swissprot DATA_DIR       Swiss-Prot models (PDB format, ~27 GB)
#   afdb_flexdb.sh fetch proteome NAME DATA_DIR   one proteome, e.g. UP000005640 (human)
#   afdb_flexdb.sh fetch proteomes DATA_DIR       all proteome archives, into one directory
#   afdb_flexdb.sh run STRUCT_DIR OUT_DIR         preprocess, prefilter, align, merge
#
# Environment (defaults in brackets):
#   THREADS [nproc]        threads for ICARUS and Foldseek
#   PLDDT [70]             residues with lower pLDDT are ignored (0 = keep all)
#   MAX_EVALUE [1e-3]      Foldseek E-value cut-off for candidate pairs
#   MAX_SEQS [1000]        Foldseek hits kept per query
#   CHUNK [2000000]        candidate pairs aligned per `icarus search` call
#   MIN_RIGID [0]          only report pairs with tm_rigid_max >= MIN_RIGID
#   ICARUS_ARGS [--fast --hinge-penalty 0.02]   database mode (see README)
#   FOLDSEEK_ARGS []       extra Foldseek search arguments, e.g. "--gpu 1"
#   ICARUS [icarus], FOLDSEEK [foldseek], AFDB [https://ftp.ebi.ac.uk/pub/databases/alphafold]
#
# Output in OUT_DIR: structures.icdb (preprocessed structures), candidates.m8
# (Foldseek hits), pairs.tsv (unordered candidate pairs), results/ (one table per
# chunk) and flexible_alignments.tsv.gz (all results, one header).
set -euo pipefail

THREADS=${THREADS:-$(nproc)}
PLDDT=${PLDDT:-70}
MAX_EVALUE=${MAX_EVALUE:-1e-3}
MAX_SEQS=${MAX_SEQS:-1000}
CHUNK=${CHUNK:-2000000}
MIN_RIGID=${MIN_RIGID:-0}
ICARUS=${ICARUS:-icarus}
FOLDSEEK=${FOLDSEEK:-foldseek}
ICARUS_ARGS=${ICARUS_ARGS:---fast --hinge-penalty 0.02}
FOLDSEEK_ARGS=${FOLDSEEK_ARGS:-}
AFDB=${AFDB:-https://ftp.ebi.ac.uk/pub/databases/alphafold}

log() { echo "[$(date '+%F %T')] $*" >&2; }
die() { echo "error: $*" >&2; exit 1; }

# archives of the AFDB bulk-download index: "swissprot", "proteome" (model organisms and
# global-health proteomes) or "all"; prints name, structures, size, species
archives() {
    curl -fsSL "$AFDB/download_metadata.json" | python3 -c '
import json, sys
kind = sys.argv[1]
for a in json.load(sys.stdin):
    name = a["archive_name"]
    t = a.get("type", "")
    if kind == "all" or (kind == "swissprot" and t == "swissprot") or (kind == "proteome" and t in ("proteome", "global_health")):
        label = a.get("species") or a.get("label") or ""
        print("%s\t%d\t%.1f GB\t%s" % (name, a["num_predicted_structures"], a["size_bytes"] / 1e9, label))
' "$1"
}

# download one archive (resumable) and unpack its PDB files into DEST
fetch_archive() {
    local name=$1 dest=$2 dl=$3
    mkdir -p "$dest" "$dl"
    if [ -e "$dl/$name.unpacked" ]; then log "$name already unpacked"; return; fi
    log "downloading $name"
    curl -fL -# --retry 10 --retry-delay 30 -C - -o "$dl/$name" "$AFDB/latest/$name"
    log "unpacking $name into $dest"
    # proteome archives hold each model as .cif.gz and .pdb.gz: keep one copy
    case "$name" in
        swissprot_*) tar -xf "$dl/$name" -C "$dest" ;;
        *) tar -xf "$dl/$name" -C "$dest" --wildcards '*.pdb.gz' ;;
    esac
    touch "$dl/$name.unpacked"
    [ -n "${KEEP_TAR:-}" ] || rm -f "$dl/$name"
}

cmd_fetch() {
    local what=${1:-}
    case "$what" in
        swissprot)
            local data=${2:?DATA_DIR}
            local name
            name=$(archives swissprot | awk -F'\t' '$1 ~ /pdb/ && !f { print $1; f = 1 }')
            [ -n "$name" ] || die "no Swiss-Prot PDB archive listed in $AFDB/download_metadata.json"
            fetch_archive "$name" "$data/swissprot" "$data/downloads"
            log "structures in $data/swissprot" ;;
        proteome)
            local id=${2:?NAME (e.g. UP000005640)} data=${3:?DATA_DIR}
            local name
            name=$(archives proteome | awk -F'\t' -v id="$id" 'index($1, id) == 1 && !f { print $1; f = 1 }')
            [ -n "$name" ] || die "no proteome archive matches $id (see: $0 list)"
            fetch_archive "$name" "$data/${name%.tar}" "$data/downloads"
            log "structures in $data/${name%.tar}" ;;
        proteomes)
            local data=${2:?DATA_DIR}
            archives proteome | cut -f1 | while read -r name; do
                fetch_archive "$name" "$data/proteomes" "$data/downloads"
            done
            log "structures in $data/proteomes" ;;
        *) die "fetch swissprot|proteome|proteomes ..." ;;
    esac
}

cmd_run() {
    local in=${1:?STRUCT_DIR} out=${2:?OUT_DIR}
    [ -d "$in" ] || die "$in is not a directory"
    command -v "$ICARUS" > /dev/null || die "icarus not found (set ICARUS=/path/to/icarus)"
    command -v "$FOLDSEEK" > /dev/null || die "foldseek not found (set FOLDSEEK=/path/to/foldseek)"
    mkdir -p "$out/results" "$out/tmp"
    local nfiles
    nfiles=$(find "$in" -type f \( -name '*.pdb*' -o -name '*.cif*' -o -name '*.ent*' \) | wc -l)
    log "$nfiles structure files in $in; $THREADS threads; pLDDT >= $PLDDT; E-value <= $MAX_EVALUE"

    if [ ! -e "$out/.createdb.done" ]; then
        log "1/4 preprocessing (DSSP, Protein Peeling)"
        "$ICARUS" createdb "$in" "$out/structures.icdb" --min-plddt "$PLDDT" -t "$THREADS" \
            2> "$out/createdb.log" || { tail "$out/createdb.log"; die "icarus createdb failed"; }
        tail -1 "$out/createdb.log" >&2
        touch "$out/.createdb.done"
    fi

    if [ ! -e "$out/.foldseek.done" ]; then
        log "2/4 candidate pairs with Foldseek"
        local include='.*'
        if find "$in" -name '*.cif.gz' -print -quit | grep -q . && find "$in" -name '*.pdb.gz' -print -quit | grep -q .; then
            include='pdb.gz$'
        fi
        [ -e "$out/fsdb.dbtype" ] || "$FOLDSEEK" createdb "$in" "$out/fsdb" --threads "$THREADS" \
            --mask-bfactor-threshold "$PLDDT" --file-include "$include" -v 1
        # shellcheck disable=SC2086
        "$FOLDSEEK" search "$out/fsdb" "$out/fsdb" "$out/fsaln" "$out/tmp/foldseek" --threads "$THREADS" \
            -e "$MAX_EVALUE" --max-seqs "$MAX_SEQS" $FOLDSEEK_ARGS -v 1 2> "$out/foldseek.log"
        "$FOLDSEEK" convertalis "$out/fsdb" "$out/fsdb" "$out/fsaln" "$out/candidates.m8" \
            --format-output query,target,evalue,bits --threads "$THREADS" -v 1
        rm -rf "$out/tmp/foldseek"
        log "   $(wc -l < "$out/candidates.m8") Foldseek hits"
        touch "$out/.foldseek.done"
    fi

    if [ ! -e "$out/.pairs.done" ]; then
        log "3/4 unordered candidate pairs, split into chunks of $CHUNK"
        awk -F'\t' '$1 != $2 { if ($1 < $2) print $1 "\t" $2; else print $2 "\t" $1 }' "$out/candidates.m8" \
            | LC_ALL=C sort -u -S 25% --parallel="$THREADS" -T "$out/tmp" > "$out/pairs.tsv"
        rm -rf "$out/chunks"
        mkdir -p "$out/chunks"
        split -l "$CHUNK" -d -a 5 "$out/pairs.tsv" "$out/chunks/pairs."
        touch "$out/.pairs.done"
    fi
    local npairs
    npairs=$(wc -l < "$out/pairs.tsv")
    [ "$npairs" -gt 0 ] || { log "no candidate pairs: nothing to align"; return; }
    log "   $npairs pairs; at ~45 ms per pair per core: ~$(awk -v n="$npairs" -v t="$THREADS" 'BEGIN { printf "%.2f h on %d threads (%.1f core-hours)", n * 0.045 / 3600 / t, t, n * 0.045 / 3600 }')"

    log "4/4 flexible alignment"
    shopt -s nullglob
    local c o
    for c in "$out"/chunks/pairs.*; do
        o="$out/results/$(basename "$c").tsv"
        [ -s "$o" ] && continue
        log "   chunk $(basename "$c") ($(wc -l < "$c") pairs)"
        # shellcheck disable=SC2086
        "$ICARUS" search "$out/structures.icdb" "$out/structures.icdb" --pairs "$c" \
            --min-rigid "$MIN_RIGID" -t "$THREADS" --transforms $ICARUS_ARGS -o "$o.tmp" \
            2>> "$out/search.log" || die "icarus search failed on $c (see $out/search.log)"
        mv "$o.tmp" "$o"
    done

    log "merging results"
    local first
    first=$(ls "$out"/results/*.tsv | head -1)
    {
        head -1 "$first"
        for o in "$out"/results/*.tsv; do tail -n +2 "$o"; done
    } | gzip > "$out/flexible_alignments.tsv.gz.tmp"
    mv "$out/flexible_alignments.tsv.gz.tmp" "$out/flexible_alignments.tsv.gz"
    log "done: $out/flexible_alignments.tsv.gz ($(zcat "$out/flexible_alignments.tsv.gz" | tail -n +2 | wc -l) alignments)"
}

case "${1:-}" in
    list) archives all ;;
    fetch) shift; cmd_fetch "$@" ;;
    run) shift; cmd_run "$@" ;;
    *) sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
