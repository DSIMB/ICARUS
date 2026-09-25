#!/usr/bin/env python3
"""
Build the benchmark structure sets used to evaluate ICARUS.

RIPC (Mayr et al. 2007, BMC Struct Biol 7:50): 40 pairs of SCOP 1.7x domains with
Repetitions, Insertions/deletions, circular Permutations and Conformational
variability; 23 pairs carry reference alignments.
SISY (same paper): 69 pairs of PDB chains derived from the SISYPHUS database.

SCOP 1.7x domain identifiers (e.g. d1hcy_2, d1adl__) are mapped onto SCOPe
identifiers, whose residue ranges are read from the RCSB data API. Domains are cut
from RCSB mmCIF files and written as PDB files that keep the author residue
numbering, so the reference alignments (given in author numbering) stay valid.

Usage: build_datasets.py OUT_DIR
"""

import concurrent.futures as cf
import json
import os
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
GRAPHQL = "https://data.rcsb.org/graphql"
CIF_URL = "https://files.rcsb.org/download/{}.cif"

AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
       "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
# Modified residues commonly found in the benchmark entries, mapped to their parent
MODIFIED = {"MSE": "MET", "SEP": "SER", "TPO": "THR", "PTR": "TYR", "CSO": "CYS",
            "CME": "CYS", "HYP": "PRO", "MLY": "LYS", "KCX": "LYS", "CSD": "CYS",
            "OCS": "CYS", "CAS": "CYS", "LLP": "LYS", "M3L": "LYS", "SAC": "SER"}


def fetch(url, data=None, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"} if data else {})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception as e:  # network hiccups: retry with backoff
            if i == tries - 1:
                raise
            time.sleep(2 ** i)


def scop_features(pdb_id):
    """Return {sid: (asym_id, auth_asym_id, [(beg, end), ...])} for one entry."""
    q = """{ entry(entry_id: "%s") { polymer_entities { polymer_entity_instances {
             rcsb_id rcsb_polymer_entity_instance_container_identifiers { asym_id auth_asym_id }
             rcsb_polymer_instance_feature { type feature_positions { beg_seq_id end_seq_id }
               additional_properties { name values } } } } } }""" % pdb_id.upper()
    d = json.loads(fetch(GRAPHQL, json.dumps({"query": q}).encode()))
    entry = d["data"]["entry"]
    out, chains = {}, {}
    if entry is None:
        return out, chains
    for ent in entry["polymer_entities"] or []:
        for inst in ent["polymer_entity_instances"] or []:
            ids = inst["rcsb_polymer_entity_instance_container_identifiers"]
            chains[ids["auth_asym_id"]] = ids["asym_id"]
            for f in inst["rcsb_polymer_instance_feature"] or []:
                if f["type"] != "SCOP":
                    continue
                sid = None
                for p in f["additional_properties"] or []:
                    if p["name"] == "SCOP_DOMAIN_ID":
                        sid = p["values"][0]
                if sid is None:
                    continue
                rng = [(fp["beg_seq_id"], fp["end_seq_id"]) for fp in f["feature_positions"]]
                prev = out.get(sid, (ids["asym_id"], ids["auth_asym_id"], []))
                out[sid] = (ids["asym_id"], ids["auth_asym_id"], prev[2] + rng)
    return out, chains


def parse_cif_atoms(text):
    """Minimal mmCIF atom_site parser (first model only)."""
    lines = text.splitlines()
    cols, rows, i = [], [], 0
    while i < len(lines):
        if lines[i].startswith("_atom_site."):
            while i < len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].split(".", 1)[1].strip())
                i += 1
            while i < len(lines) and not lines[i].startswith("#") and not lines[i].startswith("loop_"):
                if lines[i].strip():
                    rows.append(lines[i].split())
                i += 1
            break
        i += 1
    idx = {c: k for k, c in enumerate(cols)}
    first_model = None
    atoms = []
    for r in rows:
        if len(r) != len(cols):
            continue
        model = r[idx["pdbx_PDB_model_num"]]
        if first_model is None:
            first_model = model
        if model != first_model:
            break
        atoms.append({c: r[k] for c, k in idx.items()})
    return atoms


def write_domain(atoms, asym_id, ranges, path):
    """Write the atoms of one domain (label_seq_id within ranges) as a PDB file."""
    seen_alt = {}
    n = 0
    out = []
    for a in atoms:
        if a["label_asym_id"] != asym_id:
            continue
        try:
            seq = int(a["label_seq_id"])
        except ValueError:
            continue
        if ranges and not any(b <= seq <= e for b, e in ranges):
            continue
        resn = a["label_comp_id"]
        if resn not in AA3:
            if resn not in MODIFIED:
                continue
            resn = MODIFIED[resn]
        alt = a.get("label_alt_id", ".")
        key = (seq, a["label_atom_id"])
        if alt not in (".", "?"):
            if key in seen_alt:
                continue
            seen_alt[key] = alt
        name = a["label_atom_id"].strip('"')
        if a.get("type_symbol") == "SE" and name == "SE":
            name = "SD"
        el = a.get("type_symbol", name[0])
        if el == "SE":
            el = "S"
        ins = a.get("pdbx_PDB_ins_code", "?")
        ins = " " if ins in ("?", ".") else ins
        chain = a["auth_asym_id"][0]
        n += 1
        pname = name if len(name) == 4 else " " + name
        out.append("ATOM  %5d %-4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" % (
            n % 100000, pname, resn, chain, int(a["auth_seq_id"]), ins,
            float(a["Cartn_x"]), float(a["Cartn_y"]), float(a["Cartn_z"]),
            float(a.get("occupancy", 1.0)), float(a.get("B_iso_or_equiv", 0.0)), el))
    out.append("TER\nEND\n")
    with open(path, "w") as f:
        f.writelines(out)
    return n


def scope_candidates(old_id):
    """SCOP 1.7x -> SCOPe identifier candidates (chainless '_' became chain 'a')."""
    pdb, ch, dom = old_id[1:5], old_id[5], old_id[6]
    chains = [ch] if ch != "_" else ["a", "_"]
    return pdb, ch, dom, ["d%s%s%s" % (pdb, c, dom) for c in chains]


# Entries whose SCOP 1.7x identifier has no direct SCOPe/RCSB counterpart
OVERRIDES = {
    "d1jwyb_": ("1jwy", ["d1jwya4"]),   # dynamin GTPase domain, remediated from chain B into chain A
    "d1k87a2": ("4o8a", ["d4o8aa2"]),   # 1K87 obsoleted, replaced by 4O8A (same PutA construct)
}
# Domains defined by an insertion-code block (saposin-like insert of phytepsin: 1S..104S)
INS_CODE_FILTER = {"d1qdma1": "S"}


def subtract_nested(sid, feats):
    """RCSB reports some SCOPe domains as the whole chain range; remove the other
    domains of the same chain that are nested inside it."""
    asym, auth, ranges = feats[sid]
    keep = set()
    for b, e in ranges:
        keep.update(range(b, e + 1))
    for other, (a2, _, r2) in feats.items():
        if other == sid or a2 != asym:
            continue
        inner = set()
        for b, e in r2:
            inner.update(range(b, e + 1))
        if inner < keep:
            keep -= inner
    pos = sorted(keep)
    out, start = [], pos[0]
    for x, y in zip(pos, pos[1:] + [None]):
        if y != x + 1:
            out.append((start, x))
            start = y
    return out


def build_domain(old_id, out_dir):
    path = os.path.join(out_dir, old_id + ".pdb")
    if os.path.exists(path) and os.path.getsize(path) > 100:
        return old_id, "cached"
    pdb, ch, dom, cands = scope_candidates(old_id)
    if old_id in OVERRIDES:
        pdb, cands = OVERRIDES[old_id]
    feats, chains = scop_features(pdb)
    if not chains:
        return old_id, "ENTRY NOT FOUND (obsolete?)"
    hit = next((c for c in cands if c in feats), None)
    if hit:
        asym = feats[hit][0]
        ranges = subtract_nested(hit, feats)
    elif dom == "_":
        # whole chain domain: first chain (chainless entry) or the named chain
        auth = ch.upper() if ch != "_" else sorted(chains)[0]
        if auth not in chains:
            return old_id, "chain %s not found" % auth
        asym, ranges = chains[auth], []
    else:
        return old_id, "domain not found; available: %s" % sorted(feats)
    atoms = parse_cif_atoms(fetch(CIF_URL.format(pdb.upper())).decode())
    if old_id in INS_CODE_FILTER:
        code = INS_CODE_FILTER[old_id]
        atoms = [a for a in atoms if a.get("pdbx_PDB_ins_code") == code]
    n = write_domain(atoms, asym, ranges, path)
    return old_id, "ok (%s, %d atoms)" % (hit or "whole chain", n)


def build_chain(chain_id, out_dir):
    """SISY entries: 5-character ids (pdb + auth chain, '_' = first chain)."""
    path = os.path.join(out_dir, chain_id + ".pdb")
    if os.path.exists(path) and os.path.getsize(path) > 100:
        return chain_id, "cached"
    pdb, ch = chain_id[:4], chain_id[4]
    _, chains = scop_features(pdb)
    if not chains:
        return chain_id, "ENTRY NOT FOUND (obsolete?)"
    auth = ch if ch != "_" else sorted(chains)[0]
    if auth not in chains:
        return chain_id, "chain %s not found in %s" % (auth, sorted(chains))
    atoms = parse_cif_atoms(fetch(CIF_URL.format(pdb.upper())).decode())
    n = write_domain(atoms, chains[auth], [], path)
    return chain_id, "ok (%d atoms)" % n


def read_pairs(path, c1, c2):
    pairs = []
    with open(path) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(c1, c2):
                pairs.append((p[c1], p[c2]))
    return pairs


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "structures")
    ripc_dir, sisy_dir = os.path.join(out, "ripc"), os.path.join(out, "sisy")
    os.makedirs(ripc_dir, exist_ok=True)
    os.makedirs(sisy_dir, exist_ok=True)
    ripc = read_pairs(os.path.join(HERE, "data", "ripc_pairs.tsv"), 1, 2)
    sisy = read_pairs(os.path.join(HERE, "data", "sisy_pairs.tsv"), 0, 1)
    doms = sorted({d for p in ripc for d in p})
    chs = sorted({c for p in sisy for c in p})
    with cf.ThreadPoolExecutor(8) as ex:
        for d, msg in ex.map(lambda d: build_domain(d, ripc_dir), doms):
            print("RIPC", d, msg)
        for c, msg in ex.map(lambda c: build_chain(c, sisy_dir), chs):
            print("SISY", c, msg)


if __name__ == "__main__":
    main()
