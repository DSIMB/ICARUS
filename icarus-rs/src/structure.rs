//! Protein structure input: PDB and mmCIF (optionally gzip-compressed).
//!
//! One polypeptide chain is extracted per structure (first model only). Each
//! residue with a C-alpha atom becomes one position; backbone atoms are kept
//! for DSSP and, on request, all atoms are kept to write superposed models.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::dssp::Backbone;
use crate::geom::V3;

/// Author residue identifier (number + insertion code).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ResId {
    pub num: i32,
    pub icode: u8,
}

impl std::fmt::Display for ResId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.icode == b' ' {
            write!(f, "{}", self.num)
        } else {
            write!(f, "{}{}", self.num, self.icode as char)
        }
    }
}

/// One atom kept for model output.
#[derive(Debug, Clone)]
pub struct Atom {
    pub res: u32,
    pub name: [u8; 4],
    pub element: [u8; 2],
    pub xyz: V3,
    pub bfactor: f32,
}

#[derive(Debug, Clone, Default)]
pub struct Structure {
    pub name: String,
    pub chain: String,
    /// One-letter residue codes.
    pub seq: Vec<u8>,
    /// Three-letter residue names.
    pub resn: Vec<[u8; 3]>,
    pub ca: Vec<V3>,
    pub resid: Vec<ResId>,
    /// C-alpha B-factor (pLDDT for predicted models).
    pub bfac: Vec<f32>,
    pub backbone: Vec<Backbone>,
    pub atoms: Vec<Atom>,
}

impl Structure {
    pub fn len(&self) -> usize {
        self.ca.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ca.is_empty()
    }

    /// Keep only residues whose C-alpha B-factor is >= `min_b` (e.g. pLDDT
    /// masking of low-confidence regions of AlphaFold models).
    pub fn mask_low_confidence(&mut self, min_b: f32) {
        let keep: Vec<bool> = self.bfac.iter().map(|&b| b >= min_b).collect();
        if keep.iter().all(|&k| k) {
            return;
        }
        let mut new_index = vec![u32::MAX; keep.len()];
        let mut k = 0u32;
        for (i, &kp) in keep.iter().enumerate() {
            if kp {
                new_index[i] = k;
                k += 1;
            }
        }
        fn filt<T: Clone>(v: &[T], keep: &[bool]) -> Vec<T> {
            v.iter().zip(keep).filter(|(_, &k)| k).map(|(x, _)| x.clone()).collect()
        }
        self.seq = filt(&self.seq, &keep);
        self.resn = filt(&self.resn, &keep);
        self.ca = filt(&self.ca, &keep);
        self.resid = filt(&self.resid, &keep);
        self.bfac = filt(&self.bfac, &keep);
        self.backbone = filt(&self.backbone, &keep);
        self.atoms.retain(|a| keep[a.res as usize]);
        for a in self.atoms.iter_mut() {
            a.res = new_index[a.res as usize];
        }
    }
}

pub fn three_to_one(r: &[u8]) -> Option<u8> {
    Some(match r {
        b"ALA" => b'A',
        b"ARG" => b'R',
        b"ASN" => b'N',
        b"ASP" => b'D',
        b"CYS" => b'C',
        b"GLN" => b'Q',
        b"GLU" => b'E',
        b"GLY" => b'G',
        b"HIS" => b'H',
        b"ILE" => b'I',
        b"LEU" => b'L',
        b"LYS" => b'K',
        b"MET" => b'M',
        b"PHE" => b'F',
        b"PRO" => b'P',
        b"SER" => b'S',
        b"THR" => b'T',
        b"TRP" => b'W',
        b"TYR" => b'Y',
        b"VAL" => b'V',
        // common modified residues, mapped to their parent amino acid
        b"MSE" => b'M',
        b"SEP" => b'S',
        b"TPO" => b'T',
        b"PTR" => b'Y',
        b"CSO" | b"CME" | b"CSD" | b"OCS" | b"CAS" | b"CSS" | b"CYX" => b'C',
        b"HYP" => b'P',
        b"MLY" | b"KCX" | b"LLP" | b"M3L" => b'K',
        b"SAC" => b'S',
        b"HID" | b"HIE" | b"HIP" | b"HSD" | b"HSE" => b'H',
        b"SEC" => b'U',
        b"PYL" => b'O',
        b"UNK" => b'X',
        _ => return None,
    })
}

struct RawAtom {
    chain: String,
    resnum: i32,
    icode: u8,
    resn: [u8; 3],
    name: [u8; 4],
    element: [u8; 2],
    alt: u8,
    xyz: V3,
    bfactor: f32,
    hetatm: bool,
}

fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let is_gz = path.extension().map(|e| e == "gz").unwrap_or(false);
    Ok(if is_gz {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    })
}

fn pad3(s: &[u8]) -> [u8; 3] {
    let mut o = [b' '; 3];
    for (k, &c) in s.iter().take(3).enumerate() {
        o[k] = c;
    }
    o
}

fn pad4_atom(s: &[u8]) -> [u8; 4] {
    let mut o = [b' '; 4];
    for (k, &c) in s.iter().take(4).enumerate() {
        o[k] = c;
    }
    o
}

fn parse_pdb_atoms(reader: Box<dyn BufRead>) -> Result<Vec<RawAtom>> {
    let mut atoms = Vec::new();
    let mut seen_model = false;
    for line in reader.split(b'\n') {
        let line = line?;
        if line.starts_with(b"MODEL") {
            if seen_model {
                break;
            }
            seen_model = true;
            continue;
        }
        if line.starts_with(b"ENDMDL") {
            break;
        }
        let hetatm = line.starts_with(b"HETATM");
        if !(line.starts_with(b"ATOM  ") || hetatm) || line.len() < 54 {
            continue;
        }
        let field = |a: usize, b: usize| -> &[u8] { &line[a.min(line.len())..b.min(line.len())] };
        let s = |a: usize, b: usize| std::str::from_utf8(field(a, b)).unwrap_or("").trim().to_string();
        let resnum: i32 = match s(22, 26).parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let x: f64 = s(30, 38).parse().unwrap_or(f64::NAN);
        let y: f64 = s(38, 46).parse().unwrap_or(f64::NAN);
        let z: f64 = s(46, 54).parse().unwrap_or(f64::NAN);
        if !(x.is_finite() && y.is_finite() && z.is_finite()) {
            continue;
        }
        let name_field = field(12, 16);
        let el = s(76, 78);
        let mut element = [b' '; 2];
        if !el.is_empty() {
            let eb = el.as_bytes();
            if eb.len() == 1 {
                element[1] = eb[0];
            } else {
                element = [eb[0], eb[1]];
            }
        }
        atoms.push(RawAtom {
            chain: s(21, 22),
            resnum,
            icode: *field(26, 27).first().unwrap_or(&b' '),
            resn: pad3(s(17, 20).as_bytes()),
            name: pad4_atom(name_field),
            element,
            alt: *field(16, 17).first().unwrap_or(&b' '),
            xyz: [x, y, z],
            bfactor: s(60, 66).parse().unwrap_or(0.0),
            hetatm,
        });
    }
    Ok(atoms)
}

/// Split one mmCIF data line into tokens, honouring single/double quotes.
fn cif_tokens(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let b = line.as_bytes();
    let mut i = 0;
    while i < b.len() {
        while i < b.len() && b[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= b.len() {
            break;
        }
        if b[i] == b'\'' || b[i] == b'"' {
            let q = b[i];
            let start = i + 1;
            let mut j = start;
            while j < b.len() && !(b[j] == q && (j + 1 == b.len() || b[j + 1].is_ascii_whitespace())) {
                j += 1;
            }
            out.push(line[start..j.min(b.len())].to_string());
            i = j + 1;
        } else {
            let start = i;
            while i < b.len() && !b[i].is_ascii_whitespace() {
                i += 1;
            }
            out.push(line[start..i].to_string());
        }
    }
    out
}

fn parse_cif_atoms(reader: Box<dyn BufRead>) -> Result<Vec<RawAtom>> {
    let mut atoms = Vec::new();
    let mut cols: Vec<String> = Vec::new();
    let mut in_header = false;
    let mut in_loop = false;
    let mut first_model: Option<String> = None;
    let mut pending: Vec<String> = Vec::new();
    let idx = |cols: &Vec<String>, name: &str| cols.iter().position(|c| c == name);
    let mut map: Vec<Option<usize>> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with("_atom_site.") {
            if !in_header {
                cols.clear();
            }
            in_header = true;
            cols.push(line["_atom_site.".len()..].trim().to_string());
            continue;
        }
        if in_header {
            in_header = false;
            in_loop = true;
            let names = [
                "group_PDB", "type_symbol", "label_atom_id", "auth_atom_id", "label_alt_id", "label_comp_id",
                "auth_comp_id", "label_asym_id", "auth_asym_id", "label_seq_id", "auth_seq_id",
                "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv", "pdbx_PDB_model_num",
            ];
            map = names.iter().map(|n| idx(&cols, n)).collect();
        }
        if !in_loop {
            continue;
        }
        if line.starts_with('#') || line.starts_with("loop_") || line.starts_with('_') || line.starts_with("data_") {
            if !atoms.is_empty() || line.starts_with('#') {
                break;
            }
            continue;
        }
        pending.extend(cif_tokens(&line));
        if pending.len() < cols.len() {
            continue;
        }
        let tok = std::mem::take(&mut pending);
        let get = |k: usize| -> Option<&str> { map[k].and_then(|c| tok.get(c)).map(|s| s.as_str()) };
        if let Some(model) = get(16) {
            match &first_model {
                None => first_model = Some(model.to_string()),
                Some(m) if m != model => break,
                _ => {}
            }
        }
        let group = get(0).unwrap_or("ATOM");
        let hetatm = group == "HETATM";
        if group != "ATOM" && !hetatm {
            continue;
        }
        let name = get(3).or(get(2)).unwrap_or("");
        let resn = get(6).or(get(5)).unwrap_or("");
        let chain = get(8).or(get(7)).unwrap_or("").to_string();
        let resnum: i32 = match get(10).or(get(9)).and_then(|v| v.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let icode = match get(11) {
            Some(c) if c != "?" && c != "." && !c.is_empty() => c.as_bytes()[0],
            _ => b' ',
        };
        let alt = match get(4) {
            Some(c) if c != "?" && c != "." && !c.is_empty() => c.as_bytes()[0],
            _ => b' ',
        };
        let (Some(x), Some(y), Some(z)) = (
            get(12).and_then(|v| v.parse::<f64>().ok()),
            get(13).and_then(|v| v.parse::<f64>().ok()),
            get(14).and_then(|v| v.parse::<f64>().ok()),
        ) else {
            continue;
        };
        let el = get(1).unwrap_or("").as_bytes();
        let mut element = [b' '; 2];
        if el.len() == 1 {
            element[1] = el[0];
        } else if el.len() >= 2 {
            element = [el[0], el[1].to_ascii_uppercase()];
        }
        // PDB-style atom name alignment: names shorter than 4 start at column 14
        let nb = name.as_bytes();
        let aname = if nb.len() >= 4 || (element[0] != b' ' && nb.len() > 1 && element[0] == nb[0] && element[1] == nb[1]) {
            pad4_atom(nb)
        } else {
            let mut o = [b' '; 4];
            for (k, &c) in nb.iter().take(3).enumerate() {
                o[k + 1] = c;
            }
            o
        };
        atoms.push(RawAtom {
            chain,
            resnum,
            icode,
            resn: pad3(resn.as_bytes()),
            name: aname,
            element,
            alt,
            xyz: [x, y, z],
            bfactor: get(15).and_then(|v| v.parse().ok()).unwrap_or(0.0),
            hetatm,
        });
    }
    Ok(atoms)
}

fn is_cif(path: &Path) -> bool {
    let s = path.to_string_lossy().to_ascii_lowercase();
    let s = s.strip_suffix(".gz").unwrap_or(&s);
    s.ends_with(".cif") || s.ends_with(".mmcif") || s.ends_with(".bcif")
}

fn sniff_cif(path: &Path) -> Result<bool> {
    let mut r = open_text(path)?;
    let mut buf = vec![0u8; 4096];
    let n = r.read(&mut buf)?;
    let head = String::from_utf8_lossy(&buf[..n]);
    Ok(head.trim_start().starts_with("data_") || head.contains("\n_atom_site.") || head.contains("\nloop_"))
}

/// Parse a structure file. `chain`: author chain id to select (None = first
/// protein chain). `keep_atoms`: retain all atoms (needed to write models).
pub fn read_structure(path: &Path, chain: Option<&str>, keep_atoms: bool) -> Result<Structure> {
    let cif = is_cif(path) || (!path.to_string_lossy().contains(".pdb") && sniff_cif(path)?);
    let raw = if cif { parse_cif_atoms(open_text(path)?)? } else { parse_pdb_atoms(open_text(path)?)? };
    let mut name = path.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
    for ext in [".gz", ".pdb", ".cif", ".ent", ".mmcif"] {
        if let Some(stripped) = name.strip_suffix(ext) {
            name = stripped.to_string();
        }
    }
    build(name, raw, chain, keep_atoms).with_context(|| format!("while reading {}", path.display()))
}

fn build(name: String, raw: Vec<RawAtom>, chain: Option<&str>, keep_atoms: bool) -> Result<Structure> {
    // Pick the chain: requested one, or the first chain that holds an amino acid CA.
    let chain_id = match chain {
        Some(c) => c.to_string(),
        None => raw
            .iter()
            .find(|a| &a.name == b" CA " && three_to_one(&a.resn).is_some())
            .map(|a| a.chain.clone())
            .unwrap_or_default(),
    };
    let mut s = Structure { name, chain: chain_id.clone(), ..Default::default() };
    let mut cur: Option<(i32, u8)> = None;
    let mut cur_alt: u8 = b' ';
    let mut res_atoms: Vec<&RawAtom> = Vec::new();
    let flush = |res_atoms: &mut Vec<&RawAtom>, s: &mut Structure| {
        if res_atoms.is_empty() {
            return;
        }
        let a0 = res_atoms[0];
        let aa = three_to_one(&a0.resn);
        let find = |n: &[u8; 4]| res_atoms.iter().find(|a| &a.name == n).map(|a| a.xyz);
        if let (Some(aa), Some(ca)) = (aa, find(b" CA ")) {
            let idx = s.ca.len() as u32;
            s.ca.push(ca);
            s.seq.push(aa);
            s.resn.push(if a0.hetatm && aa != b'X' { mapped_resn(aa) } else { a0.resn });
            s.resid.push(ResId { num: a0.resnum, icode: a0.icode });
            s.bfac.push(res_atoms.iter().find(|a| &a.name == b" CA ").map_or(0.0, |a| a.bfactor));
            s.backbone.push(Backbone { n: find(b" N  "), ca: Some(ca), c: find(b" C  "), o: find(b" O  ") });
            if keep_atoms {
                for a in res_atoms.iter() {
                    let mut nm = a.name;
                    if a.hetatm && &a.resn == b"MSE" && &nm == b"SE  " {
                        nm = *b" SD ";
                    }
                    s.atoms.push(Atom { res: idx, name: nm, element: a.element, xyz: a.xyz, bfactor: a.bfactor });
                }
            }
        }
        res_atoms.clear();
    };
    for a in raw.iter() {
        if a.chain != chain_id {
            continue;
        }
        if a.hetatm && three_to_one(&a.resn).is_none() {
            continue;
        }
        if &a.resn == b"HOH" {
            continue;
        }
        let key = (a.resnum, a.icode);
        if cur != Some(key) {
            flush(&mut res_atoms, &mut s);
            cur = Some(key);
            cur_alt = b' ';
        }
        // alternate locations: keep the first alternate seen for this residue
        if a.alt != b' ' {
            if cur_alt == b' ' {
                cur_alt = a.alt;
            }
            if a.alt != cur_alt {
                continue;
            }
        }
        if res_atoms.iter().any(|b| b.name == a.name) {
            continue;
        }
        res_atoms.push(a);
    }
    flush(&mut res_atoms, &mut s);
    if s.ca.is_empty() {
        bail!("no protein residue with a CA atom found (chain '{}')", chain_id);
    }
    Ok(s)
}

fn mapped_resn(aa: u8) -> [u8; 3] {
    let names: [(&[u8; 3], u8); 20] = [
        (b"ALA", b'A'), (b"ARG", b'R'), (b"ASN", b'N'), (b"ASP", b'D'), (b"CYS", b'C'),
        (b"GLN", b'Q'), (b"GLU", b'E'), (b"GLY", b'G'), (b"HIS", b'H'), (b"ILE", b'I'),
        (b"LEU", b'L'), (b"LYS", b'K'), (b"MET", b'M'), (b"PHE", b'F'), (b"PRO", b'P'),
        (b"SER", b'S'), (b"THR", b'T'), (b"TRP", b'W'), (b"TYR", b'Y'), (b"VAL", b'V'),
    ];
    names.iter().find(|(_, c)| *c == aa).map(|(n, _)| **n).unwrap_or(*b"UNK")
}
