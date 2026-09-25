//! Binary database of preprocessed structures (`.icdb`).
//!
//! Stores everything a comparison needs (sequence, residue ids, C-alpha
//! coordinates, secondary structure and the Protein Unit hierarchy) so that
//! structures are parsed, DSSP-assigned and peeled only once. Fragment
//! signatures are cheap and recomputed on load.
//!
//! Layout (little endian): magic "ICDB", u32 version, u32 count, then records.

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::prep::{fragment_signatures, Prepared, PuNode, PuTree};
use crate::structure::{ResId, Structure};

const MAGIC: &[u8; 4] = b"ICDB";
const VERSION: u32 = 1;

fn w_u8<W: Write>(w: &mut W, v: u8) -> std::io::Result<()> {
    w.write_all(&[v])
}
fn w_u16<W: Write>(w: &mut W, v: u16) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
fn w_u32<W: Write>(w: &mut W, v: u32) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
fn w_i32<W: Write>(w: &mut W, v: i32) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
fn w_f32<W: Write>(w: &mut W, v: f32) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
fn w_str<W: Write>(w: &mut W, s: &str) -> std::io::Result<()> {
    w_u16(w, s.len() as u16)?;
    w.write_all(s.as_bytes())
}

struct Rd<'a> {
    b: &'a [u8],
    p: usize,
}

impl<'a> Rd<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.p + n > self.b.len() {
            bail!("truncated database");
        }
        let s = &self.b[self.p..self.p + n];
        self.p += n;
        Ok(s)
    }
    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into()?))
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into()?))
    }
    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into()?))
    }
    fn f32(&mut self) -> Result<f32> {
        Ok(f32::from_le_bytes(self.take(4)?.try_into()?))
    }
    fn str(&mut self) -> Result<String> {
        let n = self.u16()? as usize;
        Ok(String::from_utf8_lossy(self.take(n)?).to_string())
    }
}

/// Streaming `.icdb` writer: records are appended one at a time and the
/// record count in the header is patched by `finish`, so large collections
/// never need to be held in memory at once.
pub struct DbWriter {
    w: BufWriter<File>,
    count: u32,
}

impl DbWriter {
    pub fn create(path: &Path) -> Result<Self> {
        let mut w = BufWriter::new(
            File::create(path).with_context(|| format!("cannot create {}", path.display()))?,
        );
        w.write_all(MAGIC)?;
        w_u32(&mut w, VERSION)?;
        w_u32(&mut w, 0)?;
        Ok(Self { w, count: 0 })
    }

    pub fn push(&mut self, p: &Prepared) -> Result<()> {
        write_record(&mut self.w, p)?;
        self.count += 1;
        Ok(())
    }

    /// Writes the record count into the header and flushes; returns the count.
    pub fn finish(self) -> Result<u32> {
        let mut f = self.w.into_inner().map_err(|e| e.into_error())?;
        f.seek(SeekFrom::Start(8))?;
        f.write_all(&self.count.to_le_bytes())?;
        f.flush()?;
        Ok(self.count)
    }
}

pub fn write_db(path: &Path, items: &[Prepared]) -> Result<()> {
    let mut w = DbWriter::create(path)?;
    for p in items {
        w.push(p)?;
    }
    w.finish()?;
    Ok(())
}

fn write_record<W: Write>(w: &mut W, p: &Prepared) -> Result<()> {
    let s = &p.s;
    w_str(w, &s.name)?;
    w_str(w, &s.chain)?;
    w_u32(w, s.len() as u32)?;
    w.write_all(&s.seq)?;
    for r in &s.resn {
        w.write_all(r)?;
    }
    for r in &s.resid {
        w_i32(w, r.num)?;
        w_u8(w, r.icode)?;
    }
    for c in &s.ca {
        for &v in c {
            w_f32(w, v as f32)?;
        }
    }
    for &b in &s.bfac {
        w_f32(w, b)?;
    }
    for &c in &p.ss {
        w_u8(w, c as u8)?;
    }
    let t = &p.tree;
    w_u16(w, t.nodes.len() as u16)?;
    for n in &t.nodes {
        w_u32(w, n.start as u32)?;
        w_u32(w, n.end as u32)?;
        w_u16(w, n.depth as u16)?;
        w_i32(w, n.parent.map_or(-1, |x| x as i32))?;
        w_u32(w, n.leaf_mask)?;
    }
    w_u16(w, t.leaves.len() as u16)?;
    for &l in &t.leaves {
        w_u16(w, l as u16)?;
    }
    w_u16(w, t.levels.len() as u16)?;
    for lv in &t.levels {
        w_u16(w, lv.len() as u16)?;
        for &x in lv {
            w_u16(w, x as u16)?;
        }
    }
    Ok(())
}

pub fn read_db(path: &Path) -> Result<Vec<Prepared>> {
    let mut buf = Vec::new();
    BufReader::new(File::open(path).with_context(|| format!("cannot open {}", path.display()))?)
        .read_to_end(&mut buf)?;
    let mut r = Rd { b: &buf, p: 0 };
    if r.take(4)? != MAGIC {
        bail!("{} is not an ICARUS database", path.display());
    }
    let version = r.u32()?;
    if version != VERSION {
        bail!("unsupported database version {version}");
    }
    let count = r.u32()? as usize;
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let name = r.str()?;
        let chain = r.str()?;
        let n = r.u32()? as usize;
        let seq = r.take(n)?.to_vec();
        let mut resn = Vec::with_capacity(n);
        for _ in 0..n {
            let b = r.take(3)?;
            resn.push([b[0], b[1], b[2]]);
        }
        let mut resid = Vec::with_capacity(n);
        for _ in 0..n {
            let num = r.i32()?;
            let icode = r.u8()?;
            resid.push(ResId { num, icode });
        }
        let mut ca = Vec::with_capacity(n);
        for _ in 0..n {
            ca.push([r.f32()? as f64, r.f32()? as f64, r.f32()? as f64]);
        }
        let mut bfac = Vec::with_capacity(n);
        for _ in 0..n {
            bfac.push(r.f32()?);
        }
        let ss: Vec<char> = r.take(n)?.iter().map(|&c| c as char).collect();
        let nn = r.u16()? as usize;
        let mut nodes = Vec::with_capacity(nn);
        for _ in 0..nn {
            let start = r.u32()? as usize;
            let end = r.u32()? as usize;
            let depth = r.u16()? as usize;
            let parent = r.i32()?;
            let leaf_mask = r.u32()?;
            nodes.push(PuNode {
                start,
                end,
                depth,
                parent: (parent >= 0).then_some(parent as usize),
                children: vec![],
                leaf_mask,
            });
        }
        for i in 0..nn {
            if let Some(p) = nodes[i].parent {
                nodes[p].children.push(i);
            }
        }
        let nl = r.u16()? as usize;
        let mut leaves = Vec::with_capacity(nl);
        for _ in 0..nl {
            leaves.push(r.u16()? as usize);
        }
        let nlv = r.u16()? as usize;
        let mut levels = Vec::with_capacity(nlv);
        for _ in 0..nlv {
            let k = r.u16()? as usize;
            let mut lv = Vec::with_capacity(k);
            for _ in 0..k {
                lv.push(r.u16()? as usize);
            }
            levels.push(lv);
        }
        let frags = fragment_signatures(&ca);
        let s = Structure {
            name,
            chain,
            seq,
            resn,
            ca,
            resid,
            bfac,
            backbone: vec![],
            atoms: vec![],
        };
        out.push(Prepared {
            s,
            ss,
            tree: PuTree {
                nodes,
                leaves,
                levels,
            },
            frags,
        });
    }
    Ok(out)
}
