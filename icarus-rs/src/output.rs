//! Result formatting: TSV records, text reports and superposed models.

use std::fmt::Write as _;
use std::io::Write;

use crate::align::{Alignment, PairResult};
use crate::geom::dist;
use crate::prep::Prepared;
use crate::tm;

/// Per-pair summary statistics of an alignment.
pub struct Stats {
    pub tm_min: f64,
    pub tm_q: f64,
    pub tm_t: f64,
    pub n_aligned: usize,
    pub n_core: usize,
    pub rmsd_core: f64,
    pub seq_id: f64,
    /// TM-score (normalised by the longer chain) of the best run of
    /// sequence-consecutive bodies whose junctions preserve chain connectivity.
    pub tm_conn: f64,
    /// Number of bodies in that run.
    pub n_conn: usize,
}

/// Raw kernel sum (d0 of the longer chain) of the best run of connected,
/// sequence-consecutive bodies, and the run length.
fn connected_run(aln: &Alignment, moving: &Prepared, fixed: &Prepared) -> (f64, usize) {
    let lmax = moving.len().max(fixed.len());
    crate::align::connected_run(aln, &moving.s.ca, &fixed.s.ca, 1.0 / tm::d0(lmax).powi(2))
}

/// `moving` is the segmented protein of the alignment, `fixed` the other one.
pub fn stats(aln: &Alignment, moving: &Prepared, fixed: &Prepared, lnorm: usize) -> Stats {
    let (lq, lt) = (moving.len(), fixed.len());
    let (i0, iq, it) = (
        1.0 / tm::d0(lnorm).powi(2),
        1.0 / tm::d0(lq).powi(2),
        1.0 / tm::d0(lt).powi(2),
    );
    let (mut s0, mut sq, mut st, mut rms, mut ncore, mut ident) =
        (0.0, 0.0, 0.0, 0.0, 0usize, 0usize);
    for (k, &(i, j)) in aln.pairs.iter().enumerate() {
        let seg = &aln.segs[aln.pair_seg[k] as usize];
        let x = seg.tr.apply(&moving.s.ca[i as usize]);
        let d = dist(&x, &fixed.s.ca[j as usize]);
        let d2 = d * d;
        s0 += tm::kernel(d2, i0);
        sq += tm::kernel(d2, iq);
        st += tm::kernel(d2, it);
        if d <= 5.0 {
            rms += d2;
            ncore += 1;
            if moving.s.seq[i as usize] == fixed.s.seq[j as usize] {
                ident += 1;
            }
        }
    }
    let (conn_raw, n_conn) = connected_run(aln, moving, fixed);
    Stats {
        tm_conn: conn_raw / lq.max(lt) as f64,
        n_conn,
        tm_min: s0 / lnorm as f64,
        tm_q: sq / lq as f64,
        tm_t: st / lt as f64,
        n_aligned: aln.pairs.len(),
        n_core: ncore,
        rmsd_core: if ncore > 0 {
            (rms / ncore as f64).sqrt()
        } else {
            0.0
        },
        seq_id: if ncore > 0 {
            ident as f64 / ncore as f64
        } else {
            0.0
        },
    }
}

/// Segment description "qstart-qend:tstart-tend" in author numbering.
pub fn segments_string(aln: &Alignment, moving: &Prepared, fixed: &Prepared) -> String {
    let mut out = String::new();
    for &si in &aln.order {
        let sg = &aln.segs[si];
        let tj: Vec<u32> = aln
            .pairs
            .iter()
            .zip(&aln.pair_seg)
            .filter(|(_, &s)| s as usize == si)
            .map(|(p, _)| p.1)
            .collect();
        if !out.is_empty() {
            out.push(',');
        }
        let (a, b) = (moving.s.resid[sg.qs], moving.s.resid[sg.qe]);
        if tj.is_empty() {
            let _ = write!(out, "{a}-{b}:-");
        } else {
            let (lo, hi) = (
                *tj.iter().min().unwrap() as usize,
                *tj.iter().max().unwrap() as usize,
            );
            let _ = write!(out, "{a}-{b}:{}-{}", fixed.s.resid[lo], fixed.s.resid[hi]);
        }
    }
    out
}

pub const TSV_HEADER: &str =
    "query\ttarget\tlen_q\tlen_t\ttm_flex\ttm_flex_q\ttm_flex_t\ttm_rigid\ttm_rigid_max\ttm_conn\tn_conn\tn_bodies\tn_aligned\tn_core\trmsd_core\tseq_id\tpeeled\tbodies";

pub fn tsv_line(r: &PairResult, a: &Prepared, b: &Prepared) -> String {
    let (mv, fx) = if r.reversed { (b, a) } else { (a, b) };
    let st = stats(&r.flex, mv, fx, r.lnorm);
    let (tm_q, tm_t) = if r.reversed {
        (st.tm_t, st.tm_q)
    } else {
        (st.tm_q, st.tm_t)
    };
    // rigid TM-score normalised by the longer chain
    let (rmv, rfx) = if r.rigid_reversed { (b, a) } else { (a, b) };
    let rst = stats(&r.rigid, rmv, rfx, r.lnorm);
    let rigid_max = if rmv.len() >= rfx.len() {
        rst.tm_q
    } else {
        rst.tm_t
    };
    format!(
        "{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{:.2}\t{:.3}\t{}\t{}",
        a.s.name,
        b.s.name,
        a.len(),
        b.len(),
        st.tm_min,
        tm_q,
        tm_t,
        r.tm_rigid(),
        rigid_max,
        r.tm_conn().max(rigid_max),
        r.conn_n,
        r.flex.segs.len(),
        st.n_aligned,
        st.n_core,
        st.rmsd_core,
        st.seq_id,
        if r.reversed { 2 } else { 1 },
        segments_string(&r.flex, mv, fx)
    )
}

/// Per-body superpositions in target order, `;`-separated; each body is
/// "r11,r12,r13,r21,r22,r23,r31,r32,r33,t1,t2,t3" mapping the peeled
/// structure onto the rigid one (x' = R x + t).
pub fn transforms_string(aln: &Alignment) -> String {
    aln.order
        .iter()
        .map(|&si| {
            let tr = &aln.segs[si].tr;
            let mut v: Vec<String> = Vec::with_capacity(12);
            for row in &tr.r {
                for x in row {
                    v.push(format!("{x:.5}"));
                }
            }
            for x in &tr.t {
                v.push(format!("{x:.3}"));
            }
            v.join(",")
        })
        .collect::<Vec<_>>()
        .join(";")
}

/// Human readable report with the chimera alignment.
pub fn report(r: &PairResult, a: &Prepared, b: &Prepared) -> String {
    let (mv, fx) = if r.reversed { (b, a) } else { (a, b) };
    let st = stats(&r.flex, mv, fx, r.lnorm);
    let mut s = String::new();
    let _ = writeln!(s, "ICARUS flexible alignment");
    let _ = writeln!(s, "  Structure 1 : {} ({} residues)", a.s.name, a.len());
    let _ = writeln!(s, "  Structure 2 : {} ({} residues)", b.s.name, b.len());
    let _ = writeln!(
        s,
        "  Peeled (flexible) structure: {}  /  rigid: {}",
        mv.s.name, fx.s.name
    );
    let _ = writeln!(
        s,
        "  TM-score (flexible, norm. by shortest = {}) : {:.4}",
        r.lnorm, st.tm_min
    );
    let _ = writeln!(
        s,
        "  TM-score (flexible, norm. by {} / {})       : {:.4} / {:.4}",
        mv.s.name, fx.s.name, st.tm_q, st.tm_t
    );
    let _ = writeln!(
        s,
        "  TM-score (rigid, norm. by shortest)        : {:.4}",
        r.tm_rigid()
    );
    let _ = writeln!(s, "  Rigid bodies: {}   aligned: {}   within 5 A: {}   RMSD(5 A core): {:.2}   seq. id: {:.1}%",
        r.flex.segs.len(), st.n_aligned, st.n_core, st.rmsd_core, 100.0 * st.seq_id);
    let _ = writeln!(
        s,
        "\n  Rigid bodies in target order ({} residues -> {} residues):",
        mv.s.name, fx.s.name
    );
    for (rank, &si) in r.flex.order.iter().enumerate() {
        let sg = &r.flex.segs[si];
        let tj: Vec<u32> = r
            .flex
            .pairs
            .iter()
            .zip(&r.flex.pair_seg)
            .filter(|(_, &x)| x as usize == si)
            .map(|(p, _)| p.1)
            .collect();
        let tr = if tj.is_empty() {
            "-".to_string()
        } else {
            format!(
                "{}-{}",
                fx.s.resid[*tj.iter().min().unwrap() as usize],
                fx.s.resid[*tj.iter().max().unwrap() as usize]
            )
        };
        let _ = writeln!(
            s,
            "    body {:>2}: {:>6}-{:<6} -> {}",
            rank + 1,
            mv.s.resid[sg.qs].to_string(),
            mv.s.resid[sg.qe].to_string(),
            tr
        );
    }
    // chimera alignment text
    let (mut l1, mut l2, mut l3, mut lb) =
        (String::new(), String::new(), String::new(), String::new());
    let mut pos: std::collections::HashMap<u32, (u32, f64, u16)> = Default::default();
    for (k, &(i, j)) in r.flex.pairs.iter().enumerate() {
        let sg = &r.flex.segs[r.flex.pair_seg[k] as usize];
        let d = dist(&sg.tr.apply(&mv.s.ca[i as usize]), &fx.s.ca[j as usize]);
        pos.insert(i, (j, d, r.flex.pair_seg[k]));
    }
    let mut tnext = 0u32;
    for (rank, &si) in r.flex.order.iter().enumerate() {
        let sg = &r.flex.segs[si];
        let tag = (b'1' + (rank % 9) as u8) as char;
        for i in sg.qs..=sg.qe {
            if let Some(&(j, d, _)) = pos.get(&(i as u32)) {
                while tnext < j {
                    l1.push('-');
                    l2.push(' ');
                    l3.push(fx.s.seq[tnext as usize] as char);
                    lb.push(' ');
                    tnext += 1;
                }
                l1.push(mv.s.seq[i] as char);
                l2.push(if d <= 1.0 {
                    '|'
                } else if d <= 2.0 {
                    ':'
                } else if d <= 4.0 {
                    '.'
                } else {
                    ' '
                });
                l3.push(fx.s.seq[j as usize] as char);
                lb.push(tag);
                tnext = j + 1;
            } else {
                l1.push(mv.s.seq[i] as char);
                l2.push(' ');
                l3.push('-');
                lb.push(tag);
            }
        }
    }
    while (tnext as usize) < fx.len() {
        l1.push('-');
        l2.push(' ');
        l3.push(fx.s.seq[tnext as usize] as char);
        lb.push(' ');
        tnext += 1;
    }
    let _ = writeln!(s, "\n  Alignment (body number, peeled structure, distance class |<=1A :<=2A .<=4A, rigid structure):");
    let chars: Vec<(char, char, char, char)> = lb
        .chars()
        .zip(l1.chars())
        .zip(l2.chars())
        .zip(l3.chars())
        .map(|(((a, b), c), d)| (a, b, c, d))
        .collect();
    for chunk in chars.chunks(80) {
        let _ = writeln!(
            s,
            "  body   {}",
            chunk.iter().map(|c| c.0).collect::<String>()
        );
        let _ = writeln!(
            s,
            "  {:<6} {}",
            "moved",
            chunk.iter().map(|c| c.1).collect::<String>()
        );
        let _ = writeln!(
            s,
            "         {}",
            chunk.iter().map(|c| c.2).collect::<String>()
        );
        let _ = writeln!(
            s,
            "  {:<6} {}\n",
            "fixed",
            chunk.iter().map(|c| c.3).collect::<String>()
        );
    }
    s
}

/// Write the moved (segmented) structure after flexible superposition.
/// `chimera`: residues written in target order (as ICARUS does), otherwise in
/// sequence order. Requires structures read with all atoms.
pub fn write_moved_pdb<W: Write>(
    w: &mut W,
    aln: &Alignment,
    moving: &Prepared,
    chimera: bool,
) -> std::io::Result<()> {
    let s = &moving.s;
    let mut seg_of = vec![usize::MAX; s.len()];
    for (si, sg) in aln.segs.iter().enumerate() {
        for x in seg_of.iter_mut().take(sg.qe + 1).skip(sg.qs) {
            *x = si;
        }
    }
    let mut res_order: Vec<usize> = Vec::with_capacity(s.len());
    if chimera {
        for &si in &aln.order {
            res_order.extend(aln.segs[si].qs..=aln.segs[si].qe);
        }
    } else {
        res_order.extend((0..s.len()).filter(|&i| seg_of[i] != usize::MAX));
    }
    // atom ranges per residue
    let mut first = vec![usize::MAX; s.len()];
    let mut count = vec![0usize; s.len()];
    for (k, a) in s.atoms.iter().enumerate() {
        let r = a.res as usize;
        if first[r] == usize::MAX {
            first[r] = k;
        }
        count[r] += 1;
    }
    let chain = s.chain.chars().next().unwrap_or('A');
    let mut serial = 1;
    for &r in &res_order {
        let tr = &aln.segs[seg_of[r]].tr;
        let atoms: Vec<(String, [u8; 2], [f64; 3], f32)> = if first[r] == usize::MAX {
            vec![(" CA ".to_string(), *b" C", s.ca[r], 0.0)]
        } else {
            s.atoms[first[r]..first[r] + count[r]]
                .iter()
                .map(|a| {
                    (
                        String::from_utf8_lossy(&a.name).to_string(),
                        a.element,
                        a.xyz,
                        a.bfactor,
                    )
                })
                .collect()
        };
        for (name, el, xyz, b) in atoms {
            let x = tr.apply(&xyz);
            writeln!(
                w,
                "ATOM  {:>5} {:<4} {:>3} {}{:>4}{}   {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {}",
                serial % 100000,
                name,
                String::from_utf8_lossy(&s.resn[r]),
                chain,
                s.resid[r].num,
                s.resid[r].icode as char,
                x[0],
                x[1],
                x[2],
                1.0,
                b,
                String::from_utf8_lossy(&el)
            )?;
            serial += 1;
        }
    }
    writeln!(w, "TER")?;
    writeln!(w, "END")
}

/// gdt2.pl-equivalent evaluation of two superposed coordinate sets: optimal
/// sequential alignment without gap penalties (search d0 of the shortest
/// chain), then TM-score normalised by `len` (default: shortest chain).
pub struct GdtResult {
    pub tm: f64,
    pub tm_search: f64,
    pub n_aligned: usize,
    pub pairs: Vec<(u32, u32, f64)>,
}

pub fn gdt(x: &[[f64; 3]], y: &[[f64; 3]], len: Option<usize>) -> GdtResult {
    let size = x.len().min(y.len());
    let d0s = tm::d0_search(size);
    let mut work = crate::dp::DpWork::default();
    let mut pairs = Vec::new();
    let raw_s = work.align_free(x, y, 1.0 / (d0s * d0s), &mut pairs);
    let l = len.unwrap_or(size);
    let d0 = (1.24 * ((l as f64) - 15.0).cbrt() - 1.8).max(0.5);
    let inv = 1.0 / (d0 * d0);
    let mut out = Vec::with_capacity(pairs.len());
    let mut s = 0.0;
    for &(i, j) in &pairs {
        let d = dist(&x[i as usize], &y[j as usize]);
        s += tm::kernel(d * d, inv);
        out.push((i, j, d));
    }
    GdtResult {
        tm: s / l as f64,
        tm_search: raw_s / size as f64,
        n_aligned: pairs.len(),
        pairs: out,
    }
}
