//! Seed superpositions from aligned fragment pairs (AFPs).
//!
//! Two 8-residue fragments form an AFP when their intra-fragment C-alpha
//! distance matrices agree (dRMS below a tolerance), which is invariant to
//! rigid motion and needs no superposition. AFPs on the same diagonal are
//! merged into gapless blocks; each block yields one candidate transform.

use crate::geom::{superpose, Transform, V3};
use crate::prep::{FRAG, FRAG_D};

#[derive(Debug, Clone, Copy)]
pub struct Block {
    pub qi: u32,
    pub tj: u32,
    pub len: u32,
}

/// Distances used for the quick rejection test (long-range pairs of the fragment).
/// Packed indices of pairs (0,7), (0,4), (3,7), (2,6).
const QUICK: [usize; 4] = [6, 3, 21, 16];

#[inline]
fn quick_reject(a: &[f32; FRAG_D], b: &[f32; FRAG_D], tol: f32) -> bool {
    QUICK.iter().any(|&k| (a[k] - b[k]).abs() > tol)
}

#[inline]
fn drms2(a: &[f32; FRAG_D], b: &[f32; FRAG_D]) -> f32 {
    let mut s = 0.0f32;
    for k in 0..FRAG_D {
        let d = a[k] - b[k];
        s += d * d;
    }
    s / FRAG_D as f32
}

/// Find gapless blocks of similar local structure between query and target.
pub fn find_blocks(
    qf: &[[f32; FRAG_D]],
    tf: &[[f32; FRAG_D]],
    q_stride: usize,
    tol: f32,
) -> Vec<Block> {
    let tol2 = tol * tol;
    let quick_tol = 3.0 * tol;
    let mut afps: Vec<(i32, u32)> = Vec::new(); // (diagonal j - i, i)
    let stride = q_stride.max(1);
    for i in (0..qf.len()).step_by(stride) {
        let a = &qf[i];
        for (j, b) in tf.iter().enumerate() {
            if quick_reject(a, b, quick_tol) {
                continue;
            }
            if drms2(a, b) <= tol2 {
                afps.push((j as i32 - i as i32, i as u32));
            }
        }
    }
    afps.sort_unstable();
    let mut blocks = Vec::new();
    let mut k = 0;
    while k < afps.len() {
        let (diag, i0) = afps[k];
        let mut last = i0;
        let mut k2 = k + 1;
        while k2 < afps.len() && afps[k2].0 == diag && afps[k2].1 <= last + stride as u32 {
            last = afps[k2].1;
            k2 += 1;
        }
        let len = last - i0 + FRAG as u32;
        blocks.push(Block {
            qi: i0,
            tj: (i0 as i32 + diag) as u32,
            len,
        });
        k = k2;
    }
    blocks
}

/// Keep at most `max` blocks, longest first, with a quota per query region so
/// that every part of the query keeps some seeds.
pub fn select_blocks(mut blocks: Vec<Block>, qlen: usize, max: usize) -> Vec<Block> {
    if blocks.len() <= max {
        return blocks;
    }
    blocks.sort_unstable_by(|a, b| b.len.cmp(&a.len).then(a.qi.cmp(&b.qi)));
    let bin = 10usize;
    let nbins = qlen / bin + 1;
    let quota = ((max * 2) / nbins).max(8);
    let mut count = vec![0usize; nbins];
    let mut out = Vec::with_capacity(max);
    for b in blocks {
        let c = &mut count[b.qi as usize / bin];
        if *c < quota {
            *c += 1;
            out.push(b);
            if out.len() >= max {
                break;
            }
        }
    }
    out
}

/// Superposition of the block residues (query onto target).
pub fn block_transform(q: &[V3], t: &[V3], b: &Block) -> Transform {
    let (qi, tj, l) = (b.qi as usize, b.tj as usize, b.len as usize);
    superpose(&q[qi..qi + l], &t[tj..tj + l], None).0
}

#[cfg(test)]
mod tests {
    #[test]
    fn quick_indices_are_long_range_pairs() {
        // packed index of pair (a,b), a<b, in an 8-residue fragment
        let idx =
            |a: usize, b: usize| (0..a).map(|k| super::FRAG - 1 - k).sum::<usize>() + (b - a - 1);
        assert_eq!(super::QUICK, [idx(0, 7), idx(0, 4), idx(3, 7), idx(2, 6)]);
        assert_eq!(idx(6, 7), super::FRAG_D - 1);
    }
}
