//! Flexible structural alignment based on Protein Units.
//!
//! For a query whose PU hierarchy is known and a rigid target:
//!  1. candidate superpositions come from aligned fragment blocks (seeds.rs);
//!  2. every candidate is scored on every PU of the hierarchy at once, using
//!     a distance grid of the target and prefix sums over the query;
//!  3. the candidates of each PU are polished by a cheap closest-point
//!     iteration on the grid, and the best diverse ones are refined by
//!     windowed DP + iteratively reweighted superposition (which increases the
//!     TM objective monotonically for fixed residue pairs); placements of a PU
//!     are also tried for its children;
//!  4. an exact dynamic programme picks a set of disjoint PUs (from any level
//!     of the hierarchy) with non-overlapping target regions, maximising the
//!     summed score minus a per-hinge penalty;
//!  5. the chosen rigid bodies are merged into a chimera ordered along the
//!     target, aligned with the gdt2.pl DP used by ICARUS to score solutions,
//!     and each body is re-superposed on its aligned pairs until convergence;
//!  6. each body is then re-placed against the target residues left free by
//!     the others (the ICARUS "updated target" idea) when that improves the
//!     total score.
//!
//! The procedure runs in both directions (each protein peeled in turn) and the
//! best-scoring solution is kept, as in ICARUS.

use crate::dp::{DpWork, SoA};
use crate::geom::{dist2, superpose, Transform, V3};
use crate::grid::{DistGrid, KernelLut, FAR};
use crate::prep::Prepared;
use crate::seeds::{block_transform, find_blocks, select_blocks};
use crate::tm;

/// Parameters of the pairwise alignment.
#[derive(Debug, Clone)]
pub struct AlignParams {
    /// Maximum number of rigid bodies (PUs) in a solution.
    pub max_segments: usize,
    /// Score penalty per additional rigid body, in TM-score units.
    pub hinge_penalty: f64,
    /// Stride over query fragments when searching seeds.
    pub q_stride: usize,
    /// Maximum number of seed superpositions evaluated per direction.
    pub max_seeds: usize,
    /// Base number of candidate placements refined per PU (small PUs get more).
    pub per_node: usize,
    /// Fragment dRMS tolerance (Å) for aligned fragment pairs.
    pub frag_tol: f32,
    /// Peel and align in both directions (the ICARUS default).
    pub both_directions: bool,
    /// Run the refit stage (re-placing bodies on the free target).
    pub refit: bool,
}

impl Default for AlignParams {
    fn default() -> Self {
        Self {
            max_segments: 6,
            hinge_penalty: 0.0,
            q_stride: 2,
            max_seeds: 1500,
            per_node: 3,
            frag_tol: 1.0,
            both_directions: true,
            refit: true,
        }
    }
}

/// One rigid body of a flexible alignment: query residues qs..=qe moved by `tr`.
#[derive(Debug, Clone)]
pub struct Segment {
    pub qs: usize,
    pub qe: usize,
    pub tr: Transform,
}

/// An alignment of a (possibly segmented) query onto a rigid target.
#[derive(Debug, Clone, Default)]
pub struct Alignment {
    pub segs: Vec<Segment>,
    /// Order of the segments along the target (chimera order).
    pub order: Vec<usize>,
    /// Aligned (query, target) residue indices, in chimera order.
    pub pairs: Vec<(u32, u32)>,
    /// Segment of each aligned pair.
    pub pair_seg: Vec<u16>,
    /// Sum of TM kernels over aligned pairs (d0 of the normalisation length).
    pub raw: f64,
}

/// Result of aligning two structures.
#[derive(Debug, Clone)]
pub struct PairResult {
    /// The flexible solution. If `reversed`, the segmented ("query") protein is
    /// the second input and residue indices in `flex` refer to (second, first).
    pub flex: Alignment,
    pub rigid: Alignment,
    pub reversed: bool,
    pub rigid_reversed: bool,
    pub len1: usize,
    pub len2: usize,
    pub lnorm: usize,
    /// Best connected run over all candidate solutions (see `ConnRun`).
    pub conn: ConnRun,
    /// The connected run comes from the second structure being peeled.
    pub conn_reversed: bool,
}

/// The best run of sequence-consecutive, chain-connected bodies found among
/// the candidate solutions of one direction. Runs are compared by their raw
/// score minus the hinge penalty per extra body, so with a penalty the run is
/// as parsimonious as the reported flexible solution.
#[derive(Debug, Clone, Default)]
pub struct ConnRun {
    /// Raw kernel sum (d0 of the longer chain) of the run.
    pub raw: f64,
    /// `raw` minus the hinge penalty per extra body (selection criterion).
    pub score: f64,
    /// The solution the run belongs to and the indices of its bodies (in
    /// sequence order).
    pub aln: Alignment,
    pub segs: Vec<usize>,
}

impl PairResult {
    /// Connectivity-aware flexible TM-score, normalised by the longer chain
    /// (a detection score: never below the rigid TM-score of the longer chain).
    pub fn tm_conn(&self) -> f64 {
        self.conn.raw / self.len1.max(self.len2) as f64
    }
    pub fn tm_flex(&self) -> f64 {
        self.flex.raw / self.lnorm as f64
    }
    pub fn tm_rigid(&self) -> f64 {
        self.rigid.raw / self.lnorm as f64
    }
}

/// Refinement iterations (banded DP + reweighted superposition) per placement.
const REFINE_ITERS: usize = 4;
/// Stop refining a placement when the relative gain falls below this.
const REFINE_GAIN_STOP: f64 = 0.002;
/// Placements whose DP window has at most this many cells use a full DP.
const FULL_DP_CELLS: usize = 6000;
/// Half-width of the first-pass band (centred on the LIS of nearest residues).
const BAND_FIRST: usize = 20;
/// Squared distance (Å²) below which two candidate superpositions of a PU are
/// considered the same (3 reference residues all closer than this).
const SAME_CANDIDATE_D2: f64 = 9.0;
/// Placements scoring below this fraction of the best one of their PU are
/// not offered to the assembly.
const PRUNE_FRACTION: f64 = 0.5;
/// Chimera alignment / re-superposition cycles in finalisation.
const FINALIZE_ITERS: usize = 3;
/// Refit rounds and number of candidates DP-refined per body in a round.
const REFIT_ROUNDS: usize = 2;
const REFIT_CANDIDATES: usize = 3;
/// Candidate assemblies (best per body count) fully finalised / refitted.
const FINALIZED_SOLUTIONS: usize = 2;
const REFITTED_SOLUTIONS: usize = 1;

struct Ctx<'a> {
    q: &'a Prepared,
    t: &'a Prepared,
    tsoa: SoA,
    inv_d02: f64,
    inv_d02s: f64,
    grid: DistGrid,
    lut: KernelLut,
}

#[derive(Debug, Clone)]
struct Placement {
    node: usize,
    tr: Transform,
    raw: f64,
    /// Core target span (well superposed pairs).
    ta: u32,
    tb: u32,
}

#[inline]
fn cmp_desc(a: f64, b: f64) -> std::cmp::Ordering {
    b.partial_cmp(&a).unwrap_or(std::cmp::Ordering::Equal)
}

/// Iteratively reweighted superposition on fixed pairs: maximises
/// sum 1/(1+d²/d0²) (minorise-maximise with weights 1/(1+d²/d0²)²).
fn irls(
    q: &[V3],
    t: &[V3],
    pairs: &[(u32, u32)],
    tr0: Transform,
    inv_d02: f64,
    iters: usize,
) -> (Transform, f64) {
    let score = |tr: &Transform| -> f64 {
        pairs
            .iter()
            .map(|&(i, j)| tm::kernel(dist2(&tr.apply(&q[i as usize]), &t[j as usize]), inv_d02))
            .sum()
    };
    let mut best_tr = tr0;
    let mut best = score(&tr0);
    if pairs.len() < 3 {
        return (best_tr, best);
    }
    let xs: Vec<V3> = pairs.iter().map(|&(i, _)| q[i as usize]).collect();
    let ys: Vec<V3> = pairs.iter().map(|&(_, j)| t[j as usize]).collect();
    let mut w = vec![0.0; pairs.len()];
    let mut tr = tr0;
    for _ in 0..iters {
        for (k, wk) in w.iter_mut().enumerate() {
            let k1 = tm::kernel(dist2(&tr.apply(&xs[k]), &ys[k]), inv_d02);
            *wk = k1 * k1;
        }
        let (nt, _) = superpose(&xs, &ys, Some(&w));
        let s = score(&nt);
        tr = nt;
        if s > best + 1e-9 {
            best = s;
            best_tr = nt;
        } else {
            break;
        }
    }
    (best_tr, best)
}

/// Closest-point polishing of a candidate superposition of residues s..=e,
/// using the target grid (no sequence constraint). `occ` marks target
/// residues that may not be used. Returns (grid score, transform).
fn quick_refine(
    ctx: &Ctx,
    s: usize,
    e: usize,
    tr0: Transform,
    occ: Option<&[bool]>,
) -> (f32, Transform) {
    let q = &ctx.q.s.ca;
    let t = &ctx.t.s.ca;
    let score_of = |tr: &Transform| -> f32 {
        let mut sc = 0f32;
        for x in &q[s..=e] {
            let (d, j) = ctx.grid.lookup(&tr.apply(x));
            if d != FAR && occ.is_none_or(|o| !o[j as usize]) {
                sc += ctx.lut.v[d as usize];
            }
        }
        sc
    };
    let mut tr = tr0;
    let mut best = (score_of(&tr0), tr0);
    let mut xs = Vec::with_capacity(e - s + 1);
    let mut ys = Vec::with_capacity(e - s + 1);
    let mut ws = Vec::with_capacity(e - s + 1);
    for _ in 0..3 {
        xs.clear();
        ys.clear();
        ws.clear();
        for x in &q[s..=e] {
            let (d, j) = ctx.grid.lookup(&tr.apply(x));
            if d != FAR && d < 60 && occ.is_none_or(|o| !o[j as usize]) {
                let k = ctx.lut.v[d as usize] as f64;
                xs.push(*x);
                ys.push(t[j as usize]);
                ws.push(k * k);
            }
        }
        if xs.len() < 4 {
            break;
        }
        let (nt, _) = superpose(&xs, &ys, Some(&ws));
        let sc = score_of(&nt);
        tr = nt;
        if sc > best.0 + 1e-3 {
            best = (sc, nt);
        } else {
            break;
        }
    }
    best
}

fn core_span(q: &[V3], t: &[V3], pairs: &[(u32, u32)], tr: &Transform) -> (u32, u32) {
    let mut lo = u32::MAX;
    let mut hi = 0;
    for &(i, j) in pairs {
        if dist2(&tr.apply(&q[i as usize]), &t[j as usize]) < 16.0 {
            lo = lo.min(j);
            hi = hi.max(j);
        }
    }
    if lo == u32::MAX {
        return (u32::MAX, 0);
    }
    // tolerate small overlaps between neighbouring bodies
    if hi - lo >= 10 {
        (lo + 2, hi - 2)
    } else {
        (lo, hi)
    }
}

/// Fill undefined (-1) band centres by linear interpolation/extrapolation.
fn fill_centers(c: &mut [i32]) -> bool {
    let def: Vec<usize> = (0..c.len()).filter(|&i| c[i] >= 0).collect();
    if def.is_empty() {
        return false;
    }
    let (f, l) = (def[0], *def.last().unwrap());
    for i in 0..f {
        c[i] = c[f] - (f - i) as i32;
    }
    for i in l + 1..c.len() {
        c[i] = c[l] + (i - l) as i32;
    }
    for w in def.windows(2) {
        let (a, b) = (w[0], w[1]);
        for i in a + 1..b {
            c[i] = c[a] + ((c[b] - c[a]) as f64 * (i - a) as f64 / (b - a) as f64).round() as i32;
        }
    }
    true
}

/// Keep only a longest strictly increasing subsequence of the defined (>= 0)
/// centres (sequence-consistent nearest-residue map); others become -1.
fn lis_centers(c: &mut [i32]) {
    let idx: Vec<usize> = (0..c.len()).filter(|&i| c[i] >= 0).collect();
    if idx.is_empty() {
        return;
    }
    // patience sorting with predecessor links
    let mut tails: Vec<usize> = Vec::new(); // positions in idx
    let mut prev = vec![usize::MAX; idx.len()];
    for (k, &i) in idx.iter().enumerate() {
        let v = c[i];
        let pos = tails.partition_point(|&t| c[idx[t]] < v);
        if pos > 0 {
            prev[k] = tails[pos - 1];
        }
        if pos == tails.len() {
            tails.push(k);
        } else {
            tails[pos] = k;
        }
    }
    let mut keep = vec![false; c.len()];
    let mut k = *tails.last().unwrap();
    loop {
        keep[idx[k]] = true;
        if prev[k] == usize::MAX {
            break;
        }
        k = prev[k];
    }
    for (i, v) in c.iter_mut().enumerate() {
        if !keep[i] {
            *v = -1;
        }
    }
}

/// Refine a candidate superposition of residues s..=e onto the target
/// (coordinates `t`, `tsoa`: possibly masked copies of the target).
#[allow(clippy::too_many_arguments)]
fn refine_placement(
    ctx: &Ctx,
    node: usize,
    s: usize,
    e: usize,
    tr0: Transform,
    t: &[V3],
    tsoa: &SoA,
    work: &mut DpWork,
) -> Option<Placement> {
    const BAND: usize = 14;
    let q = &ctx.q.s.ca;
    let n = t.len();
    let len = e - s + 1;
    let mut tr = tr0;
    let mut best: Option<(Transform, f64, Vec<(u32, u32)>)> = None;
    let mut pairs = Vec::new();
    let mut x = Vec::with_capacity(len);
    let mut centers = vec![-1i32; len];
    for iter in 0..REFINE_ITERS {
        x.clear();
        x.extend(q[s..=e].iter().map(|p| tr.apply(p)));
        if iter == 0 || best.is_none() {
            let (mut jmin, mut jmax) = (u32::MAX, 0u32);
            for (k, p) in x.iter().enumerate() {
                let (d, j) = ctx.grid.lookup(p);
                if d != FAR && d < 70 {
                    jmin = jmin.min(j);
                    jmax = jmax.max(j);
                    centers[k] = j as i32;
                } else {
                    centers[k] = -1;
                }
            }
            if jmin == u32::MAX {
                break;
            }
            let pad = 8usize;
            let w0 = (jmin as usize).saturating_sub(pad);
            let w1 = (jmax as usize + pad).min(n - 1);
            if len * (w1 - w0 + 1) <= FULL_DP_CELLS {
                work.align_gap_open_f(&x, tsoa, w0, w1 + 1, ctx.inv_d02s as f32, -0.6, &mut pairs);
            } else {
                lis_centers(&mut centers);
                if !fill_centers(&mut centers) {
                    break;
                }
                work.align_band_f(
                    &x,
                    tsoa,
                    &centers,
                    BAND_FIRST,
                    ctx.inv_d02s as f32,
                    -0.6,
                    &mut pairs,
                );
            }
        } else {
            // band around the previous alignment path
            let Some(prev) = best.as_ref() else { break };
            centers.iter_mut().for_each(|c| *c = -1);
            for &(i, j) in &prev.2 {
                centers[i as usize - s] = j as i32;
            }
            if !fill_centers(&mut centers) {
                break;
            }
            work.align_band_f(
                &x,
                tsoa,
                &centers,
                BAND,
                ctx.inv_d02s as f32,
                -0.6,
                &mut pairs,
            );
        }
        if pairs.len() < 3 {
            break;
        }
        for p in pairs.iter_mut() {
            p.0 += s as u32;
        }
        let (ntr, raw) = irls(q, t, &pairs, tr, ctx.inv_d02, 3);
        let gain = best.as_ref().map_or(f64::INFINITY, |b| raw - b.1);
        tr = ntr;
        if gain > 1e-6 {
            best = Some((ntr, raw, pairs.clone()));
        }
        if gain < REFINE_GAIN_STOP * raw.max(1.0) {
            break;
        }
    }
    let (tr, raw, pairs) = best?;
    let (ta, tb) = core_span(q, t, &pairs, &tr);
    if ta == u32::MAX {
        return None;
    }
    Some(Placement {
        node,
        tr,
        raw,
        ta,
        tb,
    })
}

/// Exact selection of disjoint PUs with non-overlapping target spans.
/// Returns indices into `cands`.
/// Returns, for each number of bodies k = 1..=max_segments that admits a
/// solution, the best selection (indices into `cands`) with its proxy score.
fn assemble(
    cands: &[Placement],
    tree_masks: &[u32],
    n_leaves: usize,
    max_segments: usize,
    penalty: f64,
    f: &mut Vec<f32>,
) -> Vec<(f64, Vec<usize>)> {
    let nc = cands.len();
    if nc == 0 {
        return vec![];
    }
    let mut order: Vec<usize> = (0..nc).collect();
    order.sort_by_key(|&c| (cands[c].tb, cands[c].ta));
    let ends: Vec<u32> = order.iter().map(|&c| cands[c].tb).collect();
    // prev[k] = number of sorted candidates whose end is < start of candidate k
    let prev: Vec<usize> = order
        .iter()
        .map(|&c| ends.partition_point(|&b| b < cands[c].ta))
        .collect();
    let nm = 1usize << n_leaves;
    let ks = max_segments.max(1);
    let neg = f32::NEG_INFINITY;
    let stride_k = nm;
    let stride_i = (ks + 1) * nm;
    // row 0 initialised here; row i is copied from row i-1 before use
    f.clear();
    f.resize(stride_i, neg);
    f.resize((nc + 1) * stride_i, 0.0);
    f[0] = 0.0;
    let penalty = penalty as f32;
    for i in 1..=nc {
        let c = order[i - 1];
        let mc = tree_masks[cands[c].node] as usize;
        let w = cands[c].raw as f32;
        let p = prev[i - 1];
        let (before, cur) = f.split_at_mut(i * stride_i);
        cur[..stride_i].copy_from_slice(&before[(i - 1) * stride_i..i * stride_i]);
        for k in 1..=ks {
            let pen = if k >= 2 { penalty } else { 0.0 };
            let src = &before[p * stride_i + (k - 1) * stride_k..p * stride_i + k * stride_k];
            let dst = &mut cur[k * stride_k..(k + 1) * stride_k];
            // iterate over masks containing mc: sub = mask ^ mc ranges over masks disjoint from mc
            let free = (nm - 1) & !mc;
            let mut sub = free;
            loop {
                let v = src[sub];
                if v != neg {
                    let nv = v + w - pen;
                    let slot = &mut dst[sub | mc];
                    if nv > *slot {
                        *slot = nv;
                    }
                }
                if sub == 0 {
                    break;
                }
                sub = (sub - 1) & free;
            }
        }
    }
    let base = nc * stride_i;
    let mut out = Vec::new();
    for kk in 1..=ks {
        let (mut bm, mut bv) = (usize::MAX, 0.0f32);
        for mask in 0..nm {
            let v = f[base + kk * stride_k + mask];
            if v > bv {
                bv = v;
                bm = mask;
            }
        }
        if bm == usize::MAX {
            continue;
        }
        let mut chosen = Vec::new();
        let (mut i, mut k, mut mask) = (nc, kk, bm);
        while i > 0 && k > 0 {
            let v = f[i * stride_i + k * stride_k + mask];
            if v == f[(i - 1) * stride_i + k * stride_k + mask] {
                i -= 1;
                continue;
            }
            let c = order[i - 1];
            chosen.push(c);
            mask ^= tree_masks[cands[c].node] as usize;
            k -= 1;
            i = prev[i - 1];
        }
        out.push((bv as f64, chosen));
    }
    out
}

/// One pass: chimera of the bodies in target order, aligned with the gdt2 DP.
fn chimera_align(ctx: &Ctx, segs: &[Segment], keys: &[f64], work: &mut DpWork) -> Alignment {
    let q = &ctx.q.s.ca;
    let t = &ctx.t.s.ca;
    let mut order: Vec<usize> = (0..segs.len()).collect();
    order.sort_by(|&a, &b| {
        keys[a]
            .partial_cmp(&keys[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut chim: Vec<V3> = Vec::with_capacity(q.len());
    let mut origin: Vec<(u32, u16)> = Vec::with_capacity(q.len());
    for &si in &order {
        let sg = &segs[si];
        for i in sg.qs..=sg.qe {
            chim.push(sg.tr.apply(&q[i]));
            origin.push((i as u32, si as u16));
        }
    }
    let mut pairs = Vec::new();
    work.align_free_f(&chim, &ctx.tsoa, ctx.inv_d02s as f32, &mut pairs);
    let mut raw = 0.0;
    let mut out_pairs = Vec::with_capacity(pairs.len());
    let mut out_seg = Vec::with_capacity(pairs.len());
    for &(ci, tj) in &pairs {
        raw += tm::kernel(dist2(&chim[ci as usize], &t[tj as usize]), ctx.inv_d02);
        out_pairs.push((origin[ci as usize].0, tj));
        out_seg.push(origin[ci as usize].1);
    }
    Alignment {
        segs: segs.to_vec(),
        order,
        pairs: out_pairs,
        pair_seg: out_seg,
        raw,
    }
}

/// Chimera alignment, then re-superpose each body on its own pairs; iterate
/// until the score stops improving.
fn finalize(ctx: &Ctx, mut segs: Vec<Segment>, mut keys: Vec<f64>, work: &mut DpWork) -> Alignment {
    let q = &ctx.q.s.ca;
    let t = &ctx.t.s.ca;
    let mut best = Alignment::default();
    for iter in 0..FINALIZE_ITERS {
        let aln = chimera_align(ctx, &segs, &keys, work);
        if aln.raw <= best.raw * (1.0 + 1e-5) + 1e-9 && iter > 0 {
            break;
        }
        best = aln;
        for (si, sg) in segs.iter_mut().enumerate() {
            let sp: Vec<(u32, u32)> = best
                .pairs
                .iter()
                .zip(&best.pair_seg)
                .filter(|(_, &s)| s as usize == si)
                .map(|(p, _)| *p)
                .collect();
            if sp.len() >= 3 {
                let (ntr, _) = irls(q, t, &sp, sg.tr, ctx.inv_d02, 3);
                sg.tr = ntr;
            }
        }
        keys = segment_keys_from(ctx, &segs, &best);
    }
    best
}

/// Move the boundary between bodies that are consecutive in sequence when the
/// residues next to the hinge fit their aligned partners better with the
/// neighbouring body's superposition; confirmed by a full chimera alignment.
fn refine_hinges(ctx: &Ctx, aln: Alignment, work: &mut DpWork) -> Alignment {
    const DMAX: usize = 12;
    const MIN_LEN: usize = 8;
    let q = &ctx.q.s.ca;
    let t = &ctx.t.s.ca;
    if aln.segs.len() < 2 {
        return aln;
    }
    let mut partner = vec![u32::MAX; q.len()];
    for &(i, j) in &aln.pairs {
        partner[i as usize] = j;
    }
    let k = |tr: &Transform, i: usize| -> f64 {
        let j = partner[i];
        if j == u32::MAX {
            0.0
        } else {
            tm::kernel(dist2(&tr.apply(&q[i]), &t[j as usize]), ctx.inv_d02)
        }
    };
    let mut segs = aln.segs.clone();
    let mut by_start: Vec<usize> = (0..segs.len()).collect();
    by_start.sort_by_key(|&s| segs[s].qs);
    let mut changed = false;
    for w in by_start.windows(2) {
        let (a, b) = (w[0], w[1]);
        if segs[a].qe + 1 != segs[b].qs {
            continue;
        }
        let (ta, tb) = (segs[a].tr, segs[b].tr);
        let mut best = (0.3f64, 0isize);
        // shift > 0: residues move from b to a; shift < 0: from a to b
        let mut gain = 0.0;
        for d in 1..=DMAX {
            let i = segs[b].qs + d - 1;
            if segs[b].qe + 1 < segs[b].qs + d + MIN_LEN {
                break;
            }
            gain += k(&ta, i) - k(&tb, i);
            if gain > best.0 {
                best = (gain, d as isize);
            }
        }
        gain = 0.0;
        for d in 1..=DMAX {
            if segs[a].qe + 1 < segs[a].qs + d + MIN_LEN {
                break;
            }
            let i = segs[a].qe + 1 - d;
            gain += k(&tb, i) - k(&ta, i);
            if gain > best.0 {
                best = (gain, -(d as isize));
            }
        }
        if best.1 != 0 {
            let nb = (segs[b].qs as isize + best.1) as usize;
            segs[a].qe = nb - 1;
            segs[b].qs = nb;
            changed = true;
        }
    }
    if !changed {
        return aln;
    }
    let keys = segment_keys_from(ctx, &segs, &aln);
    let trial = chimera_align(ctx, &segs, &keys, work);
    if trial.raw > aln.raw {
        trial
    } else {
        aln
    }
}

/// Junction test between consecutive bodies: the last well-aligned residue of
/// the first body and the first of the next one (sequence separation `sep`)
/// must map to target positions compatible with a connected chain.
fn junction_ok(d: f64, sep: usize) -> bool {
    d <= (6.0 + 1.5 * (sep.max(1) - 1) as f64).min(25.0)
}

/// Raw kernel sum (with `inv_d02`) of the best run of sequence-consecutive
/// bodies whose junctions preserve chain connectivity, and the indices of its
/// bodies in sequence order.
/// Genuine flexibility (hinges, and also circular permutations, whose new
/// termini are spatially adjacent) keeps chains connected, whereas bodies
/// scattered over unrelated structures do not.
pub fn connected_run(aln: &Alignment, q: &[V3], t: &[V3], inv_d02: f64) -> (f64, Vec<usize>) {
    let nb = aln.segs.len();
    let mut raw = vec![0.0; nb];
    let mut first: Vec<Option<(u32, u32)>> = vec![None; nb];
    let mut last: Vec<Option<(u32, u32)>> = vec![None; nb];
    for (k, &(i, j)) in aln.pairs.iter().enumerate() {
        let b = aln.pair_seg[k] as usize;
        let d2 = dist2(&aln.segs[b].tr.apply(&q[i as usize]), &t[j as usize]);
        raw[b] += tm::kernel(d2, inv_d02);
        if d2 <= 25.0 {
            if first[b].is_none_or(|(fi, _)| i < fi) {
                first[b] = Some((i, j));
            }
            if last[b].is_none_or(|(li, _)| i > li) {
                last[b] = Some((i, j));
            }
        }
    }
    let mut by_start: Vec<usize> = (0..nb).collect();
    by_start.sort_by_key(|&b| aln.segs[b].qs);
    let (mut best, mut best_run) = (0.0f64, 0..0);
    let (mut cur, mut cur_start) = (0.0f64, 0usize);
    for (k, &b) in by_start.iter().enumerate() {
        let connected = k > 0 && {
            let a = by_start[k - 1];
            match (last[a], first[b]) {
                (Some((ai, aj)), Some((bi, bj))) if bi > ai => {
                    let d = dist2(&t[aj as usize], &t[bj as usize]).sqrt();
                    junction_ok(d, (bi - ai) as usize)
                }
                _ => false,
            }
        };
        if connected {
            cur += raw[b];
        } else {
            cur = raw[b];
            cur_start = k;
        }
        if cur > best {
            best = cur;
            best_run = cur_start..k + 1;
        }
    }
    (best, by_start[best_run].to_vec())
}

/// Median target position of the well-superposed pairs of each segment
/// (segments without pairs keep their rank in the current order).
fn segment_keys_from(ctx: &Ctx, segs: &[Segment], aln: &Alignment) -> Vec<f64> {
    let q = &ctx.q.s.ca;
    let t = &ctx.t.s.ca;
    let mut keys = vec![0.0; segs.len()];
    for (rank, &si) in aln.order.iter().enumerate() {
        let tr = &segs[si].tr;
        let mut js: Vec<u32> = aln
            .pairs
            .iter()
            .zip(&aln.pair_seg)
            .filter(|(&(i, j), &sg)| {
                sg as usize == si && dist2(&tr.apply(&q[i as usize]), &t[j as usize]) < 25.0
            })
            .map(|(&(_, j), _)| j)
            .collect();
        keys[si] = if js.is_empty() {
            let own: Vec<u32> = aln
                .pairs
                .iter()
                .zip(&aln.pair_seg)
                .filter(|(_, &sg)| sg as usize == si)
                .map(|(&(_, j), _)| j)
                .collect();
            if own.is_empty() {
                rank as f64 * t.len() as f64 / segs.len().max(1) as f64
            } else {
                own[own.len() / 2] as f64
            }
        } else {
            js.sort_unstable();
            js[js.len() / 2] as f64
        };
    }
    keys
}

/// Try to move body `b` to a better place among the target residues not used
/// by the other bodies. Returns the re-aligned solution if it scores higher.
fn refit_body(
    ctx: &Ctx,
    aln: &Alignment,
    b: usize,
    cands_all: &[(usize, f64, Transform)],
    work: &mut DpWork,
) -> Option<Alignment> {
    let q = &ctx.q.s.ca;
    let t = &ctx.t.s.ca;
    let (qs, qe) = (aln.segs[b].qs, aln.segs[b].qe);
    let mut occ = vec![false; t.len()];
    let mut own_raw = 0.0;
    for (k, &(i, j)) in aln.pairs.iter().enumerate() {
        let sg = aln.pair_seg[k] as usize;
        let d2 = dist2(&aln.segs[sg].tr.apply(&q[i as usize]), &t[j as usize]);
        if sg == b {
            own_raw += tm::kernel(d2, ctx.inv_d02);
        } else if d2 < 36.0 {
            occ[j as usize] = true;
        }
    }
    let tree = &ctx.q.tree;
    let overlaps = |node: usize| -> bool {
        let (s, e) = (tree.nodes[node].start, tree.nodes[node].end);
        let lo = s.max(qs);
        let hi = e.min(qe);
        hi >= lo && 2 * (hi - lo + 1) >= (e - s + 1).min(qe - qs + 1)
    };
    // polish all candidate transforms of overlapping PUs on the free target
    let mid = (qs + qe) / 2;
    let mut polished: Vec<(f32, Transform)> = cands_all
        .iter()
        .filter(|(node, _, _)| overlaps(*node))
        .map(|(_, _, tr)| quick_refine(ctx, qs, qe, *tr, Some(&occ)))
        .collect();
    polished.sort_by(|a, b| cmp_desc(a.0 as f64, b.0 as f64));
    let far = [1.0e4, 1.0e4, 1.0e4];
    let tmask: Vec<V3> = t
        .iter()
        .zip(&occ)
        .map(|(p, &o)| if o { far } else { *p })
        .collect();
    let tmask_soa = SoA::new(&tmask);
    let mut sigs: Vec<[V3; 3]> = vec![];
    let mut best: Option<Placement> = None;
    for (_, tr) in polished {
        if sigs.len() >= REFIT_CANDIDATES {
            break;
        }
        let sig = [tr.apply(&q[qs]), tr.apply(&q[mid]), tr.apply(&q[qe])];
        if sigs
            .iter()
            .any(|o| (0..3).all(|k| dist2(&o[k], &sig[k]) < 16.0))
        {
            continue;
        }
        sigs.push(sig);
        if let Some(pl) = refine_placement(ctx, 0, qs, qe, tr, &tmask, &tmask_soa, work) {
            if best.as_ref().is_none_or(|bp| pl.raw > bp.raw) {
                best = Some(pl);
            }
        }
    }
    let best = best?;
    if best.raw <= own_raw + 0.5 {
        return None;
    }
    let mut keys = segment_keys_from(ctx, &aln.segs, aln);
    keys[b] = (best.ta as f64 + best.tb as f64) / 2.0;
    let mut segs = aln.segs.clone();
    segs[b].tr = best.tr;
    let trial = chimera_align(ctx, &segs, &keys, work);
    (trial.raw > aln.raw + 1e-6).then_some(trial)
}

/// Align a peeled query onto a rigid target. Returns (flexible, rigid).
fn align_directional(
    q: &Prepared,
    t: &Prepared,
    lnorm: usize,
    p: &AlignParams,
    work: &mut DpWork,
) -> (Alignment, Alignment, ConnRun) {
    let prof = std::env::var_os("ICARUS_PROFILE").is_some();
    let clock = std::time::Instant::now();
    let lap = |name: &str| {
        if prof {
            eprintln!("  [{:>8.3} ms] {name}", clock.elapsed().as_secs_f64() * 1e3);
        }
    };
    let d0 = tm::d0(lnorm);
    let d0s = tm::d0_search(lnorm);
    let ctx = Ctx {
        q,
        t,
        tsoa: SoA::new(&t.s.ca),
        inv_d02: 1.0 / (d0 * d0),
        inv_d02s: 1.0 / (d0s * d0s),
        grid: DistGrid::build(&t.s.ca, 1.25, 7.5),
        lut: KernelLut::new(d0.clamp(2.0, 5.0)),
    };
    let qca = &q.s.ca;
    let tca = &t.s.ca;
    let m = qca.len();
    let tree = &q.tree;
    let nn = tree.nodes.len();

    // 1. seeds
    let blocks = select_blocks(
        find_blocks(&q.frags, &t.frags, p.q_stride, p.frag_tol),
        m,
        p.max_seeds,
    );
    let mut transforms: Vec<Transform> = blocks
        .iter()
        .map(|b| block_transform(qca, tca, b))
        .collect();
    if transforms.is_empty() {
        // no local similarity at all: fall back to centroid alignment
        let cq = centroid(qca);
        let ct = centroid(tca);
        transforms.push(Transform {
            r: Transform::identity().r,
            t: [ct[0] - cq[0], ct[1] - cq[1], ct[2] - cq[2]],
        });
    }
    lap(&format!("grid + {} seed blocks", transforms.len()));

    // 2. score every transform on every PU (prefix sums of grid kernels).
    // Small PUs fit many places equally well locally, so they get a larger
    // candidate budget (their refinement is also cheaper).
    let limits: Vec<usize> = tree
        .nodes
        .iter()
        .enumerate()
        .map(|(ni, nd)| {
            let f = ((m as f64) / (nd.len() as f64)).sqrt().clamp(1.0, 4.0);
            let per_node = p.per_node.max(1);
            let base = if ni == 0 { per_node + 1 } else { per_node };
            ((base as f64) * f).round() as usize
        })
        .collect();
    let keeps: Vec<usize> = limits.iter().map(|&l| l * 3).collect();
    let mut top: Vec<Vec<(f32, u32)>> = keeps.iter().map(|&k| Vec::with_capacity(k + 1)).collect();
    let stride = if m > 200 { 2 } else { 1 };
    let mut prefix = vec![0f32; m + 1];
    for (ti, tr) in transforms.iter().enumerate() {
        let mut acc = 0f32;
        for i in 0..m {
            if i % stride == 0 {
                let (d, _) = ctx.grid.lookup(&tr.apply(&qca[i]));
                acc += ctx.lut.v[d as usize] * stride as f32;
            }
            prefix[i + 1] = acc;
        }
        for (ni, node) in tree.nodes.iter().enumerate() {
            let sc = prefix[node.end + 1] - prefix[node.start];
            let list = &mut top[ni];
            let keep = keeps[ni];
            if list.len() < keep || sc > list[list.len() - 1].0 {
                let pos = list.partition_point(|&(v, _)| v >= sc);
                list.insert(pos, (sc, ti as u32));
                list.truncate(keep);
            }
        }
    }
    lap("seed scoring");

    // 3. polish and refine candidates of each PU (parents first; the refined
    // placements of a parent are also tried for its children)
    let mut placements: Vec<Placement> = Vec::new();
    let mut node_placements: Vec<Vec<usize>> = vec![Vec::new(); nn];
    let mut refit_pool: Vec<(usize, f64, Transform)> = Vec::new();
    for (ni, node) in tree.nodes.iter().enumerate() {
        let (s, e) = (node.start, node.end);
        let mid = (s + e) / 2;
        let mut cands: Vec<Transform> = Vec::new();
        if let Some(par) = node.parent {
            cands.extend(
                node_placements[par]
                    .iter()
                    .take(2)
                    .map(|&k| placements[k].tr),
            );
        }
        let n_inherited = cands.len();
        cands.extend(top[ni].iter().map(|&(_, ti)| transforms[ti as usize]));
        let polished: Vec<(f32, Transform)> = cands.into_iter().map(|tr| (0.0, tr)).collect();
        let mut sigs: Vec<[V3; 3]> = Vec::new();
        let mut node_best: Vec<Placement> = Vec::new();
        for (sc, tr) in polished {
            if sigs.len() >= limits[ni] + n_inherited {
                break;
            }
            let sig = [tr.apply(&qca[s]), tr.apply(&qca[mid]), tr.apply(&qca[e])];
            if sigs
                .iter()
                .any(|o| (0..3).all(|k| dist2(&o[k], &sig[k]) < SAME_CANDIDATE_D2))
            {
                continue;
            }
            sigs.push(sig);
            refit_pool.push((ni, sc as f64, tr));
            if let Some(pl) = refine_placement(&ctx, ni, s, e, tr, tca, &ctx.tsoa, work) {
                node_best.push(pl);
            }
        }
        // drop duplicates that converged to the same placement
        node_best.sort_by(|a, b| cmp_desc(a.raw, b.raw));
        let mut kept: Vec<Placement> = Vec::new();
        for pl in node_best {
            let dup = kept.iter().any(|k| {
                dist2(&k.tr.apply(&qca[mid]), &pl.tr.apply(&qca[mid])) < 4.0
                    && k.tr.rotation_angle_to(&pl.tr) < 0.2
            });
            if !dup {
                kept.push(pl);
            }
        }
        let best_raw = kept.first().map_or(0.0, |p| p.raw);
        for pl in kept {
            if pl.raw < PRUNE_FRACTION * best_raw {
                continue;
            }
            node_placements[ni].push(placements.len());
            placements.push(pl);
        }
    }
    lap(&format!(
        "refined {} placements over {} PUs",
        placements.len(),
        nn
    ));

    // rigid solution: best whole-chain placement
    let mut rigid = match node_placements[0].first() {
        Some(&k) => {
            let pl = &placements[k];
            finalize(
                &ctx,
                vec![Segment {
                    qs: 0,
                    qe: m - 1,
                    tr: pl.tr,
                }],
                vec![pl.ta as f64],
                work,
            )
        }
        None => Alignment::default(),
    };
    lap("rigid finalize");
    // connectivity-aware score, maximised over every solution examined (with
    // the same hinge penalty per extra body as the flexible solution)
    let lmax = q.len().max(t.len());
    let inv_d02_max = 1.0 / tm::d0(lmax).powi(2);
    let conn_penalty = p.hinge_penalty * lmax as f64;
    let see = |a: &Alignment, conn: &mut ConnRun| {
        let (raw, segs) = connected_run(a, qca, tca, inv_d02_max);
        let score = raw - conn_penalty * segs.len().saturating_sub(1) as f64;
        if segs.is_empty() || score <= conn.score {
            return;
        }
        *conn = ConnRun {
            raw,
            score,
            aln: a.clone(),
            segs,
        };
    };
    let mut conn = ConnRun::default();
    see(&rigid, &mut conn);
    // a one-body solution of the assembly is also a rigid superposition of the
    // whole chain (attach_gaps extends the body), found from a PU placement
    // that the whole-chain search may have missed
    let one_body = |a: &Alignment, rigid: &mut Alignment| {
        if a.segs.len() == 1 && a.segs[0].qs == 0 && a.segs[0].qe == m - 1 && a.raw > rigid.raw {
            *rigid = a.clone();
        }
    };
    if tree.nodes.len() <= 1 || placements.is_empty() {
        return (rigid.clone(), rigid, conn);
    }

    // 4. assemble rigid bodies
    let masks: Vec<u32> = tree.nodes.iter().map(|n| n.leaf_mask).collect();
    let penalty = p.hinge_penalty * lnorm as f64;
    let solutions = assemble(
        &placements,
        &masks,
        tree.leaves.len(),
        p.max_segments,
        penalty,
        &mut work.asm,
    );
    lap(&format!(
        "assembly DP ({} placements, {} leaves)",
        placements.len(),
        tree.leaves.len()
    ));
    if solutions.is_empty() {
        return (rigid.clone(), rigid, conn);
    }
    // finalise the best selection of every body count, keep the best two
    let mut finals: Vec<Alignment> = Vec::new();
    for (_, chosen) in &solutions {
        let mut segs: Vec<(Segment, f64)> = chosen
            .iter()
            .map(|&c| {
                let pl = &placements[c];
                let nd = &tree.nodes[pl.node];
                (
                    Segment {
                        qs: nd.start,
                        qe: nd.end,
                        tr: pl.tr,
                    },
                    (pl.ta as f64 + pl.tb as f64) / 2.0,
                )
            })
            .collect();
        segs.sort_by_key(|(s, _)| s.qs);
        // attach unplaced residues (PUs without a placement) to a neighbouring body
        attach_gaps(&ctx, &mut segs, m);
        let keys: Vec<f64> = segs.iter().map(|(_, k)| *k).collect();
        let segs: Vec<Segment> = segs.into_iter().map(|(s, _)| s).collect();
        // one chimera pass to rank the candidate solutions
        let aln = chimera_align(&ctx, &segs, &keys, work);
        see(&aln, &mut conn);
        one_body(&aln, &mut rigid);
        finals.push(aln);
    }
    let obj = |a: &Alignment| a.raw - penalty * (a.segs.len().max(1) - 1) as f64;
    finals.sort_by(|a, b| cmp_desc(obj(a), obj(b)));
    finals.truncate(FINALIZED_SOLUTIONS);
    let mut finals: Vec<Alignment> = finals
        .into_iter()
        .map(|a| {
            let keys = segment_keys_from(&ctx, &a.segs, &a);
            let f = finalize(&ctx, a.segs.clone(), keys, work);
            see(&f, &mut conn);
            one_body(&f, &mut rigid);
            if f.raw >= a.raw {
                f
            } else {
                a
            }
        })
        .collect();
    finals.sort_by(|a, b| cmp_desc(obj(a), obj(b)));
    finals.truncate(REFITTED_SOLUTIONS);
    lap(&format!(
        "assembly + finalize ({} candidate solutions)",
        solutions.len()
    ));
    let mut best_flex: Option<Alignment> = None;
    for mut flex in finals {
        // 5. refit bodies against the free target
        if p.refit && flex.segs.len() > 1 {
            let mut changed = false;
            for _round in 0..REFIT_ROUNDS {
                let mut improved = false;
                for b in 0..flex.segs.len() {
                    if let Some(trial) = refit_body(&ctx, &flex, b, &refit_pool, work) {
                        flex = trial;
                        improved = true;
                        changed = true;
                    }
                }
                if !improved {
                    break;
                }
            }
            if changed {
                let keys = segment_keys_from(&ctx, &flex.segs, &flex);
                let refined = finalize(&ctx, flex.segs.clone(), keys, work);
                if refined.raw > flex.raw {
                    flex = refined;
                }
            }
            lap(&format!("refit raw={:.2}", flex.raw));
        }
        {
            for _ in 0..2 {
                let before = flex.raw;
                flex = refine_hinges(&ctx, flex, work);
                if flex.raw <= before + 1e-6 {
                    break;
                }
            }
            {
                let keys = segment_keys_from(&ctx, &flex.segs, &flex);
                let f2 = finalize(&ctx, flex.segs.clone(), keys, work);
                if f2.raw > flex.raw {
                    flex = f2;
                }
            }
            lap(&format!("hinges raw={:.2}", flex.raw));
        }
        if best_flex.as_ref().is_none_or(|b| obj(&flex) > obj(b)) {
            best_flex = Some(flex);
        }
    }
    let flex = best_flex.unwrap();
    see(&flex, &mut conn);
    one_body(&flex, &mut rigid);
    if obj(&flex) >= rigid.raw {
        (flex, rigid, conn)
    } else {
        (rigid.clone(), rigid, conn)
    }
}

fn centroid(x: &[V3]) -> V3 {
    let mut c = [0.0; 3];
    for p in x {
        for k in 0..3 {
            c[k] += p[k];
        }
    }
    let n = x.len().max(1) as f64;
    [c[0] / n, c[1] / n, c[2] / n]
}

/// Extend bodies over query residues not covered by any chosen PU, choosing
/// for each gap the split between neighbouring bodies that superposes best.
fn attach_gaps(ctx: &Ctx, segs: &mut [(Segment, f64)], m: usize) {
    let q = &ctx.q.s.ca;
    let score = |tr: &Transform, a: usize, b: usize| -> f32 {
        (a..=b)
            .map(|i| ctx.lut.v[ctx.grid.lookup(&tr.apply(&q[i])).0 as usize])
            .sum()
    };
    let k = segs.len();
    if k == 0 {
        return;
    }
    if segs[0].0.qs > 0 {
        segs[0].0.qs = 0;
    }
    for i in 0..k - 1 {
        let gap_s = segs[i].0.qe + 1;
        let gap_e = segs[i + 1].0.qs;
        if gap_e <= gap_s {
            continue;
        }
        let gap_e = gap_e - 1;
        let mut best = (f32::MIN, gap_s);
        for cut in gap_s..=gap_e + 1 {
            let left = if cut > gap_s {
                score(&segs[i].0.tr, gap_s, cut - 1)
            } else {
                0.0
            };
            let right = if cut <= gap_e {
                score(&segs[i + 1].0.tr, cut, gap_e)
            } else {
                0.0
            };
            if left + right > best.0 {
                best = (left + right, cut);
            }
        }
        segs[i].0.qe = best.1 - 1;
        segs[i + 1].0.qs = best.1;
    }
    if segs[k - 1].0.qe < m - 1 {
        segs[k - 1].0.qe = m - 1;
    }
}

/// Align two prepared structures (ICARUS flexible alignment).
pub fn align_pair(a: &Prepared, b: &Prepared, p: &AlignParams, work: &mut DpWork) -> PairResult {
    let lnorm = a.len().min(b.len());
    let (f1, r1, c1) = align_directional(a, b, lnorm, p, work);
    let mut res = PairResult {
        flex: f1,
        rigid: r1,
        reversed: false,
        rigid_reversed: false,
        len1: a.len(),
        len2: b.len(),
        lnorm,
        conn: c1,
        conn_reversed: false,
    };
    if p.both_directions {
        let (f2, r2, c2) = align_directional(b, a, lnorm, p, work);
        if c2.score > res.conn.score {
            res.conn = c2;
            res.conn_reversed = true;
        }
        if f2.raw > res.flex.raw {
            res.flex = f2;
            res.reversed = true;
        }
        if r2.raw > res.rigid.raw {
            res.rigid = r2;
            res.rigid_reversed = true;
        }
    }
    res
}
