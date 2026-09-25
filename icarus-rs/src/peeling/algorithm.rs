//! Protein Peeling (Gelly et al., 2006): hierarchical cutting of a protein into
//! compact Protein Units (PUs).
//!
//! Vendored from SWORD3 (DSIMB/sword3 @ e47f837,
//! `sword3-lib/src/peeling/algorithm.rs`, CeCILL-2.1), itself a port of
//! `peeling_omp.c`. Numerics of the cutting criterion, cutting mask and stop
//! criteria are unchanged. ICARUS-specific changes: the double cut search is
//! serial (ties resolve like the original serial C code), log/matrix file
//! writers and per-iteration contact-ratio metrics are removed, and the output
//! is the list of PU delineations per iteration.
//!
//! 1. Contact probability matrix from C-alpha coordinates
//! 2. Cutting mask: no cut inside secondary structure segments shorter than
//!    `min_ss_size`
//! 3. Each iteration applies the best single or double cut (over all current
//!    PUs) maximising the Matthews-like coefficient (ab - c²)/((a+c)(b+c))
//! 4. Stop when the Compaction Index exceeds `max_r2`, no cut is possible, or
//!    the maximum number of PUs is reached

use super::contact_matrix::ContactMatrix;
use super::interval::PuSpan;

const MAX_ITERATION: usize = 64;

/// Configuration parameters for the peeling algorithm.
/// Defaults match the ICARUS v1 invocation of the Peeling binary:
/// `-R2 98 -ss2 8 -lspu 15 -mspu 0 -d0 6.0 -delta 1.5 -oss 0 -p 0 -cp 0 -npu 30`.
#[derive(Debug, Clone)]
pub struct PeelingConfig {
    pub max_r2: i32,
    pub min_ss_size: usize,
    pub min_pu_size: usize,
    pub max_pu_size: usize,
    pub d0: f64,
    pub delta: f64,
    pub pruning: bool,
    pub cutoff_pruning: f64,
    pub max_pu_number: usize,
}

impl Default for PeelingConfig {
    fn default() -> Self {
        Self {
            max_r2: 98,
            min_ss_size: 8,
            min_pu_size: 15,
            max_pu_size: 0,
            d0: 6.0,
            delta: 1.5,
            pruning: false,
            cutoff_pruning: 0.0,
            max_pu_number: 30,
        }
    }
}

/// Secondary structure classification for a residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SsType {
    Coil = 0,
    Helix = 1,
    Sheet = 2,
}

impl SsType {
    pub fn from_dssp_char(value: char) -> Self {
        match value {
            'H' | 'G' => Self::Helix,
            'E' | 'B' => Self::Sheet,
            _ => Self::Coil,
        }
    }
}

/// `mask[i + 1]` says whether a cut is allowed immediately after residue `i`:
/// boundaries strictly inside secondary-structure segments shorter than
/// `min_ss_size` are not cuttable.
pub fn build_cutting_mask(ss_types: &[SsType], min_ss_size: usize) -> Vec<bool> {
    let mut cutting_mask = vec![true; ss_types.len() + 1];
    let mut small: Vec<(usize, usize)> = Vec::new();
    let mut current_ss = SsType::Coil;
    let mut segment_start = 0usize;
    let close = |end: usize, start: usize, out: &mut Vec<(usize, usize)>| {
        if let Ok(span) = PuSpan::new(start, end) {
            if span.len() < min_ss_size {
                out.push((span.start(), span.end()));
            }
        }
    };
    for i in 0..ss_types.len() {
        if ss_types[i] != SsType::Coil && current_ss == SsType::Coil {
            current_ss = ss_types[i];
            segment_start = i;
        } else if ss_types[i] != current_ss && ss_types[i] != SsType::Coil && current_ss != SsType::Coil {
            close(i - 1, segment_start, &mut small);
            current_ss = ss_types[i];
            segment_start = i;
        } else if ss_types[i] == SsType::Coil && current_ss != SsType::Coil {
            close(i - 1, segment_start, &mut small);
            current_ss = SsType::Coil;
        }
    }
    if current_ss != SsType::Coil && !ss_types.is_empty() {
        close(ss_types.len() - 1, segment_start, &mut small);
    }
    for (s, e) in small {
        for k in s..e {
            cutting_mask[k + 1] = false;
        }
    }
    cutting_mask
}

#[derive(Debug, Clone, Copy)]
struct CutResult {
    coeff: f64,
    num_cuts: usize,
    pu_index: usize,
    start: usize,
    i1: usize,
    i2: usize,
    j1: usize,
    j2: usize,
    end: usize,
}

fn simple_cutting(
    m: &ContactMatrix,
    pu_start: usize,
    pu_end: usize,
    mask: &[bool],
    min_pu_size: usize,
    current_best: f64,
    pu_index: usize,
) -> Option<CutResult> {
    let mut best: Option<CutResult> = None;
    let mut best_coeff = current_best;
    for i in pu_start..pu_end {
        if i + 1 < mask.len() && !mask[i + 1] {
            continue;
        }
        if i - pu_start + 1 < min_pu_size || pu_end - i < min_pu_size {
            continue;
        }
        let a = m.rectangle_sum(pu_start, pu_start, i, i);
        let b = m.rectangle_sum(i + 1, i + 1, pu_end, pu_end);
        let c = m.rectangle_sum(pu_start, i + 1, i, pu_end);
        let denom = (a + c) * (b + c);
        if denom == 0.0 {
            continue;
        }
        let coeff = (a * b - c * c) / denom;
        if coeff > best_coeff {
            best_coeff = coeff;
            best = Some(CutResult { coeff, num_cuts: 1, pu_index, start: pu_start, i1: i, i2: i + 1, j1: 0, j2: 0, end: pu_end });
        }
    }
    best
}

fn double_cutting(
    m: &ContactMatrix,
    pu_start: usize,
    pu_end: usize,
    mask: &[bool],
    min_pu_size: usize,
    current_best: f64,
    pu_index: usize,
) -> Option<CutResult> {
    let min_seg = min_pu_size;
    let coo_max_i = pu_end.saturating_sub(min_seg);
    let coo_max_j = pu_end.saturating_sub(min_seg / 2);
    let coo_min_i = pu_start + min_seg - 1;
    if coo_min_i >= coo_max_i {
        return None;
    }
    let mut best_coeff = current_best;
    let mut best: Option<CutResult> = None;
    for i in coo_min_i..coo_max_i {
        if i + 1 < mask.len() && !mask[i + 1] {
            continue;
        }
        let coo_min_j = i + min_seg;
        if coo_min_j + min_seg >= coo_max_j {
            continue;
        }
        let i2 = i + 1;
        // Terms that depend on i only
        let b1 = m.rectangle_sum(pu_start, pu_start, i, i);
        for j in coo_min_j..=coo_max_j {
            if j + 1 < mask.len() && !mask[j + 1] {
                continue;
            }
            let j2 = j + 1;
            if j2 > pu_end {
                continue;
            }
            if i - pu_start + 1 < min_seg || j + 1 - i2 < min_seg {
                continue;
            }
            let a = m.rectangle_sum(i2, i2, j, j);
            let b2 = m.rectangle_sum(j2, j2, pu_end, pu_end);
            let b3 = m.rectangle_sum(pu_start, j2, i, pu_end);
            let b = b1 + b2 + 2.0 * b3;
            let c1 = m.rectangle_sum(i2, pu_start, j, i);
            let c2 = m.rectangle_sum(j2, i2, pu_end, j);
            let c = c1 + c2;
            let denom = (a + c) * (b + c);
            if denom == 0.0 {
                continue;
            }
            let coeff = (a * b - c * c) / denom;
            if coeff > best_coeff {
                best_coeff = coeff;
                best = Some(CutResult { coeff, num_cuts: 2, pu_index, start: pu_start, i1: i, i2, j1: j, j2, end: pu_end });
            }
        }
    }
    best
}

/// Compaction Index (CI) of a PU set, via the mutual information of the
/// PU-PU contact distribution.
fn compaction_index(m: &ContactMatrix, pus: &[[usize; 2]]) -> f64 {
    let n_pus = pus.len();
    let mut zone = vec![0.0f64; n_pus * n_pus];
    for x in 0..n_pus {
        for y in 0..n_pus {
            zone[x * n_pus + y] = m.rectangle_sum(pus[x][0], pus[y][0], pus[x][1], pus[y][1]);
        }
    }
    let mut marg = vec![0.0f64; n_pus];
    let mut tot = 0.0f64;
    for x in 0..n_pus {
        for y in 0..n_pus {
            marg[x] += zone[x * n_pus + y];
        }
        tot += marg[x];
    }
    if tot > 0.0 {
        zone.iter_mut().for_each(|v| *v /= tot);
        marg.iter_mut().for_each(|v| *v /= tot);
    }
    let mut entropy = 0.0f64;
    for x in 0..n_pus {
        for y in 0..n_pus {
            let pxy = zone[x * n_pus + y];
            if pxy > 1e-5 && marg[x] > 1e-5 && marg[y] > 1e-5 {
                entropy += pxy * (pxy / (marg[x] * marg[y])).ln();
            }
        }
    }
    100.0 * (1.0 - (-2.0 * entropy).exp()).sqrt()
}

fn homogeneity(m: &ContactMatrix, start: usize, end: usize) -> f64 {
    let threshold = 0.5;
    let (mut pc1, mut pc2) = (0.0f64, 0.0f64);
    for k in start..=end {
        for l in start..=end {
            let p = m.get(k, l);
            if p > threshold {
                pc1 += p;
                if (k as isize - l as isize).unsigned_abs() < 6 {
                    pc2 += p;
                }
            }
        }
    }
    let (mut h1, mut h2) = (0.0f64, 0.0f64);
    for i in start..=end {
        for j in start..=end {
            let p = m.get(i, j);
            if p < 0.0001 || p <= threshold {
                continue;
            }
            if pc1 > 0.0 {
                let pn1 = p / pc1;
                h1 += pn1 * pn1.ln();
            }
            if (i as isize - j as isize).unsigned_abs() < 6 && pc2 > 0.0 {
                let pn2 = p / pc2;
                h2 += pn2 * pn2.ln();
            }
        }
    }
    ((-h1).exp() - (-h2).exp()) / (end - start + 1) as f64
}

/// One peeling iteration: PU spans (0-based, inclusive) and its CI.
#[derive(Debug, Clone)]
pub struct Iteration {
    pub ci: f64,
    pub pus: Vec<[usize; 2]>,
}

/// Run Protein Peeling on C-alpha coordinates and per-residue secondary structure.
/// Returns the successive iterations (each one refines the previous one by
/// splitting exactly one PU into two or three).
pub fn run_peeling(ca: &[[f64; 3]], ss: &[SsType], config: &PeelingConfig) -> Vec<Iteration> {
    let n = ca.len();
    let mut iterations = Vec::new();
    if n == 0 || config.min_pu_size == 0 {
        return iterations;
    }
    let matrix = ContactMatrix::from_ca_coords(ca, config.d0, config.delta);
    let mask = build_cutting_mask(ss, config.min_ss_size);
    let mut current: Vec<[usize; 2]> = vec![[0, n - 1]];
    for _ in 1..MAX_ITERATION {
        if current.len() > config.max_pu_number {
            break;
        }
        let mut best_cut: Option<CutResult> = None;
        let mut best_coeff = 0.0f64;
        for (x, &[start, end]) in current.iter().enumerate() {
            if end - start + 1 < config.min_pu_size {
                continue;
            }
            if let Some(cut) = simple_cutting(&matrix, start, end, &mask, config.min_pu_size, best_coeff, x) {
                best_coeff = cut.coeff;
                best_cut = Some(cut);
            }
            if let Some(cut) = double_cutting(&matrix, start, end, &mask, config.min_pu_size, best_coeff, x) {
                best_coeff = cut.coeff;
                best_cut = Some(cut);
            }
        }
        let Some(cut) = best_cut else { break };
        let mut new_pus: Vec<[usize; 2]> = Vec::with_capacity(current.len() + 2);
        if cut.num_cuts == 1 {
            new_pus.push([cut.start, cut.i1]);
            new_pus.push([cut.i2, cut.end]);
        } else {
            new_pus.push([cut.start, cut.i1]);
            new_pus.push([cut.i2, cut.j1]);
            new_pus.push([cut.j2, cut.end]);
        }
        for (x, pu) in current.iter().enumerate() {
            if x != cut.pu_index {
                new_pus.push(*pu);
            }
        }
        if new_pus.len() > config.max_pu_number {
            break;
        }
        let reached_max_size = config.max_pu_size > 0
            && new_pus.iter().all(|&[s, e]| e - s + 1 <= config.max_pu_size);
        if config.pruning {
            let h = |s, e| homogeneity(&matrix, s, e) >= config.cutoff_pruning;
            let ok = if cut.num_cuts == 1 {
                h(cut.start, cut.i1) || h(cut.i2, cut.end)
            } else {
                h(cut.start, cut.i1) || h(cut.i2, cut.j1) || h(cut.j2, cut.end)
            };
            if !ok {
                break;
            }
        }
        let ci = compaction_index(&matrix, &new_pus);
        iterations.push(Iteration { ci, pus: new_pus.clone() });
        current = new_pus;
        if reached_max_size || ci > config.max_r2 as f64 {
            break;
        }
    }
    iterations
}
