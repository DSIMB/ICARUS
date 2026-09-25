//! Dynamic programming alignment of two coordinate sets in a common frame,
//! scored with the TM kernel 1/(1 + d²/d0²).

use crate::geom::V3;

const DIAG: u8 = 1;
const UP: u8 = 2; // consume x (x_i against a gap)
const LEFT: u8 = 3; // consume y

/// Target coordinates in structure-of-arrays f32 layout (SIMD-friendly).
#[derive(Default, Clone)]
pub struct SoA {
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Vec<f32>,
}

impl SoA {
    pub fn new(p: &[V3]) -> Self {
        Self {
            x: p.iter().map(|v| v[0] as f32).collect(),
            y: p.iter().map(|v| v[1] as f32).collect(),
            z: p.iter().map(|v| v[2] as f32).collect(),
        }
    }
    pub fn len(&self) -> usize {
        self.x.len()
    }
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}

#[inline]
fn kernel_row(xi: &V3, y: &SoA, lo: usize, hi: usize, inv_d02: f32, row: &mut [f32]) {
    let (x0, x1, x2) = (xi[0] as f32, xi[1] as f32, xi[2] as f32);
    let (ys0, ys1, ys2) = (&y.x[lo..hi], &y.y[lo..hi], &y.z[lo..hi]);
    for (k, r) in row.iter_mut().enumerate().take(hi - lo) {
        let dx = x0 - ys0[k];
        let dy = x1 - ys1[k];
        let dz = x2 - ys2[k];
        *r = 1.0 / (1.0 + (dx * dx + dy * dy + dz * dz) * inv_d02);
    }
}

/// Reusable DP buffers.
#[derive(Default)]
pub struct DpWork {
    /// Scratch table reused by the rigid-body assembly DP.
    pub asm: Vec<f32>,
    rowf: Vec<f32>,
    prevf: Vec<f32>,
    curf: Vec<f32>,
    prev: Vec<f64>,
    cur: Vec<f64>,
    diag_prev: Vec<bool>,
    diag_cur: Vec<bool>,
    trace: Vec<u8>,
}

impl DpWork {
    fn prepare_f(&mut self, m: usize, n: usize) {
        self.rowf.clear();
        self.rowf.resize(n, 0.0);
        self.prevf.clear();
        self.prevf.resize(n + 1, 0.0);
        self.curf.clear();
        self.curf.resize(n + 1, 0.0);
        self.diag_prev.clear();
        self.diag_prev.resize(n + 1, false);
        self.diag_cur.clear();
        self.diag_cur.resize(n + 1, false);
        self.trace.clear();
        self.trace.resize((m + 1) * (n + 1), 0);
    }

    /// f32 version of `align_gap_open` against target residues lo..hi of `y`.
    /// Returned pairs index the full target.
    pub fn align_gap_open_f(
        &mut self,
        x: &[V3],
        y: &SoA,
        lo: usize,
        hi: usize,
        inv_d02: f32,
        gap_open: f32,
        out: &mut Vec<(u32, u32)>,
    ) -> f32 {
        let (m, n) = (x.len(), hi - lo);
        self.prepare_f(m, n);
        let w = n + 1;
        for i in 1..=m {
            kernel_row(&x[i - 1], y, lo, hi, inv_d02, &mut self.rowf);
            self.curf[0] = 0.0;
            self.diag_cur[0] = false;
            let tr = &mut self.trace[i * w..(i + 1) * w];
            for j in 1..=n {
                let d = self.prevf[j - 1] + self.rowf[j - 1];
                let mut h = self.prevf[j];
                if self.diag_prev[j] && j < n {
                    h += gap_open;
                }
                let mut v = self.curf[j - 1];
                if self.diag_cur[j - 1] && i < m {
                    v += gap_open;
                }
                tr[j] = if d >= h && d >= v {
                    self.diag_cur[j] = true;
                    self.curf[j] = d;
                    DIAG
                } else if h >= v {
                    self.diag_cur[j] = false;
                    self.curf[j] = h;
                    UP
                } else {
                    self.diag_cur[j] = false;
                    self.curf[j] = v;
                    LEFT
                };
            }
            std::mem::swap(&mut self.prevf, &mut self.curf);
            std::mem::swap(&mut self.diag_prev, &mut self.diag_cur);
        }
        let score = self.prevf[n];
        self.traceback(m, n, out);
        for p in out.iter_mut() {
            p.1 += lo as u32;
        }
        score
    }

    /// f32 version of `align_free` (gdt2 DP) against the full target.
    pub fn align_free_f(
        &mut self,
        x: &[V3],
        y: &SoA,
        inv_d02: f32,
        out: &mut Vec<(u32, u32)>,
    ) -> f32 {
        let (m, n) = (x.len(), y.len());
        self.prepare_f(m, n);
        let w = n + 1;
        for i in 1..=m {
            kernel_row(&x[i - 1], y, 0, n, inv_d02, &mut self.rowf);
            self.curf[0] = 0.0;
            let tr = &mut self.trace[i * w..(i + 1) * w];
            for j in 1..=n {
                let d = self.prevf[j - 1] + self.rowf[j - 1];
                let up = self.prevf[j];
                let left = self.curf[j - 1];
                tr[j] = if d >= up && d >= left {
                    self.curf[j] = d;
                    DIAG
                } else if up <= left {
                    self.curf[j] = left;
                    LEFT
                } else {
                    self.curf[j] = up;
                    UP
                };
            }
            std::mem::swap(&mut self.prevf, &mut self.curf);
        }
        let score = self.prevf[n];
        self.traceback(m, n, out);
        score
    }

    /// Banded variant of `align_gap_open_f`: row i may only use target
    /// columns [lo_i, hi_i) with lo_i = clamp(center_i - band), both bounds
    /// non-decreasing. Local semantics (free end gaps anywhere). Pairs index
    /// the full target.
    #[allow(clippy::too_many_arguments)]
    pub fn align_band_f(
        &mut self,
        x: &[V3],
        y: &SoA,
        centers: &[i32],
        band: usize,
        inv_d02: f32,
        gap_open: f32,
        out: &mut Vec<(u32, u32)>,
    ) -> f32 {
        let m = x.len();
        let n = y.len() as i32;
        let w = 2 * band + 1;
        // monotone band bounds
        let mut lo = vec![0i32; m];
        let mut hi = vec![0i32; m];
        let mut run_lo = 0i32;
        for i in 0..m {
            let c = centers[i];
            run_lo = run_lo.max((c - band as i32).clamp(0, n - 1));
            lo[i] = run_lo;
        }
        let mut run_hi = n;
        for i in (0..m).rev() {
            let c = centers[i];
            run_hi = run_hi.min((c + band as i32 + 1).clamp(1, n));
            hi[i] = run_hi.max(lo[i] + 1);
        }
        for i in 1..m {
            hi[i] = hi[i].max(hi[i - 1]);
        }
        let neg = f32::NEG_INFINITY;
        self.prevf.clear();
        self.prevf.resize(n as usize + 1, neg);
        self.curf.clear();
        self.curf.resize(n as usize + 1, neg);
        self.diag_prev.clear();
        self.diag_prev.resize(n as usize + 1, false);
        self.diag_cur.clear();
        self.diag_cur.resize(n as usize + 1, false);
        let wmax = hi
            .iter()
            .zip(&lo)
            .map(|(h, l)| (h - l) as usize)
            .max()
            .unwrap_or(0)
            .max(w);
        self.trace.clear();
        self.trace.resize(m * wmax, 0);
        self.rowf.clear();
        self.rowf.resize(wmax, 0.0);
        let (mut best, mut bi, mut bj) = (0f32, usize::MAX, 0usize);
        let (mut plo, mut phi) = (0i32, 0i32);
        for i in 0..m {
            let (l, h) = (lo[i], hi[i]);
            kernel_row(&x[i], y, l as usize, h as usize, inv_d02, &mut self.rowf);
            let tr = &mut self.trace[i * wmax..(i + 1) * wmax];
            for j in l..h {
                let ju = j as usize;
                let s = self.rowf[(j - l) as usize];
                // predecessors (row i-1 valid in [plo, phi))
                let dprev = if i > 0 && j > plo && j - 1 < phi {
                    self.prevf[ju - 1]
                } else {
                    neg
                };
                let d = dprev.max(0.0) + s;
                let mut hup = if i > 0 && j >= plo && j < phi {
                    self.prevf[ju]
                } else {
                    neg
                };
                if hup != neg && self.diag_prev[ju] {
                    hup += gap_open;
                }
                let mut vl = if j > l { self.curf[ju - 1] } else { neg };
                if vl != neg && self.diag_cur[ju - 1] {
                    vl += gap_open;
                }
                let t = if d >= hup && d >= vl {
                    self.diag_cur[ju] = true;
                    self.curf[ju] = d;
                    if dprev > 0.0 {
                        DIAG
                    } else {
                        4
                    } // 4 = alignment start
                } else if hup >= vl {
                    self.diag_cur[ju] = false;
                    self.curf[ju] = hup;
                    UP
                } else {
                    self.diag_cur[ju] = false;
                    self.curf[ju] = vl;
                    LEFT
                };
                tr[(j - l) as usize] = t;
                if self.curf[ju] > best {
                    best = self.curf[ju];
                    bi = i;
                    bj = ju;
                }
            }
            // clear previous row outside the new band to keep invariants
            for j in plo..phi {
                self.prevf[j as usize] = neg;
            }
            std::mem::swap(&mut self.prevf, &mut self.curf);
            std::mem::swap(&mut self.diag_prev, &mut self.diag_cur);
            plo = l;
            phi = h;
        }
        out.clear();
        if bi == usize::MAX {
            return 0.0;
        }
        let (mut i, mut j) = (bi as i64, bj as i64);
        while i >= 0 && j >= lo[i as usize] as i64 && j < hi[i as usize] as i64 {
            let t = self.trace[i as usize * wmax + (j - lo[i as usize] as i64) as usize];
            match t {
                DIAG => {
                    out.push((i as u32, j as u32));
                    i -= 1;
                    j -= 1;
                }
                4 => {
                    out.push((i as u32, j as u32));
                    break;
                }
                UP => i -= 1,
                _ => j -= 1,
            }
        }
        out.reverse();
        best
    }

    fn prepare(&mut self, m: usize, n: usize) {
        self.prev.clear();
        self.prev.resize(n + 1, 0.0);
        self.cur.clear();
        self.cur.resize(n + 1, 0.0);
        self.diag_prev.clear();
        self.diag_prev.resize(n + 1, false);
        self.diag_cur.clear();
        self.diag_cur.resize(n + 1, false);
        self.trace.clear();
        self.trace.resize((m + 1) * (n + 1), 0);
    }

    fn traceback(&self, m: usize, n: usize, out: &mut Vec<(u32, u32)>) {
        out.clear();
        let w = n + 1;
        let (mut i, mut j) = (m, n);
        while i > 0 && j > 0 {
            match self.trace[i * w + j] {
                DIAG => {
                    out.push(((i - 1) as u32, (j - 1) as u32));
                    i -= 1;
                    j -= 1;
                }
                UP => i -= 1,
                _ => j -= 1,
            }
        }
        out.reverse();
    }

    /// TM-align style DP: a gap costs `gap_open` when it opens after a match,
    /// extensions are free, and end gaps are free. Returns aligned (i, j) pairs.
    pub fn align_gap_open(
        &mut self,
        x: &[V3],
        y: &[V3],
        inv_d02: f64,
        gap_open: f64,
        out: &mut Vec<(u32, u32)>,
    ) -> f64 {
        let (m, n) = (x.len(), y.len());
        self.prepare(m, n);
        let w = n + 1;
        for i in 1..=m {
            let xi = x[i - 1];
            self.cur[0] = 0.0;
            self.diag_cur[0] = false;
            for j in 1..=n {
                let yj = &y[j - 1];
                let dx = xi[0] - yj[0];
                let dy = xi[1] - yj[1];
                let dz = xi[2] - yj[2];
                let s = 1.0 / (1.0 + (dx * dx + dy * dy + dz * dz) * inv_d02);
                let d = self.prev[j - 1] + s;
                let mut h = self.prev[j];
                if self.diag_prev[j] && j < n {
                    h += gap_open;
                }
                let mut v = self.cur[j - 1];
                if self.diag_cur[j - 1] && i < m {
                    v += gap_open;
                }
                let t = if d >= h && d >= v {
                    self.diag_cur[j] = true;
                    self.cur[j] = d;
                    DIAG
                } else if h >= v {
                    self.diag_cur[j] = false;
                    self.cur[j] = h;
                    UP
                } else {
                    self.diag_cur[j] = false;
                    self.cur[j] = v;
                    LEFT
                };
                self.trace[i * w + j] = t;
            }
            std::mem::swap(&mut self.prev, &mut self.cur);
            std::mem::swap(&mut self.diag_prev, &mut self.diag_cur);
        }
        let score = self.prev[n];
        self.traceback(m, n, out);
        score
    }

    /// Gap-penalty-free DP maximising the kernel sum (the alignment step of
    /// gdt2.pl, used by ICARUS to score flexible alignments). Tie-breaking
    /// follows gdt2.pl: diagonal first, then left if up <= left.
    pub fn align_free(
        &mut self,
        x: &[V3],
        y: &[V3],
        inv_d02: f64,
        out: &mut Vec<(u32, u32)>,
    ) -> f64 {
        let (m, n) = (x.len(), y.len());
        self.prepare(m, n);
        let w = n + 1;
        for i in 1..=m {
            let xi = x[i - 1];
            self.cur[0] = 0.0;
            for j in 1..=n {
                let yj = &y[j - 1];
                let dx = xi[0] - yj[0];
                let dy = xi[1] - yj[1];
                let dz = xi[2] - yj[2];
                let s = 1.0 / (1.0 + (dx * dx + dy * dy + dz * dz) * inv_d02);
                let d = self.prev[j - 1] + s;
                let up = self.prev[j];
                let left = self.cur[j - 1];
                let t = if d >= up && d >= left {
                    self.cur[j] = d;
                    DIAG
                } else if up <= left {
                    self.cur[j] = left;
                    LEFT
                } else {
                    self.cur[j] = up;
                    UP
                };
                self.trace[i * w + j] = t;
            }
            std::mem::swap(&mut self.prev, &mut self.cur);
        }
        let score = self.prev[n];
        self.traceback(m, n, out);
        score
    }
}
