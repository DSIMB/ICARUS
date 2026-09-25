//! Distance-transform grid of a target structure.
//!
//! Every grid cell stores the distance from its centre to the nearest target
//! C-alpha (up to `rmax`) and the index of that residue. Scoring a candidate
//! superposition then costs one memory access per query residue, which makes
//! it possible to rank thousands of candidate transforms per pair.

use crate::geom::{dist2, V3};

pub struct DistGrid {
    origin: V3,
    inv_h: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    /// Nearest distance in units of 0.1 Å (255 = beyond rmax).
    dist: Vec<u8>,
    nearest: Vec<u32>,
}

pub const FAR: u8 = 255;

impl DistGrid {
    pub fn build(ca: &[V3], h: f64, rmax: f64) -> Self {
        let mut lo = [f64::MAX; 3];
        let mut hi = [f64::MIN; 3];
        for p in ca {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let pad = rmax + h;
        let origin = [lo[0] - pad, lo[1] - pad, lo[2] - pad];
        let nx = ((hi[0] - lo[0] + 2.0 * pad) / h).ceil() as usize + 1;
        let ny = ((hi[1] - lo[1] + 2.0 * pad) / h).ceil() as usize + 1;
        let nz = ((hi[2] - lo[2] + 2.0 * pad) / h).ceil() as usize + 1;
        let total = nx * ny * nz;
        let mut d2best = vec![f32::MAX; total];
        let mut nearest = vec![u32::MAX; total];
        let r = (rmax / h).ceil() as isize;
        let rmax2 = rmax * rmax;
        for (idx, p) in ca.iter().enumerate() {
            let cx = ((p[0] - origin[0]) / h).round() as isize;
            let cy = ((p[1] - origin[1]) / h).round() as isize;
            let cz = ((p[2] - origin[2]) / h).round() as isize;
            for gx in (cx - r).max(0)..=(cx + r).min(nx as isize - 1) {
                let x = origin[0] + gx as f64 * h;
                let dx2 = (x - p[0]) * (x - p[0]);
                if dx2 > rmax2 {
                    continue;
                }
                for gy in (cy - r).max(0)..=(cy + r).min(ny as isize - 1) {
                    let y = origin[1] + gy as f64 * h;
                    let dxy2 = dx2 + (y - p[1]) * (y - p[1]);
                    if dxy2 > rmax2 {
                        continue;
                    }
                    let base = (gx as usize * ny + gy as usize) * nz;
                    for gz in (cz - r).max(0)..=(cz + r).min(nz as isize - 1) {
                        let z = origin[2] + gz as f64 * h;
                        let d2 = dxy2 + (z - p[2]) * (z - p[2]);
                        if d2 > rmax2 {
                            continue;
                        }
                        let c = base + gz as usize;
                        if (d2 as f32) < d2best[c] {
                            d2best[c] = d2 as f32;
                            nearest[c] = idx as u32;
                        }
                    }
                }
            }
        }
        let dist = d2best
            .iter()
            .map(|&d2| {
                if d2 == f32::MAX {
                    FAR
                } else {
                    ((d2.sqrt() * 10.0).round() as u32).min(254) as u8
                }
            })
            .collect();
        Self {
            origin,
            inv_h: 1.0 / h,
            nx,
            ny,
            nz,
            dist,
            nearest,
        }
    }

    /// Nearest target residue to point `x`: (quantised distance in 0.1 Å, index).
    /// Returns (FAR, u32::MAX) outside the grid or beyond rmax.
    #[inline]
    pub fn lookup(&self, x: &V3) -> (u8, u32) {
        let gx = ((x[0] - self.origin[0]) * self.inv_h).round();
        let gy = ((x[1] - self.origin[1]) * self.inv_h).round();
        let gz = ((x[2] - self.origin[2]) * self.inv_h).round();
        if gx < 0.0 || gy < 0.0 || gz < 0.0 {
            return (FAR, u32::MAX);
        }
        let (gx, gy, gz) = (gx as usize, gy as usize, gz as usize);
        if gx >= self.nx || gy >= self.ny || gz >= self.nz {
            return (FAR, u32::MAX);
        }
        let c = (gx * self.ny + gy) * self.nz + gz;
        (self.dist[c], self.nearest[c])
    }
}

/// Kernel lookup table indexed by quantised distance (0.1 Å units).
pub struct KernelLut {
    pub v: [f32; 256],
}

impl KernelLut {
    pub fn new(d0: f64) -> Self {
        let mut v = [0.0f32; 256];
        for (q, item) in v.iter_mut().enumerate().take(255) {
            let d = q as f64 / 10.0;
            *item = (1.0 / (1.0 + d * d / (d0 * d0))) as f32;
        }
        v[255] = 0.0;
        Self { v }
    }
}

#[allow(dead_code)]
pub fn brute_nearest(ca: &[V3], x: &V3) -> (f64, usize) {
    let mut best = (f64::MAX, 0);
    for (i, p) in ca.iter().enumerate() {
        let d = dist2(p, x);
        if d < best.0 {
            best = (d, i);
        }
    }
    (best.0.sqrt(), best.1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_lookup_is_close_to_brute_force() {
        let ca: Vec<V3> = (0..50)
            .map(|i| {
                let t = i as f64 * 0.5;
                [10.0 * t.cos(), 10.0 * t.sin(), 1.5 * i as f64]
            })
            .collect();
        let g = DistGrid::build(&ca, 1.0, 8.0);
        for k in 0..200 {
            let x = [
                ((k * 37) % 23) as f64 - 11.0,
                ((k * 53) % 21) as f64 - 10.0,
                (k % 70) as f64,
            ];
            let (d, _) = brute_nearest(&ca, &x);
            let (q, _) = g.lookup(&x);
            if d < 7.0 {
                assert!(q != FAR);
                assert!(((q as f64) / 10.0 - d).abs() < 0.9, "d={d} q={q}");
            }
        }
    }
}
