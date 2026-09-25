//! Contact probability matrix with 2D cumulative sum table.
//!
//! Vendored from SWORD3 (DSIMB/sword3 @ e47f837,
//! `sword3-lib/src/peeling/contact_matrix.rs`, CeCILL-2.1). ICARUS keeps only
//! the prefix-sum table (the raw matrix is never read by the cutting code),
//! halving memory for long chains, and builds it serially: structures are
//! processed in parallel at a higher level.
//!
//! The contact probability between two C-alpha atoms is defined by the
//! logistic function:
//!
//! $$p(d) = \frac{1}{1 + e^{(d - D_0) / \Delta}}$$

/// Contact probability matrix represented by its (N+1)×(N+1) prefix sums.
#[derive(Debug, Clone)]
pub struct ContactMatrix {
    n: usize,
    /// cum[i][j] = sum of p over rows 0..i-1 and columns 0..j-1.
    cum: Vec<f64>,
}

impl ContactMatrix {
    pub fn from_ca_coords(ca_coords: &[[f64; 3]], d0: f64, delta: f64) -> Self {
        let n = ca_coords.len();
        let cn = n + 1;
        let mut cum = vec![0.0f64; cn * cn];
        let inv_delta = 1.0 / delta;
        let mut row = vec![0.0f64; n];
        for i in 0..n {
            let (xi, yi, zi) = (ca_coords[i][0], ca_coords[i][1], ca_coords[i][2]);
            for (j, r) in row.iter_mut().enumerate() {
                let dx = xi - ca_coords[j][0];
                let dy = yi - ca_coords[j][1];
                let dz = zi - ca_coords[j][2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                *r = 1.0 / (1.0 + ((dist - d0) * inv_delta).exp());
            }
            // cum[i+1][j+1] = row prefix + cum[i][j+1]
            let (prev, cur) = cum.split_at_mut((i + 1) * cn);
            let prev = &prev[i * cn..];
            let mut acc = 0.0;
            for j in 0..n {
                acc += row[j];
                cur[j + 1] = acc + prev[j + 1];
            }
        }
        Self { n, cum }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.n
    }

    /// Contact probability between residues i and j (0-indexed).
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.rectangle_sum(i, j, i, j)
    }

    /// Sum of contact probabilities over [row_start..=row_end, col_start..=col_end].
    #[inline]
    pub fn rectangle_sum(&self, row_start: usize, col_start: usize, row_end: usize, col_end: usize) -> f64 {
        let (r1, c1, r2, c2) = (row_start, col_start, row_end + 1, col_end + 1);
        let cn = self.n + 1;
        self.cum[r2 * cn + c2] - self.cum[r1 * cn + c2] - self.cum[r2 * cn + c1] + self.cum[r1 * cn + c1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rectangle_sums_match_brute_force() {
        let coords = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [5.0, 5.0, 0.0]];
        let mat = ContactMatrix::from_ca_coords(&coords, 6.0, 1.5);
        let p = |i: usize, j: usize| {
            let d = ((coords[i][0] - coords[j][0]).powi(2) + (coords[i][1] - coords[j][1]).powi(2)).sqrt();
            1.0 / (1.0 + ((d - 6.0) / 1.5f64).exp())
        };
        assert!((mat.get(1, 2) - p(1, 2)).abs() < 1e-12);
        let brute = p(1, 0) + p(1, 1) + p(2, 0) + p(2, 1);
        assert!((mat.rectangle_sum(1, 0, 2, 1) - brute).abs() < 1e-12);
    }
}
