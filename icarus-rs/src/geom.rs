//! Geometry: vectors, rigid transforms and weighted optimal superposition.

pub type V3 = [f64; 3];

#[inline]
pub fn sub(a: &V3, b: &V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
pub fn dist2(a: &V3, b: &V3) -> f64 {
    let d = sub(a, b);
    d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
}

#[inline]
pub fn dist(a: &V3, b: &V3) -> f64 {
    dist2(a, b).sqrt()
}

/// Rigid transform x -> R x + t.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform {
    pub r: [[f64; 3]; 3],
    pub t: V3,
}

impl Default for Transform {
    fn default() -> Self {
        Self::identity()
    }
}

impl Transform {
    pub fn identity() -> Self {
        Self { r: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], t: [0.0; 3] }
    }

    #[inline]
    pub fn apply(&self, x: &V3) -> V3 {
        let r = &self.r;
        [
            r[0][0] * x[0] + r[0][1] * x[1] + r[0][2] * x[2] + self.t[0],
            r[1][0] * x[0] + r[1][1] * x[1] + r[1][2] * x[2] + self.t[1],
            r[2][0] * x[0] + r[2][1] * x[1] + r[2][2] * x[2] + self.t[2],
        ]
    }

    pub fn apply_all(&self, xs: &[V3]) -> Vec<V3> {
        xs.iter().map(|x| self.apply(x)).collect()
    }

    /// Inverse transform.
    pub fn inverse(&self) -> Self {
        let r = &self.r;
        let rt = [[r[0][0], r[1][0], r[2][0]], [r[0][1], r[1][1], r[2][1]], [r[0][2], r[1][2], r[2][2]]];
        let t = [
            -(rt[0][0] * self.t[0] + rt[0][1] * self.t[1] + rt[0][2] * self.t[2]),
            -(rt[1][0] * self.t[0] + rt[1][1] * self.t[1] + rt[1][2] * self.t[2]),
            -(rt[2][0] * self.t[0] + rt[2][1] * self.t[1] + rt[2][2] * self.t[2]),
        ];
        Self { r: rt, t }
    }

    /// Rotation angle (radians) of R * other.R^T, a distance between rotations.
    pub fn rotation_angle_to(&self, other: &Transform) -> f64 {
        let a = &self.r;
        let b = &other.r;
        let mut tr = 0.0;
        for i in 0..3 {
            for k in 0..3 {
                tr += a[i][k] * b[i][k];
            }
        }
        ((tr - 1.0) / 2.0).clamp(-1.0, 1.0).acos()
    }
}

/// Largest eigenpair of a symmetric 4x4 matrix (cyclic Jacobi).
fn sym4_max_eigen(mut a: [[f64; 4]; 4]) -> (f64, [f64; 4]) {
    let mut v = [[0.0f64; 4]; 4];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _sweep in 0..32 {
        let mut off = 0.0;
        for p in 0..4 {
            for q in (p + 1)..4 {
                off += a[p][q] * a[p][q];
            }
        }
        let scale = a[0][0].abs() + a[1][1].abs() + a[2][2].abs() + a[3][3].abs();
        if off <= 1e-22 * (scale * scale + 1e-300) {
            break;
        }
        for p in 0..4 {
            for q in (p + 1)..4 {
                let apq = a[p][q];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..4 {
                    let akp = a[k][p];
                    let akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..4 {
                    let apk = a[p][k];
                    let aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let vkp = row[p];
                    let vkq = row[q];
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut best = 0;
    for i in 1..4 {
        if a[i][i] > a[best][best] {
            best = i;
        }
    }
    (a[best][best], [v[0][best], v[1][best], v[2][best], v[3][best]])
}

/// Optimal (weighted) superposition of `x` onto `y`: returns the transform T
/// minimising sum_i w_i |T x_i - y_i|^2, and the weighted RMSD.
/// Horn's quaternion method.
pub fn superpose(x: &[V3], y: &[V3], w: Option<&[f64]>) -> (Transform, f64) {
    let n = x.len().min(y.len());
    if n == 0 {
        return (Transform::identity(), 0.0);
    }
    let mut sw = 0.0;
    let mut cx = [0.0; 3];
    let mut cy = [0.0; 3];
    for i in 0..n {
        let wi = w.map_or(1.0, |w| w[i]);
        sw += wi;
        for k in 0..3 {
            cx[k] += wi * x[i][k];
            cy[k] += wi * y[i][k];
        }
    }
    if sw <= 0.0 {
        return (Transform::identity(), 0.0);
    }
    for k in 0..3 {
        cx[k] /= sw;
        cy[k] /= sw;
    }
    let mut s = [[0.0f64; 3]; 3];
    let mut ex = 0.0;
    for i in 0..n {
        let wi = w.map_or(1.0, |w| w[i]);
        let a = sub(&x[i], &cx);
        let b = sub(&y[i], &cy);
        ex += wi * (a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
        for p in 0..3 {
            for q in 0..3 {
                s[p][q] += wi * a[p] * b[q];
            }
        }
    }
    let (sxx, sxy, sxz) = (s[0][0], s[0][1], s[0][2]);
    let (syx, syy, syz) = (s[1][0], s[1][1], s[1][2]);
    let (szx, szy, szz) = (s[2][0], s[2][1], s[2][2]);
    let nmat = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];
    let (lambda, q) = sym4_max_eigen(nmat);
    let (q0, q1, q2, q3) = (q[0], q[1], q[2], q[3]);
    let r = [
        [q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2.0 * (q1 * q2 - q0 * q3), 2.0 * (q1 * q3 + q0 * q2)],
        [2.0 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2.0 * (q2 * q3 - q0 * q1)],
        [2.0 * (q1 * q3 - q0 * q2), 2.0 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3],
    ];
    let rc = [
        r[0][0] * cx[0] + r[0][1] * cx[1] + r[0][2] * cx[2],
        r[1][0] * cx[0] + r[1][1] * cx[1] + r[1][2] * cx[2],
        r[2][0] * cx[0] + r[2][1] * cx[1] + r[2][2] * cx[2],
    ];
    let t = [cy[0] - rc[0], cy[1] - rc[1], cy[2] - rc[2]];
    let msd = ((ex - 2.0 * lambda) / sw).max(0.0);
    (Transform { r, t }, msd.sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn superpose_recovers_known_transform() {
        let x: Vec<V3> = vec![[0.0, 0.0, 0.0], [1.5, 0.2, 0.0], [2.0, 1.4, 0.7], [0.3, 2.2, 1.9], [-1.0, 0.5, 3.0]];
        let (c, s) = (0.6f64.cos(), 0.6f64.sin());
        let t = Transform { r: [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], t: [3.0, -2.0, 5.0] };
        let y: Vec<V3> = x.iter().map(|p| t.apply(p)).collect();
        let (est, rmsd) = superpose(&x, &y, None);
        assert!(rmsd < 1e-6, "rmsd {rmsd}");
        for p in &x {
            assert!(dist(&est.apply(p), &t.apply(p)) < 1e-6);
        }
        let inv = est.inverse();
        assert!(dist(&inv.apply(&est.apply(&x[2])), &x[2]) < 1e-9);
    }

    #[test]
    fn superpose_handles_reflection_free_degenerate_input() {
        let x: Vec<V3> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let y: Vec<V3> = vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]];
        let (_, rmsd) = superpose(&x, &y, None);
        assert!(rmsd < 1e-6);
    }
}
