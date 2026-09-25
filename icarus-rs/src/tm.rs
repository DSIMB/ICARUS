//! TM-score scale parameters.

/// TM-score d0 for a normalisation length (Zhang & Skolnick 2004, as in TM-align).
pub fn d0(l: usize) -> f64 {
    if l <= 21 {
        0.5
    } else {
        (1.24 * ((l as f64) - 15.0).cbrt() - 1.8).max(0.5)
    }
}

/// d0 used to *search* alignments (TM-align / gdt2.pl convention):
/// d0 + 0.8 clamped to [4.5, 8] Å.
pub fn d0_search(l: usize) -> f64 {
    let base = if l <= 19 {
        0.168
    } else {
        (1.24 * ((l as f64) - 15.0).cbrt() - 1.8).max(0.5)
    };
    (base + 0.8).clamp(4.5, 8.0)
}

/// TM kernel 1 / (1 + d²/d0²) from a squared distance.
#[inline]
pub fn kernel(d2: f64, inv_d02: f64) -> f64 {
    1.0 / (1.0 + d2 * inv_d02)
}
