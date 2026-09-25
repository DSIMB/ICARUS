//! One inclusive interval convention for Peeling/CA projection positions.
//!
//! Vendored from SWORD3 (DSIMB/sword3 @ e47f837, `sword3-lib/src/peeling/interval.rs`, CeCILL-2.1).

use anyhow::{ensure, Result};

/// A closed interval of zero-based Peeling/CA projection positions.
///
/// Both endpoints are inclusive, so the number of residues covered is
/// `end - start + 1`. This is the only size convention in the Peeling module.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PuSpan {
    start: usize,
    end: usize,
}

impl PuSpan {
    pub fn new(start: usize, end: usize) -> Result<Self> {
        ensure!(start <= end, "PU span start {start} exceeds end {end}");
        Ok(Self { start, end })
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn end(&self) -> usize {
        self.end
    }

    /// Number of residues covered, counting both endpoints.
    pub fn len(&self) -> usize {
        self.end - self.start + 1
    }

    pub fn is_empty(&self) -> bool {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn span_len_counts_both_endpoints() {
        assert_eq!(PuSpan::new(0, 0).unwrap().len(), 1);
        assert_eq!(PuSpan::new(3, 7).unwrap().len(), 5);
    }

    #[test]
    fn span_rejects_reversed_endpoints() {
        let error = PuSpan::new(7, 3).unwrap_err();
        assert!(format!("{error:#}").contains("exceeds end"));
    }
}
