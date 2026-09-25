//! Protein Peeling: hierarchical decomposition of a structure into Protein Units.
//!
//! Vendored from SWORD3's native Rust port (see the individual files).

mod algorithm;
mod contact_matrix;
mod interval;

pub use algorithm::{build_cutting_mask, run_peeling, Iteration, PeelingConfig, SsType};
