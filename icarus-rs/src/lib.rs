//! ICARUS 2: fast flexible protein structural alignment based on Protein Units.

// Numerical kernels index several parallel arrays in the same loop.
#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::type_complexity
)]

pub mod align;
pub mod db;
pub mod dp;
pub mod dssp;
pub mod geom;
pub mod grid;
pub mod output;
pub mod peeling;
pub mod prep;
pub mod seeds;
pub mod structure;
pub mod tm;
