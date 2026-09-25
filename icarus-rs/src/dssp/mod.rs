//! DSSP secondary structure assignment (Kabsch & Sander, 1983).
//!
//! The H-bond, bridge and helix code is vendored from SWORD3 (pure Rust port of
//! `dsspcmbi`). This adapter builds the DSSP backbone chain from ICARUS residues
//! and returns one DSSP symbol per residue. Solvent accessibility and bend
//! angles are not computed: Protein Peeling only needs helix/strand segments.

mod bridge;
mod hbond;
mod helix;
mod types;

use types::{BackboneResidue, DsspChain, Point3D, BREAKDIST};

/// Backbone atoms of one residue (any may be missing).
#[derive(Debug, Clone, Copy, Default)]
pub struct Backbone {
    pub n: Option<[f64; 3]>,
    pub ca: Option<[f64; 3]>,
    pub c: Option<[f64; 3]>,
    pub o: Option<[f64; 3]>,
}

fn p(v: [f64; 3]) -> Point3D {
    Point3D::new(v[0], v[1], v[2])
}

/// Assign DSSP symbols ('H','G','I','E','B','T','S',' ') to each residue.
/// Residues with an incomplete backbone get ' ' and introduce a chain break.
pub fn assign(backbone: &[Backbone], aa: &[u8]) -> Vec<char> {
    let mut chain = DsspChain::new();
    let mut skipped = false;
    for (idx, bb) in backbone.iter().enumerate() {
        let (Some(n), Some(ca), Some(c), Some(o)) = (bb.n, bb.ca, bb.c, bb.o) else {
            skipped = true;
            continue;
        };
        if skipped && chain.len > 0 && chain.get(chain.len).aa != '!' {
            chain.push(BackboneResidue::chain_break());
        }
        skipped = false;
        let mut res = BackboneResidue {
            aa: aa.get(idx).map(|&b| b as char).unwrap_or('X'),
            source_index: Some(idx),
            n: p(n),
            ca: p(ca),
            c: p(c),
            o: p(o),
            ..BackboneResidue::default()
        };
        // Synthesize the amide H as in DSSP (none for Pro or after a break).
        res.h = res.n;
        res.has_h = false;
        if chain.len > 0 {
            let prev = chain.get(chain.len);
            if prev.aa != '!' && prev.c.distance_to(&res.n) > BREAKDIST {
                chain.push(BackboneResidue::chain_break());
            }
            let prev = chain.get(chain.len);
            if prev.aa != '!' && res.aa != 'P' {
                let co = prev.c - prev.o;
                let l = (co.x * co.x + co.y * co.y + co.z * co.z).sqrt();
                if l > 0.0 {
                    res.h =
                        Point3D::new(res.n.x + co.x / l, res.n.y + co.y / l, res.n.z + co.z / l);
                    res.has_h = true;
                }
            }
        }
        chain.push(res);
    }
    let mut out = vec![' '; backbone.len()];
    if chain.len == 0 {
        return out;
    }
    hbond::flag_hydrogen_bonds(&mut chain);
    bridge::flag_bridges(&mut chain);
    helix::flag_turns(&mut chain);
    for i in 1..=chain.len {
        let r = chain.get(i);
        if let Some(src) = r.source_index {
            out[src] = r.ss[0];
        }
    }
    out
}
