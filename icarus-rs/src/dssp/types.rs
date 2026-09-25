//! Core data structures for DSSP secondary structure assignment.
//!
//! Vendored from SWORD3 (DSIMB/sword3 @ e47f837, `sword3-lib/src/dssp/types.rs`,
//! CeCILL-2.1) and reduced to the fields used for secondary-structure
//! assignment (accessibility and legacy formatting fields removed).

/// Physical constants from the original DSSP algorithm.
pub const BREAKDIST: f64 = 2.5; // Max peptide bond C-N distance (Angstrom)
pub const CADIST: f64 = 9.0; // CA distance cutoff for H-bond candidates
pub const Q: f64 = -27888.0; // Electrostatic coupling constant (cal/mol)
pub const HBLOW: i64 = -9900; // Min H-bond energy (cal/mol)
pub const HBHIGH: i64 = -500; // Max H-bond energy (cal/mol)
pub const DIST_MIN: f64 = 0.5; // Smallest allowed distance between atoms

/// Maximum number of bridges.
pub const MAXBRIDGE: usize = 2000;

/// A 3D point in Angstroms.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3D {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn distance_to(&self, other: &Point3D) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

impl std::ops::Sub for Point3D {
    type Output = Point3D;
    fn sub(self, rhs: Point3D) -> Point3D {
        Point3D::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

/// A backbone residue with coordinates and DSSP state.
#[derive(Debug, Clone)]
pub struct BackboneResidue {
    /// One-letter amino acid code ('!' for chain break).
    pub aa: char,
    /// Index of this residue in the caller's residue array (None for breaks).
    pub source_index: Option<usize>,
    /// Backbone atom coordinates.
    pub n: Point3D,
    pub ca: Point3D,
    pub c: Point3D,
    pub o: Point3D,
    /// Synthesized backbone H position.
    pub h: Point3D,
    /// Whether this residue has a valid amide H (false for Pro, chain start).
    pub has_h: bool,
    /// Secondary structure columns:
    /// [0]=symbol, [1]=turn3, [2]=turn4, [3]=turn5, [4]=bend, [5]=chirality, [6]=beta1, [7]=beta2
    pub ss: [char; 8],
    /// Bridge partners for beta1 and beta2.
    pub partner: [usize; 2],
    /// Sheet label character.
    pub sheet_label: char,
    /// Best 2 acceptor H-bonds (this residue's NH donates to acceptor's CO).
    pub acceptor: [HydrogenBond; 2],
    /// Best 2 donor H-bonds (donor's NH donates to this residue's CO).
    pub donor: [HydrogenBond; 2],
    /// Virtual bend angle kappa (not computed here: bends are not needed).
    pub kappa: f64,
}

impl BackboneResidue {
    /// Create a chain break marker residue.
    pub fn chain_break() -> Self {
        Self {
            aa: '!',
            ..Self::default()
        }
    }
}

impl Default for BackboneResidue {
    fn default() -> Self {
        Self {
            aa: ' ',
            source_index: None,
            n: Point3D::default(),
            ca: Point3D::default(),
            c: Point3D::default(),
            o: Point3D::default(),
            h: Point3D::default(),
            has_h: false,
            ss: [' '; 8],
            partner: [0; 2],
            sheet_label: ' ',
            acceptor: [HydrogenBond::default(); 2],
            donor: [HydrogenBond::default(); 2],
            kappa: 360.0,
        }
    }
}

/// A hydrogen bond with residue index and energy.
#[derive(Debug, Clone, Copy, Default)]
pub struct HydrogenBond {
    /// 1-based residue index (0 = no bond).
    pub residue: usize,
    /// Energy in cal/mol.
    pub energy: i64,
}

/// Bridge type between two residues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BridgeType {
    Parallel,
    Antiparallel,
    NoBridge,
}

/// A beta-bridge/ladder entry in the bridge table.
#[derive(Debug, Clone)]
pub struct Bridge {
    pub sheet_name: char,
    pub ladder_name: char,
    pub btype: BridgeType,
    /// Set of linked ladder indices.
    pub link_set: Vec<bool>,
    /// i-strand begin/end (1-based).
    pub ib: usize,
    pub ie: usize,
    /// j-strand begin/end (1-based).
    pub jb: usize,
    pub je: usize,
    /// Link to previous/next ladder in the chain.
    pub from: usize,
    pub towards: usize,
}

impl Bridge {
    pub fn new(max_bridges: usize) -> Self {
        Self {
            sheet_name: ' ',
            ladder_name: ' ',
            btype: BridgeType::NoBridge,
            link_set: vec![false; max_bridges + 1],
            ib: 0,
            ie: 0,
            jb: 0,
            je: 0,
            from: 0,
            towards: 0,
        }
    }
}

/// The DSSP chain: 1-indexed array of backbone residues.
pub struct DsspChain {
    /// Residues, 1-indexed. Index 0 is unused (sentinel).
    pub residues: Vec<BackboneResidue>,
    /// Number of actual residues (length of chain, excluding index 0).
    pub len: usize,
}

impl Default for DsspChain {
    fn default() -> Self {
        Self::new()
    }
}

impl DsspChain {
    pub fn new() -> Self {
        Self {
            residues: vec![BackboneResidue::default()], // index 0 sentinel
            len: 0,
        }
    }

    /// Get residue at 1-based index.
    pub fn get(&self, i: usize) -> &BackboneResidue {
        &self.residues[i]
    }

    /// Get mutable residue at 1-based index.
    pub fn get_mut(&mut self, i: usize) -> &mut BackboneResidue {
        &mut self.residues[i]
    }

    /// Push a residue onto the chain. Returns its 1-based index.
    pub fn push(&mut self, res: BackboneResidue) -> usize {
        self.residues.push(res);
        self.len = self.residues.len() - 1;
        self.len
    }

    /// Check if there's no chain break between residues i and j (inclusive, 1-based).
    pub fn no_chain_break(&self, i: usize, j: usize) -> bool {
        if i < 1 || j > self.len || i > j {
            return false;
        }
        for k in i..=j {
            if self.residues[k].aa == '!' {
                return false;
            }
        }
        true
    }
}
