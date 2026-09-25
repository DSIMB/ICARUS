//! Per-structure preprocessing: secondary structure, Protein Unit hierarchy and
//! fragment signatures. Everything here depends on one structure only, so it
//! is computed once and reused for all comparisons involving that structure.

use crate::dssp;
use crate::geom::{dist, V3};
use crate::peeling::{run_peeling, PeelingConfig, SsType};
use crate::structure::Structure;

/// Fragment length used for seeds (aligned fragment pairs).
pub const FRAG: usize = 8;
/// Number of intra-fragment C-alpha distances (FRAG choose 2).
pub const FRAG_D: usize = FRAG * (FRAG - 1) / 2;

/// One node of the Protein Unit hierarchy: a contiguous residue range.
#[derive(Debug, Clone)]
pub struct PuNode {
    pub start: usize,
    pub end: usize,
    /// Peeling iteration at which the PU appeared (0 = whole chain).
    pub depth: usize,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    /// Bit k set if the node contains finest-level PU k.
    pub leaf_mask: u32,
}

impl PuNode {
    pub fn len(&self) -> usize {
        self.end - self.start + 1
    }
}

#[derive(Debug, Clone, Default)]
pub struct PuTree {
    pub nodes: Vec<PuNode>,
    /// Finest-level PUs (node indices), in sequence order.
    pub leaves: Vec<usize>,
    /// Node indices present at each Peeling iteration (index 0 = whole chain).
    pub levels: Vec<Vec<usize>>,
}

impl PuTree {
    /// Build the hierarchy from Peeling iterations, keeping iterations with at
    /// most `max_leaves` PUs. Every iteration splits exactly one PU.
    pub fn from_iterations(n: usize, iterations: &[crate::peeling::Iteration], max_leaves: usize) -> Self {
        let mut nodes = vec![PuNode { start: 0, end: n - 1, depth: 0, parent: None, children: vec![], leaf_mask: 0 }];
        let mut levels = vec![vec![0usize]];
        let mut current: Vec<usize> = vec![0];
        for (k, it) in iterations.iter().enumerate() {
            if it.pus.len() > max_leaves.max(1) {
                break;
            }
            let mut next = Vec::with_capacity(it.pus.len());
            for &[s, e] in &it.pus {
                if let Some(&id) = current.iter().find(|&&id| nodes[id].start == s && nodes[id].end == e) {
                    next.push(id);
                    continue;
                }
                let parent = current.iter().copied().find(|&id| nodes[id].start <= s && e <= nodes[id].end);
                let id = nodes.len();
                nodes.push(PuNode { start: s, end: e, depth: k + 1, parent, children: vec![], leaf_mask: 0 });
                if let Some(p) = parent {
                    nodes[p].children.push(id);
                }
                next.push(id);
            }
            next.sort_by_key(|&id| nodes[id].start);
            current = next;
            levels.push(current.clone());
        }
        let leaves = current;
        for (bit, &leaf) in leaves.iter().enumerate() {
            let (ls, le) = (nodes[leaf].start, nodes[leaf].end);
            for node in nodes.iter_mut() {
                if node.start <= ls && le <= node.end {
                    node.leaf_mask |= 1 << bit;
                }
            }
        }
        Self { nodes, leaves, levels }
    }
}

/// Parameters of the per-structure preprocessing.
#[derive(Debug, Clone)]
pub struct PrepParams {
    pub min_pu_size: usize,
    pub max_leaves: usize,
}

impl Default for PrepParams {
    fn default() -> Self {
        Self { min_pu_size: 15, max_leaves: 10 }
    }
}

/// A structure ready for comparison.
pub struct Prepared {
    pub s: Structure,
    pub ss: Vec<char>,
    pub tree: PuTree,
    /// Intra-fragment distance signatures, one per fragment start.
    pub frags: Vec<[f32; FRAG_D]>,
}

impl Prepared {
    pub fn len(&self) -> usize {
        self.s.len()
    }

    pub fn is_empty(&self) -> bool {
        self.s.is_empty()
    }
}

/// C-alpha based secondary structure (TM-align rules), used when backbone
/// atoms are missing (CA-only models).
fn ca_secondary_structure(ca: &[V3]) -> Vec<char> {
    let n = ca.len();
    let mut ss = vec![' '; n];
    for i in 2..n.saturating_sub(2) {
        let d13 = dist(&ca[i - 2], &ca[i]);
        let d14 = dist(&ca[i - 2], &ca[i + 1]);
        let d15 = dist(&ca[i - 2], &ca[i + 2]);
        let d24 = dist(&ca[i - 1], &ca[i + 1]);
        let d25 = dist(&ca[i - 1], &ca[i + 2]);
        let d35 = dist(&ca[i], &ca[i + 2]);
        let helix = (d15 - 6.37).abs() < 2.1
            && (d14 - 5.18).abs() < 2.1
            && (d25 - 5.18).abs() < 2.1
            && (d13 - 5.45).abs() < 2.1
            && (d24 - 5.45).abs() < 2.1
            && (d35 - 5.45).abs() < 2.1;
        let strand = (d15 - 13.0).abs() < 1.42
            && (d14 - 10.4).abs() < 1.42
            && (d25 - 10.4).abs() < 1.42
            && (d13 - 6.1).abs() < 1.42
            && (d24 - 6.1).abs() < 1.42
            && (d35 - 6.1).abs() < 1.42;
        ss[i] = if helix { 'H' } else if strand { 'E' } else { ' ' };
    }
    ss
}

pub fn fragment_signatures(ca: &[V3]) -> Vec<[f32; FRAG_D]> {
    if ca.len() < FRAG {
        return Vec::new();
    }
    (0..=ca.len() - FRAG)
        .map(|i| {
            let mut sig = [0f32; FRAG_D];
            let mut k = 0;
            for a in 0..FRAG {
                for b in (a + 1)..FRAG {
                    sig[k] = dist(&ca[i + a], &ca[i + b]) as f32;
                    k += 1;
                }
            }
            sig
        })
        .collect()
}

pub fn prepare(s: Structure, params: &PrepParams) -> Prepared {
    let n = s.len();
    let has_backbone = s.backbone.iter().filter(|b| b.n.is_some() && b.c.is_some() && b.o.is_some()).count() * 2 > n;
    let ss = if has_backbone { dssp::assign(&s.backbone, &s.seq) } else { ca_secondary_structure(&s.ca) };
    let ss_types: Vec<SsType> = ss.iter().map(|&c| SsType::from_dssp_char(c)).collect();
    let cfg = PeelingConfig {
        min_pu_size: params.min_pu_size,
        // Only the first iterations are used: stop once the leaf budget is exceeded.
        max_pu_number: params.max_leaves + 2,
        ..PeelingConfig::default()
    };
    let iterations = if n >= 2 * params.min_pu_size { run_peeling(&s.ca, &ss_types, &cfg) } else { Vec::new() };
    let tree = PuTree::from_iterations(n, &iterations, params.max_leaves);
    let frags = fragment_signatures(&s.ca);
    Prepared { s, ss, tree, frags }
}
