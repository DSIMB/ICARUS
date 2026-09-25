//! End-to-end checks of the flexible alignment on cases with a known answer.

use std::path::Path;

use icarus::align::{align_pair, AlignParams};
use icarus::dp::DpWork;
use icarus::geom::{Transform, V3};
use icarus::prep::{prepare, PrepParams, Prepared};
use icarus::structure::{read_structure, Structure};

fn load(name: &str) -> Structure {
    read_structure(
        &Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name),
        None,
        false,
    )
    .unwrap()
}

fn prep(s: Structure) -> Prepared {
    prepare(s, &PrepParams::default())
}

fn run(a: Structure, b: Structure) -> icarus::align::PairResult {
    let mut w = DpWork::default();
    align_pair(&prep(a), &prep(b), &AlignParams::default(), &mut w)
}

/// Rotation by `deg` degrees around the axis through `p` along `u` (unit).
fn rotation_about(p: V3, u: V3, deg: f64) -> Transform {
    let (s, c) = deg.to_radians().sin_cos();
    let (x, y, z) = (u[0], u[1], u[2]);
    let r = [
        [
            c + x * x * (1.0 - c),
            x * y * (1.0 - c) - z * s,
            x * z * (1.0 - c) + y * s,
        ],
        [
            y * x * (1.0 - c) + z * s,
            c + y * y * (1.0 - c),
            y * z * (1.0 - c) - x * s,
        ],
        [
            z * x * (1.0 - c) - y * s,
            z * y * (1.0 - c) + x * s,
            c + z * z * (1.0 - c),
        ],
    ];
    let rp = [
        r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2],
        r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2],
        r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2],
    ];
    Transform {
        r,
        t: [p[0] - rp[0], p[1] - rp[1], p[2] - rp[2]],
    }
}

fn transform_residues(s: &mut Structure, from: usize, tr: &Transform) {
    for i in from..s.len() {
        s.ca[i] = tr.apply(&s.ca[i]);
        let b = &mut s.backbone[i];
        for a in [&mut b.n, &mut b.ca, &mut b.c, &mut b.o]
            .into_iter()
            .flatten()
        {
            *a = tr.apply(a);
        }
    }
}

#[test]
fn self_alignment_is_perfect() {
    let r = run(load("d1nls__.pdb"), load("d1nls__.pdb"));
    assert!(r.tm_flex() > 0.99, "flexible TM {}", r.tm_flex());
    assert!(r.tm_rigid() > 0.99, "rigid TM {}", r.tm_rigid());
}

#[test]
fn hinge_motion_is_recovered() {
    let a = load("d1nls__.pdb");
    let mut b = a.clone();
    let h = a.len() / 2;
    let axis = {
        let d = [
            a.ca[h + 1][0] - a.ca[h - 1][0],
            a.ca[h + 1][1] - a.ca[h - 1][1],
            a.ca[h + 1][2] - a.ca[h - 1][2],
        ];
        let n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        // an axis roughly perpendicular to the chain at the hinge
        let v = [d[1] / n, -d[0] / n, 0.0];
        let nv = (v[0] * v[0] + v[1] * v[1]).sqrt();
        [v[0] / nv, v[1] / nv, 0.0]
    };
    transform_residues(&mut b, h, &rotation_about(a.ca[h], axis, 50.0));
    let r = run(a, b);
    assert!(
        r.tm_rigid() < 0.8,
        "rigid TM should drop, got {}",
        r.tm_rigid()
    );
    assert!(
        r.tm_flex() > 0.9,
        "flexible TM {} (rigid {})",
        r.tm_flex(),
        r.tm_rigid()
    );
    assert!(r.flex.segs.len() >= 2);
}

#[test]
fn circular_permutation_is_recovered() {
    let a = load("d1nls__.pdb");
    let k = 100;
    let mut b = a.clone();
    let rot = |v: &mut Vec<_>| v.rotate_left(k);
    rot(&mut b.ca);
    b.seq.rotate_left(k);
    b.resn.rotate_left(k);
    b.resid.rotate_left(k);
    b.bfac.rotate_left(k);
    b.backbone.rotate_left(k);
    let r = run(a, b);
    assert!(
        r.tm_flex() > 0.9,
        "flexible TM {} (rigid {})",
        r.tm_flex(),
        r.tm_rigid()
    );
}

#[test]
fn lectin_circular_permutation_pair() {
    // concanavalin A vs pea lectin (RIPC): related by a circular permutation
    let r = run(load("d1nls__.pdb"), load("d2bqpa_.pdb"));
    assert!(r.tm_rigid() < 0.6, "rigid TM {}", r.tm_rigid());
    assert!(r.tm_flex() > 0.85, "flexible TM {}", r.tm_flex());
}

#[test]
fn zero_budgets_do_not_panic() {
    let p = AlignParams {
        per_node: 0,
        max_seeds: 0,
        ..AlignParams::default()
    };
    let mut w = DpWork::default();
    let r = align_pair(
        &prep(load("d1nls__.pdb")),
        &prep(load("d2bqpa_.pdb")),
        &p,
        &mut w,
    );
    assert!(
        r.tm_flex() >= 0.0 && r.tm_flex() <= 1.0,
        "flexible TM {}",
        r.tm_flex()
    );
}
