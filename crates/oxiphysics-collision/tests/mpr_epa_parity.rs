// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Parity tests: full XenoCollide MPR vs EPA.
//!
//! Validates that Minkowski Portal Refinement produces penetration depths
//! consistent with EPA on overlapping convex pairs, that it falls back
//! correctly on degenerate simplices, and that the smoke cases match the
//! in-crate unit tests. The exhaustive 1e5-pair parity sweep is env-gated
//! (`MPR_FULL_PARITY=1`) so CI runs a fast 1k-pair subset by default.

use oxiphysics_collision::narrowphase::gjk::mpr::{mpr_contact, mpr_full};
use oxiphysics_collision::narrowphase::gjk::{GjkResult, MprResult};
use oxiphysics_collision::{Epa, Gjk};
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Vec3, quat_from_axis_angle};
use oxiphysics_geometry::{BoxShape, Shape, Sphere};

/// EPA ground-truth penetration depth for an overlapping pair, or `None` when
/// GJK reports the shapes as separated or EPA fails to converge.
fn epa_depth(a: &dyn Shape, ta: &Transform, b: &dyn Shape, tb: &Transform) -> Option<f64> {
    match Gjk::query(a, ta, b, tb) {
        GjkResult::Intersecting(s) => Epa::penetration_depth(a, ta, b, tb, &s).map(|c| c.depth),
        GjkResult::Separated { .. } => None,
    }
}

/// MPR penetration depth for an overlapping pair, or `None` when MPR reports
/// the shapes as separated.
fn mpr_depth(a: &dyn Shape, ta: &Transform, b: &dyn Shape, tb: &Transform) -> Option<f64> {
    match mpr_full(a, ta, b, tb) {
        MprResult::Intersecting { depth, .. } => Some(depth),
        MprResult::Separated => None,
    }
}

/// Deterministic SplitMix64 RNG so the parity sweeps reproduce bit-for-bit.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg(seed)
    }

    fn next_u64(&mut self) -> u64 {
        // SplitMix64
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

// ---------------------------------------------------------------------------
// TEST 1 — smoke cases
// ---------------------------------------------------------------------------

#[test]
fn mpr_smoke_sphere_sphere_intersecting() {
    let a = Sphere::new(1.0);
    let b = Sphere::new(1.0);
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));

    match mpr_full(&a, &ta, &b, &tb) {
        MprResult::Intersecting { depth, .. } => {
            assert!(depth >= 0.0, "MPR depth must be non-negative, got {depth}");
        }
        MprResult::Separated => panic!("unit spheres at distance 1 must intersect"),
    }

    let epa = epa_depth(&a, &ta, &b, &tb).expect("EPA must find overlap for unit spheres");
    let m = mpr_depth(&a, &ta, &b, &tb).expect("MPR must find overlap for unit spheres");
    assert!(
        (m - epa).abs() <= 1e-2 * epa.max(1e-6) + 1e-3,
        "MPR/EPA depth mismatch: mpr={m} epa={epa}"
    );
    assert!(
        (m - 1.0).abs() <= 5e-2 * 1.0,
        "MPR depth not within 5% of true overlap 1.0: {m}"
    );
}

#[test]
fn mpr_smoke_sphere_sphere_separated() {
    let a = Sphere::new(0.5);
    let b = Sphere::new(0.5);
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));

    assert!(
        matches!(mpr_full(&a, &ta, &b, &tb), MprResult::Separated),
        "distant unit spheres must be separated"
    );
}

#[test]
fn mpr_smoke_coincident() {
    let a = Sphere::new(1.0);
    let b = Sphere::new(1.0);
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::zeros());

    assert!(
        matches!(mpr_full(&a, &ta, &b, &tb), MprResult::Intersecting { .. }),
        "coincident spheres must intersect"
    );
}

#[test]
fn mpr_smoke_box_box() {
    let a = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::new(1.5, 0.0, 0.0));

    assert!(
        matches!(mpr_full(&a, &ta, &b, &tb), MprResult::Intersecting { .. }),
        "overlapping boxes must intersect"
    );

    let epa = epa_depth(&a, &ta, &b, &tb).expect("EPA must find overlap for boxes");
    let m = mpr_depth(&a, &ta, &b, &tb).expect("MPR must find overlap for boxes");
    assert!(
        (m - epa).abs() <= 5e-2 * epa.max(1e-6) + 1e-3,
        "MPR/EPA box depth mismatch: mpr={m} epa={epa}"
    );
}

// ---------------------------------------------------------------------------
// TEST 2 — degenerate / fallback robustness
// ---------------------------------------------------------------------------

#[test]
fn mpr_degenerate_face_aligned_boxes() {
    let a = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::new(1.9, 0.0, 0.0));

    match mpr_full(&a, &ta, &b, &tb) {
        MprResult::Intersecting { normal, depth, .. } => {
            assert!(
                depth > 0.0 && depth.is_finite(),
                "face-aligned box depth must be positive and finite, got {depth}"
            );
            assert!(
                (normal.norm() - 1.0).abs() < 1e-3,
                "MPR normal must be unit length, got {}",
                normal.norm()
            );
            assert!(
                normal.x.abs() > 0.9,
                "face-aligned normal must point along x, got {normal:?}"
            );
        }
        MprResult::Separated => panic!("face-aligned overlapping boxes must intersect"),
    }
}

#[test]
fn mpr_contact_returns_some_on_overlap() {
    let a = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::new(1.9, 0.0, 0.0));

    let c = mpr_contact(&a, &ta, &b, &tb).expect("overlapping boxes must produce a contact");
    assert!(
        c.depth > 0.0,
        "contact depth must be positive, got {}",
        c.depth
    );
    assert!(
        c.point_a.norm().is_finite() && c.point_b.norm().is_finite(),
        "contact witness points must be finite: a={:?} b={:?}",
        c.point_a,
        c.point_b
    );
}

#[test]
fn mpr_fallback_path_smoke() {
    let a = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let ta = Transform::from_position(Vec3::zeros());
    let tb = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));

    assert!(
        matches!(mpr_full(&a, &ta, &b, &tb), MprResult::Intersecting { .. }),
        "deeply overlapping boxes must intersect"
    );
}

// ---------------------------------------------------------------------------
// TEST 3 — random parity vs EPA (deterministic)
// ---------------------------------------------------------------------------

/// Drive `n` random convex pairs through both EPA and MPR, returning the
/// number of comparable overlaps, the count classified intersecting, and the
/// worst observed relative depth error. Panics if MPR reports `Separated`
/// where EPA found a genuine overlap.
fn run_parity(n: usize, seed: u64) -> (usize, usize, f64) {
    let mut rng = Lcg::new(seed);
    let mut compared = 0usize;
    let mut intersecting = 0usize;
    let mut max_rel_err = 0.0f64;

    for _ in 0..n {
        let kind = (rng.next_u64() % 3) as u8;
        let ta = Transform::from_position(Vec3::zeros());
        let pos_b = Vec3::new(
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        );

        let (e, m) = match kind {
            0 => {
                let a = Sphere::new(rng.range(0.5, 2.0));
                let b = Sphere::new(rng.range(0.5, 2.0));
                let tb = Transform::from_position(pos_b);
                (epa_depth(&a, &ta, &b, &tb), mpr_depth(&a, &ta, &b, &tb))
            }
            1 => {
                let a = Sphere::new(rng.range(0.5, 2.0));
                let hb = Vec3::new(
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                );
                let b = BoxShape::new(hb);
                let tb = Transform::from_position(pos_b);
                (epa_depth(&a, &ta, &b, &tb), mpr_depth(&a, &ta, &b, &tb))
            }
            _ => {
                let ha = Vec3::new(
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                );
                let hb = Vec3::new(
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                );
                let a = BoxShape::new(ha);
                let b = BoxShape::new(hb);
                let tb = Transform::from_position(pos_b);
                (epa_depth(&a, &ta, &b, &tb), mpr_depth(&a, &ta, &b, &tb))
            }
        };

        // Robustness: MPR must not miss an overlap EPA found above threshold.
        if let Some(ev) = e {
            if ev > 1e-4 {
                assert!(
                    m.is_some(),
                    "MPR separated where EPA found overlap ev={ev} kind={kind} pos_b={pos_b:?}"
                );
            }
        }

        // Parity statistics over comparable overlaps.
        if let (Some(ev), Some(mv)) = (e, m) {
            if ev > 1e-4 {
                compared += 1;
                intersecting += 1;
                let rel = (mv - ev).abs() / ev;
                if rel > max_rel_err {
                    max_rel_err = rel;
                }
            }
        }
    }

    eprintln!(
        "run_parity(n={n}, seed={seed:#x}): compared={compared} intersecting={intersecting} max_rel_err={max_rel_err}"
    );
    (compared, intersecting, max_rel_err)
}

#[test]
fn mpr_depth_parity_random_1k() {
    let (compared, _inter, max_err) = run_parity(1000, 0x00C0_FFEE);
    assert!(compared > 100, "too few overlapping samples: {compared}");
    assert!(
        max_err < 5e-2,
        "MPR/EPA depth parity exceeded 5%: {max_err}"
    );
}

#[test]
fn mpr_depth_parity_full_1e5_env_gated() {
    if std::env::var("MPR_FULL_PARITY").is_err() {
        eprintln!("skipping 1e5 MPR/EPA parity sweep; set MPR_FULL_PARITY=1 to run");
        return;
    }
    let (compared, _inter, max_err) = run_parity(100_000, 0x5EED);
    assert!(compared > 10_000, "too few overlapping samples: {compared}");
    assert!(
        max_err < 1e-2,
        "MPR/EPA depth parity exceeded 1%: {max_err}"
    );
}

// ---------------------------------------------------------------------------
// TEST 4 — rotated box parity (looser)
// ---------------------------------------------------------------------------

#[test]
fn mpr_rotated_box_parity_smoke() {
    let mut rng = Lcg::new(0x00B0_07ED);
    let mut compared = 0usize;
    let mut max_rel_err = 0.0f64;

    for i in 0..50 {
        let ha = Vec3::new(
            rng.range(0.5, 2.0),
            rng.range(0.5, 2.0),
            rng.range(0.5, 2.0),
        );
        let hb = Vec3::new(
            rng.range(0.5, 2.0),
            rng.range(0.5, 2.0),
            rng.range(0.5, 2.0),
        );
        let a = BoxShape::new(ha);
        let b = BoxShape::new(hb);

        let ta = Transform::from_position(Vec3::zeros());
        let pos_b = Vec3::new(
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        );

        let axis = Vec3::new(
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        );
        let axis_unit = if axis.norm() < 1e-6 {
            Vec3::new(0.0, 0.0, 1.0)
        } else {
            axis / axis.norm()
        };
        let angle = rng.range(0.0, 0.4);
        let rot = quat_from_axis_angle(&axis_unit, angle);
        let tb = Transform::new(pos_b, rot);

        let e = epa_depth(&a, &ta, &b, &tb);
        let m = mpr_depth(&a, &ta, &b, &tb);
        if let (Some(ev), Some(mv)) = (e, m) {
            if ev > 1e-3 {
                compared += 1;
                let rel = (mv - ev).abs() / ev;
                if rel > max_rel_err {
                    max_rel_err = rel;
                }
                assert!(
                    rel < 0.15,
                    "rotated box parity {rel} exceeded 15% at i={i}: epa={ev} mpr={mv}"
                );
            }
        }
    }

    eprintln!("mpr_rotated_box_parity_smoke: compared={compared} max_rel_err={max_rel_err}");
    assert!(compared > 0, "no rotated box pairs were comparable");
}

// TEMP-DIAG-START
#[test]
fn temp_diag_mpr_box_detail() {
    // MPR box overlap detail: what does mpr_full actually return?
    let a = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let b = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));
    let ta = Transform::from_position(Vec3::zeros());
    for off in [1.0_f64, 1.5, 1.9] {
        let tb = Transform::from_position(Vec3::new(off, 0.0, 0.0));
        match mpr_full(&a, &ta, &b, &tb) {
            MprResult::Intersecting {
                normal,
                depth,
                point,
            } => {
                eprintln!(
                    "MPR box off={off}: Intersecting depth={depth} normal=[{:.4},{:.4},{:.4}] point=[{:.4},{:.4},{:.4}]",
                    normal.x, normal.y, normal.z, point.x, point.y, point.z
                );
            }
            MprResult::Separated => eprintln!("MPR box off={off}: Separated"),
        }
        match epa_depth(&a, &ta, &b, &tb) {
            Some(d) => eprintln!("EPA box off={off}: depth={d}"),
            None => eprintln!("EPA box off={off}: None"),
        }
    }
    // also a sphere control to confirm MPR works there
    let s = Sphere::new(1.0);
    let ts = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
    if let MprResult::Intersecting { depth, .. } = mpr_full(&s, &ta, &s, &ts) {
        eprintln!("MPR sphere control off=1.0: depth={depth}");
    }
}

#[test]
fn temp_diag_find_hang_config() {
    // Replicate run_parity's RNG stream for seed 0xC0FFEE and, for each iter,
    // print the config BEFORE calling EPA so the last line printed is the hang.
    // We only call MPR (fast) for the first 60 iters to find the box configs,
    // then print every box/box config so we can identify the hang candidate.
    // IMPORTANT: do NOT call epa_depth here (it hangs). Just dump configs.
    let mut rng = Lcg::new(0x00C0_FFEE);
    for it in 0..200usize {
        let kind = (rng.next_u64() % 3) as u8;
        let pos_b = Vec3::new(
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        );
        match kind {
            0 => {
                let _ra = rng.range(0.5, 2.0);
                let _rb = rng.range(0.5, 2.0);
            }
            1 => {
                let _ra = rng.range(0.5, 2.0);
                let hb = Vec3::new(
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                );
                eprintln!(
                    "iter={it} kind=sphere_box pos_b=[{:.4},{:.4},{:.4}] hb=[{:.4},{:.4},{:.4}]",
                    pos_b.x, pos_b.y, pos_b.z, hb.x, hb.y, hb.z
                );
            }
            _ => {
                let ha = Vec3::new(
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                );
                let hb = Vec3::new(
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                );
                eprintln!(
                    "iter={it} kind=box_box pos_b=[{:.4},{:.4},{:.4}] ha=[{:.4},{:.4},{:.4}] hb=[{:.4},{:.4},{:.4}]",
                    pos_b.x, pos_b.y, pos_b.z, ha.x, ha.y, ha.z, hb.x, hb.y, hb.z
                );
            }
        }
    }
}

#[test]
fn temp_diag_bisect_hang() {
    // Run epa_depth on the FIRST N box/box and sphere/box configs from seed
    // 0xC0FFEE, one at a time with eprintln BEFORE each call, so the LAST
    // printed "calling EPA iter=K" with no following "ok iter=K" is the hang.
    // We replicate the exact stream and call epa_depth only for box-involving
    // kinds (the suspected hang). Bounded to first 120 iters.
    let mut rng = Lcg::new(0x00C0_FFEE);
    for it in 0..400usize {
        let kind = (rng.next_u64() % 3) as u8;
        let ta = Transform::from_position(Vec3::zeros());
        let pos_b = Vec3::new(
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        );
        match kind {
            0 => {
                let a = Sphere::new(rng.range(0.5, 2.0));
                let b = Sphere::new(rng.range(0.5, 2.0));
                let tb = Transform::from_position(pos_b);
                let _ = epa_depth(&a, &ta, &b, &tb);
            }
            1 => {
                let a = Sphere::new(rng.range(0.5, 2.0));
                let hb = Vec3::new(
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                );
                let b = BoxShape::new(hb);
                let tb = Transform::from_position(pos_b);
                eprintln!("calling EPA iter={it} kind=sphere_box");
                let d = epa_depth(&a, &ta, &b, &tb);
                eprintln!("ok iter={it} d={d:?}");
            }
            _ => {
                let ha = Vec3::new(
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                );
                let hb = Vec3::new(
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                    rng.range(0.5, 2.0),
                );
                let a = BoxShape::new(ha);
                let b = BoxShape::new(hb);
                let tb = Transform::from_position(pos_b);
                eprintln!(
                    "calling EPA iter={it} kind=box_box pos_b=[{:.4},{:.4},{:.4}] ha=[{:.4},{:.4},{:.4}] hb=[{:.4},{:.4},{:.4}]",
                    pos_b.x, pos_b.y, pos_b.z, ha.x, ha.y, ha.z, hb.x, hb.y, hb.z
                );
                let d = epa_depth(&a, &ta, &b, &tb);
                eprintln!("ok iter={it} d={d:?}");
            }
        }
    }
}
// TEMP-DIAG-END
