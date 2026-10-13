// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact-value regression tests for [`gjk_distance`], [`gjk_distance_with_witnesses`]
//! and [`closest_points_between_shapes`].
//!
//! The previous implementation returned the norm of the first Minkowski support
//! vertex, which is only correct when the shapes happen to be aligned with the
//! initial search direction. Every expectation below is a closed-form distance.

use super::*;

const TOL: f64 = 1e-6;

fn sphere(cx: f64, cy: f64, cz: f64, r: f64) -> GjkSphere {
    GjkSphere {
        center: [cx, cy, cz],
        radius: r,
    }
}

fn cuboid(c: [f64; 3], h: [f64; 3]) -> GjkBox {
    GjkBox {
        center: c,
        half_extents: h,
    }
}

fn distance_between(a: [f64; 3], b: [f64; 3]) -> f64 {
    len3(sub3(a, b))
}

// ── the examples reported by the audit ───────────────────────────────────────

#[test]
fn spheres_separated_along_y_report_the_surface_gap() {
    // F0123: unit spheres at (0,0,0) and (0,4,0) returned ~4.47 instead of 2.
    let d = gjk_distance(&sphere(0.0, 0.0, 0.0, 1.0), &sphere(0.0, 4.0, 0.0, 1.0));
    let d = d.expect("disjoint spheres");
    assert!((d - 2.0).abs() < TOL, "distance {d}");
}

#[test]
fn spheres_separated_along_y_by_five() {
    // F0124: unit spheres at (0,0,0) and (0,5,0) returned sqrt(29) ~ 5.385 instead of 3.
    let d = gjk_distance(&sphere(0.0, 0.0, 0.0, 1.0), &sphere(0.0, 5.0, 0.0, 1.0));
    let d = d.expect("disjoint spheres");
    assert!((d - 3.0).abs() < TOL, "distance {d}");
}

#[test]
fn unit_boxes_offset_in_y_report_the_x_gap() {
    // Boxes of half extent 1 at the origin and (4, 0.5, 0): gap along x is 2.
    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([4.0, 0.5, 0.0], [1.0; 3]);
    let d = gjk_distance(&a, &b).expect("disjoint boxes");
    assert!((d - 2.0).abs() < TOL, "distance {d}");
}

#[test]
fn box_witnesses_lie_on_the_facing_faces() {
    // F0124: box half 1 at the origin vs a box at (3, 0.5, 0): distance 1,
    // witness on A at x = 1 (the old vertex walk reported 2.69).
    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([3.0, 0.5, 0.0], [1.0; 3]);
    let (d, pa, pb) = gjk_distance_with_witnesses(&a, &b).expect("disjoint boxes");
    assert!((d - 1.0).abs() < TOL, "distance {d}");
    assert!((pa[0] - 1.0).abs() < 1e-9, "witness_a {pa:?}");
    assert!((pb[0] - 2.0).abs() < 1e-9, "witness_b {pb:?}");
    assert!((distance_between(pa, pb) - d).abs() < 1e-9);
    // Both witnesses lie inside their boxes.
    for k in 0..3 {
        assert!(pa[k].abs() <= 1.0 + 1e-9, "witness_a {pa:?}");
        assert!(
            (pb[k] - b.center[k]).abs() <= 1.0 + 1e-9,
            "witness_b {pb:?}"
        );
    }
}

#[test]
fn sphere_witnesses_are_the_facing_surface_points() {
    let a = sphere(0.0, 0.0, 0.0, 1.0);
    let b = sphere(0.0, 5.0, 0.0, 1.0);
    let (d, pa, pb) = gjk_distance_with_witnesses(&a, &b).expect("disjoint spheres");
    assert!((d - 3.0).abs() < TOL);
    // Witness positions on smooth shapes are accurate to ~sqrt(tolerance).
    assert!(
        distance_between(pa, [0.0, 1.0, 0.0]) < 1e-5,
        "witness_a {pa:?}"
    );
    assert!(
        distance_between(pb, [0.0, 4.0, 0.0]) < 1e-5,
        "witness_b {pb:?}"
    );
    // Both witnesses are exact surface points and realise the distance.
    assert!((distance_between(pa, a.center) - 1.0).abs() < 1e-9);
    assert!((distance_between(pb, b.center) - 1.0).abs() < 1e-9);
    assert!((distance_between(pa, pb) - d).abs() < 1e-9);
}

#[test]
fn closest_points_between_boxes_are_exact() {
    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([5.0, 0.0, 0.0], [1.0; 3]);
    let cp = closest_points_between_shapes(&a, &b).expect("disjoint boxes");
    assert!((cp.distance - 3.0).abs() < TOL, "distance {}", cp.distance);
    assert!((cp.point_a[0] - 1.0).abs() < 1e-9, "{:?}", cp.point_a);
    assert!((cp.point_b[0] - 4.0).abs() < 1e-9, "{:?}", cp.point_b);
    assert!((distance_between(cp.point_a, cp.point_b) - cp.distance).abs() < 1e-9);
}

#[test]
fn axis_aligned_pairs_keep_their_exact_distance() {
    let d = gjk_distance(&sphere(0.0, 0.0, 0.0, 1.0), &sphere(4.0, 0.0, 0.0, 1.0));
    assert!((d.expect("disjoint") - 2.0).abs() < TOL);

    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([4.0, 0.0, 0.0], [1.0; 3]);
    assert!((gjk_distance(&a, &b).expect("disjoint") - 2.0).abs() < TOL);

    // Sphere radius 0.5 at the origin vs a half-1 box at x = 4: gap 4 - 0.5 - 1.
    let s = sphere(0.0, 0.0, 0.0, 0.5);
    let d = gjk_distance(&s, &b).expect("disjoint");
    assert!((d - 2.5).abs() < TOL, "distance {d}");
}

// ── overlap and touching ─────────────────────────────────────────────────────

#[test]
fn overlapping_shapes_report_none() {
    let a = cuboid([0.0; 3], [1.0; 3]);
    assert!(gjk_distance(&a, &cuboid([1.5, 0.5, -0.2], [1.0; 3])).is_none());
    assert!(
        gjk_distance(&a, &cuboid([0.0; 3], [1.0; 3])).is_none(),
        "concentric"
    );
    assert!(
        gjk_distance(&a, &cuboid([0.0; 3], [0.25; 3])).is_none(),
        "contained"
    );
    assert!(gjk_distance(&sphere(0.0, 0.0, 0.0, 1.0), &sphere(0.5, 0.0, 0.0, 1.0)).is_none());
    assert!(
        gjk_distance(&sphere(0.1, 0.2, 0.0, 0.3), &a).is_none(),
        "sphere in box"
    );
    assert!(gjk_distance_with_witnesses(&a, &cuboid([0.9, 0.9, 0.9], [1.0; 3])).is_none());
    assert!(closest_points_between_shapes(&a, &cuboid([0.5, 0.0, 0.0], [1.0; 3])).is_none());
}

#[test]
fn face_touching_boxes_count_as_overlapping() {
    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([2.0, 0.0, 0.0], [1.0; 3]);
    assert!(gjk_distance(&a, &b).is_none());
}

#[test]
fn touching_spheres_are_either_overlapping_or_at_zero_distance() {
    let a = sphere(0.0, 0.0, 0.0, 1.0);
    let b = sphere(0.0, 2.0, 0.0, 1.0);
    match gjk_distance(&a, &b) {
        None => {}
        Some(d) => assert!(d < 1e-7, "tangent spheres reported distance {d}"),
    }
}

#[test]
fn point_shapes_report_the_point_distance_and_the_points() {
    let a = sphere(1.0, 2.0, 3.0, 0.0);
    let b = sphere(4.0, 6.0, 3.0, 0.0);
    let (d, pa, pb) = gjk_distance_with_witnesses(&a, &b).expect("disjoint points");
    assert!((d - 5.0).abs() < 1e-12);
    assert!(distance_between(pa, [1.0, 2.0, 3.0]) < 1e-12);
    assert!(distance_between(pb, [4.0, 6.0, 3.0]) < 1e-12);
}

// ── other shape types ────────────────────────────────────────────────────────

#[test]
fn capsule_and_sphere_distance() {
    // Capsule axis from (0,-1,0) to (0,1,0), r = 0.5; sphere r = 0.25 at (3, 4, 0).
    // Nearest axis point is (0,1,0): distance sqrt(9 + 9) - 0.75.
    let cap = GjkCapsule {
        p1: [0.0, -1.0, 0.0],
        p2: [0.0, 1.0, 0.0],
        radius: 0.5,
    };
    let s = sphere(3.0, 4.0, 0.0, 0.25);
    let expected = 18.0_f64.sqrt() - 0.75;
    let d = gjk_distance(&cap, &s).expect("disjoint");
    assert!((d - expected).abs() < TOL, "distance {d} vs {expected}");
}

#[test]
fn cylinder_and_sphere_distance() {
    let cyl = GjkCylinder {
        base: [0.0, 0.0, 0.0],
        axis: [0.0, 1.0, 0.0],
        height: 2.0,
        radius: 1.0,
    };
    let s = sphere(3.0, 1.0, 0.0, 0.5);
    let d = gjk_distance(&cyl, &s).expect("disjoint");
    assert!((d - 1.5).abs() < TOL, "distance {d}");
}

#[test]
fn convex_hull_and_sphere_distance() {
    let hull = GjkConvexHull::new(vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]);
    let s = sphere(4.0, 0.0, 0.0, 1.0);
    let d = gjk_distance(&hull, &s).expect("disjoint");
    assert!((d - 2.0).abs() < TOL, "distance {d}");
}

#[test]
fn rotated_box_distance_uses_the_rotated_vertex() {
    // Unit-half box rotated 45 degrees about z, centre (3,0,0): its nearest
    // vertex sits at x = 3 - sqrt(2); the fixed box has its face at x = 1.
    let c = std::f64::consts::FRAC_1_SQRT_2;
    let rot = [[c, -c, 0.0], [c, c, 0.0], [0.0, 0.0, 1.0]];
    let unit = cuboid([0.0; 3], [1.0; 3]);
    let moving = TransformedShape {
        shape: &unit,
        rotation: rot,
        translation: [3.0, 0.0, 0.0],
    };
    let fixed = cuboid([0.0; 3], [1.0; 3]);
    let d = gjk_distance(&fixed, &moving).expect("disjoint");
    let expected = 3.0 - 2.0_f64.sqrt() - 1.0;
    assert!((d - expected).abs() < TOL, "distance {d} vs {expected}");
}

#[test]
fn distance_is_symmetric_in_the_arguments() {
    let a = sphere(0.3, -0.2, 0.1, 0.7);
    let b = cuboid([3.0, 1.5, -0.5], [0.6, 0.9, 0.4]);
    let ab = gjk_distance(&a, &b).expect("disjoint");
    let ba = gjk_distance(&b, &a).expect("disjoint");
    assert!((ab - ba).abs() < 1e-8, "{ab} vs {ba}");
}

// ── randomised comparison with closed-form distances ─────────────────────────

/// Deterministic xorshift generator.
struct Rng(u64);

impl Rng {
    fn unit(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    fn vec(&mut self, lo: f64, hi: f64) -> [f64; 3] {
        [self.range(lo, hi), self.range(lo, hi), self.range(lo, hi)]
    }
}

/// Distance from point `p` to segment `a`-`b`.
fn point_segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let t = (dot3(sub3(p, a), ab) / dot3(ab, ab)).clamp(0.0, 1.0);
    len3(sub3(p, add3(a, scale3(ab, t))))
}

/// Closest-point distance between segments `p1`-`q1` and `p2`-`q2`.
fn segment_segment_distance(p1: [f64; 3], q1: [f64; 3], p2: [f64; 3], q2: [f64; 3]) -> f64 {
    let d1 = sub3(q1, p1);
    let d2 = sub3(q2, p2);
    let r = sub3(p1, p2);
    let a = dot3(d1, d1);
    let e = dot3(d2, d2);
    let f = dot3(d2, r);
    let c = dot3(d1, r);
    let b = dot3(d1, d2);
    let denom = a * e - b * b;
    let mut s = if denom > 1e-12 {
        ((b * f - c * e) / denom).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let mut t = (b * s + f) / e;
    if t < 0.0 {
        t = 0.0;
        s = (-c / a).clamp(0.0, 1.0);
    } else if t > 1.0 {
        t = 1.0;
        s = ((b - c) / a).clamp(0.0, 1.0);
    }
    len3(sub3(add3(p1, scale3(d1, s)), add3(p2, scale3(d2, t))))
}

/// Signed gap of two axis-aligned boxes: positive distance when disjoint,
/// negative (the least overlap) when they intersect.
fn box_box_signed_gap(ca: [f64; 3], ha: [f64; 3], cb: [f64; 3], hb: [f64; 3]) -> f64 {
    let gaps: Vec<f64> = (0..3)
        .map(|k| (ca[k] - cb[k]).abs() - (ha[k] + hb[k]))
        .collect();
    if gaps.iter().all(|g| *g < 0.0) {
        gaps.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    } else {
        gaps.iter()
            .map(|g| g.max(0.0) * g.max(0.0))
            .sum::<f64>()
            .sqrt()
    }
}

/// Check one random pair against its closed-form signed distance.
fn check_against_oracle(
    a: &dyn ConvexShape,
    b: &dyn ConvexShape,
    oracle: f64,
    label: &str,
    checked: &mut usize,
) {
    if oracle.abs() < 1e-3 {
        return; // borderline: either classification is acceptable
    }
    *checked += 1;
    let result = gjk_distance_with_witnesses(a, b);
    if oracle < 0.0 {
        assert!(
            result.is_none(),
            "{label}: expected overlap (oracle {oracle}), got {result:?}"
        );
    } else {
        let (d, pa, pb) = result.unwrap_or_else(|| panic!("{label}: expected distance {oracle}"));
        assert!(
            (d - oracle).abs() < 1e-7,
            "{label}: distance {d} vs oracle {oracle}"
        );
        assert!(
            (distance_between(pa, pb) - d).abs() < 1e-7,
            "{label}: witnesses {pa:?} {pb:?} do not realise distance {d}"
        );
    }
}

#[test]
fn matches_closed_form_distances_on_random_configurations() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    let mut checked = 0;
    for _ in 0..400 {
        // Sphere - sphere.
        let (ca, cb) = (rng.vec(-4.0, 4.0), rng.vec(-4.0, 4.0));
        let (ra, rb) = (rng.range(0.2, 1.5), rng.range(0.2, 1.5));
        let oracle = distance_between(ca, cb) - ra - rb;
        let sa = GjkSphere {
            center: ca,
            radius: ra,
        };
        let sb = GjkSphere {
            center: cb,
            radius: rb,
        };
        check_against_oracle(&sa, &sb, oracle, "sphere-sphere", &mut checked);

        // Box - box (axis aligned).
        let (ha, hb) = (rng.vec(0.2, 1.5), rng.vec(0.2, 1.5));
        let oracle = box_box_signed_gap(ca, ha, cb, hb);
        let ba = cuboid(ca, ha);
        let bb = cuboid(cb, hb);
        check_against_oracle(&ba, &bb, oracle, "box-box", &mut checked);

        // Sphere - box.
        let q = [
            (ca[0] - cb[0]).clamp(-hb[0], hb[0]) + cb[0],
            (ca[1] - cb[1]).clamp(-hb[1], hb[1]) + cb[1],
            (ca[2] - cb[2]).clamp(-hb[2], hb[2]) + cb[2],
        ];
        let inside = (0..3).all(|k| (ca[k] - cb[k]).abs() <= hb[k]);
        let oracle = if inside {
            -1.0
        } else {
            distance_between(ca, q) - ra
        };
        check_against_oracle(&sa, &bb, oracle, "sphere-box", &mut checked);

        // Sphere - capsule and capsule - capsule.
        let (p1, p2) = (rng.vec(-3.0, 3.0), rng.vec(-3.0, 3.0));
        if distance_between(p1, p2) > 0.3 {
            let rc = rng.range(0.1, 1.0);
            let cap = GjkCapsule { p1, p2, radius: rc };
            let oracle = point_segment_distance(ca, p1, p2) - ra - rc;
            check_against_oracle(&sa, &cap, oracle, "sphere-capsule", &mut checked);

            let (r1, r2) = (rng.vec(-3.0, 3.0), rng.vec(-3.0, 3.0));
            if distance_between(r1, r2) > 0.3 {
                let rd = rng.range(0.1, 1.0);
                let cap2 = GjkCapsule {
                    p1: r1,
                    p2: r2,
                    radius: rd,
                };
                let oracle = segment_segment_distance(p1, p2, r1, r2) - rc - rd;
                check_against_oracle(&cap, &cap2, oracle, "capsule-capsule", &mut checked);
            }
        }
    }
    assert!(
        checked > 1200,
        "only {checked} configurations were decisive"
    );
}

// ── exact oracle for convex polyhedra ────────────────────────────────────────

/// Distance from `p` to the triangle `abc`: the plane distance when the foot of
/// the perpendicular lies inside the triangle, otherwise the distance to the
/// nearest edge.
fn point_triangle_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let edges = point_segment_distance(p, a, b)
        .min(point_segment_distance(p, b, c))
        .min(point_segment_distance(p, c, a));
    let normal = cross3(sub3(b, a), sub3(c, a));
    let normal_sq = dot3(normal, normal);
    if normal_sq < 1e-24 {
        return edges;
    }
    let offset = dot3(sub3(p, a), normal) / normal_sq;
    let foot = sub3(p, scale3(normal, offset));
    let inside = [(a, b), (b, c), (c, a)]
        .iter()
        .all(|&(from, to)| dot3(cross3(sub3(to, from), sub3(foot, from)), normal) >= 0.0);
    if inside {
        offset.abs() * normal_sq.sqrt()
    } else {
        edges
    }
}

/// Exact distance between the convex hulls of two disjoint point sets.
///
/// The closest pair of two disjoint polyhedra is a vertex and a face or two
/// edges. Every face is covered by the triangles of its vertices and every edge
/// is a vertex pair, so the minimum over all vertex/triangle and edge/edge
/// combinations (chords through the interior can only be farther away) is exact.
fn polyhedron_distance(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let mut best = f64::INFINITY;
    for (own, other) in [(a, b), (b, a)] {
        for &p in own {
            for i in 0..other.len() {
                for j in (i + 1)..other.len() {
                    for k in (j + 1)..other.len() {
                        best = best.min(point_triangle_distance(p, other[i], other[j], other[k]));
                    }
                }
            }
        }
    }
    for i in 0..a.len() {
        for j in (i + 1)..a.len() {
            for k in 0..b.len() {
                for l in (k + 1)..b.len() {
                    best = best.min(segment_segment_distance(a[i], a[j], b[k], b[l]));
                }
            }
        }
    }
    best
}

/// Separating-axis gap of two convex hulls: the largest of
/// `max(lo_b - hi_a, lo_a - hi_b)` over the face normals and edge-edge cross
/// products of the point sets. Positive exactly when the hulls are disjoint;
/// otherwise its magnitude is the penetration depth.
fn polyhedron_axis_gap(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let mut axes = Vec::new();
    for cloud in [a, b] {
        for i in 0..cloud.len() {
            for j in (i + 1)..cloud.len() {
                for k in (j + 1)..cloud.len() {
                    axes.push(cross3(sub3(cloud[j], cloud[i]), sub3(cloud[k], cloud[i])));
                }
            }
        }
    }
    for i in 0..a.len() {
        for j in (i + 1)..a.len() {
            for k in 0..b.len() {
                for l in (k + 1)..b.len() {
                    axes.push(cross3(sub3(a[j], a[i]), sub3(b[l], b[k])));
                }
            }
        }
    }
    let interval = |cloud: &[[f64; 3]], u: [f64; 3]| {
        cloud
            .iter()
            .map(|p| dot3(*p, u))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                (lo.min(x), hi.max(x))
            })
    };
    let mut best = f64::NEG_INFINITY;
    for axis in axes {
        let length = len3(axis);
        if length < 1e-9 {
            continue;
        }
        let u = scale3(axis, 1.0 / length);
        let (a_lo, a_hi) = interval(a, u);
        let (b_lo, b_hi) = interval(b, u);
        best = best.max((b_lo - a_hi).max(a_lo - b_hi));
    }
    best
}

#[test]
fn matches_the_exact_distance_between_random_convex_polyhedra() {
    let mut rng = Rng(0x7E57_AB1E_5EED_0001);
    let (mut separated, mut overlapping) = (0, 0);
    for case in 0..400 {
        let count_a = 4 + (rng.unit() * 5.0) as usize;
        let count_b = 4 + (rng.unit() * 5.0) as usize;
        let center_b = rng.vec(-1.5, 1.5);
        let cloud_a: Vec<[f64; 3]> = (0..count_a).map(|_| rng.vec(-1.0, 1.0)).collect();
        let cloud_b: Vec<[f64; 3]> = (0..count_b)
            .map(|_| add3(center_b, rng.vec(-1.0, 1.0)))
            .collect();
        let gap = polyhedron_axis_gap(&cloud_a, &cloud_b);
        if gap.abs() < 1e-3 {
            continue; // borderline: either classification is acceptable
        }
        let hull_a = GjkConvexHull::new(cloud_a.clone());
        let hull_b = GjkConvexHull::new(cloud_b.clone());
        if gap < 0.0 {
            assert!(
                gjk_distance(&hull_a, &hull_b).is_none(),
                "case {case}: hulls overlap by {} but a distance was reported",
                -gap
            );
            overlapping += 1;
            continue;
        }
        let expected = polyhedron_distance(&cloud_a, &cloud_b);
        let (d, pa, pb) = gjk_distance_with_witnesses(&hull_a, &hull_b)
            .unwrap_or_else(|| panic!("case {case}: disjoint hulls (distance {expected})"));
        assert!(
            (d - expected).abs() < 1e-7,
            "case {case}: distance {d} vs exact {expected}"
        );
        assert!(
            (distance_between(pa, pb) - d).abs() < 1e-7,
            "case {case}: witnesses do not realise the distance"
        );
        // The witness pair is optimal: each hull lies behind the plane through
        // its witness that is orthogonal to the connecting segment.
        let direction = scale3(sub3(pb, pa), 1.0 / d);
        for p in &cloud_a {
            assert!(dot3(sub3(*p, pa), direction) <= 1e-7, "case {case}: A");
        }
        for p in &cloud_b {
            assert!(dot3(sub3(*p, pb), direction) >= -1e-7, "case {case}: B");
        }
        separated += 1;
    }
    assert!(
        separated > 100 && overlapping > 60,
        "{separated} separated, {overlapping} overlapping configurations were decisive"
    );
}

#[test]
fn distance_is_accurate_at_every_scale() {
    // Every tolerance of the solver is relative to the size of the problem, so
    // the answer must not depend on the units the shapes are expressed in.
    for scale in [1e-9, 1e-3, 1.0, 1e3, 1e9] {
        let a = cuboid([0.0; 3], [scale; 3]);
        let b = cuboid([4.0 * scale, 0.5 * scale, 0.0], [scale; 3]);
        let (d, pa, pb) = gjk_distance_with_witnesses(&a, &b).expect("disjoint boxes");
        assert!(
            (d - 2.0 * scale).abs() <= 1e-9 * scale,
            "scale {scale}: distance {d}"
        );
        assert!(
            (pa[0] - scale).abs() <= 1e-9 * scale,
            "scale {scale}: witness_a {pa:?}"
        );
        assert!(
            (pb[0] - 3.0 * scale).abs() <= 1e-9 * scale,
            "scale {scale}: witness_b {pb:?}"
        );
        let touching = cuboid([2.0 * scale, 0.0, 0.0], [scale; 3]);
        assert!(
            gjk_distance(&a, &touching).is_none(),
            "scale {scale}: face contact must count as overlap"
        );
    }
}
