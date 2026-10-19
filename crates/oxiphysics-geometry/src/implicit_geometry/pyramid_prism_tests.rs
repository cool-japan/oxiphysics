// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the exact signed distance functions of [`SdfPyramid`]
//! and [`SdfTriangularPrism`].
//!
//! Both solids are convex polytopes, so the tests compare against an independent
//! exact reference: inside, the largest signed face-plane distance; outside, the
//! minimum distance to the triangulated boundary (closest-point-on-triangle).

use super::functions::Sdf;
use super::types::{SdfPyramid, SdfTriangularPrism};

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// Closest point of the triangle `abc` to `p` (Ericson, *Real-Time Collision
/// Detection*, section 5.1.5).
fn closest_point_on_triangle(p: V3, a: V3, b: V3, c: V3) -> V3 {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return add(a, scale(ab, d1 / (d1 - d3)));
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return add(a, scale(ac, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add(b, scale(sub(c, b), w));
    }
    let denom = 1.0 / (va + vb + vc);
    add(a, add(scale(ab, vb * denom), scale(ac, vc * denom)))
}

/// A convex polytope given by its triangulated boundary.
struct ConvexPolytope {
    triangles: Vec<[V3; 3]>,
    /// Outward unit normal and plane offset `n . x = offset` of every triangle.
    planes: Vec<(V3, f64)>,
}

impl ConvexPolytope {
    fn new(triangles: Vec<[V3; 3]>) -> Self {
        let corners: Vec<V3> = triangles.iter().flatten().copied().collect();
        let centroid = scale(
            corners.iter().fold([0.0; 3], |acc, &c| add(acc, c)),
            1.0 / corners.len() as f64,
        );
        let planes = triangles
            .iter()
            .map(|t| {
                let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
                let unit = scale(n, 1.0 / norm(n));
                let outward = if dot(unit, sub(t[0], centroid)) < 0.0 {
                    scale(unit, -1.0)
                } else {
                    unit
                };
                (outward, dot(outward, t[0]))
            })
            .collect();
        Self { triangles, planes }
    }

    /// Exact signed distance: the largest (negative) plane distance inside, the
    /// distance to the closest boundary triangle outside.
    fn signed_distance(&self, p: V3) -> f64 {
        let deepest = self
            .planes
            .iter()
            .map(|&(n, offset)| dot(n, p) - offset)
            .fold(f64::NEG_INFINITY, f64::max);
        if deepest <= 0.0 {
            return deepest;
        }
        self.triangles
            .iter()
            .map(|t| norm(sub(p, closest_point_on_triangle(p, t[0], t[1], t[2]))))
            .fold(f64::INFINITY, f64::min)
    }
}

/// Reference polytope of the pyramid: base square on `y = 0`, apex `(0, h, 0)`.
fn pyramid_reference(half_base: f64, height: f64) -> ConvexPolytope {
    let c = [
        [-half_base, 0.0, -half_base],
        [half_base, 0.0, -half_base],
        [half_base, 0.0, half_base],
        [-half_base, 0.0, half_base],
    ];
    let apex = [0.0, height, 0.0];
    let mut triangles = vec![[c[0], c[1], c[2]], [c[0], c[2], c[3]]];
    for i in 0..4 {
        triangles.push([c[i], c[(i + 1) % 4], apex]);
    }
    ConvexPolytope::new(triangles)
}

/// Reference polytope of the prism: equilateral triangle in `(x, z)` with its
/// apex towards `+z`, extruded along `y`.
fn prism_reference(side: f64, half_height: f64) -> ConvexPolytope {
    let r = 0.5 * side;
    let k = 3.0_f64.sqrt();
    let corners_xz = [[0.0, 2.0 * r / k], [-r, -r / k], [r, -r / k]];
    let at = |i: usize, y: f64| [corners_xz[i][0], y, corners_xz[i][1]];
    let mut triangles = vec![
        [at(0, half_height), at(1, half_height), at(2, half_height)],
        [
            at(0, -half_height),
            at(1, -half_height),
            at(2, -half_height),
        ],
    ];
    for i in 0..3 {
        let j = (i + 1) % 3;
        triangles.push([at(i, -half_height), at(j, -half_height), at(j, half_height)]);
        triangles.push([at(i, -half_height), at(j, half_height), at(i, half_height)]);
    }
    ConvexPolytope::new(triangles)
}

/// Deterministic pseudo-random sequence (64-bit LCG, Knuth constants).
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn point_in(&mut self, lo: V3, hi: V3) -> V3 {
        [
            lo[0] + (hi[0] - lo[0]) * self.next_unit(),
            lo[1] + (hi[1] - lo[1]) * self.next_unit(),
            lo[2] + (hi[2] - lo[2]) * self.next_unit(),
        ]
    }
}

fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{what}: got {actual}, expected {expected}"
    );
}

// ── SdfPyramid ───────────────────────────────────────────────────────────────

#[test]
fn pyramid_reference_points_from_the_audit() {
    let pyramid = SdfPyramid::new(1.0, 1.0);
    // Above the apex: outside, at the distance to the apex.
    assert_close(pyramid.dist([0.0, 2.0, 0.0]), 1.0, 1e-12, "above the apex");
    // An interior point: nearest to a lateral face, -sqrt(2) / 8.
    assert_close(
        pyramid.dist([0.5, 0.25, 0.5]),
        -(2.0_f64.sqrt()) / 8.0,
        1e-12,
        "interior point",
    );
    // Straight below the base centre.
    assert_close(pyramid.dist([0.0, -1.0, 0.0]), 1.0, 1e-12, "below the base");
}

#[test]
fn pyramid_interior_near_the_base_measures_the_base_distance() {
    let pyramid = SdfPyramid::new(1.0, 1.0);
    // The lateral faces are 0.67 away, the base only 0.05.
    assert_close(
        pyramid.dist([0.0, 0.05, 0.0]),
        -0.05,
        1e-12,
        "near base centre",
    );
    assert_close(pyramid.dist([0.3, 0.02, -0.2]), -0.02, 1e-12, "off-centre");
    // Higher up, the lateral faces (x + y = 1, distance (1 - y) / sqrt 2) win.
    assert_close(
        pyramid.dist([0.0, 0.5, 0.0]),
        -0.5 / 2.0_f64.sqrt(),
        1e-12,
        "lateral face is the nearest boundary",
    );
}

#[test]
fn pyramid_surface_points_read_zero() {
    let pyramid = SdfPyramid::new(1.0, 1.0);
    for p in [
        [0.0, 1.0, 0.0],  // apex
        [0.0, 0.0, 0.0],  // base centre
        [1.0, 0.0, 1.0],  // base corner
        [1.0, 0.0, -0.4], // base edge
        [0.5, 0.5, 0.0],  // midpoint of the +x face's central line
        [-0.25, 0.75, 0.0],
        [0.0, 0.5, -0.5],
    ] {
        assert_close(pyramid.dist(p), 0.0, 1e-12, &format!("surface point {p:?}"));
    }
}

#[test]
fn pyramid_matches_the_exact_polytope_distance() {
    let mut rng = Lcg(0x5eed_1234_abcd_0001);
    for &(half_base, height) in &[(1.0, 1.0), (0.3, 2.0), (2.0, 0.5), (0.75, 1.9)] {
        let pyramid = SdfPyramid::new(half_base, height);
        let reference = pyramid_reference(half_base, height);
        let extent = 2.0 * half_base.max(height);
        let (lo, hi) = (
            [-extent, -0.8 * extent, -extent],
            [extent, 1.2 * extent, extent],
        );
        for _ in 0..3000 {
            let p = rng.point_in(lo, hi);
            let expected = reference.signed_distance(p);
            let actual = pyramid.dist(p);
            assert!(
                (actual - expected).abs() <= 1e-9 * (1.0 + expected.abs()),
                "hb={half_base} h={height} p={p:?}: got {actual}, exact {expected}"
            );
        }
    }
}

#[test]
fn pyramid_is_symmetric_and_continuous_across_the_base_plane() {
    let pyramid = SdfPyramid::new(0.8, 1.3);
    let mut rng = Lcg(77);
    for _ in 0..500 {
        let p = rng.point_in([-2.0, 0.0, -2.0], [2.0, 2.0, 2.0]);
        let reference = pyramid.dist(p);
        for q in [
            [-p[0], p[1], p[2]],
            [p[0], p[1], -p[2]],
            [p[2], p[1], p[0]],
            [-p[2], p[1], -p[0]],
        ] {
            assert_close(pyramid.dist(q), reference, 1e-12, "symmetry");
        }
        let above = pyramid.dist([p[0], 1e-9, p[2]]);
        let below = pyramid.dist([p[0], -1e-9, p[2]]);
        assert!(
            (above - below).abs() < 5e-9,
            "discontinuous across y = 0 at {p:?}: {above} vs {below}"
        );
    }
}

#[test]
fn pyramid_distance_is_one_lipschitz() {
    let pyramid = SdfPyramid::new(1.0, 1.5);
    let mut rng = Lcg(2024);
    for _ in 0..2000 {
        let a = rng.point_in([-3.0, -2.0, -3.0], [3.0, 3.0, 3.0]);
        let b = add(a, rng.point_in([-0.2, -0.2, -0.2], [0.2, 0.2, 0.2]));
        let gap = (pyramid.dist(a) - pyramid.dist(b)).abs();
        assert!(
            gap <= norm(sub(a, b)) * (1.0 + 1e-9),
            "{a:?} -> {b:?}: {gap}"
        );
    }
}

#[test]
fn degenerate_pyramids_collapse_to_a_segment_a_square_or_a_point() {
    // Zero base: the segment from the origin to the apex.
    let needle = SdfPyramid::new(0.0, 2.0);
    assert_close(
        needle.dist([3.0, 1.0, 4.0]),
        5.0,
        1e-12,
        "beside the needle",
    );
    assert_close(needle.dist([0.0, 5.0, 0.0]), 3.0, 1e-12, "above the needle");
    assert_close(
        needle.dist([0.0, -2.0, 0.0]),
        2.0,
        1e-12,
        "below the needle",
    );
    // Zero height: the flat base square.
    let sheet = SdfPyramid::new(1.0, 0.0);
    assert_close(sheet.dist([0.5, 3.0, -0.5]), 3.0, 1e-12, "above the sheet");
    assert_close(sheet.dist([4.0, 0.0, 0.0]), 3.0, 1e-12, "beside the sheet");
    // Both zero: a point.
    let point = SdfPyramid::new(0.0, 0.0);
    assert_close(point.dist([1.0, 2.0, 2.0]), 3.0, 1e-12, "point");
    // Negative extents behave like zero extents and never produce NaN.
    assert!(
        SdfPyramid::new(-1.0, -1.0)
            .dist([1.0, 0.0, 0.0])
            .is_finite()
    );
}

// ── SdfTriangularPrism ───────────────────────────────────────────────────────

#[test]
fn prism_reference_points_from_the_audit() {
    let (side, half_height) = (1.0, 1.0);
    let prism = SdfTriangularPrism::new(side, half_height);
    let inradius = side / (2.0 * 3.0_f64.sqrt());
    assert_close(
        prism.dist([0.0, 0.0, 0.0]),
        -inradius.min(half_height),
        1e-12,
        "centre",
    );
    assert_close(prism.dist([0.0, 5.0, 0.0]), 4.0, 1e-12, "above the prism");
    assert_close(prism.dist([0.0, -5.0, 0.0]), 4.0, 1e-12, "below the prism");
    assert!(prism.dist([5.0, 0.0, 0.0]) > 0.0, "beside the prism");
    // Nearest boundary point is the base corner (0.5, -inradius) of the section.
    assert_close(
        prism.dist([5.0, 0.0, 0.0]),
        (4.5_f64.powi(2) + inradius.powi(2)).sqrt(),
        1e-12,
        "beside the prism, corner distance",
    );
}

#[test]
fn prism_known_exterior_and_interior_distances() {
    let prism = SdfTriangularPrism::new(1.0, 1.0);
    let inradius = 1.0 / (2.0 * 3.0_f64.sqrt());
    // Straight out of the base face (the face at z = -inradius).
    assert_close(
        prism.dist([0.0, 0.0, -inradius - 1.0]),
        1.0,
        1e-12,
        "below the base face",
    );
    // Straight out of the slanted face on the +x side: outward normal (cos 30, sin 30).
    let (c30, s30) = (3.0_f64.sqrt() / 2.0, 0.5);
    let midpoint = [0.25, 0.0, 0.5 * (-inradius + 2.0 * inradius)];
    let outside = [midpoint[0] + 1.5 * c30, 0.3, midpoint[2] + 1.5 * s30];
    assert_close(prism.dist(outside), 1.5, 1e-12, "slanted face");
    // Interior points close to a cap: the cap is the nearest boundary.
    assert_close(prism.dist([0.0, 0.9, 0.0]), -0.1, 1e-12, "near the top cap");
    assert_close(
        prism.dist([0.0, -0.95, 0.05]),
        -0.05,
        1e-12,
        "near the bottom",
    );
    // The apex points towards +z.
    assert!(prism.dist([0.0, 0.0, 0.5]) < 0.0, "inside near the apex");
    assert_close(
        prism.dist([0.0, 0.0, -0.5]),
        0.5 - inradius,
        1e-12,
        "outside beyond the base edge",
    );
    assert!(prism.dist([0.0, 0.0, 0.7]) > 0.0, "beyond the apex");
}

#[test]
fn prism_surface_points_read_zero() {
    let prism = SdfTriangularPrism::new(2.0, 0.75);
    // Section corners: base (+-1, -1/sqrt 3) and apex (0, 2/sqrt 3).
    let base_z = -1.0 / 3.0_f64.sqrt();
    let apex_z = 2.0 / 3.0_f64.sqrt();
    for p in [
        [0.0, 0.0, apex_z],                   // apex edge
        [1.0, 0.3, base_z],                   // base corner edge
        [0.0, 0.75, 0.0],                     // top cap
        [0.2, -0.75, 0.1],                    // bottom cap
        [0.0, 0.1, base_z],                   // middle of the base face
        [0.5, -0.4, 0.5 * (base_z + apex_z)], // middle of the +x slanted face
        [-0.5, 0.6, 0.5 * (base_z + apex_z)], // middle of the -x slanted face
    ] {
        assert_close(prism.dist(p), 0.0, 1e-12, &format!("surface point {p:?}"));
    }
}

#[test]
fn prism_matches_the_exact_polytope_distance() {
    let mut rng = Lcg(0xfeed_beef_0000_0042);
    for &(side, half_height) in &[(1.0, 1.0), (2.0, 0.3), (0.4, 2.5), (3.0, 1.2)] {
        let prism = SdfTriangularPrism::new(side, half_height);
        let reference = prism_reference(side, half_height);
        let extent = 1.5 * side.max(half_height);
        let (lo, hi) = ([-extent; 3], [extent; 3]);
        for _ in 0..3000 {
            let p = rng.point_in(lo, hi);
            let expected = reference.signed_distance(p);
            let actual = prism.dist(p);
            assert!(
                (actual - expected).abs() <= 1e-9 * (1.0 + expected.abs()),
                "side={side} hh={half_height} p={p:?}: got {actual}, exact {expected}"
            );
        }
    }
}

#[test]
fn prism_is_mirror_symmetric_and_one_lipschitz() {
    let prism = SdfTriangularPrism::new(1.4, 0.9);
    let mut rng = Lcg(31337);
    for _ in 0..2000 {
        let a = rng.point_in([-2.0; 3], [2.0; 3]);
        let mirrored_x = [-a[0], a[1], a[2]];
        let mirrored_y = [a[0], -a[1], a[2]];
        assert_close(prism.dist(mirrored_x), prism.dist(a), 1e-12, "x mirror");
        assert_close(prism.dist(mirrored_y), prism.dist(a), 1e-12, "y mirror");
        let b = add(a, rng.point_in([-0.2; 3], [0.2; 3]));
        let gap = (prism.dist(a) - prism.dist(b)).abs();
        assert!(
            gap <= norm(sub(a, b)) * (1.0 + 1e-9),
            "{a:?} -> {b:?}: {gap}"
        );
    }
}

#[test]
fn degenerate_prisms_collapse_to_an_axis_segment_or_a_triangle() {
    // Zero side: the segment on the y axis.
    let needle = SdfTriangularPrism::new(0.0, 1.0);
    assert_close(
        needle.dist([3.0, 0.0, 4.0]),
        5.0,
        1e-12,
        "beside the segment",
    );
    assert_close(
        needle.dist([0.0, 3.0, 0.0]),
        2.0,
        1e-12,
        "above the segment",
    );
    // Zero height: the flat triangle on y = 0.
    let flat = SdfTriangularPrism::new(1.0, 0.0);
    assert_close(flat.dist([0.0, 2.0, 0.0]), 2.0, 1e-12, "above the triangle");
    assert_close(
        flat.dist([0.0, 0.0, 0.0]),
        0.0,
        1e-12,
        "inside the flat triangle",
    );
    // Negative extents never produce NaN.
    assert!(
        SdfTriangularPrism::new(-1.0, -1.0)
            .dist([1.0, 1.0, 1.0])
            .is_finite()
    );
}
