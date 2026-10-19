// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the exact signed distance functions of [`SdfCone`] and
//! [`SdfHexagonalPrism`].
//!
//! Besides the closed-form spot values, the tests compare against references that
//! do not share the production formulas: the distance to the boundary of the
//! meridian profile (cone) or of the hexagon (prism) computed segment by segment,
//! and, for the cone, a brute-force search over a densely sampled boundary.

use super::functions::Sdf;
use super::types::{SdfCone, SdfHexagonalPrism};
use std::f64::consts::{FRAC_PI_6, PI};

/// Deterministic pseudo-random numbers (64-bit linear congruential generator).
struct Lcg(u64);

impl Lcg {
    /// Uniform number in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform number in `[lo, hi)`.
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// Uniform point in the axis-aligned box `[-extent, extent]^3`.
    fn point(&mut self, extent: [f64; 3]) -> [f64; 3] {
        [
            self.range(-extent[0], extent[0]),
            self.range(-extent[1], extent[1]),
            self.range(-extent[2], extent[2]),
        ]
    }
}

/// Distance from `p` to the segment `ab` in the plane.
fn segment_distance(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ap = [p[0] - a[0], p[1] - a[1]];
    let len_sq = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if len_sq > 0.0 {
        ((ap[0] * ab[0] + ap[1] * ab[1]) / len_sq).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p[0] - a[0] - t * ab[0]).hypot(p[1] - a[1] - t * ab[1])
}

/// Signed distance to the cone from the meridian profile: the triangle with the
/// corners apex `(0, 0)`, rim `(h tan(a), h)` and base centre `(0, h)` revolved
/// around the Y axis.  The axis side of the triangle is interior to the solid, so
/// only the slant and the base segment are boundary.
fn cone_reference(half_angle: f64, height: f64, p: [f64; 3]) -> f64 {
    let slope = half_angle.tan();
    let rho = p[0].hypot(p[2]);
    let q = [rho, p[1]];
    let rim = [height * slope, height];
    let boundary =
        segment_distance(q, [0.0, 0.0], rim).min(segment_distance(q, rim, [0.0, height]));
    let inside = p[1] >= 0.0 && p[1] <= height && rho <= p[1] * slope;
    if inside { -boundary } else { boundary }
}

/// Smallest distance from `p` to a dense sampling of the boundary of the cone:
/// the lateral surface and the base disc.
fn cone_sampled_boundary_distance(half_angle: f64, height: f64, p: [f64; 3]) -> f64 {
    let radius = height * half_angle.tan();
    let (angles, rings) = (720usize, 360usize);
    let mut best = f64::INFINITY;
    for i in 0..angles {
        let phi = 2.0 * PI * i as f64 / angles as f64;
        let (s, c) = phi.sin_cos();
        for j in 0..=rings {
            let t = j as f64 / rings as f64;
            // Lateral surface point at height t * h, base disc point at radius t * R.
            let lateral = [t * radius * c, t * height, t * radius * s];
            let base = [t * radius * c, height, t * radius * s];
            for q in [lateral, base] {
                let d =
                    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
                best = best.min(d);
            }
        }
    }
    best
}

fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tol,
        "{what}: got {actual}, expected {expected}"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// SdfCone
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn cone_spot_values_with_the_apex_at_the_origin_and_the_base_at_plus_height() {
    let cone = SdfCone::new(FRAC_PI_6, 1.0);
    // Above the base disc: the distance to the base plane.
    assert_close(cone.dist([0.0, 2.0, 0.0]), 1.0, 1e-12, "above the base");
    // On the axis: nearer to the slant (0.5 * sin 30 deg) than to the base (0.5).
    assert_close(
        cone.dist([0.0, 0.5, 0.0]),
        -0.25,
        1e-12,
        "inside on the axis",
    );
    // Below the apex: the distance to the apex.
    assert_close(cone.dist([0.0, -1.0, 0.0]), 1.0, 1e-12, "below the apex");
    // Just above the base disc.
    assert_close(
        cone.dist([0.0, 1.1, 0.0]),
        0.1,
        1e-12,
        "just above the base",
    );
}

#[test]
fn points_above_the_base_are_outside_the_cone() {
    // The previous formula measured the angle from the radial plane and united the
    // cone with the base slab, which made everything above the base "inside".
    let cone = SdfCone::new(FRAC_PI_6, 1.0);
    for y in [1.01, 1.1, 1.5, 2.0, 10.0] {
        let d = cone.dist([0.0, y, 0.0]);
        assert_close(d, y - 1.0, 1e-12, "axis point above the base");
        assert!(d > 0.0, "(0, {y}, 0) must be outside, got {d}");
    }
    // A wide cone (half-angle 60 degrees) is the same story off the axis.
    let wide = SdfCone::new(PI / 3.0, 2.0);
    assert!(wide.dist([1.0, 2.5, 1.0]) > 0.0);
    assert!(wide.dist([0.0, 2.0001, 0.0]) > 0.0);
}

#[test]
fn cone_boundary_points_have_zero_distance() {
    let (alpha, h) = (0.5_f64, 2.0_f64);
    let cone = SdfCone::new(alpha, h);
    let r = h * alpha.tan();
    // Apex, centre of the base, rim and a slant point in an arbitrary meridian.
    let phi = 0.7_f64;
    let slant = |t: f64| [t * r * phi.cos(), t * h, t * r * phi.sin()];
    for (name, p) in [
        ("apex", [0.0, 0.0, 0.0]),
        ("base centre", [0.0, h, 0.0]),
        ("rim", [r, h, 0.0]),
        ("slant midpoint", slant(0.5)),
        ("slant quarter", slant(0.25)),
        (
            "base disc point",
            [0.5 * r * phi.cos(), h, 0.5 * r * phi.sin()],
        ),
    ] {
        assert_close(cone.dist(p), 0.0, 1e-12, name);
    }
}

#[test]
fn cone_distance_of_slant_and_rim_neighbours() {
    let (alpha, h) = (FRAC_PI_6, 1.0);
    let cone = SdfCone::new(alpha, h);
    let r = h * alpha.tan();
    // Straight out from the slant: the unit normal of the generator in the meridian
    // plane is (cos a, -sin a); the slant midpoint is (r/2, h/2).
    let mid = [0.5 * r, 0.5 * h];
    let normal = [alpha.cos(), -alpha.sin()];
    for offset in [0.05, 0.4, 3.0] {
        let p = [
            mid[0] + offset * normal[0],
            mid[1] + offset * normal[1],
            0.0,
        ];
        assert_close(cone.dist(p), offset, 1e-12, "outside the slant");
    }
    // The same line inside the cone: the slant is the nearest boundary as long as
    // the base is farther away.
    for depth in [0.02, 0.1, 0.2] {
        let q = [mid[0] - depth * normal[0], mid[1] - depth * normal[1], 0.0];
        assert_close(cone.dist(q), -depth, 1e-12, "inside the slant");
    }
    // Radially beyond the rim, in the plane of the base: nearest to the rim.
    assert_close(cone.dist([r + 0.3, h, 0.0]), 0.3, 1e-12, "beyond the rim");
    // Diagonally above the rim: distance to the rim point.
    assert_close(
        cone.dist([r + 0.3, h + 0.4, 0.0]),
        0.5,
        1e-12,
        "above and beyond the rim",
    );
}

#[test]
fn cone_matches_the_meridian_profile_reference() {
    let mut rng = Lcg(0x5eed_c0de);
    for &(alpha, h) in &[(FRAC_PI_6, 1.0), (0.2, 3.0), (1.1, 0.7), (PI / 4.0, 2.0)] {
        let cone = SdfCone::new(alpha, h);
        let extent = 1.5 * (h + h * alpha.tan());
        for _ in 0..4000 {
            let p = rng.point([extent, extent, extent]);
            let expected = cone_reference(alpha, h, p);
            assert_close(
                cone.dist(p),
                expected,
                1e-12 * (1.0 + expected.abs()),
                &format!("alpha {alpha}, h {h}, p {p:?}"),
            );
        }
    }
}

#[test]
fn cone_matches_brute_force_over_the_sampled_boundary() {
    let (alpha, h) = (FRAC_PI_6, 1.5);
    let cone = SdfCone::new(alpha, h);
    let mut rng = Lcg(42);
    for _ in 0..12 {
        let p = rng.point([2.0, 2.5, 2.0]);
        let d = cone.dist(p);
        let sampled = cone_sampled_boundary_distance(alpha, h, p);
        // The sampling is exact to about the sample spacing (< 0.01).
        assert_close(d.abs(), sampled, 0.012, &format!("p {p:?}"));
    }
}

#[test]
fn cone_is_a_body_of_revolution_and_has_unit_gradient() {
    let cone = SdfCone::new(0.6, 1.8);
    let mut rng = Lcg(7);
    for _ in 0..200 {
        let p = rng.point([2.5, 3.0, 2.5]);
        let angle = rng.range(0.0, 2.0 * PI);
        let (s, c) = angle.sin_cos();
        let rotated = [c * p[0] + s * p[2], p[1], -s * p[0] + c * p[2]];
        assert_close(cone.dist(p), cone.dist(rotated), 1e-12, "rotation about Y");
    }
    // An exact distance field has |grad| = 1 away from the medial axis; sample
    // points well outside the solid, where the nearest feature is unique.
    let eps = 1e-6;
    let mut checked = 0;
    for _ in 0..400 {
        let p = rng.point([3.0, 3.5, 3.0]);
        if cone.dist(p) < 0.3 {
            continue;
        }
        let g = [
            cone.dist([p[0] + eps, p[1], p[2]]) - cone.dist([p[0] - eps, p[1], p[2]]),
            cone.dist([p[0], p[1] + eps, p[2]]) - cone.dist([p[0], p[1] - eps, p[2]]),
            cone.dist([p[0], p[1], p[2] + eps]) - cone.dist([p[0], p[1], p[2] - eps]),
        ];
        let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt() / (2.0 * eps);
        assert_close(norm, 1.0, 1e-5, &format!("|grad| at {p:?}"));
        checked += 1;
    }
    assert!(checked > 100, "only {checked} gradient samples");
}

#[test]
fn cone_with_degenerate_parameters_stays_finite() {
    let probes = [[0.0, 0.5, 0.0], [1.0, 0.2, -0.3], [0.0, -2.0, 0.0]];
    for cone in [
        SdfCone::new(FRAC_PI_6, 0.0),
        SdfCone::new(FRAC_PI_6, -1.0),
        SdfCone::new(0.0, 1.0),
        SdfCone::new(PI / 2.0, 1.0),
        SdfCone::new(2.0, 1.0),
        SdfCone::new(-0.3, 1.0),
        SdfCone::new(f64::NAN, 1.0),
        SdfCone::new(0.5, f64::NAN),
    ] {
        for p in probes {
            let d = cone.dist(p);
            assert!(d.is_finite(), "{cone:?} at {p:?} gave {d}");
        }
    }
    // A zero-height cone is the apex point.
    assert_close(
        SdfCone::new(FRAC_PI_6, 0.0).dist([3.0, 4.0, 0.0]),
        5.0,
        1e-12,
        "point cone",
    );
    // A zero half-angle is the segment from the apex to (0, h, 0).
    assert_close(
        SdfCone::new(0.0, 2.0).dist([1.0, 1.0, 0.0]),
        1.0,
        1e-12,
        "segment cone",
    );
}

#[test]
fn cone_without_interior_has_a_positive_distance_off_the_solid() {
    // A zero-height cone is the apex point: the distance is the distance to the
    // origin, also below the apex where the slant and base conditions both vanish.
    let point = SdfCone::new(FRAC_PI_6, 0.0);
    assert_close(point.dist([0.0, -2.0, 0.0]), 2.0, 1e-12, "below the point");
    assert_close(point.dist([0.0, 2.0, 0.0]), 2.0, 1e-12, "above the point");
    assert_close(point.dist([1.0, 0.0, 2.0]), 5.0_f64.sqrt(), 1e-12, "beside");
    assert_close(point.dist([0.0, 0.0, 0.0]), 0.0, 1e-12, "at the apex");

    // A zero half-angle is the segment from the apex to (0, h, 0): points of its
    // axis beyond the end points used to read zero.
    let segment = SdfCone::new(0.0, 2.0);
    assert_close(segment.dist([0.0, -1.0, 0.0]), 1.0, 1e-12, "below the apex");
    assert_close(segment.dist([0.0, 3.5, 0.0]), 1.5, 1e-12, "above the end");
    assert_close(segment.dist([0.0, 1.0, 0.0]), 0.0, 1e-12, "on the segment");
    assert_close(segment.dist([0.0, 0.0, 0.0]), 0.0, 1e-12, "at the apex");
    assert_close(segment.dist([0.0, 2.0, 0.0]), 0.0, 1e-12, "at the end");
    let mut rng = Lcg(31);
    for _ in 0..200 {
        let p = rng.point([2.0, 3.0, 2.0]);
        let nearest = [0.0, p[1].clamp(0.0, 2.0), 0.0];
        let expected = ((p[0] - nearest[0]).powi(2)
            + (p[1] - nearest[1]).powi(2)
            + (p[2] - nearest[2]).powi(2))
        .sqrt();
        assert_close(segment.dist(p), expected, 1e-12, &format!("p {p:?}"));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SdfHexagonalPrism
// ─────────────────────────────────────────────────────────────────────────────

/// Signed distance to the regular hexagonal prism with circumradius `radius`
/// whose vertices are `(radius cos(60 k deg), radius sin(60 k deg))` in the
/// `(x, z)` plane, computed edge by edge.
fn hex_prism_reference(radius: f64, half_height: f64, p: [f64; 3]) -> f64 {
    let corner = |k: usize| {
        let a = PI / 3.0 * (k % 6) as f64;
        [radius * a.cos(), radius * a.sin()]
    };
    let q = [p[0], p[2]];
    let mut edge_distance = f64::INFINITY;
    let mut inside = true;
    for k in 0..6 {
        let (a, b) = (corner(k), corner(k + 1));
        edge_distance = edge_distance.min(segment_distance(q, a, b));
        // The corners run counter-clockwise: an interior point is to the left of
        // every directed side.
        let side = (b[0] - a[0]) * (q[1] - a[1]) - (b[1] - a[1]) * (q[0] - a[0]);
        if side < 0.0 {
            inside = false;
        }
    }
    let planar = if inside {
        -edge_distance
    } else {
        edge_distance
    };
    let axial = p[1].abs() - half_height;
    planar.max(0.0).hypot(axial.max(0.0)) + planar.max(axial).min(0.0)
}

#[test]
fn hex_prism_spot_values_for_the_circumscribed_radius() {
    let prism = SdfHexagonalPrism::new(1.0, 1.0);
    let apothem = 3.0_f64.sqrt() / 2.0;
    // A vertex lies on the circumscribed circle...
    assert_close(prism.dist([1.0, 0.0, 0.0]), 0.0, 1e-12, "vertex");
    // ...and the flat sides are at the apothem.
    assert_close(prism.dist([0.0, 0.0, apothem]), 0.0, 1e-12, "flat side");
    assert_close(prism.dist([1.5, 0.0, 0.0]), 0.5, 1e-12, "beyond a vertex");
    assert_close(prism.dist([0.0, 0.0, 0.0]), -apothem, 1e-12, "centre");
    assert_close(prism.dist([0.0, 2.0, 0.0]), 1.0, 1e-12, "above the top");
}

#[test]
fn hex_prism_points_outside_the_hexagon_are_outside() {
    // The reflection clamp used to be [0, 2] instead of (-inf, 0]: (1.5, 0, 0)
    // read -1.0 and (0, 0, 1.5) read -0.76, both "inside".
    let prism = SdfHexagonalPrism::new(1.0, 1.0);
    let apothem = 3.0_f64.sqrt() / 2.0;
    assert!(prism.dist([1.5, 0.0, 0.0]) > 0.0);
    assert_close(
        prism.dist([0.0, 0.0, 1.5]),
        1.5 - apothem,
        1e-12,
        "beyond a flat",
    );
    assert_close(
        prism.dist([0.0, 0.0, -1.5]),
        1.5 - apothem,
        1e-12,
        "beyond a flat",
    );
    assert_close(prism.dist([-1.5, 0.0, 0.0]), 0.5, 1e-12, "beyond a vertex");
    for k in 0..6 {
        let a = PI / 3.0 * k as f64;
        // 1.3 along a vertex direction and 1.3 along a flat-normal direction.
        let vertex_dir = [a.cos(), 0.0, a.sin()];
        let flat_dir = [(a + PI / 6.0).cos(), 0.0, (a + PI / 6.0).sin()];
        assert_close(
            prism.dist([1.3 * vertex_dir[0], 0.0, 1.3 * vertex_dir[2]]),
            0.3,
            1e-12,
            "vertex direction",
        );
        assert_close(
            prism.dist([1.3 * flat_dir[0], 0.0, 1.3 * flat_dir[2]]),
            1.3 - apothem,
            1e-12,
            "flat direction",
        );
    }
}

#[test]
fn hex_prism_radius_is_the_circumradius_for_any_size() {
    let (radius, half_height) = (2.5, 0.75);
    let prism = SdfHexagonalPrism::new(radius, half_height);
    let apothem = radius * 3.0_f64.sqrt() / 2.0;
    for k in 0..6 {
        let a = PI / 3.0 * k as f64;
        // Every vertex of the hexagon is on the surface, at both caps.
        assert_close(
            prism.dist([radius * a.cos(), 0.0, radius * a.sin()]),
            0.0,
            1e-12,
            "vertex",
        );
        // The midpoint of the side between two vertices.
        let b = PI / 3.0 * (k + 1) as f64;
        let mid = [
            0.5 * radius * (a.cos() + b.cos()),
            0.0,
            0.5 * radius * (a.sin() + b.sin()),
        ];
        assert_close(prism.dist(mid), 0.0, 1e-12, "side midpoint");
        assert_close(
            (mid[0] * mid[0] + mid[2] * mid[2]).sqrt(),
            apothem,
            1e-12,
            "apothem",
        );
    }
    assert_close(
        prism.dist([0.0, half_height, 0.0]),
        0.0,
        1e-12,
        "cap centre",
    );
    assert_close(
        prism.dist([0.0, 0.0, 0.0]),
        -(apothem.min(half_height)),
        1e-12,
        "centre",
    );
    assert_close(
        prism.dist([0.0, 3.0, 0.0]),
        3.0 - half_height,
        1e-12,
        "above",
    );
    // The corner of a cap: planar distance and axial distance combine.
    assert_close(
        prism.dist([radius + 0.3, half_height + 0.4, 0.0]),
        0.5,
        1e-12,
        "cap corner",
    );
}

#[test]
fn hex_prism_matches_the_edge_by_edge_reference() {
    let mut rng = Lcg(0xfeed_f00d);
    for &(radius, half_height) in &[(1.0, 1.0), (2.5, 0.75), (0.4, 3.0)] {
        let prism = SdfHexagonalPrism::new(radius, half_height);
        let extent = 1.6 * radius.max(half_height);
        for _ in 0..5000 {
            let p = rng.point([extent, extent, extent]);
            let expected = hex_prism_reference(radius, half_height, p);
            assert_close(
                prism.dist(p),
                expected,
                1e-12 * (1.0 + expected.abs()),
                &format!("radius {radius}, half_height {half_height}, p {p:?}"),
            );
        }
    }
}

#[test]
fn hex_prism_has_the_symmetry_of_a_regular_hexagon() {
    let prism = SdfHexagonalPrism::new(1.3, 0.9);
    let mut rng = Lcg(99);
    for _ in 0..300 {
        let p = rng.point([2.0, 2.0, 2.0]);
        let d = prism.dist(p);
        for (name, q) in [
            ("mirror x", [-p[0], p[1], p[2]]),
            ("mirror z", [p[0], p[1], -p[2]]),
            ("mirror y", [p[0], -p[1], p[2]]),
        ] {
            assert_close(prism.dist(q), d, 1e-12, name);
        }
        // Rotation by 60 degrees about the Y axis.
        let (s, c) = (PI / 3.0).sin_cos();
        let rotated = [c * p[0] + s * p[2], p[1], -s * p[0] + c * p[2]];
        assert_close(prism.dist(rotated), d, 1e-12, "rotation by 60 degrees");
    }
}

#[test]
fn hex_prism_has_unit_gradient_outside() {
    let prism = SdfHexagonalPrism::new(1.0, 1.0);
    let mut rng = Lcg(2024);
    let eps = 1e-6;
    let mut checked = 0;
    for _ in 0..400 {
        let p = rng.point([3.0, 3.0, 3.0]);
        if prism.dist(p) < 0.3 {
            continue;
        }
        let g = [
            prism.dist([p[0] + eps, p[1], p[2]]) - prism.dist([p[0] - eps, p[1], p[2]]),
            prism.dist([p[0], p[1] + eps, p[2]]) - prism.dist([p[0], p[1] - eps, p[2]]),
            prism.dist([p[0], p[1], p[2] + eps]) - prism.dist([p[0], p[1], p[2] - eps]),
        ];
        let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt() / (2.0 * eps);
        assert_close(norm, 1.0, 1e-5, &format!("|grad| at {p:?}"));
        checked += 1;
    }
    assert!(checked > 100, "only {checked} gradient samples");
}

#[test]
fn hex_prism_with_degenerate_parameters_stays_finite() {
    let probes = [[0.0, 0.5, 0.0], [1.0, 0.2, -0.3], [0.0, -2.0, 0.0]];
    for prism in [
        SdfHexagonalPrism::new(0.0, 1.0),
        SdfHexagonalPrism::new(1.0, 0.0),
        SdfHexagonalPrism::new(-1.0, 1.0),
        SdfHexagonalPrism::new(f64::NAN, 1.0),
        SdfHexagonalPrism::new(1.0, f64::NAN),
    ] {
        for p in probes {
            let d = prism.dist(p);
            assert!(d.is_finite(), "{prism:?} at {p:?} gave {d}");
        }
    }
    // A zero-radius prism is the vertical segment on the axis.
    assert_close(
        SdfHexagonalPrism::new(0.0, 1.0).dist([3.0, 0.0, 4.0]),
        5.0,
        1e-12,
        "segment prism",
    );
}
