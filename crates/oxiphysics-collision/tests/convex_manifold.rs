// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for one-shot full convex-convex contact manifolds.
//!
//! Covers the Pass-3 Item B work: general convex polyhedra produce 4-point
//! face-clip manifolds (or edge-edge witnesses) in a single dispatch query,
//! matching the box-box reference behaviour. All assertions go through the
//! public [`NarrowPhaseDispatcher`] API.

use oxiphysics_collision::narrowphase::ShapeType;
use oxiphysics_collision::{CollisionPair, NarrowPhaseDispatcher};
use oxiphysics_core::Transform;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::{BoxShape, ConvexHull};

/// Vertices of an axis-aligned cube of half-extent `h` centred at `c`.
fn cube_hull(c: [f64; 3], h: f64) -> ConvexHull {
    let mut v = Vec::new();
    for sx in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            for sz in [-1.0, 1.0] {
                v.push(Vec3::new(c[0] + sx * h, c[1] + sy * h, c[2] + sz * h));
            }
        }
    }
    ConvexHull::new(v)
}

/// Vertices of a hexagonal prism (6-sided) centred at `c`, radius `r`, the prism
/// axis along y with half-height `hh`. Produces a 6-sided convex hull.
fn hex_prism(c: [f64; 3], r: f64, hh: f64) -> ConvexHull {
    let mut v = Vec::new();
    for k in 0..6 {
        let ang = std::f64::consts::TAU * (k as f64) / 6.0;
        let x = c[0] + r * ang.cos();
        let z = c[2] + r * ang.sin();
        v.push(Vec3::new(x, c[1] + hh, z));
        v.push(Vec3::new(x, c[1] - hh, z));
    }
    ConvexHull::new(v)
}

/// A regular tetrahedron with one face flat on the plane `y = base_y`, apex up.
///
/// The three base vertices sit at `y = base_y`, the apex above the centroid.
fn tetra_face_down(center_xz: [f64; 2], base_y: f64, scale: f64) -> ConvexHull {
    let (cx, cz) = (center_xz[0], center_xz[1]);
    // Equilateral base triangle in the y=base_y plane.
    let r = scale; // circumradius
    let mut v = Vec::new();
    for k in 0..3 {
        let ang = std::f64::consts::TAU * (k as f64) / 3.0 + std::f64::consts::FRAC_PI_2;
        v.push(Vec3::new(cx + r * ang.cos(), base_y, cz + r * ang.sin()));
    }
    // Apex above the centroid.
    let apex_h = scale * 1.6;
    v.push(Vec3::new(cx, base_y + apex_h, cz));
    ConvexHull::new(v)
}

#[test]
fn test_box_box_regression_four_contacts() {
    // Box resting flat on a wide ground box: the box-box manifold (now routed
    // through the shared convex face-clip helper) must still yield exactly 4
    // corner contacts with a shared +y normal and equal depth — the no-rock
    // guarantee, unchanged by the helper refactor.
    let dispatcher = NarrowPhaseDispatcher::default();
    let cube = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let ground = BoxShape::new(Vec3::new(10.0, 0.5, 10.0));
    let t_ground = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
    let t_cube = Transform::from_position(Vec3::new(0.0, 0.5 - 0.005, 0.0));

    let result = dispatcher.dispatch(
        &cube,
        ShapeType::Box,
        &t_cube,
        &ground,
        ShapeType::Box,
        &t_ground,
        CollisionPair::new(0, 1),
    );
    let m = result.manifold.expect("overlapping boxes -> manifold");
    assert_eq!(m.contacts.len(), 4, "box-on-ground must be a 4-point manifold");

    let first_depth = m.contacts[0].depth;
    for c in &m.contacts {
        assert!(c.normal.y > 0.999, "normal points +y (B->A)");
        assert!(c.normal.x.abs() < 1e-3 && c.normal.z.abs() < 1e-3);
        assert!(
            (c.depth - first_depth).abs() < 1e-6,
            "all corner depths equal within 1e-6"
        );
    }
}

#[test]
fn test_tetra_on_box_three_point_manifold() {
    // Tetrahedron face-down resting on a box: face-face contact -> 3-point
    // manifold (the tetra base is a triangle), all depths >= 0.
    let dispatcher = NarrowPhaseDispatcher::default();
    let ground = BoxShape::new(Vec3::new(5.0, 0.5, 5.0));
    let t_ground = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
    // Base sits just below y=0 so it overlaps the ground top (y=0) slightly.
    let tetra = tetra_face_down([0.0, 0.0], -0.01, 0.8);
    let t_tetra = Transform::from_position(Vec3::zeros());

    let result = dispatcher.dispatch(
        &tetra,
        ShapeType::ConvexHull,
        &t_tetra,
        &ground,
        ShapeType::Box,
        &t_ground,
        CollisionPair::new(0, 1),
    );
    let m = result.manifold.expect("overlapping tetra/box -> manifold");
    assert!(
        m.contacts.len() >= 3,
        "tetra face-down should give a 3-point manifold, got {}",
        m.contacts.len()
    );
    for c in &m.contacts {
        assert!(c.depth >= 0.0, "all depths non-negative");
        assert!(c.normal.y.abs() > 0.9, "normal roughly vertical");
    }
}

#[test]
fn test_tetra_on_box_com_stationary() {
    // A tetra resting on a static box: a sequential-impulse solver driven with a
    // downward velocity each step must converge to zero angular velocity and zero
    // horizontal drift over 30 iterations (the manifold holds the COM stationary).
    let dispatcher = NarrowPhaseDispatcher::default();
    let ground = BoxShape::new(Vec3::new(5.0, 0.5, 5.0));
    let t_ground = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
    let tetra = tetra_face_down([0.0, 0.0], -0.01, 0.8);
    let t_tetra = Transform::from_position(Vec3::zeros());

    let m = dispatcher
        .dispatch(
            &tetra,
            ShapeType::ConvexHull,
            &t_tetra,
            &ground,
            ShapeType::Box,
            &t_ground,
            CollisionPair::new(0, 1),
        )
        .manifold
        .expect("manifold");
    assert!(m.contacts.len() >= 3);

    // Tetra is body A (dynamic), ground is body B (static).
    let com_a = [0.0_f64, 0.2, 0.0]; // approximate tetra centroid
    let inv_mass_a = 1.0_f64;
    let inv_inertia_a = [5.0_f64, 5.0, 5.0];

    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let mul_diag = |d: [f64; 3], v: [f64; 3]| [d[0] * v[0], d[1] * v[1], d[2] * v[2]];

    let mut v_a = [0.0_f64; 3];
    let mut w_a = [0.0_f64; 3];
    let mut acc = vec![0.0_f64; m.contacts.len()];

    for _ in 0..30 {
        v_a[1] = -1.0; // re-drive downward each step (gravity proxy)
        for _ in 0..12 {
            for (i, c) in m.contacts.iter().enumerate() {
                let r_a = [
                    c.point_a.x - com_a[0],
                    c.point_a.y - com_a[1],
                    c.point_a.z - com_a[2],
                ];
                let n = [c.normal.x, c.normal.y, c.normal.z];
                let wxr = cross(w_a, r_a);
                let v_point = [v_a[0] + wxr[0], v_a[1] + wxr[1], v_a[2] + wxr[2]];
                let vn = dot(v_point, n);
                let rxn = cross(r_a, n);
                let denom = inv_mass_a + dot(n, cross(mul_diag(inv_inertia_a, rxn), r_a));
                if denom <= 0.0 {
                    continue;
                }
                let d_lambda = -vn / denom;
                let old = acc[i];
                let new = (old + d_lambda).max(0.0);
                let applied = new - old;
                acc[i] = new;
                for k in 0..3 {
                    v_a[k] += n[k] * inv_mass_a * applied;
                }
                let dl = cross(r_a, [n[0] * applied, n[1] * applied, n[2] * applied]);
                let dw = mul_diag(inv_inertia_a, dl);
                for k in 0..3 {
                    w_a[k] += dw[k];
                }
            }
        }
    }

    let w_norm = (w_a[0] * w_a[0] + w_a[1] * w_a[1] + w_a[2] * w_a[2]).sqrt();
    assert!(w_norm < 1e-2, "resting tetra should not spin: |w| = {w_norm}");
    assert!(v_a[0].abs() < 1e-2, "no x drift: {}", v_a[0]);
    assert!(v_a[2].abs() < 1e-2, "no z drift: {}", v_a[2]);
}

#[test]
fn test_hull_hull_face_to_face() {
    // Two aligned hexagonal prisms stacked face-to-face (flat hexagon caps in
    // contact) must produce >= 3 contact points.
    let dispatcher = NarrowPhaseDispatcher::default();
    let lower = hex_prism([0.0, 0.0, 0.0], 1.0, 0.5);
    let upper = hex_prism([0.0, 0.99, 0.0], 1.0, 0.5);
    let t = Transform::from_position(Vec3::zeros());

    let result = dispatcher.dispatch(
        &upper,
        ShapeType::ConvexHull,
        &t,
        &lower,
        ShapeType::ConvexHull,
        &t,
        CollisionPair::new(0, 1),
    );
    let m = result.manifold.expect("overlapping hexes -> manifold");
    assert!(
        m.contacts.len() >= 3,
        "hex face-to-face should give >= 3 contacts, got {}",
        m.contacts.len()
    );
    for c in &m.contacts {
        assert!(c.normal.y.abs() > 0.9, "normal roughly vertical");
    }
}

/// A unit cube hull rotated 45 deg about the given principal axis, centred at `c`.
///
/// `axis`: 0 = x, 1 = y, 2 = z. Rotating a cube 45 deg about x (resp. z) turns
/// its top into a horizontal *edge* running along x (resp. z) rather than a face,
/// which is what sets up a genuine edge-edge contact.
fn cube_hull_rot45(c: [f64; 3], h: f64, axis: usize) -> ConvexHull {
    let ang = std::f64::consts::FRAC_PI_4;
    let (s, co) = (ang.sin(), ang.cos());
    let mut v = Vec::new();
    for sx in [-h, h] {
        for sy in [-h, h] {
            for sz in [-h, h] {
                let (rx, ry, rz) = match axis {
                    0 => (sx, sy * co - sz * s, sy * s + sz * co),
                    2 => (sx * co - sy * s, sx * s + sy * co, sz),
                    _ => (sx * co - sz * s, sy, sx * s + sz * co),
                };
                v.push(Vec3::new(c[0] + rx, c[1] + ry, c[2] + rz));
            }
        }
    }
    ConvexHull::new(v)
}

#[test]
fn test_edge_edge_two_boxes_at_45() {
    // Genuine edge-edge: cube A rotated 45 deg about x (top is a horizontal edge
    // along x); cube B rotated 45 deg about z (bottom is a horizontal edge along
    // z), placed just above A. The two top/bottom edges are perpendicular and
    // cross, so the contact is edge-edge with a vertical (+y) normal that is
    // perpendicular to BOTH edges. Neither body presents a face to +y.
    let dispatcher = NarrowPhaseDispatcher::default();
    let half = 0.5;
    // Rotated-cube half-height along y is h*sqrt(2) (corner-to-centre). Place B so
    // its bottom edge dips just below A's top edge for a shallow overlap.
    let top_a = half * std::f64::consts::SQRT_2;
    let a = cube_hull_rot45([0.0, 0.0, 0.0], half, 0);
    let b = cube_hull_rot45([0.0, 2.0 * top_a - 0.02, 0.0], half, 2);
    let t = Transform::from_position(Vec3::zeros());

    let result = dispatcher.dispatch(
        &a,
        ShapeType::ConvexHull,
        &t,
        &b,
        ShapeType::ConvexHull,
        &t,
        CollisionPair::new(0, 1),
    );
    let m = result.manifold.expect("edge-edge overlap -> manifold");
    assert!(
        (1..=2).contains(&m.contacts.len()),
        "edge-edge should give 1-2 contacts, got {}",
        m.contacts.len()
    );
    for c in &m.contacts {
        // A's top edge runs along x, B's bottom edge along z; the separation
        // normal is perpendicular to both -> vertical (~+y).
        assert!(
            c.normal.y.abs() > 0.7,
            "edge-edge normal should be ~vertical (perp to both edges), got y={}",
            c.normal.y
        );
        assert!(c.normal.x.abs() < 0.5 && c.normal.z.abs() < 0.5);
    }
}

#[test]
fn test_both_dispatch_orders_consistent() {
    // hull x box vs box x hull: opposite normals, geometrically identical
    // contacts (witness points swapped). Validates the dispatcher flip path for
    // the convex manifold arm.
    let dispatcher = NarrowPhaseDispatcher::default();
    let hull = cube_hull([0.0, 0.99, 0.0], 0.5);
    let box_b = BoxShape::new(Vec3::new(0.5, 0.5, 0.5));
    let t_hull = Transform::from_position(Vec3::zeros());
    let t_box = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));

    let m1 = dispatcher
        .dispatch(
            &hull,
            ShapeType::ConvexHull,
            &t_hull,
            &box_b,
            ShapeType::Box,
            &t_box,
            CollisionPair::new(0, 1),
        )
        .manifold
        .expect("hull x box");
    let m2 = dispatcher
        .dispatch(
            &box_b,
            ShapeType::Box,
            &t_box,
            &hull,
            ShapeType::ConvexHull,
            &t_hull,
            CollisionPair::new(1, 0),
        )
        .manifold
        .expect("box x hull");

    assert!(!m1.contacts.is_empty());
    assert_eq!(m1.contacts.len(), m2.contacts.len());

    for c1 in &m1.contacts {
        let p1 = c1.point();
        let best = m2
            .contacts
            .iter()
            .min_by(|a, b| {
                (p1 - a.point())
                    .norm()
                    .partial_cmp(&(p1 - b.point()).norm())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("nonempty");
        assert!(
            (best.normal + c1.normal).norm() < 1e-6,
            "normals should be opposite between dispatch orders"
        );
        assert!((c1.point_a - best.point_b).norm() < 1e-6, "A<->B witness swap");
        assert!((c1.point_b - best.point_a).norm() < 1e-6, "A<->B witness swap");
    }
}

#[test]
fn test_separated_hulls_no_contact() {
    let dispatcher = NarrowPhaseDispatcher::default();
    let a = cube_hull([0.0, 0.0, 0.0], 0.5);
    let b = cube_hull([10.0, 0.0, 0.0], 0.5);
    let t = Transform::from_position(Vec3::zeros());
    let result = dispatcher.dispatch(
        &a,
        ShapeType::ConvexHull,
        &t,
        &b,
        ShapeType::ConvexHull,
        &t,
        CollisionPair::new(0, 1),
    );
    assert!(
        result.manifold.is_none(),
        "well-separated hulls produce no contact"
    );
}
