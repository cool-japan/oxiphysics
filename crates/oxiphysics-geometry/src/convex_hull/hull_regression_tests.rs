// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the `ConvexHull` shape (volume, centre of mass,
//! inertia tensor, ray casts) and for the quickhull builder it relies on.

use std::collections::HashMap;

use oxiphysics_core::math::{Mat3, Vec3, quat_from_axis_angle};

use super::functions::{cross, dot, normalize, sub};
use super::types::{ConvexHull, ConvexHull3D};
use crate::box_shape::BoxShape;
use crate::shape::Shape;

/// Deterministic xorshift generator returning values in `[0, 1)`.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = (self.next() * (i as f64 + 1.0)) as usize;
            items.swap(i, j.min(i));
        }
    }
}

fn cube_points(lo: f64, hi: f64) -> Vec<[f64; 3]> {
    let mut pts = Vec::new();
    for &x in &[lo, hi] {
        for &y in &[lo, hi] {
            for &z in &[lo, hi] {
                pts.push([x, y, z]);
            }
        }
    }
    pts
}

fn to_vec3(points: &[[f64; 3]]) -> Vec<Vec3> {
    points.iter().map(|p| Vec3::new(p[0], p[1], p[2])).collect()
}

fn extent_of(points: &[[f64; 3]]) -> f64 {
    (0..3)
        .map(|axis| {
            let lo = points.iter().map(|p| p[axis]).fold(f64::INFINITY, f64::min);
            let hi = points
                .iter()
                .map(|p| p[axis])
                .fold(f64::NEG_INFINITY, f64::max);
            hi - lo
        })
        .fold(0.0, f64::max)
}

/// The hull surface must be a closed, consistently oriented 2-manifold
/// (every directed edge used once and matched by its reverse, Euler
/// characteristic 2) whose facets all face outward: no input point may lie
/// in front of a facet plane.
fn assert_closed_outward_hull(hull: &ConvexHull3D, points: &[[f64; 3]], label: &str) {
    let mut directed: HashMap<(usize, usize), usize> = HashMap::new();
    for face in &hull.faces {
        for k in 0..3 {
            *directed.entry((face[k], face[(k + 1) % 3])).or_insert(0) += 1;
        }
    }
    for (&(a, b), &count) in &directed {
        assert_eq!(
            count, 1,
            "{label}: directed edge ({a},{b}) used {count} times"
        );
        assert_eq!(
            directed.get(&(b, a)).copied(),
            Some(1),
            "{label}: edge ({a},{b}) has no opposite edge"
        );
    }
    let v = hull.vertices.len() as i64;
    let e = (directed.len() / 2) as i64;
    let f = hull.faces.len() as i64;
    assert_eq!(
        v - e + f,
        2,
        "{label}: Euler characteristic (V={v}, E={e}, F={f})"
    );
    let tolerance = 1e-9 * extent_of(points);
    for face in &hull.faces {
        let a = hull.vertices[face[0]];
        let normal = normalize(cross(
            sub(hull.vertices[face[1]], a),
            sub(hull.vertices[face[2]], a),
        ));
        assert!(
            dot(normal, normal) > 0.5,
            "{label}: degenerate facet {face:?}"
        );
        for &p in points {
            assert!(
                dot(normal, sub(p, a)) <= tolerance,
                "{label}: point {p:?} lies outside facet {face:?}"
            );
        }
    }
}

// ── Builder (ConvexHull3D) ────────────────────────────────────────────────

#[test]
fn test_builder_cube_in_any_vertex_order_is_a_closed_outward_hull() {
    let mut rng = Rng(0x1234_5678_9ABC_DEF1);
    for round in 0..40 {
        let mut pts = cube_points(0.0, 1.0);
        rng.shuffle(&mut pts);
        let hull = ConvexHull3D::build(&pts).expect("cube hull");
        assert_closed_outward_hull(&hull, &pts, &format!("cube round {round}"));
        assert_eq!(hull.faces.len(), 12, "round {round}");
        assert_eq!(hull.vertices.len(), 8, "round {round}");
        assert!((hull.compute_volume() - 1.0).abs() < 1e-12, "round {round}");
    }
}

#[test]
fn test_builder_regular_solids() {
    let tetra = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let hull = ConvexHull3D::build(&tetra).expect("tetrahedron");
    assert_closed_outward_hull(&hull, &tetra, "tetrahedron");
    assert!((hull.compute_volume() - 1.0 / 6.0).abs() < 1e-14);

    let octahedron = vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let hull = ConvexHull3D::build(&octahedron).expect("octahedron");
    assert_closed_outward_hull(&hull, &octahedron, "octahedron");
    assert_eq!(hull.faces.len(), 8);
    assert!((hull.compute_volume() - 4.0 / 3.0).abs() < 1e-13);
}

#[test]
fn test_builder_coplanar_rich_grid_points() {
    // 3x3x3 lattice: every face of the cube carries nine coplanar points and
    // every edge three collinear ones.
    let mut pts = Vec::new();
    for x in 0..3 {
        for y in 0..3 {
            for z in 0..3 {
                pts.push([x as f64, y as f64, z as f64]);
            }
        }
    }
    let mut rng = Rng(0xDEAD_BEEF_CAFE_F00D);
    for round in 0..20 {
        rng.shuffle(&mut pts);
        let hull = ConvexHull3D::build(&pts).expect("grid hull");
        assert_closed_outward_hull(&hull, &pts, &format!("grid round {round}"));
        assert!((hull.compute_volume() - 8.0).abs() < 1e-10, "round {round}");
    }
}

#[test]
fn test_builder_rotated_lattice_with_jittered_duplicates() {
    // The 3x3x3 lattice rotated by a generic rotation is only coplanar up to
    // rounding error, and every point is repeated with a 1e-14 jitter: the
    // sign of a distance to a "coplanar" vertex is noise.  The hull must
    // still be a closed outward surface of volume 8.
    let rotation = quat_from_axis_angle(&Vec3::new(1.0, 2.0, 3.0), 0.7)
        .to_rotation_matrix()
        .into_inner();
    let mut base = Vec::new();
    for x in 0..3 {
        for y in 0..3 {
            for z in 0..3 {
                let v = rotation * Vec3::new(x as f64, y as f64, z as f64);
                base.push([v.x, v.y, v.z]);
            }
        }
    }
    let mut rng = Rng(0xC0FF_EE00_5EED_1234);
    for round in 0..25 {
        let mut pts = base.clone();
        let jittered: Vec<[f64; 3]> = base
            .iter()
            .map(|p| {
                [
                    p[0] + 1e-14 * (rng.next() - 0.5),
                    p[1] + 1e-14 * (rng.next() - 0.5),
                    p[2] + 1e-14 * (rng.next() - 0.5),
                ]
            })
            .collect();
        pts.extend(jittered);
        rng.shuffle(&mut pts);
        let hull = ConvexHull3D::build(&pts).expect("rotated lattice hull");
        assert_closed_outward_hull(&hull, &pts, &format!("rotated lattice round {round}"));
        assert!(
            (hull.compute_volume() - 8.0).abs() < 1e-8,
            "round {round}: volume {}",
            hull.compute_volume()
        );
    }
}

#[test]
fn test_builder_thin_sliver_facets_do_not_hide_outside_points() {
    // A long thin needle: its facets have tiny area, so a threshold on the
    // area-weighted distance would swallow points that stick out of them.
    let pts = vec![
        [0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [0.0, 1e-3, 0.0],
        [0.0, 0.0, 1e-3],
        [5.0, 2e-3, 2e-3],
        [10.0, 1e-3, 1e-3],
    ];
    let hull = ConvexHull3D::build(&pts).expect("needle hull");
    assert_closed_outward_hull(&hull, &pts, "needle");
    for p in &pts {
        assert!(
            hull.contains_point(*p) || hull.vertices.contains(p),
            "input point {p:?} missing from the hull"
        );
    }
}

#[test]
fn test_builder_random_clouds_and_sphere_points() {
    let mut rng = Rng(0x0BAD_C0FF_EE12_3457);
    for &n in &[12_usize, 40, 150] {
        let pts: Vec<[f64; 3]> = (0..n)
            .map(|_| {
                [
                    rng.next() * 4.0 - 2.0,
                    rng.next() * 2.0 - 1.0,
                    rng.next() * 6.0,
                ]
            })
            .collect();
        let hull = ConvexHull3D::build(&pts).expect("random cloud hull");
        assert_closed_outward_hull(&hull, &pts, &format!("random cloud n={n}"));
    }
    // Fibonacci lattice on the unit sphere: all points are hull vertices.
    let n = 120;
    let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    let sphere: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let y = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
            let r = (1.0 - y * y).sqrt();
            [
                r * (golden * i as f64).cos(),
                y,
                r * (golden * i as f64).sin(),
            ]
        })
        .collect();
    let hull = ConvexHull3D::build(&sphere).expect("sphere hull");
    assert_closed_outward_hull(&hull, &sphere, "sphere");
    assert_eq!(hull.vertices.len(), n);
    let volume = hull.compute_volume();
    let ball = 4.0 / 3.0 * std::f64::consts::PI;
    assert!(
        volume < ball && volume > 0.95 * ball,
        "sphere hull volume {volume}"
    );
}

#[test]
fn test_builder_is_scale_invariant() {
    let reference = {
        let pts = cube_points(0.0, 1.0);
        ConvexHull3D::build(&pts).expect("unit cube").faces.len()
    };
    for &scale in &[1e-6, 1e-4, 1.0, 1e5] {
        let pts: Vec<[f64; 3]> = cube_points(0.0, 1.0)
            .iter()
            .map(|p| [p[0] * scale, p[1] * scale, p[2] * scale])
            .collect();
        let hull = ConvexHull3D::build(&pts).expect("scaled cube");
        assert_eq!(hull.faces.len(), reference, "scale {scale}");
        let want = scale * scale * scale;
        assert!(
            (hull.compute_volume() - want).abs() <= 1e-10 * want,
            "scale {scale}: {}",
            hull.compute_volume()
        );
    }
}

#[test]
fn test_builder_rejects_degenerate_and_non_finite_input() {
    let flat: Vec<[f64; 3]> = cube_points(0.0, 1.0)
        .iter()
        .map(|p| [p[0], p[1], 0.5])
        .collect();
    assert!(ConvexHull3D::build(&flat).is_none());
    let mut bad = cube_points(0.0, 1.0);
    bad[5][2] = f64::NAN;
    assert!(ConvexHull3D::build(&bad).is_none());
    bad[5][2] = f64::NEG_INFINITY;
    assert!(ConvexHull3D::build(&bad).is_none());
    let same = vec![[1.0, 2.0, 3.0]; 6];
    assert!(ConvexHull3D::build(&same).is_none());
}

// ── Shape implementation (ConvexHull) ─────────────────────────────────────

fn box_hull(half: [f64; 3], offset: Vec3) -> ConvexHull {
    let mut vertices = Vec::new();
    for &sx in &[-1.0, 1.0] {
        for &sy in &[-1.0, 1.0] {
            for &sz in &[-1.0, 1.0] {
                vertices.push(offset + Vec3::new(sx * half[0], sy * half[1], sz * half[2]));
            }
        }
    }
    ConvexHull::new(vertices)
}

fn assert_mat_close(got: &Mat3, want: &Mat3, tol: f64) {
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (got[(i, j)] - want[(i, j)]).abs() <= tol,
                "[{i}][{j}] = {} expected {}",
                got[(i, j)],
                want[(i, j)]
            );
        }
    }
}

#[test]
fn test_unit_cube_shape_mass_properties_ignore_vertex_order() {
    let mut rng = Rng(0x5EED_0000_0000_0001);
    for round in 0..30 {
        let mut pts = cube_points(0.0, 1.0);
        rng.shuffle(&mut pts);
        let hull = ConvexHull::new(to_vec3(&pts));
        assert!((hull.volume() - 1.0).abs() < 1e-12, "round {round}");
        assert!(
            (hull.center_of_mass() - Vec3::new(0.5, 0.5, 0.5)).norm() < 1e-12,
            "round {round}"
        );
        // Mass 6 turns I = m/6 into the identity.
        assert_mat_close(&hull.inertia_tensor(6.0), &Mat3::identity(), 1e-11);
    }
}

#[test]
fn test_unit_tetrahedron_shape_volume() {
    let hull = ConvexHull::new(vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]);
    assert!((hull.volume() - 1.0 / 6.0).abs() < 1e-14);
    assert!((hull.center_of_mass() - Vec3::new(0.25, 0.25, 0.25)).norm() < 1e-14);
}

#[test]
fn test_box_hull_agrees_with_box_shape() {
    let half = [1.0, 2.0, 3.0];
    let hull = box_hull(half, Vec3::zeros());
    let reference = BoxShape::new(Vec3::new(half[0], half[1], half[2]));
    assert!((hull.volume() - reference.volume()).abs() < 1e-10);
    assert!((hull.center_of_mass() - reference.center_of_mass()).norm() < 1e-12);
    assert_mat_close(
        &hull.inertia_tensor(9.0),
        &reference.inertia_tensor(9.0),
        1e-9,
    );
    // Default mass properties derive from the same three methods.
    let props = hull.mass_properties(100.0);
    assert!((props.mass - 100.0 * 48.0).abs() < 1e-8);
    assert!((props.center_of_mass - reference.center_of_mass()).norm() < 1e-12);
    assert_mat_close(
        &props.local_inertia,
        &reference.inertia_tensor(props.mass),
        1e-6,
    );
}

#[test]
fn test_rotated_box_hull_inertia_has_cross_terms() {
    let half = [1.0, 2.0, 3.0];
    let rotation = quat_from_axis_angle(&Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_4)
        .to_rotation_matrix()
        .into_inner();
    let straight = box_hull(half, Vec3::zeros());
    let turned = ConvexHull::new(straight.vertices.iter().map(|v| rotation * v).collect());
    let mass = 5.0;
    let local = BoxShape::new(Vec3::new(half[0], half[1], half[2])).inertia_tensor(mass);
    let want = rotation * local * rotation.transpose();
    let got = turned.inertia_tensor(mass);
    assert_mat_close(&got, &want, 1e-9);
    assert!(got[(0, 1)].abs() > 1.0, "I_xy = {}", got[(0, 1)]);
    assert!((got[(0, 1)] - got[(1, 0)]).abs() < 1e-12);
}

#[test]
fn test_translated_hull_moves_com_but_not_inertia() {
    let half = [0.5, 1.0, 1.5];
    let at_origin = box_hull(half, Vec3::zeros());
    let moved = box_hull(half, Vec3::new(10.0, -4.0, 7.0));
    assert!((moved.volume() - at_origin.volume()).abs() < 1e-10);
    assert!((moved.center_of_mass() - Vec3::new(10.0, -4.0, 7.0)).norm() < 1e-11);
    assert_mat_close(
        &moved.inertia_tensor(3.0),
        &at_origin.inertia_tensor(3.0),
        1e-9,
    );
}

#[test]
fn test_random_cloud_agrees_with_independent_polyhedral_integrals() {
    // `ConvexHull3D::compute_*` integrate with an independent formulation
    // (origin-based tetrahedra, Eberly's second moments).
    let mut rng = Rng(0xA11C_E000_0000_0042);
    for case in 0..12 {
        let n = 8 + (rng.next() * 40.0) as usize;
        let pts: Vec<[f64; 3]> = (0..n)
            .map(|_| {
                [
                    1.0 + rng.next() * 3.0,
                    -2.0 + rng.next() * 2.0,
                    0.5 + rng.next() * 1.5,
                ]
            })
            .collect();
        let reference = ConvexHull3D::build(&pts).expect("cloud hull");
        let hull = ConvexHull::new(to_vec3(&pts));
        let volume = reference.compute_volume();
        assert!(
            (hull.volume() - volume).abs() <= 1e-9 * volume,
            "case {case}: volume {} vs {volume}",
            hull.volume()
        );
        let com = reference.compute_center_of_mass();
        assert!(
            (hull.center_of_mass() - Vec3::new(com[0], com[1], com[2])).norm() < 1e-9,
            "case {case}: com"
        );
        let mass = 3.5;
        let inertia = reference.compute_moment_of_inertia(mass);
        let got = hull.inertia_tensor(mass);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (got[(i, j)] - inertia[i][j]).abs() <= 1e-8 * (1.0 + inertia[i][i].abs()),
                    "case {case}: I[{i}][{j}] = {} vs {}",
                    got[(i, j)],
                    inertia[i][j]
                );
            }
        }
        // A physical inertia tensor is symmetric positive definite and obeys
        // the triangle inequality of principal moments.
        assert!((got[(0, 1)] - got[(1, 0)]).abs() < 1e-12);
        assert!(got.symmetric_eigenvalues().iter().all(|&l| l > 0.0));
    }
}

#[test]
fn test_ray_cast_axis_ray_hits_cube_face() {
    let hull = box_hull([0.5, 0.5, 0.5], Vec3::zeros());
    let hit = hull
        .ray_cast(&Vec3::new(-3.0, 0.0, 0.0), &Vec3::new(1.0, 0.0, 0.0), 100.0)
        .expect("axis ray hits the cube");
    assert!(
        (hit.point.x + 0.5).abs() < 1e-12,
        "hit at x = {}",
        hit.point.x
    );
    assert!((hit.toi - 2.5).abs() < 1e-12);
    assert!((hit.normal - Vec3::new(-1.0, 0.0, 0.0)).norm() < 1e-12);
}

#[test]
fn test_ray_cast_matches_box_shape_for_outside_origins() {
    let half = [0.8, 1.1, 0.6];
    let hull = box_hull(half, Vec3::zeros());
    let reference = BoxShape::new(Vec3::new(half[0], half[1], half[2]));
    let mut rng = Rng(0xFEED_FACE_0000_0007);
    let mut compared = 0;
    for _ in 0..400 {
        let origin = Vec3::new(
            rng.next() * 8.0 - 4.0,
            rng.next() * 8.0 - 4.0,
            rng.next() * 8.0 - 4.0,
        );
        if origin.x.abs() < half[0] && origin.y.abs() < half[1] && origin.z.abs() < half[2] {
            continue;
        }
        // Aim roughly at the box so that many rays hit.
        let target = Vec3::new(
            (rng.next() - 0.5) * 2.0 * half[0],
            (rng.next() - 0.5) * 2.0 * half[1],
            (rng.next() - 0.5) * 2.0 * half[2],
        );
        let direction = (target - origin) * (0.5 + rng.next());
        let got = hull.ray_cast(&origin, &direction, 100.0);
        let want = reference.ray_cast(&origin, &direction, 100.0);
        match (got, want) {
            (Some(g), Some(w)) => {
                compared += 1;
                assert!((g.toi - w.toi).abs() < 1e-9, "toi {} vs {}", g.toi, w.toi);
                assert!((g.point - w.point).norm() < 1e-9);
                assert!(
                    (g.normal - w.normal).norm() < 1e-9,
                    "{:?} vs {:?}",
                    g.normal,
                    w.normal
                );
            }
            (None, None) => {}
            (g, w) => panic!("hull {g:?} vs box {w:?} for origin {origin:?} dir {direction:?}"),
        }
    }
    assert!(compared > 100, "only {compared} rays hit the box");
}

#[test]
fn test_ray_cast_misses_inside_and_limits() {
    let hull = box_hull([0.5, 0.5, 0.5], Vec3::zeros());
    let along_x = Vec3::new(1.0, 0.0, 0.0);
    assert!(
        hull.ray_cast(&Vec3::new(-3.0, 0.8, 0.0), &along_x, 100.0)
            .is_none()
    );
    assert!(
        hull.ray_cast(&Vec3::new(3.0, 0.0, 0.0), &along_x, 100.0)
            .is_none()
    );
    assert!(
        hull.ray_cast(&Vec3::new(-3.0, 0.0, 0.0), &along_x, 2.0)
            .is_none()
    );
    // Starting inside: reports the exit through the +y face.
    let exit = hull
        .ray_cast(&Vec3::new(0.0, 0.2, 0.0), &Vec3::new(0.0, 1.0, 0.0), 10.0)
        .expect("exit hit");
    assert!((exit.toi - 0.3).abs() < 1e-12);
    assert!((exit.normal - Vec3::new(0.0, 1.0, 0.0)).norm() < 1e-12);
    // Grazing along a face counts as a miss only if outside the slab.
    assert!(
        hull.ray_cast(&Vec3::new(-3.0, 0.5, 0.0), &along_x, 100.0)
            .is_some()
    );
}

#[test]
fn test_degenerate_hulls_have_no_volume_and_no_ray_hits() {
    let flat = ConvexHull::new(vec![
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 1.0),
        Vec3::new(1.0, 1.0, 1.0),
        Vec3::new(0.0, 1.0, 1.0),
    ]);
    assert_eq!(flat.volume(), 0.0);
    assert!((flat.center_of_mass() - Vec3::new(0.5, 0.5, 1.0)).norm() < 1e-14);
    assert_mat_close(&flat.inertia_tensor(2.0), &Mat3::zeros(), 0.0);
    assert!(
        flat.ray_cast(&Vec3::new(0.5, 0.5, 5.0), &Vec3::new(0.0, 0.0, -1.0), 100.0)
            .is_none()
    );
    let empty = ConvexHull::new(Vec::new());
    assert_eq!(empty.volume(), 0.0);
    assert_eq!(empty.center_of_mass(), Vec3::zeros());
    let pair = ConvexHull::new(vec![Vec3::new(-1.0, -2.0, -3.0), Vec3::new(1.0, 2.0, 3.0)]);
    assert_eq!(pair.volume(), 0.0);
    assert_eq!(pair.center_of_mass(), Vec3::zeros());
    let nan = ConvexHull::new(vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, f64::NAN),
    ]);
    assert_eq!(nan.volume(), 0.0);
}

#[test]
fn test_vertex_edits_after_first_query_are_respected() {
    let mut hull = box_hull([0.5, 0.5, 0.5], Vec3::zeros());
    assert!((hull.volume() - 1.0).abs() < 1e-12);
    // `vertices` is public: stretch the box along x after the facets were cached.
    for v in &mut hull.vertices {
        v.x *= 4.0;
    }
    assert!((hull.volume() - 4.0).abs() < 1e-12);
    let hit = hull
        .ray_cast(&Vec3::new(-5.0, 0.0, 0.0), &Vec3::new(1.0, 0.0, 0.0), 100.0)
        .expect("stretched box is hit at x = -2");
    assert!((hit.point.x + 2.0).abs() < 1e-12);
    // A flat vertex set that later gains volume.
    let mut growing = ConvexHull::new(vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(1.0, 1.0, 0.0),
    ]);
    assert_eq!(growing.volume(), 0.0);
    growing.vertices.push(Vec3::new(0.5, 0.5, 1.0));
    assert!(
        growing.volume() > 0.3,
        "pyramid volume {}",
        growing.volume()
    );
}

#[test]
fn test_clone_and_debug() {
    let hull = box_hull([1.0, 1.0, 1.0], Vec3::zeros());
    let _ = hull.volume();
    let copy = hull.clone();
    assert!((copy.volume() - 8.0).abs() < 1e-12);
    let text = format!("{hull:?}");
    assert!(text.starts_with("ConvexHull"), "{text}");
    assert!(text.contains("vertices"), "{text}");
}
