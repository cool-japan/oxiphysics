// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the general convex routes of the `ShapeKind` layer.
//!
//! Every pair without a closed-form routine (box/box, capsule/box, convex hulls,
//! sphere/convex ...) is resolved by GJK + EPA over the shapes' support
//! mappings. It used to be a three-axis world SAT whose normal always pointed
//! along `+axis`, so rotated shapes that were separated along a diagonal
//! reported contacts and the normal ignored where the other shape lay.
//! `convex_vs_mesh` used to replace each triangle by its bounding sphere.
//!
//! Contacts follow the `NarrowPhaseContact` convention: the normal points from
//! shape A towards shape B, `point_a` lies on A, `point_b` on B and
//! `point_a - point_b == normal * depth`.

use oxiphysics_collision::narrowphase::{TriangleMesh, convex_vs_mesh};
use oxiphysics_collision::{NarrowPhaseContact, ShapeKind, shape_shape_contact};

type Vec3 = [f64; 3];

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

/// A cube of half-extent `half`, turned by `yaw` about the vertical axis, as a
/// convex vertex cloud.
fn cube(center: Vec3, half: f64, yaw: f64) -> ShapeKind {
    let (s, c) = yaw.sin_cos();
    let mut vertices = Vec::new();
    for sx in [-half, half] {
        for sy in [-half, half] {
            for sz in [-half, half] {
                vertices.push([
                    center[0] + c * sx + s * sz,
                    center[1] + sy,
                    center[2] - s * sx + c * sz,
                ]);
            }
        }
    }
    ShapeKind::Convex { vertices }
}

/// The tetrahedron `x, y, z >= 0, x + y + z <= 1`, translated by `shift`.
fn tetrahedron(shift: Vec3) -> ShapeKind {
    let corners: [Vec3; 4] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    ShapeKind::Convex {
        vertices: corners
            .iter()
            .map(|c| [c[0] + shift[0], c[1] + shift[1], c[2] + shift[2]])
            .collect(),
    }
}

fn aabbs_overlap(a: &ShapeKind, b: &ShapeKind) -> bool {
    let (a_min, a_max) = a.aabb();
    let (b_min, b_max) = b.aabb();
    (0..3).all(|i| a_max[i] >= b_min[i] && b_max[i] >= a_min[i])
}

/// The contact must be internally consistent: unit normal, and the witness
/// points are `normal * depth` apart with `point_a` on the A side.
fn assert_consistent(contact: &NarrowPhaseContact, tolerance: f64) {
    assert!(
        (norm(contact.normal) - 1.0).abs() < tolerance,
        "normal {:?}",
        contact.normal
    );
    assert!(contact.depth >= 0.0, "depth {}", contact.depth);
    let gap = sub(contact.point_a, contact.point_b);
    let expected = [
        contact.normal[0] * contact.depth,
        contact.normal[1] * contact.depth,
        contact.normal[2] * contact.depth,
    ];
    assert!(
        norm(sub(gap, expected)) < tolerance,
        "point_a - point_b = {gap:?}, normal * depth = {expected:?}"
    );
}

fn assert_close(actual: Vec3, expected: Vec3, tolerance: f64) {
    assert!(
        norm(sub(actual, expected)) < tolerance,
        "{actual:?} differs from {expected:?}"
    );
}

#[test]
fn turned_cubes_separated_along_a_diagonal_do_not_touch_despite_overlapping_aabbs() {
    // Two cubes turned by 45 degrees: their faces face (1, 0, 1). Along that
    // axis each extends one unit, but the centres are 1.5 * sqrt(2) = 2.12 apart.
    let a = cube([0.0; 3], 1.0, std::f64::consts::FRAC_PI_4);
    let b = cube([1.5, 0.0, 1.5], 1.0, std::f64::consts::FRAC_PI_4);
    assert!(aabbs_overlap(&a, &b), "the bounding boxes must overlap");
    assert!(shape_shape_contact(&a, &b).is_none());
    assert!(shape_shape_contact(&b, &a).is_none());
}

#[test]
fn turned_cubes_overlapping_along_the_diagonal_report_the_face_normal_towards_b() {
    let a = cube([0.0; 3], 1.0, std::f64::consts::FRAC_PI_4);
    let b = cube([1.3, 0.0, 1.3], 1.0, std::f64::consts::FRAC_PI_4);
    let contact = shape_shape_contact(&a, &b).expect("the cubes overlap");
    let axis = std::f64::consts::FRAC_1_SQRT_2;
    let expected_depth = 2.0 - 2.6 * axis;
    assert!(
        (contact.depth - expected_depth).abs() < 1e-6,
        "depth {} vs {expected_depth}",
        contact.depth
    );
    assert_close(contact.normal, [axis, 0.0, axis], 1e-6);
    assert_consistent(&contact, 1e-6);
    // Mirrored: B first, so the normal flips and the witnesses swap.
    let mirrored = shape_shape_contact(&b, &a).expect("the cubes overlap");
    assert_close(mirrored.normal, [-axis, 0.0, -axis], 1e-6);
    assert!((mirrored.depth - contact.depth).abs() < 1e-6);
    assert_consistent(&mirrored, 1e-6);
}

#[test]
fn the_normal_points_towards_b_on_either_side() {
    let a = cube([0.0; 3], 1.0, 0.0);
    let right = shape_shape_contact(&a, &cube([1.5, 0.2, -0.1], 1.0, 0.0)).expect("overlap");
    assert!((right.depth - 0.5).abs() < 1e-9, "depth {}", right.depth);
    assert_close(right.normal, [1.0, 0.0, 0.0], 1e-9);
    assert_consistent(&right, 1e-9);
    let left = shape_shape_contact(&a, &cube([-1.5, 0.2, -0.1], 1.0, 0.0)).expect("overlap");
    assert!((left.depth - 0.5).abs() < 1e-9, "depth {}", left.depth);
    assert_close(left.normal, [-1.0, 0.0, 0.0], 1e-9);
    assert_consistent(&left, 1e-9);
    // The same holds along the other two axes.
    let above = shape_shape_contact(&a, &cube([0.1, 1.7, 0.0], 1.0, 0.0)).expect("overlap");
    assert!((above.depth - 0.3).abs() < 1e-9);
    assert_close(above.normal, [0.0, 1.0, 0.0], 1e-9);
    let behind = shape_shape_contact(&a, &cube([0.0, 0.1, -1.8], 1.0, 0.0)).expect("overlap");
    assert!((behind.depth - 0.2).abs() < 1e-9);
    assert_close(behind.normal, [0.0, 0.0, -1.0], 1e-9);
}

#[test]
fn box_box_contact_has_the_axis_depth_and_witnesses_on_the_faces() {
    let a = ShapeKind::Box {
        center: [0.0; 3],
        half_extents: [1.0, 1.0, 1.0],
    };
    let b = ShapeKind::Box {
        center: [1.5, 0.0, 0.0],
        half_extents: [1.0, 1.0, 1.0],
    };
    let contact = shape_shape_contact(&a, &b).expect("the boxes overlap");
    assert!(
        (contact.depth - 0.5).abs() < 1e-9,
        "depth {}",
        contact.depth
    );
    assert_close(contact.normal, [1.0, 0.0, 0.0], 1e-9);
    assert!(
        (contact.point_a[0] - 1.0).abs() < 1e-9,
        "{:?}",
        contact.point_a
    );
    assert!(
        (contact.point_b[0] - 0.5).abs() < 1e-9,
        "{:?}",
        contact.point_b
    );
    assert_consistent(&contact, 1e-9);
    let apart = ShapeKind::Box {
        center: [2.5, 0.0, 0.0],
        half_extents: [1.0, 1.0, 1.0],
    };
    assert!(shape_shape_contact(&a, &apart).is_none());
}

#[test]
fn capsule_against_box_is_resolved_along_the_closest_face() {
    let capsule = ShapeKind::Capsule {
        p0: [0.0, -1.0, 0.0],
        p1: [0.0, 1.0, 0.0],
        radius: 0.5,
    };
    let block = ShapeKind::Box {
        center: [1.2, 0.0, 0.0],
        half_extents: [1.0, 1.0, 1.0],
    };
    let contact = shape_shape_contact(&capsule, &block).expect("the capsule touches the box");
    // The capsule surface is curved, so EPA pins the depth to ~1e-9 but the
    // normal and the witnesses only to about its square root.
    assert!(
        (contact.depth - 0.3).abs() < 1e-6,
        "depth {}",
        contact.depth
    );
    assert_close(contact.normal, [1.0, 0.0, 0.0], 1e-3);
    assert!(
        (contact.point_a[0] - 0.5).abs() < 1e-3,
        "{:?}",
        contact.point_a
    );
    assert!(
        (contact.point_b[0] - 0.2).abs() < 1e-3,
        "{:?}",
        contact.point_b
    );
    assert_consistent(&contact, 1e-6);
    let flipped = shape_shape_contact(&block, &capsule).expect("symmetric");
    assert_close(flipped.normal, [-1.0, 0.0, 0.0], 1e-3);
    assert!((flipped.depth - 0.3).abs() < 1e-6);
}

#[test]
fn sphere_against_a_convex_cube_is_resolved_by_the_general_route() {
    let sphere = ShapeKind::Sphere {
        center: [0.0; 3],
        radius: 1.0,
    };
    let contact = shape_shape_contact(&sphere, &cube([1.5, 0.0, 0.0], 1.0, 0.0))
        .expect("the sphere reaches into the cube");
    assert!(
        (contact.depth - 0.5).abs() < 1e-6,
        "depth {}",
        contact.depth
    );
    assert_close(contact.normal, [1.0, 0.0, 0.0], 1e-3);
    assert_consistent(&contact, 1e-6);
}

#[test]
fn tetrahedra_separated_along_the_body_diagonal_do_not_touch() {
    // x + y + z <= 1 against x + y + z >= 1.8: separated, yet the bounding boxes
    // [0, 1]^3 and [0.6, 1.6]^3 overlap on every world axis.
    let a = tetrahedron([0.0; 3]);
    let b = tetrahedron([0.6; 3]);
    assert!(aabbs_overlap(&a, &b));
    assert!(shape_shape_contact(&a, &b).is_none());
    assert!(shape_shape_contact(&b, &a).is_none());
}

#[test]
fn overlapping_tetrahedra_push_apart_along_the_body_diagonal() {
    let a = tetrahedron([0.0; 3]);
    let b = tetrahedron([0.2; 3]);
    let contact = shape_shape_contact(&a, &b).expect("the tetrahedra overlap");
    let diagonal = 1.0 / 3.0_f64.sqrt();
    // A moves out along -normal, so the normal points from A's bulk to B's: it
    // has a positive component along the body diagonal.
    assert!(
        dot(contact.normal, [diagonal; 3]) > 0.5,
        "normal {:?}",
        contact.normal
    );
    // Pushing B along (1, 1, 1) by 0.4 / sqrt(3) separates the tetrahedra, so
    // the minimum translation cannot be longer.
    assert!(
        contact.depth > 0.0 && contact.depth <= 0.4 * diagonal + 1e-9,
        "depth {}",
        contact.depth
    );
    assert_consistent(&contact, 1e-6);
    let mirrored = shape_shape_contact(&b, &a).expect("symmetric");
    assert!(dot(mirrored.normal, [diagonal; 3]) < -0.5);
    assert!((mirrored.depth - contact.depth).abs() < 1e-6);
}

/// Tiny deterministic generator for the randomised checks.
struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }

    fn cloud(&mut self, center: Vec3, spread: f64, count: usize) -> ShapeKind {
        ShapeKind::Convex {
            vertices: (0..count)
                .map(|_| {
                    [
                        center[0] + self.range(-spread, spread),
                        center[1] + self.range(-spread, spread),
                        center[2] + self.range(-spread, spread),
                    ]
                })
                .collect(),
        }
    }
}

#[test]
fn random_hull_pairs_give_consistent_antisymmetric_contacts() {
    let mut rng = Lcg(0x165_166);
    let (mut touching, mut apart) = (0, 0);
    for _ in 0..160 {
        let a = rng.cloud([0.0; 3], 1.0, 10);
        let offset = [
            rng.range(-1.3, 1.3),
            rng.range(-1.3, 1.3),
            rng.range(-1.3, 1.3),
        ];
        let b = rng.cloud(offset, 1.0, 10);
        let forward = shape_shape_contact(&a, &b);
        let backward = shape_shape_contact(&b, &a);
        match (forward, backward) {
            (Some(f), Some(r)) => {
                touching += 1;
                assert_consistent(&f, 1e-6);
                assert_consistent(&r, 1e-6);
                assert!(
                    (f.depth - r.depth).abs() < 1e-6,
                    "{} vs {}",
                    f.depth,
                    r.depth
                );
                assert!(
                    norm([
                        f.normal[0] + r.normal[0],
                        f.normal[1] + r.normal[1],
                        f.normal[2] + r.normal[2]
                    ]) < 1e-5,
                    "normals {:?} and {:?} are not opposite",
                    f.normal,
                    r.normal
                );
                // `normal * depth` is the minimum translation of B: moving B any
                // further along the normal separates the pair, moving it a little
                // less does not.
                let beyond = translated(&b, f.normal, f.depth + 1e-6);
                assert!(
                    shape_shape_contact(&a, &beyond).is_none(),
                    "B moved by depth {} along {:?} still touches A",
                    f.depth,
                    f.normal
                );
                if f.depth > 1e-4 {
                    let short = translated(&b, f.normal, f.depth - 1e-6);
                    assert!(
                        shape_shape_contact(&a, &short).is_some(),
                        "a shorter move along {:?} already separates the pair",
                        f.normal
                    );
                }
            }
            (None, None) => apart += 1,
            (f, r) => panic!("asymmetric result: {f:?} vs {r:?}"),
        }
    }
    assert!(
        touching > 20 && apart > 20,
        "{touching} touching, {apart} apart"
    );
}

/// `shape` (a convex vertex cloud) moved by `distance` along `direction`.
fn translated(shape: &ShapeKind, direction: Vec3, distance: f64) -> ShapeKind {
    let ShapeKind::Convex { vertices } = shape else {
        return shape.clone();
    };
    ShapeKind::Convex {
        vertices: vertices
            .iter()
            .map(|v| {
                [
                    v[0] + direction[0] * distance,
                    v[1] + direction[1] * distance,
                    v[2] + direction[2] * distance,
                ]
            })
            .collect(),
    }
}

/// The ground triangle used by the mesh tests: it lies in the plane `y = 0`, is
/// 20 units wide and its bounding sphere reaches more than 13 units above it.
fn ground() -> TriangleMesh {
    TriangleMesh::new(vec![
        [-10.0, 0.0, -10.0],
        [10.0, 0.0, -10.0],
        [0.0, 0.0, 10.0],
    ])
}

#[test]
fn a_sphere_over_a_flat_triangle_only_touches_the_actual_surface() {
    let hovering = ShapeKind::Sphere {
        center: [0.0, 0.6, 0.0],
        radius: 0.5,
    };
    // The old bounding-sphere stand-in made this a contact.
    assert!(convex_vs_mesh(&hovering, &ground(), 0.0).is_empty());

    let sunk = ShapeKind::Sphere {
        center: [1.0, 0.4, 0.5],
        radius: 0.5,
    };
    let contacts = convex_vs_mesh(&sunk, &ground(), 0.0);
    assert_eq!(contacts.len(), 1);
    let contact = &contacts[0];
    assert!(
        (contact.depth - 0.1).abs() < 1e-12,
        "depth {}",
        contact.depth
    );
    // From the sphere (A) towards the triangle (B): straight down.
    assert_close(contact.normal, [0.0, -1.0, 0.0], 1e-12);
    assert_close(contact.point_a, [1.0, -0.1, 0.5], 1e-12);
    assert_close(contact.point_b, [1.0, 0.0, 0.5], 1e-12);
    assert_consistent(contact, 1e-12);
}

#[test]
fn a_small_triangle_does_not_bulge_into_its_bounding_sphere() {
    // The bounding sphere of this unit right triangle has radius 0.82 around its
    // centroid; the sphere hovers at 0.6 above the plane with radius 0.5.
    let mesh = TriangleMesh::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
    let hovering = ShapeKind::Sphere {
        center: [0.3, 0.6, 0.3],
        radius: 0.5,
    };
    assert!(convex_vs_mesh(&hovering, &mesh, 0.0).is_empty());
    // Beside the triangle in its own plane, 0.3 from the edge x = 0.
    let beside = ShapeKind::Sphere {
        center: [-0.3, 0.0, 0.4],
        radius: 0.5,
    };
    let contacts = convex_vs_mesh(&beside, &mesh, 0.0);
    assert_eq!(contacts.len(), 1);
    assert!((contacts[0].depth - 0.2).abs() < 1e-12);
    assert_close(contacts[0].normal, [1.0, 0.0, 0.0], 1e-12);
}

#[test]
fn contacts_shallower_than_the_minimum_depth_are_dropped() {
    let sunk = ShapeKind::Sphere {
        center: [0.0, 0.4, 0.0],
        radius: 0.5,
    };
    assert_eq!(convex_vs_mesh(&sunk, &ground(), 0.05).len(), 1);
    assert!(convex_vs_mesh(&sunk, &ground(), 0.2).is_empty());
}

#[test]
fn a_sphere_over_a_seam_touches_both_triangles() {
    // A square floor split along the diagonal through the origin.
    let floor = TriangleMesh::new(vec![
        [-1.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [1.0, 0.0, 1.0],
        [-1.0, 0.0, -1.0],
        [1.0, 0.0, 1.0],
        [-1.0, 0.0, 1.0],
    ]);
    let sphere = ShapeKind::Sphere {
        center: [0.0, 0.4, 0.0],
        radius: 0.5,
    };
    let contacts = convex_vs_mesh(&sphere, &floor, 0.0);
    assert_eq!(contacts.len(), 2);
    for contact in &contacts {
        assert!((contact.depth - 0.1).abs() < 1e-12);
        assert_close(contact.normal, [0.0, -1.0, 0.0], 1e-12);
    }
}

#[test]
fn boxes_and_capsules_rest_on_the_mesh_surface() {
    let cube_on_ground = ShapeKind::Box {
        center: [0.5, 0.4, 0.0],
        half_extents: [0.5, 0.5, 0.5],
    };
    let contacts = convex_vs_mesh(&cube_on_ground, &ground(), 0.0);
    assert_eq!(contacts.len(), 1);
    assert!(
        (contacts[0].depth - 0.1).abs() < 1e-9,
        "{}",
        contacts[0].depth
    );
    assert!(contacts[0].normal[1] < -0.999, "{:?}", contacts[0].normal);
    let lifted = ShapeKind::Box {
        center: [0.5, 0.6, 0.0],
        half_extents: [0.5, 0.5, 0.5],
    };
    assert!(convex_vs_mesh(&lifted, &ground(), 0.0).is_empty());

    let lying = ShapeKind::Capsule {
        p0: [-1.0, 0.4, 0.0],
        p1: [1.0, 0.4, 0.0],
        radius: 0.5,
    };
    let contacts = convex_vs_mesh(&lying, &ground(), 0.0);
    assert_eq!(contacts.len(), 1);
    assert!((contacts[0].depth - 0.1).abs() < 1e-12);
    assert_close(contacts[0].normal, [0.0, -1.0, 0.0], 1e-12);
}

#[test]
fn a_turned_convex_cube_over_the_mesh_touches_with_its_lowest_corner() {
    // A cube standing on an edge (turned 45 degrees about z) sinks 0.05 into
    // the ground with that edge.
    let (s, c) = std::f64::consts::FRAC_PI_4.sin_cos();
    // Half the face diagonal: how far the turned cube's edge is below its centre.
    let half_diagonal = std::f64::consts::FRAC_1_SQRT_2;
    let mut vertices = Vec::new();
    for sx in [-0.5, 0.5] {
        for sy in [-0.5, 0.5] {
            for sz in [-0.5, 0.5] {
                vertices.push([c * sx - s * sy, s * sx + c * sy + half_diagonal - 0.05, sz]);
            }
        }
    }
    let lowest = vertices.iter().map(|v| v[1]).fold(f64::INFINITY, f64::min);
    let hull = ShapeKind::Convex { vertices };
    let contacts = convex_vs_mesh(&hull, &ground(), 0.0);
    assert_eq!(contacts.len(), 1);
    assert!(
        (contacts[0].depth + lowest).abs() < 1e-6,
        "depth {} vs lowest corner {lowest}",
        contacts[0].depth
    );
    assert!(contacts[0].normal[1] < -0.999);
}

#[test]
fn zero_area_triangles_and_distant_meshes_never_collide() {
    let sliver = TriangleMesh::new(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let sphere = ShapeKind::Sphere {
        center: [1.0, 0.0, 0.0],
        radius: 1.0,
    };
    assert!(convex_vs_mesh(&sphere, &sliver, 0.0).is_empty());
    let far = ShapeKind::Sphere {
        center: [100.0, 100.0, 100.0],
        radius: 1.0,
    };
    assert!(convex_vs_mesh(&far, &ground(), 0.0).is_empty());
}
