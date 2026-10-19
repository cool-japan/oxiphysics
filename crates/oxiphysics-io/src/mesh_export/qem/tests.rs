// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for the quadric-error-metric simplifier.

use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap};

use super::*;
use crate::mesh_export::{MeshData, MeshSimplification};
use crate::mesh_orientation::count_inverted_faces;

// ---------------------------------------------------------------------------
// Mesh builders and measurements
// ---------------------------------------------------------------------------

fn unit(p: Vec3) -> Vec3 {
    let n = norm(p);
    [p[0] / n, p[1] / n, p[2] / n]
}

/// Icosphere of radius 1 with `20 * 4^level` outward-facing triangles.
fn icosphere(level: usize) -> MeshData {
    let t = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut mesh = MeshData::new();
    for p in [
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ] {
        mesh.add_vertex(unit(p));
    }
    let faces = [
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for [a, b, c] in faces {
        let (pa, pb, pc) = (mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]);
        let normal = cross(sub(pb, pa), sub(pc, pa));
        if dot(normal, pa) > 0.0 {
            mesh.add_face(a, b, c);
        } else {
            mesh.add_face(a, c, b);
        }
    }
    for _ in 0..level {
        let mut midpoints: HashMap<(usize, usize), usize> = HashMap::new();
        let faces = std::mem::take(&mut mesh.faces);
        for [a, b, c] in faces {
            let mut middle = |mesh: &mut MeshData, u: usize, v: usize| -> usize {
                *midpoints.entry((u.min(v), u.max(v))).or_insert_with(|| {
                    let point = lerp(mesh.vertices[u], mesh.vertices[v], 0.5);
                    mesh.add_vertex(unit(point))
                })
            };
            let ab = middle(&mut mesh, a, b);
            let bc = middle(&mut mesh, b, c);
            let ca = middle(&mut mesh, c, a);
            mesh.add_face(a, ab, ca);
            mesh.add_face(b, bc, ab);
            mesh.add_face(c, ca, bc);
            mesh.add_face(ab, bc, ca);
        }
    }
    mesh
}

/// `nx x ny` cells in the z = 0 plane, two counter-clockwise triangles per cell.
fn grid(nx: usize, ny: usize) -> MeshData {
    let mut mesh = MeshData::new();
    for j in 0..=ny {
        for i in 0..=nx {
            mesh.add_vertex([i as f64, j as f64, 0.0]);
        }
    }
    let stride = nx + 1;
    for j in 0..ny {
        for i in 0..nx {
            let a = j * stride + i;
            let (b, c, d) = (a + 1, a + stride, a + stride + 1);
            mesh.add_face(a, b, d);
            mesh.add_face(a, d, c);
        }
    }
    mesh
}

fn tetrahedron() -> MeshData {
    let mut mesh = MeshData::new();
    mesh.add_vertex([1.0, 1.0, 1.0]);
    mesh.add_vertex([1.0, -1.0, -1.0]);
    mesh.add_vertex([-1.0, 1.0, -1.0]);
    mesh.add_vertex([-1.0, -1.0, 1.0]);
    for [a, b, c] in [[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]] {
        mesh.add_face(a, b, c);
    }
    mesh
}

/// Octahedron with the poles at vertices 4 (top) and 5 (bottom).
fn octahedron() -> MeshData {
    let mut mesh = MeshData::new();
    for p in [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ] {
        mesh.add_vertex(p);
    }
    for i in 0..4 {
        let j = (i + 1) % 4;
        mesh.add_face(i, j, 4);
        mesh.add_face(j, i, 5);
    }
    mesh
}

/// The unit cube whose six faces are `n x n` grids of outward-facing triangles; the
/// vertices on the cube edges are shared between faces.
fn subdivided_cube(n: usize) -> MeshData {
    // (origin, u axis, v axis) of every face; `u x v` points out of the cube.
    let faces: [(Vec3, Vec3, Vec3); 6] = [
        ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]),
        ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]),
        ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
    ];
    let mut mesh = MeshData::new();
    let mut known: HashMap<[u64; 3], usize> = HashMap::new();
    for (origin, u, v) in faces {
        let mut vertex = |mesh: &mut MeshData, i: usize, j: usize| -> usize {
            let (s, t) = (i as f64 / n as f64, j as f64 / n as f64);
            let p = [
                origin[0] + s * u[0] + t * v[0],
                origin[1] + s * u[1] + t * v[1],
                origin[2] + s * u[2] + t * v[2],
            ];
            *known
                .entry(p.map(f64::to_bits))
                .or_insert_with(|| mesh.add_vertex(p))
        };
        for j in 0..n {
            for i in 0..n {
                let a = vertex(&mut mesh, i, j);
                let b = vertex(&mut mesh, i + 1, j);
                let c = vertex(&mut mesh, i + 1, j + 1);
                let d = vertex(&mut mesh, i, j + 1);
                mesh.add_face(a, b, c);
                mesh.add_face(a, c, d);
            }
        }
    }
    mesh
}

fn torus(segments: usize, sides: usize, major: f64, minor: f64) -> MeshData {
    let tau = 2.0 * std::f64::consts::PI;
    let mut mesh = MeshData::new();
    for i in 0..segments {
        let u = tau * i as f64 / segments as f64;
        for j in 0..sides {
            let v = tau * j as f64 / sides as f64;
            let ring = major + minor * v.cos();
            mesh.add_vertex([ring * u.cos(), ring * u.sin(), minor * v.sin()]);
        }
    }
    for i in 0..segments {
        for j in 0..sides {
            let a = i * sides + j;
            let b = ((i + 1) % segments) * sides + j;
            let c = ((i + 1) % segments) * sides + (j + 1) % sides;
            let d = i * sides + (j + 1) % sides;
            mesh.add_face(a, b, c);
            mesh.add_face(a, c, d);
        }
    }
    mesh
}

/// A flat ring of `n` outer and `n` inner vertices, two triangles per segment.
fn annulus(n: usize) -> MeshData {
    let tau = 2.0 * std::f64::consts::PI;
    let mut mesh = MeshData::new();
    for radius in [2.0, 1.0] {
        for i in 0..n {
            let angle = tau * i as f64 / n as f64;
            mesh.add_vertex([radius * angle.cos(), radius * angle.sin(), 0.0]);
        }
    }
    for i in 0..n {
        let next = (i + 1) % n;
        mesh.add_face(i, next, n + next);
        mesh.add_face(i, n + next, n + i);
    }
    mesh
}

fn simplified(mesh: &MeshData, target: usize) -> MeshData {
    MeshSimplification::new().simplify(mesh, target)
}

/// How many triangles use each undirected edge.
fn edge_uses(mesh: &MeshData) -> BTreeMap<(usize, usize), usize> {
    let mut uses = BTreeMap::new();
    for face in &mesh.faces {
        for (slot, &a) in face.iter().enumerate() {
            let b = face[(slot + 1) % 3];
            *uses.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    uses
}

fn is_closed_manifold(mesh: &MeshData) -> bool {
    edge_uses(mesh).values().all(|&count| count == 2)
}

fn boundary_edges(mesh: &MeshData) -> Vec<(usize, usize)> {
    edge_uses(mesh)
        .into_iter()
        .filter(|&(_, count)| count == 1)
        .map(|(edge, _)| edge)
        .collect()
}

/// `V - E + F` over the vertices that triangles use.
fn euler_characteristic(mesh: &MeshData) -> i64 {
    let used: BTreeSet<usize> = mesh.faces.iter().flatten().copied().collect();
    used.len() as i64 - edge_uses(mesh).len() as i64 + mesh.faces.len() as i64
}

fn inverted_faces(mesh: &MeshData) -> usize {
    let faces: Vec<[u32; 3]> = mesh.faces.iter().map(|f| f.map(|v| v as u32)).collect();
    count_inverted_faces(&mesh.vertices, &faces)
}

fn enclosed_volume(mesh: &MeshData) -> f64 {
    mesh.faces
        .iter()
        .map(|f| {
            let [a, b, c] = f.map(|v| mesh.vertices[v]);
            dot(a, cross(b, c)) / 6.0
        })
        .sum()
}

/// Signed area of the triangles in the z = 0 plane; positive for counter-clockwise ones.
fn signed_areas_xy(mesh: &MeshData) -> Vec<f64> {
    mesh.faces
        .iter()
        .map(|f| {
            let [a, b, c] = f.map(|v| mesh.vertices[v]);
            0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
        })
        .collect()
}

/// Every face has three distinct, in-range corners and every vertex is used.
fn assert_well_formed(mesh: &MeshData) {
    let mut used = vec![false; mesh.vertices.len()];
    for face in &mesh.faces {
        assert!(
            face.iter().all(|&v| v < mesh.vertices.len()),
            "face {face:?} indexes past {} vertices",
            mesh.vertices.len()
        );
        assert!(
            face[0] != face[1] && face[1] != face[2] && face[0] != face[2],
            "degenerate face {face:?}"
        );
        for &v in face {
            used[v] = true;
        }
    }
    assert!(used.iter().all(|&u| u), "unreferenced output vertex");
    for v in &mesh.vertices {
        assert!(v.iter().all(|c| c.is_finite()), "non-finite vertex {v:?}");
    }
}

fn close(a: Vec3, b: Vec3, tolerance: f64) -> bool {
    norm(sub(a, b)) <= tolerance
}

// ---------------------------------------------------------------------------
// Quadrics
// ---------------------------------------------------------------------------

/// Quadric of the plane with the given normal through `point`.
fn plane_quadric(weight: f64, normal: Vec3, point: Vec3) -> Quadric {
    let normal = unit(normal);
    Quadric::from_plane(weight, normal, -dot(normal, point))
}

fn total(quadrics: &[Quadric]) -> Quadric {
    let mut sum = Quadric::default();
    for q in quadrics {
        sum.add(q);
    }
    sum
}

#[test]
fn three_orthogonal_planes_meet_at_their_common_point() {
    let q = total(&[
        plane_quadric(1.0, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        plane_quadric(1.0, [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]),
        plane_quadric(1.0, [0.0, 0.0, 1.0], [0.0, 0.0, 3.0]),
    ]);
    let optimum = q.minimizer().expect("three planes meet in a point");
    assert!(close(optimum, [1.0, 2.0, 3.0], 1e-12), "{optimum:?}");
    assert!(q.error_at([1.0, 2.0, 3.0]).abs() < 1e-12);
    assert!((q.error_at([0.0, 0.0, 0.0]) - 14.0).abs() < 1e-12);
    assert!((q.error_at([2.0, 2.0, 3.0]) - 1.0).abs() < 1e-12);
}

#[test]
fn skewed_planes_meet_at_their_intersection() {
    let point = [0.3, -1.2, 2.5];
    let q = total(&[
        plane_quadric(2.0, [1.0, 1.0, 0.0], point),
        plane_quadric(0.5, [0.0, 1.0, 1.0], point),
        plane_quadric(3.0, [1.0, 0.0, 2.0], point),
    ]);
    let optimum = q.minimizer().expect("three independent planes");
    assert!(close(optimum, point, 1e-10), "{optimum:?}");
}

#[test]
fn heavier_planes_pull_the_optimum_towards_themselves() {
    // x = 0 with weight 1 against x = 1 with weight 3 gives x = 3 / 4.
    let q = total(&[
        plane_quadric(1.0, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
        plane_quadric(3.0, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        plane_quadric(1.0, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
        plane_quadric(1.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
    ]);
    let optimum = q.minimizer().expect("full rank");
    assert!(close(optimum, [0.75, 0.0, 0.0], 1e-12), "{optimum:?}");
}

#[test]
fn quadrics_add_coefficient_wise() {
    let a = plane_quadric(1.5, [1.0, 2.0, 3.0], [0.0, 1.0, 0.0]);
    let b = plane_quadric(0.25, [0.0, 1.0, -1.0], [2.0, 0.0, 0.0]);
    let sum = a.sum(&b);
    for ((s, x), y) in sum.q.iter().zip(&a.q).zip(&b.q) {
        assert_eq!(*s, x + y);
    }
    let probe = [0.4, -0.7, 1.9];
    assert!((sum.error_at(probe) - (a.error_at(probe) + b.error_at(probe))).abs() < 1e-12);
}

#[test]
fn the_error_is_the_weighted_squared_distance_to_the_plane() {
    let q = plane_quadric(2.0, [0.0, 0.0, 1.0], [5.0, 5.0, 4.0]);
    assert!((q.error_at([0.0, 0.0, 6.5]) - 2.0 * 2.5 * 2.5).abs() < 1e-12);
    assert!(q.error_at([123.0, -45.0, 4.0]).abs() < 1e-9);
}

#[test]
fn rank_deficient_quadrics_have_no_minimizer() {
    let two = total(&[
        plane_quadric(1.0, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        plane_quadric(1.0, [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]),
    ]);
    assert_eq!(two.minimizer(), None);
    let one = plane_quadric(1.0, [0.0, 0.0, 1.0], [0.0, 0.0, 3.0]);
    assert_eq!(one.minimizer(), None);
    assert_eq!(Quadric::default().minimizer(), None);
}

#[test]
fn nearly_parallel_planes_count_as_singular() {
    let tilted = |angle: f64| {
        total(&[
            plane_quadric(1.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
            plane_quadric(1.0, [angle.sin(), 0.0, angle.cos()], [0.0, 0.0, 0.01]),
            plane_quadric(1.0, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
        ])
    };
    // sin^2 of the angle is the determinant: 1e-14 is rounding noise, 1e-6 is a solution.
    assert_eq!(tilted(1e-7).minimizer(), None);
    let optimum = tilted(1e-3).minimizer().expect("well enough conditioned");
    assert!(
        (optimum[0] - 10.0).abs() < 0.1 && optimum[1].abs() < 1e-9 && optimum[2].abs() < 1e-3,
        "{optimum:?}"
    );
}

#[test]
fn non_finite_quadrics_have_no_minimizer() {
    let mut q = total(&[
        plane_quadric(1.0, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        plane_quadric(1.0, [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]),
        plane_quadric(1.0, [0.0, 0.0, 1.0], [0.0, 0.0, 3.0]),
    ]);
    q.q[3] = f64::NAN;
    assert_eq!(q.minimizer(), None);
    q.q[0] = f64::NAN;
    assert_eq!(q.minimizer(), None);
    assert_eq!(q.minimum_on_segment([0.0; 3], [1.0, 0.0, 0.0]), None);
}

#[test]
fn the_segment_minimum_is_exact_and_clamped_to_the_edge() {
    let q = total(&[
        plane_quadric(1.0, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        plane_quadric(1.0, [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]),
    ]);
    // error(t) = (2t - 1)^2 on the segment from (0, 2, 0) to (2, 2, 0).
    assert_eq!(
        q.minimum_on_segment([0.0, 2.0, 0.0], [2.0, 2.0, 0.0]),
        Some(0.5)
    );
    // The minimum x = 1 is before the segment starts at x = 3.
    assert_eq!(
        q.minimum_on_segment([3.0, 2.0, 0.0], [5.0, 2.0, 0.0]),
        Some(0.0)
    );
    assert_eq!(
        q.minimum_on_segment([-4.0, 2.0, 0.0], [-2.0, 2.0, 0.0]),
        Some(1.0)
    );
    // The segment runs along the line where both planes meet: the error is constant.
    assert_eq!(q.minimum_on_segment([1.0, 2.0, 0.0], [1.0, 2.0, 7.0]), None);
    // Off-axis segment: the minimum of (4t - 1)^2 + (4t - 2)^2 is at t = 3 / 8.
    let t = q
        .minimum_on_segment([0.0, 0.0, 0.0], [4.0, 4.0, 0.0])
        .expect("curved along the segment");
    assert!((t - 0.375).abs() < 1e-12, "{t}");
}

// ---------------------------------------------------------------------------
// Placement of the merged vertex
// ---------------------------------------------------------------------------

/// A simplifier holding one triangle whose first two vertices and quadrics are replaced,
/// to cost the edge between them.
fn edge_between(from: Vec3, to: Vec3, quadric_a: Quadric, quadric_b: Quadric) -> Simplifier {
    let mut simplifier =
        Simplifier::new(vec![from, to, [0.0, 0.0, 100.0]], [0.0; 3], vec![[0, 1, 2]]);
    simplifier.quadrics[0] = quadric_a;
    simplifier.quadrics[1] = quadric_b;
    simplifier
}

#[test]
fn placement_uses_the_optimum_when_it_is_well_defined() {
    let simplifier = edge_between(
        [0.5, 1.5, 2.5],
        [1.5, 2.5, 3.5],
        total(&[
            plane_quadric(1.0, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            plane_quadric(1.0, [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]),
        ]),
        plane_quadric(1.0, [0.0, 0.0, 1.0], [0.0, 0.0, 3.0]),
    );
    let (target, cost) = simplifier.best_target(0, 1);
    assert!(close(target, [1.0, 2.0, 3.0], 1e-12), "{target:?}");
    assert!(cost.abs() < 1e-12);
}

#[test]
fn placement_falls_back_to_the_minimum_on_the_edge_when_the_system_is_singular() {
    let simplifier = edge_between(
        [0.0, 2.0, 0.0],
        [2.0, 2.0, 0.0],
        plane_quadric(1.0, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        plane_quadric(1.0, [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]),
    );
    assert_eq!(
        simplifier.quadrics[0]
            .sum(&simplifier.quadrics[1])
            .minimizer(),
        None
    );
    let (target, cost) = simplifier.best_target(0, 1);
    assert!(close(target, [1.0, 2.0, 0.0], 1e-12), "{target:?}");
    assert!(cost.abs() < 1e-12);
}

#[test]
fn placement_on_a_flat_edge_keeps_the_first_endpoint() {
    let plane = plane_quadric(1.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
    let simplifier = edge_between([0.25, 0.5, 0.0], [3.0, -2.0, 0.0], plane, plane);
    let (target, cost) = simplifier.best_target(0, 1);
    assert_eq!(target, [0.25, 0.5, 0.0]);
    assert_eq!(cost, 0.0);
}

#[test]
fn a_distant_optimum_is_replaced_by_a_point_on_the_edge() {
    let quadric = total(&[
        plane_quadric(1.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
        plane_quadric(1.0, [1e-3_f64.sin(), 0.0, 1e-3_f64.cos()], [0.0, 0.0, 0.01]),
        plane_quadric(1.0, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
    ]);
    let optimum = quadric.minimizer().expect("conditioned well enough");
    assert!((optimum[0] - 10.0).abs() < 0.1, "{optimum:?}");
    let simplifier = edge_between(
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        quadric,
        Quadric::default(),
    );
    let (target, cost) = simplifier.best_target(0, 1);
    assert!(
        norm(sub(target, [0.5, 0.0, 0.0])) <= 0.5 + 1e-12,
        "target {target:?} leaves the edge"
    );
    assert!(close(target, [1.0, 0.0, 0.0], 1e-9), "{target:?}");
    assert!((cost - quadric.error_at(target)).abs() < 1e-15);
}

#[test]
fn points_on_a_segment_snap_to_its_endpoints() {
    let (from, to) = ([0.1, -2.7, 3.3], [4.9, 0.3, -1.1]);
    assert_eq!(point_on_segment(from, to, 0.0), from);
    assert_eq!(point_on_segment(from, to, 1e-13), from);
    assert_eq!(point_on_segment(from, to, 1.0), to);
    assert_eq!(point_on_segment(from, to, 1.0 - 1e-13), to);
    assert_eq!(point_on_segment(from, to, 0.5), lerp(from, to, 0.5));
    assert_eq!(point_on_segment(from, to, 0.25), lerp(from, to, 0.25));
}

#[test]
fn an_optimum_that_is_an_endpoint_up_to_rounding_is_snapped_to_it() {
    let (from, to) = ([0.1, 0.2, 0.3], [1.1, 0.2, 0.3]);
    assert_eq!(
        snap_to_endpoint(from, to, [0.1 + 1e-14, 0.2, 0.3 - 1e-15], 1.0),
        from
    );
    assert_eq!(
        snap_to_endpoint(from, to, [1.1 - 1e-14, 0.2 + 1e-15, 0.3], 1.0),
        to
    );
    let inside = [0.6, 0.2, 0.3];
    assert_eq!(snap_to_endpoint(from, to, inside, 1.0), inside);
    let beside = [0.1, 0.2 + 1e-6, 0.3];
    assert_eq!(snap_to_endpoint(from, to, beside, 1.0), beside);
}

#[test]
fn a_corner_where_three_planes_meet_keeps_its_exact_coordinates() {
    // The solution of the 3x3 system is the corner only up to rounding; the placement
    // must still return the corner vertex itself.
    let corner = [0.3, -1.2, 2.5];
    let quadric = total(&[
        plane_quadric(2.0, [1.0, 1.0, 0.0], corner),
        plane_quadric(0.5, [0.0, 1.0, 1.0], corner),
        plane_quadric(3.0, [1.0, 0.0, 2.0], corner),
    ]);
    let found = quadric.minimizer().expect("three independent planes");
    assert!(close(found, corner, 1e-10));
    let simplifier = edge_between(corner, [1.3, -1.2, 2.5], quadric, Quadric::default());
    let (target, _) = simplifier.best_target(0, 1);
    assert_eq!(target, corner);
}

// ---------------------------------------------------------------------------
// Priority queue
// ---------------------------------------------------------------------------

fn candidate(cost: f64, length_sq: f64, a: usize, b: usize) -> Candidate {
    Candidate {
        cost,
        length_sq,
        a,
        b,
        version_a: 0,
        version_b: 0,
        target: [0.0; 3],
    }
}

#[test]
fn costs_are_never_nan_or_negative() {
    assert_eq!(finite_cost(f64::NAN), f64::INFINITY);
    assert_eq!(finite_cost(-f64::NAN), f64::INFINITY);
    assert_eq!(finite_cost(f64::INFINITY), f64::INFINITY);
    assert_eq!(finite_cost(2.5), 2.5);
    assert_eq!(finite_cost(-1e-30).to_bits(), 0.0_f64.to_bits());
    assert_eq!(finite_cost(-0.0).to_bits(), 0.0_f64.to_bits());
}

#[test]
fn the_heap_pops_the_cheapest_then_shortest_then_lowest_edge() {
    let mut heap = BinaryHeap::new();
    heap.extend([
        candidate(finite_cost(f64::NAN), 1.0, 0, 2),
        candidate(3.0, 1.0, 0, 1),
        candidate(1.0, 4.0, 5, 6),
        candidate(1.0, 4.0, 2, 3),
        candidate(1.0, 1.0, 9, 10),
        candidate(finite_cost(-0.0), 2.0, 7, 8),
    ]);
    let order: Vec<(usize, usize)> = std::iter::from_fn(|| heap.pop())
        .map(|c| (c.a, c.b))
        .collect();
    assert_eq!(order, vec![(7, 8), (9, 10), (2, 3), (5, 6), (0, 1), (0, 2)]);
}

// ---------------------------------------------------------------------------
// Orientation test
// ---------------------------------------------------------------------------

#[test]
fn small_moves_keep_the_orientation() {
    let before = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let after = [[0.1, 0.1, 0.05], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert!(keeps_orientation(before, after));
}

#[test]
fn a_triangle_turned_over_does_not_keep_the_orientation() {
    let before = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let turned = [[0.5, 2.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert!(!keeps_orientation(before, turned));
}

#[test]
fn a_right_angle_turn_is_not_allowed() {
    let before = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let upright = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    assert!(!keeps_orientation(before, upright));
    let tilted = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.2]];
    assert!(keeps_orientation(before, tilted));
}

#[test]
fn collapsing_a_triangle_to_a_line_is_not_allowed() {
    let before = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let line = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    assert!(!keeps_orientation(before, line));
    let point = [[1.0, 1.0, 1.0]; 3];
    assert!(!keeps_orientation(before, point));
}

#[test]
fn a_degenerate_triangle_has_no_orientation_to_keep() {
    let flat = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let after = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]];
    assert!(keeps_orientation(flat, after));
}

// ---------------------------------------------------------------------------
// Collapse bookkeeping
// ---------------------------------------------------------------------------

/// The adjacency of `simplifier` agrees with its triangle list.
fn assert_adjacency_is_consistent(simplifier: &Simplifier) {
    let alive_faces = simplifier.face_alive.iter().filter(|&&alive| alive).count();
    assert_eq!(simplifier.live_faces, alive_faces);
    for (index, face) in simplifier.faces.iter().enumerate() {
        if !simplifier.face_alive[index] {
            continue;
        }
        assert!(
            face[0] != face[1] && face[1] != face[2] && face[0] != face[2],
            "live face {index} is degenerate: {face:?}"
        );
        for &v in face {
            assert!(
                simplifier.vertex_alive[v],
                "live face {index} uses dead vertex {v}"
            );
            assert_eq!(
                simplifier.vertex_faces[v]
                    .iter()
                    .filter(|&&f| f == index)
                    .count(),
                1,
                "face {index} is not listed exactly once at vertex {v}"
            );
        }
    }
    for (vertex, faces) in simplifier.vertex_faces.iter().enumerate() {
        if !simplifier.vertex_alive[vertex] {
            assert!(faces.is_empty(), "dead vertex {vertex} still has faces");
        }
        for &f in faces {
            assert!(
                simplifier.face_alive[f],
                "vertex {vertex} lists dead face {f}"
            );
            assert!(
                simplifier.faces[f].contains(&vertex),
                "vertex {vertex} lists face {f} that does not use it"
            );
        }
    }
}

fn simplifier_for(mesh: &MeshData) -> Simplifier {
    Simplifier::new(mesh.vertices.clone(), [0.0; 3], mesh.faces.clone())
}

#[test]
fn the_survivor_takes_over_the_sum_of_both_quadrics() {
    let mesh = grid(2, 2);
    let mut simplifier = simplifier_for(&mesh);
    let (survivor, absorbed) = (4, 5);
    let before_survivor = simplifier.quadrics[survivor];
    let before_absorbed = simplifier.quadrics[absorbed];
    let neighbour_before = simplifier.quadrics[1];
    simplifier.collapse(survivor, absorbed, [1.25, 1.0, 0.0]);
    for ((after, x), y) in simplifier.quadrics[survivor]
        .q
        .iter()
        .zip(&before_survivor.q)
        .zip(&before_absorbed.q)
    {
        assert_eq!(*after, x + y);
    }
    assert_eq!(
        simplifier.quadrics[1], neighbour_before,
        "bystanders keep their quadric"
    );
    assert_eq!(simplifier.positions[survivor], [1.25, 1.0, 0.0]);
    assert_eq!(simplifier.absolute[survivor], [1.25, 1.0, 0.0]);
    assert_eq!(simplifier.absolute[absorbed], mesh.vertices[absorbed]);
    assert!(!simplifier.vertex_alive[absorbed]);
    assert!(simplifier.vertex_faces[absorbed].is_empty());
    assert_eq!(simplifier.versions[survivor], 1);
    assert_eq!(simplifier.versions[0], 0);
    // The edge (4, 5) is interior, so its two triangles vanish.
    assert_eq!(simplifier.live_faces, 6);
    assert_adjacency_is_consistent(&simplifier);
}

#[test]
fn a_boundary_edge_collapse_removes_a_single_triangle() {
    let mesh = grid(2, 2);
    let mut simplifier = simplifier_for(&mesh);
    simplifier.collapse(1, 2, [1.5, 0.0, 0.0]);
    assert_eq!(simplifier.live_faces, 7);
    assert_adjacency_is_consistent(&simplifier);
}

#[test]
fn candidates_that_mention_a_collapsed_vertex_become_stale() {
    let mesh = grid(3, 3);
    let mut simplifier = simplifier_for(&mesh);
    simplifier.seed_all_edges();
    let candidates: Vec<Candidate> = simplifier.heap.clone().into_vec();
    assert!(candidates.iter().all(|c| simplifier.is_current(c)));
    simplifier.collapse(5, 6, [1.5, 1.0, 0.0]);
    let (touching, bystanders): (Vec<&Candidate>, Vec<&Candidate>) = candidates
        .iter()
        .partition(|c| [c.a, c.b].iter().any(|&v| v == 5 || v == 6));
    assert!(!touching.is_empty() && !bystanders.is_empty());
    assert!(touching.iter().all(|c| !simplifier.is_current(c)));
    assert!(bystanders.iter().all(|c| simplifier.is_current(c)));
}

#[test]
fn only_the_edges_of_the_survivor_are_costed_again() {
    let mesh = grid(3, 3);
    let mut simplifier = simplifier_for(&mesh);
    simplifier.seed_all_edges();
    simplifier.collapse(5, 6, [1.5, 1.0, 0.0]);
    let before = simplifier.heap.len();
    simplifier.refresh_edges(5);
    let neighbours: Vec<usize> = simplifier
        .neighbour_counts(5)
        .into_iter()
        .map(|(v, _)| v)
        .collect();
    assert_eq!(simplifier.heap.len() - before, neighbours.len());
    let fresh: Vec<Candidate> = simplifier
        .heap
        .iter()
        .copied()
        .filter(|c| simplifier.is_current(c) && (c.a == 5 || c.b == 5))
        .collect();
    assert_eq!(fresh.len(), neighbours.len());
    for c in fresh {
        let other = if c.a == 5 { c.b } else { c.a };
        assert!(neighbours.contains(&other));
        assert_eq!(c.version_a, simplifier.versions[c.a]);
        assert_eq!(c.version_b, simplifier.versions[c.b]);
    }
}

#[test]
fn refused_edges_are_reseeded_only_while_they_still_exist() {
    let mesh = grid(2, 2);
    let mut simplifier = simplifier_for(&mesh);
    simplifier.deferred = vec![(4, 5), (4, 5), (0, 8), (1, 3)];
    simplifier.collapse(0, 1, [0.5, 0.0, 0.0]);
    simplifier.deferred.push((0, 1));
    simplifier.reseed_deferred();
    // (4, 5) is an edge (once, whatever the duplicates), (0, 8) and (1, 3) were never
    // edges, and vertex 1 is gone.
    let pushed: Vec<(usize, usize)> = simplifier.heap.iter().map(|c| (c.a, c.b)).collect();
    assert_eq!(pushed, vec![(4, 5)]);
    assert!(simplifier.deferred.is_empty());
}

// ---------------------------------------------------------------------------
// Topology rules
// ---------------------------------------------------------------------------

#[test]
fn a_tetrahedron_cannot_lose_a_vertex() {
    let simplifier = simplifier_for(&tetrahedron());
    for a in 0..4 {
        for b in a + 1..4 {
            assert!(
                !simplifier.preserves_topology(a, b),
                "edge ({a}, {b}) would fold the tetrahedron flat"
            );
        }
    }
}

#[test]
fn an_octahedron_edge_may_collapse() {
    let simplifier = simplifier_for(&octahedron());
    assert!(simplifier.preserves_topology(0, 4));
    assert!(!simplifier.preserves_topology(4, 5), "not an edge");
}

#[test]
fn an_interior_edge_between_two_boundary_vertices_is_not_collapsed() {
    // Two triangles forming a quad: the diagonal joins two boundary vertices.
    let mut mesh = MeshData::new();
    for p in [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ] {
        mesh.add_vertex(p);
    }
    mesh.add_face(0, 1, 2);
    mesh.add_face(0, 2, 3);
    let simplifier = simplifier_for(&mesh);
    assert!(!simplifier.preserves_topology(0, 2));
    assert!(simplifier.preserves_topology(0, 1));
    let out = simplified(&mesh, 1);
    assert_eq!(out.faces.len(), 1);
    assert_well_formed(&out);
}

#[test]
fn a_vertex_adjacent_to_both_ends_beyond_the_edge_triangles_blocks_the_collapse() {
    // Vertices 0 and 1 share the edge (0, 1), whose only triangle has the apex 2, but
    // 3 is adjacent to both of them as well: collapsing would glue two triangles.
    let mut mesh = MeshData::new();
    for p in [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0, -1.0, 0.0],
        [1.0, -3.0, 0.0],
        [-2.0, -3.0, 0.0],
    ] {
        mesh.add_vertex(p);
    }
    mesh.add_face(0, 1, 2);
    mesh.add_face(0, 3, 4);
    mesh.add_face(1, 4, 3);
    let simplifier = simplifier_for(&mesh);
    assert!(
        simplifier.neighbour_counts(0).iter().any(|&(v, _)| v == 3)
            && simplifier.neighbour_counts(1).iter().any(|&(v, _)| v == 3)
    );
    assert!(!simplifier.preserves_topology(0, 1));
}

#[test]
fn non_manifold_edges_are_left_alone() {
    // Three triangles share the edge (0, 1).
    let mut mesh = MeshData::new();
    for p in [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, -1.0, 0.0],
        [0.5, 0.0, 1.0],
    ] {
        mesh.add_vertex(p);
    }
    for apex in [2, 3, 4] {
        mesh.add_face(0, 1, apex);
    }
    let out = simplified(&mesh, 1);
    assert_eq!(out.faces.len(), 3);
    assert_eq!(out.vertices, mesh.vertices);
}

#[test]
fn an_isolated_triangle_cannot_be_collapsed() {
    let mut mesh = MeshData::new();
    for p in [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]] {
        mesh.add_vertex(p);
    }
    mesh.add_face(0, 1, 2);
    let simplifier = simplifier_for(&mesh);
    for (a, b) in [(0, 1), (1, 2), (0, 2)] {
        assert!(!simplifier.preserves_topology(a, b), "edge ({a}, {b})");
    }
    let out = simplified(&mesh, 0);
    assert_eq!(out.faces, vec![[0, 1, 2]]);
}

/// Two triangles sharing an edge plus a third one that touches them only at the vertices
/// 0 and 4: the three form a ring around a hole.
fn ring_of_triangles() -> MeshData {
    let mut mesh = MeshData::new();
    for p in [
        [0.0, 0.0, 0.0],
        [1.0, 0.5, 0.0],
        [-0.5, 2.0, 0.0],
        [-1.0, 1.0, 0.0],
        [2.0, 2.0, 0.0],
    ] {
        mesh.add_vertex(p);
    }
    for face in [[3, 1, 0], [1, 4, 3], [4, 2, 0]] {
        mesh.add_face(face[0], face[1], face[2]);
    }
    mesh
}

#[test]
fn an_ear_that_touches_the_surface_at_vertices_only_is_not_collapsed() {
    let mesh = ring_of_triangles();
    assert_eq!(euler_characteristic(&mesh), 0);
    let simplifier = simplifier_for(&mesh);
    // The triangle (4, 2, 0) has no neighbour across an edge: collapsing one of its edges
    // would delete it and open the ring.
    for (a, b) in [(2, 4), (0, 2), (0, 4)] {
        assert!(!simplifier.preserves_topology(a, b), "edge ({a}, {b})");
    }
    let out = simplified(&mesh, 0);
    assert_well_formed(&out);
    assert_eq!(euler_characteristic(&out), 0);
}

/// Number of connected pieces, joined through shared vertices.
fn component_count(mesh: &MeshData) -> usize {
    let mut parent: Vec<usize> = (0..mesh.vertices.len()).collect();
    fn root(parent: &mut [usize], mut v: usize) -> usize {
        while parent[v] != v {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        v
    }
    for face in &mesh.faces {
        let first = root(&mut parent, face[0]);
        for &v in &face[1..] {
            let other = root(&mut parent, v);
            parent[other] = first;
        }
    }
    let mut used: Vec<usize> = mesh
        .faces
        .iter()
        .flatten()
        .map(|&v| root(&mut parent, v))
        .collect();
    used.sort_unstable();
    used.dedup();
    used.len()
}

#[test]
fn collapses_never_change_the_topology_of_small_random_meshes() {
    // Random triangle soups on few vertices produce non-manifold vertices, rings, ears
    // and several components; however far they are reduced, the number of connected
    // pieces and the Euler characteristic must stay the same.
    let mut state = 0x5151_7777_1234_9999_u64;
    let mut below = move |n: usize| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % n as u64) as usize
    };
    let mut checked = 0;
    for _ in 0..6000 {
        let vertex_count = 4 + below(5);
        let mut mesh = MeshData::new();
        for _ in 0..vertex_count {
            let point = [
                below(1000) as f64 / 1000.0,
                below(1000) as f64 / 1000.0,
                below(300) as f64 / 1000.0,
            ];
            mesh.add_vertex(point);
        }
        for _ in 0..2 + below(5) {
            let (a, b, c) = (
                below(vertex_count),
                below(vertex_count),
                below(vertex_count),
            );
            let mut key = [a, b, c];
            key.sort_unstable();
            let repeated = key[0] == key[1] || key[1] == key[2];
            let known = mesh.faces.iter().any(|f| {
                let mut other = *f;
                other.sort_unstable();
                other == key
            });
            if !repeated && !known {
                mesh.add_face(a, b, c);
            }
        }
        if mesh.faces.len() < 2 || edge_uses(&mesh).values().any(|&uses| uses > 2) {
            continue;
        }
        for drop in [1, 2, 3, 100] {
            let out = simplified(&mesh, mesh.faces.len().saturating_sub(drop));
            assert_eq!(
                component_count(&out),
                component_count(&mesh),
                "faces {:?} -> {:?}",
                mesh.faces,
                out.faces
            );
            assert_eq!(
                euler_characteristic(&out),
                euler_characteristic(&mesh),
                "faces {:?} -> {:?}",
                mesh.faces,
                out.faces
            );
            checked += 1;
        }
    }
    assert!(checked > 2000, "only {checked} meshes were checked");
}

#[test]
fn a_fold_over_blocks_the_collapse() {
    // A fan around vertex 0 whose ring is not star-shaped from vertex 3.
    let mut mesh = MeshData::new();
    for p in [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ] {
        mesh.add_vertex(p);
    }
    for (a, b) in [(1, 2), (2, 3), (3, 4), (4, 1)] {
        mesh.add_face(0, a, b);
    }
    let simplifier = simplifier_for(&mesh);
    // Moving the centre next to vertex 3 and across the ring would flip triangles.
    assert!(simplifier.preserves_orientation(0, 1, [0.5, 0.0, 0.0]));
    assert!(!simplifier.preserves_orientation(0, 1, [-3.0, 0.0, 0.0]));
    assert!(simplifier.allows_collapse(0, 1, [0.5, 0.0, 0.0]));
    assert!(!simplifier.allows_collapse(0, 1, [-3.0, 0.0, 0.0]));
}

// ---------------------------------------------------------------------------
// Whole simplifications
// ---------------------------------------------------------------------------

#[test]
fn a_sphere_stays_a_closed_manifold_close_to_the_sphere() {
    let sphere = icosphere(2);
    assert_eq!(sphere.faces.len(), 320);
    let out = simplified(&sphere, 80);
    assert_well_formed(&out);
    assert_eq!(out.faces.len(), 80);
    assert!(
        is_closed_manifold(&out),
        "an edge is not shared by two triangles"
    );
    assert_eq!(euler_characteristic(&out), 2);
    assert_eq!(inverted_faces(&out), 0);
    for v in &out.vertices {
        let radius = norm(*v);
        assert!((0.9..=1.1).contains(&radius), "radius {radius}");
    }
    let (reference, volume) = (enclosed_volume(&sphere), enclosed_volume(&out));
    assert!(
        (volume - reference).abs() < 0.1 * reference,
        "volume {volume} against {reference}"
    );
}

#[test]
fn sphere_simplification_is_deterministic() {
    let sphere = icosphere(2);
    for target in [200, 80, 12] {
        let first = simplified(&sphere, target);
        let second = simplified(&sphere, target);
        assert_eq!(first.faces, second.faces, "target {target}");
        assert_eq!(first.vertices, second.vertices, "target {target}");
    }
}

#[test]
fn every_target_down_to_a_tetrahedron_keeps_the_sphere_closed() {
    let sphere = icosphere(1);
    for target in (4..=80).rev().step_by(2) {
        let out = simplified(&sphere, target);
        assert_well_formed(&out);
        assert!(
            out.faces.len() <= target,
            "target {target}: {}",
            out.faces.len()
        );
        assert!(is_closed_manifold(&out), "target {target}");
        assert_eq!(euler_characteristic(&out), 2, "target {target}");
        assert_eq!(inverted_faces(&out), 0, "target {target}");
    }
}

#[test]
fn a_tetrahedron_stays_a_tetrahedron() {
    let tetra = tetrahedron();
    for target in [3, 2, 0] {
        let out = simplified(&tetra, target);
        assert_eq!(out.faces.len(), 4, "target {target}");
        assert_eq!(out.vertices.len(), 4);
        assert!(is_closed_manifold(&out));
        assert_eq!(euler_characteristic(&out), 2);
        let mut seen: Vec<[usize; 3]> = out
            .faces
            .iter()
            .map(|f| {
                let mut sorted = *f;
                sorted.sort_unstable();
                sorted
            })
            .collect();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), 4, "back-to-back or duplicate triangles");
    }
}

#[test]
fn an_octahedron_loses_one_vertex_and_stays_closed() {
    let out = simplified(&octahedron(), 6);
    assert_well_formed(&out);
    assert_eq!(out.faces.len(), 6);
    assert_eq!(out.vertices.len(), 5);
    assert!(is_closed_manifold(&out));
    assert_eq!(euler_characteristic(&out), 2);
    assert_eq!(inverted_faces(&out), 0);
}

#[test]
fn a_subdivided_cube_reduces_to_the_twelve_triangle_cube() {
    let cube = subdivided_cube(6);
    assert_eq!(cube.faces.len(), 432);
    assert!(is_closed_manifold(&cube));
    assert_eq!(euler_characteristic(&cube), 2);
    // The quadrics of the flat sides pin every vertex to the cube: the corners are the
    // exact optimum of the three planes that meet there, and no collapse may cut them.
    let out = simplified(&cube, 12);
    assert_well_formed(&out);
    assert_eq!(out.faces.len(), 12);
    assert_eq!(out.vertices.len(), 8);
    for v in &out.vertices {
        assert!(
            v.iter().all(|&c| c == 0.0 || c == 1.0),
            "{v:?} is not a corner"
        );
    }
    assert!(is_closed_manifold(&out));
    assert_eq!(inverted_faces(&out), 0);
    assert!((enclosed_volume(&out) - 1.0).abs() < 1e-12);
}

#[test]
fn a_subdivided_cube_keeps_its_corners_and_stays_on_its_sides() {
    let cube = subdivided_cube(6);
    for target in [200, 100, 50, 24] {
        let out = simplified(&cube, target);
        assert_well_formed(&out);
        assert!(
            out.faces.len() <= target,
            "target {target}: {}",
            out.faces.len()
        );
        assert!(is_closed_manifold(&out), "target {target}");
        assert_eq!(inverted_faces(&out), 0, "target {target}");
        for v in &out.vertices {
            assert!(
                v.iter().all(|c| (0.0..=1.0).contains(c))
                    && v.iter().any(|&c| c == 0.0 || c == 1.0),
                "target {target}: {v:?} left the surface of the cube"
            );
        }
        let corners = out
            .vertices
            .iter()
            .filter(|v| v.iter().all(|&c| c == 0.0 || c == 1.0))
            .count();
        assert_eq!(corners, 8, "target {target}");
        assert!(
            (enclosed_volume(&out) - 1.0).abs() < 1e-9,
            "target {target}"
        );
    }
}

#[test]
fn a_torus_keeps_its_genus() {
    let donut = torus(24, 12, 3.0, 1.0);
    assert_eq!(euler_characteristic(&donut), 0);
    for target in [300, 150, 80, 40] {
        let out = simplified(&donut, target);
        assert_well_formed(&out);
        assert_eq!(out.faces.len(), target, "target {target}");
        assert!(is_closed_manifold(&out), "target {target}");
        assert_eq!(euler_characteristic(&out), 0, "target {target}");
        assert_eq!(inverted_faces(&out), 0, "target {target}");
    }
}

#[test]
fn a_planar_grid_stays_planar_and_keeps_its_outline() {
    let n = 12;
    let mesh = grid(n, n);
    let out = simplified(&mesh, 60);
    assert_well_formed(&out);
    assert!(
        out.faces.len() <= 60 && out.faces.len() >= 59,
        "{}",
        out.faces.len()
    );
    let side = n as f64;
    // Flat and boundary-preserving: no vertex leaves the plane, and the quadrics of the
    // flat faces never move one off its input position either.
    for v in &out.vertices {
        assert_eq!(v[2], 0.0);
        assert!(
            v[0].fract() == 0.0 && v[1].fract() == 0.0,
            "{v:?} is not an input vertex"
        );
        assert!((0.0..=side).contains(&v[0]) && (0.0..=side).contains(&v[1]));
    }
    // Every triangle is counter-clockwise and together they tile the square exactly.
    let areas = signed_areas_xy(&out);
    assert!(areas.iter().all(|&a| a > 0.0), "{areas:?}");
    assert!((areas.iter().sum::<f64>() - side * side).abs() < 1e-9);
    // The boundary runs along the sides of the square and the corners survive.
    let on_side = |(a, b): (usize, usize)| {
        let (p, q) = (out.vertices[a], out.vertices[b]);
        (p[0] == 0.0 && q[0] == 0.0)
            || (p[0] == side && q[0] == side)
            || (p[1] == 0.0 && q[1] == 0.0)
            || (p[1] == side && q[1] == side)
    };
    let boundary = boundary_edges(&out);
    assert!(!boundary.is_empty() && boundary.iter().all(|&e| on_side(e)));
    for corner in [[0.0, 0.0], [side, 0.0], [0.0, side], [side, side]] {
        assert!(
            out.vertices
                .iter()
                .any(|v| v[0] == corner[0] && v[1] == corner[1]),
            "corner {corner:?} was lost"
        );
    }
}

#[test]
fn a_large_planar_grid_is_reduced_consistently() {
    let n = 40;
    let out = simplified(&grid(n, n), 320);
    assert_well_formed(&out);
    assert!(
        out.faces.len() <= 320 && out.faces.len() >= 319,
        "{}",
        out.faces.len()
    );
    let areas = signed_areas_xy(&out);
    assert!(areas.iter().all(|&a| a > 0.0));
    assert!((areas.iter().sum::<f64>() - (n * n) as f64).abs() < 1e-8);
    assert_eq!(euler_characteristic(&out), 1);
}

#[test]
fn an_annulus_keeps_both_of_its_boundary_loops() {
    let ring = annulus(12);
    assert_eq!(euler_characteristic(&ring), 0);
    for target in [20, 12, 8, 4] {
        let out = simplified(&ring, target);
        assert_well_formed(&out);
        assert!(out.faces.len() >= 6, "an annulus needs two loops of three");
        assert_eq!(euler_characteristic(&out), 0, "target {target}");
        // Every vertex of a thin annulus lies on a boundary: edges to the other loop
        // are interior edges between boundary vertices, which are never collapsed.
        let boundary = boundary_edges(&out);
        let mut degree: BTreeMap<usize, usize> = BTreeMap::new();
        for &(a, b) in &boundary {
            *degree.entry(a).or_insert(0) += 1;
            *degree.entry(b).or_insert(0) += 1;
        }
        assert!(
            degree.values().all(|&d| d == 2),
            "target {target}: {degree:?}"
        );
        let radii: Vec<f64> = out.vertices.iter().map(|v| norm(*v)).collect();
        let inner = radii.iter().filter(|&&r| r < 1.5).count();
        assert!(
            inner >= 3 && radii.len() - inner >= 3,
            "target {target}: {radii:?}"
        );
    }
}

#[test]
fn separate_components_are_simplified_independently() {
    let first = icosphere(2);
    let mut both = first.clone();
    let offset = both.vertices.len();
    for v in &first.vertices {
        both.add_vertex([v[0] + 10.0, v[1], v[2]]);
    }
    for f in &first.faces {
        both.add_face(f[0] + offset, f[1] + offset, f[2] + offset);
    }
    let out = simplified(&both, 160);
    assert_well_formed(&out);
    assert_eq!(out.faces.len(), 160);
    assert!(is_closed_manifold(&out));
    assert_eq!(euler_characteristic(&out), 4);
    let (left, right): (Vec<&[usize; 3]>, Vec<&[usize; 3]>) =
        out.faces.iter().partition(|f| out.vertices[f[0]][0] < 5.0);
    assert_eq!((left.len(), right.len()), (80, 80));
}

#[test]
fn the_result_does_not_depend_on_where_the_mesh_is() {
    let mesh = grid(10, 10);
    let shift = [1.0e7, -3.0e7, 5.0e6];
    let mut moved = mesh.clone();
    for v in &mut moved.vertices {
        *v = add(*v, shift);
    }
    let near = simplified(&mesh, 40);
    let far = simplified(&moved, 40);
    assert_eq!(near.faces, far.faces);
    for (a, b) in near.vertices.iter().zip(&far.vertices) {
        assert_eq!(add(*a, shift), *b);
    }
    // Curved input far from the origin is simplified just as well: rounding of the
    // translated coordinates may reorder equal-cost collapses, but the result is still a
    // closed sphere-like mesh around the translated centre.
    let mut moved_sphere = icosphere(2);
    for v in &mut moved_sphere.vertices {
        *v = add(*v, shift);
    }
    let far = simplified(&moved_sphere, 100);
    assert_well_formed(&far);
    assert_eq!(far.faces.len(), 100);
    assert!(is_closed_manifold(&far));
    assert_eq!(euler_characteristic(&far), 2);
    assert_eq!(inverted_faces(&far), 0);
    for v in &far.vertices {
        let radius = norm(sub(*v, shift));
        assert!((0.95..=1.1).contains(&radius), "radius {radius}");
    }
}

#[test]
fn unusable_triangles_are_dropped_without_a_panic() {
    let mut mesh = grid(4, 4);
    let valid = mesh.faces.len();
    let nan_vertex = mesh.add_vertex([f64::NAN, 0.0, 0.0]);
    let inf_vertex = mesh.add_vertex([0.0, f64::INFINITY, 0.0]);
    mesh.add_face(0, 1, 999);
    mesh.add_face(usize::MAX, 1, 2);
    mesh.add_face(3, 3, 4);
    mesh.add_face(0, 0, 0);
    mesh.add_face(0, 1, nan_vertex);
    mesh.add_face(2, inf_vertex, 3);
    let out = simplified(&mesh, 10);
    assert_well_formed(&out);
    assert!(out.faces.len() <= 10);
    assert!(out.faces.len() < valid);
    // Nothing else was lost: the planar grid still tiles the square.
    let areas = signed_areas_xy(&out);
    assert!(areas.iter().all(|&a| a > 0.0));
    assert!((areas.iter().sum::<f64>() - 16.0).abs() < 1e-9);
}

#[test]
fn a_mesh_of_only_unusable_triangles_simplifies_to_nothing() {
    let mut mesh = MeshData::new();
    mesh.add_vertex([0.0; 3]);
    mesh.add_vertex([1.0, 0.0, 0.0]);
    for face in [[0, 0, 1], [0, 1, 7], [5, 6, 7]] {
        mesh.add_face(face[0], face[1], face[2]);
    }
    let out = simplified(&mesh, 1);
    assert!(out.faces.is_empty() && out.vertices.is_empty());
}

#[test]
fn material_ids_follow_the_surviving_triangles() {
    // Two flat patches far apart; a triangle never changes patch, so its material
    // identifies the patch it came from.
    let patch = grid(6, 6);
    let mut mesh = patch.clone();
    let offset = mesh.vertices.len();
    for v in &patch.vertices {
        mesh.add_vertex([v[0] + 100.0, v[1], v[2]]);
    }
    for f in &patch.faces {
        mesh.add_face(f[0] + offset, f[1] + offset, f[2] + offset);
    }
    mesh.material_ids = (0..mesh.faces.len())
        .map(|i| if i < patch.faces.len() { 7 } else { 9 })
        .collect();
    let out = simplified(&mesh, 30);
    assert_well_formed(&out);
    assert_eq!(out.material_ids.len(), out.faces.len());
    for (face, &material) in out.faces.iter().zip(&out.material_ids) {
        let x = out.vertices[face[0]][0];
        assert_eq!(material, if x < 50.0 { 7 } else { 9 }, "face {face:?}");
    }
    assert!(out.material_ids.contains(&7) && out.material_ids.contains(&9));
    // Materials that do not describe every triangle are not carried over.
    let mut partial = mesh.clone();
    partial.material_ids.pop();
    assert!(simplified(&partial, 30).material_ids.is_empty());
}

#[test]
fn a_target_at_or_above_the_triangle_count_returns_the_mesh_unchanged() {
    let mesh = grid(3, 3);
    for target in [18, 19, 1000, usize::MAX] {
        let out = simplified(&mesh, target);
        assert_eq!(out.faces, mesh.faces);
        assert_eq!(out.vertices, mesh.vertices);
    }
}

#[test]
fn target_zero_reduces_as_far_as_the_surface_allows() {
    let out = simplified(&icosphere(1), 0);
    assert_eq!(out.faces.len(), 4);
    assert!(is_closed_manifold(&out));
    // A disc is reduced to a single triangle; no piece of the surface is ever deleted.
    let out = simplified(&grid(3, 3), 0);
    assert_well_formed(&out);
    assert_eq!(out.faces.len(), 1);
    assert_eq!(out.vertices.len(), 3);
}

#[test]
fn flat_regions_keep_the_input_coordinates_bit_for_bit() {
    // Irregular spacing and coordinates of very different sizes: the simplifier works on
    // shifted copies of the coordinates, yet every vertex that a collapse leaves on an
    // endpoint must come out with exactly the coordinates it went in with.
    let mut mesh = grid(9, 9);
    for v in &mut mesh.vertices {
        *v = [v[0] * 0.37 + 0.001, v[1] * v[1] * 0.113 + 0.0007, 3.3];
    }
    let input: BTreeSet<[u64; 3]> = mesh.vertices.iter().map(|v| v.map(f64::to_bits)).collect();
    let out = simplified(&mesh, 40);
    assert_well_formed(&out);
    assert!(out.faces.len() <= 40, "{}", out.faces.len());
    for v in &out.vertices {
        assert!(
            input.contains(&v.map(f64::to_bits)),
            "{v:?} is not one of the input vertices"
        );
    }
}

#[test]
fn rough_surfaces_keep_a_consistent_adjacency_at_every_target() {
    // A deterministic pseudo-random height field.
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut mesh = grid(14, 11);
    for v in &mut mesh.vertices {
        v[2] = 0.8 * next();
    }
    for target in [250, 100, 40, 12, 1, 0] {
        let mut simplifier = simplifier_for(&mesh);
        simplifier.run(target);
        assert_adjacency_is_consistent(&simplifier);
        let input_indices: Vec<usize> = (0..mesh.faces.len()).collect();
        let (faces, indices) = simplifier.surviving_faces(&input_indices);
        assert_eq!(faces.len(), simplifier.live_faces);
        assert!(
            indices.windows(2).all(|pair| pair[0] < pair[1]),
            "input order is kept"
        );
        let out = simplified(&mesh, target);
        assert_well_formed(&out);
        assert_eq!(out.faces.len(), simplifier.live_faces);
        assert!(
            out.faces.len() <= target.max(1),
            "target {target}: {}",
            out.faces.len()
        );
        assert!(
            edge_uses(&out).values().all(|&c| c <= 2),
            "non-manifold edge"
        );
        // A disc stays a disc.
        assert_eq!(euler_characteristic(&out), 1, "target {target}");
    }
}
