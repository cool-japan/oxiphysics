// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`super::mesh_boolean`] on generic (mid-triangle)
//! intersections: results must be split along the intersection curves, closed
//! with consistent orientation, and have the exact boolean volume.

use std::collections::HashMap;

use super::{
    MeshBooleanOp, SimpleMesh, V3, add, cross, dot, euler_characteristic, is_watertight,
    mesh_boolean, normalize, sub,
};

/// Every directed edge appears exactly once and its reverse exactly once:
/// a closed, consistently oriented 2-manifold surface.
fn is_oriented_closed(mesh: &SimpleMesh) -> bool {
    let mut directed: HashMap<(usize, usize), u32> = HashMap::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            *directed.entry((t[k], t[(k + 1) % 3])).or_insert(0) += 1;
        }
    }
    directed
        .iter()
        .all(|(&(a, b), &c)| c == 1 && directed.get(&(b, a)) == Some(&1))
}

fn volume(m: &SimpleMesh) -> f64 {
    m.signed_volume()
}

/// Signed volume measured about `o` (well conditioned far from the origin).
fn volume_about(m: &SimpleMesh, o: V3) -> f64 {
    m.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| sub(m.vertices[i], o));
            dot(a, cross(b, c))
        })
        .sum::<f64>()
        / 6.0
}

fn assert_closed(m: &SimpleMesh, what: &str) {
    assert!(!m.triangles.is_empty(), "{what}: empty result");
    assert!(is_watertight(m), "{what}: not watertight");
    assert!(is_oriented_closed(m), "{what}: inconsistent orientation");
}

fn rotate(v: V3, axis: V3, angle: f64) -> V3 {
    let k = normalize(axis);
    let (s, c) = angle.sin_cos();
    let kxv = cross(k, v);
    let kdv = dot(k, v);
    [
        v[0] * c + kxv[0] * s + k[0] * kdv * (1.0 - c),
        v[1] * c + kxv[1] * s + k[1] * kdv * (1.0 - c),
        v[2] * c + kxv[2] * s + k[2] * kdv * (1.0 - c),
    ]
}

fn rotated_cube(centre: V3, half: f64, axis: V3, angle: f64) -> SimpleMesh {
    let mut m = SimpleMesh::unit_cube([0.0; 3], half);
    for v in &mut m.vertices {
        *v = add(rotate(*v, axis, angle), centre);
    }
    m
}

/// Closed octahedron `|x| + |y| + |z| <= r` about `c`, outward winding.
fn octahedron(c: V3, r: f64) -> SimpleMesh {
    let vertices = vec![
        [c[0] + r, c[1], c[2]],
        [c[0] - r, c[1], c[2]],
        [c[0], c[1] + r, c[2]],
        [c[0], c[1] - r, c[2]],
        [c[0], c[1], c[2] + r],
        [c[0], c[1], c[2] - r],
    ];
    let triangles = vec![
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    SimpleMesh::from_data(vertices, triangles)
}

#[test]
fn offset_unit_cubes_give_exact_watertight_booleans() {
    // A = [-0.5, 0.5]^3, B = A + (0.5, 0, 0): mid-face cuts on the ±X faces and
    // overlapping coplanar ±Y / ±Z faces.
    let a = SimpleMesh::unit_cube([0.0; 3], 0.5);
    let b = SimpleMesh::unit_cube([0.5, 0.0, 0.0], 0.5);

    let union = mesh_boolean(&a, &b, MeshBooleanOp::Union);
    assert_closed(&union, "union");
    assert!(
        (volume(&union) - 1.5).abs() < 1e-9,
        "union volume {}",
        volume(&union)
    );
    assert_eq!(euler_characteristic(&union), 2);

    let diff = mesh_boolean(&a, &b, MeshBooleanOp::Difference);
    assert_closed(&diff, "difference");
    assert!(
        (volume(&diff) - 0.5).abs() < 1e-9,
        "difference volume {}",
        volume(&diff)
    );
    assert_eq!(euler_characteristic(&diff), 2);

    let inter = mesh_boolean(&a, &b, MeshBooleanOp::Intersection);
    assert_closed(&inter, "intersection");
    assert!(
        (volume(&inter) - 0.5).abs() < 1e-9,
        "intersection volume {}",
        volume(&inter)
    );
    assert_eq!(euler_characteristic(&inter), 2);

    // The result really is the box [-0.5, 1] x [-0.5, 0.5]^2: every vertex lies
    // on its boundary and the bounding box matches.
    let (lo, hi) = union.aabb();
    for (got, want) in lo.iter().chain(&hi).zip([-0.5, -0.5, -0.5, 1.0, 0.5, 0.5]) {
        assert!((got - want).abs() < 1e-12, "union bbox {lo:?} {hi:?}");
    }
    // The split faces must not keep any surface strictly inside the union.
    for t in &union.triangles {
        let c = super::scale(
            add(
                add(union.vertices[t[0]], union.vertices[t[1]]),
                union.vertices[t[2]],
            ),
            1.0 / 3.0,
        );
        let on_boundary = (c[0] + 0.5).abs() < 1e-12
            || (c[0] - 1.0).abs() < 1e-12
            || (c[1].abs() - 0.5).abs() < 1e-12
            || (c[2].abs() - 0.5).abs() < 1e-12;
        assert!(on_boundary, "interior face at centroid {c:?}");
    }
}

#[test]
fn cube_and_octahedron_match_analytic_volumes() {
    // The octahedron's six tips poke through the centres of the cube faces, so
    // every cube face is cut by a closed diamond-shaped seam.
    let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
    let octa = octahedron([0.0; 3], 1.4);
    let v_cube = 8.0;
    let v_octa = 4.0 / 3.0 * 1.4_f64.powi(3);
    // Each tip beyond a face is a pyramid with a 0.4-diamond base and height 0.4.
    let tip = (2.0 * 0.4 * 0.4) * 0.4 / 3.0;
    let v_inter = v_octa - 6.0 * tip;

    let inter = mesh_boolean(&cube, &octa, MeshBooleanOp::Intersection);
    assert_closed(&inter, "cube ∩ octahedron");
    assert!(
        (volume(&inter) - v_inter).abs() < 1e-9,
        "∩ {}",
        volume(&inter)
    );

    let union = mesh_boolean(&cube, &octa, MeshBooleanOp::Union);
    assert_closed(&union, "cube ∪ octahedron");
    assert!(
        (volume(&union) - (v_cube + v_octa - v_inter)).abs() < 1e-9,
        "∪ {}",
        volume(&union)
    );

    let diff = mesh_boolean(&octa, &cube, MeshBooleanOp::Difference);
    assert_closed(&diff, "octahedron \\ cube");
    assert!(
        (volume(&diff) - 6.0 * tip).abs() < 1e-9,
        "\\ {}",
        volume(&diff)
    );
    // Six disjoint tips: Euler characteristic 2 per component.
    assert_eq!(euler_characteristic(&diff), 12);
}

#[test]
fn rotated_cubes_satisfy_inclusion_exclusion() {
    let a = SimpleMesh::unit_cube([0.0; 3], 1.0);
    let b = rotated_cube([0.7, 0.4, -0.3], 0.9, [0.3, -0.5, 0.8], 0.61);
    let (va, vb) = (volume(&a), volume(&b));
    assert!((vb - 1.8_f64.powi(3)).abs() < 1e-12);

    let union = mesh_boolean(&a, &b, MeshBooleanOp::Union);
    let inter = mesh_boolean(&a, &b, MeshBooleanOp::Intersection);
    let diff_ab = mesh_boolean(&a, &b, MeshBooleanOp::Difference);
    let diff_ba = mesh_boolean(&b, &a, MeshBooleanOp::Difference);
    for (m, what) in [
        (&union, "union"),
        (&inter, "intersection"),
        (&diff_ab, "A \\ B"),
        (&diff_ba, "B \\ A"),
    ] {
        assert_closed(m, what);
    }
    let vi = volume(&inter);
    assert!(vi > 0.5 && vi < vb, "implausible intersection volume {vi}");
    assert!((volume(&union) + vi - (va + vb)).abs() < 1e-9);
    assert!((volume(&diff_ab) - (va - vi)).abs() < 1e-9);
    assert!((volume(&diff_ba) - (vb - vi)).abs() < 1e-9);
}

#[test]
fn seam_loop_inside_a_single_face_triangle() {
    // A small cube straddling the big cube's +X face away from the face
    // diagonal: the whole seam is a closed square inside one triangle, so the
    // re-triangulated triangle has a hole-like inner boundary.
    let big = SimpleMesh::unit_cube([0.0; 3], 1.0);
    let small = SimpleMesh::unit_cube([1.0, 0.5, -0.5], 0.2);
    let half_small = 0.4_f64.powi(3) / 2.0;

    let union = mesh_boolean(&big, &small, MeshBooleanOp::Union);
    assert_closed(&union, "union");
    assert!((volume(&union) - (8.0 + half_small)).abs() < 1e-9);

    let diff = mesh_boolean(&big, &small, MeshBooleanOp::Difference);
    assert_closed(&diff, "difference");
    assert!((volume(&diff) - (8.0 - half_small)).abs() < 1e-9);

    let inter = mesh_boolean(&big, &small, MeshBooleanOp::Intersection);
    assert_closed(&inter, "intersection");
    assert!((volume(&inter) - half_small).abs() < 1e-9);
}

#[test]
fn split_pieces_keep_parent_orientation() {
    // Union of cube and octahedron: every output triangle must face away from
    // the (star-shaped about the origin) solid.
    let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
    let octa = octahedron([0.0; 3], 1.4);
    let union = mesh_boolean(&cube, &octa, MeshBooleanOp::Union);
    assert!(union.n_triangles() > 20, "faces were not split");
    for t in &union.triangles {
        let [p, q, r] = t.map(|i| union.vertices[i]);
        let n = cross(sub(q, p), sub(r, p));
        let c = add(add(p, q), r);
        assert!(dot(n, c) > 0.0, "inward-facing piece {p:?} {q:?} {r:?}");
    }
}

#[test]
fn invalid_triangles_are_ignored() {
    let a = SimpleMesh::unit_cube([0.0; 3], 0.5);
    let mut b = SimpleMesh::unit_cube([0.5, 0.0, 0.0], 0.5);
    // Out-of-range index, repeated index and a non-finite vertex.
    b.triangles.push([0, 1, 999]);
    b.triangles.push([2, 2, 3]);
    let bad = b.add_vertex([f64::NAN, 0.0, 0.0]);
    b.triangles.push([0, 1, bad]);
    let union = mesh_boolean(&a, &b, MeshBooleanOp::Union);
    assert_closed(&union, "union with sanitised operand");
    assert!((volume(&union) - 1.5).abs() < 1e-9);
}

#[test]
fn translated_far_from_origin_still_exact() {
    // Large absolute coordinates: tolerances scale with the coordinates.
    let o = [1.0e4, -2.0e4, 3.0e4];
    let a = SimpleMesh::unit_cube(o, 0.5);
    let b = SimpleMesh::unit_cube([o[0] + 0.5, o[1] + 0.25, o[2] - 0.125], 0.5);
    let union = mesh_boolean(&a, &b, MeshBooleanOp::Union);
    let inter = mesh_boolean(&a, &b, MeshBooleanOp::Intersection);
    assert_closed(&union, "union");
    assert_closed(&inter, "intersection");
    let vi = 0.5 * 0.75 * 0.875;
    let (v_inter, v_union) = (volume_about(&inter, o), volume_about(&union, o));
    assert!((v_inter - vi).abs() < 1e-6, "∩ {v_inter}");
    assert!((v_union - (2.0 - vi)).abs() < 1e-6, "∪ {v_union}");
}

/// Deterministic LCG in `[0, 1]` for reproducible random configurations.
struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 33) as f64) / f64::from(u32::MAX)
    }
}

/// Icosphere of radius `r` about `c` (`levels` midpoint subdivisions of an
/// icosahedron), rotated by `angle` about a fixed oblique axis, outward winding.
fn icosphere(c: V3, r: f64, levels: usize, angle: f64) -> SimpleMesh {
    let g = (1.0 + 5f64.sqrt()) / 2.0;
    let mut v: Vec<V3> = vec![
        [-1.0, g, 0.0],
        [1.0, g, 0.0],
        [-1.0, -g, 0.0],
        [1.0, -g, 0.0],
        [0.0, -1.0, g],
        [0.0, 1.0, g],
        [0.0, -1.0, -g],
        [0.0, 1.0, -g],
        [g, 0.0, -1.0],
        [g, 0.0, 1.0],
        [-g, 0.0, -1.0],
        [-g, 0.0, 1.0],
    ];
    let mut f: Vec<[usize; 3]> = vec![
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
    for _ in 0..levels {
        let mut cache: HashMap<(usize, usize), usize> = HashMap::new();
        let mut next = Vec::with_capacity(f.len() * 4);
        for tri in &f {
            let mut mid = [0usize; 3];
            for k in 0..3 {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                let key = if a < b { (a, b) } else { (b, a) };
                mid[k] = *cache.entry(key).or_insert_with(|| {
                    v.push(super::scale(add(v[a], v[b]), 0.5));
                    v.len() - 1
                });
            }
            next.push([tri[0], mid[0], mid[2]]);
            next.push([tri[1], mid[1], mid[0]]);
            next.push([tri[2], mid[2], mid[1]]);
            next.push(mid);
        }
        f = next;
    }
    let verts = v
        .iter()
        .map(|p| {
            add(
                rotate(super::scale(normalize(*p), r), [0.2, 0.7, -0.4], angle),
                c,
            )
        })
        .collect();
    SimpleMesh::from_data(verts, f)
}

/// Union, intersection and difference of `a` and `b` are closed, consistently
/// oriented, and satisfy `|A ∪ B| + |A ∩ B| = |A| + |B|` and
/// `|A \ B| = |A| - |A ∩ B|`.
fn assert_boolean_identities(a: &SimpleMesh, b: &SimpleMesh, what: &str) {
    let union = mesh_boolean(a, b, MeshBooleanOp::Union);
    let inter = mesh_boolean(a, b, MeshBooleanOp::Intersection);
    let diff = mesh_boolean(a, b, MeshBooleanOp::Difference);
    for (m, op) in [
        (&union, "union"),
        (&inter, "intersection"),
        (&diff, "difference"),
    ] {
        assert!(
            is_watertight(m) && is_oriented_closed(m),
            "{what}: {op} is not a closed oriented surface"
        );
    }
    let (va, vb, vi) = (volume(a), volume(b), volume(&inter));
    let scale = va.abs() + vb.abs();
    assert!(vi > 0.0, "{what}: empty intersection");
    assert!(
        (volume(&union) + vi - (va + vb)).abs() < 1e-9 * scale,
        "{what}: union {} intersection {vi} operands {va} {vb}",
        volume(&union)
    );
    assert!(
        (volume(&diff) - (va - vi)).abs() < 1e-9 * scale,
        "{what}: difference {} intersection {vi} operand {va}",
        volume(&diff)
    );
}

#[test]
fn randomized_overlaps_are_closed_and_volume_exact() {
    // Small random rotations / offsets of a cube and an octahedron around
    // (1, 0, 0) produce many nearly aligned (but resolvable) cuts; random
    // icosphere pairs produce long seams across many small triangles.
    let mut rng = Lcg(0x9E37_79B9_7F4A_7C15);
    let cube = SimpleMesh::unit_cube([0.0; 3], 1.0);
    for it in 0..12 {
        let tiny = 1e-8;
        let o = [
            1.0 + (rng.next_f64() - 0.5) * tiny,
            (rng.next_f64() - 0.5) * tiny,
            (rng.next_f64() - 0.5) * tiny,
        ];
        let angle = (rng.next_f64() - 0.5) * 0.1;
        let axis = [
            rng.next_f64() - 0.5,
            rng.next_f64() - 0.5,
            rng.next_f64() - 0.5,
        ];
        let rotated = rotated_cube(o, 1.0, axis, angle);
        assert_boolean_identities(&cube, &rotated, &format!("cube/rotated cube #{it}"));
        let octa = octahedron(o, 1.0);
        assert_boolean_identities(&cube, &octa, &format!("cube/octahedron #{it}"));
    }
    for it in 0..3 {
        let s1 = icosphere([0.3 * rng.next_f64(), 0.2, 0.1], 1.0, 2, 0.0);
        let s2 = icosphere([0.7, 0.3 * rng.next_f64(), -0.2], 0.8, 2, rng.next_f64());
        assert_boolean_identities(&s1, &s2, &format!("icosphere pair #{it}"));
        assert_boolean_identities(&cube, &s2, &format!("cube/icosphere #{it}"));
    }
}
