// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Winding-order analysis of triangle meshes.
//!
//! A surface is *consistently oriented* when every two triangles that share an edge
//! traverse it in opposite directions.  Triangles that break this rule have their
//! winding reversed relative to their neighbours; they are the "inverted" triangles
//! of a mesh.  The analysis is purely topological, so it does not depend on the
//! direction in which the mesh happens to face.
//!
//! The mesh is split into connected components (triangles joined through edges that
//! exactly two triangles share).  The orientation of a seed triangle is propagated
//! breadth-first, which divides the component into triangles that agree with the seed
//! and triangles that must be flipped to agree.  Which of the two groups is the
//! inverted one is decided as follows:
//!
//! * A **closed** component (every edge shared by two triangles) bounds a volume, so
//!   the sign of its signed volume says which orientation points outwards: the group
//!   that gives a negative volume is the inverted one.  A closed component that is
//!   consistently oriented but inside out is therefore reported as entirely inverted.
//! * An **open** component has no inside; the smaller group is reported.

use std::collections::VecDeque;

/// The connection of a triangle to its neighbour across one edge.
#[derive(Clone, Copy, Debug)]
struct Link {
    /// The neighbouring triangle.
    face: usize,
    /// Both triangles traverse the shared edge in the same direction, so exactly
    /// one of them must be flipped for the pair to be consistent.
    same_direction: bool,
}

/// Relative size below which a component's signed volume is treated as zero.
const VOLUME_TOLERANCE: f64 = 1e-9;

/// Number of triangles in `faces` whose winding disagrees with the surface they belong to.
///
/// Triangles with a repeated vertex index or an index outside `vertices` have no
/// orientation and are ignored.
pub(crate) fn count_inverted_faces(vertices: &[[f64; 3]], faces: &[[u32; 3]]) -> usize {
    let usable = |face: &[u32; 3]| {
        face[0] != face[1]
            && face[1] != face[2]
            && face[0] != face[2]
            && face.iter().all(|&v| (v as usize) < vertices.len())
    };

    // One record per triangle edge: (low vertex, high vertex, triangle, edge slot,
    // traversed from the low to the high vertex).
    let mut half_edges: Vec<(u32, u32, usize, usize, bool)> = Vec::with_capacity(3 * faces.len());
    for (index, face) in faces.iter().enumerate() {
        if !usable(face) {
            continue;
        }
        for slot in 0..3 {
            let from = face[slot];
            let to = face[(slot + 1) % 3];
            half_edges.push((from.min(to), from.max(to), index, slot, from < to));
        }
    }
    half_edges.sort_unstable();

    // Link the two triangles of every edge that exactly two triangles share.
    let mut links: Vec<[Option<Link>; 3]> = vec![[None; 3]; faces.len()];
    let mut start = 0;
    while start < half_edges.len() {
        let (low, high, ..) = half_edges[start];
        let mut end = start + 1;
        while end < half_edges.len() && half_edges[end].0 == low && half_edges[end].1 == high {
            end += 1;
        }
        if end - start == 2 {
            let (_, _, face_a, slot_a, forward_a) = half_edges[start];
            let (_, _, face_b, slot_b, forward_b) = half_edges[start + 1];
            let same_direction = forward_a == forward_b;
            links[face_a][slot_a] = Some(Link {
                face: face_b,
                same_direction,
            });
            links[face_b][slot_b] = Some(Link {
                face: face_a,
                same_direction,
            });
        }
        start = end;
    }

    // Propagate the orientation of a seed triangle through every component.
    // `flipped[f]` tells whether triangle `f` must be flipped to match its seed.
    let mut flipped: Vec<Option<bool>> = vec![None; faces.len()];
    let mut queue: VecDeque<usize> = VecDeque::new();
    let mut component: Vec<usize> = Vec::new();
    let mut inverted = 0usize;
    for seed in 0..faces.len() {
        if flipped[seed].is_some() || !usable(&faces[seed]) {
            continue;
        }
        component.clear();
        flipped[seed] = Some(false);
        queue.push_back(seed);
        while let Some(face) = queue.pop_front() {
            component.push(face);
            let own = flipped[face].unwrap_or(false);
            for link in links[face].iter().flatten() {
                if flipped[link.face].is_none() {
                    flipped[link.face] = Some(own != link.same_direction);
                    queue.push_back(link.face);
                }
            }
        }
        inverted += inverted_in_component(vertices, faces, &links, &flipped, &component);
    }
    inverted
}

/// Inverted triangles of one connected component (see the module documentation).
fn inverted_in_component(
    vertices: &[[f64; 3]],
    faces: &[[u32; 3]],
    links: &[[Option<Link>; 3]],
    flipped: &[Option<bool>],
    component: &[usize],
) -> usize {
    // Measure volumes from the centre of the component to keep them well conditioned.
    let mut centre = [0.0f64; 3];
    let mut corners = 0usize;
    for &face in component {
        for &vertex in &faces[face] {
            for (sum, coordinate) in centre.iter_mut().zip(vertices[vertex as usize]) {
                *sum += coordinate;
            }
            corners += 1;
        }
    }
    for value in &mut centre {
        *value /= corners.max(1) as f64;
    }

    let mut with_seed = 0usize;
    let mut against_seed = 0usize;
    let mut closed = true;
    let mut volume = 0.0f64;
    let mut absolute_volume = 0.0f64;
    for &face in component {
        let against = flipped[face] == Some(true);
        if against {
            against_seed += 1;
        } else {
            with_seed += 1;
        }
        if links[face].iter().any(Option::is_none) {
            closed = false;
        }
        let [a, b, c] = faces[face].map(|v| {
            let p = vertices[v as usize];
            [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]]
        });
        let triple = (a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            / 6.0;
        volume += if against { -triple } else { triple };
        absolute_volume += triple.abs();
    }

    if closed && volume.abs() > VOLUME_TOLERANCE * absolute_volume {
        // With the seed's orientation the surface encloses `volume`: positive means
        // the seed's side points outwards, so the triangles against it are inverted.
        if volume > 0.0 {
            against_seed
        } else {
            with_seed
        }
    } else {
        with_seed.min(against_seed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Regular icosahedron with outward-facing triangles.
    fn icosahedron() -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
        let t = (1.0 + 5.0f64.sqrt()) / 2.0;
        let vertices = vec![
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
        ];
        let faces: Vec<[u32; 3]> = vec![
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
        let faces = faces
            .into_iter()
            .map(|[a, b, c]| {
                let (pa, pb, pc) = (
                    vertices[a as usize],
                    vertices[b as usize],
                    vertices[c as usize],
                );
                let ab = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
                let ac = [pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]];
                let n = [
                    ab[1] * ac[2] - ab[2] * ac[1],
                    ab[2] * ac[0] - ab[0] * ac[2],
                    ab[0] * ac[1] - ab[1] * ac[0],
                ];
                let outward = n[0] * pa[0] + n[1] * pa[1] + n[2] * pa[2];
                if outward > 0.0 { [a, b, c] } else { [a, c, b] }
            })
            .collect();
        (vertices, faces)
    }

    /// Icosphere with `20 * 4^level` outward-facing triangles.
    fn icosphere(level: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
        let (mut vertices, mut faces) = icosahedron();
        for p in &mut vertices {
            let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            *p = [p[0] / n, p[1] / n, p[2] / n];
        }
        for _ in 0..level {
            let mut midpoints: HashMap<(u32, u32), u32> = HashMap::new();
            let mut next = Vec::with_capacity(faces.len() * 4);
            for &[a, b, c] in &faces {
                let mut middle = |u: u32, v: u32| -> u32 {
                    *midpoints.entry((u.min(v), u.max(v))).or_insert_with(|| {
                        let (pu, pv) = (vertices[u as usize], vertices[v as usize]);
                        let m = [
                            (pu[0] + pv[0]) / 2.0,
                            (pu[1] + pv[1]) / 2.0,
                            (pu[2] + pv[2]) / 2.0,
                        ];
                        let n = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
                        vertices.push([m[0] / n, m[1] / n, m[2] / n]);
                        (vertices.len() - 1) as u32
                    })
                };
                let ab = middle(a, b);
                let bc = middle(b, c);
                let ca = middle(c, a);
                next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
            }
            faces = next;
        }
        (vertices, faces)
    }

    /// `n x n` cells in the z = 0 plane, two triangles per cell, normals along +z.
    fn grid(n: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
        let stride = n + 1;
        let mut vertices = Vec::new();
        for j in 0..=n {
            for i in 0..=n {
                vertices.push([i as f64, j as f64, 0.0]);
            }
        }
        let mut faces = Vec::new();
        for j in 0..n {
            for i in 0..n {
                let a = (j * stride + i) as u32;
                let b = a + 1;
                let c = a + stride as u32;
                let d = c + 1;
                faces.push([a, b, d]);
                faces.push([a, d, c]);
            }
        }
        (vertices, faces)
    }

    fn flip(face: &mut [u32; 3]) {
        face.swap(1, 2);
    }

    #[test]
    fn the_test_icosahedron_is_a_closed_manifold() {
        let (vertices, faces) = icosahedron();
        assert_eq!((vertices.len(), faces.len()), (12, 20));
        let mut uses: HashMap<(u32, u32), Vec<bool>> = HashMap::new();
        for f in &faces {
            for k in 0..3 {
                let (u, v) = (f[k], f[(k + 1) % 3]);
                uses.entry((u.min(v), u.max(v))).or_default().push(u < v);
            }
        }
        assert_eq!(uses.len(), 30);
        for (edge, dirs) in uses {
            assert_eq!(dirs.len(), 2, "edge {edge:?}");
            assert_ne!(
                dirs[0], dirs[1],
                "edge {edge:?} is traversed twice the same way"
            );
        }
    }

    #[test]
    fn consistently_oriented_closed_meshes_have_no_inverted_faces() {
        for level in 0..=3 {
            let (vertices, faces) = icosphere(level);
            assert_eq!(faces.len(), 20 * 4usize.pow(level as u32));
            assert_eq!(count_inverted_faces(&vertices, &faces), 0, "level {level}");
        }
    }

    #[test]
    fn a_single_flipped_face_is_the_only_inverted_face() {
        let (vertices, faces) = icosphere(2);
        for index in [0, 1, 17, 100, faces.len() - 1] {
            let mut broken = faces.clone();
            flip(&mut broken[index]);
            assert_eq!(
                count_inverted_faces(&vertices, &broken),
                1,
                "face {index} flipped"
            );
        }
    }

    #[test]
    fn several_flipped_faces_are_counted_exactly() {
        let (vertices, faces) = icosphere(2);
        let mut broken = faces.clone();
        for index in [3, 4, 40, 41, 200, 300, 301] {
            flip(&mut broken[index]);
        }
        assert_eq!(count_inverted_faces(&vertices, &broken), 7);
    }

    #[test]
    fn an_inside_out_closed_mesh_is_inverted_everywhere() {
        let (vertices, mut faces) = icosphere(1);
        for face in &mut faces {
            flip(face);
        }
        assert_eq!(count_inverted_faces(&vertices, &faces), 80);
    }

    #[test]
    fn the_verdict_does_not_depend_on_where_the_mesh_is() {
        let (mut vertices, mut faces) = icosphere(1);
        for p in &mut vertices {
            *p = [p[0] * 3.0 + 1.0e4, p[1] * 3.0 - 2.5e4, p[2] * 3.0 + 7.0];
        }
        assert_eq!(count_inverted_faces(&vertices, &faces), 0);
        flip(&mut faces[5]);
        assert_eq!(count_inverted_faces(&vertices, &faces), 1);
    }

    #[test]
    fn components_are_judged_independently() {
        let (vertices_a, faces_a) = icosphere(1);
        let mut vertices = vertices_a.clone();
        let mut faces = faces_a.clone();
        let offset = vertices.len() as u32;
        vertices.extend(vertices_a.iter().map(|p| [p[0] + 10.0, p[1], p[2]]));
        faces.extend(
            faces_a
                .iter()
                .map(|f| [f[0] + offset, f[1] + offset, f[2] + offset]),
        );
        assert_eq!(count_inverted_faces(&vertices, &faces), 0);

        // Turn the second sphere inside out: its 80 faces are inverted, the first sphere's are not.
        let mut inside_out = faces.clone();
        for face in inside_out.iter_mut().skip(faces_a.len()) {
            flip(face);
        }
        assert_eq!(count_inverted_faces(&vertices, &inside_out), 80);

        // One bad face in each sphere.
        let mut two = faces.clone();
        flip(&mut two[2]);
        flip(&mut two[faces_a.len() + 9]);
        assert_eq!(count_inverted_faces(&vertices, &two), 2);
    }

    #[test]
    fn open_surfaces_report_the_minority_orientation() {
        let (vertices, faces) = grid(4);
        assert_eq!(count_inverted_faces(&vertices, &faces), 0);
        for index in [0, 5, faces.len() - 1] {
            let mut broken = faces.clone();
            flip(&mut broken[index]);
            assert_eq!(count_inverted_faces(&vertices, &broken), 1, "face {index}");
        }
        let mut broken = faces.clone();
        for index in [6, 7, 8] {
            flip(&mut broken[index]);
        }
        assert_eq!(count_inverted_faces(&vertices, &broken), 3);
        // A consistently reversed open sheet cannot be told apart from a sheet seen from below.
        let mut reversed = faces.clone();
        for face in &mut reversed {
            flip(face);
        }
        assert_eq!(count_inverted_faces(&vertices, &reversed), 0);
    }

    #[test]
    fn faces_without_an_orientation_are_ignored() {
        let (vertices, mut faces) = grid(2);
        faces.push([0, 0, 1]);
        faces.push([3, 3, 3]);
        faces.push([0, 1, 99]);
        faces.push([u32::MAX, 1, 2]);
        assert_eq!(count_inverted_faces(&vertices, &faces), 0);
        assert_eq!(count_inverted_faces(&[], &[]), 0);
        assert_eq!(count_inverted_faces(&vertices, &[]), 0);
    }

    #[test]
    fn non_manifold_edges_do_not_link_faces() {
        // Three triangles share the edge 0-1: no unique neighbour, so no verdict.
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let faces = vec![[0, 1, 2], [0, 1, 3], [0, 1, 4]];
        assert_eq!(count_inverted_faces(&vertices, &faces), 0);
    }

    #[test]
    fn duplicated_coincident_faces_are_inconsistent() {
        let (vertices, mut faces) = grid(1);
        // The same triangle twice with the same winding shares its edges in the same direction.
        faces.push(faces[0]);
        assert_eq!(count_inverted_faces(&vertices, &faces), 1);
    }

    #[test]
    fn non_finite_coordinates_fall_back_to_the_minority_rule() {
        let (mut vertices, mut faces) = icosphere(1);
        vertices[3] = [f64::NAN, 0.0, 0.0];
        assert_eq!(count_inverted_faces(&vertices, &faces), 0);
        flip(&mut faces[10]);
        assert_eq!(count_inverted_faces(&vertices, &faces), 1);
    }
}
