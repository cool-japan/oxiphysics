// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Edge-point stencils of the butterfly subdivision scheme.
//!
//! Every stencil is a set of weights that sum to one, so each new vertex is an
//! affine combination of existing vertices and the scheme commutes with
//! translations, rotations and uniform scalings of the input mesh.
//!
//! The adjacency is undirected (edge to incident faces), so the stencils do not
//! depend on the winding of the input faces.

use std::collections::HashMap;

use super::{add3, midpoint3, scale3};

/// Weight of each end point of the edge in the regular eight-point stencil.
const REGULAR_END_WEIGHT: f64 = 0.5;
/// Weight of the two vertices opposite to the edge in the regular stencil.
const REGULAR_OPPOSITE_WEIGHT: f64 = 1.0 / 8.0;
/// Weight of each of the four wing vertices in the regular stencil.
const REGULAR_WING_WEIGHT: f64 = -1.0 / 16.0;
/// Weight of each end point of a boundary edge in the 4-point boundary rule.
const BOUNDARY_END_WEIGHT: f64 = 9.0 / 16.0;
/// Weight of the two outer boundary neighbours in the 4-point boundary rule.
const BOUNDARY_OUTER_WEIGHT: f64 = -1.0 / 16.0;
/// Weight of the extraordinary vertex itself in the modified butterfly stencil.
const EXTRAORDINARY_CENTRE_WEIGHT: f64 = 0.75;
/// Valence of a regular interior vertex of a triangle mesh.
const REGULAR_VALENCE: usize = 6;

/// Undirected edge key with the smaller vertex index first.
fn edge_key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// Edge-face adjacency of a triangle mesh.
pub(super) struct EdgeAdjacency<'a> {
    faces: &'a [[usize; 3]],
    /// Undirected edge to the faces that contain it, in face order.
    edge_faces: HashMap<(usize, usize), Vec<usize>>,
    /// Sorted, duplicate-free neighbour vertices of every vertex.
    neighbours: Vec<Vec<usize>>,
}

impl<'a> EdgeAdjacency<'a> {
    /// Build the adjacency of `faces` over `vertex_count` vertices.
    pub(super) fn new(vertex_count: usize, faces: &'a [[usize; 3]]) -> Self {
        let mut edge_faces: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); vertex_count];
        for (fi, &[a, b, c]) in faces.iter().enumerate() {
            for (u, v) in [(a, b), (b, c), (c, a)] {
                edge_faces.entry(edge_key(u, v)).or_default().push(fi);
                neighbours[u].push(v);
                neighbours[v].push(u);
            }
        }
        for list in &mut neighbours {
            list.sort_unstable();
            list.dedup();
        }
        Self {
            faces,
            edge_faces,
            neighbours,
        }
    }

    /// Faces that contain the edge `{u, v}`.
    fn faces_of(&self, u: usize, v: usize) -> &[usize] {
        self.edge_faces
            .get(&edge_key(u, v))
            .map_or(&[][..], Vec::as_slice)
    }

    /// Vertex of `face` that is neither `u` nor `v`.
    fn third_vertex(&self, face: usize, u: usize, v: usize) -> Option<usize> {
        self.faces[face].iter().copied().find(|&w| w != u && w != v)
    }

    /// Vertex opposite to the edge `{u, v}` in the face on the other side of
    /// `face`. `None` unless the edge is shared by exactly two faces, one of them
    /// `face`.
    fn across(&self, u: usize, v: usize, face: usize) -> Option<usize> {
        let &[first, second] = self.faces_of(u, v) else {
            return None;
        };
        let other = if first == face {
            second
        } else if second == face {
            first
        } else {
            return None;
        };
        self.third_vertex(other, u, v)
    }

    /// Neighbours of the interior vertex `v` in cyclic order, starting with
    /// `start`. `None` when `v` is not an interior vertex with a single closed
    /// ring of faces (boundary vertex, non-manifold edge, pinched vertex).
    fn interior_ring(&self, v: usize, start: usize) -> Option<Vec<usize>> {
        let &[first_face, _] = self.faces_of(v, start) else {
            return None;
        };
        let neighbour_count = self.neighbours[v].len();
        let mut ring = vec![start];
        let mut face = first_face;
        let mut current = self.third_vertex(face, v, start)?;
        while current != start {
            if ring.len() >= neighbour_count {
                return None;
            }
            ring.push(current);
            let &[f0, f1] = self.faces_of(v, current) else {
                return None;
            };
            let next_face = if f0 == face {
                f1
            } else if f1 == face {
                f0
            } else {
                return None;
            };
            current = self.third_vertex(next_face, v, current)?;
            face = next_face;
        }
        (ring.len() == neighbour_count).then_some(ring)
    }

    /// The second boundary neighbour of the boundary vertex `v`, given the first
    /// one (`known`). `None` unless `v` has exactly two boundary edges.
    fn other_boundary_neighbour(&self, v: usize, known: usize) -> Option<usize> {
        let mut others = self.neighbours[v]
            .iter()
            .copied()
            .filter(|&u| u != known && self.faces_of(v, u).len() == 1);
        let found = others.next()?;
        others.next().is_none().then_some(found)
    }

    /// Position of the new vertex on the edge `{a, b}`.
    ///
    /// * Interior edge: the regular eight-point stencil, or the modified
    ///   butterfly stencil when an end point is an extraordinary vertex.
    /// * Boundary edge: the four-point rule along the boundary curve.
    /// * Anything else (non-manifold edge, missing wing vertices): the midpoint.
    pub(super) fn edge_point(&self, vertices: &[[f64; 3]], a: usize, b: usize) -> [f64; 3] {
        let stencil = match *self.faces_of(a, b) {
            [f1, f2] => self.interior_edge_point(vertices, a, b, f1, f2),
            [_] => self.boundary_edge_point(vertices, a, b),
            _ => None,
        };
        stencil.unwrap_or_else(|| midpoint3(vertices[a], vertices[b]))
    }

    fn interior_edge_point(
        &self,
        vertices: &[[f64; 3]],
        a: usize,
        b: usize,
        f1: usize,
        f2: usize,
    ) -> Option<[f64; 3]> {
        let c = self.third_vertex(f1, a, b)?;
        let d = self.third_vertex(f2, a, b)?;
        if c == d {
            return None;
        }
        let extraordinary: Vec<[f64; 3]> = [(a, b), (b, a)]
            .into_iter()
            .filter_map(|(vertex, other)| {
                let ring = self.interior_ring(vertex, other)?;
                (ring.len() != REGULAR_VALENCE)
                    .then(|| extraordinary_point(vertices, vertex, &ring))
                    .flatten()
            })
            .collect();
        if !extraordinary.is_empty() {
            let sum = extraordinary
                .iter()
                .fold([0.0; 3], |acc, &point| add3(acc, point));
            return Some(scale3(sum, 1.0 / extraordinary.len() as f64));
        }

        // Wings: the vertices opposite to the four outer edges of the two faces
        // `abc` and `abd`.
        let wings = [
            self.across(c, b, f1)?,
            self.across(a, c, f1)?,
            self.across(d, a, f2)?,
            self.across(b, d, f2)?,
        ];
        let ends = scale3(add3(vertices[a], vertices[b]), REGULAR_END_WEIGHT);
        let opposite = scale3(add3(vertices[c], vertices[d]), REGULAR_OPPOSITE_WEIGHT);
        let wing_sum = wings
            .iter()
            .fold([0.0; 3], |acc, &w| add3(acc, vertices[w]));
        Some(add3(
            add3(ends, opposite),
            scale3(wing_sum, REGULAR_WING_WEIGHT),
        ))
    }

    fn boundary_edge_point(&self, vertices: &[[f64; 3]], a: usize, b: usize) -> Option<[f64; 3]> {
        let before = self.other_boundary_neighbour(a, b)?;
        let after = self.other_boundary_neighbour(b, a)?;
        let inner = scale3(add3(vertices[a], vertices[b]), BOUNDARY_END_WEIGHT);
        let outer = scale3(
            add3(vertices[before], vertices[after]),
            BOUNDARY_OUTER_WEIGHT,
        );
        Some(add3(inner, outer))
    }
}

/// Modified butterfly weights of the ring vertices of an extraordinary vertex
/// of the given valence (Zorin, Schroeder & Sweldens 1996); the first weight
/// belongs to the other end point of the edge. They sum to `1/4`, and the
/// vertex itself carries [`EXTRAORDINARY_CENTRE_WEIGHT`]. `None` for valences
/// below three.
fn extraordinary_weights(valence: usize) -> Option<Vec<f64>> {
    match valence {
        0..=2 => None,
        3 => Some(vec![5.0 / 12.0, -1.0 / 12.0, -1.0 / 12.0]),
        4 => Some(vec![3.0 / 8.0, 0.0, -1.0 / 8.0, 0.0]),
        k => {
            let count = k as f64;
            Some(
                (0..k)
                    .map(|j| {
                        let angle = 2.0 * std::f64::consts::PI * j as f64 / count;
                        (0.25 + angle.cos() + 0.5 * (2.0 * angle).cos()) / count
                    })
                    .collect(),
            )
        }
    }
}

/// Edge point from the one-ring of the extraordinary vertex `vertex`.
fn extraordinary_point(vertices: &[[f64; 3]], vertex: usize, ring: &[usize]) -> Option<[f64; 3]> {
    let weights = extraordinary_weights(ring.len())?;
    let mut point = scale3(vertices[vertex], EXTRAORDINARY_CENTRE_WEIGHT);
    for (&neighbour, &weight) in ring.iter().zip(&weights) {
        point = add3(point, scale3(vertices[neighbour], weight));
    }
    Some(point)
}
