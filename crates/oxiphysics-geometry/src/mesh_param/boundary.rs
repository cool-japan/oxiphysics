// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ordered boundary loops and connectivity checks for the parameterization
//! algorithms.
//!
//! A parameterization that pins the boundary of a patch (Tutte, harmonic map) must
//! visit the boundary vertices in the order in which they follow each other around
//! the patch; vertex-id order is meaningless.  [`boundary_loops`] extracts that
//! order deterministically from the triangle list.

use crate::error::{Error, Result};
use std::collections::{BTreeMap, HashMap};

/// Closed boundary loops of a triangle mesh, longest first.
///
/// A boundary edge is an undirected edge that belongs to exactly one triangle.  It
/// is directed as in that triangle, so for consistently wound triangles every loop
/// is traversed with the mesh on its left (counter-clockwise around the outer
/// boundary of a patch seen from the side its triangles face).  Each loop lists
/// every vertex once, starts at its smallest vertex id and is closed implicitly
/// (the last vertex connects back to the first).  Loops are ordered by decreasing
/// vertex count, ties by their smallest vertex id.
///
/// The result does not depend on hash-map iteration order.  At a boundary vertex
/// that is pinched between several fans of triangles (a bow-tie) the walk follows
/// the boundary edge with the smallest head first and splits the trail into one
/// loop per closed cycle, so no loop repeats a vertex.  Boundary edges that do
/// not close into a cycle (inconsistently wound input) are ignored.  A closed
/// mesh has no loops.
pub fn boundary_loops(triangles: &[[usize; 3]]) -> Vec<Vec<usize>> {
    // Triangles incident to every undirected edge.
    let mut incidence: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in triangles {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            if a != b {
                *incidence.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
    }

    // Boundary edges grouped by tail vertex.  Each list is sorted in descending
    // order so that `pop` yields the smallest head first.
    let mut outgoing: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for tri in triangles {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            if a != b && incidence.get(&(a.min(b), a.max(b))) == Some(&1) {
                outgoing.entry(a).or_default().push(b);
            }
        }
    }
    for heads in outgoing.values_mut() {
        heads.sort_unstable_by(|x, y| y.cmp(x));
    }

    let tails: Vec<usize> = outgoing.keys().copied().collect();
    let mut loops: Vec<Vec<usize>> = Vec::new();
    for tail in tails {
        while outgoing.get(&tail).is_some_and(|heads| !heads.is_empty()) {
            extract_cycles(tail, &mut outgoing, &mut loops);
        }
    }

    let mut loops: Vec<Vec<usize>> = loops
        .into_iter()
        .filter(|cycle| cycle.len() >= 3)
        .map(|mut cycle| {
            let first = cycle
                .iter()
                .enumerate()
                .min_by_key(|&(_, &v)| v)
                .map_or(0, |(at, _)| at);
            cycle.rotate_left(first);
            cycle
        })
        .collect();
    loops.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a[0].cmp(&b[0])));
    loops
}

/// Follow unused boundary edges from `start` and append every cycle that closes
/// along the way to `loops`.
///
/// The trail of visited vertices is kept without repeats: as soon as an edge
/// leads back to a vertex of the trail, the cycle it closes is cut out and the
/// walk resumes from that vertex.  A vertex without unused outgoing edges is a
/// dead end (open or inconsistently wound boundary) and is dropped from the
/// trail.  Every step consumes an edge or removes a vertex, so the walk ends.
fn extract_cycles(
    start: usize,
    outgoing: &mut BTreeMap<usize, Vec<usize>>,
    loops: &mut Vec<Vec<usize>>,
) {
    let mut trail = vec![start];
    let mut position: HashMap<usize, usize> = HashMap::from([(start, 0)]);
    while let Some(&current) = trail.last() {
        match outgoing.get_mut(&current).and_then(Vec::pop) {
            Some(next) => match position.get(&next).copied() {
                Some(at) => {
                    loops.push(trail[at..].to_vec());
                    for cut in trail.drain(at + 1..) {
                        position.remove(&cut);
                    }
                }
                None => {
                    position.insert(next, trail.len());
                    trail.push(next);
                }
            },
            None => {
                if let Some(dead_end) = trail.pop() {
                    position.remove(&dead_end);
                }
            }
        }
    }
}

/// The longest boundary loop of a triangle mesh (see [`boundary_loops`]).
///
/// Returns `None` for a closed mesh.  Walking the loop keeps the mesh on the
/// left, so mapping its vertices to increasing angles of a circle preserves the
/// orientation of the triangles.
pub fn boundary_loop(triangles: &[[usize; 3]]) -> Option<Vec<usize>> {
    boundary_loops(triangles).into_iter().next()
}

/// Which vertices a triangle list uses and how many connected pieces it forms.
pub(super) struct Connectivity {
    /// `true` for every vertex that belongs to at least one triangle.
    pub(super) referenced: Vec<bool>,
    /// Number of connected components of the triangle mesh.
    pub(super) components: usize,
}

/// Root of `v` in the union-find forest `parent` (with path halving).
fn find_root(parent: &mut [usize], mut v: usize) -> usize {
    while parent[v] != v {
        parent[v] = parent[parent[v]];
        v = parent[v];
    }
    v
}

/// Validate the triangle indices against `n_verts` and measure the connectivity.
///
/// # Errors
/// [`Error::IndexOutOfBounds`] when a triangle refers to a vertex `>= n_verts`.
pub(super) fn analyse_connectivity(
    n_verts: usize,
    triangles: &[[usize; 3]],
) -> Result<Connectivity> {
    let mut referenced = vec![false; n_verts];
    let mut parent: Vec<usize> = (0..n_verts).collect();
    for tri in triangles {
        for &v in tri {
            if v >= n_verts {
                return Err(Error::index_out_of_bounds(v, n_verts));
            }
            referenced[v] = true;
        }
        let root = find_root(&mut parent, tri[0]);
        for &v in &tri[1..] {
            let other = find_root(&mut parent, v);
            if other != root {
                parent[other] = root;
            }
        }
    }
    let components = (0..n_verts)
        .filter(|&v| referenced[v] && find_root(&mut parent, v) == v)
        .count();
    Ok(Connectivity {
        referenced,
        components,
    })
}

/// Require a non-empty, single-component triangle mesh.
///
/// # Errors
/// [`Error::InvalidMesh`] when there are no triangles or more than one connected
/// component (parameterize each chart separately and pack the results).
pub(super) fn require_single_component(connectivity: &Connectivity) -> Result<()> {
    match connectivity.components {
        0 => Err(Error::invalid_mesh("the mesh has no triangles")),
        1 => Ok(()),
        count => Err(Error::invalid_mesh(format!(
            "the mesh has {count} connected components; parameterize each chart separately"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Triangles of an `n x n` vertex grid, counter-clockwise in the xy plane,
    /// vertex `(i, j)` has id `j * n + i`.
    fn grid(n: usize) -> Vec<[usize; 3]> {
        let id = |i: usize, j: usize| j * n + i;
        let mut triangles = Vec::new();
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                triangles.push([id(i, j), id(i + 1, j), id(i + 1, j + 1)]);
                triangles.push([id(i, j), id(i + 1, j + 1), id(i, j + 1)]);
            }
        }
        triangles
    }

    #[test]
    fn grid_boundary_follows_the_loop_not_the_vertex_ids() {
        // Ascending ids would read 0,1,2,3,5,6,7,8; the loop is 0,1,2,5,8,7,6,3.
        assert_eq!(boundary_loop(&grid(3)), Some(vec![0, 1, 2, 5, 8, 7, 6, 3]));
        let loop5 = boundary_loop(&grid(5)).unwrap_or_default();
        assert_eq!(loop5.len(), 16);
        assert_eq!(&loop5[..5], &[0, 1, 2, 3, 4]);
        assert_eq!(&loop5[5..9], &[9, 14, 19, 24]);
        assert_eq!(&loop5[9..13], &[23, 22, 21, 20]);
        assert_eq!(&loop5[13..], &[15, 10, 5]);
    }

    #[test]
    fn reversing_the_winding_reverses_the_loop() {
        let flipped: Vec<[usize; 3]> = grid(3).iter().map(|&[a, b, c]| [a, c, b]).collect();
        // Same cycle, opposite direction, still starting at the smallest vertex.
        assert_eq!(boundary_loop(&flipped), Some(vec![0, 3, 6, 7, 8, 5, 2, 1]));
    }

    #[test]
    fn closed_meshes_have_no_boundary_loop() {
        // Tetrahedron with outward winding.
        let tetrahedron = [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];
        assert!(boundary_loops(&tetrahedron).is_empty());
        assert_eq!(boundary_loop(&tetrahedron), None);
        assert!(boundary_loops(&[]).is_empty());
    }

    #[test]
    fn annulus_reports_both_loops_longest_first() {
        // 4x4 vertex grid without the two triangles of the cell (1, 1); cell
        // (i, j) owns the triangles 2c and 2c + 1 with c = j * 3 + i.
        let n = 4;
        let hole_cell = (n - 1) + 1;
        let hole: Vec<[usize; 3]> = grid(n)
            .into_iter()
            .enumerate()
            .filter(|&(t, _)| t / 2 != hole_cell)
            .map(|(_, tri)| tri)
            .collect();
        let loops = boundary_loops(&hole);
        assert_eq!(loops.len(), 2);
        assert_eq!(loops[0].len(), 12);
        assert_eq!(loops[1].len(), 4);
        assert_eq!(loops[0][0], 0);
        assert_eq!(loops[1][0], 5);
        // The inner loop runs the other way round (the mesh stays on its left).
        assert_eq!(loops[1], vec![5, 9, 10, 6]);
    }

    #[test]
    fn bow_tie_splits_into_loops_without_repeated_vertices() {
        // Two triangles that touch only at vertex 0.
        let bow_tie = [[0, 1, 2], [0, 3, 4]];
        let loops = boundary_loops(&bow_tie);
        assert_eq!(loops, vec![vec![0, 1, 2], vec![0, 3, 4]]);
        for cycle in &loops {
            let mut sorted = cycle.clone();
            sorted.sort_unstable();
            sorted.dedup();
            assert_eq!(sorted.len(), cycle.len());
        }
    }

    #[test]
    fn inconsistent_or_non_manifold_input_terminates() {
        // One triangle flipped against its neighbour, plus an edge shared by three
        // triangles.
        let flipped = [[0, 1, 2], [1, 2, 3]];
        let _ = boundary_loops(&flipped);
        let fin = [[0, 1, 2], [0, 1, 3], [0, 1, 4]];
        let _ = boundary_loops(&fin);
        // Degenerate triangles with repeated vertices are ignored.
        assert!(boundary_loops(&[[0, 0, 1], [2, 2, 2]]).is_empty());
    }

    #[test]
    fn result_is_deterministic() {
        let reference = boundary_loops(&grid(7));
        for _ in 0..20 {
            assert_eq!(boundary_loops(&grid(7)), reference);
        }
    }

    #[test]
    fn connectivity_counts_components_and_validates_indices() {
        let two_pieces = [[0, 1, 2], [3, 4, 5]];
        let connectivity = analyse_connectivity(8, &two_pieces).unwrap();
        assert_eq!(connectivity.components, 2);
        assert_eq!(
            connectivity.referenced,
            vec![true, true, true, true, true, true, false, false]
        );
        assert!(require_single_component(&connectivity).is_err());

        let joined = [[0, 1, 2], [2, 3, 4]];
        let connectivity = analyse_connectivity(5, &joined).unwrap();
        assert_eq!(connectivity.components, 1);
        assert!(require_single_component(&connectivity).is_ok());

        let none = analyse_connectivity(3, &[]).unwrap();
        assert!(require_single_component(&none).is_err());
        assert!(matches!(
            analyse_connectivity(3, &[[0, 1, 3]]),
            Err(Error::IndexOutOfBounds { index: 3, len: 3 })
        ));
    }
}
