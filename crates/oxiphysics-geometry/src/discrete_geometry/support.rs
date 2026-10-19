// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Building blocks shared by the discrete-geometry solvers.
//!
//! * [`EdgeTable`]: the unique edges of a triangle list with the edge opposite every
//!   triangle corner, the number of triangles at every edge, and manifold /
//!   orientation checks;
//! * [`Components`]: connected components of the vertices through the triangles;
//! * [`ScaledSpd`]: a Jacobi-scaled conjugate-gradient solver for symmetric
//!   positive definite systems in which some unknowns are pinned to zero.

use crate::error::{Error, Result};
use oxiphysics_core::sparse::{CsrMatrix, pcg_solve};

/// Relative residual (against the right-hand side) at which the refinement stops.
const RELATIVE_RESIDUAL: f64 = 1e-12;
/// Largest relative residual that is still accepted as a solution.
const ACCEPTED_RESIDUAL: f64 = 1e-8;
/// Tolerance handed to every conjugate-gradient pass; each pass works on a
/// right-hand side of unit norm.
const PASS_TOLERANCE: f64 = 1e-10;
/// Largest number of conjugate-gradient passes (one solve plus refinements).
const MAX_PASSES: usize = 8;

/// Euclidean norm of a slice.
fn l2_norm(values: &[f64]) -> f64 {
    values.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// The edges of a triangle list.
///
/// Edges are numbered in lexicographic order of their vertex pair `(low, high)`.
#[derive(Debug, Clone)]
pub(super) struct EdgeTable {
    /// Unique undirected edges as `(low, high)` vertex pairs, sorted.
    pub(super) edges: Vec<(usize, usize)>,
    /// For every triangle the edge *opposite* each of its corners: entry `c` joins
    /// corners `c + 1` and `c + 2` (modulo three).
    pub(super) face_edges: Vec<[usize; 3]>,
    /// Number of triangles incident to every edge (one for a boundary edge).
    pub(super) face_count: Vec<usize>,
    /// `true` if every edge shared by two triangles is traversed in opposite
    /// directions by them, i.e. the triangle orientations are consistent.
    pub(super) consistently_oriented: bool,
}

impl EdgeTable {
    /// Build the table.
    ///
    /// # Errors
    /// [`Error::InvalidMesh`] if a triangle references a vertex `>= vertex_count`
    /// or if an edge is shared by more than two triangles.
    pub(super) fn build(vertex_count: usize, faces: &[[usize; 3]]) -> Result<Self> {
        for (f, tri) in faces.iter().enumerate() {
            if let Some(&bad) = tri.iter().find(|&&v| v >= vertex_count) {
                return Err(Error::invalid_mesh(format!(
                    "triangle {f} references vertex {bad} but the mesh has {vertex_count} vertices"
                )));
            }
        }
        // One record per triangle side: (low, high, slot, forward) where
        // `slot = 3 * triangle + corner` and `forward` tells whether the counter-
        // clockwise traversal of the side runs from the low to the high vertex.
        let mut sides: Vec<(usize, usize, usize, bool)> = Vec::with_capacity(3 * faces.len());
        for (f, tri) in faces.iter().enumerate() {
            for c in 0..3 {
                let from = tri[(c + 1) % 3];
                let to = tri[(c + 2) % 3];
                sides.push((from.min(to), from.max(to), 3 * f + c, from <= to));
            }
        }
        sides.sort_unstable_by_key(|&(low, high, slot, _)| (low, high, slot));

        let mut edges: Vec<(usize, usize)> = Vec::new();
        let mut face_count: Vec<usize> = Vec::new();
        let mut face_edges = vec![[0usize; 3]; faces.len()];
        let mut consistently_oriented = true;
        let mut start = 0;
        while start < sides.len() {
            let (low, high, ..) = sides[start];
            let mut end = start + 1;
            while end < sides.len() && sides[end].0 == low && sides[end].1 == high {
                end += 1;
            }
            let count = end - start;
            if count > 2 {
                return Err(Error::invalid_mesh(format!(
                    "edge ({low}, {high}) is shared by {count} triangles; \
                     only manifold meshes are supported"
                )));
            }
            if count == 2 && sides[start].3 == sides[start + 1].3 {
                consistently_oriented = false;
            }
            let id = edges.len();
            edges.push((low, high));
            face_count.push(count);
            for side in &sides[start..end] {
                face_edges[side.2 / 3][side.2 % 3] = id;
            }
            start = end;
        }
        Ok(Self {
            edges,
            face_edges,
            face_count,
            consistently_oriented,
        })
    }

    /// Vertices that lie on an edge with a single incident triangle.
    pub(super) fn boundary_vertices(&self, vertex_count: usize) -> Vec<bool> {
        let mut boundary = vec![false; vertex_count];
        for (&(low, high), &count) in self.edges.iter().zip(self.face_count.iter()) {
            if count == 1 {
                boundary[low] = true;
                boundary[high] = true;
            }
        }
        boundary
    }
}

/// Connected components of the vertices through a triangle list.
#[derive(Debug, Clone)]
pub(super) struct Components {
    /// Component of every vertex, numbered by first appearance (a vertex that is in
    /// no triangle forms a component of its own).
    pub(super) id: Vec<usize>,
    /// Number of components.
    pub(super) count: usize,
}

impl Components {
    /// Label the components; every index in `faces` must be `< vertex_count`.
    pub(super) fn label(vertex_count: usize, faces: &[[usize; 3]]) -> Self {
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        let mut parent: Vec<usize> = (0..vertex_count).collect();
        for tri in faces {
            for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2])] {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    parent[ra.max(rb)] = ra.min(rb);
                }
            }
        }
        let mut label_of_root = vec![usize::MAX; vertex_count];
        let mut id = vec![0usize; vertex_count];
        let mut count = 0usize;
        for (v, slot) in id.iter_mut().enumerate() {
            let root = find(&mut parent, v);
            if label_of_root[root] == usize::MAX {
                label_of_root[root] = count;
                count += 1;
            }
            *slot = label_of_root[root];
        }
        Self { id, count }
    }
}

/// A symmetric positive semi-definite system scaled to a unit diagonal, with some
/// unknowns pinned to zero.
///
/// The system is given as a list of `(row, column, value)` entries (repeated
/// positions are summed).  Pinned unknowns, and unknowns whose diagonal entry is
/// not positive (isolated vertices), have their row and column replaced by those
/// of the identity and a zero right-hand side, so the remaining matrix is
/// symmetric positive definite as soon as every connected block contains a pinned
/// unknown.  The matrix is scaled symmetrically by `S = diag(A)^-1/2`; conjugate
/// gradients on `S A S` generate exactly the iterates of Jacobi-preconditioned
/// conjugate gradients on `A`, but the absolute thresholds of the routine of
/// `oxiphysics-core` then act on a well-scaled problem.  The first solve is
/// followed by iterative refinement on the true residual until it is below
/// `1e-12` of the right-hand side.
pub(super) struct ScaledSpd {
    /// `S A S` with identity rows for the pinned unknowns.
    matrix: CsrMatrix,
    /// The diagonal of `S` (one for pinned unknowns).
    scale: Vec<f64>,
    /// Pinned unknowns.
    fixed: Vec<bool>,
}

impl ScaledSpd {
    /// Assemble the system of `n` unknowns from `entries`; `pinned[i]` pins unknown
    /// `i` to zero.  Every index in `entries` must be `< n` and `pinned` must have
    /// length `n`.
    pub(super) fn new(n: usize, entries: &[(usize, usize, f64)], pinned: &[bool]) -> Self {
        let mut diagonal = vec![0.0_f64; n];
        for &(row, col, value) in entries {
            if row == col && !pinned[row] {
                diagonal[row] += value;
            }
        }
        let fixed: Vec<bool> = diagonal
            .iter()
            .zip(pinned.iter())
            .map(|(&d, &p)| p || !(d.is_finite() && d > 0.0))
            .collect();
        let scale: Vec<f64> = diagonal
            .iter()
            .zip(fixed.iter())
            .map(|(&d, &is_fixed)| if is_fixed { 1.0 } else { 1.0 / d.sqrt() })
            .collect();
        let mut rows = Vec::with_capacity(entries.len() + n);
        let mut cols = Vec::with_capacity(entries.len() + n);
        let mut vals = Vec::with_capacity(entries.len() + n);
        for &(row, col, value) in entries {
            if !fixed[row] && !fixed[col] {
                rows.push(row);
                cols.push(col);
                vals.push(value * scale[row] * scale[col]);
            }
        }
        for (i, &is_fixed) in fixed.iter().enumerate() {
            if is_fixed {
                rows.push(i);
                cols.push(i);
                vals.push(1.0);
            }
        }
        Self {
            matrix: CsrMatrix::from_triplets(n, n, &rows, &cols, &vals),
            scale,
            fixed,
        }
    }

    /// Solve `A x = rhs`; the entries of `rhs` at pinned unknowns are ignored and
    /// the corresponding entries of `x` are zero.
    ///
    /// # Errors
    /// * [`Error::LengthMismatch`] if `rhs` has the wrong length;
    /// * [`Error::General`] for a non-finite right-hand side;
    /// * [`Error::ConvergenceFailure`] if the relative residual stays above `1e-8`
    ///   (a singular system, i.e. a block without a pinned unknown and an
    ///   inconsistent right-hand side).
    pub(super) fn solve(&self, rhs: &[f64], operation: &'static str) -> Result<Vec<f64>> {
        let n = self.scale.len();
        if rhs.len() != n {
            return Err(Error::length_mismatch(n, rhs.len()));
        }
        let scaled_rhs: Vec<f64> = (0..n)
            .map(|i| {
                if self.fixed[i] {
                    0.0
                } else {
                    rhs[i] * self.scale[i]
                }
            })
            .collect();
        let norm = l2_norm(&scaled_rhs);
        if !norm.is_finite() {
            return Err(Error::general(format!(
                "{operation}: non-finite right-hand side"
            )));
        }
        if norm == 0.0 {
            return Ok(vec![0.0; n]);
        }
        let unit_preconditioner = vec![1.0_f64; n];
        let max_iter = 2 * n + 1000;
        let mut x = vec![0.0_f64; n];
        let mut residual = scaled_rhs.clone();
        let mut residual_norm = norm;
        let mut iterations = 0usize;
        for _ in 0..MAX_PASSES {
            let unit_rhs: Vec<f64> = residual.iter().map(|r| r / residual_norm).collect();
            let (step, passes, _) = pcg_solve(
                &self.matrix,
                &unit_rhs,
                &unit_preconditioner,
                PASS_TOLERANCE,
                max_iter,
            );
            iterations += passes as usize;
            for (xi, di) in x.iter_mut().zip(step.iter()) {
                *xi += residual_norm * di;
            }
            let applied = self.matrix.matvec(&x);
            residual = scaled_rhs
                .iter()
                .zip(applied.iter())
                .map(|(b, a)| b - a)
                .collect();
            let new_norm = l2_norm(&residual);
            if !new_norm.is_finite() {
                return Err(Error::convergence_failure(
                    operation,
                    iterations,
                    f64::INFINITY,
                ));
            }
            let stagnated = new_norm > 0.5 * residual_norm;
            residual_norm = new_norm;
            if residual_norm <= RELATIVE_RESIDUAL * norm || stagnated {
                break;
            }
        }
        if residual_norm > ACCEPTED_RESIDUAL * norm {
            return Err(Error::convergence_failure(
                operation,
                iterations,
                residual_norm / norm,
            ));
        }
        Ok(x.iter()
            .zip(self.scale.iter())
            .map(|(v, s)| v * s)
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two triangles forming a square, counter-clockwise.
    fn two_triangles() -> Vec<[usize; 3]> {
        vec![[0, 1, 2], [0, 2, 3]]
    }

    /// Boundary of a tetrahedron, outward counter-clockwise.
    fn tetrahedron_faces() -> Vec<[usize; 3]> {
        vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
    }

    #[test]
    fn edge_table_of_a_closed_surface() {
        let table = EdgeTable::build(4, &tetrahedron_faces()).expect("valid mesh");
        assert_eq!(table.edges.len(), 6);
        assert!(table.face_count.iter().all(|&c| c == 2));
        assert!(table.consistently_oriented);
        assert!(table.boundary_vertices(4).iter().all(|&b| !b));
        // The edges are sorted and every triangle lists the edge opposite each corner.
        assert!(table.edges.windows(2).all(|w| w[0] < w[1]));
        for (tri, ids) in tetrahedron_faces().iter().zip(table.face_edges.iter()) {
            for c in 0..3 {
                let (a, b) = (tri[(c + 1) % 3], tri[(c + 2) % 3]);
                assert_eq!(table.edges[ids[c]], (a.min(b), a.max(b)));
            }
        }
    }

    #[test]
    fn edge_table_finds_the_boundary() {
        let table = EdgeTable::build(4, &two_triangles()).expect("valid mesh");
        assert_eq!(table.edges, vec![(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)]);
        assert_eq!(table.face_count, vec![1, 2, 1, 1, 1]);
        assert!(table.consistently_oriented);
        assert_eq!(
            table.boundary_vertices(5),
            vec![true, true, true, true, false]
        );
    }

    #[test]
    fn edge_table_detects_inconsistent_orientation() {
        let table = EdgeTable::build(4, &[[0, 1, 2], [0, 3, 2]]).expect("valid mesh");
        assert!(!table.consistently_oriented);
    }

    #[test]
    fn edge_table_rejects_non_manifold_edges_and_bad_indices() {
        let fan = [[0, 1, 2], [0, 1, 3], [0, 1, 4]];
        assert!(matches!(
            EdgeTable::build(5, &fan),
            Err(Error::InvalidMesh { .. })
        ));
        assert!(matches!(
            EdgeTable::build(3, &[[0, 1, 3]]),
            Err(Error::InvalidMesh { .. })
        ));
    }

    #[test]
    fn components_are_numbered_by_first_appearance() {
        let faces = [[0, 1, 2], [3, 4, 5], [2, 1, 6]];
        let components = Components::label(8, &faces);
        assert_eq!(components.count, 3);
        assert_eq!(components.id, vec![0, 0, 0, 1, 1, 1, 0, 2]);
    }

    /// The path Laplacian with a Dirichlet end: `A = tridiag(-1, 2, -1)` with the
    /// last unknown eliminated, whose solution for `b = e` is known in closed form.
    #[test]
    fn solves_a_dirichlet_path_system() {
        let n = 200;
        let mut entries = Vec::new();
        for i in 0..n {
            entries.push((i, i, if i == 0 || i == n - 1 { 1.0 } else { 2.0 }));
            if i + 1 < n {
                entries.push((i, i + 1, -1.0));
                entries.push((i + 1, i, -1.0));
            }
        }
        let mut pinned = vec![false; n];
        pinned[n - 1] = true;
        let system = ScaledSpd::new(n, &entries, &pinned);
        // Unit load at the free end: x_i = n - 1 - i (a straight ramp to the pin).
        let mut rhs = vec![0.0; n];
        rhs[0] = 1.0;
        let x = system.solve(&rhs, "test").expect("converges");
        for (i, xi) in x.iter().enumerate() {
            let exact = (n - 1 - i) as f64;
            assert!((xi - exact).abs() < 1e-8 * (1.0 + exact), "x[{i}] = {xi}");
        }
    }

    #[test]
    fn pinned_and_isolated_unknowns_are_zero() {
        // Unknown 2 has no entries at all (an isolated vertex), unknown 0 is pinned.
        let entries = [
            (0, 0, 1.0),
            (0, 1, -1.0),
            (1, 0, -1.0),
            (1, 1, 1.0),
            (3, 3, 4.0),
        ];
        let system = ScaledSpd::new(4, &entries, &[true, false, false, false]);
        let x = system
            .solve(&[5.0, 3.0, 7.0, 8.0], "test")
            .expect("converges");
        assert_eq!(x[0], 0.0);
        assert!((x[1] - 3.0).abs() < 1e-10, "x[1] = {}", x[1]);
        assert_eq!(x[2], 0.0);
        assert!((x[3] - 2.0).abs() < 1e-10, "x[3] = {}", x[3]);
    }

    #[test]
    fn refinement_reaches_high_accuracy() {
        // A 2-D grid Laplacian with Dirichlet rows at one corner.
        let m = 24;
        let n = m * m;
        let mut entries = Vec::new();
        for i in 0..m {
            for j in 0..m {
                let v = i * m + j;
                entries.push((v, v, 4.0));
                if i + 1 < m {
                    entries.push((v, v + m, -1.0));
                    entries.push((v + m, v, -1.0));
                }
                if j + 1 < m {
                    entries.push((v, v + 1, -1.0));
                    entries.push((v + 1, v, -1.0));
                }
            }
        }
        let mut pinned = vec![false; n];
        pinned[0] = true;
        let system = ScaledSpd::new(n, &entries, &pinned);
        let rhs: Vec<f64> = (0..n).map(|v| ((v * 7919) % 101) as f64 - 50.0).collect();
        let x = system.solve(&rhs, "test").expect("converges");
        let mut ax = vec![0.0; n];
        for &(row, col, value) in &entries {
            if row != 0 && col != 0 {
                ax[row] += value * x[col];
            }
        }
        for v in 1..n {
            assert!(
                (ax[v] - rhs[v]).abs() < 1e-9,
                "row {v}: {} vs {}",
                ax[v],
                rhs[v]
            );
        }
    }

    #[test]
    fn singular_inconsistent_system_is_reported() {
        // The pure Neumann path Laplacian with a right-hand side that does not sum to zero.
        let n = 10;
        let mut entries = Vec::new();
        for i in 0..n {
            entries.push((i, i, if i == 0 || i == n - 1 { 1.0 } else { 2.0 }));
            if i + 1 < n {
                entries.push((i, i + 1, -1.0));
                entries.push((i + 1, i, -1.0));
            }
        }
        let system = ScaledSpd::new(n, &entries, &vec![false; n]);
        let rhs = vec![1.0; n];
        assert!(matches!(
            system.solve(&rhs, "singular"),
            Err(Error::ConvergenceFailure { .. })
        ));
        assert!(matches!(
            system.solve(&rhs[..3], "short"),
            Err(Error::LengthMismatch { .. })
        ));
        let mut nan_rhs = rhs.clone();
        nan_rhs[2] = f64::NAN;
        assert!(matches!(
            system.solve(&nan_rhs, "nan"),
            Err(Error::General(_))
        ));
        assert_eq!(
            system.solve(&vec![0.0; n], "zero").expect("zero"),
            vec![0.0; n]
        );
    }
}
