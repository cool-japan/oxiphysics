// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Convex-combination (Tutte / harmonic) embeddings of a disc-like patch.
//!
//! The longest boundary loop is pinned to the unit circle in proportion to its arc
//! length and every other vertex is the weighted barycentre of its neighbours.
//! With strictly positive weights and a convex boundary this is a bijection onto
//! the disc (Tutte 1963, Floater 2003), so no UV triangle folds.  The linear
//! system is assembled once as a sparse matrix and solved with preconditioned
//! conjugate gradients.

use super::ParamTriMesh;
use super::boundary::{analyse_connectivity, boundary_loops, require_single_component};
use super::solve::solve_spd;
use crate::error::{Error, Result};
use oxiphysics_core::sparse::CsrMatrix;
use std::collections::BTreeMap;

/// Smallest edge weight of the cotangent scheme.
///
/// Cotangent weights are negative across obtuse angles, which would break the
/// no-fold guarantee.  Clamping to a tiny positive value keeps every weight
/// positive; the weights are dimensionless, so the floor is meaningful for any
/// mesh scale.
const MIN_COTANGENT_WEIGHT: f64 = 1.0e-8;

/// How the weight of an edge is derived from the mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum EdgeWeights {
    /// Every edge weighs 1: Tutte's barycentric embedding.
    Uniform,
    /// `cot(alpha) + cot(beta)` of the two angles opposite the edge, clamped to a
    /// small positive value: the discrete harmonic map.
    Cotangent,
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross_length(a: [f64; 3], b: [f64; 3]) -> f64 {
    let c = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    dot3(c, c).sqrt()
}

/// Weight of every undirected mesh edge `(low, high)`, in ascending edge order.
fn edge_weights(mesh: &ParamTriMesh, scheme: EdgeWeights) -> BTreeMap<(usize, usize), f64> {
    let mut weights: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    for tri in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            if a != b {
                let weight = weights.entry((a.min(b), a.max(b))).or_insert(0.0);
                if scheme == EdgeWeights::Uniform {
                    *weight = 1.0;
                }
            }
        }
        if scheme != EdgeWeights::Cotangent {
            continue;
        }
        // The cotangent of the angle at each corner belongs to the opposite edge.
        for k in 0..3 {
            let (corner, p, q) = (tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]);
            if corner == p || corner == q || p == q {
                continue;
            }
            let u = sub3(mesh.positions[p], mesh.positions[corner]);
            let v = sub3(mesh.positions[q], mesh.positions[corner]);
            let area2 = cross_length(u, v);
            let lengths = dot3(u, u).sqrt() * dot3(v, v).sqrt();
            // Skip slivers whose angle is numerically 0 or pi: its cotangent is
            // meaningless and would swamp the system.
            if area2 > 1.0e-12 * lengths
                && let Some(weight) = weights.get_mut(&(p.min(q), p.max(q)))
            {
                *weight += dot3(u, v) / area2;
            }
        }
    }
    if scheme == EdgeWeights::Cotangent {
        for weight in weights.values_mut() {
            *weight = weight.max(MIN_COTANGENT_WEIGHT);
        }
    }
    weights
}

/// Embed the patch `mesh` in the unit disc with the given edge weights.
///
/// Vertices that no triangle uses receive `[0, 0]`.
///
/// # Errors
/// * [`Error::IndexOutOfBounds`] for a triangle that refers to a missing vertex;
/// * [`Error::InvalidMesh`] for a mesh without triangles, with several connected
///   components, or without a boundary (a closed surface cannot be mapped into
///   the disc);
/// * [`Error::ConvergenceFailure`] if the linear solver stalls.
pub(super) fn circle_embedding(mesh: &ParamTriMesh, scheme: EdgeWeights) -> Result<Vec<[f64; 2]>> {
    let n = mesh.vertex_count();
    if n == 0 {
        return Ok(Vec::new());
    }
    let connectivity = analyse_connectivity(n, &mesh.triangles)?;
    require_single_component(&connectivity)?;
    let loops = boundary_loops(&mesh.triangles);
    let Some(outer) = loops.first() else {
        return Err(Error::invalid_mesh(
            "a closed mesh has no boundary to map onto the circle; cut a seam first",
        ));
    };
    let mut uvs =
        super::TutteParameterization::map_boundary_loop_to_circle(&mesh.positions, outer, n);

    let mut pinned = vec![false; n];
    for &v in outer {
        pinned[v] = true;
    }
    // Free vertices in ascending id order; `free_of[v]` is the row of vertex `v`.
    let mut free_of = vec![usize::MAX; n];
    let mut free_vertices: Vec<usize> = Vec::new();
    for (v, slot) in free_of.iter_mut().enumerate() {
        if connectivity.referenced[v] && !pinned[v] {
            *slot = free_vertices.len();
            free_vertices.push(v);
        }
    }
    if free_vertices.is_empty() {
        return Ok(uvs);
    }

    // Sorted neighbour lists with their weights.
    let mut neighbours: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for (&(a, b), &weight) in &edge_weights(mesh, scheme) {
        neighbours[a].push((b, weight));
        neighbours[b].push((a, weight));
    }

    // Row `r` of the Laplacian:  (sum_j w_ij) x_i - sum_{j free} w_ij x_j  =
    //                            sum_{j pinned} w_ij x_j.
    let count = free_vertices.len();
    let mut row_ptr = Vec::with_capacity(count + 1);
    let mut col_idx: Vec<usize> = Vec::new();
    let mut values: Vec<f64> = Vec::new();
    let mut rhs = [vec![0.0; count], vec![0.0; count]];
    row_ptr.push(0);
    for (row, &v) in free_vertices.iter().enumerate() {
        let mut entries: Vec<(usize, f64)> = Vec::with_capacity(neighbours[v].len() + 1);
        let mut diagonal = 0.0;
        for &(other, weight) in &neighbours[v] {
            diagonal += weight;
            if free_of[other] == usize::MAX {
                rhs[0][row] += weight * uvs[other][0];
                rhs[1][row] += weight * uvs[other][1];
            } else {
                entries.push((free_of[other], -weight));
            }
        }
        entries.push((row, diagonal));
        entries.sort_unstable_by_key(|&(col, _)| col);
        for (col, value) in entries {
            col_idx.push(col);
            values.push(value);
        }
        row_ptr.push(col_idx.len());
    }
    let laplacian = CsrMatrix {
        nrows: count,
        ncols: count,
        row_ptr,
        col_idx,
        values,
    };
    for (axis, b) in rhs.iter().enumerate() {
        let solution = solve_spd(&laplacian, b, "circle embedding")?;
        for (&v, &x) in free_vertices.iter().zip(solution.iter()) {
            uvs[v][axis] = x;
        }
    }
    Ok(uvs)
}
