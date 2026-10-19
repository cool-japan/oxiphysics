// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Least Squares Conformal Maps (Lévy, Petitjean, Ray, Maillot, *Least squares
//! conformal maps for automatic texture atlas generation*, SIGGRAPH 2002).
//!
//! # Formulation
//!
//! Every triangle gets an orthonormal local frame `(e1, e2)` with the vertices at
//! `(0, 0)`, `(x1, 0)` and `(x2, y2)` (`y2 > 0`).  A map `f = u + i v` that is
//! linear on the triangle is conformal there exactly when it satisfies the
//! Cauchy–Riemann equations `(∂x + i ∂y) f = 0`, which in terms of the three
//! vertex values reads
//!
//! ```text
//!   ((x2 − x1) + i y2) f0  −  (x2 + i y2) f1  +  x1 f2  =  0 .
//! ```
//!
//! LSCM minimises the sum of the squared residuals of these complex rows, each
//! weighted by `1 / (2 · 2·area)`, which is the area-weighted conformal energy
//! `Σ_T A_T |(∂x + i ∂y) f|²`.  Two vertices are pinned to remove the similarity
//! ambiguity; the remaining free vertices solve the sparse symmetric positive
//! definite normal equations, assembled directly in CSR form and solved with
//! Jacobi-preconditioned conjugate gradients.

use super::ParamTriMesh;
use super::boundary::{analyse_connectivity, boundary_loops, require_single_component};
use super::solve::solve_spd;
use crate::error::{Error, Result};
use oxiphysics_core::sparse::CsrMatrix;

/// A triangle whose height over its first edge is below this fraction of the
/// product of its two edge lengths (an angle of about `1e-12` rad) is degenerate.
const DEGENERATE_SINE: f64 = 1.0e-12;

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Complex conformality row of a triangle: `Σ_j c_j f_j = 0` for a conformal map.
///
/// Each coefficient is `(re, im)`, already scaled by `1 / sqrt(2 · 2A)` so that
/// the squared row residual is the area-weighted energy of the triangle.  `None`
/// for a degenerate triangle, which constrains nothing.
pub(super) fn conformal_row(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> Option<[(f64, f64); 3]> {
    let e1 = sub3(p1, p0);
    let d2 = sub3(p2, p0);
    let x1 = dot3(e1, e1).sqrt();
    let doubled_area = dot3(cross3(e1, d2), cross3(e1, d2)).sqrt();
    let side_product = x1 * dot3(d2, d2).sqrt();
    let well_formed = x1 > 0.0 && doubled_area > DEGENERATE_SINE * side_product;
    if !well_formed {
        return None;
    }
    let x2 = dot3(d2, e1) / x1;
    let y2 = doubled_area / x1;
    let weight = 1.0 / (2.0 * doubled_area).sqrt();
    Some([
        (weight * (x2 - x1), weight * y2),
        (weight * -x2, weight * -y2),
        (weight * x1, 0.0),
    ])
}

/// The normal equations `M x = b` of the LSCM least-squares problem over the free
/// (not pinned, referenced) vertices; the unknowns are interleaved `(u, v)` pairs.
pub(super) struct NormalEquations {
    /// Free vertices in ascending id order.
    pub(super) free_vertices: Vec<usize>,
    /// Symmetric positive semi-definite matrix of size `2 * free_vertices.len()`.
    pub(super) matrix: CsrMatrix,
    /// Right-hand side: minus the coupling of the free vertices to the pinned ones.
    pub(super) rhs: Vec<f64>,
}

/// Assemble the normal equations of LSCM for the given pins.
///
/// The Hermitian form `f^H conj(c) c^T f` of every triangle row is added block by
/// block: the complex entry `α + iβ` of vertex pair `(j, k)` becomes the real
/// `2 x 2` block `[[α, −β], [β, α]]`.  Entries that couple to a pinned vertex move
/// to the right-hand side.
///
/// # Errors
/// * [`Error::IndexOutOfBounds`] for a triangle or pin outside `positions`;
/// * [`Error::InvalidMesh`] for equal pins, a pin that no triangle uses, or a mesh
///   that is empty or has several connected components;
/// * [`Error::DegenerateGeometry`] for non-finite positions, a mesh whose triangles
///   are all degenerate, non-degenerate triangles that do not form a single patch
///   with both pins in it, or a free vertex that only degenerate triangles use.
pub(super) fn assemble_normal_equations(
    positions: &[[f64; 3]],
    triangles: &[[usize; 3]],
    pins: [(usize, [f64; 2]); 2],
) -> Result<NormalEquations> {
    let n = positions.len();
    let connectivity = analyse_connectivity(n, triangles)?;
    require_single_component(&connectivity)?;
    for &(pin, _) in &pins {
        if pin >= n {
            return Err(Error::index_out_of_bounds(pin, n));
        }
        if !connectivity.referenced[pin] {
            return Err(Error::invalid_mesh(format!(
                "pinned vertex {pin} is not used by any triangle"
            )));
        }
    }
    if pins[0].0 == pins[1].0 {
        return Err(Error::invalid_mesh(
            "LSCM needs two distinct pinned vertices",
        ));
    }
    let finite = |p: &[f64; 3]| p.iter().all(|c| c.is_finite());
    if (0..n).any(|v| connectivity.referenced[v] && !finite(&positions[v])) {
        return Err(Error::degenerate_geometry(
            "a vertex position is not finite",
        ));
    }

    let mut free_of = vec![usize::MAX; n];
    let mut free_vertices: Vec<usize> = Vec::new();
    for (v, slot) in free_of.iter_mut().enumerate() {
        if connectivity.referenced[v] && v != pins[0].0 && v != pins[1].0 {
            *slot = free_vertices.len();
            free_vertices.push(v);
        }
    }
    let free_count = free_vertices.len();
    let pinned_uv = |v: usize| -> [f64; 2] { if v == pins[0].0 { pins[0].1 } else { pins[1].1 } };

    // Conformality rows of the non-degenerate triangles.
    let rows: Vec<Option<[(f64, f64); 3]>> = triangles
        .iter()
        .map(|t| conformal_row(positions[t[0]], positions[t[1]], positions[t[2]]))
        .collect();
    if rows.iter().all(Option::is_none) {
        return Err(Error::degenerate_geometry("every triangle is degenerate"));
    }
    // The non-degenerate triangles alone must hold the patch together, pins
    // included: a part that is joined to the rest only through degenerate
    // triangles constrains nothing, so nothing would fix its position, rotation
    // and scale.
    let regular: Vec<[usize; 3]> = triangles
        .iter()
        .zip(rows.iter())
        .filter(|(_, row)| row.is_some())
        .map(|(tri, _)| *tri)
        .collect();
    let regular_connectivity = analyse_connectivity(n, &regular)?;
    if regular_connectivity.components != 1
        || pins
            .iter()
            .any(|&(pin, _)| !regular_connectivity.referenced[pin])
    {
        return Err(Error::degenerate_geometry(
            "the non-degenerate triangles do not form one patch that contains both pinned vertices",
        ));
    }

    // Block pattern: free vertices that share a non-degenerate triangle, including
    // each vertex with itself.
    let mut pattern: Vec<Vec<usize>> = vec![Vec::new(); free_count];
    for (tri, row) in triangles.iter().zip(rows.iter()) {
        if row.is_none() {
            continue;
        }
        for &j in tri {
            for &k in tri {
                if free_of[j] != usize::MAX && free_of[k] != usize::MAX {
                    pattern[free_of[j]].push(free_of[k]);
                }
            }
        }
    }
    for (f, columns) in pattern.iter_mut().enumerate() {
        columns.push(f);
        columns.sort_unstable();
        columns.dedup();
    }

    // CSR layout: row 2f and row 2f + 1 both hold the columns (2k, 2k + 1) of every
    // block column k of the pattern.
    let mut row_ptr = Vec::with_capacity(2 * free_count + 1);
    let mut col_idx: Vec<usize> = Vec::new();
    row_ptr.push(0);
    for columns in &pattern {
        for _ in 0..2 {
            for &k in columns {
                col_idx.push(2 * k);
                col_idx.push(2 * k + 1);
            }
            row_ptr.push(col_idx.len());
        }
    }
    let mut values = vec![0.0; col_idx.len()];
    let mut rhs = vec![0.0; 2 * free_count];
    let mut touched = vec![false; free_count];

    for (tri, row) in triangles.iter().zip(rows.iter()) {
        let Some(c) = row else { continue };
        for j in 0..3 {
            let fj = free_of[tri[j]];
            if fj == usize::MAX {
                continue;
            }
            touched[fj] = true;
            for k in 0..3 {
                // conj(c_j) * c_k = alpha + i beta.
                let alpha = c[j].0 * c[k].0 + c[j].1 * c[k].1;
                let beta = c[j].0 * c[k].1 - c[j].1 * c[k].0;
                let fk = free_of[tri[k]];
                if fk == usize::MAX {
                    let [uk, vk] = pinned_uv(tri[k]);
                    rhs[2 * fj] -= alpha * uk - beta * vk;
                    rhs[2 * fj + 1] -= beta * uk + alpha * vk;
                } else if let Ok(slot) = pattern[fj].binary_search(&fk) {
                    let upper = row_ptr[2 * fj] + 2 * slot;
                    let lower = row_ptr[2 * fj + 1] + 2 * slot;
                    values[upper] += alpha;
                    values[upper + 1] -= beta;
                    values[lower] += beta;
                    values[lower + 1] += alpha;
                }
            }
        }
    }
    if touched.iter().any(|&t| !t) {
        return Err(Error::degenerate_geometry(
            "a vertex is used only by degenerate triangles",
        ));
    }
    Ok(NormalEquations {
        free_vertices,
        matrix: CsrMatrix {
            nrows: 2 * free_count,
            ncols: 2 * free_count,
            row_ptr,
            col_idx,
            values,
        },
        rhs,
    })
}

/// Least-squares conformal map of a triangle mesh with two pinned vertices.
///
/// `pins` gives two distinct vertices and their UV positions; every other vertex
/// is placed to minimise the conformal energy.  A developable patch is unfolded
/// exactly (zero energy) up to the similarity that the pins impose, so with pins at
/// their true 3-D distance a planar patch is reproduced isometrically.  Vertices
/// that no triangle uses receive `[0, 0]`.
///
/// # Errors
/// * [`Error::IndexOutOfBounds`] for a triangle or pin outside `positions`;
/// * [`Error::InvalidMesh`] for equal pins, a pin that no triangle uses, or a mesh
///   that is empty or has several connected components;
/// * [`Error::DegenerateGeometry`] for non-finite positions, a mesh whose triangles
///   are all degenerate, non-degenerate triangles that do not form a single patch
///   with both pins in it, or a vertex that only degenerate triangles use;
/// * [`Error::ConvergenceFailure`] if the conjugate-gradient solver stalls.
pub fn lscm_solve(
    positions: &[[f64; 3]],
    triangles: &[[usize; 3]],
    pins: [(usize, [f64; 2]); 2],
) -> Result<Vec<[f64; 2]>> {
    let system = assemble_normal_equations(positions, triangles, pins)?;
    let mut uvs = vec![[0.0f64; 2]; positions.len()];
    for &(pin, uv) in &pins {
        uvs[pin] = uv;
    }
    if system.free_vertices.is_empty() {
        return Ok(uvs);
    }
    let solution = solve_spd(&system.matrix, &system.rhs, "least squares conformal map")?;
    for (slot, &v) in system.free_vertices.iter().enumerate() {
        uvs[v] = [solution[2 * slot], solution[2 * slot + 1]];
    }
    Ok(uvs)
}

/// Largest number of pin candidates for which every pair is compared.
///
/// `4096` candidates are about 8.4 million distance evaluations per pass, a few
/// tens of milliseconds; larger boundaries fall back to a double sweep.
const EXACT_PAIR_LIMIT: usize = 4096;

/// Relative slack within which two squared distances count as equal.
///
/// Geometrically tied candidates (the two diagonals of a square patch, the
/// vertices of a regular ring) differ by rounding noise of about `1e-16` once the
/// mesh is moved or rotated; if that noise decided the pins, the UV layout would
/// depend on the pose of the mesh.  Among candidates whose distance is within this
/// slack of the largest one, the smallest vertex ids win.
const TIE_TOLERANCE: f64 = 1.0e-9;

/// Smallest threshold a squared distance must reach to tie with `largest`.
fn tie_threshold(largest: f64) -> f64 {
    largest * (1.0 - TIE_TOLERANCE)
}

/// Pick the two pinned vertices among `candidates`: the farthest-apart pair.
///
/// With at most [`EXACT_PAIR_LIMIT`] distinct candidates every pair is compared,
/// so the pair realises the exact diameter of the candidate set.  Larger sets use
/// a double sweep (the candidate farthest from the smallest one, then the
/// candidate farthest from that one), which finds a pair at least half as far
/// apart as the diameter and, for the elongated boundaries of typical patches,
/// the diameter itself.  Ties, up to rounding noise ([`TIE_TOLERANCE`]), go to the
/// smaller vertex ids.  The vertex with the smaller id is pinned to `(0, 0)` and
/// the other to `(d, 0)`, where `d` is their 3-D distance, so that the UV layout
/// has the scale of the mesh.
///
/// # Errors
/// [`Error::TooFewPoints`] for fewer than two distinct candidates and
/// [`Error::IndexOutOfBounds`] for a candidate outside `positions`.
pub(super) fn farthest_pair_pins(
    positions: &[[f64; 3]],
    mut candidates: Vec<usize>,
) -> Result<[(usize, [f64; 2]); 2]> {
    candidates.sort_unstable();
    candidates.dedup();
    if let Some(&last) = candidates.last()
        && last >= positions.len()
    {
        return Err(Error::index_out_of_bounds(last, positions.len()));
    }
    let (Some(&smallest), Some(&second_smallest)) = (candidates.first(), candidates.get(1)) else {
        return Err(Error::too_few_points(2, candidates.len()));
    };
    let squared_distance = |a: usize, b: usize| {
        let d = sub3(positions[a], positions[b]);
        dot3(d, d)
    };
    let fallback = (smallest, second_smallest);
    let (first, second) = if candidates.len() <= EXACT_PAIR_LIMIT {
        // Pass 1: the largest squared distance.  Pass 2: the first pair, in
        // ascending (a, b) order, that ties with it.
        let mut largest = f64::NEG_INFINITY;
        for (at, &a) in candidates.iter().enumerate() {
            for &b in &candidates[at + 1..] {
                largest = largest.max(squared_distance(a, b));
            }
        }
        let threshold = tie_threshold(largest);
        candidates
            .iter()
            .enumerate()
            .find_map(|(at, &a)| {
                candidates[at + 1..]
                    .iter()
                    .find(|&&b| squared_distance(a, b) >= threshold)
                    .map(|&b| (a, b))
            })
            .unwrap_or(fallback)
    } else {
        let farthest_from = |from: usize| -> usize {
            let largest = candidates
                .iter()
                .map(|&v| squared_distance(from, v))
                .fold(f64::NEG_INFINITY, f64::max);
            let threshold = tie_threshold(largest);
            candidates
                .iter()
                .copied()
                .find(|&v| squared_distance(from, v) >= threshold)
                .unwrap_or(smallest)
        };
        let b = farthest_from(smallest);
        let c = farthest_from(b);
        if b == c {
            // All candidates coincide: any two distinct ones will do.
            fallback
        } else {
            (b.min(c), b.max(c))
        }
    };
    let separation = squared_distance(first, second).sqrt();
    let span = if separation.is_finite() && separation > 0.0 {
        separation
    } else {
        1.0
    };
    Ok([(first, [0.0, 0.0]), (second, [span, 0.0])])
}

/// LSCM parameterization of a disc-like mesh with automatically chosen pins: the
/// farthest-apart pair of boundary vertices (see [`farthest_pair_pins`]).
///
/// # Errors
/// [`Error::InvalidMesh`] for a mesh without triangles, with several connected
/// components or without a boundary (a closed surface cannot be unfolded), plus
/// the errors of [`lscm_solve`].
pub(super) fn parameterize(mesh: &ParamTriMesh) -> Result<Vec<[f64; 2]>> {
    let n = mesh.vertex_count();
    if n == 0 {
        return Ok(Vec::new());
    }
    let connectivity = analyse_connectivity(n, &mesh.triangles)?;
    require_single_component(&connectivity)?;
    let boundary: Vec<usize> = boundary_loops(&mesh.triangles)
        .into_iter()
        .flatten()
        .collect();
    if boundary.is_empty() {
        return Err(Error::invalid_mesh(
            "a closed mesh cannot be unfolded by LSCM; cut a seam first",
        ));
    }
    let pins = farthest_pair_pins(&mesh.positions, boundary)?;
    lscm_solve(&mesh.positions, &mesh.triangles, pins)
}
