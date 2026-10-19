// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Voronoi-pole analysis of an unorganised boundary sample (Amenta, Choi &
//! Kolluri, "The Power Crust", 2001; Amenta & Bern, "Surface reconstruction by
//! Voronoi filtering", 1999).
//!
//! The Voronoi vertices of a dense surface sample are the circumcentres of its
//! Delaunay tetrahedra. For every sample `s` the *poles* are the farthest
//! Voronoi vertices of its Voronoi cell on either side of the surface; for an
//! ε-sample they converge to the medial axis, and the distance from a pole to
//! `s` converges to the radius of the corresponding medial ball.
//!
//! This module builds the exact Delaunay tetrahedralisation, derives an
//! *oriented* outward normal for every sample, selects the inner pole of every
//! sample and classifies every Voronoi vertex (Delaunay tet) as lying inside or
//! outside the sampled solid:
//!
//! * samples on the convex hull get their outward direction from the
//!   area-weighted normals of their incident hull faces (Amenta's convention of
//!   a pole "at infinity");
//! * every other sample takes the direction towards its farthest Voronoi vertex
//!   as an unoriented normal axis, and is oriented by a minimum-spanning-tree
//!   propagation (Hoppe et al. 1992) seeded at the hull samples over the
//!   Delaunay graph, preferring short, normal-consistent edges;
//! * the inner pole of a sample `s` with outward normal `n` is the farthest
//!   circumcentre `c` of a tet incident to `s` with `(c − s)·n < 0`;
//! * the normal is finally refined to the pole vector `s − c`, which is a far
//!   better normal estimate than hull faces or the propagation axis.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use super::{dist3, dot3, sub3};
use crate::computational_geometry::delaunay_3d::{Tetrahedralization, delaunay_3d};
use crate::error::{Error, Result};

/// Circumsphere of one Delaunay tetrahedron (= one Voronoi vertex).
#[derive(Debug, Clone, Copy)]
pub(super) struct Circumsphere {
    /// Circumcentre (the Voronoi vertex).
    pub(super) center: [f64; 3],
    /// Circumradius (distance from the Voronoi vertex to its nearest samples).
    pub(super) radius: f64,
}

/// Result of the Voronoi-pole analysis of a boundary sample.
#[derive(Debug, Clone)]
pub(super) struct PoleAnalysis {
    /// Delaunay tetrahedralisation of the samples.
    pub(super) tetra: Tetrahedralization,
    /// Finite, well-conditioned circumsphere of every tet (`None` for flat or
    /// unboundedly large tets, which cannot carry an interior medial ball).
    pub(super) spheres: Vec<Option<Circumsphere>>,
    /// Whether the Voronoi vertex of each tet lies inside the sampled solid.
    pub(super) inside: Vec<bool>,
    /// Oriented outward unit normal of every sample (`None` for samples that
    /// are not referenced by any tet, i.e. exact duplicates).
    pub(super) normals: Vec<Option<[f64; 3]>>,
    /// Distance from every sample to its nearest Delaunay neighbour (the local
    /// sampling density); `f64::INFINITY` for unreferenced samples.
    pub(super) spacing: Vec<f64>,
    /// Delaunay-graph neighbours of every sample (sorted, deduplicated).
    pub(super) neighbors: Vec<Vec<usize>>,
    /// Inner pole of every sample as a tet index, after discarding poles whose
    /// radius is below `min_radius_factor × spacing`.
    pub(super) inner_pole: Vec<Option<usize>>,
    /// Length of the diagonal of the samples' bounding box.
    pub(super) diagonal: f64,
}

/// Total-order wrapper so `f64` keys can live in a [`BinaryHeap`] (min-heap).
#[derive(Debug, Clone, Copy, PartialEq)]
struct MinKey {
    key: f64,
    from: usize,
    to: usize,
}

impl Eq for MinKey {}

impl PartialOrd for MinKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MinKey {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse for a min-heap; tie-break on indices for determinism.
        other
            .key
            .total_cmp(&self.key)
            .then_with(|| other.to.cmp(&self.to))
            .then_with(|| other.from.cmp(&self.from))
    }
}

/// Run the full Voronoi-pole analysis on `boundary`.
///
/// `min_radius_factor` (λ) drops the inner pole of a sample when its radius is
/// below `λ × spacing`, i.e. when the pole is not resolved by the sampling
/// (sharp creases, noise).
///
/// # Errors
///
/// * [`Error::TooFewPoints`] when fewer than four samples are given;
/// * [`Error::DegenerateGeometry`] for non-finite coordinates or when all
///   samples are coplanar (no 3-D Delaunay complex exists).
pub(super) fn analyze(boundary: &[[f64; 3]], min_radius_factor: f64) -> Result<PoleAnalysis> {
    let n = boundary.len();
    if n < 4 {
        return Err(Error::too_few_points(4, n));
    }
    if boundary.iter().any(|p| p.iter().any(|c| !c.is_finite())) {
        return Err(Error::degenerate_geometry(
            "medial axis: boundary sample has a non-finite coordinate",
        ));
    }
    let tetra = delaunay_3d(boundary).ok_or_else(|| {
        Error::degenerate_geometry(
            "medial axis: boundary samples are coplanar, so no 3-D Delaunay complex exists",
        )
    })?;
    let diagonal = bbox_diagonal(boundary);

    let spheres: Vec<Option<Circumsphere>> = tetra
        .tets
        .iter()
        .map(|t| {
            circumsphere(
                boundary[t.v[0]],
                boundary[t.v[1]],
                boundary[t.v[2]],
                boundary[t.v[3]],
            )
            .filter(|s| s.radius <= diagonal)
        })
        .collect();

    // Vertex → incident tets, and the Delaunay graph.
    let mut vertex_tets: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (ti, t) in tetra.tets.iter().enumerate() {
        for a in 0..4 {
            vertex_tets[t.v[a]].push(ti);
            for b in 0..4 {
                if a != b {
                    neighbors[t.v[a]].push(t.v[b]);
                }
            }
        }
    }
    for list in &mut neighbors {
        list.sort_unstable();
        list.dedup();
    }

    let spacing: Vec<f64> = (0..n)
        .map(|i| {
            neighbors[i]
                .iter()
                .map(|&j| dist3(boundary[i], boundary[j]))
                .fold(f64::INFINITY, f64::min)
        })
        .collect();

    let hull_normals = hull_vertex_normals(&tetra, boundary);
    let axes = unoriented_axes(boundary, &neighbors, &vertex_tets, &spheres, &hull_normals);
    let coarse_normals = orient_normals(boundary, &neighbors, &spacing, &axes, &hull_normals);

    // Raw inner poles (before the λ filter) drive the inside classification.
    let raw_pole: Vec<Option<usize>> = (0..n)
        .map(|i| {
            let ni = coarse_normals[i]?;
            let mut best: Option<(usize, f64)> = None;
            for &ti in &vertex_tets[i] {
                if let Some(s) = spheres[ti]
                    && dot3(sub3(s.center, boundary[i]), ni) < 0.0
                    && best.is_none_or(|(_, r)| s.radius > r)
                {
                    best = Some((ti, s.radius));
                }
            }
            best.map(|(ti, _)| ti)
        })
        .collect();

    // Refine every normal with its pole vector (sample minus inner pole),
    // which converges to the true surface normal for dense samples
    // (Amenta & Bern 1999); it keeps the propagated orientation because the
    // inner pole lies on the inner side by construction.
    let normals: Vec<Option<[f64; 3]>> = (0..n)
        .map(|i| {
            let coarse = coarse_normals[i]?;
            Some(
                raw_pole[i]
                    .and_then(|ti| spheres[ti])
                    .and_then(|s| unit_or_none(sub3(boundary[i], s.center)))
                    .unwrap_or(coarse),
            )
        })
        .collect();

    let mut inside = vec![false; tetra.tets.len()];
    for (ti, t) in tetra.tets.iter().enumerate() {
        let Some(s) = spheres[ti] else { continue };
        if s.radius <= 0.0 {
            continue;
        }
        let mut vote = 0.0;
        for &v in &t.v {
            if let Some(nv) = normals[v] {
                vote += dot3(sub3(s.center, boundary[v]), nv) / s.radius;
            }
        }
        inside[ti] = vote < 0.0;
    }
    for &ti in raw_pole.iter().flatten() {
        inside[ti] = true;
    }

    let inner_pole: Vec<Option<usize>> = raw_pole
        .iter()
        .enumerate()
        .map(|(i, pole)| {
            pole.filter(|&ti| {
                spheres[ti].is_some_and(|s| s.radius >= min_radius_factor * spacing[i])
            })
        })
        .collect();

    Ok(PoleAnalysis {
        tetra,
        spheres,
        inside,
        normals,
        spacing,
        neighbors,
        inner_pole,
        diagonal,
    })
}

/// Bounding-box diagonal of a point set (0 for an empty set).
pub(super) fn bbox_diagonal(points: &[[f64; 3]]) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    if points.is_empty() {
        return 0.0;
    }
    dist3(lo, hi)
}

/// Circumsphere of the tetrahedron `(a, b, c, d)`.
///
/// Solves `2 (p_k − a) · x = |p_k − a|²` for `k ∈ {b, c, d}` (with `x` the
/// centre relative to `a`) by Cramer's rule. Returns `None` when the tet is
/// numerically flat: the determinant is compared against the product of the
/// three edge lengths, so the test is scale invariant.
pub(super) fn circumsphere(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
) -> Option<Circumsphere> {
    let r1 = sub3(b, a);
    let r2 = sub3(c, a);
    let r3 = sub3(d, a);
    let det = r1[0] * (r2[1] * r3[2] - r2[2] * r3[1]) - r1[1] * (r2[0] * r3[2] - r2[2] * r3[0])
        + r1[2] * (r2[0] * r3[1] - r2[1] * r3[0]);
    let scale = dot3(r1, r1).sqrt() * dot3(r2, r2).sqrt() * dot3(r3, r3).sqrt();
    if !det.is_finite() || scale <= 0.0 || det.abs() <= 1e-10 * scale {
        return None;
    }
    let h1 = 0.5 * dot3(r1, r1);
    let h2 = 0.5 * dot3(r2, r2);
    let h3 = 0.5 * dot3(r3, r3);
    // Cramer's rule with the rhs substituted column-wise.
    let x = (h1 * (r2[1] * r3[2] - r2[2] * r3[1]) - r1[1] * (h2 * r3[2] - r2[2] * h3)
        + r1[2] * (h2 * r3[1] - r2[1] * h3))
        / det;
    let y = (r1[0] * (h2 * r3[2] - r2[2] * h3) - h1 * (r2[0] * r3[2] - r2[2] * r3[0])
        + r1[2] * (r2[0] * h3 - h2 * r3[0]))
        / det;
    let z = (r1[0] * (r2[1] * h3 - h2 * r3[1]) - r1[1] * (r2[0] * h3 - h2 * r3[0])
        + h1 * (r2[0] * r3[1] - r2[1] * r3[0]))
        / det;
    let center = [a[0] + x, a[1] + y, a[2] + z];
    let radius = (x * x + y * y + z * z).sqrt();
    if !radius.is_finite() {
        return None;
    }
    Some(Circumsphere { center, radius })
}

/// Area-weighted outward normal of every convex-hull vertex (`None` for
/// vertices not on the hull).
///
/// A hull face is a tet face without a neighbour (`neigh == -1`). Its outward
/// sense is fixed against the tet's opposite vertex rather than trusting the
/// face winding.
fn hull_vertex_normals(tetra: &Tetrahedralization, points: &[[f64; 3]]) -> Vec<Option<[f64; 3]>> {
    let mut acc: Vec<Option<[f64; 3]>> = vec![None; points.len()];
    for t in &tetra.tets {
        for i in 0..4 {
            if t.neigh[i] >= 0 {
                continue;
            }
            let f = t.face(i);
            let a = points[f[0]];
            let e1 = sub3(points[f[1]], a);
            let e2 = sub3(points[f[2]], a);
            let mut nrm = cross3(e1, e2);
            if dot3(nrm, sub3(points[t.v[i]], a)) > 0.0 {
                nrm = [-nrm[0], -nrm[1], -nrm[2]];
            }
            for &v in &f {
                let cur = acc[v].unwrap_or([0.0; 3]);
                acc[v] = Some([cur[0] + nrm[0], cur[1] + nrm[1], cur[2] + nrm[2]]);
            }
        }
    }
    acc.into_iter().map(|o| o.and_then(unit_or_none)).collect()
}

/// Unoriented normal axis of every sample.
///
/// Hull samples use their (already oriented) hull normal. Interior samples use
/// the direction towards their farthest finite Voronoi vertex (Amenta's first
/// pole). A sample without any well-conditioned incident tet falls back to its
/// umbrella vector (sample minus the centroid of its Delaunay neighbours).
fn unoriented_axes(
    points: &[[f64; 3]],
    neighbors: &[Vec<usize>],
    vertex_tets: &[Vec<usize>],
    spheres: &[Option<Circumsphere>],
    hull_normals: &[Option<[f64; 3]>],
) -> Vec<Option<[f64; 3]>> {
    (0..points.len())
        .map(|i| {
            if vertex_tets[i].is_empty() {
                return None;
            }
            if let Some(h) = hull_normals[i] {
                return Some(h);
            }
            let mut best: Option<Circumsphere> = None;
            for &ti in &vertex_tets[i] {
                if let Some(s) = spheres[ti]
                    && best.is_none_or(|b| s.radius > b.radius)
                {
                    best = Some(s);
                }
            }
            if let Some(b) = best
                && let Some(u) = unit_or_none(sub3(b.center, points[i]))
            {
                return Some(u);
            }
            // Umbrella fallback: sample minus the centroid of its neighbours.
            let nb = &neighbors[i];
            if nb.is_empty() {
                return None;
            }
            let inv = 1.0 / nb.len() as f64;
            let mut centroid = [0.0; 3];
            for &j in nb {
                for k in 0..3 {
                    centroid[k] += points[j][k] * inv;
                }
            }
            unit_or_none(sub3(points[i], centroid))
        })
        .collect()
}

/// Orient the unoriented normal axes consistently by a minimum-spanning-tree
/// propagation over the Delaunay graph, seeded at the hull samples.
///
/// The edge key `|p_i − p_j| / min(h_i, h_j) + (1 − |a_i · a_j|)` (with `h`
/// the local spacing and `a` the axes) makes the propagation travel along
/// short surface-neighbour edges between nearly parallel normals first, so it
/// does not jump across thin parts of the solid.
fn orient_normals(
    points: &[[f64; 3]],
    neighbors: &[Vec<usize>],
    spacing: &[f64],
    axes: &[Option<[f64; 3]>],
    hull_normals: &[Option<[f64; 3]>],
) -> Vec<Option<[f64; 3]>> {
    let n = points.len();
    let mut oriented: Vec<Option<[f64; 3]>> = vec![None; n];
    let mut heap = BinaryHeap::new();
    let push_edges = |heap: &mut BinaryHeap<MinKey>, i: usize, ni: [f64; 3]| {
        for &j in &neighbors[i] {
            let Some(aj) = axes[j] else { continue };
            let h = spacing[i].min(spacing[j]).max(f64::MIN_POSITIVE);
            let key = dist3(points[i], points[j]) / h + (1.0 - dot3(ni, aj).abs());
            heap.push(MinKey {
                key,
                from: i,
                to: j,
            });
        }
    };
    for ((slot, &hull), axis) in oriented.iter_mut().zip(hull_normals).zip(axes) {
        if let Some(h) = hull
            && axis.is_some()
        {
            *slot = Some(h);
        }
    }
    for (i, &slot) in oriented.iter().enumerate() {
        if let Some(ni) = slot {
            push_edges(&mut heap, i, ni);
        }
    }
    while let Some(MinKey { from, to, .. }) = heap.pop() {
        if oriented[to].is_some() {
            continue;
        }
        let (Some(nf), Some(at)) = (oriented[from], axes[to]) else {
            continue;
        };
        let nt = if dot3(nf, at) < 0.0 {
            [-at[0], -at[1], -at[2]]
        } else {
            at
        };
        oriented[to] = Some(nt);
        push_edges(&mut heap, to, nt);
    }
    oriented
}

/// Cross product.
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Normalise `v`, or `None` when it has (numerically) zero length.
pub(super) fn unit_or_none(v: [f64; 3]) -> Option<[f64; 3]> {
    let len = dot3(v, v).sqrt();
    if len.is_finite() && len > 1e-300 {
        Some([v[0] / len, v[1] / len, v[2] / len])
    } else {
        None
    }
}
