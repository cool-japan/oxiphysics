// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Curve-skeleton extraction by Laplacian contraction of a point cloud
//! (Au, Tai, Chu, Cohen-Or & Lee, "Skeleton extraction by mesh contraction",
//! SIGGRAPH 2008; point-cloud variant of Cao, Tagliasacchi, Olson, Zhang &
//! Su, "Point cloud skeletons via Laplacian-based contraction", SMI 2010).
//!
//! **Neighbourhood.** Each sample gets its `k` nearest neighbours; they are
//! projected onto the PCA tangent plane and the sample's *one-ring* is the set
//! of neighbours whose Voronoi edge borders the sample's 2-D Voronoi cell
//! (i.e. its Delaunay neighbours in the tangent plane), in angular order.
//! Consecutive ring members span the one-ring triangles.
//!
//! **Contraction.** Each iteration solves, per coordinate, the least-squares
//! system
//!
//! ```text
//! [ W_L · L ]          [   0     ]
//! [  W_H    ] · P'  =  [ W_H · P ]
//! ```
//!
//! through its normal equations `(W_L² LᵀL + W_H²) P' = W_H² P` with a
//! Jacobi-preconditioned conjugate-gradient solver. `L` is the cotangent
//! Laplacian of the one-ring triangles of the *input* cloud (negative
//! cotangents clamped to zero; a sample without a non-degenerate triangle
//! falls back to uniform weights over its ring). The weights are kept fixed:
//! re-evaluating them on the contracting cloud would feed the exploding
//! cotangents of collapsing (degenerate) rings into the solve.
//! Afterwards `W_L ← s_L · W_L` (`s_L = 3`, capped at 2048) and
//! `W_H,i ← W_H,i⁰ · √(A_i⁰ / A_i)` where `A_i` is the one-ring area, so
//! already-collapsed regions are held in place. A sample whose one-ring area
//! has dropped below `10⁻³` of its initial value has collapsed onto the
//! skeleton: its Laplacian row is dropped and its attraction weight is raised
//! to the maximum, pinning it (otherwise continued smoothing would shrink the
//! skeleton along its own length and clump it). The iteration stops once the
//! total one-ring area is below `10⁻²` of its initial value, i.e. when the
//! cross-sections have collapsed (cf. the termination test of Cao et al.), or
//! when, below `10⁻¹`, an iteration no longer halves the area (the collapse
//! has stalled; that iterate only shrank the skeleton and is discarded).
//!
//! **Extraction.** The contracted cloud is sampled by farthest-point sampling
//! with a radius of twice the mean sample spacing; every sample is assigned to
//! its nearest node, nodes whose samples are neighbours are connected, each
//! node is moved to the centre of the sphere fitted (algebraic least squares)
//! to the *original* positions of its samples — a local medial-ball estimate
//! that stays centred even when a cluster covers only part of a
//! cross-section, falling back to their centroid (Au's embedding refinement)
//! for flat or ill-conditioned clusters — with the mean distance to them as its
//! radius,
//! and the node graph is skeletonised (spanning tree, twig pruning, loop
//! closure).

use super::kdtree::KdTree;
use super::poles::{bbox_diagonal, cross3, unit_or_none};
use super::skeleton_graph::{SkelNode, SkeletonParams, farthest_point_clusters, skeletonize};
use super::{CurveSkeleton, CurveSkeletonNode, dist3, dot3, sub3};
use crate::error::{Error, Result};

/// Number of nearest neighbours used to build the one-rings.
const KNN: usize = 12;
/// Laplacian-weight growth factor per iteration.
const S_L: f64 = 3.0;
/// Upper bound of the Laplacian weight.
const W_L_MAX: f64 = 2048.0;
/// Upper bound of the attraction weight.
const W_H_MAX: f64 = 1e6;
/// Upper bound of a single cotangent.
const COT_MAX: f64 = 1e4;
/// A sample whose one-ring area has shrunk below this fraction of its initial
/// area has collapsed onto the skeleton and is frozen.
const FREEZE_RATIO: f64 = 1e-3;
/// Contraction terminates once the total one-ring area has shrunk below this
/// fraction of its initial value: the cross-sections have collapsed and further
/// iterations would only shrink the skeleton itself along its slow modes.
const TERMINATION_RATIO: f64 = 1e-2;
/// Once the total one-ring area is below this fraction of its initial value,
/// an iteration that fails to halve it means the cross-sections have stopped
/// collapsing; continuing would only shrink the skeleton along its own slow
/// modes (e.g. close the ring of a torus), so the contraction stops.
const STALL_ACTIVE_RATIO: f64 = 0.1;

/// One-ring neighbourhood of a sample.
#[derive(Debug, Clone)]
struct Ring {
    /// Ring members in angular order around the sample.
    members: Vec<usize>,
    /// `closed[k]`: members `k` and `k + 1 (mod m)` span a one-ring triangle.
    closed: Vec<bool>,
}

/// Validate a boundary sample for contraction.
fn validate(points: &[[f64; 3]], steps: usize) -> Result<()> {
    if points.len() < 4 {
        return Err(Error::too_few_points(4, points.len()));
    }
    if points.iter().any(|p| p.iter().any(|c| !c.is_finite())) {
        return Err(Error::degenerate_geometry(
            "laplacian contraction: point has a non-finite coordinate",
        ));
    }
    if steps == 0 {
        return Err(Error::out_of_range("contraction_steps", 0.0, 1.0, f64::MAX));
    }
    if bbox_diagonal(points) <= 0.0 {
        return Err(Error::degenerate_geometry(
            "laplacian contraction: all points coincide",
        ));
    }
    Ok(())
}

/// Contract `points` with `steps` iterations of Laplacian contraction and
/// return the contracted positions (same order as the input).
pub(super) fn contract(points: &[[f64; 3]], steps: usize) -> Result<Vec<[f64; 3]>> {
    validate(points, steps)?;
    let n = points.len();
    let all: Vec<usize> = (0..n).collect();
    let tree = KdTree::build(points, &all);
    let rings: Vec<Ring> = (0..n).map(|i| one_ring(points, &tree, i)).collect();

    let area0: Vec<f64> = (0..n).map(|i| ring_area(points, i, &rings[i])).collect();
    let total0: f64 = area0.iter().sum();
    let mut w_l = 1.0_f64;
    let mut w_h = vec![1.0_f64; n];
    let mut frozen = vec![false; n];
    let mut current = points.to_vec();
    let mut previous_total = total0;
    for _ in 0..steps {
        let lap = cotangent_laplacian(points, &rings, &frozen);
        let next = solve_contraction(&lap, w_l, &w_h, &current)?;
        let before = std::mem::replace(&mut current, next);
        let area: Vec<f64> = (0..n).map(|i| ring_area(&current, i, &rings[i])).collect();
        let total: f64 = area.iter().sum();
        w_l = (w_l * S_L).min(W_L_MAX);
        for i in 0..n {
            w_h[i] = if area[i] > 0.0 && area0[i] > 0.0 {
                (area0[i] / area[i]).sqrt().min(W_H_MAX)
            } else {
                W_H_MAX
            };
            // Once a one-ring has collapsed the sample lies on the skeleton;
            // further smoothing would only shorten or clump the skeleton, so
            // the sample is pinned where it is.
            if area0[i] > 0.0 && area[i] < FREEZE_RATIO * area0[i] {
                frozen[i] = true;
            }
            if frozen[i] {
                w_h[i] = W_H_MAX;
            }
        }
        let collapsed = total0 > 0.0 && total / total0 < TERMINATION_RATIO;
        let stalled =
            total0 > 0.0 && total / total0 < STALL_ACTIVE_RATIO && total > 0.5 * previous_total;
        if stalled {
            // The stalled iterate only shrank the skeleton: discard it.
            current = before;
            break;
        }
        if collapsed || frozen.iter().all(|&f| f) {
            break;
        }
        previous_total = total;
    }
    Ok(current)
}

/// Full curve-skeleton extraction (contraction + node extraction).
pub(super) fn curve_skeleton(points: &[[f64; 3]], steps: usize) -> Result<CurveSkeleton> {
    let contracted = contract(points, steps)?;
    let n = points.len();
    let all: Vec<usize> = (0..n).collect();
    let tree = KdTree::build(points, &all);
    let spacing: f64 = (0..n)
        .map(|i| {
            tree.nearest_filtered(points[i], |j| j != i)
                .map_or(0.0, |(_, d2)| d2.sqrt())
        })
        .sum::<f64>()
        / n as f64;
    let diagonal = bbox_diagonal(points);
    let sample_radius = (2.0 * spacing).max(1e-9 * diagonal);

    // Farthest-point sampling of the contracted cloud.
    let (assign, m) = farthest_point_clusters(&contracted, sample_radius);

    // Embedding refinement: each node is placed at the centre of the medial
    // ball fitted to the original samples of its cluster (centroid fallback),
    // with the mean distance to those samples as its radius.
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); m];
    for (i, &a) in assign.iter().enumerate() {
        members[a].push(i);
    }
    let nodes: Vec<SkelNode> = members
        .iter()
        .map(|idx| {
            let center = refined_center(points, idx);
            let radius = if idx.is_empty() {
                0.0
            } else {
                idx.iter().map(|&i| dist3(points[i], center)).sum::<f64>() / idx.len() as f64
            };
            SkelNode {
                position: center,
                radius,
                anchor: true,
            }
        })
        .collect();

    // Connect nodes whose samples are neighbours in the original cloud.
    let mut edges = Vec::new();
    for i in 0..n {
        for (j, _) in tree.k_nearest(points[i], KNN.min(n - 1), |j| j != i) {
            let (a, b) = (assign[i], assign[j]);
            if a != b {
                edges.push((a.min(b), a.max(b)));
            }
        }
    }
    edges.sort_unstable();
    edges.dedup();

    let ma = skeletonize(
        &nodes,
        &edges,
        SkeletonParams {
            merge_tolerance: 1e-9 * diagonal,
            prune_scale: 0.5,
            loop_factor: 1.0,
        },
    );
    let mut skeleton = CurveSkeleton {
        nodes: ma
            .points
            .iter()
            .map(|p| CurveSkeletonNode {
                position: p.position,
                radius: p.radius,
                neighbors: Vec::new(),
            })
            .collect(),
    };
    for branch in &ma.branches {
        for w in branch.windows(2) {
            let (a, b) = (w[0], w[1]);
            if a == b {
                continue;
            }
            if !skeleton.nodes[a].neighbors.contains(&b) {
                skeleton.nodes[a].neighbors.push(b);
            }
            if !skeleton.nodes[b].neighbors.contains(&a) {
                skeleton.nodes[b].neighbors.push(a);
            }
        }
    }
    Ok(skeleton)
}

/// Centre of the medial ball of a cluster of original samples.
///
/// Fits the sphere `|p − c|² = ρ²` to the samples by the algebraic
/// (Kåsa) least-squares system `2 p·c + d = |p|²` (`d = ρ² − |c|²`), solved
/// relative to the centroid for conditioning. The fit is rejected in favour
/// of the centroid when the system is singular (fewer than four samples,
/// coplanar samples) or the centre falls outside the cluster's extent (a flat
/// patch, whose best sphere is unboundedly large).
fn refined_center(points: &[[f64; 3]], idx: &[usize]) -> [f64; 3] {
    if idx.is_empty() {
        return [0.0; 3];
    }
    let inv = 1.0 / idx.len() as f64;
    let mut centroid = [0.0_f64; 3];
    for &i in idx {
        for k in 0..3 {
            centroid[k] += points[i][k] * inv;
        }
    }
    if idx.len() < 4 {
        return centroid;
    }
    let extent = idx
        .iter()
        .map(|&i| dist3(points[i], centroid))
        .fold(0.0_f64, f64::max);
    // Normal equations of rows [2x, 2y, 2z, 1] · (c, d) = |x|², x = p − centroid.
    let mut ata = [[0.0_f64; 4]; 4];
    let mut atb = [0.0_f64; 4];
    for &i in idx {
        let x = sub3(points[i], centroid);
        let row = [2.0 * x[0], 2.0 * x[1], 2.0 * x[2], 1.0];
        let rhs = dot3(x, x);
        for r in 0..4 {
            atb[r] += row[r] * rhs;
            for c in 0..4 {
                ata[r][c] += row[r] * row[c];
            }
        }
    }
    let Some(sol) = solve4(ata, atb) else {
        return centroid;
    };
    let offset = [sol[0], sol[1], sol[2]];
    if !(dot3(offset, offset).sqrt() <= extent) {
        return centroid;
    }
    [
        centroid[0] + offset[0],
        centroid[1] + offset[1],
        centroid[2] + offset[2],
    ]
}

/// Gaussian elimination with partial pivoting for a 4×4 system; `None` when
/// the matrix is numerically singular.
fn solve4(mut a: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    let scale = a
        .iter()
        .flat_map(|r| r.iter())
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    if !(scale > 0.0) {
        return None;
    }
    for col in 0..4 {
        let pivot = (col..4).max_by(|&x, &y| a[x][col].abs().total_cmp(&a[y][col].abs()))?;
        if a[pivot][col].abs() <= 1e-12 * scale {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for r in (col + 1)..4 {
            let f = a[r][col] / a[col][col];
            for c in col..4 {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = [0.0_f64; 4];
    for r in (0..4).rev() {
        let mut acc = b[r];
        for c in (r + 1)..4 {
            acc -= a[r][c] * x[c];
        }
        x[r] = acc / a[r][r];
    }
    if x.iter().all(|v| v.is_finite()) {
        Some(x)
    } else {
        None
    }
}

/// One-ring of sample `i`: Delaunay neighbours in the PCA tangent plane among
/// the `k` nearest neighbours, in angular order.
fn one_ring(points: &[[f64; 3]], tree: &KdTree<'_>, i: usize) -> Ring {
    let p = points[i];
    let knn: Vec<usize> = tree
        .k_nearest(p, KNN.min(points.len() - 1), |j| {
            j != i && dist3(points[j], p) > 0.0
        })
        .into_iter()
        .map(|(j, _)| j)
        .collect();
    if knn.is_empty() {
        return Ring {
            members: Vec::new(),
            closed: Vec::new(),
        };
    }
    let normal = pca_normal(points, p, &knn);
    let (e1, e2) = tangent_basis(normal);
    let q: Vec<[f64; 2]> = knn
        .iter()
        .map(|&j| {
            let d = sub3(points[j], p);
            [dot3(d, e1), dot3(d, e2)]
        })
        .collect();

    // Neighbour j borders the Voronoi cell of the origin iff the part of its
    // bisector line x·q_j = |q_j|²/2 satisfying every other half-plane
    // x·q_k ≤ |q_k|²/2 has positive length.
    let mut members: Vec<(f64, usize)> = Vec::new();
    for (a, &qa) in q.iter().enumerate() {
        let len2 = qa[0] * qa[0] + qa[1] * qa[1];
        if len2 <= 0.0 {
            continue;
        }
        let mid = [0.5 * qa[0], 0.5 * qa[1]];
        let dir = [-qa[1], qa[0]];
        let (mut t_lo, mut t_hi) = (f64::NEG_INFINITY, f64::INFINITY);
        for (b, &qb) in q.iter().enumerate() {
            if a == b {
                continue;
            }
            let rhs = 0.5 * (qb[0] * qb[0] + qb[1] * qb[1]) - (mid[0] * qb[0] + mid[1] * qb[1]);
            let coef = dir[0] * qb[0] + dir[1] * qb[1];
            if coef.abs() <= 1e-300 {
                if rhs < 0.0 {
                    t_lo = 1.0;
                    t_hi = 0.0;
                }
                continue;
            }
            let t = rhs / coef;
            if coef > 0.0 {
                t_hi = t_hi.min(t);
            } else {
                t_lo = t_lo.max(t);
            }
        }
        if t_hi - t_lo > 1e-12 * len2.sqrt() {
            members.push((qa[1].atan2(qa[0]), knn[a]));
        }
    }
    members.sort_by(|x, y| x.0.total_cmp(&y.0));
    let m = members.len();
    let closed: Vec<bool> = (0..m)
        .map(|k| {
            if m < 2 {
                return false;
            }
            let a0 = members[k].0;
            let a1 = members[(k + 1) % m].0;
            let mut gap = a1 - a0;
            if gap <= 0.0 {
                gap += 2.0 * std::f64::consts::PI;
            }
            gap < std::f64::consts::PI
        })
        .collect();
    Ring {
        members: members.into_iter().map(|(_, j)| j).collect(),
        closed,
    }
}

/// Smallest-eigenvalue eigenvector of the neighbourhood covariance.
fn pca_normal(points: &[[f64; 3]], p: [f64; 3], knn: &[usize]) -> [f64; 3] {
    let mut mean = p;
    for &j in knn {
        for k in 0..3 {
            mean[k] += points[j][k];
        }
    }
    let cnt = (knn.len() + 1) as f64;
    for v in &mut mean {
        *v /= cnt;
    }
    let mut cov = [[0.0_f64; 3]; 3];
    let mut add = |x: [f64; 3]| {
        let d = sub3(x, mean);
        for r in 0..3 {
            for c in 0..3 {
                cov[r][c] += d[r] * d[c];
            }
        }
    };
    add(p);
    for &j in knn {
        add(points[j]);
    }
    smallest_eigenvector(cov)
}

/// Cyclic Jacobi eigen-decomposition of a symmetric 3×3 matrix; returns the
/// unit eigenvector of the smallest eigenvalue.
fn smallest_eigenvector(mut a: [[f64; 3]; 3]) -> [f64; 3] {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..64 {
        let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
        let scale = a[0][0].abs() + a[1][1].abs() + a[2][2].abs();
        if off <= 1e-15 * scale.max(f64::MIN_POSITIVE) {
            break;
        }
        for (p, q) in [(0usize, 1usize), (0, 2), (1, 2)] {
            if a[p][q].abs() <= f64::MIN_POSITIVE {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            // A ← Jᵀ A J
            for k in 0..3 {
                let akp = a[k][p];
                let akq = a[k][q];
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for k in 0..3 {
                let apk = a[p][k];
                let aqk = a[q][k];
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            for row in &mut v {
                let vp = row[p];
                let vq = row[q];
                row[p] = c * vp - s * vq;
                row[q] = s * vp + c * vq;
            }
        }
    }
    let mut best = 0;
    for k in 1..3 {
        if a[k][k] < a[best][best] {
            best = k;
        }
    }
    unit_or_none([v[0][best], v[1][best], v[2][best]]).unwrap_or([0.0, 0.0, 1.0])
}

/// Orthonormal tangent basis perpendicular to `n`.
fn tangent_basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let helper = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let e1 = unit_or_none(cross3(n, helper)).unwrap_or([0.0, 1.0, 0.0]);
    let e2 = cross3(n, e1);
    (e1, e2)
}

/// Sum of the one-ring triangle areas of sample `i` at positions `x`.
fn ring_area(x: &[[f64; 3]], i: usize, ring: &Ring) -> f64 {
    let m = ring.members.len();
    let mut area = 0.0;
    for k in 0..m {
        if !ring.closed[k] {
            continue;
        }
        let a = ring.members[k];
        let b = ring.members[(k + 1) % m];
        let c = cross3(sub3(x[a], x[i]), sub3(x[b], x[i]));
        area += 0.5 * dot3(c, c).sqrt();
    }
    area
}

/// Clamped cotangent of the angle at `at` in triangle `(at, u, v)`.
fn cot_at(at: [f64; 3], u: [f64; 3], v: [f64; 3]) -> f64 {
    let a = sub3(u, at);
    let b = sub3(v, at);
    let c = cross3(a, b);
    let s = dot3(c, c).sqrt();
    let d = dot3(a, b);
    if s <= 1e-300 {
        return if d > 0.0 { COT_MAX } else { 0.0 };
    }
    (d / s).clamp(0.0, COT_MAX)
}

/// Sparse Laplacian in row form: `rows[i]` lists `(j, w_ij)`; `L_ii = −Σ w_ij`.
struct Laplacian {
    rows: Vec<Vec<(usize, f64)>>,
}

/// Cotangent Laplacian of the one-ring triangles at positions `x`; frozen
/// samples get an empty row (they no longer move under smoothing).
fn cotangent_laplacian(x: &[[f64; 3]], rings: &[Ring], frozen: &[bool]) -> Laplacian {
    let rows = rings
        .iter()
        .enumerate()
        .map(|(i, ring)| {
            if frozen[i] {
                return Vec::new();
            }
            let m = ring.members.len();
            let mut w = vec![0.0_f64; m];
            for k in 0..m {
                if !ring.closed[k] {
                    continue;
                }
                let k1 = (k + 1) % m;
                let (a, b) = (ring.members[k], ring.members[k1]);
                // Angle at b is opposite edge (i, a); angle at a opposite (i, b).
                w[k] += 0.5 * cot_at(x[b], x[i], x[a]);
                w[k1] += 0.5 * cot_at(x[a], x[i], x[b]);
            }
            if w.iter().all(|&v| v <= 0.0) {
                w.iter_mut().for_each(|v| *v = 1.0);
            }
            ring.members.iter().copied().zip(w).collect()
        })
        .collect();
    Laplacian { rows }
}

impl Laplacian {
    /// `y = L x`.
    fn apply(&self, x: &[f64]) -> Vec<f64> {
        self.rows
            .iter()
            .enumerate()
            .map(|(i, row)| row.iter().map(|&(j, w)| w * (x[j] - x[i])).sum())
            .collect()
    }

    /// `y = Lᵀ x`.
    fn apply_transpose(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; x.len()];
        for (i, row) in self.rows.iter().enumerate() {
            for &(j, w) in row {
                y[j] += w * x[i];
                y[i] -= w * x[i];
            }
        }
        y
    }

    /// Diagonal of `LᵀL`.
    fn gram_diagonal(&self, n: usize) -> Vec<f64> {
        let mut d = vec![0.0; n];
        for (i, row) in self.rows.iter().enumerate() {
            let diag: f64 = row.iter().map(|&(_, w)| w).sum();
            d[i] += diag * diag;
            for &(j, w) in row {
                d[j] += w * w;
            }
        }
        d
    }
}

/// Solve `(W_L² LᵀL + W_H²) P' = W_H² P` per coordinate with Jacobi-PCG.
fn solve_contraction(
    lap: &Laplacian,
    w_l: f64,
    w_h: &[f64],
    current: &[[f64; 3]],
) -> Result<Vec<[f64; 3]>> {
    let n = current.len();
    let wl2 = w_l * w_l;
    let gram = lap.gram_diagonal(n);
    let precond: Vec<f64> = (0..n)
        .map(|i| {
            let d = wl2 * gram[i] + w_h[i] * w_h[i];
            if d > 0.0 { 1.0 / d } else { 1.0 }
        })
        .collect();
    let apply = |x: &[f64]| -> Vec<f64> {
        let lx = lap.apply(x);
        let ltlx = lap.apply_transpose(&lx);
        (0..n)
            .map(|i| wl2 * ltlx[i] + w_h[i] * w_h[i] * x[i])
            .collect()
    };
    let mut out = current.to_vec();
    for k in 0..3 {
        let b: Vec<f64> = (0..n).map(|i| w_h[i] * w_h[i] * current[i][k]).collect();
        let x0: Vec<f64> = (0..n).map(|i| current[i][k]).collect();
        let x = pcg(&apply, &b, x0, &precond)?;
        for i in 0..n {
            out[i][k] = x[i];
        }
    }
    Ok(out)
}

/// Preconditioned conjugate gradients for an SPD operator.
fn pcg<F: Fn(&[f64]) -> Vec<f64>>(
    apply: &F,
    b: &[f64],
    mut x: Vec<f64>,
    precond: &[f64],
) -> Result<Vec<f64>> {
    let n = b.len();
    let max_iter = (10 * n).clamp(200, 20_000);
    let b_norm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
    let tol = 1e-10 * b_norm.max(f64::MIN_POSITIVE);
    let ax = apply(&x);
    let mut r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
    let mut z: Vec<f64> = (0..n).map(|i| precond[i] * r[i]).collect();
    let mut p = z.clone();
    let mut rz: f64 = (0..n).map(|i| r[i] * z[i]).sum();
    let mut res = r.iter().map(|v| v * v).sum::<f64>().sqrt();
    for _ in 0..max_iter {
        if res <= tol {
            return Ok(x);
        }
        let ap = apply(&p);
        let pap: f64 = (0..n).map(|i| p[i] * ap[i]).sum();
        if pap <= 0.0 || !pap.is_finite() {
            break;
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        res = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        for i in 0..n {
            z[i] = precond[i] * r[i];
        }
        let rz_new: f64 = (0..n).map(|i| r[i] * z[i]).sum();
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    if res <= 1e-6 * b_norm.max(f64::MIN_POSITIVE) {
        Ok(x)
    } else {
        Err(Error::convergence_failure(
            "laplacian contraction (conjugate gradient)",
            max_iter,
            res,
        ))
    }
}
