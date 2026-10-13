// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Medial-axis construction from unorganised boundary samples, one routine per
//! [`SkeletonMethod`](super::SkeletonMethod).

use std::collections::HashMap;

use super::kdtree::KdTree;
use super::poles::{PoleAnalysis, analyze};
use super::skeleton_graph::{SkelNode, SkeletonParams, cluster_balls, skeletonize};
use super::{MedialAxis, MedialAxisPoint, dist3, dot3, sub3};
use crate::error::{Error, Result};

/// Tuning parameters for [`MedialAxis::from_points_with_config`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MedialAxisConfig {
    /// λ: a medial ball whose radius is below `λ ×` the local sample spacing is
    /// not resolved by the sampling (sharp creases, noise) and is discarded.
    pub min_radius_factor: f64,
    /// Scale-axis factor `s − 1` of the twig pruning: a leaf branch is removed
    /// when all of its balls lie inside some ball of the remaining skeleton
    /// scaled by `1 + prune_scale`. `0` keeps every twig that is not strictly
    /// covered.
    pub prune_scale: f64,
    /// A loop of the skeleton is kept only when it is longer than
    /// `2π · loop_factor ·` the local medial radius (smaller cycles are
    /// sampling artefacts inside one tube).
    pub loop_factor: f64,
    /// Number of voxels along the longest bounding-box side used by
    /// [`SkeletonMethod::Thinning`](super::SkeletonMethod::Thinning).
    pub thinning_resolution: usize,
}

impl Default for MedialAxisConfig {
    fn default() -> Self {
        Self {
            min_radius_factor: 0.5,
            prune_scale: 0.5,
            loop_factor: 1.0,
            thinning_resolution: 64,
        }
    }
}

impl MedialAxisConfig {
    /// Validate the parameters.
    pub(super) fn validate(&self) -> Result<()> {
        if !(self.min_radius_factor.is_finite() && self.min_radius_factor >= 0.0) {
            return Err(Error::out_of_range(
                "min_radius_factor",
                self.min_radius_factor,
                0.0,
                f64::MAX,
            ));
        }
        if !(self.prune_scale.is_finite() && self.prune_scale >= 0.0) {
            return Err(Error::out_of_range(
                "prune_scale",
                self.prune_scale,
                0.0,
                f64::MAX,
            ));
        }
        if !(self.loop_factor.is_finite() && self.loop_factor > 0.0) {
            return Err(Error::out_of_range(
                "loop_factor",
                self.loop_factor,
                f64::MIN_POSITIVE,
                f64::MAX,
            ));
        }
        if !(4..=1024).contains(&self.thinning_resolution) {
            return Err(Error::out_of_range(
                "thinning_resolution",
                self.thinning_resolution as f64,
                4.0,
                1024.0,
            ));
        }
        Ok(())
    }

    /// Skeleton-graph parameters for a shape of the given diagonal.
    pub(super) fn skeleton_params(&self, diagonal: f64) -> SkeletonParams {
        SkeletonParams {
            merge_tolerance: 1e-6 * diagonal,
            prune_scale: self.prune_scale,
            loop_factor: self.loop_factor,
        }
    }
}

/// Error returned when no medial ball survives the resolution filter.
pub(super) fn unresolved() -> Error {
    Error::degenerate_geometry(
        "medial axis: no interior medial ball is resolved by the boundary sampling",
    )
}

/// Voronoi-pole curve skeleton: inner Voronoi vertices linked by Voronoi edges,
/// anchored at the inner poles of the samples.
pub(super) fn voronoi_skeleton(
    boundary: &[[f64; 3]],
    cfg: &MedialAxisConfig,
) -> Result<MedialAxis> {
    let pa = analyze(boundary, cfg.min_radius_factor)?;
    let mut is_anchor = vec![false; pa.tetra.tets.len()];
    for &ti in pa.inner_pole.iter().flatten() {
        is_anchor[ti] = true;
    }
    let mut node_of = vec![usize::MAX; pa.tetra.tets.len()];
    let mut nodes = Vec::new();
    for (ti, sphere) in pa.spheres.iter().enumerate() {
        if let Some(s) = sphere
            && pa.inside[ti]
        {
            node_of[ti] = nodes.len();
            nodes.push(SkelNode {
                position: s.center,
                radius: s.radius,
                anchor: is_anchor[ti],
            });
        }
    }
    if !nodes.iter().any(|n| n.anchor) {
        return Err(unresolved());
    }
    let mut edges = Vec::new();
    for (ti, t) in pa.tetra.tets.iter().enumerate() {
        if node_of[ti] == usize::MAX {
            continue;
        }
        for &nb in &t.neigh {
            if nb < 0 {
                continue;
            }
            let nb = nb as usize;
            if nb > ti && node_of[nb] != usize::MAX {
                edges.push((node_of[ti], node_of[nb]));
            }
        }
    }
    Ok(skeletonize(
        &nodes,
        &edges,
        cfg.skeleton_params(pa.diagonal),
    ))
}

/// Unstructured medial-surface sample: one point per distinct inner pole.
pub(super) fn medial_surface(boundary: &[[f64; 3]], cfg: &MedialAxisConfig) -> Result<MedialAxis> {
    let pa = analyze(boundary, cfg.min_radius_factor)?;
    let mut pole_tets: Vec<usize> = pa.inner_pole.iter().flatten().copied().collect();
    pole_tets.sort_unstable();
    pole_tets.dedup();
    if pole_tets.is_empty() {
        return Err(unresolved());
    }
    // Co-spherical samples give many tets with one circumcentre: merge
    // centres closer than `tol`, found through a hash grid of cell size `tol`
    // anchored at the bounding-box minimum.
    let tol = (1e-6 * pa.diagonal).max(f64::MIN_POSITIVE);
    let mut lo = [f64::INFINITY; 3];
    for p in boundary {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
        }
    }
    let cell_of = |c: [f64; 3]| -> [i64; 3] {
        let mut key = [0_i64; 3];
        for k in 0..3 {
            // Inner poles lie inside the convex hull, so the quotient is
            // bounded by about `1e6`; the clamp only guards pathological input.
            key[k] = ((c[k] - lo[k]) / tol).floor().clamp(-4e15, 4e15) as i64;
        }
        key
    };
    let mut grid: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    let mut points: Vec<MedialAxisPoint> = Vec::new();
    let mut counts: Vec<f64> = Vec::new();
    let mut seeds: Vec<[f64; 3]> = Vec::new();
    for ti in pole_tets {
        let Some(s) = pa.spheres[ti] else { continue };
        let key = cell_of(s.center);
        let mut existing = None;
        'search: for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let probe = [key[0] + dx, key[1] + dy, key[2] + dz];
                    if let Some(list) = grid.get(&probe)
                        && let Some(&k) = list.iter().find(|&&k| dist3(seeds[k], s.center) <= tol)
                    {
                        existing = Some(k);
                        break 'search;
                    }
                }
            }
        }
        match existing {
            Some(k) => {
                let c = counts[k];
                let p = &mut points[k];
                for d in 0..3 {
                    p.position[d] = (p.position[d] * c + s.center[d]) / (c + 1.0);
                }
                p.radius = (p.radius * c + s.radius) / (c + 1.0);
                p.feature_size = p.radius;
                counts[k] = c + 1.0;
            }
            None => {
                grid.entry(key).or_default().push(seeds.len());
                seeds.push(s.center);
                counts.push(1.0);
                points.push(MedialAxisPoint::new(s.center, s.radius, [0.0; 3], s.radius));
            }
        }
    }
    Ok(MedialAxis {
        points,
        branches: Vec::new(),
        junctions: Vec::new(),
    })
}

/// Shrinking-ball medial ball of every sample (Ma, Bae & Choi, "3D medial axis
/// point approximation using nearest neighbors and the normal field", 2012),
/// linked through surface-neighbour relations and skeletonised.
pub(super) fn ball_rolling(boundary: &[[f64; 3]], cfg: &MedialAxisConfig) -> Result<MedialAxis> {
    let pa = analyze(boundary, cfg.min_radius_factor)?;
    let balls = shrinking_balls(boundary, &pa, cfg.min_radius_factor);
    let mut node_of = vec![usize::MAX; boundary.len()];
    let mut nodes = Vec::new();
    for (i, ball) in balls.iter().enumerate() {
        if let Some((center, radius)) = *ball {
            node_of[i] = nodes.len();
            nodes.push(SkelNode {
                position: center,
                radius,
                anchor: true,
            });
        }
    }
    if nodes.is_empty() {
        return Err(unresolved());
    }
    // Surface-neighbour relation: short Delaunay edges between samples whose
    // normals agree. This never links samples facing each other across a thin
    // gap or across the solid.
    let mut edges = Vec::new();
    for i in 0..boundary.len() {
        if node_of[i] == usize::MAX {
            continue;
        }
        let Some(ni) = pa.normals[i] else { continue };
        for &j in &pa.neighbors[i] {
            if j <= i || node_of[j] == usize::MAX {
                continue;
            }
            let Some(nj) = pa.normals[j] else { continue };
            let reach = 4.0 * pa.spacing[i].max(pa.spacing[j]);
            if dot3(ni, nj) > 0.0 && dist3(boundary[i], boundary[j]) <= reach {
                edges.push((node_of[i], node_of[j]));
            }
        }
    }
    // Shrinking-ball centres scatter around the medial axis at the scale of the
    // sampling; cluster them at twice the mean sample spacing so the curve
    // skeleton follows the axis instead of zig-zagging through the scatter.
    let finite: Vec<f64> = pa
        .spacing
        .iter()
        .copied()
        .filter(|s| s.is_finite())
        .collect();
    let mean_spacing = finite.iter().sum::<f64>() / finite.len().max(1) as f64;
    let (nodes, edges) = cluster_balls(&nodes, &edges, 2.0 * mean_spacing);
    Ok(skeletonize(
        &nodes,
        &edges,
        cfg.skeleton_params(pa.diagonal),
    ))
}

/// Separation angle below which a first shrinking step is attributed to a
/// locally planar, under-resolved neighbourhood (Ma et al. 2012 denoising).
const DENOISE_PLANAR: f64 = 32.0 * std::f64::consts::PI / 180.0;
/// Separation angle below which a later shrinking step is attributed to noise
/// and the previous (larger) ball is kept.
const DENOISE_PRESERVE: f64 = 20.0 * std::f64::consts::PI / 180.0;

/// Shrinking-ball iteration for every sample with an oriented normal.
///
/// The ball touching sample `p` with inward normal `−n` is
/// `c(r) = p − n r`. Starting from a radius equal to the bounding-box
/// diagonal, the nearest sample `q` to the centre is found; if `q` lies inside
/// the ball, the radius shrinks to that of the ball through `p` and `q`,
/// `r = |p − q|² / (2 n·(p − q))`, until the ball is empty.
///
/// Denoising follows Ma et al.: the *separation angle* `∠(p, c, q)` of a new
/// ball measures how far apart its two touching samples are. A first step
/// with an angle below 32° means the sample's neighbourhood is locally planar
/// at the sampling scale and no medial ball is resolved; a later step with an
/// angle below 20° is caused by noise, so the previous ball is kept.
fn shrinking_balls(
    boundary: &[[f64; 3]],
    pa: &PoleAnalysis,
    min_radius_factor: f64,
) -> Vec<Option<([f64; 3], f64)>> {
    let used: Vec<usize> = (0..boundary.len())
        .filter(|&i| pa.normals[i].is_some())
        .collect();
    let tree = KdTree::build(boundary, &used);
    let eps = 1e-12 * pa.diagonal.max(f64::MIN_POSITIVE);
    (0..boundary.len())
        .map(|i| {
            let n = pa.normals[i]?;
            let p = boundary[i];
            let initial = pa.diagonal;
            let mut r = initial;
            for iteration in 0..256 {
                let c = [p[0] - n[0] * r, p[1] - n[1] * r, p[2] - n[2] * r];
                let (q_idx, d2) = tree.nearest_filtered(c, |j| dist3(boundary[j], p) > eps)?;
                if d2.sqrt() >= r - eps {
                    break;
                }
                let q = boundary[q_idx];
                let pq = sub3(p, q);
                let denom = 2.0 * dot3(n, pq);
                if denom <= 0.0 {
                    break;
                }
                let r_new = dot3(pq, pq) / denom;
                if !(r_new.is_finite() && r_new < r) {
                    break;
                }
                let c_new = [
                    p[0] - n[0] * r_new,
                    p[1] - n[1] * r_new,
                    p[2] - n[2] * r_new,
                ];
                let separation = angle_between(sub3(p, c_new), sub3(q, c_new));
                if iteration == 0 && separation < DENOISE_PLANAR {
                    return None;
                }
                if iteration > 0 && separation < DENOISE_PRESERVE {
                    break;
                }
                r = r_new;
            }
            if r >= initial || r < min_radius_factor * pa.spacing[i] {
                return None;
            }
            Some(([p[0] - n[0] * r, p[1] - n[1] * r, p[2] - n[2] * r], r))
        })
        .collect()
}

/// Angle between two vectors in radians (0 when either is zero).
fn angle_between(a: [f64; 3], b: [f64; 3]) -> f64 {
    let na = dot3(a, a).sqrt();
    let nb = dot3(b, b).sqrt();
    if na <= 0.0 || nb <= 0.0 {
        return 0.0;
    }
    (dot3(a, b) / (na * nb)).clamp(-1.0, 1.0).acos()
}
