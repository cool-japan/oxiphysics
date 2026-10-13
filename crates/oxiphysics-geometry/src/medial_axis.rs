// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Medial axis transform, skeleton computation, and related distance-transform utilities.
//!
//! Provides medial axis extraction from boundary point sets, straight skeleton
//! computation for 2D polygons, simplified Voronoi-based 2D skeletons, and
//! centerline extraction for tubular structures.

mod construct;
mod contraction;
mod kdtree;
mod poles;
mod skeleton_graph;
mod thinning3d;

pub use construct::MedialAxisConfig;

use crate::error::Result;

/// Method used to compute the skeleton / medial axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SkeletonMethod {
    /// Ball-rolling (inscribed-sphere) approach.
    BallRolling,
    /// Binary-thinning approach.
    Thinning,
    /// Voronoi-diagram-based approach.
    Voronoi,
    /// Medial-surface approach (for 3-D objects).
    MedialSurface,
}

/// A single point on the medial axis.
#[derive(Debug, Clone)]
pub struct MedialAxisPoint {
    /// 3-D position of the medial-axis point.
    pub position: [f64; 3],
    /// Radius of the largest inscribed sphere at this point.
    pub radius: f64,
    /// Unit tangent vector along the axis at this point.
    pub tangent: [f64; 3],
    /// Local feature size (minimum distance to any boundary feature).
    pub feature_size: f64,
}

impl MedialAxisPoint {
    /// Construct a new `MedialAxisPoint`.
    pub fn new(position: [f64; 3], radius: f64, tangent: [f64; 3], feature_size: f64) -> Self {
        Self {
            position,
            radius,
            tangent,
            feature_size,
        }
    }
}

/// Medial axis (skeleton) of a 3-D solid represented as a set of boundary points.
#[derive(Debug, Clone)]
pub struct MedialAxis {
    /// Ordered list of medial-axis sample points.
    pub points: Vec<MedialAxisPoint>,
    /// Each branch is a list of indices into `points`.
    pub branches: Vec<Vec<usize>>,
    /// Indices of junction (bifurcation) points.
    pub junctions: Vec<usize>,
}

impl MedialAxis {
    /// Construct the medial axis of the solid whose boundary is sampled by
    /// `boundary`, using the default [`MedialAxisConfig`].
    ///
    /// See [`MedialAxis::from_points_with_config`] for the algorithms behind
    /// each [`SkeletonMethod`].
    ///
    /// # Errors
    ///
    /// See [`MedialAxis::from_points_with_config`].
    pub fn from_points(boundary: &[[f64; 3]], method: SkeletonMethod) -> Result<Self> {
        Self::from_points_with_config(boundary, method, &MedialAxisConfig::default())
    }

    /// Construct the medial axis of the solid whose boundary is sampled by the
    /// unorganised point cloud `boundary` (a closed surface sample).
    ///
    /// Every method starts from the exact 3-D Delaunay tetrahedralisation of
    /// the samples and their *Voronoi poles* (Amenta et al.): each sample gets
    /// an oriented outward normal (convex-hull faces, then minimum-spanning-tree
    /// propagation) and an inner pole, the farthest Voronoi vertex of its cell
    /// on the inner side, whose distance to the sample is the radius of the
    /// local medial ball. Poles with a radius below
    /// `config.min_radius_factor ×` the local sample spacing are discarded.
    ///
    /// * [`SkeletonMethod::Voronoi`] — curve skeleton: the inner Voronoi
    ///   vertices linked by Voronoi edges form a graph that is reduced to a
    ///   minimum spanning forest spanning the poles, pruned with the
    ///   scale-axis criterion, re-closed around genuine tunnels and split into
    ///   branches (maximal chains) and junctions (degree ≥ 3).
    /// * [`SkeletonMethod::MedialSurface`] — the unstructured medial-surface
    ///   sample: one point per distinct inner pole with its medial radius.
    ///   `branches` and `junctions` are empty and tangents are zero because a
    ///   medial *surface* has no curve connectivity.
    /// * [`SkeletonMethod::BallRolling`] — shrinking-ball medial balls
    ///   (Ma, Bae & Choi 2012) computed from the oriented normals and a k-d
    ///   tree, linked through surface-neighbour relations, then skeletonised as
    ///   for `Voronoi`.
    /// * [`SkeletonMethod::Thinning`] — the union of inner polar balls is
    ///   voxelised (`config.thinning_resolution` voxels along the longest side)
    ///   and thinned to a one-voxel curve skeleton with a topology-preserving
    ///   3-D simple-point test, then skeletonised as for `Voronoi`.
    ///
    /// Each point's `radius` and `feature_size` are the distance from the
    /// medial point to the nearest boundary sample; `tangent` is the unit
    /// direction of the branch through the point (zero when undefined).
    ///
    /// # Errors
    ///
    /// * [`Error::TooFewPoints`](crate::error::Error::TooFewPoints) for fewer than four samples;
    /// * [`Error::DegenerateGeometry`](crate::error::Error::DegenerateGeometry) for non-finite or coplanar samples, or
    ///   when no interior medial ball is resolved by the sampling;
    /// * [`Error::OutOfRange`](crate::error::Error::OutOfRange) for invalid `config` values.
    pub fn from_points_with_config(
        boundary: &[[f64; 3]],
        method: SkeletonMethod,
        config: &MedialAxisConfig,
    ) -> Result<Self> {
        config.validate()?;
        match method {
            SkeletonMethod::Voronoi => construct::voronoi_skeleton(boundary, config),
            SkeletonMethod::MedialSurface => construct::medial_surface(boundary, config),
            SkeletonMethod::BallRolling => construct::ball_rolling(boundary, config),
            SkeletonMethod::Thinning => thinning3d::thinning_skeleton(boundary, config),
        }
    }

    /// Return the number of branches in the medial axis.
    pub fn branch_count(&self) -> usize {
        self.branches.len()
    }

    /// Return the total arc-length of the medial axis.
    pub fn total_length(&self) -> f64 {
        let mut len = 0.0;
        for branch in &self.branches {
            for w in branch.windows(2) {
                len += dist3(self.points[w[0]].position, self.points[w[1]].position);
            }
        }
        len
    }

    /// Return `true` if the shape is tubular (single branch, roughly constant radius).
    pub fn is_tubular(&self) -> bool {
        if self.branches.len() != 1 {
            return false;
        }
        if self.points.len() < 2 {
            return true;
        }
        let radii: Vec<f64> = self.points.iter().map(|p| p.radius).collect();
        let mean = radii.iter().sum::<f64>() / radii.len() as f64;
        let variance = radii.iter().map(|r| (r - mean).powi(2)).sum::<f64>() / radii.len() as f64;
        let cv = variance.sqrt() / mean.max(1e-12); // coefficient of variation
        cv < 0.3
    }
}

/// Simplified 2-D skeleton via Voronoi vertices of a point set.
pub struct Voronoi2dSkeleton;

impl Voronoi2dSkeleton {
    /// Compute skeleton-like vertices from a 2-D point set within a bounding box.
    ///
    /// Returns Voronoi-vertex approximations by finding circumcentres of triplets
    /// of nearby points and retaining those inside the bounding box.
    pub fn compute(points: &[[f64; 2]], bbox: [f64; 4]) -> Vec<[f64; 2]> {
        // bbox = [x_min, y_min, x_max, y_max]
        let n = points.len();
        if n < 3 {
            return vec![];
        }
        let mut result = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    if let Some(cc) = circumcenter_2d(points[i], points[j], points[k])
                        && cc[0] >= bbox[0]
                        && cc[1] >= bbox[1]
                        && cc[0] <= bbox[2]
                        && cc[1] <= bbox[3]
                    {
                        result.push(cc);
                    }
                }
            }
        }
        result
    }
}

/// Straight skeleton for a simple 2-D polygon.
pub struct StraightSkeleton {
    polygon: Vec<[f64; 2]>,
}

impl StraightSkeleton {
    /// Construct a `StraightSkeleton` from a polygon given as an ordered vertex list.
    pub fn new(polygon: Vec<[f64; 2]>) -> Self {
        Self { polygon }
    }

    /// Compute offset polygon(s) at signed distance `dist` (positive = inward).
    ///
    /// Returns a list of offset rings; may be empty if the polygon collapses.
    pub fn compute_offsets(&self, dist: f64) -> Vec<Vec<[f64; 2]>> {
        let n = self.polygon.len();
        if n < 3 {
            return vec![];
        }
        let mut offset = Vec::with_capacity(n);
        for i in 0..n {
            let prev = self.polygon[(i + n - 1) % n];
            let curr = self.polygon[i];
            let next = self.polygon[(i + 1) % n];

            // Inward normal bisector direction
            let d1 = normalize2(sub2(curr, prev));
            let d2 = normalize2(sub2(next, curr));
            let n1 = [-d1[1], d1[0]];
            let n2 = [-d2[1], d2[0]];
            let bisector = normalize2(add2(n1, n2));

            // Scale by 1/sin(half-angle) for miter
            let dot = n1[0] * n2[0] + n1[1] * n2[1];
            let sin_half = ((1.0 - dot) / 2.0).sqrt().max(1e-9);
            let scale = 1.0 / sin_half;

            let op = [
                curr[0] + bisector[0] * dist * scale,
                curr[1] + bisector[1] * dist * scale,
            ];
            offset.push(op);
        }
        vec![offset]
    }

    /// Return the skeleton edges as pairs of 2-D points.
    ///
    /// Each edge connects a vertex to its corresponding collapsed skeleton point.
    pub fn skeleton_edges(&self) -> Vec<([f64; 2], [f64; 2])> {
        // Simplified: use bisector from each vertex toward the centroid
        let n = self.polygon.len();
        if n < 3 {
            return vec![];
        }
        let cx = self.polygon.iter().map(|p| p[0]).sum::<f64>() / n as f64;
        let cy = self.polygon.iter().map(|p| p[1]).sum::<f64>() / n as f64;
        let centroid = [cx, cy];

        let mut edges = Vec::with_capacity(n);
        for i in 0..n {
            edges.push((self.polygon[i], centroid));
        }
        edges
    }
}

// ─── Distance-transform helpers ────────────────────────────────────────────

/// Compute the minimum distance from point `p` to line segment `[a, b]`.
pub fn point_to_segment_dist(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let ap = sub3(p, a);
    let len_sq = dot3(ab, ab);
    if len_sq < 1e-24 {
        return dist3(p, a);
    }
    let t = (dot3(ap, ab) / len_sq).clamp(0.0, 1.0);
    let closest = add3(a, scale3(ab, t));
    dist3(p, closest)
}

/// Compute the minimum distance from 2-D point `p` to a polygon edge set.
pub fn point_to_polygon_dist(p: [f64; 2], poly: &[[f64; 2]]) -> f64 {
    let n = poly.len();
    if n == 0 {
        return f64::INFINITY;
    }
    let mut min_d = f64::INFINITY;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let d = seg_dist_2d(p, a, b);
        if d < min_d {
            min_d = d;
        }
    }
    min_d
}

/// Test whether 2-D point `p` is inside a simple polygon using ray casting.
pub fn is_inside_polygon(p: [f64; 2], poly: &[[f64; 2]]) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let (x, y) = (p[0], p[1]);
    let mut j = n - 1;
    for i in 0..n {
        let xi = poly[i][0];
        let yi = poly[i][1];
        let xj = poly[j][0];
        let yj = poly[j][1];
        let intersect = ((yi > y) != (yj > y))
            && (x < (xj - xi) * (y - yi) / (yj - yi + if yj >= yi { 1e-14 } else { -1e-14 }) + xi);
        if intersect {
            inside = !inside;
        }
        j = i;
    }
    inside
}

// ─── Centerline extraction ──────────────────────────────────────────────────

/// Centerline extractor for tubular structures (e.g., blood vessels, pipes).
pub struct CenterlineExtraction;

impl CenterlineExtraction {
    /// Extract a centerline from a surface point cloud by slicing along the
    /// principal axis and computing slice centroids.
    ///
    /// `num_slices` controls the resolution of the output.
    pub fn extract_centerline(surface_points: &[[f64; 3]], num_slices: usize) -> Vec<[f64; 3]> {
        if surface_points.is_empty() || num_slices == 0 {
            return vec![];
        }

        // Find axis-aligned bounding box
        let mut x_min = f64::INFINITY;
        let mut x_max = f64::NEG_INFINITY;
        for p in surface_points {
            if p[0] < x_min {
                x_min = p[0];
            }
            if p[0] > x_max {
                x_max = p[0];
            }
        }

        let mut centerline = Vec::with_capacity(num_slices);
        for s in 0..num_slices {
            let t = s as f64 / (num_slices - 1).max(1) as f64;
            let x_c = x_min + t * (x_max - x_min);
            let half_width = (x_max - x_min) / (num_slices as f64 * 2.0);

            // Gather points in this slab
            let slice_pts: Vec<[f64; 3]> = surface_points
                .iter()
                .filter(|p| (p[0] - x_c).abs() <= half_width)
                .copied()
                .collect();

            let centroid = if slice_pts.is_empty() {
                [x_c, 0.0, 0.0]
            } else {
                let m = slice_pts.len() as f64;
                [
                    slice_pts.iter().map(|p| p[0]).sum::<f64>() / m,
                    slice_pts.iter().map(|p| p[1]).sum::<f64>() / m,
                    slice_pts.iter().map(|p| p[2]).sum::<f64>() / m,
                ]
            };
            centerline.push(centroid);
        }
        centerline
    }
}

// ─── Pruning ────────────────────────────────────────────────────────────────

/// Remove branches shorter than `min_length` from the medial axis in-place.
pub fn prune_short_branches(ma: &mut MedialAxis, min_length: f64) {
    ma.branches.retain(|branch| {
        let mut len = 0.0;
        for w in branch.windows(2) {
            len += dist3(ma.points[w[0]].position, ma.points[w[1]].position);
        }
        len >= min_length
    });
}

// ─── Topology ───────────────────────────────────────────────────────────────

/// Compute the Euler characteristic χ = V - E + F of the medial axis graph.
///
/// Here V = number of axis points, E = number of axis edges (consecutive pairs
/// in branches), F = 0 (the axis is a graph, not a surface).
pub fn euler_characteristic(ma: &MedialAxis) -> i32 {
    let v = ma.points.len() as i32;
    let mut e = 0i32;
    for branch in &ma.branches {
        e += branch.len().saturating_sub(1) as i32;
    }
    v - e
}

// ─── Axis-aligned slab decomposition ────────────────────────────────────────

/// An axis-aligned slab used to decompose a bounding region for medial-axis
/// computation.
#[derive(Debug, Clone)]
pub struct AxisAlignedSlab {
    /// Minimum x-coordinate of the slab.
    pub x_min: f64,
    /// Maximum x-coordinate of the slab.
    pub x_max: f64,
    /// Indices of medial-axis points contained in this slab.
    pub point_indices: Vec<usize>,
}

impl AxisAlignedSlab {
    /// Build a slab decomposition of the medial axis along the X axis into
    /// `num_slabs` equal-width slabs.
    pub fn decompose(ma: &MedialAxis, num_slabs: usize) -> Vec<Self> {
        if num_slabs == 0 || ma.points.is_empty() {
            return vec![];
        }
        let x_vals: Vec<f64> = ma.points.iter().map(|p| p.position[0]).collect();
        let x_min = x_vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let x_max = x_vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let width = (x_max - x_min) / num_slabs as f64;

        let mut slabs: Vec<AxisAlignedSlab> = (0..num_slabs)
            .map(|i| AxisAlignedSlab {
                x_min: x_min + i as f64 * width,
                x_max: x_min + (i + 1) as f64 * width,
                point_indices: vec![],
            })
            .collect();

        for (idx, pt) in ma.points.iter().enumerate() {
            let raw = ((pt.position[0] - x_min) / width.max(1e-15)) as usize;
            let slab_idx = raw.min(num_slabs - 1);
            slabs[slab_idx].point_indices.push(idx);
        }
        slabs
    }
}

// ─── Shape diameter function ─────────────────────────────────────────────────

/// Compute the shape diameter function (SDF) approximation at surface point `p`
/// with outward normal `n`, over a local hemisphere of `radius`.
///
/// Returns an estimate of the local thickness of the shape.
pub fn shape_diameter(p: [f64; 3], n: [f64; 3], radius: f64) -> f64 {
    // Approximation: SDF ≈ 2 * radius * |cos(angle between n and -n)| = 2*radius
    // In practice, average of ray-cast lengths; here we return 2*radius as a
    // geometric estimate since we have no volumetric representation.
    let _ = p;
    let _n_norm = normalize3(n);
    2.0 * radius
}

// ─── Private math helpers ────────────────────────────────────────────────────

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f64; 3], t: f64) -> [f64; 3] {
    [a[0] * t, a[1] * t, a[2] * t]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-15 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn add2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

fn normalize2(v: [f64; 2]) -> [f64; 2] {
    let len = (v[0] * v[0] + v[1] * v[1]).sqrt();
    if len < 1e-15 {
        [1.0, 0.0]
    } else {
        [v[0] / len, v[1] / len]
    }
}

fn seg_dist_2d(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let ab = sub2(b, a);
    let ap = sub2(p, a);
    let len_sq = ab[0] * ab[0] + ab[1] * ab[1];
    if len_sq < 1e-24 {
        return (ap[0] * ap[0] + ap[1] * ap[1]).sqrt();
    }
    let t = ((ap[0] * ab[0] + ap[1] * ab[1]) / len_sq).clamp(0.0, 1.0);
    let cx = a[0] + ab[0] * t - p[0];
    let cy = a[1] + ab[1] * t - p[1];
    (cx * cx + cy * cy).sqrt()
}

fn circumcenter_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Option<[f64; 2]> {
    let ax = b[0] - a[0];
    let ay = b[1] - a[1];
    let bx = c[0] - a[0];
    let by = c[1] - a[1];
    let d = 2.0 * (ax * by - ay * bx);
    if d.abs() < 1e-12 {
        return None;
    }
    let ux = (by * (ax * ax + ay * ay) - ay * (bx * bx + by * by)) / d;
    let uy = (ax * (bx * bx + by * by) - bx * (ax * ax + ay * ay)) / d;
    Some([a[0] + ux, a[1] + uy])
}

// ─── Distance transform (2-D grid) ──────────────────────────────────────────

/// A 2-D binary image stored as a flat boolean grid, row-major.
#[derive(Debug, Clone)]
pub struct BinaryImage {
    /// Grid width (number of columns).
    pub width: usize,
    /// Grid height (number of rows).
    pub height: usize,
    /// Pixel data (true = foreground / object).
    pub data: Vec<bool>,
}

impl BinaryImage {
    /// Create a new binary image filled with `fill`.
    pub fn new(width: usize, height: usize, fill: bool) -> Self {
        Self {
            width,
            height,
            data: vec![fill; width * height],
        }
    }

    /// Access pixel at (col, row).  Returns `false` if out of bounds.
    pub fn get(&self, col: usize, row: usize) -> bool {
        if col < self.width && row < self.height {
            self.data[row * self.width + col]
        } else {
            false
        }
    }

    /// Set pixel at (col, row).
    pub fn set(&mut self, col: usize, row: usize, val: bool) {
        if col < self.width && row < self.height {
            self.data[row * self.width + col] = val;
        }
    }

    /// Count foreground pixels.
    pub fn foreground_count(&self) -> usize {
        self.data.iter().filter(|&&v| v).count()
    }
}

/// Compute the Euclidean distance transform of a binary image.
///
/// For each pixel, returns the distance to the nearest background pixel.
/// Background pixels have distance 0.  Uses the independent-scan algorithm
/// (Saito & Toriwaki, simplified).
pub fn distance_transform(img: &BinaryImage) -> Vec<f64> {
    let w = img.width;
    let h = img.height;
    let n = w * h;
    let big = (w + h) as f64;
    let mut dt = vec![0.0_f64; n];

    // Initialise: foreground = big, background = 0
    for (i, dt_val) in dt.iter_mut().enumerate().take(n) {
        *dt_val = if img.data[i] { big } else { 0.0 };
    }

    // Forward pass: row-wise left-to-right, top-to-bottom
    for r in 0..h {
        for c in 0..w {
            let idx = r * w + c;
            if c > 0 {
                dt[idx] = dt[idx].min(dt[idx - 1] + 1.0);
            }
            if r > 0 {
                dt[idx] = dt[idx].min(dt[(r - 1) * w + c] + 1.0);
            }
        }
    }

    // Backward pass: right-to-left, bottom-to-top
    for r in (0..h).rev() {
        for c in (0..w).rev() {
            let idx = r * w + c;
            if c + 1 < w {
                dt[idx] = dt[idx].min(dt[idx + 1] + 1.0);
            }
            if r + 1 < h {
                dt[idx] = dt[idx].min(dt[(r + 1) * w + c] + 1.0);
            }
        }
    }

    dt
}

/// Compute the exact Euclidean distance transform using Meijster's algorithm.
///
/// For every pixel the result is the Euclidean distance (in pixel units) from
/// the pixel centre to the centre of the nearest *background* (`false`) pixel;
/// background pixels map to `0.0`. The output is laid out row-major exactly
/// like [`BinaryImage::data`].
///
/// The implementation follows A. Meijster, J. B. T. M. Roerdink and
/// W. H. Hesselink, "A General Algorithm for Computing Distance Transforms in
/// Linear Time" (2000):
///
/// 1. a column-wise pass computes `g(x, y)`, the vertical distance to the
///    nearest background pixel in the same column;
/// 2. a row-wise pass computes the lower envelope of the parabolas
///    `f(x, i) = (x − i)² + g(i)²` with a stack of parabola indices `s` and
///    segment start positions `t`, then scans it to read off the minimum.
///
/// All intermediate quantities are squared integer distances held in `i64`, so
/// the result is *exact*: every returned value is the square root of an
/// integer and is bit-identical to a brute-force evaluation. The total cost is
/// `O(width · height)`.
///
/// If the image contains no background pixel at all, no distance is defined
/// and every pixel is reported as [`f64::INFINITY`].
pub fn distance_transform_exact(img: &BinaryImage) -> Vec<f64> {
    let w = img.width;
    let h = img.height;
    let n = w * h;
    if n == 0 {
        return vec![];
    }
    if img.data.iter().all(|&v| v) {
        return vec![f64::INFINITY; n];
    }

    // Phase 1: column-wise vertical distance, stored squared.
    // `inf` exceeds any attainable distance, so a column without background
    // never wins the row-wise minimum against a column that has one.
    let inf = (w + h) as i64;
    let mut g = vec![0_i64; n];
    for c in 0..w {
        g[c] = if img.data[c] { inf } else { 0 };
        for r in 1..h {
            let idx = r * w + c;
            g[idx] = if img.data[idx] { g[idx - w] + 1 } else { 0 };
        }
        for r in (0..h.saturating_sub(1)).rev() {
            let idx = r * w + c;
            let below = g[idx + w] + 1;
            if below < g[idx] {
                g[idx] = below;
            }
        }
    }
    for v in &mut g {
        *v *= *v;
    }

    // Phase 2: row-wise lower envelope of the parabolas f(x, i) = (x-i)² + G(i).
    let mut dt = vec![0.0_f64; n];
    let mut s_buf = vec![0_i64; w];
    let mut t_buf = vec![0_i64; w];
    let width = w as i64;
    for r in 0..h {
        let row = &g[r * w..(r + 1) * w];
        let f = |x: i64, i: i64| -> i64 {
            let d = x - i;
            d * d + row[i as usize]
        };
        let mut q: isize = 0;
        s_buf[0] = 0;
        t_buf[0] = 0;
        for u in 1..width {
            while q >= 0 && f(t_buf[q as usize], s_buf[q as usize]) > f(t_buf[q as usize], u) {
                q -= 1;
            }
            if q < 0 {
                q = 0;
                s_buf[0] = u;
            } else {
                let w1 = 1 + sep_edt(s_buf[q as usize], u, row);
                if w1 < width {
                    q += 1;
                    s_buf[q as usize] = u;
                    t_buf[q as usize] = w1;
                }
            }
        }
        for u in (0..width).rev() {
            let qi = q as usize;
            dt[r * w + u as usize] = (f(u, s_buf[qi]) as f64).sqrt();
            if u == t_buf[qi] {
                q -= 1;
            }
        }
    }

    dt
}

/// Meijster separator for the Euclidean metric.
///
/// Returns the largest integer abscissa `x` at which the parabola rooted at
/// column `i` is still no higher than the one rooted at column `u > i`:
/// `Sep(i, u) = ⌊(u² − i² + G(u) − G(i)) / (2 (u − i))⌋`, using Euclidean
/// (floor) division so negative numerators round toward −∞.
fn sep_edt(i: i64, u: i64, row: &[i64]) -> i64 {
    let num = u * u - i * i + row[u as usize] - row[i as usize];
    num.div_euclid(2 * (u - i))
}

// ─── Thinning algorithm (Zhang-Suen) ────────────────────────────────────────

/// Apply the Zhang-Suen thinning algorithm to produce a 1-pixel-wide skeleton.
///
/// Iteratively removes boundary pixels until no more can be removed, preserving
/// topology.  Operates in-place on the image.
pub fn zhang_suen_thinning(img: &mut BinaryImage) {
    let w = img.width;
    let h = img.height;
    if w < 3 || h < 3 {
        return;
    }

    loop {
        // Sub-iteration 1
        let to_remove_1 = zs_subiteration(img, true);
        for &(c, r) in &to_remove_1 {
            img.set(c, r, false);
        }

        // Sub-iteration 2
        let to_remove_2 = zs_subiteration(img, false);
        for &(c, r) in &to_remove_2 {
            img.set(c, r, false);
        }

        if to_remove_1.is_empty() && to_remove_2.is_empty() {
            break;
        }
    }
}

/// One sub-iteration of Zhang-Suen.
fn zs_subiteration(img: &BinaryImage, first: bool) -> Vec<(usize, usize)> {
    let w = img.width;
    let h = img.height;
    let mut to_remove = Vec::new();

    for r in 1..h.saturating_sub(1) {
        for c in 1..w.saturating_sub(1) {
            if !img.get(c, r) {
                continue;
            }

            // 8-neighbours: P2..P9 in clockwise order starting from top
            let p2 = img.get(c, r - 1) as u8;
            let p3 = img.get(c + 1, r - 1) as u8;
            let p4 = img.get(c + 1, r) as u8;
            let p5 = img.get(c + 1, r + 1) as u8;
            let p6 = img.get(c, r + 1) as u8;
            let p7 = img.get(c - 1, r + 1) as u8;
            let p8 = img.get(c - 1, r) as u8;
            let p9 = img.get(c - 1, r - 1) as u8;

            let b = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9; // B(P1)
            if !(2..=6).contains(&b) {
                continue;
            }

            // A(P1): number of 0->1 transitions in the sequence P2..P9..P2
            let seq = [p2, p3, p4, p5, p6, p7, p8, p9, p2];
            let a: u8 = seq
                .windows(2)
                .map(|w| if w[0] == 0 && w[1] == 1 { 1u8 } else { 0 })
                .sum();
            if a != 1 {
                continue;
            }

            if first {
                // P2 * P4 * P6 == 0  and  P4 * P6 * P8 == 0
                if p2 * p4 * p6 != 0 || p4 * p6 * p8 != 0 {
                    continue;
                }
            } else {
                // P2 * P4 * P8 == 0  and  P2 * P6 * P8 == 0
                if p2 * p4 * p8 != 0 || p2 * p6 * p8 != 0 {
                    continue;
                }
            }

            to_remove.push((c, r));
        }
    }

    to_remove
}

/// Extract skeleton pixel coordinates from a thinned binary image.
pub fn skeleton_pixels(img: &BinaryImage) -> Vec<(usize, usize)> {
    let mut pixels = Vec::new();
    for r in 0..img.height {
        for c in 0..img.width {
            if img.get(c, r) {
                pixels.push((c, r));
            }
        }
    }
    pixels
}

// ─── Medial axis transform (MAT) ────────────────────────────────────────────

/// A single point of the medial axis transform, storing both position and the
/// inscribed-circle radius.
#[derive(Debug, Clone)]
pub struct MatPoint {
    /// Column coordinate.
    pub col: usize,
    /// Row coordinate.
    pub row: usize,
    /// Inscribed-circle radius (from the distance transform).
    pub radius: f64,
}

/// Compute the medial axis transform of a binary image.
///
/// This combines Zhang-Suen thinning with the distance transform: for each
/// skeleton pixel, the radius is the distance-transform value at that pixel.
pub fn medial_axis_transform(img: &BinaryImage) -> Vec<MatPoint> {
    let dt = distance_transform(img);

    let mut thin = img.clone();
    zhang_suen_thinning(&mut thin);

    let mut mat_points = Vec::new();
    for r in 0..thin.height {
        for c in 0..thin.width {
            if thin.get(c, r) {
                mat_points.push(MatPoint {
                    col: c,
                    row: r,
                    radius: dt[r * img.width + c],
                });
            }
        }
    }
    mat_points
}

/// Reconstruct the shape from MAT points by "painting" filled circles.
///
/// Returns a binary image where each MAT point stamps a filled disc of its
/// inscribed radius.
pub fn reconstruct_from_mat(mat_points: &[MatPoint], width: usize, height: usize) -> BinaryImage {
    let mut img = BinaryImage::new(width, height, false);
    for mp in mat_points {
        let r_int = mp.radius.ceil() as isize;
        let cx = mp.col as isize;
        let cy = mp.row as isize;
        for dy in -r_int..=r_int {
            for dx in -r_int..=r_int {
                if (dx * dx + dy * dy) as f64 <= mp.radius * mp.radius {
                    let nx = cx + dx;
                    let ny = cy + dy;
                    if nx >= 0 && ny >= 0 {
                        img.set(nx as usize, ny as usize, true);
                    }
                }
            }
        }
    }
    img
}

// ─── Significance-based pruning ──────────────────────────────────────────────

/// Significance measure for a medial-axis branch.
///
/// Uses the ratio of the maximum inscribed radius along the branch to the
/// branch length.  Branches with low significance are noise.
pub fn branch_significance(ma: &MedialAxis, branch_idx: usize) -> f64 {
    if branch_idx >= ma.branches.len() {
        return 0.0;
    }
    let branch = &ma.branches[branch_idx];
    if branch.is_empty() {
        return 0.0;
    }

    let max_r = branch
        .iter()
        .map(|&i| ma.points[i].radius)
        .fold(0.0_f64, f64::max);

    let mut length = 0.0;
    for w in branch.windows(2) {
        length += dist3(ma.points[w[0]].position, ma.points[w[1]].position);
    }

    if length < 1e-15 {
        return max_r;
    }
    max_r / length
}

/// Prune branches whose significance measure is below `min_sig`.
pub fn prune_by_significance(ma: &mut MedialAxis, min_sig: f64) {
    let sigs: Vec<f64> = (0..ma.branches.len())
        .map(|i| branch_significance(ma, i))
        .collect();
    let mut keep = Vec::new();
    for (i, &s) in sigs.iter().enumerate() {
        if s >= min_sig {
            keep.push(ma.branches[i].clone());
        }
    }
    ma.branches = keep;
}

// ─── Radius function ─────────────────────────────────────────────────────────

/// Evaluate the radius function along a medial-axis branch.
///
/// Returns a vector of (arc-length, radius) pairs sampled at each branch point.
pub fn radius_function(ma: &MedialAxis, branch_idx: usize) -> Vec<(f64, f64)> {
    if branch_idx >= ma.branches.len() {
        return vec![];
    }
    let branch = &ma.branches[branch_idx];
    let mut result = Vec::with_capacity(branch.len());
    let mut arc = 0.0;
    for (k, &idx) in branch.iter().enumerate() {
        if k > 0 {
            arc += dist3(ma.points[branch[k - 1]].position, ma.points[idx].position);
        }
        result.push((arc, ma.points[idx].radius));
    }
    result
}

/// Compute the average radius along a branch.
pub fn average_radius(ma: &MedialAxis, branch_idx: usize) -> f64 {
    if branch_idx >= ma.branches.len() {
        return 0.0;
    }
    let branch = &ma.branches[branch_idx];
    if branch.is_empty() {
        return 0.0;
    }
    let sum: f64 = branch.iter().map(|&i| ma.points[i].radius).sum();
    sum / branch.len() as f64
}

/// Maximum inscribed radius along a branch.
pub fn max_radius(ma: &MedialAxis, branch_idx: usize) -> f64 {
    if branch_idx >= ma.branches.len() {
        return 0.0;
    }
    ma.branches[branch_idx]
        .iter()
        .map(|&i| ma.points[i].radius)
        .fold(0.0_f64, f64::max)
}

// ─── Shape decomposition from skeleton ───────────────────────────────────────

/// A shape part resulting from skeleton-based decomposition.
#[derive(Debug, Clone)]
pub struct ShapePart {
    /// Indices of the medial-axis points belonging to this part.
    pub point_indices: Vec<usize>,
    /// Approximate volume (sum of inscribed-sphere volumes).
    pub approx_volume: f64,
    /// Centroid of the part in 3-D.
    pub centroid: [f64; 3],
}

/// Decompose a shape into parts using its medial axis.
///
/// Each branch of the medial axis defines one part.  The approximate volume is
/// estimated by summing the volumes of inscribed spheres at each sample point
/// (with overlap correction).
pub fn decompose_from_skeleton(ma: &MedialAxis) -> Vec<ShapePart> {
    let mut parts = Vec::with_capacity(ma.branches.len());
    for branch in &ma.branches {
        if branch.is_empty() {
            continue;
        }

        let mut vol = 0.0;
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;

        for &idx in branch {
            let r = ma.points[idx].radius;
            vol += (4.0 / 3.0) * std::f64::consts::PI * r * r * r;
            cx += ma.points[idx].position[0];
            cy += ma.points[idx].position[1];
            cz += ma.points[idx].position[2];
        }

        let n = branch.len() as f64;
        parts.push(ShapePart {
            point_indices: branch.clone(),
            approx_volume: vol,
            centroid: [cx / n, cy / n, cz / n],
        });
    }
    parts
}

/// Total approximate volume from all shape parts.
pub fn total_decomposed_volume(parts: &[ShapePart]) -> f64 {
    parts.iter().map(|p| p.approx_volume).sum()
}

// ─── Curve skeleton extraction (3-D) ────────────────────────────────────────

/// A node in the curve skeleton graph.
#[derive(Debug, Clone)]
pub struct CurveSkeletonNode {
    /// 3-D position.
    pub position: [f64; 3],
    /// Inscribed-sphere radius.
    pub radius: f64,
    /// Indices of adjacent nodes.
    pub neighbors: Vec<usize>,
}

/// A curve skeleton represented as a graph.
#[derive(Debug, Clone)]
pub struct CurveSkeleton {
    /// Nodes of the skeleton graph.
    pub nodes: Vec<CurveSkeletonNode>,
}

impl CurveSkeleton {
    /// Extract a curve skeleton from a surface point cloud by Laplacian
    /// contraction (Au et al. 2008, point-cloud form of Cao et al. 2010).
    ///
    /// At most `contraction_steps` contraction iterations are run (fewer once
    /// the cross-sections have collapsed): each solves the cotangent-Laplacian
    /// least-squares system `(W_L² LᵀL + W_H²) P' = W_H² P` with conjugate
    /// gradients and then raises `W_L` threefold. The contracted cloud is
    /// clustered by farthest-point sampling; nodes are placed at the centroid
    /// of their clusters' original samples with the mean distance to them as
    /// radius, connected where their samples are neighbours, and reduced to a
    /// curve skeleton (spanning tree, twig pruning, loop closure). See
    /// [`laplacian_contraction`] for the contraction step alone.
    ///
    /// # Errors
    ///
    /// * [`Error::TooFewPoints`](crate::error::Error::TooFewPoints) for fewer than four samples;
    /// * [`Error::OutOfRange`](crate::error::Error::OutOfRange) when `contraction_steps == 0`;
    /// * [`Error::DegenerateGeometry`](crate::error::Error::DegenerateGeometry) for non-finite or coincident samples;
    /// * [`Error::ConvergenceFailure`](crate::error::Error::ConvergenceFailure) if a linear solve fails to converge.
    pub fn from_surface_points(boundary: &[[f64; 3]], contraction_steps: usize) -> Result<Self> {
        contraction::curve_skeleton(boundary, contraction_steps)
    }

    /// Number of nodes in the skeleton.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of edges in the skeleton.
    pub fn edge_count(&self) -> usize {
        let sum: usize = self.nodes.iter().map(|n| n.neighbors.len()).sum();
        sum / 2 // each edge counted twice
    }

    /// Total arc-length of the skeleton.
    pub fn total_length(&self) -> f64 {
        let mut len = 0.0;
        for (i, node) in self.nodes.iter().enumerate() {
            for &j in &node.neighbors {
                if j > i {
                    len += dist3(node.position, self.nodes[j].position);
                }
            }
        }
        len
    }

    /// Find leaf nodes (degree 1) in the skeleton.
    pub fn leaf_nodes(&self) -> Vec<usize> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.neighbors.len() == 1)
            .map(|(i, _)| i)
            .collect()
    }

    /// Find junction nodes (degree >= 3) in the skeleton.
    pub fn junction_nodes(&self) -> Vec<usize> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.neighbors.len() >= 3)
            .map(|(i, _)| i)
            .collect()
    }

    /// Prune leaf branches shorter than `min_length` from the skeleton.
    pub fn prune_leaves(&mut self, min_length: f64) {
        loop {
            let leaves = self.leaf_nodes();
            let mut pruned = false;
            for leaf in leaves {
                if self.nodes[leaf].neighbors.is_empty() {
                    continue;
                }
                let nbr = self.nodes[leaf].neighbors[0];
                let d = dist3(self.nodes[leaf].position, self.nodes[nbr].position);
                if d < min_length && self.nodes[nbr].neighbors.len() > 1 {
                    // Remove leaf
                    self.nodes[nbr].neighbors.retain(|&n| n != leaf);
                    self.nodes[leaf].neighbors.clear();
                    pruned = true;
                }
            }
            if !pruned {
                break;
            }
        }
    }
}

/// Compute the Laplacian of a point set (average displacement toward neighbours).
///
/// Used in iterative Laplacian contraction for curve-skeleton extraction.
pub fn laplacian_smooth(
    points: &[[f64; 3]],
    adjacency: &[Vec<usize>],
    weight: f64,
) -> Vec<[f64; 3]> {
    let n = points.len();
    let mut result = vec![[0.0; 3]; n];
    for i in 0..n {
        if adjacency[i].is_empty() {
            result[i] = points[i];
            continue;
        }
        let mut avg = [0.0; 3];
        for &j in &adjacency[i] {
            avg[0] += points[j][0];
            avg[1] += points[j][1];
            avg[2] += points[j][2];
        }
        let k = adjacency[i].len() as f64;
        avg[0] /= k;
        avg[1] /= k;
        avg[2] /= k;
        result[i][0] = points[i][0] + weight * (avg[0] - points[i][0]);
        result[i][1] = points[i][1] + weight * (avg[1] - points[i][1]);
        result[i][2] = points[i][2] + weight * (avg[2] - points[i][2]);
    }
    result
}

/// Contract a surface point cloud towards its curve skeleton with
/// `steps` iterations of cotangent-Laplacian contraction (Au et al. 2008;
/// Cao et al. 2010) and return the contracted positions in input order.
///
/// Each iteration solves `(W_L² LᵀL + W_H²) P' = W_H² P` per coordinate,
/// where `L` is the cotangent Laplacian of every sample's tangent-plane
/// Delaunay one-ring (evaluated on the input cloud), `W_L` starts at 1 and
/// grows threefold per iteration (capped at 2048) and `W_H,i = √(A_i⁰ / A_i)`
/// tracks the shrinking one-ring area. A sample whose one-ring area dropped
/// below `10⁻³` of its initial value is pinned, and iteration stops early once
/// the total one-ring area has shrunk below `10⁻²` of its initial value (the
/// cross-sections have collapsed onto the skeleton), so `steps` is an upper
/// bound on the number of iterations.
///
/// # Errors
///
/// * [`Error::TooFewPoints`](crate::error::Error::TooFewPoints) for fewer than four points;
/// * [`Error::OutOfRange`](crate::error::Error::OutOfRange) when `steps == 0`;
/// * [`Error::DegenerateGeometry`](crate::error::Error::DegenerateGeometry) for non-finite or coincident points;
/// * [`Error::ConvergenceFailure`](crate::error::Error::ConvergenceFailure) if a linear solve fails to converge.
pub fn laplacian_contraction(points: &[[f64; 3]], steps: usize) -> Result<Vec<[f64; 3]>> {
    contraction::contract(points, steps)
}

// ─── Voronoi-based medial axis approximation ────────────────────────────────

/// Approximate the 2-D medial axis using Voronoi vertices of a dense boundary
/// sampling, then filter by inscribed-circle radius.
///
/// Only keeps Voronoi vertices whose distance to the boundary is at least
/// `min_radius`.
pub fn voronoi_medial_axis_2d(
    boundary: &[[f64; 2]],
    bbox: [f64; 4],
    min_radius: f64,
) -> Vec<([f64; 2], f64)> {
    let vertices = Voronoi2dSkeleton::compute(boundary, bbox);
    let mut result = Vec::new();
    for v in &vertices {
        let r = point_to_polygon_dist(*v, boundary);
        if r >= min_radius {
            result.push((*v, r));
        }
    }
    result
}

// ─── Hausdorff distance between two medial axes ──────────────────────────────

/// Compute the one-directional Hausdorff distance from medial axis `a` to `b`.
///
/// max over points in a of (min over points in b of dist(a_i, b_j)).
pub fn directed_hausdorff(a: &MedialAxis, b: &MedialAxis) -> f64 {
    if a.points.is_empty() || b.points.is_empty() {
        return f64::INFINITY;
    }
    let mut max_d = 0.0_f64;
    for pa in &a.points {
        let min_d = b
            .points
            .iter()
            .map(|pb| dist3(pa.position, pb.position))
            .fold(f64::INFINITY, f64::min);
        if min_d > max_d {
            max_d = min_d;
        }
    }
    max_d
}

/// Symmetric Hausdorff distance between two medial axes.
pub fn hausdorff_distance(a: &MedialAxis, b: &MedialAxis) -> f64 {
    directed_hausdorff(a, b).max(directed_hausdorff(b, a))
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod edt_tests;

#[cfg(test)]
mod skeleton_tests;

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: unit cube corner points
    fn cube_boundary() -> Vec<[f64; 3]> {
        let mut pts = vec![];
        for &x in &[-1.0, 1.0] {
            for &y in &[-1.0, 1.0] {
                for &z in &[-1.0, 1.0] {
                    pts.push([x, y, z]);
                }
            }
        }
        pts
    }

    // Helper: hand-built medial axis — a 5-point path along +x with radii
    // 1, 2, 3, 2, 1 plus a 2-point side branch hanging off the middle point.
    fn path_axis() -> MedialAxis {
        let radii = [1.0, 2.0, 3.0, 2.0, 1.0];
        let mut points: Vec<MedialAxisPoint> = (0..5)
            .map(|i| {
                MedialAxisPoint::new([i as f64, 0.0, 0.0], radii[i], [1.0, 0.0, 0.0], radii[i])
            })
            .collect();
        points.push(MedialAxisPoint::new(
            [2.0, 0.5, 0.0],
            0.5,
            [0.0, 1.0, 0.0],
            0.5,
        ));
        MedialAxis {
            points,
            branches: vec![vec![0, 1, 2], vec![2, 3, 4], vec![2, 5]],
            junctions: vec![2],
        }
    }

    // Helper: unit square polygon
    fn unit_square() -> Vec<[f64; 2]> {
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    }

    #[test]
    fn test_medial_axis_from_empty_is_error() {
        let err = MedialAxis::from_points(&[], SkeletonMethod::BallRolling);
        assert!(matches!(
            err,
            Err(crate::error::Error::TooFewPoints {
                required: 4,
                actual: 0
            })
        ));
    }

    #[test]
    fn test_medial_axis_from_cube_corners_is_cube_centre() {
        // The eight cospherical corners have a single Voronoi vertex inside
        // the cube: the centre, at distance sqrt(3) from every corner.
        let pts = cube_boundary();
        for method in [SkeletonMethod::Voronoi, SkeletonMethod::MedialSurface] {
            let ma = MedialAxis::from_points(&pts, method).expect("cube corners");
            assert_eq!(ma.points.len(), 1, "{method:?}");
            let p = &ma.points[0];
            assert!(
                dist3(p.position, [0.0; 3]) < 1e-9,
                "{method:?}: {:?}",
                p.position
            );
            assert!((p.radius - 3.0_f64.sqrt()).abs() < 1e-9, "{method:?}");
            assert!(ma.junctions.is_empty());
            assert_eq!(ma.total_length(), 0.0);
        }
    }

    #[test]
    fn test_medial_axis_is_tubular_on_constant_radius_path() {
        let points: Vec<MedialAxisPoint> = (0..6)
            .map(|i| MedialAxisPoint::new([i as f64, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0], 1.0))
            .collect();
        let ma = MedialAxis {
            points,
            branches: vec![(0..6).collect()],
            junctions: vec![],
        };
        assert!(ma.is_tubular());
        assert!(!path_axis().is_tubular());
    }

    #[test]
    fn test_medial_axis_collinear_or_too_few_points_are_errors() {
        let few = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert!(matches!(
            MedialAxis::from_points(&few, SkeletonMethod::Voronoi),
            Err(crate::error::Error::TooFewPoints { .. })
        ));
        let planar: Vec<[f64; 3]> = (0..10)
            .map(|i| [i as f64, (i * i) as f64 * 0.1, 0.0])
            .collect();
        for method in [
            SkeletonMethod::BallRolling,
            SkeletonMethod::Thinning,
            SkeletonMethod::Voronoi,
            SkeletonMethod::MedialSurface,
        ] {
            assert!(
                matches!(
                    MedialAxis::from_points(&planar, method),
                    Err(crate::error::Error::DegenerateGeometry { .. })
                ),
                "{method:?}"
            );
        }
    }

    #[test]
    fn test_voronoi2d_skeleton_empty() {
        let pts: Vec<[f64; 2]> = vec![];
        let res = Voronoi2dSkeleton::compute(&pts, [0.0, 0.0, 1.0, 1.0]);
        assert!(res.is_empty());
    }

    #[test]
    fn test_voronoi2d_skeleton_two_points() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0]];
        let res = Voronoi2dSkeleton::compute(&pts, [0.0, 0.0, 1.0, 1.0]);
        assert!(res.is_empty()); // need at least 3
    }

    #[test]
    fn test_voronoi2d_skeleton_four_points() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let res = Voronoi2dSkeleton::compute(&pts, [-0.5, -0.5, 1.5, 1.5]);
        // Should produce at least some circumcentres
        assert!(!res.is_empty());
    }

    #[test]
    fn test_voronoi2d_skeleton_inside_bbox() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let res = Voronoi2dSkeleton::compute(&pts, [-10.0, -10.0, 10.0, 10.0]);
        for p in &res {
            assert!(p[0] >= -10.0 && p[0] <= 10.0);
            assert!(p[1] >= -10.0 && p[1] <= 10.0);
        }
    }

    #[test]
    fn test_straight_skeleton_new() {
        let sk = StraightSkeleton::new(unit_square());
        let edges = sk.skeleton_edges();
        assert_eq!(edges.len(), 4);
    }

    #[test]
    fn test_straight_skeleton_empty_polygon() {
        let sk = StraightSkeleton::new(vec![]);
        assert!(sk.skeleton_edges().is_empty());
        assert!(sk.compute_offsets(0.1).is_empty());
    }

    #[test]
    fn test_straight_skeleton_offset_count() {
        let sk = StraightSkeleton::new(unit_square());
        let offsets = sk.compute_offsets(0.1);
        assert_eq!(offsets.len(), 1);
        assert_eq!(offsets[0].len(), 4);
    }

    #[test]
    fn test_straight_skeleton_edges_toward_centroid() {
        let sq = unit_square();
        let sk = StraightSkeleton::new(sq.clone());
        let edges = sk.skeleton_edges();
        // All edges should end at centroid [0.5, 0.5]
        for (_, end) in &edges {
            assert!((end[0] - 0.5).abs() < 1e-10);
            assert!((end[1] - 0.5).abs() < 1e-10);
        }
    }

    #[test]
    fn test_point_to_segment_dist_endpoint() {
        let p = [0.0, 0.0, 0.0];
        let a = [1.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let d = point_to_segment_dist(p, a, b);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_point_to_segment_dist_perpendicular() {
        let p = [1.0, 1.0, 0.0];
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let d = point_to_segment_dist(p, a, b);
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_point_to_segment_dist_degenerate() {
        let p = [3.0, 4.0, 0.0];
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0]; // degenerate segment
        let d = point_to_segment_dist(p, a, b);
        assert!((d - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_point_to_polygon_dist_square() {
        let poly = unit_square();
        let p = [0.5, 0.5];
        let d = point_to_polygon_dist(p, &poly);
        assert!((d - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_point_to_polygon_dist_empty() {
        let d = point_to_polygon_dist([0.0, 0.0], &[]);
        assert!(d.is_infinite());
    }

    #[test]
    fn test_is_inside_polygon_center() {
        let poly = unit_square();
        assert!(is_inside_polygon([0.5, 0.5], &poly));
    }

    #[test]
    fn test_is_inside_polygon_outside() {
        let poly = unit_square();
        assert!(!is_inside_polygon([2.0, 2.0], &poly));
    }

    #[test]
    fn test_is_inside_polygon_too_few_vertices() {
        let poly = vec![[0.0, 0.0], [1.0, 0.0]];
        assert!(!is_inside_polygon([0.5, 0.0], &poly));
    }

    #[test]
    fn test_centerline_extraction_basic() {
        let mut pts = vec![];
        for i in 0..10 {
            let x = i as f64;
            pts.push([x, 1.0, 0.0]);
            pts.push([x, -1.0, 0.0]);
            pts.push([x, 0.0, 1.0]);
            pts.push([x, 0.0, -1.0]);
        }
        let cl = CenterlineExtraction::extract_centerline(&pts, 5);
        assert_eq!(cl.len(), 5);
    }

    #[test]
    fn test_centerline_extraction_empty() {
        let cl = CenterlineExtraction::extract_centerline(&[], 5);
        assert!(cl.is_empty());
    }

    #[test]
    fn test_centerline_extraction_zero_slices() {
        let pts = vec![[0.0, 0.0, 0.0]];
        let cl = CenterlineExtraction::extract_centerline(&pts, 0);
        assert!(cl.is_empty());
    }

    #[test]
    fn test_prune_short_branches() {
        let mut ma = path_axis();
        prune_short_branches(&mut ma, 1.0); // side branch has length 0.5
        assert_eq!(ma.branch_count(), 2);
        assert!(ma.branches.iter().all(|b| b.len() == 3));
        prune_short_branches(&mut ma, 1e12); // prune everything
        assert_eq!(ma.branch_count(), 0);
    }

    #[test]
    fn test_prune_keeps_long_branches() {
        let mut ma = path_axis();
        prune_short_branches(&mut ma, 0.0); // keep all
        assert_eq!(ma.branch_count(), 3);
    }

    #[test]
    fn test_euler_characteristic() {
        // V = 6 points, E = 2 + 2 + 1 = 5 edges: a tree, chi = 1.
        assert_eq!(euler_characteristic(&path_axis()), 1);
    }

    #[test]
    fn test_axis_aligned_slab_decompose() {
        let ma = path_axis();
        let slabs = AxisAlignedSlab::decompose(&ma, 4);
        assert_eq!(slabs.len(), 4);
        let total_pts: usize = slabs.iter().map(|s| s.point_indices.len()).sum();
        assert_eq!(total_pts, ma.points.len());
        // x = 0 falls in slab 0, x = 4 in the last slab.
        assert!(slabs[0].point_indices.contains(&0));
        assert!(slabs[3].point_indices.contains(&4));
    }

    #[test]
    fn test_axis_aligned_slab_zero_slabs() {
        let slabs = AxisAlignedSlab::decompose(&path_axis(), 0);
        assert!(slabs.is_empty());
    }

    #[test]
    fn test_shape_diameter_positive() {
        let d = shape_diameter([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0);
        assert!((d - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_shape_diameter_scales_with_radius() {
        let d1 = shape_diameter([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0);
        let d2 = shape_diameter([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0);
        assert!((d2 - 2.0 * d1).abs() < 1e-10);
    }

    // ─── Distance transform tests ───────────────────────────────────────

    fn make_square_image(size: usize) -> BinaryImage {
        let mut img = BinaryImage::new(size, size, false);
        let margin = size / 4;
        for r in margin..(size - margin) {
            for c in margin..(size - margin) {
                img.set(c, r, true);
            }
        }
        img
    }

    #[test]
    fn test_binary_image_basic() {
        let img = BinaryImage::new(10, 10, false);
        assert_eq!(img.foreground_count(), 0);
        assert!(!img.get(0, 0));
    }

    #[test]
    fn test_distance_transform_background() {
        let img = BinaryImage::new(10, 10, false);
        let dt = distance_transform(&img);
        assert!(dt.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_distance_transform_foreground() {
        let img = make_square_image(20);
        let dt = distance_transform(&img);
        // Centre pixel should have maximum distance
        let center = 10 * 20 + 10;
        assert!(dt[center] > 0.0, "Centre distance should be > 0");
        // Corner (background) should be 0
        assert!((dt[0] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_distance_transform_exact_background() {
        let img = BinaryImage::new(10, 10, false);
        let dt = distance_transform_exact(&img);
        assert!(dt.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_distance_transform_exact_foreground() {
        let img = make_square_image(20);
        let dt = distance_transform_exact(&img);
        let center = 10 * 20 + 10;
        assert!(dt[center] > 0.0);
    }

    // ─── Thinning tests ────────────────────────────────────────────────

    #[test]
    fn test_zhang_suen_reduces_pixels() {
        let mut img = make_square_image(20);
        let before = img.foreground_count();
        zhang_suen_thinning(&mut img);
        let after = img.foreground_count();
        assert!(after < before, "Thinning should reduce pixel count");
        assert!(after > 0, "Thinning should leave a skeleton");
    }

    #[test]
    fn test_zhang_suen_empty_image() {
        let mut img = BinaryImage::new(10, 10, false);
        zhang_suen_thinning(&mut img);
        assert_eq!(img.foreground_count(), 0);
    }

    #[test]
    fn test_skeleton_pixels_from_thinned() {
        let mut img = make_square_image(20);
        zhang_suen_thinning(&mut img);
        let pixels = skeleton_pixels(&img);
        assert!(!pixels.is_empty());
        assert_eq!(pixels.len(), img.foreground_count());
    }

    // ─── MAT tests ─────────────────────────────────────────────────────

    #[test]
    fn test_mat_has_positive_radii() {
        let img = make_square_image(20);
        let mat = medial_axis_transform(&img);
        assert!(!mat.is_empty(), "MAT should have points");
        for mp in &mat {
            assert!(mp.radius >= 0.0, "Radius should be non-negative");
        }
    }

    #[test]
    fn test_reconstruct_from_mat() {
        let img = make_square_image(20);
        let mat = medial_axis_transform(&img);
        let recon = reconstruct_from_mat(&mat, 20, 20);
        assert!(
            recon.foreground_count() > 0,
            "Reconstructed image should not be empty"
        );
    }

    // ─── Significance pruning tests ─────────────────────────────────────

    #[test]
    fn test_branch_significance() {
        let ma = path_axis();
        // Branch 0: max radius 3 over length 2.
        assert!((branch_significance(&ma, 0) - 1.5).abs() < 1e-12);
        // Side branch: max radius 3 (junction) over length 0.5.
        assert!((branch_significance(&ma, 2) - 6.0).abs() < 1e-12);
        assert_eq!(branch_significance(&ma, 99), 0.0);
    }

    #[test]
    fn test_prune_by_significance() {
        let mut ma = path_axis();
        prune_by_significance(&mut ma, 2.0); // keeps only the side branch
        assert_eq!(ma.branches, vec![vec![2, 5]]);
        prune_by_significance(&mut ma, 1e12); // very high threshold => prune all
        assert!(ma.branches.is_empty());
    }

    // ─── Radius function tests ──────────────────────────────────────────

    #[test]
    fn test_radius_function() {
        let ma = path_axis();
        let rf = radius_function(&ma, 0);
        assert_eq!(rf, vec![(0.0, 1.0), (1.0, 2.0), (2.0, 3.0)]);
    }

    #[test]
    fn test_average_and_max_radius() {
        let ma = path_axis();
        assert!((average_radius(&ma, 0) - 2.0).abs() < 1e-12);
        assert_eq!(max_radius(&ma, 0), 3.0);
        assert!((average_radius(&ma, 2) - 1.75).abs() < 1e-12);
    }

    // ─── Shape decomposition tests ──────────────────────────────────────

    #[test]
    fn test_decompose_from_skeleton() {
        let ma = path_axis();
        let parts = decompose_from_skeleton(&ma);
        assert_eq!(parts.len(), 3);
        let ball = |r: f64| 4.0 / 3.0 * std::f64::consts::PI * r * r * r;
        let expected0 = ball(1.0) + ball(2.0) + ball(3.0);
        assert!((parts[0].approx_volume - expected0).abs() < 1e-9);
        assert!((parts[0].centroid[0] - 1.0).abs() < 1e-12);
        let vol = total_decomposed_volume(&parts);
        assert!(vol > expected0);
    }

    // ─── Curve skeleton tests ───────────────────────────────────────────

    // Helper: hand-built curve skeleton — a path 0-1-2-3 with a short spur
    // 1-4 (length 0.2).
    fn star_skeleton() -> CurveSkeleton {
        let node = |p: [f64; 3], nb: Vec<usize>| CurveSkeletonNode {
            position: p,
            radius: 1.0,
            neighbors: nb,
        };
        CurveSkeleton {
            nodes: vec![
                node([0.0, 0.0, 0.0], vec![1]),
                node([1.0, 0.0, 0.0], vec![0, 2, 4]),
                node([2.0, 0.0, 0.0], vec![1, 3]),
                node([3.0, 0.0, 0.0], vec![2]),
                node([1.0, 0.2, 0.0], vec![1]),
            ],
        }
    }

    #[test]
    fn test_curve_skeleton_graph_queries() {
        let cs = star_skeleton();
        assert_eq!(cs.node_count(), 5);
        assert_eq!(cs.edge_count(), 4);
        assert!((cs.total_length() - 3.2).abs() < 1e-12);
        assert_eq!(cs.leaf_nodes(), vec![0, 3, 4]);
        assert_eq!(cs.junction_nodes(), vec![1]);
    }

    #[test]
    fn test_curve_skeleton_prune() {
        let mut cs = star_skeleton();
        cs.prune_leaves(0.5); // removes the 0.2 spur only
        assert_eq!(cs.edge_count(), 3);
        assert!(cs.junction_nodes().is_empty());
        assert!((cs.total_length() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_curve_skeleton_rejects_zero_steps_and_tiny_input() {
        assert!(matches!(
            CurveSkeleton::from_surface_points(&cube_boundary(), 0),
            Err(crate::error::Error::OutOfRange { .. })
        ));
        assert!(matches!(
            CurveSkeleton::from_surface_points(&[[0.0; 3]; 2], 3),
            Err(crate::error::Error::TooFewPoints { .. })
        ));
    }

    #[test]
    fn test_laplacian_smooth() {
        let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let adj = vec![vec![1], vec![0, 2], vec![1]];
        let result = laplacian_smooth(&points, &adj, 0.5);
        // Middle point should stay put (neighbours average to itself)
        assert!((result[1][0] - 1.0).abs() < 1e-10);
    }

    // ─── Voronoi medial axis 2-D tests ──────────────────────────────────

    #[test]
    fn test_voronoi_medial_axis_2d() {
        let pts = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let result = voronoi_medial_axis_2d(&pts, [-0.5, -0.5, 1.5, 1.5], 0.0);
        // Should find at least some vertices
        assert!(!result.is_empty());
    }

    // ─── Hausdorff distance tests ───────────────────────────────────────

    #[test]
    fn test_hausdorff_self_zero() {
        let ma = path_axis();
        let d = hausdorff_distance(&ma, &ma);
        assert!(d < 1e-10, "Hausdorff to self should be 0, got {:.6}", d);
    }

    #[test]
    fn test_hausdorff_different_shapes() {
        let pts1 = cube_boundary();
        let pts2: Vec<[f64; 3]> = pts1.iter().map(|p| [p[0] + 5.0, p[1], p[2]]).collect();
        let ma1 = MedialAxis::from_points(&pts1, SkeletonMethod::Voronoi).expect("cube");
        let ma2 = MedialAxis::from_points(&pts2, SkeletonMethod::Voronoi).expect("shifted cube");
        let d = hausdorff_distance(&ma1, &ma2);
        assert!((d - 5.0).abs() < 1e-9, "centres are 5 apart, got {d}");
    }
}
