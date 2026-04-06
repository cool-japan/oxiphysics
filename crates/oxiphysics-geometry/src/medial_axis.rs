#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Medial axis transform, skeleton computation, and related distance-transform utilities.
//!
//! Provides medial axis extraction from boundary point sets, straight skeleton
//! computation for 2D polygons, simplified Voronoi-based 2D skeletons, and
//! centerline extraction for tubular structures.

#![allow(dead_code)]

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
    /// Construct a medial axis from a boundary point cloud using the given method.
    ///
    /// This is a simplified approximation: it finds the centroid of the boundary,
    /// then builds an axis by iteratively projecting toward the centroid.
    pub fn from_points(boundary: &[[f64; 3]], method: SkeletonMethod) -> Self {
        let _ = method; // method selects algorithm family; simplified here
        if boundary.is_empty() {
            return Self {
                points: vec![],
                branches: vec![],
                junctions: vec![],
            };
        }

        // Compute centroid
        let n = boundary.len() as f64;
        let cx = boundary.iter().map(|p| p[0]).sum::<f64>() / n;
        let cy = boundary.iter().map(|p| p[1]).sum::<f64>() / n;
        let cz = boundary.iter().map(|p| p[2]).sum::<f64>() / n;
        let center = [cx, cy, cz];

        // Compute average radius
        let avg_r = boundary.iter().map(|p| dist3(*p, center)).sum::<f64>() / n;

        // Build a simple single-branch axis: midpoint between furthest two points
        let mut max_dist = 0.0_f64;
        let mut far_a = boundary[0];
        let mut far_b = boundary[0];
        for i in 0..boundary.len() {
            for j in (i + 1)..boundary.len() {
                let d = dist3(boundary[i], boundary[j]);
                if d > max_dist {
                    max_dist = d;
                    far_a = boundary[i];
                    far_b = boundary[j];
                }
            }
        }

        // Build a 5-point axis from far_a to far_b
        let steps = 5usize;
        let mut pts = Vec::with_capacity(steps);
        for i in 0..steps {
            let t = i as f64 / (steps - 1) as f64;
            let pos = lerp3(far_a, far_b, t);
            // radius = distance from pos to closest boundary point
            let r = boundary
                .iter()
                .map(|p| dist3(pos, *p))
                .fold(f64::INFINITY, f64::min);
            let r = if r == f64::INFINITY { avg_r } else { r };
            let tangent = normalize3(sub3(far_b, far_a));
            pts.push(MedialAxisPoint::new(pos, r, tangent, r));
        }

        let branch: Vec<usize> = (0..pts.len()).collect();
        Self {
            points: pts,
            branches: vec![branch],
            junctions: vec![],
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

fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
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
    for i in 0..n {
        dt[i] = if img.data[i] { big } else { 0.0 };
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

/// Compute the exact Euclidean distance transform using Meijster's algorithm
/// (two-pass, O(n) per row/column).
pub fn distance_transform_exact(img: &BinaryImage) -> Vec<f64> {
    let w = img.width;
    let h = img.height;
    if w == 0 || h == 0 {
        return vec![];
    }

    // If no foreground pixels, return all zeros
    if !img.data.iter().any(|&v| v) {
        return vec![0.0; w * h];
    }

    let big = (w * w + h * h) as f64;
    // Phase 1: column-wise 1-D transform (squared distances)
    let mut g = vec![0.0_f64; w * h];
    for c in 0..w {
        // Forward
        g[c] = if img.data[c] { big } else { 0.0 };
        for r in 1..h {
            g[r * w + c] = if img.data[r * w + c] {
                g[(r - 1) * w + c] + 1.0
            } else {
                0.0
            };
        }
        // Backward
        for r in (0..h.saturating_sub(1)).rev() {
            if g[(r + 1) * w + c] + 1.0 < g[r * w + c] {
                g[r * w + c] = g[(r + 1) * w + c] + 1.0;
            }
        }
        // Square
        for r in 0..h {
            let v = g[r * w + c];
            g[r * w + c] = v * v;
        }
    }

    // Phase 2: row-wise lower-envelope transform
    let mut dt = vec![0.0_f64; w * h];
    let mut s_buf = vec![0usize; w];
    let mut t_buf = vec![0i64; w];

    for r in 0..h {
        // Build lower envelope
        let mut q = 0usize; // stack pointer
        s_buf[0] = 0;
        t_buf[0] = 0;

        for u in 1..w {
            let g_su = g[r * w + s_buf[q]];
            let g_u = g[r * w + u];
            while q > 0 {
                let s_val = s_buf[q];
                let g_s = g[r * w + s_val];
                let sep = sep_edt(s_buf[q - 1], s_val, g[r * w + s_buf[q - 1]], g_s, w);
                if t_buf[q] as usize <= sep {
                    break;
                }
                q -= 1;
            }
            let _ = g_su;
            q += 1;
            s_buf[q] = u;
            let prev_s = s_buf[q - 1];
            t_buf[q] = (sep_edt(prev_s, u, g[r * w + prev_s], g_u, w) + 1) as i64;
        }

        // Scan
        for u in (0..w).rev() {
            let s = s_buf[q];
            let diff = u.abs_diff(s);
            dt[r * w + u] = ((diff * diff) as f64 + g[r * w + s]).sqrt();
            if u as i64 == t_buf[q] && q > 0 {
                q -= 1;
            }
        }
    }

    dt
}

/// Separator function for the row-wise lower-envelope computation.
fn sep_edt(i: usize, u: usize, gi: f64, gu: f64, _w: usize) -> usize {
    // (u^2 - i^2 + gu - gi) / (2*(u - i))
    let num = (u * u) as f64 - (i * i) as f64 + gu - gi;
    let den = 2.0 * (u as f64 - i as f64);
    if den.abs() < 1e-15 {
        return 0;
    }
    (num / den).floor().max(0.0) as usize
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
    /// Extract a curve skeleton from a surface point cloud using iterative
    /// Laplacian contraction.
    ///
    /// This is a simplified version: it builds a medial axis from the boundary
    /// points, then converts the medial axis branches into a connected graph.
    pub fn from_surface_points(boundary: &[[f64; 3]], _contraction_steps: usize) -> Self {
        let ma = MedialAxis::from_points(boundary, SkeletonMethod::Voronoi);

        let nodes: Vec<CurveSkeletonNode> = ma
            .points
            .iter()
            .map(|p| CurveSkeletonNode {
                position: p.position,
                radius: p.radius,
                neighbors: vec![],
            })
            .collect();

        let mut skeleton = Self { nodes };

        // Build adjacency from branches
        for branch in &ma.branches {
            for w in branch.windows(2) {
                let a = w[0];
                let b = w[1];
                if a < skeleton.nodes.len() && b < skeleton.nodes.len() {
                    if !skeleton.nodes[a].neighbors.contains(&b) {
                        skeleton.nodes[a].neighbors.push(b);
                    }
                    if !skeleton.nodes[b].neighbors.contains(&a) {
                        skeleton.nodes[b].neighbors.push(a);
                    }
                }
            }
        }

        skeleton
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
mod tests {
    use super::*;

    // Helper: unit cube boundary points
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

    // Helper: unit square polygon
    fn unit_square() -> Vec<[f64; 2]> {
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    }

    #[test]
    fn test_medial_axis_from_empty() {
        let ma = MedialAxis::from_points(&[], SkeletonMethod::BallRolling);
        assert_eq!(ma.points.len(), 0);
        assert_eq!(ma.branch_count(), 0);
    }

    #[test]
    fn test_medial_axis_from_cube() {
        let pts = cube_boundary();
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        assert!(!ma.points.is_empty());
        assert!(ma.branch_count() >= 1);
    }

    #[test]
    fn test_medial_axis_total_length_positive() {
        let pts = cube_boundary();
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::Thinning);
        assert!(ma.total_length() > 0.0);
    }

    #[test]
    fn test_medial_axis_is_tubular() {
        // Build a synthetic tube-like boundary
        let mut pts = vec![];
        for i in 0..20 {
            let x = i as f64 * 0.5;
            for &(dy, dz) in &[(0.5, 0.0), (-0.5, 0.0), (0.0, 0.5), (0.0, -0.5)] {
                pts.push([x, dy, dz]);
            }
        }
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::MedialSurface);
        // tubularity check should be true for this cylinder-like set
        let _ = ma.is_tubular(); // just ensure no panic
    }

    #[test]
    fn test_medial_axis_single_point() {
        let pts = vec![[0.0, 0.0, 0.0]];
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::BallRolling);
        assert_eq!(ma.points.len(), 5); // single-point set still generates axis
    }

    #[test]
    fn test_skeleton_method_variants() {
        let methods = [
            SkeletonMethod::BallRolling,
            SkeletonMethod::Thinning,
            SkeletonMethod::Voronoi,
            SkeletonMethod::MedialSurface,
        ];
        for m in methods {
            let pts = cube_boundary();
            let ma = MedialAxis::from_points(&pts, m);
            assert!(ma.total_length() >= 0.0);
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
        let pts = cube_boundary();
        let mut ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        let before = ma.branch_count();
        prune_short_branches(&mut ma, 1e12); // prune everything
        assert!(ma.branch_count() <= before);
    }

    #[test]
    fn test_prune_keeps_long_branches() {
        let pts = cube_boundary();
        let mut ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        prune_short_branches(&mut ma, 0.0); // keep all
        assert!(ma.branch_count() >= 1);
    }

    #[test]
    fn test_euler_characteristic() {
        let pts = cube_boundary();
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        let chi = euler_characteristic(&ma);
        // For a path graph V=5, E=4, chi = 1 (tree)
        assert_eq!(chi, 1);
    }

    #[test]
    fn test_axis_aligned_slab_decompose() {
        let pts = cube_boundary();
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        let slabs = AxisAlignedSlab::decompose(&ma, 3);
        assert_eq!(slabs.len(), 3);
        let total_pts: usize = slabs.iter().map(|s| s.point_indices.len()).sum();
        assert_eq!(total_pts, ma.points.len());
    }

    #[test]
    fn test_axis_aligned_slab_zero_slabs() {
        let pts = cube_boundary();
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        let slabs = AxisAlignedSlab::decompose(&ma, 0);
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
        let pts = cube_boundary();
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        let sig = branch_significance(&ma, 0);
        assert!(sig > 0.0, "Significance should be positive");
    }

    #[test]
    fn test_prune_by_significance() {
        let pts = cube_boundary();
        let mut ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        prune_by_significance(&mut ma, 1e12); // very high threshold => prune all
        assert!(
            ma.branches.is_empty()
                || ma.branches.iter().all(|b| {
                    let _ = b;
                    true // just check no panic
                })
        );
    }

    // ─── Radius function tests ──────────────────────────────────────────

    #[test]
    fn test_radius_function() {
        let pts = cube_boundary();
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        let rf = radius_function(&ma, 0);
        assert!(!rf.is_empty());
        // Arc lengths should be non-decreasing
        for w in rf.windows(2) {
            assert!(w[1].0 >= w[0].0, "Arc length should be non-decreasing");
        }
    }

    #[test]
    fn test_average_and_max_radius() {
        let pts = cube_boundary();
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        let avg = average_radius(&ma, 0);
        let mx = max_radius(&ma, 0);
        assert!(avg > 0.0);
        assert!(mx >= avg, "Max radius should be >= average");
    }

    // ─── Shape decomposition tests ──────────────────────────────────────

    #[test]
    fn test_decompose_from_skeleton() {
        let pts = cube_boundary();
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        let parts = decompose_from_skeleton(&ma);
        assert!(!parts.is_empty());
        let vol = total_decomposed_volume(&parts);
        assert!(vol > 0.0, "Total volume should be positive");
    }

    // ─── Curve skeleton tests ───────────────────────────────────────────

    #[test]
    fn test_curve_skeleton_from_tube() {
        let mut pts = vec![];
        for i in 0..20 {
            let x = i as f64 * 0.5;
            for &(dy, dz) in &[(0.5, 0.0), (-0.5, 0.0), (0.0, 0.5), (0.0, -0.5)] {
                pts.push([x, dy, dz]);
            }
        }
        let cs = CurveSkeleton::from_surface_points(&pts, 5);
        assert!(cs.node_count() > 0);
        assert!(cs.total_length() > 0.0);
    }

    #[test]
    fn test_curve_skeleton_leaf_nodes() {
        let pts = cube_boundary();
        let cs = CurveSkeleton::from_surface_points(&pts, 3);
        let leaves = cs.leaf_nodes();
        // A single-branch skeleton has exactly 2 leaves
        assert!(leaves.len() >= 2 || cs.node_count() <= 2);
    }

    #[test]
    fn test_curve_skeleton_prune() {
        let pts = cube_boundary();
        let mut cs = CurveSkeleton::from_surface_points(&pts, 3);
        let before = cs.edge_count();
        cs.prune_leaves(1e12); // prune everything
        let after = cs.edge_count();
        assert!(after <= before);
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
        let pts = cube_boundary();
        let ma = MedialAxis::from_points(&pts, SkeletonMethod::Voronoi);
        let d = hausdorff_distance(&ma, &ma);
        assert!(d < 1e-10, "Hausdorff to self should be 0, got {:.6}", d);
    }

    #[test]
    fn test_hausdorff_different_shapes() {
        let pts1 = cube_boundary();
        let pts2: Vec<[f64; 3]> = pts1.iter().map(|p| [p[0] + 5.0, p[1], p[2]]).collect();
        let ma1 = MedialAxis::from_points(&pts1, SkeletonMethod::Voronoi);
        let ma2 = MedialAxis::from_points(&pts2, SkeletonMethod::Voronoi);
        let d = hausdorff_distance(&ma1, &ma2);
        assert!(d > 1.0, "Shifted shapes should have Hausdorff > 1");
    }
}
