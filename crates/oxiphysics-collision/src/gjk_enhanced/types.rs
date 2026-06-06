//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    ConvexShape, johnson_segment, johnson_tetrahedron, johnson_triangle, vlen, vlen_sq, vneg, vsub,
};

/// A sphere support shape.
#[derive(Debug, Clone, Copy)]
pub struct EnhSphere {
    /// Centre of the sphere.
    pub centre: [f64; 3],
    /// Radius.
    pub radius: f64,
}
impl EnhSphere {
    /// Create a new sphere.
    pub fn new(centre: [f64; 3], radius: f64) -> Self {
        Self { centre, radius }
    }
}
/// Helper: a shape translated by an offset.
pub(super) struct TranslatedShape<'a> {
    pub(super) shape: &'a dyn ConvexShape,
    pub(super) offset: [f64; 3],
}
/// A Minkowski difference support point carrying witness points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SupportPoint {
    /// Point in the Minkowski difference: `w = support_a(d) - support_b(-d)`.
    pub w: [f64; 3],
    /// Witness point on shape A.
    pub a: [f64; 3],
    /// Witness point on shape B.
    pub b: [f64; 3],
}
impl SupportPoint {
    /// Create a new support point.
    pub fn new(w: [f64; 3], a: [f64; 3], b: [f64; 3]) -> Self {
        Self { w, a, b }
    }
}
/// A cached simplex used to warm-start the next GJK query.
///
/// Stores up to 4 support points from the previous frame.
#[derive(Debug, Clone, Default)]
pub struct WarmStartCache {
    /// Cached simplex vertices.
    pub points: Vec<SupportPoint>,
    /// Last known search direction.
    pub direction: [f64; 3],
    /// Whether the cache is valid for use.
    pub valid: bool,
}
impl WarmStartCache {
    /// Create an empty warm-start cache.
    pub fn new() -> Self {
        Self {
            points: Vec::with_capacity(4),
            direction: [1.0, 0.0, 0.0],
            valid: false,
        }
    }
    /// Invalidate the cache (call after large body motions).
    pub fn invalidate(&mut self) {
        self.valid = false;
        self.points.clear();
    }
    /// Update the cache after a successful GJK run.
    pub fn update(&mut self, points: &[SupportPoint], direction: [f64; 3]) {
        self.points.clear();
        self.points.extend_from_slice(points);
        self.direction = direction;
        self.valid = true;
    }
}
/// Result of the EPA penetration-depth query.
#[derive(Debug, Clone)]
pub struct EpaEnhancedResult {
    /// Penetration depth.
    pub depth: f64,
    /// Contact normal pointing from B into A.
    pub normal: [f64; 3],
    /// Contact point (midpoint on both shapes).
    pub contact_point: [f64; 3],
    /// Number of EPA iterations performed.
    pub iterations: u32,
}
/// Unified result from the combined GJK-EPA pipeline.
#[derive(Debug, Clone)]
pub struct ContactResult {
    /// Whether the shapes are in contact (intersecting).
    pub in_contact: bool,
    /// Separation distance (positive if separated, negative if penetrating).
    pub signed_distance: f64,
    /// Contact or closest point on shape A.
    pub point_a: [f64; 3],
    /// Contact or closest point on shape B.
    pub point_b: [f64; 3],
    /// Normal pointing from B to A (only meaningful on contact).
    pub normal: [f64; 3],
    /// GJK metrics.
    pub gjk_metrics: GjkMetrics,
    /// EPA iterations (0 if no EPA was run).
    pub epa_iterations: u32,
}
/// Result of an enhanced GJK query.
#[derive(Debug, Clone)]
pub struct GjkResult {
    /// Whether the two shapes are intersecting.
    pub intersecting: bool,
    /// Closest point on shape A (only valid if not intersecting).
    pub closest_a: [f64; 3],
    /// Closest point on shape B (only valid if not intersecting).
    pub closest_b: [f64; 3],
    /// Squared distance between closest points (0 if intersecting).
    pub dist_sq: f64,
    /// The final simplex.
    pub simplex: Vec<SupportPoint>,
    /// Performance metrics.
    pub metrics: GjkMetrics,
}
impl GjkResult {
    /// Distance between the closest points (0 if intersecting).
    pub fn distance(&self) -> f64 {
        self.dist_sq.sqrt()
    }
}
/// Approximate a torus as a convex shape for GJK by treating it as a
/// rounded-rectangle cross-section sweep (approximate support).
///
/// The torus has major radius `R` (from centre to tube centre) and minor
/// radius `r` (tube radius).
#[derive(Debug, Clone, Copy)]
pub struct TorusApprox {
    /// Major radius.
    pub major_r: f64,
    /// Minor radius (tube radius).
    pub minor_r: f64,
    /// Centre of the torus.
    pub centre: [f64; 3],
}
impl TorusApprox {
    /// Create a new approximate torus shape.
    pub fn new(centre: [f64; 3], major_r: f64, minor_r: f64) -> Self {
        Self {
            major_r,
            minor_r,
            centre,
        }
    }
}
/// Cached support query result for reducing redundant support evaluations.
#[derive(Debug, Clone)]
pub struct SupportCache {
    /// Cached direction (last query direction).
    pub last_dir: [f64; 3],
    /// Cached support point for shape A.
    pub support_a: [f64; 3],
    /// Cached support point for shape B.
    pub support_b: [f64; 3],
    /// Number of cache hits.
    pub hits: u64,
    /// Total number of queries.
    pub total: u64,
}
impl SupportCache {
    /// Create a new empty support cache.
    pub fn new() -> Self {
        Self {
            last_dir: [f64::NAN; 3],
            support_a: [0.0; 3],
            support_b: [0.0; 3],
            hits: 0,
            total: 0,
        }
    }
    /// Attempt to retrieve cached results for the given direction.
    ///
    /// Returns `Some((sa, sb))` if the direction matches within tolerance.
    pub fn get(&mut self, dir: [f64; 3], tol: f64) -> Option<([f64; 3], [f64; 3])> {
        self.total += 1;
        if self.last_dir[0].is_nan() {
            return None;
        }
        let diff = vsub(dir, self.last_dir);
        if vlen(diff) < tol {
            self.hits += 1;
            Some((self.support_a, self.support_b))
        } else {
            None
        }
    }
    /// Store a new support result for the given direction.
    pub fn put(&mut self, dir: [f64; 3], sa: [f64; 3], sb: [f64; 3]) {
        self.last_dir = dir;
        self.support_a = sa;
        self.support_b = sb;
    }
    /// Cache hit rate as a fraction in \[0, 1\].
    pub fn hit_rate(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        self.hits as f64 / self.total as f64
    }
}
/// Result of a time-of-impact (TOI) query.
#[derive(Debug, Clone)]
pub struct ToiResult {
    /// Whether a collision occurs during `[0, 1]`.
    pub hit: bool,
    /// Time of impact in `[0, 1]` (1.0 if no hit).
    pub toi: f64,
    /// Contact normal at impact.
    pub normal: [f64; 3],
    /// Number of GJK iterations used.
    pub iterations: u32,
}
/// A GJK simplex of dimension 0–3 (point, line, triangle, tetrahedron).
#[derive(Debug, Clone)]
pub struct EnhSimplex {
    /// Simplex vertices (up to 4).
    pub verts: Vec<SupportPoint>,
}
impl EnhSimplex {
    /// Create an empty simplex.
    pub fn new() -> Self {
        Self {
            verts: Vec::with_capacity(4),
        }
    }
    /// Create a simplex pre-populated from a warm-start cache.
    pub fn from_cache(cache: &WarmStartCache) -> Self {
        Self {
            verts: cache.points.clone(),
        }
    }
    /// Current dimension (number of vertices − 1; −1 if empty).
    pub fn dim(&self) -> i32 {
        self.verts.len() as i32 - 1
    }
    /// Add a new vertex to the simplex.
    pub fn add(&mut self, p: SupportPoint) {
        self.verts.push(p);
    }
    /// Reduce the simplex to the sub-simplex closest to the origin.
    ///
    /// Returns the new search direction and whether the origin is contained.
    pub fn reduce(&mut self) -> ([f64; 3], bool) {
        match self.verts.len() {
            0 => ([1.0, 0.0, 0.0], false),
            1 => {
                let w = vneg(self.verts[0].w);
                (w, vlen_sq(self.verts[0].w) < 1e-18)
            }
            2 => self.reduce_line(),
            3 => self.reduce_triangle(),
            4 => self.reduce_tetrahedron(),
            _ => ([1.0, 0.0, 0.0], false),
        }
    }
    fn reduce_line(&mut self) -> ([f64; 3], bool) {
        let a = self.verts[1].w;
        let b = self.verts[0].w;
        let (_, _, pt) = johnson_segment(a, b);
        let dist = vlen_sq(pt);
        if dist < 1e-18 {
            return ([1.0, 0.0, 0.0], true);
        }
        let (la, lb, _) = johnson_segment(a, b);
        if lb < 1e-10 {
            self.verts.remove(0);
        } else if la < 1e-10 {
            self.verts.remove(1);
        }
        (vneg(pt), false)
    }
    fn reduce_triangle(&mut self) -> ([f64; 3], bool) {
        let a = self.verts[2].w;
        let b = self.verts[1].w;
        let c = self.verts[0].w;
        let (la, lb, lc, pt) = johnson_triangle(a, b, c);
        let dist = vlen_sq(pt);
        if dist < 1e-18 {
            return ([1.0, 0.0, 0.0], true);
        }
        let mut new_verts = Vec::with_capacity(3);
        let weights = [lc, lb, la];
        for (i, &w) in weights.iter().enumerate() {
            if w > 1e-10 {
                new_verts.push(self.verts[i]);
            }
        }
        self.verts = new_verts;
        (vneg(pt), false)
    }
    fn reduce_tetrahedron(&mut self) -> ([f64; 3], bool) {
        let a = self.verts[3].w;
        let b = self.verts[2].w;
        let c = self.verts[1].w;
        let d = self.verts[0].w;
        let (pt, inside) = johnson_tetrahedron(a, b, c, d);
        if inside {
            return ([1.0, 0.0, 0.0], true);
        }
        let faces = [(3usize, 2usize, 1usize), (3, 2, 0), (3, 1, 0), (2, 1, 0)];
        let _ = pt;
        let mut best_dist = f64::INFINITY;
        let mut best_face = (3, 2, 1);
        for &(i, j, k) in &faces {
            let wi = self.verts[i].w;
            let wj = self.verts[j].w;
            let wk = self.verts[k].w;
            let (_, _, _, fpt) = johnson_triangle(wi, wj, wk);
            let d = vlen_sq(fpt);
            if d < best_dist {
                best_dist = d;
                best_face = (i, j, k);
            }
        }
        let (fi, fj, fk) = best_face;
        self.verts = vec![self.verts[fk], self.verts[fj], self.verts[fi]];
        let (_, _, _, closest) =
            johnson_triangle(self.verts[2].w, self.verts[1].w, self.verts[0].w);
        (vneg(closest), false)
    }
}
/// A capsule support shape (cylinder with hemispherical caps).
#[derive(Debug, Clone, Copy)]
pub struct EnhCapsule {
    /// Start point of the capsule axis.
    pub a: [f64; 3],
    /// End point of the capsule axis.
    pub b: [f64; 3],
    /// Radius.
    pub radius: f64,
}
impl EnhCapsule {
    /// Create a new capsule.
    pub fn new(a: [f64; 3], b: [f64; 3], radius: f64) -> Self {
        Self { a, b, radius }
    }
}
/// Minkowski sum of two convex shapes (for rounded shapes).
///
/// The support of `A ⊕ B` in direction `d` is `support_A(d) + support_B(d)`.
pub struct MinkowskiSum<'a> {
    /// First shape.
    pub a: &'a dyn ConvexShape,
    /// Second shape (e.g. a sphere for rounding).
    pub b: &'a dyn ConvexShape,
}
impl<'a> MinkowskiSum<'a> {
    /// Create a Minkowski sum from two shapes.
    pub fn new(a: &'a dyn ConvexShape, b: &'a dyn ConvexShape) -> Self {
        Self { a, b }
    }
}
/// An axis-aligned box support shape.
#[derive(Debug, Clone, Copy)]
pub struct EnhBox {
    /// Centre of the box.
    pub centre: [f64; 3],
    /// Half-extents.
    pub half: [f64; 3],
}
impl EnhBox {
    /// Create a new axis-aligned box.
    pub fn new(centre: [f64; 3], half: [f64; 3]) -> Self {
        Self { centre, half }
    }
}
/// A convex hull support shape (arbitrary point cloud).
#[derive(Debug, Clone)]
pub struct ConvexHull {
    /// Vertices of the convex hull.
    pub vertices: Vec<[f64; 3]>,
}
impl ConvexHull {
    /// Create a convex hull from a list of vertices.
    pub fn new(vertices: Vec<[f64; 3]>) -> Self {
        Self { vertices }
    }
}
/// An ellipsoid support shape.
#[derive(Debug, Clone, Copy)]
pub struct EnhEllipsoid {
    /// Centre.
    pub centre: [f64; 3],
    /// Semi-axes (a, b, c).
    pub semi_axes: [f64; 3],
}
impl EnhEllipsoid {
    /// Create a new ellipsoid.
    pub fn new(centre: [f64; 3], semi_axes: [f64; 3]) -> Self {
        Self { centre, semi_axes }
    }
}
/// Performance metrics collected during a GJK run.
#[derive(Debug, Clone, Default)]
pub struct GjkMetrics {
    /// Number of iterations performed.
    pub iterations: u32,
    /// Number of support function calls.
    pub support_calls: u32,
    /// Number of simplex reductions performed.
    pub simplex_reductions: u32,
    /// Whether warm-starting was used.
    pub warm_started: bool,
    /// Whether the algorithm terminated early (convergence before max iter).
    pub converged_early: bool,
    /// Final distance squared (or 0 if intersecting).
    pub final_dist_sq: f64,
}
impl GjkMetrics {
    /// Create a new zeroed metrics structure.
    pub fn new() -> Self {
        Self::default()
    }
}
/// An EPA polytope face.
#[derive(Debug, Clone, Copy)]
pub(super) struct EpaFaceE {
    pub(super) indices: [usize; 3],
    pub(super) normal: [f64; 3],
    pub(super) dist: f64,
}
