//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::Transform;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_geometry::Shape;

use super::functions::*;

/// A single iteration record in a GJK debug trace.
#[derive(Debug, Clone)]
pub struct GjkIterationRecord {
    /// Simplex state at start of this iteration.
    pub simplex_size: usize,
    /// Search direction used.
    pub direction: Vec3,
    /// New support point found.
    pub new_support: Vec3,
    /// Whether the dot-product test passed (new_support · direction > 0).
    pub dot_positive: bool,
}
/// The GJK simplex (up to 4 vertices in 3D).
#[derive(Debug, Clone)]
pub struct Simplex {
    /// Vertices of the simplex.
    pub points: Vec<SupportPoint>,
}
impl Simplex {
    /// Create a new empty simplex.
    pub fn new() -> Self {
        Self {
            points: Vec::with_capacity(4),
        }
    }
    /// Push a new point onto the simplex.
    pub fn push(&mut self, point: SupportPoint) {
        self.points.push(point);
    }
    /// Number of points in the simplex.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Whether the simplex is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
    /// Compute the closest point on the simplex to the origin.
    /// Returns (closest_point, barycentric_coordinates).
    pub fn closest_point_to_origin(&self) -> (Vec3, Vec<Real>) {
        match self.points.len() {
            1 => (self.points[0].point, vec![1.0]),
            2 => closest_point_line(&self.points[0].point, &self.points[1].point),
            3 => closest_point_triangle(
                &self.points[0].point,
                &self.points[1].point,
                &self.points[2].point,
            ),
            _ => (Vec3::zeros(), vec![0.25; 4]),
        }
    }
}
/// Result of a GJK ray-cast query.
#[derive(Debug, Clone, Copy)]
pub struct RayCastResult {
    /// Hit parameter `t` along the ray: `origin + t * direction`.
    pub t: f64,
    /// World-space hit point on the shape surface.
    pub point: Vec3,
    /// Outward surface normal at the hit point.
    pub normal: Vec3,
}
/// A frame-coherent GJK solver that carries the previous simplex forward.
///
/// In physics simulation the relative motion of two bodies between frames is
/// small.  Seeding GJK with the simplex from the previous collision query
/// (warm-starting) typically halves the iteration count.
pub struct WarmStartGjk {
    /// Simplex retained from the previous frame.
    pub(super) prev_simplex: Option<Simplex>,
    /// Previous search direction retained for re-use.
    pub(super) prev_direction: Option<Vec3>,
    /// Running average of iterations (for profiling).
    pub(super) avg_iterations: f64,
    /// Number of queries performed.
    pub(super) query_count: usize,
}
impl WarmStartGjk {
    /// Create a new warm-start solver with no cached state.
    pub fn new() -> Self {
        Self {
            prev_simplex: None,
            prev_direction: None,
            avg_iterations: 0.0,
            query_count: 0,
        }
    }
    /// Reset cached state (call when pairing changes).
    pub fn reset(&mut self) {
        self.prev_simplex = None;
        self.prev_direction = None;
    }
    /// Running-average iteration count across all queries.
    pub fn avg_iterations(&self) -> f64 {
        self.avg_iterations
    }
    /// Total number of queries performed.
    pub fn query_count(&self) -> usize {
        self.query_count
    }
    /// Perform a GJK intersection + closest-point query with warm starting.
    ///
    /// If a previous simplex is available it is used as the starting point,
    /// reducing the number of iterations required.
    pub fn query(
        &mut self,
        shape_a: &dyn Shape,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        transform_b: &Transform,
    ) -> GjkResult {
        let mut direction = if let Some(d) = self.prev_direction {
            d
        } else {
            let d = transform_b.position - transform_a.position;
            if d.norm_squared() < TOLERANCE {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                d
            }
        };
        let mut simplex = Simplex::new();
        if let Some(ref ps) = self.prev_simplex {
            for sp in &ps.points {
                if sp.point.dot(&direction) > 0.0 {
                    simplex.push(*sp);
                }
            }
        }
        if simplex.is_empty() {
            let first = support(shape_a, transform_a, shape_b, transform_b, &direction);
            simplex.push(first);
            direction = -first.point;
        } else if let Some(new_dir) = do_simplex(&mut simplex) {
            direction = new_dir;
        }
        let mut iters = 0usize;
        let result = loop {
            iters += 1;
            if iters > MAX_ITERATIONS {
                break GjkResult::Intersecting(simplex.clone());
            }
            let new_pt = support(shape_a, transform_a, shape_b, transform_b, &direction);
            if new_pt.point.dot(&direction) < -TOLERANCE {
                let (closest, bary) = simplex.closest_point_to_origin();
                let dist = closest.norm();
                let (pa, pb) = barycentric_witnesses(&simplex, &bary);
                break GjkResult::Separated {
                    point_a: pa,
                    point_b: pb,
                    distance: dist,
                };
            }
            simplex.push(new_pt);
            if let Some(new_dir) = do_simplex(&mut simplex) {
                direction = new_dir;
            } else {
                break GjkResult::Intersecting(simplex.clone());
            }
        };
        self.prev_direction = Some(direction);
        match &result {
            GjkResult::Intersecting(s) => {
                self.prev_simplex = Some(s.clone());
            }
            GjkResult::Separated {
                point_b, point_a, ..
            } => {
                let d = *point_b - *point_a;
                if d.norm_squared() > TOLERANCE {
                    self.prev_direction = Some(d.normalize());
                }
                self.prev_simplex = None;
            }
        }
        self.query_count += 1;
        let n = self.query_count as f64;
        self.avg_iterations = self.avg_iterations * (n - 1.0) / n + iters as f64 / n;
        result
    }
}
/// Stateful GJK solver with warm-starting and signed-distance capability.
///
/// Stores the result of the previous frame's query so that the next
/// call can reuse the cached simplex as its initial starting point,
/// typically converging in far fewer iterations.
pub struct GjkSolver {
    /// Cached simplex from the last successful query.
    pub(super) cached_simplex: Option<Simplex>,
    /// Last known search direction (used for warm-starting).
    pub(super) last_direction: Option<Vec3>,
}
impl GjkSolver {
    /// Create a new solver with no cached state.
    pub fn new() -> Self {
        Self {
            cached_simplex: None,
            last_direction: None,
        }
    }
    /// Reset cached warm-start state (call when shapes change identity).
    pub fn reset(&mut self) {
        self.cached_simplex = None;
        self.last_direction = None;
    }
    /// Compute the *signed* distance between two convex shapes.
    ///
    /// - Positive: the shapes are separated; value equals the gap.
    /// - Zero: the shapes are exactly touching.
    /// - Negative: the shapes are penetrating; magnitude is an approximation
    ///   of the overlap (for an exact depth, follow up with EPA).
    ///
    /// The solver's internal state is updated with the result.
    pub fn compute_signed_distance(
        &mut self,
        shape_a: &dyn Shape,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        transform_b: &Transform,
    ) -> f64 {
        match Gjk::query(shape_a, transform_a, shape_b, transform_b) {
            GjkResult::Separated {
                distance,
                point_a,
                point_b,
            } => {
                let dir = point_b - point_a;
                self.last_direction = Some(if dir.norm_squared() > 1e-14 {
                    dir.normalize()
                } else {
                    Vec3::new(1.0, 0.0, 0.0)
                });
                self.cached_simplex = None;
                distance
            }
            GjkResult::Intersecting(simplex) => {
                let depth = simplex
                    .points
                    .iter()
                    .map(|p| p.point.norm())
                    .fold(f64::MAX, f64::min);
                self.cached_simplex = Some(simplex);
                self.last_direction = None;
                -(depth.max(0.0))
            }
        }
    }
    /// Compute the closest (witness) points on each shape.
    ///
    /// Returns `Some((point_on_a, point_on_b))` when the shapes are
    /// separated, and `None` when they are intersecting (use EPA for
    /// penetration contacts in that case).
    ///
    /// The returned points are in world space.
    pub fn compute_witness_points(
        &mut self,
        shape_a: &dyn Shape,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        transform_b: &Transform,
    ) -> Option<(Vec3, Vec3)> {
        match Gjk::query(shape_a, transform_a, shape_b, transform_b) {
            GjkResult::Separated {
                point_a,
                point_b,
                distance,
            } => {
                let _ = distance;
                self.cached_simplex = None;
                Some((point_a, point_b))
            }
            GjkResult::Intersecting(simplex) => {
                self.cached_simplex = Some(simplex);
                None
            }
        }
    }
    /// Seed the solver with an initial direction from a previous frame.
    ///
    /// When the relative configuration of the two shapes changes slowly
    /// between frames (the common case in physics simulation), starting GJK
    /// from the previous separation direction reduces the average iteration
    /// count substantially.
    ///
    /// Call this once per frame *before* any query methods.
    pub fn warm_start(&mut self, previous_direction: Vec3) {
        if previous_direction.norm_squared() > 1e-14 {
            self.last_direction = Some(previous_direction.normalize());
        }
    }
    /// Return the warm-start direction stored by the solver.
    ///
    /// Useful for passing across frames or serialising simulation state.
    pub fn warm_direction(&self) -> Option<Vec3> {
        self.last_direction
    }
    /// Whether the solver currently has a cached simplex.
    pub fn has_cached_simplex(&self) -> bool {
        self.cached_simplex.is_some()
    }
    /// Consume and return the cached intersecting simplex (if any).
    ///
    /// This hands the simplex off to an EPA solver and clears the cache.
    pub fn take_simplex(&mut self) -> Option<Simplex> {
        self.cached_simplex.take()
    }
}
/// Generate a speculative contact for CCD.
///
/// When two shapes are separated by less than `threshold`, a speculative contact
/// is generated that allows the constraint solver to keep them apart over the
/// next timestep, preventing tunnelling without requiring CCD.
#[derive(Debug, Clone)]
pub struct SpeculativeContact {
    /// Closest point on shape A.
    pub point_a: Vec3,
    /// Closest point on shape B.
    pub point_b: Vec3,
    /// Contact normal (from B toward A).
    pub normal: Vec3,
    /// Signed gap (negative = penetrating).
    pub gap: f64,
}
/// Voronoi simplex solver for determining the closest feature of a simplex
/// to the origin.
///
/// This is a more robust implementation that explicitly identifies which
/// feature (vertex, edge, or face) of the simplex is closest.
pub struct VoronoiSimplexSolver {
    /// Points in the simplex.
    pub(super) points: Vec<Vec3>,
    /// Barycentric coordinates of the closest point.
    pub(super) bary: Vec<f64>,
    /// Which vertices are active in the closest feature.
    pub(super) active: Vec<bool>,
}
impl Default for VoronoiSimplexSolver {
    fn default() -> Self {
        Self::new()
    }
}
impl VoronoiSimplexSolver {
    /// Create a new empty solver.
    pub fn new() -> Self {
        Self {
            points: Vec::with_capacity(4),
            bary: Vec::with_capacity(4),
            active: Vec::with_capacity(4),
        }
    }
    /// Reset the solver.
    pub fn reset(&mut self) {
        self.points.clear();
        self.bary.clear();
        self.active.clear();
    }
    /// Add a point to the simplex.
    pub fn add_point(&mut self, p: Vec3) {
        self.points.push(p);
        self.bary.push(0.0);
        self.active.push(true);
    }
    /// Number of points in the simplex.
    pub fn num_points(&self) -> usize {
        self.points.len()
    }
    /// Compute the closest point to the origin and identify the closest feature.
    ///
    /// Returns `(closest_point, feature_type)` where feature_type is:
    /// - 0 = vertex
    /// - 1 = edge
    /// - 2 = face
    /// - 3 = interior (origin inside tetrahedron)
    pub fn solve(&mut self) -> (Vec3, ClosestFeature) {
        match self.points.len() {
            1 => {
                self.bary = vec![1.0];
                (self.points[0], ClosestFeature::Vertex(0))
            }
            2 => self.solve_line(),
            3 => self.solve_triangle(),
            4 => self.solve_tetrahedron(),
            _ => (Vec3::zeros(), ClosestFeature::Vertex(0)),
        }
    }
    fn solve_line(&mut self) -> (Vec3, ClosestFeature) {
        let a = &self.points[0];
        let b = &self.points[1];
        let ab = b - a;
        let ao = -a;
        let t = ao.dot(&ab) / ab.dot(&ab);
        if t <= 0.0 {
            self.bary = vec![1.0, 0.0];
            (*a, ClosestFeature::Vertex(0))
        } else if t >= 1.0 {
            self.bary = vec![0.0, 1.0];
            (*b, ClosestFeature::Vertex(1))
        } else {
            self.bary = vec![1.0 - t, t];
            (a + ab * t, ClosestFeature::Edge(0, 1))
        }
    }
    fn solve_triangle(&mut self) -> (Vec3, ClosestFeature) {
        let (pt, bary) = closest_point_triangle(&self.points[0], &self.points[1], &self.points[2]);
        self.bary = bary.clone();
        let num_nonzero = bary.iter().filter(|&&b| b.abs() > 1e-10).count();
        let feature = match num_nonzero {
            1 => {
                let idx = bary.iter().position(|&b| b.abs() > 1e-10).unwrap_or(0);
                ClosestFeature::Vertex(idx)
            }
            2 => {
                let indices: Vec<usize> = bary
                    .iter()
                    .enumerate()
                    .filter(|&(_, &b)| b.abs() > 1e-10)
                    .map(|(i, _)| i)
                    .collect();
                ClosestFeature::Edge(indices[0], indices[1])
            }
            _ => ClosestFeature::Face(0, 1, 2),
        };
        (pt, feature)
    }
    fn solve_tetrahedron(&mut self) -> (Vec3, ClosestFeature) {
        let a = self.points[0];
        let b = self.points[1];
        let c = self.points[2];
        let d = self.points[3];
        let ab = b - a;
        let ac = c - a;
        let ad = d - a;
        let det = ab.dot(&ac.cross(&ad));
        if det.abs() < 1e-14 {
            return self.solve_triangle();
        }
        let ao = -a;
        let abc = ab.cross(&ac);
        let acd = ac.cross(&ad);
        let adb = ad.cross(&ab);
        let inside_abc = abc.dot(&ao) * abc.dot(&(d - a)) >= 0.0;
        let inside_acd = acd.dot(&ao) * acd.dot(&(b - a)) >= 0.0;
        let inside_adb = adb.dot(&ao) * adb.dot(&(c - a)) >= 0.0;
        let bc = c - b;
        let bd = d - b;
        let bo = -b;
        let bcd = bc.cross(&bd);
        let inside_bcd = bcd.dot(&bo) * bcd.dot(&(a - b)) >= 0.0;
        if inside_abc && inside_acd && inside_adb && inside_bcd {
            self.bary = vec![0.25; 4];
            return (Vec3::zeros(), ClosestFeature::Interior);
        }
        let mut best_dist = f64::MAX;
        let mut best_point = Vec3::zeros();
        let mut best_feature = ClosestFeature::Vertex(0);
        let faces: [(usize, usize, usize); 4] = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)];
        for &(i, j, k) in &faces {
            let (pt, _bary) =
                closest_point_triangle(&self.points[i], &self.points[j], &self.points[k]);
            let dist = pt.norm_squared();
            if dist < best_dist {
                best_dist = dist;
                best_point = pt;
                best_feature = ClosestFeature::Face(i, j, k);
            }
        }
        (best_point, best_feature)
    }
}
/// GJK algorithm for convex shape intersection and closest-point queries.
pub struct Gjk;
impl Gjk {
    /// Test whether two convex shapes intersect.
    pub fn intersect(
        shape_a: &dyn Shape,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        transform_b: &Transform,
    ) -> bool {
        matches!(
            Self::query(shape_a, transform_a, shape_b, transform_b),
            GjkResult::Intersecting(_)
        )
    }
    /// Compute closest points between two convex shapes.
    pub fn closest_points(
        shape_a: &dyn Shape,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        transform_b: &Transform,
    ) -> Option<(Vec3, Vec3)> {
        match Self::query(shape_a, transform_a, shape_b, transform_b) {
            GjkResult::Separated {
                point_a, point_b, ..
            } => Some((point_a, point_b)),
            GjkResult::Intersecting(_) => None,
        }
    }
    /// Full GJK query returning either intersection or closest points.
    pub fn query(
        shape_a: &dyn Shape,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        transform_b: &Transform,
    ) -> GjkResult {
        let mut simplex = Simplex::new();
        let mut direction = transform_b.position - transform_a.position;
        if direction.norm_squared() < TOLERANCE {
            direction = Vec3::new(1.0, 0.0, 0.0);
        }
        let first = support(shape_a, transform_a, shape_b, transform_b, &direction);
        simplex.push(first);
        direction = -first.point;
        for _ in 0..MAX_ITERATIONS {
            let new_point = support(shape_a, transform_a, shape_b, transform_b, &direction);
            if new_point.point.dot(&direction) < -TOLERANCE {
                let (closest, bary) = simplex.closest_point_to_origin();
                let _ = bary;
                let dist = closest.norm();
                let point_a = transform_a.position;
                let point_b = transform_b.position;
                return GjkResult::Separated {
                    point_a,
                    point_b,
                    distance: dist,
                };
            }
            simplex.push(new_point);
            if let Some(new_dir) = do_simplex(&mut simplex) {
                direction = new_dir;
            } else {
                return GjkResult::Intersecting(simplex);
            }
        }
        GjkResult::Intersecting(simplex)
    }
}
impl Gjk {
    /// Test intersection between two uniformly-scaled convex shapes.
    ///
    /// `scale_a` and `scale_b` uniformly inflate/deflate the respective shapes.
    pub fn intersect_scaled(
        shape_a: &dyn Shape,
        transform_a: &Transform,
        scale_a: f64,
        shape_b: &dyn Shape,
        transform_b: &Transform,
        scale_b: f64,
    ) -> bool {
        let mut simplex = Simplex::new();
        let mut direction = transform_b.position - transform_a.position;
        if direction.norm_squared() < TOLERANCE {
            direction = Vec3::new(1.0, 0.0, 0.0);
        }
        let first = scaled_support_minkowski(
            shape_a,
            transform_a,
            scale_a,
            shape_b,
            transform_b,
            scale_b,
            &direction,
        );
        simplex.push(first);
        direction = -first.point;
        for _ in 0..MAX_ITERATIONS {
            let new_pt = scaled_support_minkowski(
                shape_a,
                transform_a,
                scale_a,
                shape_b,
                transform_b,
                scale_b,
                &direction,
            );
            if new_pt.point.dot(&direction) < -TOLERANCE {
                return false;
            }
            simplex.push(new_pt);
            if do_simplex(&mut simplex).is_none() {
                return true;
            }
        }
        true
    }
}
/// Result of a GJK query.
pub enum GjkResult {
    /// Shapes are intersecting; simplex is passed to EPA.
    Intersecting(Simplex),
    /// Shapes are separated; closest points returned.
    Separated {
        /// Closest point on shape A.
        point_a: Vec3,
        /// Closest point on shape B.
        point_b: Vec3,
        /// Distance between closest points.
        distance: Real,
    },
}
/// Full debug trace of a GJK run.
#[derive(Debug, Clone)]
pub struct GjkDebugTrace {
    pub(super) iterations: Vec<GjkIterationRecord>,
    pub(super) terminated_intersecting: bool,
}
impl GjkDebugTrace {
    /// Run GJK while recording each iteration.
    pub fn collect(
        shape_a: &dyn Shape,
        transform_a: &Transform,
        shape_b: &dyn Shape,
        transform_b: &Transform,
    ) -> Self {
        let mut simplex = Simplex::new();
        let mut direction = transform_b.position - transform_a.position;
        if direction.norm_squared() < TOLERANCE {
            direction = Vec3::new(1.0, 0.0, 0.0);
        }
        let mut records = Vec::new();
        let first = support(shape_a, transform_a, shape_b, transform_b, &direction);
        simplex.push(first);
        direction = -first.point;
        let mut terminated_intersecting = false;
        for _ in 0..MAX_ITERATIONS {
            let new_point = support(shape_a, transform_a, shape_b, transform_b, &direction);
            let dot = new_point.point.dot(&direction);
            records.push(GjkIterationRecord {
                simplex_size: simplex.len(),
                direction,
                new_support: new_point.point,
                dot_positive: dot >= -TOLERANCE,
            });
            if dot < -TOLERANCE {
                break;
            }
            simplex.push(new_point);
            if let Some(new_dir) = do_simplex(&mut simplex) {
                direction = new_dir;
            } else {
                terminated_intersecting = true;
                break;
            }
        }
        Self {
            iterations: records,
            terminated_intersecting,
        }
    }
    /// Number of GJK iterations performed.
    pub fn n_iterations(&self) -> usize {
        self.iterations.len()
    }
    /// Whether GJK terminated because the origin was found inside the simplex.
    pub fn terminated_as_intersecting(&self) -> bool {
        self.terminated_intersecting
    }
    /// All iteration records.
    pub fn records(&self) -> &[GjkIterationRecord] {
        &self.iterations
    }
    /// All search directions used.
    pub fn directions(&self) -> Vec<Vec3> {
        self.iterations.iter().map(|r| r.direction).collect()
    }
    /// All support points found.
    pub fn support_points(&self) -> Vec<Vec3> {
        self.iterations.iter().map(|r| r.new_support).collect()
    }
    /// Render the trace as a human-readable string (for debugging).
    pub fn summary(&self) -> String {
        let mut s = format!(
            "GJK trace: {} iterations, intersecting={}\n",
            self.iterations.len(),
            self.terminated_intersecting
        );
        for (i, rec) in self.iterations.iter().enumerate() {
            s.push_str(
                &format!(
                    "  iter {i}: simplex_size={}, dir=({:.3},{:.3},{:.3}), support=({:.3},{:.3},{:.3}), dot_ok={}\n",
                    rec.simplex_size, rec.direction.x, rec.direction.y, rec.direction.z,
                    rec.new_support.x, rec.new_support.y, rec.new_support.z, rec
                    .dot_positive,
                ),
            );
        }
        s
    }
}
/// Cache for support function evaluations to avoid redundant computation.
pub struct SupportCache {
    /// Cached directions (unit vectors).
    pub(super) directions: Vec<Vec3>,
    /// Corresponding support points.
    pub(super) results: Vec<SupportPoint>,
    /// Maximum cache size.
    pub(super) max_size: usize,
}
impl SupportCache {
    /// Create a new support cache with the given maximum size.
    pub fn new(max_size: usize) -> Self {
        Self {
            directions: Vec::with_capacity(max_size),
            results: Vec::with_capacity(max_size),
            max_size,
        }
    }
    /// Look up a cached result for a direction.
    ///
    /// Returns `Some(support_point)` if a sufficiently similar direction
    /// was previously queried (dot product > 0.999).
    pub fn lookup(&self, direction: &Vec3) -> Option<&SupportPoint> {
        let dir_norm = direction.norm();
        if dir_norm < 1e-14 {
            return None;
        }
        let unit_dir = direction * (1.0 / dir_norm);
        for (i, cached_dir) in self.directions.iter().enumerate() {
            if unit_dir.dot(cached_dir) > 0.999 {
                return Some(&self.results[i]);
            }
        }
        None
    }
    /// Insert a direction-result pair into the cache.
    pub fn insert(&mut self, direction: Vec3, result: SupportPoint) {
        let dir_norm = direction.norm();
        if dir_norm < 1e-14 {
            return;
        }
        let unit_dir = direction * (1.0 / dir_norm);
        if self.directions.len() >= self.max_size {
            self.directions.remove(0);
            self.results.remove(0);
        }
        self.directions.push(unit_dir);
        self.results.push(result);
    }
    /// Clear the cache.
    pub fn clear(&mut self) {
        self.directions.clear();
        self.results.clear();
    }
    /// Number of cached entries.
    pub fn len(&self) -> usize {
        self.directions.len()
    }
    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.directions.is_empty()
    }
}
/// Identifies the closest feature of a simplex to a query point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ClosestFeature {
    /// Closest feature is a vertex.
    Vertex(usize),
    /// Closest feature is an edge between two vertices.
    Edge(usize, usize),
    /// Closest feature is a triangular face.
    Face(usize, usize, usize),
    /// The query point is inside the tetrahedron.
    Interior,
}
/// Result of an MPR query.
#[derive(Debug, Clone)]
pub enum MprResult {
    /// Shapes are intersecting; contact normal and depth provided.
    Intersecting {
        /// Contact normal (pointing from B to A).
        normal: Vec3,
        /// Penetration depth.
        depth: f64,
        /// Contact point (on boundary between shapes).
        point: Vec3,
    },
    /// Shapes are separated.
    Separated,
}
/// A support point in Minkowski difference space.
#[derive(Debug, Clone, Copy)]
pub struct SupportPoint {
    /// Point in Minkowski difference (A - B).
    pub point: Vec3,
    /// Support point on shape A.
    pub support_a: Vec3,
    /// Support point on shape B.
    pub support_b: Vec3,
}
/// Result of a continuous-collision detection (CCD) sweep.
#[derive(Debug, Clone)]
pub struct CcdResult {
    /// Time of impact in \[0, 1\].
    pub toi: f64,
    /// Contact normal at the time of impact.
    pub normal: Vec3,
    /// Contact point at the time of impact.
    pub point: Vec3,
}
/// Result of a GJK distance query.
#[derive(Debug, Clone)]
pub struct GjkDistanceResult {
    /// Signed distance: positive = separated, negative = penetrating (approx).
    pub distance: f64,
    /// Closest point on shape A (world space).
    pub closest_a: Vec3,
    /// Closest point on shape B (world space).
    pub closest_b: Vec3,
    /// Number of GJK iterations performed.
    pub iterations: usize,
}
