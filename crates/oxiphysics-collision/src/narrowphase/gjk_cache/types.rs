//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions_2::*;
use oxiphysics_core::math::Vec3;
use std::collections::HashMap;

#[allow(unused_imports)]
use super::functions::*;

/// Result of a GJK distance query.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GjkDistanceResult {
    /// Signed distance between the shapes (negative = penetrating).
    pub distance: f64,
    /// Closest point on shape A (world space).
    pub point_a: [f64; 3],
    /// Closest point on shape B (world space).
    pub point_b: [f64; 3],
    /// Number of GJK iterations performed.
    pub iterations: usize,
    /// Whether the query was warm-started.
    pub warm_started: bool,
    /// Whether the shapes are intersecting.
    pub intersecting: bool,
}
#[allow(dead_code)]
impl GjkDistanceResult {
    /// Separation vector from B to A (pointing from B's closest point to A's).
    pub fn separation_vector(&self) -> [f64; 3] {
        sub3_arr(self.point_a, self.point_b)
    }
}
/// A GJK cache entry with a frame-access timestamp for eviction.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TimestampedCacheEntry {
    /// The underlying cache.
    pub cache: GjkCache,
    /// Frame number of last access.
    pub last_accessed: u32,
}
#[allow(dead_code)]
impl TimestampedCacheEntry {
    /// Create a new entry for the current frame.
    pub fn new(frame: u32) -> Self {
        Self {
            cache: GjkCache::new(),
            last_accessed: frame,
        }
    }
    /// Age in frames since last access.
    pub fn age(&self, current_frame: u32) -> u32 {
        current_frame.saturating_sub(self.last_accessed)
    }
}
/// A point in the Configuration Space Obstacle (Minkowski difference) of two shapes.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct CsoPoint {
    /// The Minkowski-difference point: `a_support - b_support`.
    pub p: [f64; 3],
    /// Support point on shape A.
    pub a_support: [f64; 3],
    /// Support point on shape B.
    pub b_support: [f64; 3],
}
/// Barycentric coordinates for the closest point on a simplex to the origin,
/// computed via Johnson's algorithm.
///
/// This is the core of the GJK distance algorithm: given a simplex (up to 4
/// vertices in the Minkowski difference), find the closest point to the origin
/// and the corresponding barycentric weights.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct JohnsonSubResult {
    /// Closest point to origin on the simplex.
    pub closest: [f64; 3],
    /// Barycentric weights for each simplex vertex (sums to 1).
    pub bary: [f64; 4],
    /// Number of active vertices (1–4).
    pub num_active: usize,
    /// Whether the origin is inside the simplex (distance ≈ 0).
    pub origin_inside: bool,
}
#[allow(dead_code)]
impl JohnsonSubResult {
    /// Distance from origin to the closest simplex point.
    pub fn distance(&self) -> f64 {
        len3_arr(self.closest)
    }
    /// Distance squared.
    pub fn distance_sq(&self) -> f64 {
        dot3_arr(self.closest, self.closest)
    }
}
/// Statistics for monitoring the effectiveness of GJK warm starting.
#[derive(Debug, Default, Clone)]
pub struct GjkStats {
    /// Total number of GJK queries processed.
    pub total_queries: u64,
    /// Number of those queries that used a warm-start hint.
    pub warm_started: u64,
}
impl GjkStats {
    /// Create zeroed statistics.
    pub fn new() -> Self {
        Self::default()
    }
    /// Record a query; set `was_warm` to `true` if a cache hint was used.
    pub fn record(&mut self, was_warm: bool) {
        self.total_queries += 1;
        if was_warm {
            self.warm_started += 1;
        }
    }
    /// Fraction of queries that used a warm-start hint (0.0 – 1.0).
    pub fn warm_start_ratio(&self) -> f64 {
        if self.total_queries == 0 {
            0.0
        } else {
            self.warm_started as f64 / self.total_queries as f64
        }
    }
}
/// A registry that associates (body_a_id, body_b_id) pairs with their own
/// `GjkCache` instance for persistent warm starting across frames.
#[allow(dead_code)]
pub struct GjkCacheRegistry {
    pub(super) caches: HashMap<(u64, u64), GjkCache>,
    pub(super) termination: GjkTermination,
}
#[allow(dead_code)]
impl GjkCacheRegistry {
    /// Create a new registry with default termination criteria.
    pub fn new() -> Self {
        Self {
            caches: HashMap::new(),
            termination: GjkTermination::default_criteria(),
        }
    }
    /// Create a registry with custom termination criteria.
    pub fn with_termination(termination: GjkTermination) -> Self {
        Self {
            caches: HashMap::new(),
            termination,
        }
    }
    fn key(id_a: u64, id_b: u64) -> (u64, u64) {
        if id_a <= id_b {
            (id_a, id_b)
        } else {
            (id_b, id_a)
        }
    }
    /// Get or create the GJK cache for the given pair.
    pub fn get_or_create(&mut self, id_a: u64, id_b: u64) -> &mut GjkCache {
        self.caches.entry(Self::key(id_a, id_b)).or_default()
    }
    /// Remove the cache for a given pair (e.g., when a body is destroyed).
    pub fn remove(&mut self, id_a: u64, id_b: u64) {
        self.caches.remove(&Self::key(id_a, id_b));
    }
    /// Remove all caches that reference a specific body ID.
    pub fn remove_body(&mut self, id: u64) {
        self.caches.retain(|&(a, b), _| a != id && b != id);
    }
    /// Run a proximity query with persistent warm starting for a given pair.
    pub fn query<FA, FB>(
        &mut self,
        id_a: u64,
        id_b: u64,
        support_a: &mut FA,
        support_b: &mut FB,
    ) -> ProximityResult
    where
        FA: FnMut([f64; 3]) -> [f64; 3],
        FB: FnMut([f64; 3]) -> [f64; 3],
    {
        let termination = self.termination.clone();
        let cache = self.get_or_create(id_a, id_b);
        gjk_proximity(support_a, support_b, cache, &termination)
    }
    /// Number of cached pairs.
    pub fn len(&self) -> usize {
        self.caches.len()
    }
    /// Whether the registry is empty.
    pub fn is_empty(&self) -> bool {
        self.caches.is_empty()
    }
    /// Clear all caches (e.g., on scene reset).
    pub fn clear(&mut self) {
        self.caches.clear();
    }
    /// Aggregate statistics across all caches.
    pub fn aggregate_stats(&self) -> GjkCacheStats {
        let mut agg = GjkCacheStats::new();
        for cache in self.caches.values() {
            agg.queries += cache.stats.queries;
            agg.iterations += cache.stats.iterations;
            agg.warm_start_saved += cache.stats.warm_start_saved;
        }
        agg
    }
}
/// Shape type tags used by the narrowphase dispatcher.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ShapeType {
    /// Sphere defined by centre + radius.
    Sphere,
    /// Axis-aligned box defined by centre + half-extents.
    Aabb,
    /// Convex hull (support function only).
    ConvexHull,
    /// Capsule (segment + radius).
    Capsule,
}
/// A support point record holding points from both shapes and the
/// Minkowski-difference point.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct CachedSupport {
    /// Support point on shape A in world space.
    pub point_a: [f64; 3],
    /// Support point on shape B in world space.
    pub point_b: [f64; 3],
    /// Minkowski-difference point: `point_a - point_b`.
    pub minkowski_point: [f64; 3],
}
#[allow(dead_code)]
impl CachedSupport {
    /// Create a new cached support from raw arrays.
    pub fn new(point_a: [f64; 3], point_b: [f64; 3]) -> Self {
        Self {
            point_a,
            point_b,
            minkowski_point: [
                point_a[0] - point_b[0],
                point_a[1] - point_b[1],
                point_a[2] - point_b[2],
            ],
        }
    }
    /// Return the zero support.
    pub fn zero() -> Self {
        Self {
            point_a: [0.0; 3],
            point_b: [0.0; 3],
            minkowski_point: [0.0; 3],
        }
    }
}
/// A fixed-capacity cache for support function evaluations.
///
/// Avoids redundant calls to expensive support functions by caching up to
/// `CAP` direction→support mappings for a single shape.
#[allow(dead_code)]
pub struct SupportCache<const CAP: usize> {
    /// Cached directions (unit vectors).
    pub(super) directions: [[f64; 3]; CAP],
    /// Cached support points corresponding to `directions`.
    pub(super) supports: [[f64; 3]; CAP],
    /// Number of valid entries.
    pub(super) count: usize,
    /// Cache hits since last reset.
    pub(super) hits: u32,
    /// Cache misses since last reset.
    pub(super) misses: u32,
}
#[allow(dead_code)]
impl<const CAP: usize> SupportCache<CAP> {
    /// Create an empty support cache.
    pub fn new() -> Self {
        Self {
            directions: [[0.0; 3]; CAP],
            supports: [[0.0; 3]; CAP],
            count: 0,
            hits: 0,
            misses: 0,
        }
    }
    /// Look up a cached support point for the given direction.
    ///
    /// `tol`: cosine similarity threshold (e.g., 0.9999) for a cache hit.
    pub fn lookup(&mut self, dir: [f64; 3], tol: f64) -> Option<[f64; 3]> {
        for i in 0..self.count {
            let cos_sim = dot3_arr(self.directions[i], dir);
            if cos_sim >= tol {
                self.hits += 1;
                return Some(self.supports[i]);
            }
        }
        self.misses += 1;
        None
    }
    /// Insert a direction–support pair into the cache.
    ///
    /// If the cache is full, the oldest entry is evicted (ring buffer).
    pub fn insert(&mut self, dir: [f64; 3], support: [f64; 3]) {
        if CAP == 0 {
            return;
        }
        if self.count < CAP {
            self.directions[self.count] = dir;
            self.supports[self.count] = support;
            self.count += 1;
        } else {
            for i in 0..CAP - 1 {
                self.directions[i] = self.directions[i + 1];
                self.supports[i] = self.supports[i + 1];
            }
            self.directions[CAP - 1] = dir;
            self.supports[CAP - 1] = support;
        }
    }
    /// Reset the cache, clearing all entries and statistics.
    pub fn reset(&mut self) {
        self.count = 0;
        self.hits = 0;
        self.misses = 0;
    }
    /// Cache hit rate (0.0–1.0).
    pub fn hit_rate(&self) -> f64 {
        let total = self.hits + self.misses;
        if total == 0 {
            0.0
        } else {
            self.hits as f64 / total as f64
        }
    }
    /// Number of valid entries.
    pub fn len(&self) -> usize {
        self.count
    }
    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
}
/// Descriptor for a single collision shape with position info.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ShapeDesc {
    /// Unique body/shape id.
    pub id: u64,
    /// Shape type.
    pub shape_type: ShapeType,
    /// Centre/origin in world space.
    pub position: [f64; 3],
    /// Primary shape parameter: radius for spheres/capsules, ignored for others.
    pub radius: f64,
    /// Secondary shape parameter: half-extents for AABB, segment half-length for capsule.
    pub half_extents: [f64; 3],
}
/// Narrowphase dispatch result incorporating the GJK cache.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NarrowphaseResult {
    /// Whether the shapes are overlapping.
    pub overlapping: bool,
    /// Signed separation distance (positive = gap, negative = penetration).
    pub distance: f64,
    /// Closest point on shape A.
    pub point_a: [f64; 3],
    /// Closest point on shape B.
    pub point_b: [f64; 3],
    /// Contact normal (from B toward A).
    pub normal: [f64; 3],
    /// Number of GJK iterations.
    pub iterations: usize,
    /// Whether the cache was used.
    pub warm_started: bool,
}
impl NarrowphaseResult {
    /// Create a result indicating full separation without running GJK.
    #[allow(dead_code)]
    pub fn separated(point_a: [f64; 3], point_b: [f64; 3], distance: f64) -> Self {
        let d = sub3_arr(point_a, point_b);
        let len = len3_arr(d);
        let normal = if len > 1e-10 {
            scale3_arr(d, 1.0 / len)
        } else {
            [0.0, 1.0, 0.0]
        };
        Self {
            overlapping: false,
            distance,
            point_a,
            point_b,
            normal,
            iterations: 0,
            warm_started: false,
        }
    }
}
/// A GJK cache registry with frame-based eviction.
#[allow(dead_code)]
pub struct TimestampedGjkRegistry {
    pub(super) entries: HashMap<(u64, u64), TimestampedCacheEntry>,
    pub(super) current_frame: u32,
    pub(super) policy: EvictionPolicy,
}
#[allow(dead_code)]
impl TimestampedGjkRegistry {
    /// Create a new registry.
    pub fn new(policy: EvictionPolicy) -> Self {
        Self {
            entries: HashMap::new(),
            current_frame: 0,
            policy,
        }
    }
    fn key(id_a: u64, id_b: u64) -> (u64, u64) {
        if id_a <= id_b {
            (id_a, id_b)
        } else {
            (id_b, id_a)
        }
    }
    /// Advance to the next frame, evicting stale entries.
    pub fn advance_frame(&mut self) {
        self.current_frame += 1;
        let frame = self.current_frame;
        let policy = &self.policy;
        self.entries
            .retain(|_, entry| !policy.should_evict(entry.age(frame)));
    }
    /// Get or create a cache entry for the given pair.
    pub fn get_or_create(&mut self, id_a: u64, id_b: u64) -> &mut GjkCache {
        let frame = self.current_frame;
        let entry = self
            .entries
            .entry(Self::key(id_a, id_b))
            .or_insert_with(|| TimestampedCacheEntry::new(frame));
        entry.last_accessed = frame;
        &mut entry.cache
    }
    /// Remove all entries referencing body `id`.
    pub fn remove_body(&mut self, id: u64) {
        self.entries.retain(|&(a, b), _| a != id && b != id);
    }
    /// Number of active entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Whether the registry has no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
/// Result of a proximity query: distance and witness points (or overlap info).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ProximityResult {
    /// Distance between closest features (0 if overlapping).
    pub distance: f64,
    /// Closest point on shape A.
    pub point_a: [f64; 3],
    /// Closest point on shape B.
    pub point_b: [f64; 3],
    /// Contact normal (from B toward A), unit vector. Valid only when `distance < tolerance`.
    pub normal: [f64; 3],
    /// Whether the shapes overlap.
    pub overlapping: bool,
}
/// Result of a batch proximity test.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BatchProximityEntry {
    /// Index of shape A in the batch.
    pub idx_a: usize,
    /// Index of shape B in the batch.
    pub idx_b: usize,
    /// Distance between the shapes.
    pub distance: f64,
    /// Whether the shapes are intersecting.
    pub intersecting: bool,
}
/// GJK overlap tester that warm-starts from a per-pair `GjkPairCache`.
///
/// Tracks `cache_hits` and `cache_misses` to expose a `hit_ratio`.
#[allow(dead_code)]
pub struct WarmStartedGjk {
    /// Per-pair cache.
    pub pair_cache: GjkPairCache,
    /// Number of times the warm-start direction was reused.
    pub cache_hits: u64,
    /// Number of cold-start queries.
    pub cache_misses: u64,
}
#[allow(dead_code)]
impl WarmStartedGjk {
    /// Create a new warm-started GJK with an empty cache.
    pub fn new() -> Self {
        Self {
            pair_cache: GjkPairCache::new(),
            cache_hits: 0,
            cache_misses: 0,
        }
    }
    /// Test whether two shapes overlap.
    ///
    /// `support_a(dir)` and `support_b(dir)` return support points on the
    /// respective shapes.  Note: the caller must handle direction negation
    /// for shape B if needed.
    pub fn test_overlap<FA, FB>(&mut self, mut support_a: FA, mut support_b: FB) -> bool
    where
        FA: FnMut([f64; 3]) -> [f64; 3],
        FB: FnMut([f64; 3]) -> [f64; 3],
    {
        let max_iters = 64;
        let init_dir = if self.pair_cache.is_valid() {
            self.cache_hits += 1;
            self.pair_cache.last_direction
        } else {
            self.cache_misses += 1;
            [1.0, 0.0, 0.0]
        };
        let mut simplex: Vec<CsoPoint> = Vec::with_capacity(4);
        let pa0 = support_a(init_dir);
        let pb0 = support_b(init_dir);
        let cso0 = CsoPoint {
            p: sub3_arr(pa0, pb0),
            a_support: pa0,
            b_support: pb0,
        };
        if dot3_arr(cso0.p, init_dir) < 0.0 {
            self.pair_cache.last_direction = init_dir;
            self.pair_cache.simplex_len = 1;
            self.pair_cache.last_simplex[0] = cso0;
            self.pair_cache.hit = false;
            return false;
        }
        simplex.push(cso0);
        let mut direction = negate3(cso0.p);
        for _ in 0..max_iters {
            let dir_len = len3_arr(direction);
            if dir_len < 1e-10 {
                self.pair_cache.hit = true;
                self.update_cache(&simplex, direction);
                return true;
            }
            let norm_dir = scale3_arr(direction, 1.0 / dir_len);
            let pa = support_a(norm_dir);
            let pb = support_b(norm_dir);
            let new_cso = CsoPoint {
                p: sub3_arr(pa, pb),
                a_support: pa,
                b_support: pb,
            };
            if dot3_arr(new_cso.p, norm_dir) < 0.0 {
                self.pair_cache.hit = false;
                self.update_cache(&simplex, norm_dir);
                return false;
            }
            simplex.push(new_cso);
            let extracted: Vec<[f64; 3]> = simplex.iter().map(|c| c.p).collect();
            let mut mink_simplex: Vec<CachedSupport> = extracted
                .iter()
                .map(|&p| CachedSupport {
                    point_a: [0.0; 3],
                    point_b: [0.0; 3],
                    minkowski_point: p,
                })
                .collect();
            if let Some(new_dir) = do_simplex_cached(&mut mink_simplex) {
                let new_len = mink_simplex.len();
                simplex.truncate(new_len);
                direction = new_dir;
            } else {
                self.pair_cache.hit = true;
                self.update_cache(&simplex, norm_dir);
                return true;
            }
        }
        self.pair_cache.hit = false;
        false
    }
    fn update_cache(&mut self, simplex: &[CsoPoint], direction: [f64; 3]) {
        let n = simplex.len().min(4);
        let zero_cso = CsoPoint {
            p: [0.0; 3],
            a_support: [0.0; 3],
            b_support: [0.0; 3],
        };
        let mut arr = [zero_cso; 4];
        arr[..n].copy_from_slice(&simplex[..n]);
        self.pair_cache.last_simplex = arr;
        self.pair_cache.simplex_len = n;
        self.pair_cache.last_direction = direction;
    }
    /// Hit ratio: fraction of queries that used the warm-start cache (0.0–1.0).
    pub fn hit_ratio(&self) -> f64 {
        let total = self.cache_hits + self.cache_misses;
        if total == 0 {
            0.0
        } else {
            self.cache_hits as f64 / total as f64
        }
    }
}
/// A positioned GJK cache that stores body positions alongside the simplex,
/// so validity can be checked automatically.
#[allow(dead_code)]
pub struct PositionedGjkCache {
    /// The GJK cache data.
    pub cache: GjkCache,
    /// Recorded position of body A when the cache was last updated.
    pub pos_a: [f64; 3],
    /// Recorded position of body B when the cache was last updated.
    pub pos_b: [f64; 3],
    /// Whether the cache has been populated at least once.
    pub populated: bool,
}
#[allow(dead_code)]
impl PositionedGjkCache {
    /// Create an empty positioned cache.
    pub fn new() -> Self {
        Self {
            cache: GjkCache::new(),
            pos_a: [0.0; 3],
            pos_b: [0.0; 3],
            populated: false,
        }
    }
    /// Check whether the cache is valid for the given current positions.
    pub fn is_valid_for(&self, current_pos_a: [f64; 3], current_pos_b: [f64; 3]) -> bool {
        self.populated && cache_still_valid(self.pos_a, self.pos_b, current_pos_a, current_pos_b)
    }
    /// Update the cache with new positions and simplex data.
    pub fn update(&mut self, pos_a: [f64; 3], pos_b: [f64; 3]) {
        self.pos_a = pos_a;
        self.pos_b = pos_b;
        self.populated = true;
    }
    /// Invalidate the cache (e.g., after a teleport).
    pub fn invalidate(&mut self) {
        self.cache.reset();
        self.populated = false;
    }
}
/// Cached state for a single GJK pair: last simplex (up to 4 vertices) and
/// last search direction.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GjkPairCache {
    /// Last simplex: up to 4 CSO vertices.
    pub last_simplex: [CsoPoint; 4],
    /// Number of valid entries in `last_simplex`.
    pub simplex_len: usize,
    /// Last search direction used.
    pub last_direction: [f64; 3],
    /// Whether the last query detected overlap.
    pub hit: bool,
}
#[allow(dead_code)]
impl GjkPairCache {
    /// Create an empty (invalid) cache.
    pub fn new() -> Self {
        let zero_cso = CsoPoint {
            p: [0.0; 3],
            a_support: [0.0; 3],
            b_support: [0.0; 3],
        };
        Self {
            last_simplex: [zero_cso; 4],
            simplex_len: 0,
            last_direction: [1.0, 0.0, 0.0],
            hit: false,
        }
    }
    /// Whether this cache has valid data.
    pub fn is_valid(&self) -> bool {
        self.simplex_len > 0
    }
}
/// Per-pair GJK warm-start manager.
///
/// Call [`GjkWarmStart::cached_query_hint`] at the start of each GJK query to
/// obtain a warm-start direction, then call [`GjkWarmStart::update_cache`]
/// after the query completes.
pub struct GjkWarmStart {
    /// Cache keyed by an ordered (smaller_id, larger_id) pair so insertion
    /// order doesn't matter.
    pub(super) cache: HashMap<(u64, u64), SimplexCache>,
}
impl GjkWarmStart {
    /// Create an empty warm-start cache.
    pub fn new() -> Self {
        Self {
            cache: HashMap::new(),
        }
    }
    /// Canonical key: always put the smaller id first.
    fn key(id_a: u64, id_b: u64) -> (u64, u64) {
        if id_a <= id_b {
            (id_a, id_b)
        } else {
            (id_b, id_a)
        }
    }
    /// Return the cached initial search direction for the pair, if any.
    ///
    /// Returns `None` on a cache miss (first query for this pair).
    pub fn get_initial_direction(&self, id_a: u64, id_b: u64) -> Option<Vec3> {
        self.cache.get(&Self::key(id_a, id_b))?.direction()
    }
    /// Store the result of a GJK query so it can warm-start the next frame.
    pub fn update_cache(&mut self, id_a: u64, id_b: u64, witness_a: Vec3, witness_b: Vec3) {
        let entry = self.cache.entry(Self::key(id_a, id_b)).or_default();
        entry.witness_a = witness_a;
        entry.witness_b = witness_b;
        entry.simplex_bits = 0b0001;
        entry.cache_valid = true;
    }
    /// Convenience wrapper: return `(witness_b − witness_a).normalize()` as
    /// the warm-start hint for the next GJK query, or `None` on cache miss.
    pub fn cached_query_hint(&self, id_a: u64, id_b: u64) -> Option<Vec3> {
        self.get_initial_direction(id_a, id_b)
    }
    /// Remove a cache entry (e.g. when a body is destroyed).
    pub fn invalidate(&mut self, id_a: u64, id_b: u64) {
        self.cache.remove(&Self::key(id_a, id_b));
    }
    /// Number of cached pairs.
    pub fn len(&self) -> usize {
        self.cache.len()
    }
    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }
}
/// Controls when the GJK distance iteration terminates.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GjkTermination {
    /// Absolute distance tolerance: stop when `|dist - prev_dist| < abs_tol`.
    pub abs_tol: f64,
    /// Relative distance tolerance: stop when `|dist - prev_dist| / dist < rel_tol`.
    pub rel_tol: f64,
    /// Maximum number of iterations.
    pub max_iters: usize,
    /// Stop immediately if `dist < proximity_tol` (near-intersection).
    pub proximity_tol: f64,
}
#[allow(dead_code)]
impl GjkTermination {
    /// Default conservative termination criteria.
    pub fn default_criteria() -> Self {
        Self {
            abs_tol: 1e-8,
            rel_tol: 1e-6,
            max_iters: 64,
            proximity_tol: 1e-10,
        }
    }
    /// Tight criteria for high-accuracy distance queries.
    pub fn tight() -> Self {
        Self {
            abs_tol: 1e-12,
            rel_tol: 1e-10,
            max_iters: 128,
            proximity_tol: 1e-14,
        }
    }
    /// Loose criteria for broad-phase proximity queries.
    pub fn loose() -> Self {
        Self {
            abs_tol: 1e-4,
            rel_tol: 1e-3,
            max_iters: 16,
            proximity_tol: 1e-6,
        }
    }
    /// Check whether the algorithm should terminate given the current and
    /// previous distance values and the number of iterations so far.
    pub fn should_terminate(&self, dist: f64, prev_dist: f64, iter: usize) -> bool {
        if iter >= self.max_iters {
            return true;
        }
        if dist < self.proximity_tol {
            return true;
        }
        let delta = (prev_dist - dist).abs();
        if delta < self.abs_tol {
            return true;
        }
        if dist > 0.0 && delta / dist < self.rel_tol {
            return true;
        }
        false
    }
}
/// Policy controlling when stale GJK cache entries are evicted.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EvictionPolicy {
    /// Maximum number of frames a cache entry can survive without being accessed.
    pub max_age_frames: u32,
    /// Maximum total entries allowed before evicting oldest.
    pub max_entries: usize,
}
#[allow(dead_code)]
impl EvictionPolicy {
    /// Conservative policy: keep entries for 60 frames, cap at 1024.
    pub fn conservative() -> Self {
        Self {
            max_age_frames: 60,
            max_entries: 1024,
        }
    }
    /// Aggressive policy: keep entries for only 10 frames, cap at 256.
    pub fn aggressive() -> Self {
        Self {
            max_age_frames: 10,
            max_entries: 256,
        }
    }
    /// Check whether a cache entry with age `frame_age` should be evicted.
    pub fn should_evict(&self, frame_age: u32) -> bool {
        frame_age > self.max_age_frames
    }
}
/// A raw-array GJK cache that stores simplex vertices and last search
/// direction for warm-starting subsequent queries.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GjkCache {
    /// Cached simplex vertices from the previous query.
    pub simplex_vertices: Vec<CachedSupport>,
    /// Last search direction used.
    pub last_direction: [f64; 3],
    /// Whether warm starting is enabled.
    pub warm_start: bool,
    /// Accumulated performance statistics.
    pub stats: GjkCacheStats,
}
#[allow(dead_code)]
impl GjkCache {
    /// Create a new cache with warm starting enabled.
    pub fn new() -> Self {
        Self {
            simplex_vertices: Vec::new(),
            last_direction: [1.0, 0.0, 0.0],
            warm_start: true,
            stats: GjkCacheStats::new(),
        }
    }
    /// Reset the cache, discarding all stored data.
    pub fn reset(&mut self) {
        self.simplex_vertices.clear();
        self.last_direction = [1.0, 0.0, 0.0];
    }
    /// Warm start from a previous frame's simplex.
    pub fn warm_start_from(&mut self, prev_simplex: &[CachedSupport]) {
        self.simplex_vertices = prev_simplex.to_vec();
        if let Some(last) = prev_simplex.last() {
            let p = last.minkowski_point;
            let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            if len > 1e-10 {
                self.last_direction = [-p[0] / len, -p[1] / len, -p[2] / len];
            }
        }
        self.warm_start = true;
    }
    /// Check intersection using cached simplex data.
    ///
    /// `support_fn(direction) -> CachedSupport` returns a support point
    /// in the Minkowski difference for the given search direction.
    pub fn intersect_cached<F>(&mut self, mut support_fn: F) -> bool
    where
        F: FnMut([f64; 3]) -> CachedSupport,
    {
        let max_iters = 64;
        let mut simplex: Vec<CachedSupport> = Vec::new();
        let init_dir = if self.warm_start && !self.simplex_vertices.is_empty() {
            self.last_direction
        } else {
            [1.0, 0.0, 0.0]
        };
        let first = support_fn(init_dir);
        if dot3_arr(first.minkowski_point, init_dir) < 0.0 {
            self.stats.record_query(1, self.warm_start, 10);
            return false;
        }
        simplex.push(first);
        let mut direction = negate3(first.minkowski_point);
        for iter in 0..max_iters {
            let dir_len = len3_arr(direction);
            if dir_len < 1e-10 {
                self.simplex_vertices = simplex;
                self.stats
                    .record_query(iter as u64 + 1, self.warm_start, 10);
                return true;
            }
            let norm_dir = scale3_arr(direction, 1.0 / dir_len);
            let new_pt = support_fn(norm_dir);
            if dot3_arr(new_pt.minkowski_point, norm_dir) < 0.0 {
                self.stats
                    .record_query(iter as u64 + 1, self.warm_start, 10);
                return false;
            }
            simplex.push(new_pt);
            if let Some(new_dir) = do_simplex_cached(&mut simplex) {
                direction = new_dir;
            } else {
                self.simplex_vertices = simplex;
                self.last_direction = norm_dir;
                self.stats
                    .record_query(iter as u64 + 1, self.warm_start, 10);
                return true;
            }
        }
        self.simplex_vertices = simplex;
        false
    }
    /// Compute closest points between two shapes using cached data.
    ///
    /// Returns `Some((closest_on_a, closest_on_b))` if the shapes are separated,
    /// or `None` if they are intersecting.
    pub fn closest_points_cached<F>(&mut self, mut support_fn: F) -> Option<([f64; 3], [f64; 3])>
    where
        F: FnMut([f64; 3]) -> CachedSupport,
    {
        let max_iters = 64;
        let tol = 1e-8;
        let init_dir = if self.warm_start && !self.simplex_vertices.is_empty() {
            self.last_direction
        } else {
            [1.0, 0.0, 0.0]
        };
        let mut simplex: Vec<CachedSupport> = Vec::new();
        let first = support_fn(init_dir);
        simplex.push(first);
        let mut closest_a = first.point_a;
        let mut closest_b = first.point_b;
        let mut prev_dist_sq = f64::INFINITY;
        for _ in 0..max_iters {
            let closest_mk = closest_point_on_simplex_to_origin(&simplex);
            let dist_sq = dot3_arr(closest_mk, closest_mk);
            if dist_sq < tol {
                return None;
            }
            if (prev_dist_sq - dist_sq).abs() < tol * tol {
                self.simplex_vertices = simplex;
                return Some((closest_a, closest_b));
            }
            prev_dist_sq = dist_sq;
            let direction = negate3(closest_mk);
            let dir_len = len3_arr(direction);
            if dir_len < 1e-10 {
                return None;
            }
            let norm_dir = scale3_arr(direction, 1.0 / dir_len);
            let new_pt = support_fn(norm_dir);
            let proj = dot3_arr(new_pt.minkowski_point, norm_dir);
            let cur_dist = dist_sq.sqrt();
            if cur_dist - proj < tol {
                closest_a = new_pt.point_a;
                closest_b = new_pt.point_b;
                self.simplex_vertices = simplex;
                return Some((closest_a, closest_b));
            }
            simplex.push(new_pt);
            if simplex.len() > 4 {
                simplex.remove(0);
            }
            closest_a = new_pt.point_a;
            closest_b = new_pt.point_b;
        }
        self.simplex_vertices = simplex;
        Some((closest_a, closest_b))
    }
}
/// Cached simplex / witness data from the previous frame's GJK query.
#[derive(Debug, Clone)]
pub struct SimplexCache {
    /// Last cached closest point on body A (world space).
    pub witness_a: Vec3,
    /// Last cached closest point on body B (world space).
    pub witness_b: Vec3,
    /// Bitmask encoding which simplex vertices were part of the cached simplex.
    pub simplex_bits: u8,
    /// Whether this cache entry is valid and can be used as a warm-start hint.
    pub cache_valid: bool,
}
impl SimplexCache {
    /// Create a new, invalid cache entry.
    pub fn new() -> Self {
        Self {
            witness_a: Vec3::zeros(),
            witness_b: Vec3::zeros(),
            simplex_bits: 0,
            cache_valid: false,
        }
    }
    /// Return the cached search direction (B − A, normalised), if valid.
    pub fn direction(&self) -> Option<Vec3> {
        if !self.cache_valid {
            return None;
        }
        let d = self.witness_b - self.witness_a;
        let len = d.norm();
        if len > 1e-10 { Some(d / len) } else { None }
    }
}
/// Performance statistics for cached GJK queries.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct GjkCacheStats {
    /// Total iterations across all queries.
    pub iterations: u64,
    /// Number of iterations saved by warm starting.
    pub warm_start_saved: u64,
    /// Number of queries performed.
    pub queries: u64,
}
#[allow(dead_code)]
impl GjkCacheStats {
    /// Create zeroed stats.
    pub fn new() -> Self {
        Self::default()
    }
    /// Record a query with the given iteration count and whether it was warm started.
    pub fn record_query(&mut self, iters: u64, was_warm: bool, cold_estimate: u64) {
        self.queries += 1;
        self.iterations += iters;
        if was_warm && cold_estimate > iters {
            self.warm_start_saved += cold_estimate - iters;
        }
    }
    /// Average iterations per query.
    pub fn avg_iterations(&self) -> f64 {
        if self.queries == 0 {
            0.0
        } else {
            self.iterations as f64 / self.queries as f64
        }
    }
}
/// Extended hit-rate tracker that also monitors validity rejection rate.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct HitRateStats {
    /// Total warm-start attempts.
    pub attempts: u64,
    /// Attempts that passed validity check.
    pub valid_hits: u64,
    /// Attempts that failed validity check (cache was stale).
    pub stale_misses: u64,
    /// Cold-start queries (no cache at all).
    pub cold_starts: u64,
}
#[allow(dead_code)]
impl HitRateStats {
    /// Create zeroed stats.
    pub fn new() -> Self {
        Self::default()
    }
    /// Record a warm-start attempt.
    pub fn record_warm_attempt(&mut self, valid: bool) {
        self.attempts += 1;
        if valid {
            self.valid_hits += 1;
        } else {
            self.stale_misses += 1;
        }
    }
    /// Record a cold start (no cached data).
    pub fn record_cold_start(&mut self) {
        self.cold_starts += 1;
    }
    /// Effective hit rate: valid hits / (valid hits + stale misses + cold starts).
    pub fn effective_hit_rate(&self) -> f64 {
        let total = self.valid_hits + self.stale_misses + self.cold_starts;
        if total == 0 {
            0.0
        } else {
            self.valid_hits as f64 / total as f64
        }
    }
    /// Staleness rate: how often cache data is present but stale.
    pub fn staleness_rate(&self) -> f64 {
        if self.attempts == 0 {
            0.0
        } else {
            self.stale_misses as f64 / self.attempts as f64
        }
    }
}
/// Stateful incremental GJK contact detector for a persistent shape pair.
///
/// Maintains warm-start data across multiple frames.
#[allow(dead_code)]
pub struct GjkContactPair {
    /// The underlying GJK cache.
    pub cache: GjkCache,
    /// Termination criteria.
    pub termination: GjkTermination,
    /// Number of frames this pair has been tracked.
    pub age: u32,
    /// Last known result: whether the shapes were intersecting.
    pub last_intersecting: bool,
    /// Last known distance (0 if intersecting).
    pub last_distance: f64,
}
#[allow(dead_code)]
impl GjkContactPair {
    /// Create a new contact pair with default settings.
    pub fn new() -> Self {
        Self {
            cache: GjkCache::new(),
            termination: GjkTermination::default_criteria(),
            age: 0,
            last_intersecting: false,
            last_distance: f64::INFINITY,
        }
    }
    /// Run a GJK query and update the pair state.
    ///
    /// Returns the query result.
    pub fn query<F>(&mut self, support_fn: &mut F) -> GjkDistanceResult
    where
        F: FnMut([f64; 3]) -> CachedSupport,
    {
        self.age += 1;
        let result = gjk_distance(support_fn, &mut self.cache, &self.termination);
        self.last_intersecting = result.intersecting;
        self.last_distance = result.distance;
        result
    }
    /// Whether the pair is currently overlapping (based on last query).
    pub fn is_overlapping(&self) -> bool {
        self.last_intersecting
    }
    /// Distance from the last query.
    pub fn distance(&self) -> f64 {
        self.last_distance
    }
}
