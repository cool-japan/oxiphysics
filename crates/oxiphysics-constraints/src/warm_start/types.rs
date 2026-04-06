//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{cross3, dist_sq_3, dot3};

/// A single cached impulse with aging information.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct WarmStartRecord {
    /// 3-D impulse vector.
    pub impulse: [f64; 3],
    /// Number of frames since this record was refreshed.
    pub age: u32,
    /// The contact pair this record belongs to.
    pub pair_id: (u32, u32),
}
#[allow(dead_code)]
impl WarmStartRecord {
    /// Create a fresh record (age = 0).
    pub fn new(pair_id: (u32, u32), impulse: [f64; 3]) -> Self {
        Self {
            impulse,
            age: 0,
            pair_id,
        }
    }
    /// Scale the impulse by `decay_factor`, simulating energy loss over time.
    pub fn decay(&mut self, decay_factor: f64) {
        self.impulse[0] *= decay_factor;
        self.impulse[1] *= decay_factor;
        self.impulse[2] *= decay_factor;
    }
}
/// Tracks how much work warm-starting saves compared to cold-starting.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct WarmStartContribStats {
    /// Number of contact pairs where a cached impulse was available.
    pub hits: u64,
    /// Number of contact pairs with no cached impulse.
    pub misses: u64,
    /// Cumulative squared magnitude of applied warm-start impulses.
    pub total_impulse_sq: f64,
    /// Largest single warm-start impulse magnitude observed.
    pub peak_impulse: f64,
}
#[allow(dead_code)]
impl WarmStartContribStats {
    /// Create zeroed statistics.
    pub fn new() -> Self {
        Self::default()
    }
    /// Record a warm-start hit with the given impulse.
    pub fn record_hit(&mut self, impulse: [f64; 3]) {
        self.hits += 1;
        let mag_sq = impulse[0] * impulse[0] + impulse[1] * impulse[1] + impulse[2] * impulse[2];
        self.total_impulse_sq += mag_sq;
        let mag = mag_sq.sqrt();
        if mag > self.peak_impulse {
            self.peak_impulse = mag;
        }
    }
    /// Record a warm-start miss (no cached impulse).
    pub fn record_miss(&mut self) {
        self.misses += 1;
    }
    /// Hit rate in \[0, 1\] (0 when no data).
    pub fn hit_rate(&self) -> f64 {
        let total = self.hits + self.misses;
        if total == 0 {
            0.0
        } else {
            self.hits as f64 / total as f64
        }
    }
    /// Root-mean-square warm-start impulse magnitude.
    pub fn rms_impulse(&self) -> f64 {
        if self.hits == 0 {
            0.0
        } else {
            (self.total_impulse_sq / self.hits as f64).sqrt()
        }
    }
    /// Reset all counters.
    pub fn reset(&mut self) {
        *self = Self::default();
    }
}
/// Measures how effective warm starting was during a solve frame.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct WarmStartQualityMetrics {
    /// Number of constraint pairs that had a warm-start cache hit.
    pub hits: u32,
    /// Number of constraint pairs that had no cache (cold start).
    pub misses: u32,
    /// Sum of |warm_impulse - final_impulse| / |final_impulse| across hits.
    /// Lower is better (the warm-start guess was close to the solution).
    pub relative_error_sum: f64,
    /// Number of solver iterations actually used this frame.
    pub iterations_used: u32,
    /// Residual velocity error after solving.
    pub residual: f64,
}
#[allow(dead_code)]
impl WarmStartQualityMetrics {
    /// Create empty metrics for a new frame.
    pub fn new() -> Self {
        Self {
            hits: 0,
            misses: 0,
            relative_error_sum: 0.0,
            iterations_used: 0,
            residual: 0.0,
        }
    }
    /// Record a cache hit and the quality of the warm-start guess.
    ///
    /// * `warm_magnitude` — magnitude of the re-applied warm impulse.
    /// * `final_magnitude` — magnitude of the converged impulse.
    pub fn record_hit(&mut self, warm_magnitude: f64, final_magnitude: f64) {
        self.hits += 1;
        if final_magnitude.abs() > 1e-12 {
            self.relative_error_sum +=
                (warm_magnitude - final_magnitude).abs() / final_magnitude.abs();
        }
    }
    /// Record a cache miss (cold start).
    pub fn record_miss(&mut self) {
        self.misses += 1;
    }
    /// Hit rate in \[0, 1\].
    pub fn hit_rate(&self) -> f64 {
        let total = self.hits + self.misses;
        if total == 0 {
            0.0
        } else {
            self.hits as f64 / total as f64
        }
    }
    /// Average relative error for cache hits. Lower is better.
    pub fn mean_relative_error(&self) -> f64 {
        if self.hits == 0 {
            0.0
        } else {
            self.relative_error_sum / self.hits as f64
        }
    }
    /// A composite quality score in \[0, 1\] where 1 = perfect warm start.
    ///
    /// Combines hit rate and accuracy: `quality = hit_rate * (1 - clamp(mean_error, 0, 1))`.
    pub fn quality_score(&self) -> f64 {
        let accuracy = (1.0 - self.mean_relative_error().clamp(0.0, 1.0)).max(0.0);
        self.hit_rate() * accuracy
    }
}
/// Decides whether to apply a cached impulse based on how much a body's
/// velocity has changed since the impulse was stored.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AdaptiveVelocityGate {
    /// Maximum allowed squared linear-velocity change before warm-start is
    /// skipped.  Typical value: (0.5 m/s)^2 = 0.25.
    pub max_delta_sq: f64,
}
#[allow(dead_code)]
impl AdaptiveVelocityGate {
    /// Create the gate with a given velocity-change threshold.
    pub fn new(max_velocity_delta: f64) -> Self {
        Self {
            max_delta_sq: max_velocity_delta * max_velocity_delta,
        }
    }
    /// Return `true` when the cached impulse should be applied.
    ///
    /// Returns `false` (skip warm-start) when the body's current linear
    /// velocity differs from the snapshot by more than `max_velocity_delta`.
    pub fn should_apply(&self, snapshot: &BodyVelocitySnapshot, current_linear: [f64; 3]) -> bool {
        snapshot.linear_delta_sq(current_linear) <= self.max_delta_sq
    }
    /// Return the warm-start impulse if the velocity gate allows it, otherwise
    /// return `[0, 0, 0]`.
    pub fn apply_or_zero(
        &self,
        snapshot: &BodyVelocitySnapshot,
        current_linear: [f64; 3],
        cached_impulse: [f64; 3],
    ) -> [f64; 3] {
        if self.should_apply(snapshot, current_linear) {
            cached_impulse
        } else {
            [0.0, 0.0, 0.0]
        }
    }
}
/// A warm-start candidate entry for quality-based selection.
#[allow(dead_code)]
pub struct WarmStartCandidate {
    /// The cached impulse data.
    pub cache: WarmStartCache,
    /// Frame age of this entry (0 = current frame).
    pub age: u32,
    /// Quality score in \[0, 1\].
    pub quality: f64,
}
#[allow(dead_code)]
impl WarmStartCandidate {
    /// Create a new candidate.
    pub fn new(cache: WarmStartCache, age: u32, quality: f64) -> Self {
        Self {
            cache,
            age,
            quality,
        }
    }
    /// Compute a composite score: quality / (1 + age).
    pub fn composite_score(&self) -> f64 {
        self.quality / (1.0 + self.age as f64)
    }
}
/// Identifies a contact point for cross-frame matching.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ContactFingerprint {
    /// Body pair key (canonical ordering).
    pub body_pair: (u64, u64),
    /// Contact point on body A in local space.
    pub local_point_a: [f64; 3],
    /// Contact point on body B in local space.
    pub local_point_b: [f64; 3],
    /// Contact normal direction.
    pub normal: [f64; 3],
}
/// Tracks age-based decay for cached impulses.
///
/// Each entry records how many frames old the cached impulse is.
/// After `max_age` frames without renewal the entry is evicted.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ImpulseAging {
    /// Per-pair age counters (number of frames since last active contact).
    pub(super) ages: HashMap<(u64, u64), u32>,
    /// Maximum age before an entry is evicted.
    pub(super) max_age: u32,
    /// Per-frame decay multiplier applied to cached impulses.
    pub(super) decay_factor: f64,
}
#[allow(dead_code)]
impl ImpulseAging {
    /// Create a new aging tracker.
    ///
    /// * `max_age` — frames to keep a cached impulse without renewal.
    /// * `decay_factor` — multiplier applied each frame (e.g. 0.9).
    pub fn new(max_age: u32, decay_factor: f64) -> Self {
        Self {
            ages: HashMap::new(),
            max_age,
            decay_factor: decay_factor.clamp(0.0, 1.0),
        }
    }
    /// Mark a pair as active this frame (resets its age to 0).
    pub fn touch(&mut self, a: u64, b: u64) {
        let k = WarmStartMap::key(a, b);
        self.ages.insert(k, 0);
    }
    /// Advance one frame: increment all ages and remove entries that exceed max_age.
    /// Returns the list of evicted keys.
    pub fn tick(&mut self) -> Vec<(u64, u64)> {
        let mut evicted = Vec::new();
        self.ages.retain(|&k, age| {
            *age += 1;
            if *age > self.max_age {
                evicted.push(k);
                false
            } else {
                true
            }
        });
        evicted
    }
    /// Compute the decay multiplier for a given pair based on its age.
    /// Returns `decay_factor ^ age`, or 0.0 if the pair is not tracked.
    pub fn decay_multiplier(&self, a: u64, b: u64) -> f64 {
        let k = WarmStartMap::key(a, b);
        match self.ages.get(&k) {
            Some(&age) => self.decay_factor.powi(age as i32),
            None => 0.0,
        }
    }
    /// Return the current age of a pair, or `None` if not tracked.
    pub fn age_of(&self, a: u64, b: u64) -> Option<u32> {
        let k = WarmStartMap::key(a, b);
        self.ages.get(&k).copied()
    }
    /// Number of entries currently tracked.
    pub fn len(&self) -> usize {
        self.ages.len()
    }
    /// Whether the tracker is empty.
    pub fn is_empty(&self) -> bool {
        self.ages.is_empty()
    }
    /// Apply aging decay to a warm-start map, evicting expired entries.
    pub fn apply_to_map(&mut self, map: &mut WarmStartMap) {
        let evicted = self.tick();
        for k in &evicted {
            map.cache.remove(k);
        }
        let keys: Vec<(u64, u64)> = map.cache.keys().copied().collect();
        for k in keys {
            let multiplier = match self.ages.get(&k) {
                Some(&age) => self.decay_factor.powi(age as i32),
                None => 0.0,
            };
            if multiplier < 1e-12 {
                map.cache.remove(&k);
            } else if let Some(entry) = map.cache.get_mut(&k) {
                entry.lambda_n *= multiplier;
                entry.lambda_t1 *= multiplier;
                entry.lambda_t2 *= multiplier;
            }
        }
    }
}
/// Cached accumulated impulses from the previous simulation frame.
///
/// These are stored per contact pair and re-applied at the start of each frame
/// to warm-start the solver.
#[derive(Debug, Clone, Default)]
pub struct WarmStartCache {
    /// Accumulated normal impulse from previous frame.
    pub lambda_n: f64,
    /// Accumulated first tangent impulse from previous frame.
    pub lambda_t1: f64,
    /// Accumulated second tangent impulse from previous frame.
    pub lambda_t2: f64,
}
impl WarmStartCache {
    /// Create a cache entry with explicit impulse values.
    #[allow(dead_code)]
    pub fn with_impulses(lambda_n: f64, lambda_t1: f64, lambda_t2: f64) -> Self {
        Self {
            lambda_n,
            lambda_t1,
            lambda_t2,
        }
    }
    /// Compute the total impulse magnitude (L2 norm of all three components).
    #[allow(dead_code)]
    pub fn magnitude(&self) -> f64 {
        (self.lambda_n * self.lambda_n
            + self.lambda_t1 * self.lambda_t1
            + self.lambda_t2 * self.lambda_t2)
            .sqrt()
    }
    /// Scale all impulse components by a factor.
    #[allow(dead_code)]
    pub fn scale(&self, factor: f64) -> Self {
        Self {
            lambda_n: self.lambda_n * factor,
            lambda_t1: self.lambda_t1 * factor,
            lambda_t2: self.lambda_t2 * factor,
        }
    }
    /// Return true if all impulse components are effectively zero.
    #[allow(dead_code)]
    pub fn is_negligible(&self, epsilon: f64) -> bool {
        self.lambda_n.abs() < epsilon
            && self.lambda_t1.abs() < epsilon
            && self.lambda_t2.abs() < epsilon
    }
    /// Linearly interpolate between this cache and another.
    #[allow(dead_code)]
    pub fn lerp(&self, other: &Self, t: f64) -> Self {
        let t = t.clamp(0.0, 1.0);
        Self {
            lambda_n: self.lambda_n + t * (other.lambda_n - self.lambda_n),
            lambda_t1: self.lambda_t1 + t * (other.lambda_t1 - self.lambda_t1),
            lambda_t2: self.lambda_t2 + t * (other.lambda_t2 - self.lambda_t2),
        }
    }
}
/// Matches contact constraint pairs across simulation frames.
///
/// Uses spatial proximity of local-space contact points to decide whether
/// a contact from the current frame corresponds to one from the previous frame.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ConstraintPairMatcher {
    /// Position tolerance squared for matching contact points.
    pub(super) position_tolerance_sq: f64,
    /// Normal direction tolerance (dot product threshold; 1.0 = exact match).
    pub(super) normal_tolerance: f64,
    /// Previous frame fingerprints keyed by body pair.
    pub(super) prev_fingerprints: HashMap<(u64, u64), Vec<ContactFingerprint>>,
}
#[allow(dead_code)]
impl ConstraintPairMatcher {
    /// Create a new matcher with given tolerances.
    ///
    /// * `position_tolerance` — max distance between matching contact points.
    /// * `normal_tolerance` — min dot product between normals (e.g. 0.95).
    pub fn new(position_tolerance: f64, normal_tolerance: f64) -> Self {
        Self {
            position_tolerance_sq: position_tolerance * position_tolerance,
            normal_tolerance: normal_tolerance.clamp(0.0, 1.0),
            prev_fingerprints: HashMap::new(),
        }
    }
    /// Store the current frame's fingerprints for use in the next frame.
    pub fn push_frame(&mut self, fingerprints: &[ContactFingerprint]) {
        self.prev_fingerprints.clear();
        for fp in fingerprints {
            self.prev_fingerprints
                .entry(fp.body_pair)
                .or_default()
                .push(fp.clone());
        }
    }
    /// Try to find a matching previous-frame fingerprint for the given current-frame contact.
    ///
    /// Returns the index into the previous frame's fingerprint list for this body pair,
    /// or `None` if no match is found.
    pub fn find_match(&self, current: &ContactFingerprint) -> Option<usize> {
        let prev_list = self.prev_fingerprints.get(&current.body_pair)?;
        for (i, prev) in prev_list.iter().enumerate() {
            let dist_a_sq = dist_sq_3(current.local_point_a, prev.local_point_a);
            let dist_b_sq = dist_sq_3(current.local_point_b, prev.local_point_b);
            let dot_n = dot3(current.normal, prev.normal);
            if dist_a_sq <= self.position_tolerance_sq
                && dist_b_sq <= self.position_tolerance_sq
                && dot_n >= self.normal_tolerance
            {
                return Some(i);
            }
        }
        None
    }
    /// Number of body pairs stored from the previous frame.
    pub fn prev_pair_count(&self) -> usize {
        self.prev_fingerprints.len()
    }
    /// Total number of fingerprints stored from the previous frame.
    pub fn prev_fingerprint_count(&self) -> usize {
        self.prev_fingerprints.values().map(|v| v.len()).sum()
    }
    /// Clear all stored fingerprints.
    pub fn clear(&mut self) {
        self.prev_fingerprints.clear();
    }
}
/// Adaptive impulse threshold that adjusts based on recent impulse statistics.
///
/// The threshold grows when large impulses are observed (body is active) and
/// shrinks back toward a floor when everything is quiet.
#[allow(dead_code)]
pub struct AdaptiveImpulseThreshold {
    /// Current threshold value.
    pub threshold: f64,
    /// Minimum (floor) threshold.
    pub floor: f64,
    /// Maximum (ceiling) threshold.
    pub ceiling: f64,
    /// Rate at which the threshold adapts upward (per impulse observation).
    pub up_rate: f64,
    /// Rate at which the threshold decays toward the floor each frame.
    pub decay_rate: f64,
}
#[allow(dead_code)]
impl AdaptiveImpulseThreshold {
    /// Create an adaptive threshold with default parameters.
    pub fn new(floor: f64, ceiling: f64) -> Self {
        Self {
            threshold: floor,
            floor,
            ceiling,
            up_rate: 0.1,
            decay_rate: 0.02,
        }
    }
    /// Observe an impulse magnitude and update the threshold.
    pub fn observe(&mut self, impulse_magnitude: f64) {
        if impulse_magnitude > self.threshold {
            self.threshold = (self.threshold + self.up_rate * impulse_magnitude).min(self.ceiling);
        }
    }
    /// Advance one frame: decay the threshold toward the floor.
    pub fn tick(&mut self) {
        self.threshold = (self.threshold * (1.0 - self.decay_rate)).max(self.floor);
    }
    /// Return true if the impulse magnitude exceeds the threshold.
    pub fn exceeds(&self, impulse_magnitude: f64) -> bool {
        impulse_magnitude > self.threshold
    }
    /// Reset to floor threshold.
    pub fn reset(&mut self) {
        self.threshold = self.floor;
    }
}
/// Cached impulse data for an entire contact island, transferred from one
/// frame to the next.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct IslandWarmStartData {
    /// Island identifier.
    pub island_id: IslandId,
    /// Per-pair impulse cache for pairs within this island.
    pub pairs: HashMap<(u32, u32), [f64; 3]>,
    /// Frame age: 0 = freshly recorded, increases each frame.
    pub age: u32,
}
#[allow(dead_code)]
impl IslandWarmStartData {
    /// Create empty warm-start data for the given island.
    pub fn new(island_id: IslandId) -> Self {
        Self {
            island_id,
            pairs: HashMap::new(),
            age: 0,
        }
    }
    /// Store the warm-start impulse for a contact pair within the island.
    pub fn insert_pair(&mut self, pair: (u32, u32), impulse: [f64; 3]) {
        let key = if pair.0 <= pair.1 {
            pair
        } else {
            (pair.1, pair.0)
        };
        self.pairs.insert(key, impulse);
    }
    /// Retrieve the warm-start impulse for a contact pair, if present.
    pub fn get_pair(&self, pair: (u32, u32)) -> Option<[f64; 3]> {
        let key = if pair.0 <= pair.1 {
            pair
        } else {
            (pair.1, pair.0)
        };
        self.pairs.get(&key).copied()
    }
    /// Increment the age counter by one.
    pub fn increment_age(&mut self) {
        self.age += 1;
    }
    /// Scale all cached impulses by `scale` (for alpha-decay).
    pub fn decay_all(&mut self, scale: f64) {
        for v in self.pairs.values_mut() {
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;
        }
    }
    /// Number of contact pairs in this island.
    pub fn pair_count(&self) -> usize {
        self.pairs.len()
    }
}
/// HashMap-backed cache of warm-start impulses keyed by contact pair.
#[allow(dead_code)]
pub struct WarmStartImpulseCache {
    pub(super) records: HashMap<(u32, u32), WarmStartRecord>,
    /// Maximum age before a record is considered stale and evicted.
    pub max_age: u32,
}
#[allow(dead_code)]
impl WarmStartImpulseCache {
    /// Create an empty cache with the given maximum age.
    pub fn new(max_age: u32) -> Self {
        Self {
            records: HashMap::new(),
            max_age,
        }
    }
    /// Insert (or overwrite) a record for the given pair with a fresh impulse.
    pub fn insert(&mut self, pair: (u32, u32), impulse: [f64; 3]) {
        self.records
            .insert(pair, WarmStartRecord::new(pair, impulse));
    }
    /// Look up the cached impulse for a pair, if present.
    pub fn lookup(&self, pair: (u32, u32)) -> Option<[f64; 3]> {
        self.records.get(&pair).map(|r| r.impulse)
    }
    /// Increment the age of every record by 1.
    pub fn age_all(&mut self) {
        for record in self.records.values_mut() {
            record.age += 1;
        }
    }
    /// Remove all records whose age exceeds `max_age`.
    pub fn evict_aged(&mut self) {
        let max_age = self.max_age;
        self.records.retain(|_, r| r.age <= max_age);
    }
    /// Return quality classification for a cached pair.
    pub fn quality(&self, pair: (u32, u32)) -> WarmStartQuality {
        match self.records.get(&pair) {
            None => WarmStartQuality::Stale,
            Some(r) if r.age == 0 => WarmStartQuality::Fresh,
            Some(r) if r.age <= self.max_age => {
                WarmStartQuality::Aged(1.0 - r.age as f64 / self.max_age as f64)
            }
            _ => WarmStartQuality::Stale,
        }
    }
}
/// Combines a `WarmStartImpulseCache` with a quality threshold to decide
/// whether to apply cached impulses or fall back to the candidate.
#[allow(dead_code)]
pub struct AdaptiveWarmStart {
    /// Internal cache of impulse records.
    pub cache: WarmStartImpulseCache,
    /// Minimum quality score (0..=1) required to use cached data.
    pub quality_threshold: f64,
}
#[allow(dead_code)]
impl AdaptiveWarmStart {
    /// Create a new adaptive warm-start with the given max age and threshold.
    pub fn new(max_age: u32, quality_threshold: f64) -> Self {
        Self {
            cache: WarmStartImpulseCache::new(max_age),
            quality_threshold,
        }
    }
    /// Return the cached impulse for `pair` if its quality meets the threshold,
    /// otherwise return `candidate_impulse` unchanged.
    pub fn apply_if_valid(&self, pair: (u32, u32), candidate_impulse: [f64; 3]) -> [f64; 3] {
        match self.cache.quality(pair) {
            WarmStartQuality::Fresh => self.cache.lookup(pair).unwrap_or(candidate_impulse),
            WarmStartQuality::Aged(score) if score >= self.quality_threshold => {
                self.cache.lookup(pair).unwrap_or(candidate_impulse)
            }
            _ => candidate_impulse,
        }
    }
}
/// Velocity snapshot used to decide whether warm-start is still valid.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, Default)]
pub struct BodyVelocitySnapshot {
    /// Linear velocity at the time the impulse was cached.
    pub linear: [f64; 3],
    /// Angular velocity at the time the impulse was cached.
    pub angular: [f64; 3],
}
#[allow(dead_code)]
impl BodyVelocitySnapshot {
    /// Create a snapshot from linear and angular velocities.
    pub fn new(linear: [f64; 3], angular: [f64; 3]) -> Self {
        Self { linear, angular }
    }
    /// Compute the squared L2 distance between this snapshot and a current
    /// velocity reading.
    ///
    /// Only the linear component is compared — this is sufficient for most
    /// stacking scenarios.
    pub fn linear_delta_sq(&self, current: [f64; 3]) -> f64 {
        let d = [
            current[0] - self.linear[0],
            current[1] - self.linear[1],
            current[2] - self.linear[2],
        ];
        d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
    }
}
/// Circular buffer storing impulse histories over N frames.
///
/// Useful for detecting sleeping bodies (very small impulses over many frames)
/// and for providing higher-order warm-start predictions.
#[allow(dead_code)]
pub struct ImpulseHistory {
    /// Ring buffer of impulse magnitudes.
    pub(super) data: Vec<f64>,
    /// Write head index.
    pub(super) head: usize,
    /// Number of valid entries.
    pub(super) count: usize,
}
#[allow(dead_code)]
impl ImpulseHistory {
    /// Create a new impulse history buffer with the given capacity.
    pub fn new(capacity: usize) -> Self {
        let capacity = capacity.max(1);
        Self {
            data: vec![0.0; capacity],
            head: 0,
            count: 0,
        }
    }
    /// Push a new impulse magnitude, overwriting the oldest if full.
    pub fn push(&mut self, magnitude: f64) {
        let cap = self.data.len();
        self.data[self.head] = magnitude;
        self.head = (self.head + 1) % cap;
        if self.count < cap {
            self.count += 1;
        }
    }
    /// Mean of all stored values.
    pub fn mean(&self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }
        self.data[..self.count].iter().sum::<f64>() / self.count as f64
    }
    /// Maximum of all stored values.
    pub fn max_val(&self) -> f64 {
        self.data[..self.count]
            .iter()
            .copied()
            .fold(0.0_f64, f64::max)
    }
    /// True if all stored values are below `threshold`.
    pub fn all_below(&self, threshold: f64) -> bool {
        self.data[..self.count].iter().all(|&v| v < threshold)
    }
    /// Number of stored entries.
    pub fn len(&self) -> usize {
        self.count
    }
    /// True if the buffer is empty.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }
    /// Number of entries the buffer can hold.
    pub fn capacity(&self) -> usize {
        self.data.len()
    }
    /// Clear all stored values.
    pub fn clear(&mut self) {
        self.data.fill(0.0);
        self.head = 0;
        self.count = 0;
    }
}
/// Matches current-frame contact points to previous-frame contact points by
/// proximity, enabling per-point impulse warm-starting.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PersistentManifoldMatcher {
    /// Maximum distance (squared) for two contact points to be considered the
    /// same across frames.
    pub proximity_sq: f64,
}
#[allow(dead_code)]
impl PersistentManifoldMatcher {
    /// Create a matcher with the given proximity radius.
    pub fn new(proximity_radius: f64) -> Self {
        Self {
            proximity_sq: proximity_radius * proximity_radius,
        }
    }
    /// For each point in `current` find the index of the closest point in
    /// `previous` (if within proximity).  Returns `None` for unmatched points.
    pub fn match_points(
        &self,
        current: &[ContactPoint3D],
        previous: &[ContactPoint3D],
    ) -> Vec<Option<usize>> {
        current
            .iter()
            .map(|cp| {
                let mut best_idx = None;
                let mut best_dist_sq = self.proximity_sq;
                for (i, pp) in previous.iter().enumerate() {
                    let d = [
                        cp.position[0] - pp.position[0],
                        cp.position[1] - pp.position[1],
                        cp.position[2] - pp.position[2],
                    ];
                    let dist_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                    if dist_sq <= best_dist_sq {
                        best_dist_sq = dist_sq;
                        best_idx = Some(i);
                    }
                }
                best_idx
            })
            .collect()
    }
    /// Count how many current points were successfully matched.
    pub fn match_count(&self, current: &[ContactPoint3D], previous: &[ContactPoint3D]) -> usize {
        self.match_points(current, previous)
            .iter()
            .filter(|m| m.is_some())
            .count()
    }
}
/// Exponential moving average aging for impulse caches.
///
/// Each frame the stored impulse is blended toward zero:
/// `impulse = impulse * (1 - alpha)` where `alpha ∈ [0, 1]`.
///
/// This differs from `AlphaImpulseAging` in that the blend factor is applied
/// independently per entry via a per-entry weight.
#[allow(dead_code)]
pub struct ExponentialImpulseAging {
    /// Blend factor per frame (0 = no decay, 1 = instant reset).
    pub alpha: f64,
    /// Minimum impulse magnitude below which the entry is dropped entirely.
    pub prune_threshold: f64,
}
#[allow(dead_code)]
impl ExponentialImpulseAging {
    /// Create a new exponential aging with blend factor and prune threshold.
    pub fn new(alpha: f64, prune_threshold: f64) -> Self {
        let alpha = alpha.clamp(0.0, 1.0);
        Self {
            alpha,
            prune_threshold,
        }
    }
    /// Age a single impulse value.
    pub fn age_scalar(&self, impulse: f64) -> f64 {
        impulse * (1.0 - self.alpha)
    }
    /// Age a cache entry; returns `None` if the result is negligible.
    pub fn age_cache(&self, cache: &WarmStartCache) -> Option<WarmStartCache> {
        let aged = cache.scale(1.0 - self.alpha);
        if aged.is_negligible(self.prune_threshold) {
            None
        } else {
            Some(aged)
        }
    }
    /// Age all entries in a `WarmStartMap` in place, pruning negligible ones.
    pub fn apply_to_map(&self, map: &mut WarmStartMap) {
        let factor = 1.0 - self.alpha;
        map.scale_all(factor);
        map.prune(self.prune_threshold);
    }
}
/// An identifier for a contact island (a connected component of touching
/// bodies).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct IslandId(pub u32);
/// A velocity-level contact solver with warm starting.
///
/// Computes normal and friction impulses using a sequential-impulse approach.
/// Warm starting accelerates convergence by re-applying the previous frame's
/// accumulated impulses as an initial guess.
#[allow(dead_code)]
pub struct ContactVelocitySolver {
    /// Coefficient of restitution (bounciness), in \[0, 1\].
    pub restitution: f64,
    /// Coefficient of friction (Coulomb), >= 0.
    pub friction: f64,
    /// Penetration slop: allowed overlap before Baumgarte kicks in.
    pub slop: f64,
    /// Baumgarte stabilization factor, typically 0.1–0.3.
    pub baumgarte: f64,
}
#[allow(dead_code)]
impl ContactVelocitySolver {
    /// Create a new solver with the given restitution and friction coefficients.
    ///
    /// Uses default values for `slop` (0.005) and `baumgarte` (0.2).
    pub fn new(restitution: f64, friction: f64) -> Self {
        Self {
            restitution,
            friction,
            slop: 0.005,
            baumgarte: 0.2,
        }
    }
    /// Compute the normal impulse increment Δλ_n with warm starting.
    ///
    /// Returns `(delta_lambda_n, updated_cache)` where the updated cache has the
    /// new accumulated normal impulse stored (tangent impulses are carried through
    /// unchanged — update them via `solve_friction_impulse`).
    ///
    /// # Parameters
    /// - `rel_vel_n`: relative velocity along the contact normal (positive = separating)
    /// - `inv_mass_a`, `inv_mass_b`: inverse masses of bodies A and B
    /// - `inv_inertia_a`, `inv_inertia_b`: scalar inverse inertia (diagonal approximation)
    /// - `jacobian_n`: normal direction `[nx, ny, nz]`
    /// - `r_a`, `r_b`: lever arms from body CoM to contact point `[x, y, z]`
    /// - `penetration`: overlap depth (positive = penetrating)
    /// - `dt`: time step
    /// - `warm`: warm-start cache from the previous frame
    #[allow(clippy::too_many_arguments)]
    pub fn solve_normal_impulse(
        rel_vel_n: f64,
        inv_mass_a: f64,
        inv_mass_b: f64,
        inv_inertia_a: f64,
        inv_inertia_b: f64,
        jacobian_n: [f64; 3],
        r_a: [f64; 3],
        r_b: [f64; 3],
        penetration: f64,
        dt: f64,
        warm: &WarmStartCache,
    ) -> (f64, WarmStartCache) {
        let cross_a = cross3(r_a, jacobian_n);
        let cross_b = cross3(r_b, jacobian_n);
        let k = inv_mass_a
            + inv_mass_b
            + dot3(cross_a, cross_a) * inv_inertia_a
            + dot3(cross_b, cross_b) * inv_inertia_b;
        if k < 1e-12 {
            return (0.0, warm.clone());
        }
        let effective_mass = 1.0 / k;
        let baumgarte = 0.2_f64;
        let slop = 0.005_f64;
        let bias = (baumgarte / dt) * (penetration - slop).max(0.0);
        let delta_lambda = effective_mass * (-rel_vel_n + bias);
        let new_accumulated = (warm.lambda_n + delta_lambda).max(0.0);
        let applied = new_accumulated - warm.lambda_n;
        let updated_cache = WarmStartCache {
            lambda_n: new_accumulated,
            lambda_t1: warm.lambda_t1,
            lambda_t2: warm.lambda_t2,
        };
        (applied, updated_cache)
    }
    /// Compute friction impulses along the two tangent directions.
    ///
    /// Applies Coulomb friction cone clamping: √(λ_t1² + λ_t2²) ≤ μ·λ_n.
    ///
    /// # Parameters
    /// - `rel_vel_t1`, `rel_vel_t2`: relative velocities along the two tangent axes
    /// - `lambda_n`: the current (accumulated) normal impulse magnitude
    /// - `inv_mass_a`, `inv_mass_b`: inverse masses
    ///
    /// Returns `(delta_lambda_t1, delta_lambda_t2)`.
    pub fn solve_friction_impulse(
        rel_vel_t1: f64,
        rel_vel_t2: f64,
        lambda_n: f64,
        inv_mass_a: f64,
        inv_mass_b: f64,
    ) -> (f64, f64) {
        let inv_mass_sum = inv_mass_a + inv_mass_b;
        if inv_mass_sum < 1e-12 {
            return (0.0, 0.0);
        }
        let effective_mass = 1.0 / inv_mass_sum;
        let raw_t1 = effective_mass * (-rel_vel_t1);
        let raw_t2 = effective_mass * (-rel_vel_t2);
        let limit = lambda_n.abs();
        let mag = (raw_t1 * raw_t1 + raw_t2 * raw_t2).sqrt();
        if mag > limit && mag > 1e-12 {
            let scale = limit / mag;
            (raw_t1 * scale, raw_t2 * scale)
        } else {
            (raw_t1, raw_t2)
        }
    }
}
/// Automatically adjusts the warm-start scale factor based on recent quality.
///
/// If warm starting is producing good guesses (low error, high hit rate), the
/// scale factor moves toward 1.0. If quality degrades, it shrinks toward a
/// minimum floor to avoid destabilizing the solver.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AdaptiveWarmStartScaler {
    /// Current scale factor applied to warm-start impulses.
    pub scale: f64,
    /// Minimum allowed scale factor.
    pub min_scale: f64,
    /// Maximum allowed scale factor (typically 1.0).
    pub max_scale: f64,
    /// How quickly the scale adapts (learning rate).
    pub adaptation_rate: f64,
    /// Quality threshold below which the scale decreases.
    pub quality_threshold: f64,
}
#[allow(dead_code)]
impl AdaptiveWarmStartScaler {
    /// Create a new adaptive scaler with default parameters.
    pub fn new() -> Self {
        Self {
            scale: 0.85,
            min_scale: 0.3,
            max_scale: 1.0,
            adaptation_rate: 0.1,
            quality_threshold: 0.5,
        }
    }
    /// Create with custom parameters.
    pub fn with_params(
        initial_scale: f64,
        min_scale: f64,
        max_scale: f64,
        adaptation_rate: f64,
        quality_threshold: f64,
    ) -> Self {
        Self {
            scale: initial_scale.clamp(min_scale, max_scale),
            min_scale,
            max_scale,
            adaptation_rate: adaptation_rate.clamp(0.0, 1.0),
            quality_threshold,
        }
    }
    /// Update the scale factor based on the latest frame's quality metrics.
    ///
    /// Returns the new scale factor.
    pub fn update(&mut self, metrics: &WarmStartQualityMetrics) -> f64 {
        let q = metrics.quality_score();
        if q >= self.quality_threshold {
            self.scale += self.adaptation_rate * (self.max_scale - self.scale);
        } else {
            self.scale -= self.adaptation_rate * (self.scale - self.min_scale);
        }
        self.scale = self.scale.clamp(self.min_scale, self.max_scale);
        self.scale
    }
    /// Apply the current scale to a warm-start cache entry.
    pub fn apply(&self, cache: &WarmStartCache) -> WarmStartCache {
        cache.scale(self.scale)
    }
}
/// A contact point described by its world-space position and normal.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ContactPoint3D {
    /// World-space contact position.
    pub position: [f64; 3],
    /// Outward contact normal (unit length).
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
}
#[allow(dead_code)]
impl ContactPoint3D {
    /// Create a contact point.
    pub fn new(position: [f64; 3], normal: [f64; 3], depth: f64) -> Self {
        Self {
            position,
            normal,
            depth,
        }
    }
}
/// Manages warm-start data for all active contact islands.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct IslandWarmStartManager {
    /// Map from island ID to its cached warm-start data.
    pub islands: HashMap<u32, IslandWarmStartData>,
    /// Maximum age before an island's data is discarded.
    pub max_age: u32,
}
#[allow(dead_code)]
impl IslandWarmStartManager {
    /// Create a new manager with the given max age.
    pub fn new(max_age: u32) -> Self {
        Self {
            islands: HashMap::new(),
            max_age,
        }
    }
    /// Insert or replace warm-start data for an island.
    pub fn store_island(&mut self, data: IslandWarmStartData) {
        self.islands.insert(data.island_id.0, data);
    }
    /// Retrieve (immutable) warm-start data for an island, if present.
    pub fn get_island(&self, id: IslandId) -> Option<&IslandWarmStartData> {
        self.islands.get(&id.0)
    }
    /// Age all islands by one frame and evict those exceeding `max_age`.
    pub fn advance_frame(&mut self) {
        let max_age = self.max_age;
        for data in self.islands.values_mut() {
            data.increment_age();
        }
        self.islands.retain(|_, data| data.age <= max_age);
    }
    /// Apply alpha-decay to all island impulse caches.
    pub fn decay_all(&mut self, alpha: f64) {
        let keep = (1.0 - alpha).max(0.0);
        for data in self.islands.values_mut() {
            data.decay_all(keep);
        }
    }
    /// Number of active islands.
    pub fn island_count(&self) -> usize {
        self.islands.len()
    }
}
/// Per-frame impulse aging using a linear alpha factor.
///
/// Each frame the stored impulse is scaled by `(1 - alpha)`, where alpha is
/// in \[0, 1\].  Alpha = 0 means no decay (full warm-start every frame);
/// alpha = 1 means the impulse is zeroed each frame (no warm-start).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AlphaImpulseAging {
    /// Decay factor in \[0, 1\].
    pub alpha: f64,
}
#[allow(dead_code)]
impl AlphaImpulseAging {
    /// Create a new alpha-based impulse aging policy.
    pub fn new(alpha: f64) -> Self {
        Self {
            alpha: alpha.clamp(0.0, 1.0),
        }
    }
    /// Apply aging to a single impulse triple.
    pub fn age(&self, impulse: [f64; 3]) -> [f64; 3] {
        let keep = 1.0 - self.alpha;
        [impulse[0] * keep, impulse[1] * keep, impulse[2] * keep]
    }
    /// Apply aging in-place to every entry in a `WarmStartMap`.
    pub fn age_map(&self, map: &mut WarmStartMap) {
        let keep = 1.0 - self.alpha;
        for v in map.cache.values_mut() {
            v.lambda_n *= keep;
            v.lambda_t1 *= keep;
            v.lambda_t2 *= keep;
        }
    }
}
/// Running statistics for a warm-start cache.
#[allow(dead_code)]
pub struct CacheStatistics {
    /// Total number of cache lookups.
    pub total_lookups: u64,
    /// Number of successful hits (entry found and applied).
    pub hits: u64,
    /// Number of misses (no entry found).
    pub misses: u64,
    /// Number of entries pruned (below threshold or too old).
    pub prunes: u64,
    /// Cumulative sum of hit impulse magnitudes (for mean computation).
    pub impulse_sum: f64,
    /// Cumulative sum of squared hit impulse magnitudes (for RMS).
    pub impulse_sq_sum: f64,
    /// Peak impulse magnitude seen.
    pub peak_impulse: f64,
}
#[allow(dead_code)]
impl CacheStatistics {
    /// Create new empty statistics.
    pub fn new() -> Self {
        Self {
            total_lookups: 0,
            hits: 0,
            misses: 0,
            prunes: 0,
            impulse_sum: 0.0,
            impulse_sq_sum: 0.0,
            peak_impulse: 0.0,
        }
    }
    /// Record a cache hit with the applied impulse magnitude.
    pub fn record_hit(&mut self, impulse_magnitude: f64) {
        self.total_lookups += 1;
        self.hits += 1;
        self.impulse_sum += impulse_magnitude;
        self.impulse_sq_sum += impulse_magnitude * impulse_magnitude;
        if impulse_magnitude > self.peak_impulse {
            self.peak_impulse = impulse_magnitude;
        }
    }
    /// Record a cache miss.
    pub fn record_miss(&mut self) {
        self.total_lookups += 1;
        self.misses += 1;
    }
    /// Record a prune event.
    pub fn record_prune(&mut self) {
        self.prunes += 1;
    }
    /// Hit rate in \[0, 1\].
    pub fn hit_rate(&self) -> f64 {
        if self.total_lookups == 0 {
            0.0
        } else {
            self.hits as f64 / self.total_lookups as f64
        }
    }
    /// Mean impulse magnitude of hits (0 if no hits).
    pub fn mean_impulse(&self) -> f64 {
        if self.hits == 0 {
            0.0
        } else {
            self.impulse_sum / self.hits as f64
        }
    }
    /// RMS impulse magnitude of hits.
    pub fn rms_impulse(&self) -> f64 {
        if self.hits == 0 {
            0.0
        } else {
            (self.impulse_sq_sum / self.hits as f64).sqrt()
        }
    }
    /// Reset all counters.
    pub fn reset(&mut self) {
        *self = Self::new();
    }
}
/// Combines alpha decay and manifold-quality gating.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct QualityBasedAging {
    /// Frame-to-frame alpha decay (see [`AlphaImpulseAging`]).
    pub alpha: f64,
    /// Minimum quality below which the impulse is zeroed entirely.
    pub quality_floor: f64,
}
#[allow(dead_code)]
impl QualityBasedAging {
    /// Create a quality-based aging policy.
    pub fn new(alpha: f64, quality_floor: f64) -> Self {
        Self {
            alpha: alpha.clamp(0.0, 1.0),
            quality_floor: quality_floor.clamp(0.0, 1.0),
        }
    }
    /// Compute the effective scaling factor for the given manifold quality.
    ///
    /// If `quality < quality_floor`, returns 0 (impulse discarded).
    /// Otherwise returns `quality * (1 - alpha)`.
    pub fn effective_scale(&self, quality: f64) -> f64 {
        if quality < self.quality_floor {
            0.0
        } else {
            quality * (1.0 - self.alpha)
        }
    }
    /// Scale a single impulse by the effective factor derived from `quality`.
    pub fn age_impulse(&self, impulse: [f64; 3], quality: f64) -> [f64; 3] {
        let s = self.effective_scale(quality);
        [impulse[0] * s, impulse[1] * s, impulse[2] * s]
    }
}
/// Strategy for how warm-start impulses are applied from the previous frame.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[allow(dead_code)]
pub enum WarmStartStrategy {
    /// No warm starting — every frame starts from zero impulses.
    None,
    /// Full warm starting — previous impulses are re-applied at 100 %.
    #[default]
    Full,
    /// Decayed warm starting — previous impulses are scaled by a decay factor
    /// each frame (e.g. 0.85 means 85 % of the old impulse is kept).
    Decayed(f64),
    /// Adaptive warm starting — the scale factor is adjusted automatically
    /// based on convergence quality metrics.
    Adaptive,
}
/// Quality score for a contact manifold in \[0, 1\].
///
/// 1.0 = perfect match from last frame (no decay), 0.0 = completely new
/// contact (full decay).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ManifoldQuality(pub f64);
#[allow(dead_code)]
impl ManifoldQuality {
    /// Clamp the raw quality value to \[0, 1\].
    pub fn new(raw: f64) -> Self {
        Self(raw.clamp(0.0, 1.0))
    }
    /// Apply quality-based scaling to a stored impulse.
    ///
    /// When quality = 1 the impulse is unchanged; when quality = 0 the
    /// impulse is zeroed.
    pub fn scale_impulse(&self, impulse: [f64; 3]) -> [f64; 3] {
        let q = self.0;
        [impulse[0] * q, impulse[1] * q, impulse[2] * q]
    }
}
/// Quality classification for a warm-start impulse entry.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub enum WarmStartQuality {
    /// Data is recent (age == 0 last frame).
    Fresh,
    /// Data has been around for several frames but is within max_age.
    Aged(f64),
    /// Data is too old to trust.
    Stale,
}
/// A map from contact pair keys to their warm-start cache.
///
/// Keys are `(body_a_id, body_b_id)` with `body_a_id < body_b_id` to ensure
/// consistent ordering regardless of which body is "A" or "B" in the constraint.
#[allow(dead_code)]
pub struct WarmStartMap {
    pub(super) cache: HashMap<(u64, u64), WarmStartCache>,
}
#[allow(dead_code)]
impl WarmStartMap {
    /// Create a new, empty warm-start map.
    pub fn new() -> Self {
        Self {
            cache: HashMap::new(),
        }
    }
    /// Create with a pre-allocated capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            cache: HashMap::with_capacity(capacity),
        }
    }
    /// Retrieve cached warm-start data for a body pair, or return defaults if not found.
    pub fn get_or_default(&self, a: u64, b: u64) -> WarmStartCache {
        let k = Self::key(a, b);
        self.cache.get(&k).cloned().unwrap_or_default()
    }
    /// Retrieve cached warm-start data for a body pair.
    pub fn get(&self, a: u64, b: u64) -> Option<&WarmStartCache> {
        let k = Self::key(a, b);
        self.cache.get(&k)
    }
    /// Store warm-start data for a body pair.
    pub fn store(&mut self, a: u64, b: u64, ws: WarmStartCache) {
        let k = Self::key(a, b);
        self.cache.insert(k, ws);
    }
    /// Remove a specific pair's cache entry.
    pub fn remove(&mut self, a: u64, b: u64) -> Option<WarmStartCache> {
        let k = Self::key(a, b);
        self.cache.remove(&k)
    }
    /// Clear all cached warm-start data (call at simulation reset).
    pub fn clear(&mut self) {
        self.cache.clear();
    }
    /// Number of cached entries.
    pub fn len(&self) -> usize {
        self.cache.len()
    }
    /// Whether the map is empty.
    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }
    /// Compute the canonical key for a body pair, ensuring `a < b` ordering.
    pub fn key(a: u64, b: u64) -> (u64, u64) {
        if a <= b { (a, b) } else { (b, a) }
    }
    /// Apply a uniform scale factor to all cached impulses.
    pub fn scale_all(&mut self, factor: f64) {
        for entry in self.cache.values_mut() {
            entry.lambda_n *= factor;
            entry.lambda_t1 *= factor;
            entry.lambda_t2 *= factor;
        }
    }
    /// Remove entries whose impulse magnitude is below a threshold.
    pub fn prune(&mut self, threshold: f64) {
        self.cache.retain(|_, v| v.magnitude() >= threshold);
    }
    /// Apply a strategy-based transformation before warm-starting a new frame.
    pub fn apply_strategy(&mut self, strategy: WarmStartStrategy) {
        match strategy {
            WarmStartStrategy::None => self.clear(),
            WarmStartStrategy::Full => {}
            WarmStartStrategy::Decayed(factor) => self.scale_all(factor),
            WarmStartStrategy::Adaptive => {}
        }
    }
    /// Iterate over all cached pairs and their impulse data.
    pub fn iter(&self) -> impl Iterator<Item = (&(u64, u64), &WarmStartCache)> {
        self.cache.iter()
    }
    /// Compute aggregate statistics: total impulse magnitude and max single entry.
    pub fn stats(&self) -> (f64, f64) {
        let mut total = 0.0_f64;
        let mut max_mag = 0.0_f64;
        for entry in self.cache.values() {
            let m = entry.magnitude();
            total += m;
            if m > max_mag {
                max_mag = m;
            }
        }
        (total, max_mag)
    }
}
/// Blend two warm-start caches according to a blending mode.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub enum BlendMode {
    /// Linearly interpolate between prev and current.
    Lerp(f64),
    /// Take the larger-magnitude component-wise.
    Max,
    /// Take the smaller-magnitude component-wise.
    Min,
    /// Average both caches.
    Average,
}
/// Selection policy for choosing which warm-start entry to use when multiple
/// candidates exist.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub enum CacheSelectionPolicy {
    /// Use the most recently updated entry.
    MostRecent,
    /// Use the entry with the largest impulse magnitude.
    LargestMagnitude,
    /// Use the entry with the smallest impulse magnitude.
    SmallestMagnitude,
    /// Use the entry with quality score closest to a target value.
    ClosestQuality(f64),
}
