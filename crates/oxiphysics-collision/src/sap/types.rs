//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::{HashMap, HashSet};

use super::functions::{aabb3_overlaps, aabb3_point_dist_sq, axis_is_sorted};

/// An endpoint on one axis used by the Sweep-and-Prune broadphase.
#[derive(Debug, Clone)]
pub struct SapEndpoint {
    /// The coordinate value of this endpoint.
    pub value: f64,
    /// Which object this endpoint belongs to.
    pub object_id: u64,
    /// `true` for the minimum (left) endpoint, `false` for the maximum (right).
    pub is_min: bool,
}
/// Extended SAP that tracks statistics.
pub struct StatTrackingSap {
    /// Inner SAP structure.
    pub inner: IncrementalSap,
    /// Statistics from the last query.
    pub stats: SapStats,
}
impl StatTrackingSap {
    /// Create a new statistics-tracking SAP.
    pub fn new() -> Self {
        Self {
            inner: IncrementalSap::new(),
            stats: SapStats::default(),
        }
    }
    /// Insert a body.
    pub fn insert(&mut self, id: u32, aabb: Aabb3) {
        self.inner.insert(id, aabb);
    }
    /// Remove a body.
    pub fn remove(&mut self, id: u32) {
        self.inner.remove(id);
    }
    /// Update a body's AABB.
    pub fn update(&mut self, id: u32, aabb: Aabb3) {
        self.inner.update(id, aabb);
    }
    /// Query all overlapping pairs, recording statistics.
    pub fn query_pairs(&mut self) -> Vec<(u32, u32)> {
        let n_endpoints = self.inner.endpoints_x.len();
        let pairs = self.inner.query_pairs();
        self.stats = SapStats {
            pair_count: pairs.len(),
            sweep_count: n_endpoints,
            body_count: self.inner.body_count(),
        };
        pairs
    }
}
/// An axis-aligned bounding box stored inside the SAP structure.
#[derive(Debug, Clone)]
pub struct SapObject {
    /// Unique identifier for this object.
    pub id: u64,
    /// Minimum corner `[x, y, z]`.
    pub min: [f64; 3],
    /// Maximum corner `[x, y, z]`.
    pub max: [f64; 3],
}
/// Result of a single broadphase step.
#[derive(Debug, Clone)]
pub struct BroadphaseResult {
    /// All currently overlapping pairs.
    pub pairs: Vec<(u32, u32)>,
    /// New pairs this frame (begin events).
    pub new_pairs: Vec<(u32, u32)>,
    /// Pairs that ended this frame (end events).
    pub lost_pairs: Vec<(u32, u32)>,
    /// Statistics for this step.
    pub stats: SapStats,
}
/// An endpoint on one axis used by the incremental multi-axis SAP.
#[derive(Debug, Clone)]
pub struct SapEndpointU32 {
    /// The coordinate value of this endpoint.
    pub value: f64,
    /// Which body this endpoint belongs to.
    pub body_id: u32,
    /// `true` for the minimum (left) endpoint, `false` for the maximum (right).
    pub is_min: bool,
}
/// Uniform-grid broadphase that hashes objects into fixed-size cells.
pub struct GridBroadphase {
    /// Side length of each cubic cell.
    pub cell_size: f64,
    /// Map from cell coordinate to the list of object ids in that cell.
    pub cells: HashMap<(i32, i32, i32), Vec<u64>>,
}
impl GridBroadphase {
    /// Create a new grid with the given cell size.
    pub fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            cells: HashMap::new(),
        }
    }
    /// Insert an object into every grid cell that its AABB overlaps.
    pub fn insert(&mut self, id: u64, min: [f64; 3], max: [f64; 3]) {
        let min_cx = self.cell_coord(min[0]);
        let min_cy = self.cell_coord(min[1]);
        let min_cz = self.cell_coord(min[2]);
        let max_cx = self.cell_coord(max[0]);
        let max_cy = self.cell_coord(max[1]);
        let max_cz = self.cell_coord(max[2]);
        for cx in min_cx..=max_cx {
            for cy in min_cy..=max_cy {
                for cz in min_cz..=max_cz {
                    self.cells.entry((cx, cy, cz)).or_default().push(id);
                }
            }
        }
    }
    /// Convert a world-space coordinate to its cell index.
    #[inline]
    pub fn cell_coord(&self, pos: f64) -> i32 {
        (pos / self.cell_size).floor() as i32
    }
    /// Return all object-id pairs that share at least one cell (conservative).
    pub fn query_potential_pairs(&self) -> Vec<(u64, u64)> {
        let mut pairs: Vec<(u64, u64)> = Vec::new();
        for ids in self.cells.values() {
            for i in 0..ids.len() {
                for j in (i + 1)..ids.len() {
                    let (lo, hi) = if ids[i] < ids[j] {
                        (ids[i], ids[j])
                    } else {
                        (ids[j], ids[i])
                    };
                    pairs.push((lo, hi));
                }
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        pairs
    }
    /// Remove all objects from all cells.
    pub fn clear(&mut self) {
        self.cells.clear();
    }
}
impl GridBroadphase {
    /// Return all object ids in the cell that contains world-space point `p`.
    pub fn query_point(&self, p: [f64; 3]) -> Vec<u64> {
        let key = (
            self.cell_coord(p[0]),
            self.cell_coord(p[1]),
            self.cell_coord(p[2]),
        );
        self.cells.get(&key).cloned().unwrap_or_default()
    }
    /// Return total number of object-cell entries (not unique objects).
    pub fn total_entries(&self) -> usize {
        self.cells.values().map(|v| v.len()).sum()
    }
    /// Return the number of occupied cells.
    pub fn occupied_cells(&self) -> usize {
        self.cells.len()
    }
    /// Remove a specific object id from all cells it occupies.
    ///
    /// Scans all cells — `O(cells * max_per_cell)`.  For high-density scenes
    /// prefer `clear` + bulk re-insert.
    pub fn remove(&mut self, id: u64) {
        for ids in self.cells.values_mut() {
            ids.retain(|&x| x != id);
        }
        self.cells.retain(|_, v| !v.is_empty());
    }
    /// Return a conservative over-approximation of all objects that could
    /// overlap `aabb` by checking every cell the AABB touches.
    pub fn query_aabb(&self, min: [f64; 3], max: [f64; 3]) -> Vec<u64> {
        let min_cx = self.cell_coord(min[0]);
        let min_cy = self.cell_coord(min[1]);
        let min_cz = self.cell_coord(min[2]);
        let max_cx = self.cell_coord(max[0]);
        let max_cy = self.cell_coord(max[1]);
        let max_cz = self.cell_coord(max[2]);
        let mut result = Vec::new();
        for cx in min_cx..=max_cx {
            for cy in min_cy..=max_cy {
                for cz in min_cz..=max_cz {
                    if let Some(ids) = self.cells.get(&(cx, cy, cz)) {
                        for &id in ids {
                            if !result.contains(&id) {
                                result.push(id);
                            }
                        }
                    }
                }
            }
        }
        result.sort_unstable();
        result
    }
}
/// A sorted endpoint array for one axis used by the incremental SAP.
///
/// Maintains a sorted list of `SapEndpointU32` entries and provides
/// insertion-sort for small moves.
#[derive(Debug, Clone, Default)]
pub struct SapAxis {
    /// Sorted endpoints on this axis.
    pub endpoints: Vec<SapEndpointU32>,
}
impl SapAxis {
    /// Create an empty axis.
    pub fn new() -> Self {
        Self {
            endpoints: Vec::new(),
        }
    }
    /// Insert min and max endpoints for `body_id` with the given values.
    pub fn insert(&mut self, body_id: u32, min_val: f64, max_val: f64) {
        self.endpoints.push(SapEndpointU32 {
            value: min_val,
            body_id,
            is_min: true,
        });
        self.endpoints.push(SapEndpointU32 {
            value: max_val,
            body_id,
            is_min: false,
        });
        let n = self.endpoints.len();
        for start in [n - 2, n - 1] {
            let mut i = start;
            while i > 0 && self.endpoints[i - 1].value > self.endpoints[i].value {
                self.endpoints.swap(i - 1, i);
                i -= 1;
            }
        }
    }
    /// Remove all endpoints belonging to `body_id`.
    pub fn remove(&mut self, body_id: u32) {
        self.endpoints.retain(|e| e.body_id != body_id);
    }
    /// Update the endpoints for `body_id` using insertion sort (efficient for small moves).
    pub fn update(&mut self, body_id: u32, min_val: f64, max_val: f64) {
        for ep in self.endpoints.iter_mut() {
            if ep.body_id == body_id {
                ep.value = if ep.is_min { min_val } else { max_val };
            }
        }
        let n = self.endpoints.len();
        for i in 1..n {
            let mut j = i;
            while j > 0 && self.endpoints[j - 1].value > self.endpoints[j].value {
                self.endpoints.swap(j - 1, j);
                j -= 1;
            }
        }
    }
    /// Sweep the sorted endpoints and return overlapping pairs on this axis.
    pub fn overlapping_pairs(&self) -> HashSet<(u32, u32)> {
        let mut pairs = HashSet::new();
        let mut active: Vec<u32> = Vec::new();
        for ep in &self.endpoints {
            if ep.is_min {
                for &aid in &active {
                    let key = if ep.body_id < aid {
                        (ep.body_id, aid)
                    } else {
                        (aid, ep.body_id)
                    };
                    pairs.insert(key);
                }
                active.push(ep.body_id);
            } else {
                active.retain(|&id| id != ep.body_id);
            }
        }
        pairs
    }
}
impl SapAxis {
    /// Return `true` if the endpoint list is sorted.
    pub fn is_sorted(&self) -> bool {
        axis_is_sorted(&self.endpoints)
    }
    /// Return the number of distinct bodies tracked on this axis.
    pub fn body_count(&self) -> usize {
        self.endpoints.iter().filter(|e| e.is_min).count()
    }
    /// Return the minimum value across all endpoints.
    pub fn min_value(&self) -> Option<f64> {
        self.endpoints.first().map(|e| e.value)
    }
    /// Return the maximum value across all endpoints.
    pub fn max_value(&self) -> Option<f64> {
        self.endpoints.last().map(|e| e.value)
    }
    /// Return all body ids currently tracked on this axis.
    pub fn tracked_bodies(&self) -> Vec<u32> {
        let mut ids: Vec<u32> = self
            .endpoints
            .iter()
            .filter(|e| e.is_min)
            .map(|e| e.body_id)
            .collect();
        ids.sort_unstable();
        ids.dedup();
        ids
    }
}
/// A snapshot of an `IncrementalSap` that can be restored cheaply.
#[derive(Debug, Clone)]
pub struct SapSnapshot {
    pub(super) aabbs: HashMap<u32, Aabb3>,
}
impl SapSnapshot {
    /// Capture the current AABB map from `sap`.
    pub fn capture(sap: &IncrementalSap) -> Self {
        Self {
            aabbs: sap.aabbs.clone(),
        }
    }
    /// Restore `sap` to the state recorded in this snapshot.
    ///
    /// Bodies added after the snapshot are removed; bodies removed after the
    /// snapshot are re-added; changed AABBs are reverted.
    pub fn restore(self, sap: &mut IncrementalSap) {
        let current_ids: Vec<u32> = sap.aabbs.keys().copied().collect();
        for id in &current_ids {
            if !self.aabbs.contains_key(id) {
                sap.remove(*id);
            }
        }
        for (id, aabb) in self.aabbs {
            if sap.aabbs.contains_key(&id) {
                sap.update(id, aabb);
            } else {
                sap.insert(id, aabb);
            }
        }
    }
}
/// An overlap event: a pair has just started or stopped overlapping.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OverlapEvent {
    /// Bodies `(a, b)` began overlapping this frame.
    Begin(u32, u32),
    /// Bodies `(a, b)` stopped overlapping this frame.
    End(u32, u32),
}
/// A simple AABB for the incremental SAP (raw f64 arrays).
#[derive(Debug, Clone)]
pub struct Aabb3 {
    /// Minimum corner `[x, y, z]`.
    pub min: [f64; 3],
    /// Maximum corner `[x, y, z]`.
    pub max: [f64; 3],
}
/// Statistics collected during a SAP broadphase query.
#[derive(Debug, Clone, Default)]
pub struct SapStats {
    /// Number of overlapping pairs found in the last query.
    pub pair_count: usize,
    /// Number of endpoint comparisons made during the last sweep.
    pub sweep_count: usize,
    /// Number of active bodies tracked.
    pub body_count: usize,
}
/// SAP that tracks pair-begin / pair-end events across frames.
///
/// Call `step_pairs` each frame; it returns only the *delta* — new overlaps
/// and vanished overlaps — rather than the full pair list.
pub struct EventDrivenSap {
    pub(super) inner: IncrementalSap,
    /// Pairs that were overlapping in the *previous* frame.
    pub(super) prev_pairs: HashSet<(u32, u32)>,
}
impl EventDrivenSap {
    /// Create a new event-driven SAP.
    pub fn new() -> Self {
        Self {
            inner: IncrementalSap::new(),
            prev_pairs: HashSet::new(),
        }
    }
    /// Insert a body.
    pub fn insert(&mut self, id: u32, aabb: Aabb3) {
        self.inner.insert(id, aabb);
    }
    /// Remove a body.  Any pairs it was part of will appear as `End` events
    /// on the next `step_pairs` call.
    pub fn remove(&mut self, id: u32) {
        self.inner.remove(id);
    }
    /// Update a body's AABB.
    pub fn update(&mut self, id: u32, aabb: Aabb3) {
        self.inner.update(id, aabb);
    }
    /// Compute overlap events for the current frame.
    ///
    /// Returns `(events, current_pairs)` where `events` lists only the
    /// begin/end transitions relative to the previous frame.
    pub fn step_pairs(&mut self) -> (Vec<OverlapEvent>, Vec<(u32, u32)>) {
        let current_vec = self.inner.query_pairs();
        let current: HashSet<(u32, u32)> = current_vec.iter().copied().collect();
        let mut events = Vec::new();
        for &pair in &current {
            if !self.prev_pairs.contains(&pair) {
                events.push(OverlapEvent::Begin(pair.0, pair.1));
            }
        }
        for &pair in &self.prev_pairs {
            if !current.contains(&pair) {
                events.push(OverlapEvent::End(pair.0, pair.1));
            }
        }
        self.prev_pairs = current;
        events.sort_by_key(|e| match e {
            OverlapEvent::Begin(a, b) | OverlapEvent::End(a, b) => (*a, *b),
        });
        (events, current_vec)
    }
    /// Number of currently tracked bodies.
    pub fn body_count(&self) -> usize {
        self.inner.body_count()
    }
}
/// Sweep-and-Prune broadphase that operates on the X axis and filters on Y/Z.
pub struct SweepAndPrune {
    /// All registered objects.
    pub objects: Vec<SapObject>,
    /// Sorted endpoints on the X axis.
    pub sorted_x: Vec<SapEndpoint>,
}
impl SweepAndPrune {
    /// Create an empty `SweepAndPrune` instance.
    pub fn new() -> Self {
        Self {
            objects: Vec::new(),
            sorted_x: Vec::new(),
        }
    }
    /// Add an object with the given AABB.
    pub fn add_object(&mut self, id: u64, min: [f64; 3], max: [f64; 3]) {
        self.objects.push(SapObject { id, min, max });
        self.sorted_x.push(SapEndpoint {
            value: min[0],
            object_id: id,
            is_min: true,
        });
        self.sorted_x.push(SapEndpoint {
            value: max[0],
            object_id: id,
            is_min: false,
        });
        self.sort_axis();
    }
    /// Remove the object with the given `id`.
    pub fn remove_object(&mut self, id: u64) {
        self.objects.retain(|o| o.id != id);
        self.sorted_x.retain(|e| e.object_id != id);
    }
    /// Update (replace) the AABB for an existing object.
    pub fn update_object(&mut self, id: u64, min: [f64; 3], max: [f64; 3]) {
        if let Some(obj) = self.objects.iter_mut().find(|o| o.id == id) {
            obj.min = min;
            obj.max = max;
        }
        for ep in self.sorted_x.iter_mut().filter(|e| e.object_id == id) {
            ep.value = if ep.is_min { min[0] } else { max[0] };
        }
        self.sort_axis();
    }
    /// Sort the X-axis endpoint list by value.
    pub fn sort_axis(&mut self) {
        self.sorted_x.sort_by(|a, b| {
            a.value
                .partial_cmp(&b.value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    /// Return all overlapping pairs.
    ///
    /// Sweeps the sorted X endpoints to collect candidate pairs that overlap on
    /// X, then filters each candidate by checking Y and Z overlap as well.
    /// The returned list is deduplicated and has `(min_id, max_id)` ordering.
    pub fn query_overlapping_pairs(&mut self) -> Vec<(u64, u64)> {
        self.sort_axis();
        let obj_map: HashMap<u64, &SapObject> = self.objects.iter().map(|o| (o.id, o)).collect();
        let mut pairs: Vec<(u64, u64)> = Vec::new();
        let mut active: Vec<u64> = Vec::new();
        for ep in &self.sorted_x {
            if ep.is_min {
                if let Some(obj_a) = obj_map.get(&ep.object_id) {
                    for &active_id in &active {
                        if let Some(obj_b) = obj_map.get(&active_id)
                            && Self::overlaps_on_axis(
                                obj_a.min[1],
                                obj_a.max[1],
                                obj_b.min[1],
                                obj_b.max[1],
                            )
                            && Self::overlaps_on_axis(
                                obj_a.min[2],
                                obj_a.max[2],
                                obj_b.min[2],
                                obj_b.max[2],
                            )
                        {
                            let (lo, hi) = if ep.object_id < active_id {
                                (ep.object_id, active_id)
                            } else {
                                (active_id, ep.object_id)
                            };
                            pairs.push((lo, hi));
                        }
                    }
                }
                active.push(ep.object_id);
            } else {
                active.retain(|&id| id != ep.object_id);
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        pairs
    }
    /// Test whether two intervals `[min_a, max_a]` and `[min_b, max_b]` overlap.
    #[inline]
    pub fn overlaps_on_axis(min_a: f64, max_a: f64, min_b: f64, max_b: f64) -> bool {
        min_a <= max_b && min_b <= max_a
    }
    /// Return the number of registered objects.
    pub fn object_count(&self) -> usize {
        self.objects.len()
    }
}
impl SweepAndPrune {
    /// Batch-insert multiple objects at once, sorting the axis list only once.
    ///
    /// More efficient than calling `add_object` repeatedly when adding many
    /// objects at the same time.
    pub fn add_object_batch(&mut self, objects: &[(u64, [f64; 3], [f64; 3])]) {
        for &(id, min, max) in objects {
            self.objects.push(SapObject { id, min, max });
            self.sorted_x.push(SapEndpoint {
                value: min[0],
                object_id: id,
                is_min: true,
            });
            self.sorted_x.push(SapEndpoint {
                value: max[0],
                object_id: id,
                is_min: false,
            });
        }
        self.sort_axis();
    }
    /// Compute the variance of AABB centres along each axis.
    ///
    /// Returns `[var_x, var_y, var_z]`.  The axis with the highest variance is
    /// the best choice for the SAP sweep direction (maximises early-out
    /// opportunities).
    pub fn compute_axis_variance(&self) -> [f64; 3] {
        let n = self.objects.len();
        if n == 0 {
            return [0.0; 3];
        }
        let mut sum = [0.0_f64; 3];
        let mut sum_sq = [0.0_f64; 3];
        for obj in &self.objects {
            for ax in 0..3 {
                let c = (obj.min[ax] + obj.max[ax]) * 0.5;
                sum[ax] += c;
                sum_sq[ax] += c * c;
            }
        }
        let nf = n as f64;
        let mut var = [0.0_f64; 3];
        for ax in 0..3 {
            let mean = sum[ax] / nf;
            var[ax] = (sum_sq[ax] / nf) - mean * mean;
        }
        var
    }
    /// Reorder the SAP to sweep along the axis with the highest variance.
    ///
    /// Updates `self.axis` (if it were a stored field; here we rebuild
    /// `sorted_x` to reflect the chosen axis) and re-sorts the endpoint list.
    ///
    /// Returns the chosen axis index (`0`=X, `1`=Y, `2`=Z).
    pub fn reorder_axes(&mut self) -> usize {
        let var = self.compute_axis_variance();
        let best_axis = if var[0] >= var[1] && var[0] >= var[2] {
            0
        } else if var[1] >= var[2] {
            1
        } else {
            2
        };
        self.sorted_x.clear();
        for obj in &self.objects {
            self.sorted_x.push(SapEndpoint {
                value: obj.min[best_axis],
                object_id: obj.id,
                is_min: true,
            });
            self.sorted_x.push(SapEndpoint {
                value: obj.max[best_axis],
                object_id: obj.id,
                is_min: false,
            });
        }
        self.sort_axis();
        best_axis
    }
}
/// Incremental multi-axis Sweep-and-Prune broadphase.
///
/// Maintains sorted endpoint lists on X, Y, and Z axes and intersects
/// overlap sets from all three to produce the final pair set.
pub struct IncrementalSap {
    /// Sorted endpoints on the X axis.
    pub endpoints_x: Vec<SapEndpointU32>,
    /// Sorted endpoints on the Y axis.
    pub endpoints_y: Vec<SapEndpointU32>,
    /// Sorted endpoints on the Z axis.
    pub endpoints_z: Vec<SapEndpointU32>,
    /// Currently tracked AABBs keyed by body id.
    pub aabbs: HashMap<u32, Aabb3>,
    /// Active overlapping pairs from the last query.
    pub active_pairs: HashSet<(u32, u32)>,
}
impl IncrementalSap {
    /// Create an empty incremental SAP.
    pub fn new() -> Self {
        Self {
            endpoints_x: Vec::new(),
            endpoints_y: Vec::new(),
            endpoints_z: Vec::new(),
            aabbs: HashMap::new(),
            active_pairs: HashSet::new(),
        }
    }
    /// Insert a body with the given AABB.
    pub fn insert(&mut self, body_id: u32, aabb: Aabb3) {
        self.endpoints_x.push(SapEndpointU32 {
            value: aabb.min[0],
            body_id,
            is_min: true,
        });
        self.endpoints_x.push(SapEndpointU32 {
            value: aabb.max[0],
            body_id,
            is_min: false,
        });
        self.endpoints_y.push(SapEndpointU32 {
            value: aabb.min[1],
            body_id,
            is_min: true,
        });
        self.endpoints_y.push(SapEndpointU32 {
            value: aabb.max[1],
            body_id,
            is_min: false,
        });
        self.endpoints_z.push(SapEndpointU32 {
            value: aabb.min[2],
            body_id,
            is_min: true,
        });
        self.endpoints_z.push(SapEndpointU32 {
            value: aabb.max[2],
            body_id,
            is_min: false,
        });
        self.aabbs.insert(body_id, aabb);
    }
    /// Remove a body from the SAP.
    pub fn remove(&mut self, body_id: u32) {
        self.endpoints_x.retain(|e| e.body_id != body_id);
        self.endpoints_y.retain(|e| e.body_id != body_id);
        self.endpoints_z.retain(|e| e.body_id != body_id);
        self.aabbs.remove(&body_id);
        self.active_pairs
            .retain(|&(a, b)| a != body_id && b != body_id);
    }
    /// Update the AABB for an existing body.
    pub fn update(&mut self, body_id: u32, new_aabb: Aabb3) {
        for ep in self.endpoints_x.iter_mut().filter(|e| e.body_id == body_id) {
            ep.value = if ep.is_min {
                new_aabb.min[0]
            } else {
                new_aabb.max[0]
            };
        }
        for ep in self.endpoints_y.iter_mut().filter(|e| e.body_id == body_id) {
            ep.value = if ep.is_min {
                new_aabb.min[1]
            } else {
                new_aabb.max[1]
            };
        }
        for ep in self.endpoints_z.iter_mut().filter(|e| e.body_id == body_id) {
            ep.value = if ep.is_min {
                new_aabb.min[2]
            } else {
                new_aabb.max[2]
            };
        }
        self.aabbs.insert(body_id, new_aabb);
    }
    /// Batch update multiple bodies at once, then re-sort once.
    pub fn batch_update(&mut self, updates: &[(u32, Aabb3)]) {
        for (body_id, new_aabb) in updates {
            for ep in self
                .endpoints_x
                .iter_mut()
                .filter(|e| e.body_id == *body_id)
            {
                ep.value = if ep.is_min {
                    new_aabb.min[0]
                } else {
                    new_aabb.max[0]
                };
            }
            for ep in self
                .endpoints_y
                .iter_mut()
                .filter(|e| e.body_id == *body_id)
            {
                ep.value = if ep.is_min {
                    new_aabb.min[1]
                } else {
                    new_aabb.max[1]
                };
            }
            for ep in self
                .endpoints_z
                .iter_mut()
                .filter(|e| e.body_id == *body_id)
            {
                ep.value = if ep.is_min {
                    new_aabb.min[2]
                } else {
                    new_aabb.max[2]
                };
            }
            self.aabbs.insert(*body_id, new_aabb.clone());
        }
    }
    /// Query all overlapping pairs by intersecting results from all three axes.
    pub fn query_pairs(&mut self) -> Vec<(u32, u32)> {
        let pairs_x = Self::sort_and_sweep_axis(&mut self.endpoints_x);
        let pairs_y = Self::sort_and_sweep_axis(&mut self.endpoints_y);
        let pairs_z = Self::sort_and_sweep_axis(&mut self.endpoints_z);
        let result: HashSet<(u32, u32)> = pairs_x
            .intersection(&pairs_y)
            .copied()
            .collect::<HashSet<_>>()
            .intersection(&pairs_z)
            .copied()
            .collect();
        self.active_pairs = result.clone();
        let mut pairs: Vec<(u32, u32)> = result.into_iter().collect();
        pairs.sort_unstable();
        pairs
    }
    /// Sort endpoints on one axis and sweep to find overlapping pairs.
    pub fn sort_and_sweep_axis(endpoints: &mut [SapEndpointU32]) -> HashSet<(u32, u32)> {
        endpoints.sort_by(|a, b| {
            a.value
                .partial_cmp(&b.value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut pairs = HashSet::new();
        let mut active: Vec<u32> = Vec::new();
        for ep in endpoints.iter() {
            if ep.is_min {
                for &aid in &active {
                    let (lo, hi) = if ep.body_id < aid {
                        (ep.body_id, aid)
                    } else {
                        (aid, ep.body_id)
                    };
                    pairs.insert((lo, hi));
                }
                active.push(ep.body_id);
            } else {
                active.retain(|&id| id != ep.body_id);
            }
        }
        pairs
    }
    /// Return the number of tracked bodies.
    pub fn body_count(&self) -> usize {
        self.aabbs.len()
    }
    /// Return the current set of active pairs (from last query).
    pub fn current_pairs(&self) -> &HashSet<(u32, u32)> {
        &self.active_pairs
    }
    /// Insert a new AABB for `id` and update the active pair set.
    pub fn insert_aabb(&mut self, id: u32, min: [f64; 3], max: [f64; 3]) {
        self.insert(id, Aabb3 { min, max });
        self.query_pairs();
    }
    /// Remove the AABB for `id` and clear its affected pairs.
    pub fn remove_aabb(&mut self, id: u32) {
        self.remove(id);
        self.active_pairs.retain(|&(a, b)| a != id && b != id);
    }
    /// Update the AABB for `id` using incremental re-sort (insertion sort for small moves).
    pub fn update_aabb(&mut self, id: u32, min: [f64; 3], max: [f64; 3]) {
        self.update(id, Aabb3 { min, max });
        self.query_pairs();
    }
    /// Return the current overlapping pair set.
    pub fn active_pairs(&self) -> &HashSet<(u32, u32)> {
        &self.active_pairs
    }
    /// Compute broadphase pairs between two disjoint sets using 3-axis SAP.
    ///
    /// Only pairs `(a, b)` where `a ∈ set_a` and `b ∈ set_b` are returned.
    pub fn bipartite_pairs(&self, set_a: &[u32], set_b: &[u32]) -> Vec<(u32, u32)> {
        let set_a_hs: HashSet<u32> = set_a.iter().copied().collect();
        let set_b_hs: HashSet<u32> = set_b.iter().copied().collect();
        let mut temp = IncrementalSap::new();
        for &id in set_a.iter().chain(set_b.iter()) {
            if let Some(aabb) = self.aabbs.get(&id) {
                temp.insert(id, aabb.clone());
            }
        }
        let all_pairs = temp.query_pairs();
        all_pairs
            .into_iter()
            .filter(|&(a, b)| {
                (set_a_hs.contains(&a) && set_b_hs.contains(&b))
                    || (set_a_hs.contains(&b) && set_b_hs.contains(&a))
            })
            .collect()
    }
}
impl IncrementalSap {
    /// Compute the variance of AABB centres on each axis, returning
    /// `[var_x, var_y, var_z]`.  The axis with highest variance is the best
    /// sweep direction.
    pub fn compute_axis_variance(&self) -> [f64; 3] {
        let n = self.aabbs.len();
        if n == 0 {
            return [0.0; 3];
        }
        let mut sum = [0.0_f64; 3];
        let mut sum_sq = [0.0_f64; 3];
        for aabb in self.aabbs.values() {
            for ax in 0..3 {
                let c = (aabb.min[ax] + aabb.max[ax]) * 0.5;
                sum[ax] += c;
                sum_sq[ax] += c * c;
            }
        }
        let nf = n as f64;
        let mut var = [0.0_f64; 3];
        for ax in 0..3 {
            let mean = sum[ax] / nf;
            var[ax] = (sum_sq[ax] / nf) - mean * mean;
        }
        var
    }
    /// Reorder all three endpoint arrays to sweep along the axis with the
    /// highest AABB-centre variance, then re-sort.
    ///
    /// Returns the chosen axis index (`0`=X, `1`=Y, `2`=Z).
    pub fn reorder_axes(&mut self) -> usize {
        let var = self.compute_axis_variance();
        let best = if var[0] >= var[1] && var[0] >= var[2] {
            0
        } else if var[1] >= var[2] {
            1
        } else {
            2
        };
        self.endpoints_x.clear();
        for (id, aabb) in &self.aabbs {
            self.endpoints_x.push(SapEndpointU32 {
                value: aabb.min[best],
                body_id: *id,
                is_min: true,
            });
            self.endpoints_x.push(SapEndpointU32 {
                value: aabb.max[best],
                body_id: *id,
                is_min: false,
            });
        }
        self.endpoints_x.sort_by(|a, b| {
            a.value
                .partial_cmp(&b.value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        best
    }
}
impl IncrementalSap {
    /// Return `true` if any body in the SAP overlaps `aabb`.
    pub fn any_overlap(&self, aabb: &Aabb3) -> bool {
        for stored in self.aabbs.values() {
            if aabb3_overlaps(stored, aabb) {
                return true;
            }
        }
        false
    }
    /// Return the ids of all bodies whose AABB overlaps `aabb`.
    pub fn query_aabb(&self, aabb: &Aabb3) -> Vec<u32> {
        self.aabbs
            .iter()
            .filter_map(|(&id, stored)| {
                if aabb3_overlaps(stored, aabb) {
                    Some(id)
                } else {
                    None
                }
            })
            .collect()
    }
    /// Return the ids of all bodies within squared distance `radius_sq` of point `p`.
    pub fn query_sphere_sq(&self, p: [f64; 3], radius_sq: f64) -> Vec<u32> {
        self.aabbs
            .iter()
            .filter_map(|(&id, aabb)| {
                if aabb3_point_dist_sq(aabb, p) <= radius_sq {
                    Some(id)
                } else {
                    None
                }
            })
            .collect()
    }
    /// Remove all bodies and reset to empty.
    pub fn clear(&mut self) {
        self.endpoints_x.clear();
        self.endpoints_y.clear();
        self.endpoints_z.clear();
        self.aabbs.clear();
        self.active_pairs.clear();
    }
    /// Return the AABB for `id`, or `None` if not present.
    pub fn get_aabb(&self, id: u32) -> Option<&Aabb3> {
        self.aabbs.get(&id)
    }
    /// Return an iterator over all tracked body ids.
    pub fn body_ids(&self) -> impl Iterator<Item = u32> + '_ {
        self.aabbs.keys().copied()
    }
    /// Batch-insert all bodies from `other` into `self`.
    ///
    /// Bodies already in `self` are overwritten (their AABBs updated).
    pub fn merge_from(&mut self, other: &IncrementalSap) {
        for (&id, aabb) in &other.aabbs {
            if self.aabbs.contains_key(&id) {
                self.update(id, aabb.clone());
            } else {
                self.insert(id, aabb.clone());
            }
        }
    }
    /// Return `true` if `id` is currently tracked.
    pub fn contains(&self, id: u32) -> bool {
        self.aabbs.contains_key(&id)
    }
}
/// A multi-phase SAP that integrates event tracking and statistics.
pub struct MultiPhaseSap {
    pub(super) event_sap: EventDrivenSap,
    pub(super) stat_sap: StatTrackingSap,
}
impl MultiPhaseSap {
    /// Create a new multi-phase SAP.
    pub fn new() -> Self {
        Self {
            event_sap: EventDrivenSap::new(),
            stat_sap: StatTrackingSap::new(),
        }
    }
    /// Insert a body into both SAP structures.
    pub fn insert(&mut self, id: u32, aabb: Aabb3) {
        self.event_sap.insert(id, aabb.clone());
        self.stat_sap.insert(id, aabb);
    }
    /// Remove a body from both SAP structures.
    pub fn remove(&mut self, id: u32) {
        self.event_sap.remove(id);
        self.stat_sap.remove(id);
    }
    /// Update a body's AABB in both SAP structures.
    pub fn update(&mut self, id: u32, aabb: Aabb3) {
        self.event_sap.update(id, aabb.clone());
        self.stat_sap.update(id, aabb);
    }
    /// Perform a full broadphase step and return the combined result.
    pub fn step(&mut self) -> BroadphaseResult {
        let pairs = self.stat_sap.query_pairs();
        let stats = self.stat_sap.stats.clone();
        let (events, _) = self.event_sap.step_pairs();
        let mut new_pairs = Vec::new();
        let mut lost_pairs = Vec::new();
        for ev in events {
            match ev {
                OverlapEvent::Begin(a, b) => new_pairs.push((a, b)),
                OverlapEvent::End(a, b) => lost_pairs.push((a, b)),
            }
        }
        BroadphaseResult {
            pairs,
            new_pairs,
            lost_pairs,
            stats,
        }
    }
}
