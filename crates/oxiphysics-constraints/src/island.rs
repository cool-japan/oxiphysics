// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Island manager: groups bodies into connected components for per-island sleeping.
//!
//! Uses union-find (disjoint set union with path compression + union by rank)
//! to merge bodies that share a contact pair. Bodies not in any contact pair
//! form singleton islands.
//!
//! Extended features:
//! - Island splitting when contact pairs are removed between frames
//! - Island merging when new contact pairs are added
//! - Island sleep detection with configurable KE thresholds
//! - Island iteration utilities (body/constraint counts, max velocity)
//! - Island statistics for profiling and debugging

use oxiphysics_core::BodyHandle;
use std::collections::HashMap;

/// A connected component of touching/constrained bodies.
#[derive(Debug, Clone)]
pub struct Island {
    /// Body handles belonging to this island.
    pub body_handles: Vec<BodyHandle>,
    /// Whether the entire island is currently sleeping.
    pub is_sleeping: bool,
    /// Accumulated time that all bodies in the island have been below the KE threshold.
    pub sleep_timer: f64,
}

impl Island {
    fn new(body_handles: Vec<BodyHandle>) -> Self {
        Self {
            body_handles,
            is_sleeping: false,
            sleep_timer: 0.0,
        }
    }
}

/// Statistics for a single island.
#[derive(Debug, Clone, Default)]
pub struct IslandStats {
    /// Number of bodies in the island.
    pub body_count: usize,
    /// Number of contact pairs internal to the island.
    pub constraint_count: usize,
    /// Maximum linear speed among bodies in the island.
    pub max_linear_speed: f64,
    /// Maximum angular speed among bodies in the island.
    pub max_angular_speed: f64,
    /// Total kinetic energy of the island.
    pub total_kinetic_energy: f64,
    /// Whether the island is currently sleeping.
    pub is_sleeping: bool,
    /// How long the island has been below the sleep threshold.
    pub sleep_timer: f64,
}

/// Aggregate statistics across all islands.
#[derive(Debug, Clone, Default)]
pub struct IslandManagerStats {
    /// Total number of islands.
    pub island_count: usize,
    /// Number of sleeping islands.
    pub sleeping_count: usize,
    /// Number of awake islands.
    pub awake_count: usize,
    /// Total bodies across all islands.
    pub total_bodies: usize,
    /// Largest island body count.
    pub largest_island_size: usize,
    /// Smallest island body count.
    pub smallest_island_size: usize,
    /// Average island size (bodies per island).
    pub average_island_size: f64,
}

/// Union-Find (Disjoint Set Union) helper.
struct Dsu {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl Dsu {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            // Path halving
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            self.parent[ra] = rb;
        } else if self.rank[ra] > self.rank[rb] {
            self.parent[rb] = ra;
        } else {
            self.parent[rb] = ra;
            self.rank[ra] += 1;
        }
    }

    /// Check whether two elements are in the same set.
    #[cfg(test)]
    fn connected(&mut self, a: usize, b: usize) -> bool {
        self.find(a) == self.find(b)
    }

    /// Count the number of distinct sets.
    #[cfg(test)]
    fn set_count(&mut self) -> usize {
        let n = self.parent.len();
        let mut roots = std::collections::HashSet::new();
        for i in 0..n {
            roots.insert(self.find(i));
        }
        roots.len()
    }
}

/// Manages grouping of bodies into islands (connected components) for per-island sleeping.
#[derive(Debug, Default)]
pub struct IslandManager {
    islands: Vec<Island>,
    /// Map from `BodyHandle.index` to island index.
    body_to_island: HashMap<u32, usize>,
}

impl IslandManager {
    /// Create a new, empty island manager.
    pub fn new() -> Self {
        Self {
            islands: Vec::new(),
            body_to_island: HashMap::new(),
        }
    }

    /// Build islands from a set of active bodies and contact pairs.
    ///
    /// All bodies in `active_bodies` are assigned to an island. Bodies that share
    /// an entry in `contact_pairs` are merged into the same island. Bodies not
    /// appearing in any contact pair form singleton islands.
    ///
    /// Previous island data (including sleep timers) is preserved when the same
    /// body ends up in an island with the same sleep state.
    pub fn build_islands(
        &mut self,
        active_bodies: &[BodyHandle],
        contact_pairs: &[(BodyHandle, BodyHandle)],
    ) {
        if active_bodies.is_empty() {
            self.islands.clear();
            self.body_to_island.clear();
            return;
        }

        // Assign a local index to every body.
        let mut handle_to_local: HashMap<u32, usize> = HashMap::with_capacity(active_bodies.len());
        for (i, h) in active_bodies.iter().enumerate() {
            handle_to_local.insert(h.index, i);
        }

        let n = active_bodies.len();
        let mut dsu = Dsu::new(n);

        // Union bodies connected by contact pairs.
        for (ha, hb) in contact_pairs {
            if let (Some(&ia), Some(&ib)) = (
                handle_to_local.get(&ha.index),
                handle_to_local.get(&hb.index),
            ) {
                dsu.union(ia, ib);
            }
        }

        // Collect previous sleep timers keyed by the *set* of body indices (represented
        // by the canonical root's first body for simplicity — we carry the timer forward
        // on a best-effort basis using the old body_to_island mapping).
        let old_sleep_timer: HashMap<u32, f64> = self
            .islands
            .iter()
            .flat_map(|isl| {
                let timer = isl.sleep_timer;
                isl.body_handles.iter().map(move |h| (h.index, timer))
            })
            .collect();

        // Group bodies by DSU root.
        let mut root_to_bodies: HashMap<usize, Vec<BodyHandle>> = HashMap::new();
        for (i, h) in active_bodies.iter().enumerate() {
            let root = dsu.find(i);
            root_to_bodies.entry(root).or_default().push(*h);
        }

        // Rebuild islands and body_to_island map.
        let mut new_islands: Vec<Island> = Vec::with_capacity(root_to_bodies.len());
        let mut new_body_to_island: HashMap<u32, usize> = HashMap::with_capacity(n);

        for (_root, bodies) in root_to_bodies {
            let island_idx = new_islands.len();

            // Carry forward the minimum sleep timer from the old mapping (conservative).
            let min_old_timer = bodies
                .iter()
                .filter_map(|h| old_sleep_timer.get(&h.index).copied())
                .fold(f64::INFINITY, f64::min);
            let sleep_timer = if min_old_timer == f64::INFINITY {
                0.0
            } else {
                min_old_timer
            };

            let was_sleeping = bodies.iter().all(|h| {
                self.body_to_island
                    .get(&h.index)
                    .and_then(|&idx| self.islands.get(idx))
                    .map(|isl| isl.is_sleeping)
                    .unwrap_or(false)
            });

            let mut island = Island::new(bodies.clone());
            island.sleep_timer = sleep_timer;
            island.is_sleeping = was_sleeping;

            for h in &bodies {
                new_body_to_island.insert(h.index, island_idx);
            }

            new_islands.push(island);
        }

        self.islands = new_islands;
        self.body_to_island = new_body_to_island;
    }

    /// Update sleep state of each island.
    ///
    /// An island goes to sleep when ALL bodies in it have kinetic energy below
    /// `ke_threshold` for a continuous duration exceeding `sleep_time` seconds.
    /// If any body in an island exceeds the threshold the timer resets to zero
    /// and the island wakes up.
    pub fn update_sleep<F>(
        &mut self,
        dt: f64,
        ke_threshold: f64,
        sleep_time: f64,
        kinetic_energy: F,
    ) where
        F: Fn(BodyHandle) -> f64,
    {
        for island in &mut self.islands {
            let all_quiet = island
                .body_handles
                .iter()
                .all(|&h| kinetic_energy(h) < ke_threshold);

            if all_quiet {
                island.sleep_timer += dt;
                if island.sleep_timer >= sleep_time {
                    island.is_sleeping = true;
                }
            } else {
                island.sleep_timer = 0.0;
                island.is_sleeping = false;
            }
        }
    }

    /// Return the island index for a body, if any.
    pub fn island_of(&self, handle: BodyHandle) -> Option<usize> {
        self.body_to_island.get(&handle.index).copied()
    }

    /// Iterate over all islands.
    pub fn islands(&self) -> &[Island] {
        &self.islands
    }

    /// Wake up the island that contains the given body handle.
    ///
    /// Resets the sleep timer and clears the sleeping flag for the entire island.
    /// Does nothing if the handle is not tracked.
    pub fn wake_island(&mut self, handle: BodyHandle) {
        if let Some(&idx) = self.body_to_island.get(&handle.index)
            && let Some(island) = self.islands.get_mut(idx)
        {
            island.is_sleeping = false;
            island.sleep_timer = 0.0;
        }
    }

    // ── Island splitting ──────────────────────────────────────────────────

    /// Detect islands that should be split because contact pairs have been removed.
    ///
    /// Compares the current island state against a new set of contact pairs
    /// and returns the indices of islands that need to be rebuilt (i.e., they
    /// have bodies that are no longer connected by the new contact set).
    pub fn detect_splits(&self, new_contact_pairs: &[(BodyHandle, BodyHandle)]) -> Vec<usize> {
        // Build adjacency from new contacts
        let mut adjacency: HashMap<u32, Vec<u32>> = HashMap::new();
        for (ha, hb) in new_contact_pairs {
            adjacency.entry(ha.index).or_default().push(hb.index);
            adjacency.entry(hb.index).or_default().push(ha.index);
        }

        let mut split_indices = Vec::new();

        for (island_idx, island) in self.islands.iter().enumerate() {
            if island.body_handles.len() <= 1 {
                continue;
            }

            // BFS/flood from the first body; if not all bodies are reachable,
            // the island needs splitting.
            let first = island.body_handles[0].index;
            let body_set: std::collections::HashSet<u32> =
                island.body_handles.iter().map(|h| h.index).collect();

            let mut visited = std::collections::HashSet::new();
            let mut queue = std::collections::VecDeque::new();
            visited.insert(first);
            queue.push_back(first);

            while let Some(current) = queue.pop_front() {
                if let Some(neighbors) = adjacency.get(&current) {
                    for &nb in neighbors {
                        if body_set.contains(&nb) && visited.insert(nb) {
                            queue.push_back(nb);
                        }
                    }
                }
            }

            if visited.len() < body_set.len() {
                split_indices.push(island_idx);
            }
        }

        split_indices
    }

    // ── Island merging ────────────────────────────────────────────────────

    /// Detect pairs of islands that should be merged because a new contact pair
    /// connects bodies from different islands.
    ///
    /// Returns pairs `(island_a, island_b)` where `island_a < island_b`.
    pub fn detect_merges(
        &self,
        new_contact_pairs: &[(BodyHandle, BodyHandle)],
    ) -> Vec<(usize, usize)> {
        let mut merge_pairs = Vec::new();

        for (ha, hb) in new_contact_pairs {
            let ia = self.body_to_island.get(&ha.index).copied();
            let ib = self.body_to_island.get(&hb.index).copied();
            if let (Some(a), Some(b)) = (ia, ib)
                && a != b
            {
                let pair = if a < b { (a, b) } else { (b, a) };
                if !merge_pairs.contains(&pair) {
                    merge_pairs.push(pair);
                }
            }
        }

        merge_pairs
    }

    // ── Island sleep detection ────────────────────────────────────────────

    /// Check whether a specific island is a candidate for sleeping based on
    /// per-body velocity thresholds.
    ///
    /// Returns `true` if all bodies have both linear and angular speeds below
    /// the given thresholds.
    pub fn is_island_sleep_candidate<F, G>(
        &self,
        island_idx: usize,
        linear_speed: F,
        angular_speed: G,
        linear_threshold: f64,
        angular_threshold: f64,
    ) -> bool
    where
        F: Fn(BodyHandle) -> f64,
        G: Fn(BodyHandle) -> f64,
    {
        if let Some(island) = self.islands.get(island_idx) {
            island.body_handles.iter().all(|&h| {
                linear_speed(h) < linear_threshold && angular_speed(h) < angular_threshold
            })
        } else {
            false
        }
    }

    /// Force all islands to wake up. Useful after a global event (e.g. explosion).
    pub fn wake_all(&mut self) {
        for island in &mut self.islands {
            island.is_sleeping = false;
            island.sleep_timer = 0.0;
        }
    }

    /// Force a specific island to sleep immediately.
    pub fn force_sleep(&mut self, island_idx: usize) {
        if let Some(island) = self.islands.get_mut(island_idx) {
            island.is_sleeping = true;
        }
    }

    // ── Island iteration utilities ────────────────────────────────────────

    /// Return handles of all bodies in the island at the given index.
    pub fn bodies_in_island(&self, island_idx: usize) -> &[BodyHandle] {
        self.islands
            .get(island_idx)
            .map(|i| i.body_handles.as_slice())
            .unwrap_or(&[])
    }

    /// Return indices of all awake islands.
    pub fn awake_island_indices(&self) -> Vec<usize> {
        self.islands
            .iter()
            .enumerate()
            .filter(|(_, i)| !i.is_sleeping)
            .map(|(idx, _)| idx)
            .collect()
    }

    /// Return indices of all sleeping islands.
    pub fn sleeping_island_indices(&self) -> Vec<usize> {
        self.islands
            .iter()
            .enumerate()
            .filter(|(_, i)| i.is_sleeping)
            .map(|(idx, _)| idx)
            .collect()
    }

    /// Count the number of contact pairs internal to a given island.
    pub fn count_constraints_in_island(
        &self,
        island_idx: usize,
        contact_pairs: &[(BodyHandle, BodyHandle)],
    ) -> usize {
        let Some(island) = self.islands.get(island_idx) else {
            return 0;
        };
        let body_set: std::collections::HashSet<u32> =
            island.body_handles.iter().map(|h| h.index).collect();
        contact_pairs
            .iter()
            .filter(|(ha, hb)| body_set.contains(&ha.index) && body_set.contains(&hb.index))
            .count()
    }

    /// Find the island with the most bodies.
    pub fn largest_island_index(&self) -> Option<usize> {
        self.islands
            .iter()
            .enumerate()
            .max_by_key(|(_, i)| i.body_handles.len())
            .map(|(idx, _)| idx)
    }

    /// Total number of bodies across all islands.
    pub fn total_body_count(&self) -> usize {
        self.islands.iter().map(|i| i.body_handles.len()).sum()
    }

    /// Number of islands.
    pub fn island_count(&self) -> usize {
        self.islands.len()
    }

    // ── Island statistics ─────────────────────────────────────────────────

    /// Compute statistics for a single island.
    pub fn island_stats<F, G, H>(
        &self,
        island_idx: usize,
        contact_pairs: &[(BodyHandle, BodyHandle)],
        linear_speed: F,
        angular_speed: G,
        kinetic_energy: H,
    ) -> Option<IslandStats>
    where
        F: Fn(BodyHandle) -> f64,
        G: Fn(BodyHandle) -> f64,
        H: Fn(BodyHandle) -> f64,
    {
        let island = self.islands.get(island_idx)?;

        let body_count = island.body_handles.len();

        let body_set: std::collections::HashSet<u32> =
            island.body_handles.iter().map(|h| h.index).collect();
        let constraint_count = contact_pairs
            .iter()
            .filter(|(ha, hb)| body_set.contains(&ha.index) && body_set.contains(&hb.index))
            .count();

        let mut max_linear = 0.0_f64;
        let mut max_angular = 0.0_f64;
        let mut total_ke = 0.0_f64;

        for &h in &island.body_handles {
            max_linear = max_linear.max(linear_speed(h));
            max_angular = max_angular.max(angular_speed(h));
            total_ke += kinetic_energy(h);
        }

        Some(IslandStats {
            body_count,
            constraint_count,
            max_linear_speed: max_linear,
            max_angular_speed: max_angular,
            total_kinetic_energy: total_ke,
            is_sleeping: island.is_sleeping,
            sleep_timer: island.sleep_timer,
        })
    }

    /// Compute aggregate statistics across all islands.
    pub fn manager_stats(&self) -> IslandManagerStats {
        let island_count = self.islands.len();
        if island_count == 0 {
            return IslandManagerStats::default();
        }

        let sleeping_count = self.islands.iter().filter(|i| i.is_sleeping).count();
        let awake_count = island_count - sleeping_count;
        let total_bodies: usize = self.islands.iter().map(|i| i.body_handles.len()).sum();
        let largest = self
            .islands
            .iter()
            .map(|i| i.body_handles.len())
            .max()
            .unwrap_or(0);
        let smallest = self
            .islands
            .iter()
            .map(|i| i.body_handles.len())
            .min()
            .unwrap_or(0);
        let average = total_bodies as f64 / island_count as f64;

        IslandManagerStats {
            island_count,
            sleeping_count,
            awake_count,
            total_bodies,
            largest_island_size: largest,
            smallest_island_size: smallest,
            average_island_size: average,
        }
    }

    /// Collect the handles of all bodies that are currently in sleeping islands.
    pub fn sleeping_bodies(&self) -> Vec<BodyHandle> {
        self.islands
            .iter()
            .filter(|i| i.is_sleeping)
            .flat_map(|i| i.body_handles.iter().copied())
            .collect()
    }

    /// Collect the handles of all bodies that are currently in awake islands.
    pub fn awake_bodies(&self) -> Vec<BodyHandle> {
        self.islands
            .iter()
            .filter(|i| !i.is_sleeping)
            .flat_map(|i| i.body_handles.iter().copied())
            .collect()
    }

    /// Wake up all islands that contain any body from the given list.
    ///
    /// Useful when an external force or collision affects specific bodies.
    pub fn wake_bodies(&mut self, handles: &[BodyHandle]) {
        for h in handles {
            self.wake_island(*h);
        }
    }

    /// Check whether a body is in a sleeping island.
    pub fn is_body_sleeping(&self, handle: BodyHandle) -> bool {
        self.body_to_island
            .get(&handle.index)
            .and_then(|&idx| self.islands.get(idx))
            .map(|i| i.is_sleeping)
            .unwrap_or(false)
    }

    // ── Energy-based sleep with hysteresis ────────────────────────────────

    /// Update sleep state using energy-based hysteresis.
    ///
    /// An island enters sleep when its total KE stays below `sleep_threshold`
    /// for `sleep_frames` consecutive frames.  It wakes when KE exceeds
    /// `wake_threshold` (which should be ≥ `sleep_threshold`).
    ///
    /// Returns the number of islands whose sleep state changed.
    pub fn update_sleep_hysteresis<F>(
        &mut self,
        island_ke: F,
        sleep_threshold: f64,
        wake_threshold: f64,
        sleep_frames: u32,
        frame_counters: &mut Vec<u32>,
    ) -> usize
    where
        F: Fn(usize) -> f64, // argument is island index
    {
        // Ensure frame counter vector is sized correctly
        while frame_counters.len() < self.islands.len() {
            frame_counters.push(0);
        }
        frame_counters.truncate(self.islands.len());

        let mut changed = 0usize;

        for (i, island) in self.islands.iter_mut().enumerate() {
            let ke = island_ke(i);
            if island.is_sleeping {
                // Wake if energy spikes above wake_threshold
                if ke > wake_threshold {
                    island.is_sleeping = false;
                    island.sleep_timer = 0.0;
                    frame_counters[i] = 0;
                    changed += 1;
                }
            } else {
                if ke < sleep_threshold {
                    frame_counters[i] += 1;
                    if frame_counters[i] >= sleep_frames {
                        island.is_sleeping = true;
                        changed += 1;
                    }
                } else {
                    frame_counters[i] = 0;
                }
            }
        }
        changed
    }

    // ── Priority-sorted island iteration ─────────────────────────────────

    /// Return island indices sorted by priority (highest energy first).
    ///
    /// Useful for adaptive sub-stepping where the most active islands
    /// should be processed first.
    pub fn sorted_island_indices_by_energy<F>(&self, island_ke: F) -> Vec<usize>
    where
        F: Fn(usize) -> f64,
    {
        let mut indices: Vec<usize> = (0..self.islands.len())
            .filter(|&i| !self.islands[i].is_sleeping)
            .collect();
        indices.sort_by(|&a, &b| {
            let ke_a = island_ke(a);
            let ke_b = island_ke(b);
            ke_b.partial_cmp(&ke_a).unwrap_or(std::cmp::Ordering::Equal)
        });
        indices
    }

    // ── Island body-count histogram ───────────────────────────────────────

    /// Compute a histogram of island sizes.
    ///
    /// Returns a `Vec` where index `i` holds the number of islands with
    /// exactly `i + 1` bodies.  Islands larger than `max_size` are bucketed
    /// into the last slot.
    pub fn island_size_histogram(&self, max_size: usize) -> Vec<usize> {
        let max_size = max_size.max(1);
        let mut hist = vec![0usize; max_size];
        for island in &self.islands {
            let sz = island.body_handles.len().saturating_sub(1);
            let idx = sz.min(max_size - 1);
            hist[idx] += 1;
        }
        hist
    }

    // ── Cross-island contact detection ────────────────────────────────────

    /// Find all contact pairs that span two different islands.
    ///
    /// Such cross-island contacts indicate that a new constraint has been
    /// created that will merge the two islands on the next `build_islands`.
    pub fn cross_island_contacts<'a>(
        &self,
        contact_pairs: &'a [(BodyHandle, BodyHandle)],
    ) -> Vec<&'a (BodyHandle, BodyHandle)> {
        contact_pairs
            .iter()
            .filter(|(ha, hb)| {
                let ia = self.body_to_island.get(&ha.index).copied();
                let ib = self.body_to_island.get(&hb.index).copied();
                match (ia, ib) {
                    (Some(a), Some(b)) => a != b,
                    _ => false,
                }
            })
            .collect()
    }

    // ── Island energy threshold classification ────────────────────────────

    /// Classify islands into "hot", "warm", and "cold" categories based on KE.
    ///
    /// * Hot: ke ≥ `hot_threshold`
    /// * Warm: `cold_threshold` ≤ ke < `hot_threshold`
    /// * Cold: ke < `cold_threshold`
    ///
    /// Returns `(hot, warm, cold)` index lists.
    pub fn classify_islands_by_energy<F>(
        &self,
        island_ke: F,
        hot_threshold: f64,
        cold_threshold: f64,
    ) -> (Vec<usize>, Vec<usize>, Vec<usize>)
    where
        F: Fn(usize) -> f64,
    {
        let mut hot = Vec::new();
        let mut warm = Vec::new();
        let mut cold = Vec::new();
        for i in 0..self.islands.len() {
            let ke = island_ke(i);
            if ke >= hot_threshold {
                hot.push(i);
            } else if ke >= cold_threshold {
                warm.push(i);
            } else {
                cold.push(i);
            }
        }
        (hot, warm, cold)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn handle(index: u32) -> BodyHandle {
        BodyHandle::new(index, 0)
    }

    // -----------------------------------------------------------------------
    // test_island_single_body
    // -----------------------------------------------------------------------
    #[test]
    fn test_island_single_body() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        mgr.build_islands(&[h0], &[]);

        assert_eq!(mgr.islands().len(), 1, "one body -> one singleton island");
        assert_eq!(mgr.islands()[0].body_handles.len(), 1);
        assert_eq!(mgr.island_of(h0), Some(0));
    }

    // -----------------------------------------------------------------------
    // test_island_two_touching_bodies
    // -----------------------------------------------------------------------
    #[test]
    fn test_island_two_touching_bodies() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[(h0, h1)]);

        assert_eq!(
            mgr.islands().len(),
            1,
            "two bodies in contact -> one island"
        );
        assert_eq!(mgr.islands()[0].body_handles.len(), 2);

        let i0 = mgr.island_of(h0).expect("h0 should be in an island");
        let i1 = mgr.island_of(h1).expect("h1 should be in an island");
        assert_eq!(i0, i1, "both bodies should be in the same island");
    }

    // -----------------------------------------------------------------------
    // test_island_separate_groups
    // -----------------------------------------------------------------------
    #[test]
    fn test_island_separate_groups() {
        let mut mgr = IslandManager::new();
        let ha = handle(0);
        let hb = handle(1);
        let hc = handle(2);
        let hd = handle(3);
        // A-B in contact; C-D in contact; no cross contacts.
        mgr.build_islands(&[ha, hb, hc, hd], &[(ha, hb), (hc, hd)]);

        assert_eq!(mgr.islands().len(), 2, "should produce two islands");
        let mut sizes: Vec<usize> = mgr.islands().iter().map(|i| i.body_handles.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![2, 2]);

        // A and B must share an island; C and D must share an island.
        assert_eq!(mgr.island_of(ha), mgr.island_of(hb));
        assert_eq!(mgr.island_of(hc), mgr.island_of(hd));
        assert_ne!(
            mgr.island_of(ha),
            mgr.island_of(hc),
            "the two groups must be in different islands"
        );
    }

    // -----------------------------------------------------------------------
    // test_island_sleep_timer
    // -----------------------------------------------------------------------
    #[test]
    fn test_island_sleep_timer() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[(h0, h1)]);

        // All bodies have zero KE — timer should advance and island should sleep.
        let ke = |_: BodyHandle| 0.0_f64;
        let dt = 0.1_f64;
        let ke_threshold = 1.0_f64;
        let sleep_time = 0.5_f64;

        for step in 1..=6 {
            mgr.update_sleep(dt, ke_threshold, sleep_time, ke);
            let timer = mgr.islands()[0].sleep_timer;
            // Timer should grow monotonically.
            assert!(
                timer > 0.0,
                "sleep timer should be positive after step {step}"
            );
        }

        // After 6 x 0.1 s = 0.6 s > sleep_time (0.5 s), island must be asleep.
        assert!(
            mgr.islands()[0].is_sleeping,
            "island should be sleeping after timer exceeds sleep_time"
        );
    }

    // -----------------------------------------------------------------------
    // test_island_wake_on_contact
    // -----------------------------------------------------------------------
    #[test]
    fn test_island_wake_on_contact() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[(h0, h1)]);

        // Force the island to sleep.
        let ke = |_: BodyHandle| 0.0_f64;
        for _ in 0..10 {
            mgr.update_sleep(0.1, 1.0, 0.5, ke);
        }
        assert!(mgr.islands()[0].is_sleeping, "island should be sleeping");

        // Wake the island via wake_island.
        mgr.wake_island(h0);

        assert!(
            !mgr.islands()[0].is_sleeping,
            "island should be awake after wake_island"
        );
        assert_eq!(
            mgr.islands()[0].sleep_timer,
            0.0,
            "sleep timer should be reset to zero"
        );
    }

    // -----------------------------------------------------------------------
    // Island splitting tests
    // -----------------------------------------------------------------------
    #[test]
    fn test_detect_splits_no_split() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[(h0, h1)]);

        // Same contacts: no split needed
        let splits = mgr.detect_splits(&[(h0, h1)]);
        assert!(
            splits.is_empty(),
            "no split expected when contacts unchanged"
        );
    }

    #[test]
    fn test_detect_splits_contact_removed() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        let h2 = handle(2);
        // Chain: 0-1-2
        mgr.build_islands(&[h0, h1, h2], &[(h0, h1), (h1, h2)]);
        assert_eq!(mgr.islands().len(), 1, "chain should be one island");

        // Remove the 1-2 contact: 0-1 and 2 should split
        let splits = mgr.detect_splits(&[(h0, h1)]);
        assert_eq!(splits.len(), 1, "one island should need splitting");
    }

    #[test]
    fn test_detect_splits_singleton_no_split() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        mgr.build_islands(&[h0], &[]);

        let splits = mgr.detect_splits(&[]);
        assert!(splits.is_empty(), "singleton never needs splitting");
    }

    // -----------------------------------------------------------------------
    // Island merging tests
    // -----------------------------------------------------------------------
    #[test]
    fn test_detect_merges() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        let h2 = handle(2);
        // Two separate islands: {0,1} and {2}
        mgr.build_islands(&[h0, h1, h2], &[(h0, h1)]);
        assert_eq!(mgr.islands().len(), 2);

        // New contact connects 1 and 2
        let merges = mgr.detect_merges(&[(h1, h2)]);
        assert_eq!(merges.len(), 1, "one merge expected");
        let (a, b) = merges[0];
        assert!(a < b, "pair should be ordered");
    }

    #[test]
    fn test_detect_merges_same_island() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[(h0, h1)]);

        // Contact within the same island: no merge needed
        let merges = mgr.detect_merges(&[(h0, h1)]);
        assert!(merges.is_empty(), "no merge for same-island contact");
    }

    // -----------------------------------------------------------------------
    // Sleep detection tests
    // -----------------------------------------------------------------------
    #[test]
    fn test_is_island_sleep_candidate() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[(h0, h1)]);

        // All bodies below threshold
        let lin = |_: BodyHandle| 0.001_f64;
        let ang = |_: BodyHandle| 0.001_f64;
        assert!(mgr.is_island_sleep_candidate(0, lin, ang, 0.01, 0.01));

        // One body above threshold
        let lin2 = |h: BodyHandle| if h.index == 0 { 0.1 } else { 0.001 };
        assert!(!mgr.is_island_sleep_candidate(0, lin2, ang, 0.01, 0.01));
    }

    #[test]
    fn test_wake_all() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[(h0, h1)]);

        // Put to sleep
        let ke = |_: BodyHandle| 0.0_f64;
        for _ in 0..10 {
            mgr.update_sleep(0.1, 1.0, 0.5, ke);
        }
        assert!(mgr.islands()[0].is_sleeping);

        mgr.wake_all();
        assert!(!mgr.islands()[0].is_sleeping);
        assert_eq!(mgr.islands()[0].sleep_timer, 0.0);
    }

    #[test]
    fn test_force_sleep() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        mgr.build_islands(&[h0], &[]);
        assert!(!mgr.islands()[0].is_sleeping);

        mgr.force_sleep(0);
        assert!(mgr.islands()[0].is_sleeping);
    }

    // -----------------------------------------------------------------------
    // Island iteration utility tests
    // -----------------------------------------------------------------------
    #[test]
    fn test_bodies_in_island() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        let h2 = handle(2);
        mgr.build_islands(&[h0, h1, h2], &[(h0, h1)]);

        // One island has 2 bodies, another has 1
        let largest = mgr.largest_island_index().unwrap();
        assert_eq!(mgr.bodies_in_island(largest).len(), 2);

        // Out-of-bounds returns empty
        assert!(mgr.bodies_in_island(999).is_empty());
    }

    #[test]
    fn test_awake_and_sleeping_indices() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[]);
        // Two singleton islands, both awake
        assert_eq!(mgr.awake_island_indices().len(), 2);
        assert_eq!(mgr.sleeping_island_indices().len(), 0);

        // Put first island to sleep
        mgr.force_sleep(0);
        assert_eq!(mgr.awake_island_indices().len(), 1);
        assert_eq!(mgr.sleeping_island_indices().len(), 1);
    }

    #[test]
    fn test_count_constraints_in_island() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        let h2 = handle(2);
        let contacts = [(h0, h1), (h1, h0)]; // duplicate pair
        mgr.build_islands(&[h0, h1, h2], &contacts);

        let idx = mgr.island_of(h0).unwrap();
        let count = mgr.count_constraints_in_island(idx, &contacts);
        assert_eq!(count, 2, "both contact entries should be counted");
    }

    #[test]
    fn test_total_body_count() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        let h2 = handle(2);
        mgr.build_islands(&[h0, h1, h2], &[(h0, h1)]);
        assert_eq!(mgr.total_body_count(), 3);
        assert_eq!(mgr.island_count(), 2);
    }

    // -----------------------------------------------------------------------
    // Island statistics tests
    // -----------------------------------------------------------------------
    #[test]
    fn test_island_stats() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        let contacts = [(h0, h1)];
        mgr.build_islands(&[h0, h1], &contacts);

        let lin = |h: BodyHandle| if h.index == 0 { 3.0 } else { 1.0 };
        let ang = |h: BodyHandle| if h.index == 0 { 0.5 } else { 2.0 };
        let ke = |h: BodyHandle| if h.index == 0 { 10.0 } else { 5.0 };

        let stats = mgr.island_stats(0, &contacts, lin, ang, ke).unwrap();
        assert_eq!(stats.body_count, 2);
        assert_eq!(stats.constraint_count, 1);
        assert!((stats.max_linear_speed - 3.0).abs() < 1e-12);
        assert!((stats.max_angular_speed - 2.0).abs() < 1e-12);
        assert!((stats.total_kinetic_energy - 15.0).abs() < 1e-12);
    }

    #[test]
    fn test_manager_stats() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        let h2 = handle(2);
        mgr.build_islands(&[h0, h1, h2], &[(h0, h1)]);

        mgr.force_sleep(mgr.island_of(h2).unwrap());

        let stats = mgr.manager_stats();
        assert_eq!(stats.island_count, 2);
        assert_eq!(stats.sleeping_count, 1);
        assert_eq!(stats.awake_count, 1);
        assert_eq!(stats.total_bodies, 3);
        assert_eq!(stats.largest_island_size, 2);
        assert_eq!(stats.smallest_island_size, 1);
        assert!((stats.average_island_size - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_manager_stats_empty() {
        let mgr = IslandManager::new();
        let stats = mgr.manager_stats();
        assert_eq!(stats.island_count, 0);
        assert_eq!(stats.total_bodies, 0);
    }

    // -----------------------------------------------------------------------
    // Sleeping / awake body collection tests
    // -----------------------------------------------------------------------
    #[test]
    fn test_sleeping_and_awake_bodies() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        let h2 = handle(2);
        mgr.build_islands(&[h0, h1, h2], &[(h0, h1)]);

        // Put the singleton island (h2) to sleep
        mgr.force_sleep(mgr.island_of(h2).unwrap());

        let sleeping = mgr.sleeping_bodies();
        assert_eq!(sleeping.len(), 1);
        assert_eq!(sleeping[0].index, 2);

        let awake = mgr.awake_bodies();
        assert_eq!(awake.len(), 2);
    }

    #[test]
    fn test_is_body_sleeping() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[]);

        assert!(!mgr.is_body_sleeping(h0));
        // Sleep the island that contains h0
        let island_of_h0 = mgr.island_of(h0).unwrap();
        mgr.force_sleep(island_of_h0);
        assert!(mgr.is_body_sleeping(h0));
        // h1 is sleeping only if it shares the same island as h0
        let island_of_h1 = mgr.island_of(h1).unwrap();
        if island_of_h0 != island_of_h1 {
            assert!(!mgr.is_body_sleeping(h1));
        }
    }

    #[test]
    fn test_wake_bodies() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        let h2 = handle(2);
        mgr.build_islands(&[h0, h1, h2], &[]);

        // Sleep all
        for i in 0..mgr.island_count() {
            mgr.force_sleep(i);
        }
        assert!(mgr.islands().iter().all(|i| i.is_sleeping));

        // Wake specific bodies
        mgr.wake_bodies(&[h0, h2]);
        assert!(!mgr.is_body_sleeping(h0));
        assert!(!mgr.is_body_sleeping(h2));
    }

    // -----------------------------------------------------------------------
    // DSU internal tests
    // -----------------------------------------------------------------------
    #[test]
    fn test_dsu_connected() {
        let mut dsu = Dsu::new(5);
        dsu.union(0, 1);
        dsu.union(2, 3);
        assert!(dsu.connected(0, 1));
        assert!(dsu.connected(2, 3));
        assert!(!dsu.connected(0, 2));

        dsu.union(1, 3);
        assert!(dsu.connected(0, 3));
    }

    #[test]
    fn test_dsu_set_count() {
        let mut dsu = Dsu::new(6);
        assert_eq!(dsu.set_count(), 6);

        dsu.union(0, 1);
        dsu.union(2, 3);
        assert_eq!(dsu.set_count(), 4);

        dsu.union(0, 2);
        assert_eq!(dsu.set_count(), 3);
    }

    // -----------------------------------------------------------------------
    // Chain island tests
    // -----------------------------------------------------------------------
    #[test]
    fn test_chain_island() {
        let mut mgr = IslandManager::new();
        let handles: Vec<BodyHandle> = (0..5).map(handle).collect();
        // Chain: 0-1-2-3-4
        let contacts: Vec<(BodyHandle, BodyHandle)> =
            (0..4).map(|i| (handles[i], handles[i + 1])).collect();
        mgr.build_islands(&handles, &contacts);

        assert_eq!(mgr.islands().len(), 1, "chain should form one island");
        assert_eq!(mgr.islands()[0].body_handles.len(), 5);
    }

    #[test]
    fn test_rebuild_islands_preserves_sleep_timer() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[(h0, h1)]);

        // Advance sleep timer
        let ke = |_: BodyHandle| 0.0_f64;
        mgr.update_sleep(0.1, 1.0, 10.0, ke);
        let timer_before = mgr.islands()[0].sleep_timer;
        assert!(timer_before > 0.0);

        // Rebuild with same contacts — timer should be preserved
        mgr.build_islands(&[h0, h1], &[(h0, h1)]);
        let timer_after = mgr.islands()[0].sleep_timer;
        assert!(
            (timer_after - timer_before).abs() < 1e-12,
            "sleep timer should be preserved across rebuild"
        );
    }

    // ── update_sleep_hysteresis tests ─────────────────────────────────────

    #[test]
    fn test_sleep_hysteresis_stays_awake_high_energy() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        mgr.build_islands(&[h0], &[]);
        let mut counters = vec![0u32];
        // Energy always above sleep threshold → never sleeps
        mgr.update_sleep_hysteresis(|_| 10.0, 1.0, 2.0, 5, &mut counters);
        assert!(
            !mgr.islands()[0].is_sleeping,
            "high-energy island should not sleep"
        );
    }

    #[test]
    fn test_sleep_hysteresis_sleeps_after_enough_frames() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        mgr.build_islands(&[h0], &[]);
        let mut counters = vec![0u32];
        // Apply 3 frames below threshold (need 3 to sleep)
        for _ in 0..3 {
            mgr.update_sleep_hysteresis(|_| 0.0, 1.0, 2.0, 3, &mut counters);
        }
        assert!(
            mgr.islands()[0].is_sleeping,
            "island should sleep after 3 low-energy frames"
        );
    }

    #[test]
    fn test_sleep_hysteresis_wakes_on_high_energy() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        mgr.build_islands(&[h0], &[]);
        // Force island to sleep first
        mgr.force_sleep(0);
        assert!(mgr.islands()[0].is_sleeping);

        let mut counters = vec![0u32];
        // Apply high energy → should wake
        let changed = mgr.update_sleep_hysteresis(|_| 5.0, 1.0, 2.0, 3, &mut counters);
        assert!(
            !mgr.islands()[0].is_sleeping,
            "island should wake on high energy"
        );
        assert_eq!(changed, 1);
    }

    // ── sorted_island_indices_by_energy tests ─────────────────────────────

    #[test]
    fn test_sorted_islands_highest_energy_first() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        let h2 = handle(2);
        mgr.build_islands(&[h0, h1, h2], &[]);
        // ke by island index: 0→1.0, 1→3.0, 2→2.0 → sorted: [1, 2, 0]
        let ke_map = [1.0_f64, 3.0, 2.0];
        let sorted = mgr.sorted_island_indices_by_energy(|i| ke_map[i]);
        assert_eq!(sorted.len(), 3);
        // First index should have highest energy
        let first_ke = ke_map[sorted[0]];
        let last_ke = ke_map[sorted[sorted.len() - 1]];
        assert!(
            first_ke >= last_ke,
            "first island should have highest energy"
        );
    }

    // ── island_size_histogram tests ───────────────────────────────────────

    #[test]
    fn test_island_size_histogram_single_islands() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[]); // 2 singleton islands
        let hist = mgr.island_size_histogram(3);
        assert_eq!(hist[0], 2, "should have 2 singleton islands");
    }

    #[test]
    fn test_island_size_histogram_pair_island() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[(h0, h1)]); // 1 island of size 2
        let hist = mgr.island_size_histogram(3);
        assert_eq!(hist[1], 1, "should have 1 island of size 2");
    }

    // ── cross_island_contacts tests ───────────────────────────────────────

    #[test]
    fn test_cross_island_contacts_detected() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[]); // separate islands
        // Contact between different islands
        let contacts = vec![(h0, h1)];
        let cross = mgr.cross_island_contacts(&contacts);
        assert_eq!(cross.len(), 1, "should detect 1 cross-island contact");
    }

    #[test]
    fn test_cross_island_contacts_none_within_island() {
        let mut mgr = IslandManager::new();
        let h0 = handle(0);
        let h1 = handle(1);
        mgr.build_islands(&[h0, h1], &[(h0, h1)]); // same island
        let contacts = vec![(h0, h1)];
        let cross = mgr.cross_island_contacts(&contacts);
        assert!(
            cross.is_empty(),
            "within-island contacts are not cross-island"
        );
    }

    // ── classify_islands_by_energy tests ─────────────────────────────────

    #[test]
    fn test_classify_islands_hot_warm_cold() {
        let mut mgr = IslandManager::new();
        for i in 0..3u32 {
            let h = handle(i);
            mgr.build_islands(&[h], &[]);
        }
        // Rebuild all 3 separately (3 singletons)
        let h0 = handle(0);
        let h1 = handle(1);
        let h2 = handle(2);
        mgr.build_islands(&[h0, h1, h2], &[]);
        // ke: island 0 → 100 (hot), island 1 → 5 (warm), island 2 → 0 (cold)
        let ke_map = [100.0_f64, 5.0, 0.0];
        let (hot, warm, cold) = mgr.classify_islands_by_energy(|i| ke_map[i], 50.0, 1.0);
        assert_eq!(hot.len(), 1, "one hot island");
        assert_eq!(warm.len(), 1, "one warm island");
        assert_eq!(cold.len(), 1, "one cold island");
    }
}
