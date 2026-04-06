//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{BodyId, ConstraintId, SLEEP_VELOCITY_THRESHOLD};

/// Cache of contact manifold entries for all active body pairs.
#[derive(Debug, Clone, Default)]
pub struct ContactManifoldCache {
    /// All manifold entries.
    pub entries: Vec<ContactManifoldEntry>,
}
impl ContactManifoldCache {
    /// Create an empty cache.
    pub fn new() -> Self {
        Self::default()
    }
    /// Find the manifold entry for a body pair, if present.
    pub fn get(&self, body_a: BodyId, body_b: BodyId) -> Option<&ContactManifoldEntry> {
        self.entries.iter().find(|e| {
            (e.body_a == body_a && e.body_b == body_b) || (e.body_a == body_b && e.body_b == body_a)
        })
    }
    /// Find the mutable manifold entry for a body pair, if present.
    pub fn get_mut(&mut self, body_a: BodyId, body_b: BodyId) -> Option<&mut ContactManifoldEntry> {
        self.entries.iter_mut().find(|e| {
            (e.body_a == body_a && e.body_b == body_b) || (e.body_a == body_b && e.body_b == body_a)
        })
    }
    /// Insert a new manifold entry or replace the existing one for the same pair.
    pub fn insert(&mut self, entry: ContactManifoldEntry) {
        if let Some(existing) = self.get_mut(entry.body_a, entry.body_b) {
            *existing = entry;
        } else {
            self.entries.push(entry);
        }
    }
    /// Remove the manifold entry for a body pair.
    pub fn remove(&mut self, body_a: BodyId, body_b: BodyId) {
        self.entries.retain(|e| {
            !((e.body_a == body_a && e.body_b == body_b)
                || (e.body_a == body_b && e.body_b == body_a))
        });
    }
    /// Advance all entries by one frame.  Remove entries older than `max_age` frames.
    pub fn advance_all(&mut self, decay: f64, max_age: u32) {
        for e in &mut self.entries {
            e.advance_frame(decay);
        }
        self.entries.retain(|e| e.frame <= max_age);
    }
    /// Number of cached manifolds.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
}
/// Energy report for a single island.
#[derive(Debug, Clone)]
pub struct IslandEnergyReport {
    /// Island index.
    pub island_idx: usize,
    /// Total kinetic energy (translational + rotational) of all dynamic bodies.
    pub kinetic_energy: f64,
    /// Number of dynamic bodies.
    pub n_dynamic_bodies: usize,
    /// Whether the island qualifies for sleeping.
    pub can_sleep: bool,
}
impl IslandEnergyReport {
    /// Compute the energy report for a given island.
    pub fn compute(island_idx: usize, island: &Island, bodies: &[IslandBody]) -> Self {
        let mut ke = 0.0_f64;
        let mut n_dynamic = 0usize;
        for &bid in &island.bodies {
            if let Some(body) = bodies.iter().find(|b| b.id == bid) {
                if body.is_static || body.inv_mass == 0.0 {
                    continue;
                }
                n_dynamic += 1;
                let lv = body.linear_vel;
                let av = body.angular_vel;
                let mass = 1.0 / body.inv_mass;
                ke += 0.5 * mass * (lv[0] * lv[0] + lv[1] * lv[1] + lv[2] * lv[2]);
                ke += 0.5 * mass * 0.1 * (av[0] * av[0] + av[1] * av[1] + av[2] * av[2]);
            }
        }
        let can_sleep = island.can_sleep(bodies);
        IslandEnergyReport {
            island_idx,
            kinetic_energy: ke,
            n_dynamic_bodies: n_dynamic,
            can_sleep,
        }
    }
}
/// Simple PGS (Projected Gauss-Seidel) solver that iterates over one island's
/// constraints and applies corrective impulses to the participating bodies.
#[derive(Debug, Clone)]
pub struct SequentialImpulseSolver {
    /// Number of PGS iterations per call to [`Self::solve_island`].
    pub iterations: usize,
}
impl SequentialImpulseSolver {
    /// Create a solver with the given iteration count.
    pub fn new(iterations: usize) -> Self {
        Self { iterations }
    }
    /// Run PGS on the given bodies and constraints for one time step.
    ///
    /// Only `Contact` constraints are actively solved in this implementation;
    /// `Joint` and `Distance` constraints are accepted but currently act as
    /// structural hints only (they still influence island membership).
    pub fn solve_island(&self, bodies: &mut [IslandBody], constraints: &[IslandConstraint]) {
        for _ in 0..self.iterations {
            for c in constraints {
                let (ia, ib) = match (
                    bodies.iter().position(|b| b.id == c.body_a),
                    bodies.iter().position(|b| b.id == c.body_b),
                ) {
                    (Some(a), Some(b)) => (a, b),
                    _ => continue,
                };
                match &c.constraint_type {
                    ConstraintType::Contact { normal, depth } => {
                        let (ba, bb) = if ia < ib {
                            let (left, right) = bodies.split_at_mut(ib);
                            (&mut left[ia], &mut right[0])
                        } else {
                            let (left, right) = bodies.split_at_mut(ia);
                            (&mut right[0], &mut left[ib])
                        };
                        solve_contact_constraint(ba, bb, *normal, *depth, 0.0);
                    }
                    ConstraintType::Joint { .. } | ConstraintType::Distance { .. } => {}
                }
            }
        }
    }
}
/// Detects and applies island merges when new constraints bridge two previously
/// independent islands.
///
/// The merger performs a single-pass scan of new constraints and updates the
/// provided `IslandManager` in place.
pub struct IslandMerger;
impl IslandMerger {
    /// Scan `new_constraints` for cross-island constraints and merge the
    /// affected islands.
    ///
    /// Returns the number of merges performed.
    pub fn merge_from_constraints(
        manager: &mut IslandManager,
        new_constraints: &[IslandConstraint],
    ) -> usize {
        let mut merges = 0;
        for c in new_constraints {
            let ia = manager.bodies.iter().position(|b| b.id == c.body_a);
            let ib = manager.bodies.iter().position(|b| b.id == c.body_b);
            if let (Some(ia), Some(ib)) = (ia, ib) {
                let island_a = manager
                    .islands
                    .iter()
                    .position(|isl| isl.bodies.contains(&manager.bodies[ia].id));
                let island_b = manager
                    .islands
                    .iter()
                    .position(|isl| isl.bodies.contains(&manager.bodies[ib].id));
                if let (Some(ia_idx), Some(ib_idx)) = (island_a, island_b)
                    && ia_idx != ib_idx
                {
                    let bodies_b: Vec<BodyId> = manager.islands[ib_idx].bodies.clone();
                    let constraints_b: Vec<ConstraintId> =
                        manager.islands[ib_idx].constraints.clone();
                    manager.islands[ia_idx].bodies.extend(bodies_b);
                    manager.islands[ia_idx].constraints.extend(constraints_b);
                    manager.islands[ib_idx].bodies.clear();
                    manager.islands[ib_idx].constraints.clear();
                    merges += 1;
                }
            }
        }
        manager.islands.retain(|isl| !isl.bodies.is_empty());
        merges
    }
}
/// A parallel island solver that dispatches each island to an independent solve.
///
/// In a real engine these would run on separate threads; here we simulate
/// the interface sequentially.
pub struct ParallelIslandSolver {
    /// Base PGS solver.
    pub base_solver: SequentialImpulseSolver,
    /// Maximum sub-steps per island.
    pub max_substeps: usize,
}
impl ParallelIslandSolver {
    /// Create a new parallel island solver.
    pub fn new(iterations: usize, max_substeps: usize) -> Self {
        Self {
            base_solver: SequentialImpulseSolver::new(iterations),
            max_substeps,
        }
    }
    /// Solve all active (non-sleeping) islands.
    ///
    /// Returns the number of islands processed.
    pub fn solve_all_islands(&self, manager: &mut IslandManager) -> usize {
        let body_to_island = build_body_to_island_map(&manager.islands);
        let substep_plan = adaptive_substeps(&manager.islands, &manager.bodies, self.max_substeps);
        let mut islands_processed = 0;
        for plan in &substep_plan {
            if plan.n_substeps == 0 {
                continue;
            }
            let island = &manager.islands[plan.island_idx];
            let constraint_ids: Vec<ConstraintId> = island.constraints.clone();
            let body_ids: Vec<BodyId> = island.bodies.clone();
            let constraints: Vec<IslandConstraint> = constraint_ids
                .iter()
                .filter_map(|&cid| manager.constraints.iter().find(|c| c.id == cid).cloned())
                .collect();
            let body_indices: Vec<usize> = body_ids
                .iter()
                .filter_map(|&bid| manager.bodies.iter().position(|b| b.id == bid))
                .collect();
            for _ in 0..plan.n_substeps {
                for c in &constraints {
                    let ia = manager.bodies.iter().position(|b| b.id == c.body_a);
                    let ib = manager.bodies.iter().position(|b| b.id == c.body_b);
                    if let (Some(ia), Some(ib)) = (ia, ib)
                        && let ConstraintType::Contact { normal, depth } = c.constraint_type
                    {
                        let (ba, bb) = if ia < ib {
                            let (left, right) = manager.bodies.split_at_mut(ib);
                            (&mut left[ia], &mut right[0])
                        } else {
                            let (left, right) = manager.bodies.split_at_mut(ia);
                            (&mut right[0], &mut left[ib])
                        };
                        solve_contact_constraint(ba, bb, normal, depth, 0.0);
                    }
                }
            }
            islands_processed += 1;
            let _ = body_to_island.get(&0);
            let _ = body_indices;
        }
        islands_processed
    }
}
/// High-level solver that wraps `IslandManager` with energy tracking, small-island
/// merging, and energy-based sleep decisions.
#[derive(Debug, Default)]
pub struct IslandSolver {
    /// The underlying island manager.
    pub manager: IslandManager,
    /// Energy threshold below which an island is considered for sleeping.
    pub sleep_energy_threshold: f64,
    /// Minimum island size (in bodies) to keep separate; smaller islands are merged.
    pub min_island_size: usize,
}
impl IslandSolver {
    /// Create a new island solver with default thresholds.
    pub fn new() -> Self {
        Self {
            manager: IslandManager::new(),
            sleep_energy_threshold: 1e-4,
            min_island_size: 2,
        }
    }
    /// Create a new island solver with explicit thresholds.
    pub fn with_thresholds(sleep_energy_threshold: f64, min_island_size: usize) -> Self {
        Self {
            manager: IslandManager::new(),
            sleep_energy_threshold,
            min_island_size,
        }
    }
    /// Compute the total kinetic energy for a specific island by index.
    ///
    /// Returns `None` if the island index is out of range.
    ///
    /// Energy formula: KE = Σ ½ m v² + Σ ½ m * 0.1 * ω²  (rotational proxy).
    pub fn compute_island_energy(&self, island_idx: usize) -> Option<f64> {
        let island = self.manager.islands.get(island_idx)?;
        let ke = island
            .bodies
            .iter()
            .filter_map(|&bid| self.manager.bodies.iter().find(|b| b.id == bid))
            .filter(|b| !b.is_static && b.inv_mass > 0.0)
            .map(|b| {
                let mass = 1.0 / b.inv_mass;
                let lv = b.linear_vel;
                let av = b.angular_vel;
                0.5 * mass * (lv[0] * lv[0] + lv[1] * lv[1] + lv[2] * lv[2])
                    + 0.5 * mass * 0.1 * (av[0] * av[0] + av[1] * av[1] + av[2] * av[2])
            })
            .sum();
        Some(ke)
    }
    /// Compute energy for all islands and return a `Vec<(island_idx, energy)>`.
    pub fn compute_all_island_energies(&self) -> Vec<(usize, f64)> {
        (0..self.manager.islands.len())
            .filter_map(|i| self.compute_island_energy(i).map(|e| (i, e)))
            .collect()
    }
    /// Merge all islands that have fewer than `self.min_island_size` bodies into their
    /// nearest neighbour island (the one sharing the most constraints).
    ///
    /// Returns the number of islands merged.
    pub fn merge_small_islands(&mut self) -> usize {
        let threshold = self.min_island_size;
        let mut merges = 0;
        loop {
            let small_idx = self
                .manager
                .islands
                .iter()
                .position(|isl| isl.bodies.len() < threshold && !isl.bodies.is_empty());
            let small_idx = match small_idx {
                Some(i) => i,
                None => break,
            };
            let small_bodies: std::collections::HashSet<BodyId> = self.manager.islands[small_idx]
                .bodies
                .iter()
                .copied()
                .collect();
            let target_idx = (0..self.manager.islands.len())
                .filter(|&i| i != small_idx && !self.manager.islands[i].bodies.is_empty())
                .max_by_key(|&i| {
                    let other_bodies: std::collections::HashSet<BodyId> =
                        self.manager.islands[i].bodies.iter().copied().collect();
                    self.manager
                        .constraints
                        .iter()
                        .filter(|c| {
                            (small_bodies.contains(&c.body_a) && other_bodies.contains(&c.body_b))
                                || (small_bodies.contains(&c.body_b)
                                    && other_bodies.contains(&c.body_a))
                        })
                        .count()
                });
            match target_idx {
                Some(t) => {
                    let bodies: Vec<BodyId> = self.manager.islands[small_idx].bodies.clone();
                    let constraints: Vec<ConstraintId> =
                        self.manager.islands[small_idx].constraints.clone();
                    self.manager.islands[t].bodies.extend(bodies);
                    self.manager.islands[t].constraints.extend(constraints);
                    self.manager.islands[small_idx].bodies.clear();
                    self.manager.islands[small_idx].constraints.clear();
                    merges += 1;
                }
                None => break,
            }
        }
        self.manager.islands.retain(|isl| !isl.bodies.is_empty());
        merges
    }
    /// Compute which islands qualify for sleeping based on energy thresholds.
    ///
    /// An island qualifies for sleep when its total kinetic energy drops below
    /// `self.sleep_energy_threshold`.
    ///
    /// Returns a list of island indices that should be put to sleep.
    pub fn compute_sleeping_criteria(&self) -> Vec<usize> {
        (0..self.manager.islands.len())
            .filter(|&i| {
                let island = &self.manager.islands[i];
                if island.is_sleeping {
                    return false;
                }
                if island.bodies.is_empty() {
                    return false;
                }
                let can_sleep_vel = island.can_sleep(&self.manager.bodies);
                let energy = self.compute_island_energy(i).unwrap_or(0.0);
                let energy_ok = energy < self.sleep_energy_threshold;
                can_sleep_vel && energy_ok
            })
            .collect()
    }
    /// Mark the islands returned by `compute_sleeping_criteria` as sleeping.
    ///
    /// Returns the number of islands put to sleep.
    pub fn apply_sleeping(&mut self) -> usize {
        let to_sleep = self.compute_sleeping_criteria();
        let n = to_sleep.len();
        for i in to_sleep {
            if let Some(island) = self.manager.islands.get_mut(i) {
                island.is_sleeping = true;
            }
        }
        n
    }
    /// Wake all sleeping islands (e.g., after a large impulse is detected).
    pub fn wake_all_islands(&mut self) {
        for island in &mut self.manager.islands {
            island.is_sleeping = false;
        }
    }
}
/// Union-Find (Disjoint Set Union) with path compression and union by rank.
#[derive(Debug)]
pub struct UnionFind {
    /// Parent array; `parent[i] == i` for roots.
    pub parent: Vec<usize>,
    /// Rank array used for union by rank.
    pub rank: Vec<usize>,
}
impl UnionFind {
    /// Construct a new structure with `n` singleton components.
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    /// Find the representative of `x` with path compression.
    pub fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut cur = x;
        while cur != root {
            let next = self.parent[cur];
            self.parent[cur] = root;
            cur = next;
        }
        root
    }
    /// Merge the components containing `x` and `y` (union by rank).
    pub fn union(&mut self, x: usize, y: usize) {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx == ry {
            return;
        }
        match self.rank[rx].cmp(&self.rank[ry]) {
            std::cmp::Ordering::Less => self.parent[rx] = ry,
            std::cmp::Ordering::Greater => self.parent[ry] = rx,
            std::cmp::Ordering::Equal => {
                self.parent[ry] = rx;
                self.rank[rx] += 1;
            }
        }
    }
    /// Returns `true` when `x` and `y` are in the same component.
    pub fn same_component(&mut self, x: usize, y: usize) -> bool {
        self.find(x) == self.find(y)
    }
}
/// Discriminated union of supported constraint flavours.
#[derive(Debug, Clone)]
pub enum ConstraintType {
    /// Frictionless contact constraint.
    Contact {
        /// Contact normal (world-space, pointing from B to A).
        normal: [f64; 3],
        /// Penetration depth (positive = overlapping).
        depth: f64,
    },
    /// Rigid joint that enforces a fixed relative pose between two anchors.
    Joint {
        /// Anchor point in body A's local space.
        local_anchor_a: [f64; 3],
        /// Anchor point in body B's local space.
        local_anchor_b: [f64; 3],
    },
    /// Distance constraint that keeps two bodies at a target separation.
    Distance {
        /// Desired distance between the two body origins.
        target_length: f64,
    },
}
/// Detects islands that have become disconnected (need splitting).
///
/// An island needs splitting when some of its constraints have been removed
/// and the remaining connectivity does not cover all bodies.
pub struct IslandSplitDetector;
impl IslandSplitDetector {
    /// Detect islands that need splitting given an active constraint set.
    ///
    /// Returns the indices of islands where not all bodies are connected.
    pub fn detect(manager: &IslandManager, active_constraints: &[IslandConstraint]) -> Vec<usize> {
        let mut split_candidates = Vec::new();
        for (island_idx, island) in manager.islands.iter().enumerate() {
            if island.bodies.len() <= 1 {
                continue;
            }
            let body_set: std::collections::HashSet<BodyId> =
                island.bodies.iter().copied().collect();
            let mut adjacency: std::collections::HashMap<BodyId, Vec<BodyId>> =
                island.bodies.iter().map(|&b| (b, Vec::new())).collect();
            for c in active_constraints {
                if body_set.contains(&c.body_a) && body_set.contains(&c.body_b) {
                    adjacency.entry(c.body_a).or_default().push(c.body_b);
                    adjacency.entry(c.body_b).or_default().push(c.body_a);
                }
            }
            let start = island.bodies[0];
            let mut visited = std::collections::HashSet::new();
            let mut queue = std::collections::VecDeque::new();
            visited.insert(start);
            queue.push_back(start);
            while let Some(current) = queue.pop_front() {
                if let Some(neighbors) = adjacency.get(&current) {
                    for &nb in neighbors {
                        if visited.insert(nb) {
                            queue.push_back(nb);
                        }
                    }
                }
            }
            if visited.len() < island.bodies.len() {
                split_candidates.push(island_idx);
            }
        }
        split_candidates
    }
}
/// A cached contact manifold for one body-pair within an island.
///
/// Stores up to 4 contact points (classic manifold reduction limit) with
/// accumulated impulses for warm-starting.
#[derive(Debug, Clone)]
pub struct ContactManifoldEntry {
    /// Body A identifier.
    pub body_a: BodyId,
    /// Body B identifier.
    pub body_b: BodyId,
    /// Contact normal (world-space, from B toward A).
    pub normal: [f64; 3],
    /// Contact points (world-space), up to 4.
    pub points: [[f64; 3]; 4],
    /// Penetration depths at each contact point.
    pub depths: [f64; 4],
    /// Accumulated normal impulse per contact point (warm-start).
    pub normal_lambdas: [f64; 4],
    /// Number of active contact points.
    pub n_points: usize,
    /// Frame counter for freshness tracking.
    pub frame: u32,
}
impl ContactManifoldEntry {
    /// Create a new manifold entry for a body pair.
    pub fn new(body_a: BodyId, body_b: BodyId, normal: [f64; 3]) -> Self {
        Self {
            body_a,
            body_b,
            normal,
            points: [[0.0; 3]; 4],
            depths: [0.0; 4],
            normal_lambdas: [0.0; 4],
            n_points: 0,
            frame: 0,
        }
    }
    /// Add a contact point.  Returns `false` when the manifold is full (> 4 pts).
    pub fn add_point(&mut self, point: [f64; 3], depth: f64) -> bool {
        if self.n_points >= 4 {
            return false;
        }
        self.points[self.n_points] = point;
        self.depths[self.n_points] = depth;
        self.normal_lambdas[self.n_points] = 0.0;
        self.n_points += 1;
        true
    }
    /// Get the warm-start impulse for a contact point.
    pub fn warm_start_impulse(&self, point_idx: usize, factor: f64) -> f64 {
        if point_idx < self.n_points {
            self.normal_lambdas[point_idx] * factor.clamp(0.0, 1.0)
        } else {
            0.0
        }
    }
    /// Update the accumulated impulse for a contact point.
    pub fn update_lambda(&mut self, point_idx: usize, lambda: f64) {
        if point_idx < self.n_points {
            self.normal_lambdas[point_idx] = lambda;
        }
    }
    /// Total accumulated normal impulse across all contact points.
    pub fn total_normal_impulse(&self) -> f64 {
        self.normal_lambdas[..self.n_points].iter().sum()
    }
    /// Advance to a new frame (increments frame counter, optionally decays impulses).
    pub fn advance_frame(&mut self, decay: f64) {
        self.frame += 1;
        for lam in self.normal_lambdas.iter_mut().take(self.n_points) {
            *lam *= decay.clamp(0.0, 1.0);
        }
    }
}
/// Manages the complete set of bodies and constraints, and partitions them into
/// simulation islands (connected components) for parallel solving.
#[derive(Debug, Default)]
pub struct IslandManager {
    /// All bodies registered with the manager.
    pub bodies: Vec<IslandBody>,
    /// All constraints registered with the manager.
    pub constraints: Vec<IslandConstraint>,
    /// Current set of islands built by [`Self::build_islands`].
    pub islands: Vec<Island>,
}
impl IslandManager {
    /// Create a new, empty island manager.
    pub fn new() -> Self {
        Self::default()
    }
    /// Register a body and return its index in `self.bodies`.
    pub fn add_body(&mut self, body: IslandBody) -> usize {
        let idx = self.bodies.len();
        self.bodies.push(body);
        idx
    }
    /// Register a constraint and return its index in `self.constraints`.
    pub fn add_constraint(&mut self, constraint: IslandConstraint) -> usize {
        let idx = self.constraints.len();
        self.constraints.push(constraint);
        idx
    }
    /// Partition bodies and constraints into connected-component islands.
    ///
    /// Uses Union-Find over the constraint graph: each constraint merges its two
    /// participant bodies into the same component.  Bodies not referenced by any
    /// constraint become singleton islands.
    pub fn build_islands(&mut self) {
        let n = self.bodies.len();
        if n == 0 {
            self.islands.clear();
            return;
        }
        let body_index =
            |bid: BodyId| -> Option<usize> { self.bodies.iter().position(|b| b.id == bid) };
        let mut uf = UnionFind::new(n);
        for c in &self.constraints {
            if let (Some(ia), Some(ib)) = (body_index(c.body_a), body_index(c.body_b)) {
                uf.union(ia, ib);
            }
        }
        let mut root_to_local: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for i in 0..n {
            let root = uf.find(i);
            root_to_local.entry(root).or_default().push(i);
        }
        let mut new_islands: Vec<Island> = Vec::with_capacity(root_to_local.len());
        let mut root_to_island: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        for (root, local_indices) in &root_to_local {
            let island_idx = new_islands.len();
            root_to_island.insert(*root, island_idx);
            let body_ids: Vec<BodyId> = local_indices.iter().map(|&i| self.bodies[i].id).collect();
            new_islands.push(Island {
                bodies: body_ids,
                constraints: Vec::new(),
                is_sleeping: false,
            });
        }
        for c in &self.constraints {
            if let Some(ia) = body_index(c.body_a) {
                let root = uf.find(ia);
                if let Some(&island_idx) = root_to_island.get(&root) {
                    new_islands[island_idx].constraints.push(c.id);
                }
            }
        }
        for island in &mut new_islands {
            island.is_sleeping = island.can_sleep(&self.bodies);
        }
        self.islands = new_islands;
    }
    /// Total number of islands.
    pub fn island_count(&self) -> usize {
        self.islands.len()
    }
    /// Number of islands that are currently sleeping.
    pub fn sleeping_island_count(&self) -> usize {
        self.islands.iter().filter(|i| i.is_sleeping).count()
    }
    /// Total number of bodies across all islands.
    pub fn total_bodies(&self) -> usize {
        self.bodies.len()
    }
    /// Wake up an entire island: mark all of its bodies as awake.
    pub fn wake_island(&mut self, island_idx: usize) {
        if island_idx >= self.islands.len() {
            return;
        }
        let body_ids: Vec<BodyId> = self.islands[island_idx].bodies.clone();
        for bid in body_ids {
            if let Some(body) = self.bodies.iter_mut().find(|b| b.id == bid) {
                body.is_sleeping = false;
            }
        }
        self.islands[island_idx].is_sleeping = false;
    }
}
/// A body as seen by the island solver.
#[derive(Debug, Clone)]
pub struct IslandBody {
    /// Unique body identifier.
    pub id: BodyId,
    /// Inverse mass (0 for static/kinematic bodies).
    pub inv_mass: f64,
    /// Linear velocity (world-space, m/s).
    pub linear_vel: [f64; 3],
    /// Angular velocity (world-space, rad/s).
    pub angular_vel: [f64; 3],
    /// World-space position.
    pub position: [f64; 3],
    /// World-space orientation as a unit quaternion `[x, y, z, w]`.
    pub orientation: [f64; 4],
    /// True when this body has infinite mass (never moved by the solver).
    pub is_static: bool,
    /// True when the body is currently sleeping (excluded from integration).
    pub is_sleeping: bool,
}
impl IslandBody {
    /// Returns true when the body's velocity is below the sleep threshold.
    fn is_slow(&self) -> bool {
        let lv = &self.linear_vel;
        let av = &self.angular_vel;
        let lv2 = lv[0] * lv[0] + lv[1] * lv[1] + lv[2] * lv[2];
        let av2 = av[0] * av[0] + av[1] * av[1] + av[2] * av[2];
        lv2 < SLEEP_VELOCITY_THRESHOLD * SLEEP_VELOCITY_THRESHOLD
            && av2 < SLEEP_VELOCITY_THRESHOLD * SLEEP_VELOCITY_THRESHOLD
    }
}
/// An edge in the island dependency graph.
///
/// Islands are "dependent" when a static body is shared between them,
/// or when a joint constraint spans two otherwise disconnected components.
#[derive(Debug, Clone)]
pub struct IslandDependency {
    /// First island index.
    pub island_a: usize,
    /// Second island index.
    pub island_b: usize,
}
/// A warm-start cache that maps `(BodyId, BodyId)` → accumulated impulse.
///
/// Provides O(1) lookup and insertion for frequently-used constraint pairs.
#[derive(Debug, Clone, Default)]
pub struct ImpulseWarmStartCache {
    /// Key: sorted (min_id, max_id); Value: accumulated impulse magnitude.
    pub(super) map: std::collections::HashMap<(BodyId, BodyId), f64>,
    /// Maximum number of entries.
    pub capacity: usize,
}
impl ImpulseWarmStartCache {
    /// Create a new warm-start cache with a given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            map: std::collections::HashMap::with_capacity(capacity),
            capacity,
        }
    }
    fn key(a: BodyId, b: BodyId) -> (BodyId, BodyId) {
        if a <= b { (a, b) } else { (b, a) }
    }
    /// Store the accumulated impulse for a body pair.
    pub fn store(&mut self, body_a: BodyId, body_b: BodyId, lambda: f64) {
        if self.map.len() >= self.capacity && !self.map.contains_key(&Self::key(body_a, body_b)) {
            return;
        }
        self.map.insert(Self::key(body_a, body_b), lambda);
    }
    /// Retrieve the warm-start impulse for a body pair (or 0 if not cached).
    pub fn retrieve(&self, body_a: BodyId, body_b: BodyId) -> f64 {
        self.map
            .get(&Self::key(body_a, body_b))
            .copied()
            .unwrap_or(0.0)
    }
    /// Decay all stored impulses by `factor`.
    pub fn decay_all(&mut self, factor: f64) {
        for v in self.map.values_mut() {
            *v *= factor;
        }
    }
    /// Remove entries below a minimum threshold (clean up negligible impulses).
    pub fn prune(&mut self, min_lambda: f64) {
        self.map.retain(|_, v| v.abs() >= min_lambda);
    }
    /// Number of cached entries.
    pub fn len(&self) -> usize {
        self.map.len()
    }
    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
    /// Clear all entries.
    pub fn clear(&mut self) {
        self.map.clear();
    }
}
/// Residual (violation) information for a single constraint.
///
/// Useful for debugging, convergence monitoring, and adaptive sub-stepping.
#[derive(Debug, Clone)]
pub struct ConstraintResidual {
    /// Constraint identifier.
    pub constraint_id: ConstraintId,
    /// Residual magnitude (velocity-level violation).
    pub residual: f64,
    /// Whether the constraint is currently active (non-zero impulse).
    pub is_active: bool,
}
/// A connected component of bodies and the constraints linking them.
#[derive(Debug, Clone, Default)]
pub struct Island {
    /// Indices into `IslandManager::bodies` that belong to this island.
    pub bodies: Vec<BodyId>,
    /// Indices into `IslandManager::constraints` that belong to this island.
    pub constraints: Vec<ConstraintId>,
    /// Whether the entire island is currently sleeping.
    pub is_sleeping: bool,
}
impl Island {
    /// Number of bodies in this island.
    pub fn size(&self) -> usize {
        self.bodies.len()
    }
    /// Returns `true` when every dynamic body in the island is slow enough to sleep.
    ///
    /// A body is considered dynamic when `inv_mass > 0` and `is_static == false`.
    pub fn can_sleep(&self, all_bodies: &[IslandBody]) -> bool {
        for &bid in &self.bodies {
            if let Some(body) = all_bodies.iter().find(|b| b.id == bid)
                && !body.is_static
                && body.inv_mass > 0.0
                && !body.is_slow()
            {
                return false;
            }
        }
        true
    }
}
/// A constraint between two bodies, as seen by the island solver.
#[derive(Debug, Clone)]
pub struct IslandConstraint {
    /// Unique constraint identifier.
    pub id: ConstraintId,
    /// First body participating in this constraint.
    pub body_a: BodyId,
    /// Second body participating in this constraint.
    pub body_b: BodyId,
    /// Constraint flavour and its parameters.
    pub constraint_type: ConstraintType,
    /// Accumulated impulse magnitude from the previous frame (warm start).
    pub lambda: f64,
}
/// Sub-step budget for an island (number of solver iterations).
#[derive(Debug, Clone)]
pub struct IslandSubStep {
    /// Island index.
    pub island_idx: usize,
    /// Number of sub-steps (iterations) allocated.
    pub n_substeps: usize,
}
/// Aggregated statistics about the current island partition.
#[derive(Debug, Clone, Default)]
pub struct IslandStatistics {
    /// Total number of bodies across all islands.
    pub total_bodies: usize,
    /// Total number of constraints across all islands.
    pub total_constraints: usize,
    /// Total number of islands.
    pub n_islands: usize,
    /// Number of islands that are currently sleeping.
    pub n_sleeping_islands: usize,
    /// Size (body count) of the largest island.
    pub largest_island: usize,
}
impl IslandStatistics {
    /// Compute statistics from a live [`IslandManager`].
    pub fn compute(manager: &IslandManager) -> Self {
        let n_sleeping_islands = manager.sleeping_island_count();
        let largest_island = manager
            .islands
            .iter()
            .map(|i| i.bodies.len())
            .max()
            .unwrap_or(0);
        Self {
            total_bodies: manager.bodies.len(),
            total_constraints: manager.constraints.len(),
            n_islands: manager.islands.len(),
            n_sleeping_islands,
            largest_island,
        }
    }
}
/// Constraint sorting key for deterministic island processing.
#[derive(Debug, Clone)]
pub struct SortedConstraint {
    /// Island index this constraint belongs to.
    pub island_idx: usize,
    /// The constraint itself.
    pub constraint: IslandConstraint,
}
/// Per-island solver budget (iterations, time, etc.) used for adaptive control.
#[derive(Debug, Clone)]
pub struct IslandBudget {
    /// Island index.
    pub island_idx: usize,
    /// Allocated PGS iterations for this island.
    pub iterations: usize,
    /// Allocated time budget in microseconds (advisory; not enforced here).
    pub time_us: u64,
    /// Priority score (higher = solve first).
    pub priority: f64,
}
impl IslandBudget {
    /// Create a new budget entry.
    pub fn new(island_idx: usize, iterations: usize, time_us: u64, priority: f64) -> Self {
        Self {
            island_idx,
            iterations,
            time_us,
            priority,
        }
    }
}
