//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::contact::ContactConstraint;
use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::Real;
use oxiphysics_rigid::{BodyState, BodyType, RigidBodySet};
use std::collections::HashMap;

/// A work item for solving one island's constraints.
///
/// Contains all indices needed to solve the island independently.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct IslandWorkItem {
    /// Index of the island in the islands array.
    pub island_index: usize,
    /// Body handles to process.
    pub body_handles: Vec<BodyHandle>,
    /// Contact constraint indices.
    pub contact_indices: Vec<usize>,
    /// Joint constraint indices.
    pub joint_indices: Vec<usize>,
    /// Estimated cost for load balancing (bodies + constraints).
    pub estimated_cost: usize,
}
/// Accumulates per-body velocity history for sleep detection.
///
/// A body is eligible to sleep when its accumulated velocity falls below
/// the threshold for `required_steps` consecutive steps.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct VelocityAccumulator {
    /// Summed velocity magnitudes per body (body_handle key → sum).
    pub velocity_sum: std::collections::HashMap<(u32, u32), f64>,
    /// Number of steps each body has been below threshold.
    pub below_threshold_steps: std::collections::HashMap<(u32, u32), usize>,
    /// Linear velocity threshold.
    pub linear_threshold: f64,
    /// Number of steps required before a body is considered sleepable.
    pub required_steps: usize,
}
#[allow(dead_code)]
impl VelocityAccumulator {
    /// Create a new accumulator.
    pub fn new(linear_threshold: f64, required_steps: usize) -> Self {
        Self {
            velocity_sum: std::collections::HashMap::new(),
            below_threshold_steps: std::collections::HashMap::new(),
            linear_threshold,
            required_steps,
        }
    }
    /// Update the accumulator with the current velocity magnitude of a body.
    pub fn update(&mut self, index: u32, generation: u32, vel_magnitude: f64) {
        let key = (index, generation);
        *self.velocity_sum.entry(key).or_insert(0.0) += vel_magnitude;
        if vel_magnitude < self.linear_threshold {
            *self.below_threshold_steps.entry(key).or_insert(0) += 1;
        } else {
            self.below_threshold_steps.insert(key, 0);
        }
    }
    /// Check whether a body is eligible for sleeping.
    pub fn is_sleep_eligible(&self, index: u32, generation: u32) -> bool {
        let key = (index, generation);
        self.below_threshold_steps
            .get(&key)
            .map(|&steps| steps >= self.required_steps)
            .unwrap_or(false)
    }
    /// Check whether all bodies in an island are sleep-eligible.
    pub fn island_sleep_eligible(&self, island: &Island) -> bool {
        island
            .body_handles
            .iter()
            .all(|h| self.is_sleep_eligible(h.index, h.generation))
    }
    /// Reset the accumulator for a body (e.g. when it is woken up).
    pub fn reset(&mut self, index: u32, generation: u32) {
        let key = (index, generation);
        self.velocity_sum.remove(&key);
        self.below_threshold_steps.remove(&key);
    }
    /// Total number of tracked bodies.
    pub fn tracked_count(&self) -> usize {
        self.velocity_sum.len()
    }
}
/// An island of connected bodies and their associated constraints.
#[derive(Debug, Clone)]
pub struct Island {
    /// Body handles in this island.
    pub body_handles: Vec<BodyHandle>,
    /// Indices into the contact constraint array.
    pub contact_indices: Vec<usize>,
    /// Indices into the joint constraint array.
    pub joint_indices: Vec<usize>,
    /// Whether the entire island is sleeping.
    pub sleeping: bool,
}
impl Island {
    /// Create a new empty island.
    pub fn new() -> Self {
        Self {
            body_handles: Vec::new(),
            contact_indices: Vec::new(),
            joint_indices: Vec::new(),
            sleeping: false,
        }
    }
    /// Number of bodies in this island.
    #[allow(dead_code)]
    pub fn body_count(&self) -> usize {
        self.body_handles.len()
    }
    /// Total number of constraints (contacts + joints) in this island.
    #[allow(dead_code)]
    pub fn constraint_count(&self) -> usize {
        self.contact_indices.len() + self.joint_indices.len()
    }
    /// Whether the island has any constraints.
    #[allow(dead_code)]
    pub fn has_constraints(&self) -> bool {
        !self.contact_indices.is_empty() || !self.joint_indices.is_empty()
    }
    /// Check whether a body handle belongs to this island.
    #[allow(dead_code)]
    pub fn contains_body(&self, handle: BodyHandle) -> bool {
        self.body_handles
            .iter()
            .any(|&h| h.index == handle.index && h.generation == handle.generation)
    }
}
/// Per-island statistics for profiling.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct IslandSolveStats {
    /// Number of bodies.
    pub body_count: usize,
    /// Number of contact constraints.
    pub contact_count: usize,
    /// Number of joint constraints.
    pub joint_count: usize,
    /// Maximum linear speed in the island.
    pub max_linear_speed: f64,
    /// Maximum angular speed in the island.
    pub max_angular_speed: f64,
    /// Whether the island is sleeping.
    pub sleeping: bool,
}
/// Result of an island split analysis.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct IslandSplitResult {
    /// Index of the island that was analysed.
    pub island_index: usize,
    /// Whether splitting is recommended.
    pub should_split: bool,
    /// Suggested number of sub-islands.
    pub suggested_parts: usize,
    /// Reason for the recommendation.
    pub reason: String,
}
/// Aggregate statistics for the island contact graph.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct IslandContactGraphStats {
    /// Total number of islands.
    pub island_count: usize,
    /// Total awake islands.
    pub awake_island_count: usize,
    /// Total sleeping islands.
    pub sleeping_island_count: usize,
    /// Maximum contacts in any single island.
    pub max_contacts: usize,
    /// Maximum bodies in any single island.
    pub max_bodies: usize,
    /// Mean contacts per awake island.
    pub mean_contacts: f64,
    /// Mean bodies per awake island.
    pub mean_bodies: f64,
    /// Load imbalance ratio: max_cost / mean_cost (ideal = 1.0).
    pub load_imbalance: f64,
}
/// Per-island size summary.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct IslandSizeSummary {
    /// Index of the island.
    pub index: usize,
    /// Number of bodies.
    pub bodies: usize,
    /// Number of contacts.
    pub contacts: usize,
    /// Number of joints.
    pub joints: usize,
    /// Sleeping flag.
    pub sleeping: bool,
}
/// Iterator over the body handles in an island with their velocity data.
///
/// Yields `(BodyHandle, linear_speed, angular_speed)` for every dynamic body.
#[allow(dead_code)]
pub struct IslandBodyIter<'a> {
    pub(super) handles: std::slice::Iter<'a, BodyHandle>,
    pub(super) bodies: &'a RigidBodySet,
}
#[allow(dead_code)]
impl<'a> IslandBodyIter<'a> {
    /// Create a new iterator for the given island.
    pub fn new(island: &'a Island, bodies: &'a RigidBodySet) -> Self {
        Self {
            handles: island.body_handles.iter(),
            bodies,
        }
    }
}
/// A contact that lies on the boundary between two islands.
///
/// Boundary contacts arise when a newly detected contact connects bodies
/// from two previously separate islands.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct BoundaryContact {
    /// Index of the contact constraint.
    pub contact_index: usize,
    /// Island index containing body A.
    pub island_a: usize,
    /// Island index containing body B.
    pub island_b: usize,
}
/// Adjacency list representation of the body contact graph.
///
/// Each entry maps a body handle key `(index, generation)` to a list of
/// body handle keys it is in contact with.
#[allow(dead_code)]
pub struct BodyGraph {
    pub(super) adjacency: HashMap<(u32, u32), Vec<(u32, u32)>>,
}
#[allow(dead_code)]
impl BodyGraph {
    /// Build the body graph from contact constraints.
    pub fn from_contacts(contacts: &[ContactConstraint]) -> Self {
        let mut adjacency: HashMap<(u32, u32), Vec<(u32, u32)>> = HashMap::new();
        for contact in contacts {
            let ka = (
                contact.body_handle_a.index,
                contact.body_handle_a.generation,
            );
            let kb = (
                contact.body_handle_b.index,
                contact.body_handle_b.generation,
            );
            adjacency.entry(ka).or_default().push(kb);
            adjacency.entry(kb).or_default().push(ka);
        }
        Self { adjacency }
    }
    /// Build the body graph from contacts and joints.
    pub fn from_constraints(
        contacts: &[ContactConstraint],
        joints: &[Box<dyn Constraint>],
    ) -> Self {
        let mut adjacency: HashMap<(u32, u32), Vec<(u32, u32)>> = HashMap::new();
        for contact in contacts {
            let ka = (
                contact.body_handle_a.index,
                contact.body_handle_a.generation,
            );
            let kb = (
                contact.body_handle_b.index,
                contact.body_handle_b.generation,
            );
            adjacency.entry(ka).or_default().push(kb);
            adjacency.entry(kb).or_default().push(ka);
        }
        for joint in joints {
            let handles = joint.body_handles();
            if handles.len() >= 2 {
                let ka = (handles[0].index, handles[0].generation);
                let kb = (handles[1].index, handles[1].generation);
                adjacency.entry(ka).or_default().push(kb);
                adjacency.entry(kb).or_default().push(ka);
            }
        }
        Self { adjacency }
    }
    /// Get neighbors of a body.
    pub fn neighbors(&self, index: u32, generation: u32) -> &[(u32, u32)] {
        self.adjacency
            .get(&(index, generation))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }
    /// Number of bodies in the graph.
    pub fn node_count(&self) -> usize {
        self.adjacency.len()
    }
    /// Number of edges (each contact counted once per direction, so total / 2 = unique contacts).
    pub fn edge_count(&self) -> usize {
        let total: usize = self.adjacency.values().map(|v| v.len()).sum();
        total / 2
    }
    /// BFS traversal from a starting body, returning all reachable body keys.
    pub fn bfs_from(&self, start: (u32, u32)) -> Vec<(u32, u32)> {
        let mut visited = std::collections::HashSet::new();
        let mut queue = std::collections::VecDeque::new();
        let mut result = Vec::new();
        if !self.adjacency.contains_key(&start) {
            return result;
        }
        visited.insert(start);
        queue.push_back(start);
        while let Some(current) = queue.pop_front() {
            result.push(current);
            if let Some(neighbors) = self.adjacency.get(&current) {
                for &nb in neighbors {
                    if visited.insert(nb) {
                        queue.push_back(nb);
                    }
                }
            }
        }
        result
    }
    /// Find all connected components in the graph.
    ///
    /// Returns a list of components, where each component is a list of body keys.
    pub fn connected_components(&self) -> Vec<Vec<(u32, u32)>> {
        let mut visited = std::collections::HashSet::new();
        let mut components = Vec::new();
        for &key in self.adjacency.keys() {
            if visited.contains(&key) {
                continue;
            }
            let component = self.bfs_from(key);
            for &k in &component {
                visited.insert(k);
            }
            components.push(component);
        }
        components
    }
    /// Check if two bodies are connected (reachable from each other).
    pub fn are_connected(&self, a: (u32, u32), b: (u32, u32)) -> bool {
        let reachable = self.bfs_from(a);
        reachable.contains(&b)
    }
    /// Degree of a body (number of contacts/joints it participates in).
    pub fn degree(&self, index: u32, generation: u32) -> usize {
        self.adjacency
            .get(&(index, generation))
            .map(|v| v.len())
            .unwrap_or(0)
    }
}
/// Manages grouping of bodies into islands (connected components).
///
/// Bodies are connected if they share a contact or joint constraint.
/// Each island can be solved independently, enabling parallelism.
/// Islands where all bodies are below sleep thresholds are put to sleep.
#[derive(Debug)]
pub struct IslandManager {
    /// Linear velocity threshold for sleeping.
    pub linear_sleep_threshold: Real,
    /// Angular velocity threshold for sleeping.
    pub angular_sleep_threshold: Real,
    /// Time bodies must be below thresholds before sleeping.
    pub time_before_sleep: Real,
}
impl IslandManager {
    /// Create a new island manager with the given sleep parameters.
    pub fn new(
        linear_sleep_threshold: Real,
        angular_sleep_threshold: Real,
        time_before_sleep: Real,
    ) -> Self {
        Self {
            linear_sleep_threshold,
            angular_sleep_threshold,
            time_before_sleep,
        }
    }
    /// Build islands from bodies and their constraints.
    ///
    /// Uses union-find to group connected bodies. Static bodies are not
    /// included as island roots but serve as connectors.
    pub fn build_islands(
        &self,
        bodies: &RigidBodySet,
        contacts: &[ContactConstraint],
        joints: &[Box<dyn Constraint>],
    ) -> Vec<Island> {
        let all_handles: Vec<BodyHandle> = bodies.iter().map(|(h, _)| h).collect();
        if all_handles.is_empty() {
            return Vec::new();
        }
        let mut handle_to_idx: HashMap<(u32, u32), usize> = HashMap::new();
        for (i, h) in all_handles.iter().enumerate() {
            handle_to_idx.insert((h.index, h.generation), i);
        }
        let n = all_handles.len();
        let mut parent: Vec<usize> = (0..n).collect();
        let mut rank: Vec<usize> = vec![0; n];
        let find = |parent: &mut Vec<usize>, mut x: usize| -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        };
        let union = |parent: &mut Vec<usize>, rank: &mut Vec<usize>, a: usize, b: usize| {
            let ra = find(parent, a);
            let rb = find(parent, b);
            if ra == rb {
                return;
            }
            if rank[ra] < rank[rb] {
                parent[ra] = rb;
            } else if rank[ra] > rank[rb] {
                parent[rb] = ra;
            } else {
                parent[rb] = ra;
                rank[ra] += 1;
            }
        };
        for contact in contacts {
            let key_a = (
                contact.body_handle_a.index,
                contact.body_handle_a.generation,
            );
            let key_b = (
                contact.body_handle_b.index,
                contact.body_handle_b.generation,
            );
            if let (Some(&ia), Some(&ib)) = (handle_to_idx.get(&key_a), handle_to_idx.get(&key_b)) {
                union(&mut parent, &mut rank, ia, ib);
            }
        }
        for joint in joints {
            let handles = joint.body_handles();
            if handles.len() >= 2 {
                let key_a = (handles[0].index, handles[0].generation);
                let key_b = (handles[1].index, handles[1].generation);
                if let (Some(&ia), Some(&ib)) =
                    (handle_to_idx.get(&key_a), handle_to_idx.get(&key_b))
                {
                    union(&mut parent, &mut rank, ia, ib);
                }
            }
        }
        let mut root_to_island: HashMap<usize, usize> = HashMap::new();
        let mut islands: Vec<Island> = Vec::new();
        for (i, h) in all_handles.iter().enumerate() {
            if let Some(body) = bodies.get(*h)
                && body.body_type == BodyType::Static
            {
                continue;
            }
            let root = find(&mut parent, i);
            let island_idx = if let Some(&idx) = root_to_island.get(&root) {
                idx
            } else {
                let idx = islands.len();
                islands.push(Island::new());
                root_to_island.insert(root, idx);
                idx
            };
            islands[island_idx].body_handles.push(*h);
        }
        for (ci, contact) in contacts.iter().enumerate() {
            let key_a = (
                contact.body_handle_a.index,
                contact.body_handle_a.generation,
            );
            if let Some(&ia) = handle_to_idx.get(&key_a) {
                let root = find(&mut parent, ia);
                if let Some(&island_idx) = root_to_island.get(&root) {
                    islands[island_idx].contact_indices.push(ci);
                }
            }
        }
        for (ji, joint) in joints.iter().enumerate() {
            let handles = joint.body_handles();
            if !handles.is_empty() {
                let key = (handles[0].index, handles[0].generation);
                if let Some(&ia) = handle_to_idx.get(&key) {
                    let root = find(&mut parent, ia);
                    if let Some(&island_idx) = root_to_island.get(&root) {
                        islands[island_idx].joint_indices.push(ji);
                    }
                }
            }
        }
        for island in &mut islands {
            island.sleeping = self.check_island_sleeping(bodies, island);
        }
        islands
    }
    /// Check whether all dynamic bodies in an island qualify for sleeping.
    fn check_island_sleeping(&self, bodies: &RigidBodySet, island: &Island) -> bool {
        for &handle in &island.body_handles {
            if let Some(body) = bodies.get(handle) {
                if body.body_type != BodyType::Dynamic {
                    continue;
                }
                if body.state != BodyState::Sleeping
                    && (body.velocity.norm_squared()
                        >= self.linear_sleep_threshold * self.linear_sleep_threshold
                        || body.angular_velocity.norm_squared()
                            >= self.angular_sleep_threshold * self.angular_sleep_threshold)
                {
                    return false;
                }
            }
        }
        true
    }
    /// Put all bodies in a sleeping island to sleep.
    pub fn apply_sleeping(&self, bodies: &mut RigidBodySet, islands: &[Island]) {
        for island in islands {
            if island.sleeping {
                for &handle in &island.body_handles {
                    if let Some(body) = bodies.get_mut(handle)
                        && body.body_type == BodyType::Dynamic
                    {
                        body.state = BodyState::Sleeping;
                        body.velocity = oxiphysics_core::math::Vec3::zeros();
                        body.angular_velocity = oxiphysics_core::math::Vec3::zeros();
                    }
                }
            }
        }
    }
    /// Wake all bodies in a specific island.
    #[allow(dead_code)]
    pub fn wake_island(&self, bodies: &mut RigidBodySet, island: &Island) {
        for &handle in &island.body_handles {
            if let Some(body) = bodies.get_mut(handle)
                && body.state == BodyState::Sleeping
            {
                body.state = BodyState::Active;
            }
        }
    }
    /// Find the island containing a given body handle from a list of islands.
    #[allow(dead_code)]
    pub fn find_island_for_body(&self, handle: BodyHandle, islands: &[Island]) -> Option<usize> {
        islands
            .iter()
            .position(|island| island.contains_body(handle))
    }
}
/// A directed dependency between two islands.
///
/// If island A "depends on" island B it means a constraint spans both islands
/// and they must be solved jointly or in order.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct IslandDependency {
    /// Index of the source island.
    pub from: usize,
    /// Index of the target island.
    pub to: usize,
    /// Number of constraints linking these islands.
    pub constraint_count: usize,
}
/// A priority-ordered list of island indices based on estimated cost.
///
/// Higher-cost islands are placed first so a scheduler can pick them up
/// and assign them to worker threads before smaller islands.
#[allow(dead_code)]
pub struct IslandPriorityQueue {
    /// Ordered pairs of (cost, island_index).
    pub(super) items: Vec<(usize, usize)>,
}
#[allow(dead_code)]
impl IslandPriorityQueue {
    /// Build a priority queue from a list of islands (highest cost first).
    pub fn build(islands: &[Island]) -> Self {
        let mut items: Vec<(usize, usize)> = islands
            .iter()
            .enumerate()
            .filter(|(_, i)| !i.sleeping)
            .map(|(idx, i)| (i.constraint_count() + i.body_count(), idx))
            .collect();
        items.sort_by_key(|b| std::cmp::Reverse(b.0));
        Self { items }
    }
    /// Pop the highest-cost island index, or `None` if empty.
    pub fn pop(&mut self) -> Option<usize> {
        if self.items.is_empty() {
            None
        } else {
            Some(self.items.remove(0).1)
        }
    }
    /// Peek at the highest-cost island index without removing it.
    pub fn peek(&self) -> Option<usize> {
        self.items.first().map(|&(_, idx)| idx)
    }
    /// Number of islands remaining.
    pub fn len(&self) -> usize {
        self.items.len()
    }
    /// Whether the queue is empty.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
}
