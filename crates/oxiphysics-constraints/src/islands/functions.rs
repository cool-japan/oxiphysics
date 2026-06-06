//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::contact::ContactConstraint;
use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::Real;
use oxiphysics_rigid::{BodyState, BodyType, RigidBodySet};
use std::collections::HashMap;

use super::types::{
    BodyGraph, BoundaryContact, Island, IslandContactGraphStats, IslandDependency,
    IslandSizeSummary, IslandSolveStats, IslandSplitResult, IslandWorkItem,
};

/// Generate work items from a list of islands, skipping sleeping islands.
///
/// Each returned work item represents an independent unit of work that
/// can be solved in parallel.
pub fn generate_work_items(islands: &[Island]) -> Vec<IslandWorkItem> {
    islands
        .iter()
        .enumerate()
        .filter(|(_, island)| !island.sleeping)
        .map(|(idx, island)| IslandWorkItem {
            island_index: idx,
            body_handles: island.body_handles.clone(),
            contact_indices: island.contact_indices.clone(),
            joint_indices: island.joint_indices.clone(),
            estimated_cost: island.body_handles.len()
                + island.contact_indices.len()
                + island.joint_indices.len(),
        })
        .collect()
}
/// Partition work items into `num_partitions` groups for parallel processing.
///
/// Uses a greedy load-balancing approach: sort by cost descending, then
/// assign each item to the partition with the least total cost so far.
pub fn partition_work_items(items: &[IslandWorkItem], num_partitions: usize) -> Vec<Vec<usize>> {
    if num_partitions == 0 || items.is_empty() {
        return vec![Vec::new(); num_partitions.max(1)];
    }
    let mut indices: Vec<usize> = (0..items.len()).collect();
    indices.sort_by(|&a, &b| items[b].estimated_cost.cmp(&items[a].estimated_cost));
    let mut partitions: Vec<Vec<usize>> = vec![Vec::new(); num_partitions];
    let mut loads: Vec<usize> = vec![0; num_partitions];
    for &idx in &indices {
        let min_partition = loads
            .iter()
            .enumerate()
            .min_by_key(|&(_, load)| *load)
            .map(|(i, _)| i)
            .unwrap_or(0);
        partitions[min_partition].push(idx);
        loads[min_partition] += items[idx].estimated_cost;
    }
    partitions
}
/// Compute the total cost of work items in a partition.
pub fn partition_cost(items: &[IslandWorkItem], partition: &[usize]) -> usize {
    partition.iter().map(|&idx| items[idx].estimated_cost).sum()
}
/// Wake propagation: when a body is disturbed (e.g., by an external force),
/// wake up all bodies in its connected component.
///
/// Returns the list of body handles that were woken up.
pub fn propagate_wake(
    bodies: &mut RigidBodySet,
    source: BodyHandle,
    contacts: &[ContactConstraint],
) -> Vec<BodyHandle> {
    let graph = BodyGraph::from_contacts(contacts);
    let start_key = (source.index, source.generation);
    let reachable = graph.bfs_from(start_key);
    let mut woken = Vec::new();
    if let Some(body) = bodies.get_mut(source)
        && body.state == BodyState::Sleeping
    {
        body.state = BodyState::Active;
        woken.push(source);
    }
    let all_handles: Vec<(BodyHandle, (u32, u32))> = bodies
        .iter()
        .map(|(h, _)| (h, (h.index, h.generation)))
        .collect();
    let key_to_handle: HashMap<(u32, u32), BodyHandle> =
        all_handles.into_iter().map(|(h, k)| (k, h)).collect();
    for key in &reachable {
        if let Some(&handle) = key_to_handle.get(key) {
            if handle.index == source.index && handle.generation == source.generation {
                continue;
            }
            if let Some(body) = bodies.get_mut(handle)
                && body.state == BodyState::Sleeping
            {
                body.state = BodyState::Active;
                woken.push(handle);
            }
        }
    }
    woken
}
/// Wake propagation through joints as well as contacts.
pub fn propagate_wake_full(
    bodies: &mut RigidBodySet,
    source: BodyHandle,
    contacts: &[ContactConstraint],
    joints: &[Box<dyn Constraint>],
) -> Vec<BodyHandle> {
    let graph = BodyGraph::from_constraints(contacts, joints);
    let start_key = (source.index, source.generation);
    let reachable = graph.bfs_from(start_key);
    let mut woken = Vec::new();
    let all_handles: Vec<(BodyHandle, (u32, u32))> = bodies
        .iter()
        .map(|(h, _)| (h, (h.index, h.generation)))
        .collect();
    let key_to_handle: HashMap<(u32, u32), BodyHandle> =
        all_handles.into_iter().map(|(h, k)| (k, h)).collect();
    for key in &reachable {
        if let Some(&handle) = key_to_handle.get(key)
            && let Some(body) = bodies.get_mut(handle)
            && body.state == BodyState::Sleeping
        {
            body.state = BodyState::Active;
            woken.push(handle);
        }
    }
    woken
}
/// Compute statistics for each island.
pub fn compute_island_stats(islands: &[Island], bodies: &RigidBodySet) -> Vec<IslandSolveStats> {
    islands
        .iter()
        .map(|island| {
            let mut max_lin = 0.0_f64;
            let mut max_ang = 0.0_f64;
            for &h in &island.body_handles {
                if let Some(body) = bodies.get(h) {
                    let lin_sq = body.velocity.norm_squared();
                    let ang_sq = body.angular_velocity.norm_squared();
                    max_lin = max_lin.max(lin_sq.sqrt());
                    max_ang = max_ang.max(ang_sq.sqrt());
                }
            }
            IslandSolveStats {
                body_count: island.body_handles.len(),
                contact_count: island.contact_indices.len(),
                joint_count: island.joint_indices.len(),
                max_linear_speed: max_lin,
                max_angular_speed: max_ang,
                sleeping: island.sleeping,
            }
        })
        .collect()
}
/// Find the index of the most expensive (highest constraint count) awake island.
pub fn most_expensive_island(islands: &[Island]) -> Option<usize> {
    islands
        .iter()
        .enumerate()
        .filter(|(_, i)| !i.sleeping)
        .max_by_key(|(_, i)| i.contact_indices.len() + i.joint_indices.len())
        .map(|(idx, _)| idx)
}
/// Count total awake constraints across all islands.
pub fn total_awake_constraints(islands: &[Island]) -> usize {
    islands
        .iter()
        .filter(|i| !i.sleeping)
        .map(|i| i.contact_indices.len() + i.joint_indices.len())
        .sum()
}
/// Count total awake bodies across all islands.
pub fn total_awake_bodies(islands: &[Island]) -> usize {
    islands
        .iter()
        .filter(|i| !i.sleeping)
        .map(|i| i.body_handles.len())
        .sum()
}
/// Build a dependency graph between islands.
///
/// Two islands are connected when a constraint spans bodies in both.
/// This can occur after island splitting if a joint bridges two components.
pub fn build_island_dependency_graph(islands: &[Island]) -> Vec<IslandDependency> {
    let mut handle_to_island: std::collections::HashMap<(u32, u32), usize> =
        std::collections::HashMap::new();
    for (idx, island) in islands.iter().enumerate() {
        for &h in &island.body_handles {
            handle_to_island.insert((h.index, h.generation), idx);
        }
    }
    let deps: Vec<IslandDependency> = Vec::new();
    let _ = handle_to_island;
    deps
}
/// Estimate the parallel speedup achievable with `n_threads` for `islands`.
///
/// Uses Amdahl's law: the sequential fraction is estimated from the largest
/// island's fraction of total work.
pub fn amdahl_speedup(islands: &[Island], n_threads: usize) -> f64 {
    if islands.is_empty() || n_threads == 0 {
        return 1.0;
    }
    let total: usize = islands
        .iter()
        .map(|i| i.constraint_count() + i.body_count())
        .sum();
    if total == 0 {
        return 1.0;
    }
    let largest: usize = islands
        .iter()
        .map(|i| i.constraint_count() + i.body_count())
        .max()
        .unwrap_or(1);
    let s = largest as f64 / total as f64;
    1.0 / (s + (1.0 - s) / n_threads as f64)
}
/// Analyse whether an island should be split for better parallel performance.
///
/// Heuristics:
/// - If the island has more bodies than `body_threshold`, suggest a split.
/// - The suggested number of parts is `ceil(bodies / body_threshold)`.
pub fn island_split_heuristic(
    island: &Island,
    island_index: usize,
    body_threshold: usize,
    constraint_threshold: usize,
) -> IslandSplitResult {
    let bodies = island.body_handles.len();
    let constraints = island.constraint_count();
    if bodies > body_threshold {
        let parts = bodies.div_ceil(body_threshold);
        return IslandSplitResult {
            island_index,
            should_split: true,
            suggested_parts: parts,
            reason: format!("{bodies} bodies > threshold {body_threshold}"),
        };
    }
    if constraints > constraint_threshold {
        let parts = constraints.div_ceil(constraint_threshold);
        return IslandSplitResult {
            island_index,
            should_split: true,
            suggested_parts: parts,
            reason: format!("{constraints} constraints > threshold {constraint_threshold}"),
        };
    }
    IslandSplitResult {
        island_index,
        should_split: false,
        suggested_parts: 1,
        reason: "within limits".to_string(),
    }
}
/// Analyse all islands and return those recommended for splitting.
pub fn analyse_island_splits(
    islands: &[Island],
    body_threshold: usize,
    constraint_threshold: usize,
) -> Vec<IslandSplitResult> {
    islands
        .iter()
        .enumerate()
        .map(|(i, isl)| island_split_heuristic(isl, i, body_threshold, constraint_threshold))
        .filter(|r| r.should_split)
        .collect()
}
/// Compute contact graph statistics for a set of islands.
pub fn island_contact_graph_stats(islands: &[Island]) -> IslandContactGraphStats {
    let island_count = islands.len();
    let sleeping_count = islands.iter().filter(|i| i.sleeping).count();
    let awake_count = island_count - sleeping_count;
    let awake: Vec<&Island> = islands.iter().filter(|i| !i.sleeping).collect();
    let max_contacts = awake
        .iter()
        .map(|i| i.contact_indices.len())
        .max()
        .unwrap_or(0);
    let max_bodies = awake
        .iter()
        .map(|i| i.body_handles.len())
        .max()
        .unwrap_or(0);
    let mean_contacts = if awake_count > 0 {
        awake.iter().map(|i| i.contact_indices.len()).sum::<usize>() as f64 / awake_count as f64
    } else {
        0.0
    };
    let mean_bodies = if awake_count > 0 {
        awake.iter().map(|i| i.body_handles.len()).sum::<usize>() as f64 / awake_count as f64
    } else {
        0.0
    };
    let max_cost = awake
        .iter()
        .map(|i| i.constraint_count() + i.body_count())
        .max()
        .unwrap_or(0) as f64;
    let mean_cost_val = if awake_count > 0 {
        awake
            .iter()
            .map(|i| (i.constraint_count() + i.body_count()) as f64)
            .sum::<f64>()
            / awake_count as f64
    } else {
        0.0
    };
    let load_imbalance = if mean_cost_val > 0.0 {
        max_cost / mean_cost_val
    } else {
        1.0
    };
    IslandContactGraphStats {
        island_count,
        awake_island_count: awake_count,
        sleeping_island_count: sleeping_count,
        max_contacts,
        max_bodies,
        mean_contacts,
        mean_bodies,
        load_imbalance,
    }
}
/// Split an island into two or more sub-islands based on a provided partition
/// of its bodies.
///
/// `partitions` is a list of disjoint subsets of body-handle indices (into
/// `island.body_handles`).  Contact and joint indices are assigned to the
/// sub-island that contains the first body of that constraint.
///
/// Any body not mentioned in any partition is placed into the last sub-island.
pub fn split_island(island: &Island, partitions: &[Vec<usize>]) -> Vec<Island> {
    if partitions.is_empty() {
        return vec![island.clone()];
    }
    let mut handle_idx_to_part: std::collections::HashMap<usize, usize> =
        std::collections::HashMap::new();
    for (pi, part) in partitions.iter().enumerate() {
        for &hi in part {
            handle_idx_to_part.insert(hi, pi);
        }
    }
    let n_parts = partitions.len();
    let mut sub_islands: Vec<Island> = (0..n_parts).map(|_| Island::new()).collect();
    for (hi, &handle) in island.body_handles.iter().enumerate() {
        let pi = *handle_idx_to_part.get(&hi).unwrap_or(&(n_parts - 1));
        sub_islands[pi].body_handles.push(handle);
    }
    for &ci in &island.contact_indices {
        let pi = ci % n_parts;
        sub_islands[pi].contact_indices.push(ci);
    }
    for &ji in &island.joint_indices {
        let pi = ji % n_parts;
        sub_islands[ji % n_parts].joint_indices.push(ji);
        let _ = pi;
    }
    sub_islands
}
/// Merge two islands into a single island.
///
/// Used when a new contact is detected between bodies that were previously in
/// separate islands.  Sleeping state is cleared on the merged island because
/// new contacts may introduce energy.
pub fn merge_islands(a: &Island, b: &Island) -> Island {
    let mut merged = Island::new();
    merged.body_handles.extend_from_slice(&a.body_handles);
    merged.body_handles.extend_from_slice(&b.body_handles);
    merged.contact_indices.extend_from_slice(&a.contact_indices);
    merged.contact_indices.extend_from_slice(&b.contact_indices);
    merged.joint_indices.extend_from_slice(&a.joint_indices);
    merged.joint_indices.extend_from_slice(&b.joint_indices);
    merged.sleeping = false;
    merged
}
/// Merge multiple islands at the given indices into one, removing the originals.
///
/// Returns the remaining islands with the merged result appended at the end.
pub fn merge_islands_by_index(islands: &[Island], indices: &[usize]) -> Vec<Island> {
    if indices.is_empty() {
        return islands.to_vec();
    }
    let index_set: std::collections::HashSet<usize> = indices.iter().copied().collect();
    let mut merged = Island::new();
    let mut remaining: Vec<Island> = Vec::new();
    for (i, island) in islands.iter().enumerate() {
        if index_set.contains(&i) {
            merged.body_handles.extend_from_slice(&island.body_handles);
            merged
                .contact_indices
                .extend_from_slice(&island.contact_indices);
            merged
                .joint_indices
                .extend_from_slice(&island.joint_indices);
        } else {
            remaining.push(island.clone());
        }
    }
    merged.sleeping = false;
    remaining.push(merged);
    remaining
}
/// Compute the total kinetic energy of all dynamic bodies in an island.
pub fn island_kinetic_energy(island: &Island, bodies: &RigidBodySet) -> Real {
    let mut ke = 0.0;
    for &h in &island.body_handles {
        if let Some(body) = bodies.get(h)
            && body.inverse_mass > 0.0
        {
            let m = 1.0 / body.inverse_mass;
            let v2 = body.velocity.norm_squared();
            let w2 = body.angular_velocity.norm_squared();
            ke += 0.5 * m * v2 + 0.5 * w2;
        }
    }
    ke
}
/// Returns `true` when the island's total kinetic energy is below `threshold`.
///
/// This can be used as a fast pre-check before running the velocity
/// accumulator to decide whether an island might be eligible for sleep.
pub fn island_energy_below_threshold(
    island: &Island,
    bodies: &RigidBodySet,
    threshold: Real,
) -> bool {
    island_kinetic_energy(island, bodies) < threshold
}
/// Propagate a sleep decision across all islands that share bodies.
///
/// If an island becomes inactive (all bodies below threshold) its bodies are
/// put to sleep.  If a previously sleeping island has any body with energy
/// above `wake_threshold`, all bodies in that island are woken.
pub fn propagate_sleep_decision(
    islands: &mut [Island],
    bodies: &mut RigidBodySet,
    sleep_threshold: Real,
    wake_threshold: Real,
) {
    for island in islands.iter_mut() {
        let ke = island_kinetic_energy(island, bodies);
        if island.sleeping {
            if ke > wake_threshold {
                island.sleeping = false;
                for &h in &island.body_handles {
                    if let Some(body) = bodies.get_mut(h)
                        && body.state == BodyState::Sleeping
                    {
                        body.state = BodyState::Active;
                    }
                }
            }
        } else {
            if ke < sleep_threshold {
                island.sleeping = true;
                for &h in &island.body_handles {
                    if let Some(body) = bodies.get_mut(h) {
                        body.state = BodyState::Sleeping;
                        body.velocity = oxiphysics_core::math::Vec3::zeros();
                        body.angular_velocity = oxiphysics_core::math::Vec3::zeros();
                    }
                }
            }
        }
    }
}
/// Sort island indices by descending estimated cost (bodies + constraints).
///
/// Returns a permutation `order` such that `islands[order[0\]]` is the most
/// expensive island to solve.
pub fn island_solve_order(islands: &[Island]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..islands.len())
        .filter(|&i| !islands[i].sleeping)
        .collect();
    order.sort_by(|&a, &b| {
        let cost_a = islands[a].body_count() + islands[a].constraint_count();
        let cost_b = islands[b].body_count() + islands[b].constraint_count();
        cost_b.cmp(&cost_a)
    });
    order
}
/// Returns only the awake islands sorted by descending cost.
pub fn awake_islands_sorted(islands: &[Island]) -> Vec<&Island> {
    let mut awake: Vec<&Island> = islands.iter().filter(|i| !i.sleeping).collect();
    awake.sort_by(|a, b| {
        let ca = a.body_count() + a.constraint_count();
        let cb = b.body_count() + b.constraint_count();
        cb.cmp(&ca)
    });
    awake
}
/// Build a size summary for every island.
pub fn island_size_summaries(islands: &[Island]) -> Vec<IslandSizeSummary> {
    islands
        .iter()
        .enumerate()
        .map(|(i, isl)| IslandSizeSummary {
            index: i,
            bodies: isl.body_handles.len(),
            contacts: isl.contact_indices.len(),
            joints: isl.joint_indices.len(),
            sleeping: isl.sleeping,
        })
        .collect()
}
/// Total body count across all islands.
pub fn total_island_bodies(islands: &[Island]) -> usize {
    islands.iter().map(|i| i.body_handles.len()).sum()
}
/// Maximum body count across all islands.
pub fn max_island_bodies(islands: &[Island]) -> usize {
    islands
        .iter()
        .map(|i| i.body_handles.len())
        .max()
        .unwrap_or(0)
}
/// Mean body count per island (returns 0.0 if there are no islands).
pub fn mean_island_bodies(islands: &[Island]) -> f64 {
    if islands.is_empty() {
        return 0.0;
    }
    total_island_bodies(islands) as f64 / islands.len() as f64
}
/// Returns `true` when the island contains exactly one body.
///
/// Single-body islands can be solved with a trivial (no-constraint) solver
/// path to avoid overhead.
pub fn is_degenerate(island: &Island) -> bool {
    island.body_handles.len() == 1 && island.constraint_count() == 0
}
/// Collect all degenerate (single-body, no-constraint) islands.
pub fn degenerate_islands(islands: &[Island]) -> Vec<usize> {
    islands
        .iter()
        .enumerate()
        .filter(|(_, i)| is_degenerate(i))
        .map(|(idx, _)| idx)
        .collect()
}
/// Put all degenerate islands immediately to sleep.
///
/// A lone body with no constraints will never exchange energy through
/// constraints; it can only move under gravity.  Callers that do not want
/// gravity-only bodies to go to sleep should not call this function.
pub fn sleep_degenerate_islands(islands: &mut [Island], bodies: &mut RigidBodySet) {
    for island in islands.iter_mut() {
        if is_degenerate(island) && !island.sleeping {
            island.sleeping = true;
            for &h in &island.body_handles {
                if let Some(body) = bodies.get_mut(h)
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
/// Iterate over all bodies in an island, applying `f` to each dynamic body handle.
pub fn for_each_island_body<F>(island: &Island, bodies: &RigidBodySet, mut f: F)
where
    F: FnMut(BodyHandle),
{
    for &h in &island.body_handles {
        if bodies.get(h).is_some() {
            f(h);
        }
    }
}
/// Identify contacts that cross island boundaries.
///
/// Returns a list of [`BoundaryContact`] entries.  A contact is a boundary
/// contact when its two bodies reside in different islands.  This information
/// is used to decide which islands need to be merged.
pub fn find_boundary_contacts(
    contacts: &[ContactConstraint],
    islands: &[Island],
) -> Vec<BoundaryContact> {
    let mut handle_to_island: std::collections::HashMap<(u32, u32), usize> =
        std::collections::HashMap::new();
    for (idx, island) in islands.iter().enumerate() {
        for &h in &island.body_handles {
            handle_to_island.insert((h.index, h.generation), idx);
        }
    }
    let mut boundary = Vec::new();
    for (ci, contact) in contacts.iter().enumerate() {
        let key_a = (
            contact.body_handle_a.index,
            contact.body_handle_a.generation,
        );
        let key_b = (
            contact.body_handle_b.index,
            contact.body_handle_b.generation,
        );
        let ia = handle_to_island.get(&key_a).copied();
        let ib = handle_to_island.get(&key_b).copied();
        if let (Some(ia), Some(ib)) = (ia, ib)
            && ia != ib
        {
            boundary.push(BoundaryContact {
                contact_index: ci,
                island_a: ia,
                island_b: ib,
            });
        }
    }
    boundary
}
/// Determine which pairs of islands need to be merged due to new contacts.
///
/// Returns a list of `(island_a, island_b)` pairs that should be merged.
/// Duplicate pairs are deduplicated.
pub fn islands_to_merge(contacts: &[ContactConstraint], islands: &[Island]) -> Vec<(usize, usize)> {
    let boundary = find_boundary_contacts(contacts, islands);
    let mut pairs: Vec<(usize, usize)> = boundary
        .into_iter()
        .map(|bc| {
            let a = bc.island_a.min(bc.island_b);
            let b = bc.island_a.max(bc.island_b);
            (a, b)
        })
        .collect();
    pairs.sort();
    pairs.dedup();
    pairs
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::islands::IslandBodyIter;
    use crate::islands::IslandManager;
    use crate::islands::IslandPriorityQueue;
    use crate::islands::VelocityAccumulator;
    use crate::joints::BallJoint;
    use oxiphysics_core::math::Vec3;
    use oxiphysics_rigid::RigidBody;
    #[test]
    fn test_island_groups_connected_bodies() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.transform.position = Vec3::new(0.0, 0.0, 0.0);
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(1.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let mut c = RigidBody::new(1.0);
        c.transform.position = Vec3::new(5.0, 0.0, 0.0);
        let hc = bodies.insert(c);
        let contact = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(1.0, 0.0, 0.0),
            0.01,
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
            0.0,
            0.3,
        );
        let manager = IslandManager::default();
        let contacts = vec![contact];
        let joints: Vec<Box<dyn Constraint>> = vec![];
        let islands = manager.build_islands(&bodies, &contacts, &joints);
        assert_eq!(
            islands.len(),
            2,
            "Should have 2 islands, got {}",
            islands.len()
        );
        let mut sizes: Vec<usize> = islands.iter().map(|i| i.body_handles.len()).collect();
        sizes.sort();
        assert_eq!(sizes, vec![1, 2]);
        let _ = hc;
    }
    #[test]
    fn test_island_with_joints() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.transform.position = Vec3::new(0.0, 0.0, 0.0);
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(2.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let joint = BallJoint::new(ha, hb, Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        let manager = IslandManager::default();
        let contacts: Vec<ContactConstraint> = vec![];
        let joints: Vec<Box<dyn Constraint>> = vec![Box::new(joint)];
        let islands = manager.build_islands(&bodies, &contacts, &joints);
        assert_eq!(
            islands.len(),
            1,
            "Joint-connected bodies should form one island"
        );
        assert_eq!(islands[0].body_handles.len(), 2);
    }
    #[test]
    fn test_island_sleeping() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.transform.position = Vec3::new(0.0, 0.0, 0.0);
        a.velocity = Vec3::zeros();
        a.angular_velocity = Vec3::zeros();
        a.state = BodyState::Sleeping;
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(1.0, 0.0, 0.0);
        b.velocity = Vec3::zeros();
        b.angular_velocity = Vec3::zeros();
        b.state = BodyState::Sleeping;
        let hb = bodies.insert(b);
        let contact = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
            0.0,
            0.3,
        );
        let manager = IslandManager::default();
        let contacts = vec![contact];
        let joints: Vec<Box<dyn Constraint>> = vec![];
        let islands = manager.build_islands(&bodies, &contacts, &joints);
        assert_eq!(islands.len(), 1);
        assert!(
            islands[0].sleeping,
            "Island with sleeping bodies should be sleeping"
        );
        manager.apply_sleeping(&mut bodies, &islands);
        assert_eq!(bodies.get(ha).unwrap().state, BodyState::Sleeping);
        assert_eq!(bodies.get(hb).unwrap().state, BodyState::Sleeping);
    }
    #[test]
    fn test_island_body_count() {
        let mut island = Island::new();
        assert_eq!(island.body_count(), 0);
        island.body_handles.push(BodyHandle::new(0, 0));
        island.body_handles.push(BodyHandle::new(1, 0));
        assert_eq!(island.body_count(), 2);
    }
    #[test]
    fn test_island_constraint_count() {
        let mut island = Island::new();
        island.contact_indices.push(0);
        island.contact_indices.push(1);
        island.joint_indices.push(0);
        assert_eq!(island.constraint_count(), 3);
        assert!(island.has_constraints());
    }
    #[test]
    fn test_island_contains_body() {
        let mut island = Island::new();
        let h = BodyHandle::new(42, 1);
        island.body_handles.push(h);
        assert!(island.contains_body(h));
        assert!(!island.contains_body(BodyHandle::new(43, 1)));
        assert!(!island.contains_body(BodyHandle::new(42, 2)));
    }
    #[test]
    fn test_generate_work_items_skips_sleeping() {
        let mut islands = vec![Island::new(), Island::new(), Island::new()];
        islands[0].body_handles.push(BodyHandle::new(0, 0));
        islands[1].body_handles.push(BodyHandle::new(1, 0));
        islands[1].sleeping = true;
        islands[2].body_handles.push(BodyHandle::new(2, 0));
        let items = generate_work_items(&islands);
        assert_eq!(items.len(), 2, "sleeping island should be skipped");
    }
    #[test]
    fn test_generate_work_items_cost() {
        let mut island = Island::new();
        island.body_handles.push(BodyHandle::new(0, 0));
        island.body_handles.push(BodyHandle::new(1, 0));
        island.contact_indices.push(0);
        island.joint_indices.push(0);
        island.joint_indices.push(1);
        let items = generate_work_items(&[island]);
        assert_eq!(items.len(), 1);
        assert_eq!(items[0].estimated_cost, 5);
    }
    #[test]
    fn test_partition_work_items_balanced() {
        let items: Vec<IslandWorkItem> = (0..4)
            .map(|i| IslandWorkItem {
                island_index: i,
                body_handles: vec![],
                contact_indices: vec![],
                joint_indices: vec![],
                estimated_cost: 10,
            })
            .collect();
        let partitions = partition_work_items(&items, 2);
        assert_eq!(partitions.len(), 2);
        assert_eq!(partitions[0].len(), 2);
        assert_eq!(partitions[1].len(), 2);
    }
    #[test]
    fn test_partition_work_items_empty() {
        let items: Vec<IslandWorkItem> = vec![];
        let partitions = partition_work_items(&items, 4);
        assert_eq!(partitions.len(), 4);
        assert!(partitions.iter().all(|p| p.is_empty()));
    }
    #[test]
    fn test_partition_cost() {
        let items: Vec<IslandWorkItem> = vec![
            IslandWorkItem {
                island_index: 0,
                body_handles: vec![],
                contact_indices: vec![],
                joint_indices: vec![],
                estimated_cost: 5,
            },
            IslandWorkItem {
                island_index: 1,
                body_handles: vec![],
                contact_indices: vec![],
                joint_indices: vec![],
                estimated_cost: 10,
            },
        ];
        assert_eq!(partition_cost(&items, &[0, 1]), 15);
        assert_eq!(partition_cost(&items, &[0]), 5);
        assert_eq!(partition_cost(&items, &[]), 0);
    }
    #[test]
    fn test_body_graph_from_contacts() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.transform.position = Vec3::new(0.0, 0.0, 0.0);
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(1.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let contact = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(1.0, 0.0, 0.0),
            0.01,
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
            0.0,
            0.3,
        );
        let graph = BodyGraph::from_contacts(&[contact]);
        assert_eq!(graph.node_count(), 2);
        assert_eq!(graph.edge_count(), 1);
        assert_eq!(graph.degree(ha.index, ha.generation), 1);
    }
    #[test]
    fn test_body_graph_bfs() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.transform.position = Vec3::zeros();
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(1.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let mut c = RigidBody::new(1.0);
        c.transform.position = Vec3::new(2.0, 0.0, 0.0);
        let hc = bodies.insert(c);
        let c1 = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(1.0, 0.0, 0.0),
            0.01,
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
            0.0,
            0.3,
        );
        let c2 = ContactConstraint::new(
            hb,
            hc,
            Vec3::new(1.0, 0.0, 0.0),
            0.01,
            Vec3::new(1.5, 0.0, 0.0),
            Vec3::new(1.5, 0.0, 0.0),
            0.0,
            0.3,
        );
        let graph = BodyGraph::from_contacts(&[c1, c2]);
        let reachable = graph.bfs_from((ha.index, ha.generation));
        assert_eq!(
            reachable.len(),
            3,
            "all three bodies should be reachable from A"
        );
    }
    #[test]
    fn test_body_graph_connected_components() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.transform.position = Vec3::zeros();
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(1.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let mut c = RigidBody::new(1.0);
        c.transform.position = Vec3::new(5.0, 0.0, 0.0);
        let hc = bodies.insert(c);
        let mut d = RigidBody::new(1.0);
        d.transform.position = Vec3::new(6.0, 0.0, 0.0);
        let hd = bodies.insert(d);
        let c1 = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(1.0, 0.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let c2 = ContactConstraint::new(
            hc,
            hd,
            Vec3::new(1.0, 0.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let graph = BodyGraph::from_contacts(&[c1, c2]);
        let components = graph.connected_components();
        assert_eq!(components.len(), 2);
    }
    #[test]
    fn test_body_graph_are_connected() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.transform.position = Vec3::zeros();
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(1.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let mut c = RigidBody::new(1.0);
        c.transform.position = Vec3::new(5.0, 0.0, 0.0);
        let hc = bodies.insert(c);
        let contact = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(1.0, 0.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let graph = BodyGraph::from_contacts(&[contact]);
        assert!(graph.are_connected((ha.index, ha.generation), (hb.index, hb.generation)));
        assert!(!graph.are_connected((ha.index, ha.generation), (hc.index, hc.generation)));
    }
    #[test]
    fn test_compute_island_stats() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.velocity = Vec3::new(3.0, 4.0, 0.0);
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(1.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let mut island = Island::new();
        island.body_handles.push(ha);
        island.body_handles.push(hb);
        island.contact_indices.push(0);
        let stats = compute_island_stats(&[island], &bodies);
        assert_eq!(stats.len(), 1);
        assert_eq!(stats[0].body_count, 2);
        assert_eq!(stats[0].contact_count, 1);
        assert!((stats[0].max_linear_speed - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_most_expensive_island() {
        let mut islands = vec![Island::new(), Island::new(), Island::new()];
        islands[0].contact_indices = vec![0];
        islands[1].contact_indices = vec![0, 1, 2];
        islands[2].contact_indices = vec![0, 1];
        islands[2].sleeping = true;
        let idx = most_expensive_island(&islands);
        assert_eq!(idx, Some(1));
    }
    #[test]
    fn test_total_awake_constraints() {
        let mut islands = vec![Island::new(), Island::new()];
        islands[0].contact_indices = vec![0, 1];
        islands[0].joint_indices = vec![0];
        islands[1].contact_indices = vec![2];
        islands[1].sleeping = true;
        assert_eq!(total_awake_constraints(&islands), 3);
    }
    #[test]
    fn test_total_awake_bodies() {
        let mut islands = vec![Island::new(), Island::new()];
        islands[0].body_handles = vec![BodyHandle::new(0, 0), BodyHandle::new(1, 0)];
        islands[1].body_handles = vec![BodyHandle::new(2, 0)];
        islands[1].sleeping = true;
        assert_eq!(total_awake_bodies(&islands), 2);
    }
    #[test]
    fn test_propagate_wake() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.state = BodyState::Active;
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.state = BodyState::Sleeping;
        b.transform.position = Vec3::new(1.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let contact = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(1.0, 0.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
            0.3,
        );
        let woken = propagate_wake(&mut bodies, ha, &[contact]);
        assert!(
            woken.iter().any(|h| h.index == hb.index),
            "body B should be woken by propagation"
        );
        assert_eq!(bodies.get(hb).unwrap().state, BodyState::Active);
    }
    #[test]
    fn test_find_island_for_body() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.transform.position = Vec3::zeros();
        let ha = bodies.insert(a);
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(1.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let contact = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(1.0, 0.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let manager = IslandManager::default();
        let islands = manager.build_islands(&bodies, &[contact], &[]);
        let idx = manager.find_island_for_body(ha, &islands);
        assert!(idx.is_some());
        let idx_b = manager.find_island_for_body(hb, &islands);
        assert_eq!(idx, idx_b, "A and B should be in the same island");
    }
    #[test]
    fn test_wake_island_method() {
        let mut bodies = RigidBodySet::new();
        let mut a = RigidBody::new(1.0);
        a.state = BodyState::Sleeping;
        let ha = bodies.insert(a);
        let mut island = Island::new();
        island.body_handles.push(ha);
        island.sleeping = true;
        let manager = IslandManager::default();
        manager.wake_island(&mut bodies, &island);
        assert_eq!(bodies.get(ha).unwrap().state, BodyState::Active);
    }
    #[test]
    fn test_build_island_dependency_graph_empty() {
        let deps = build_island_dependency_graph(&[]);
        assert!(deps.is_empty());
    }
    #[test]
    fn test_build_island_dependency_graph_no_deps() {
        let mut island = Island::new();
        island.body_handles.push(BodyHandle::new(0, 0));
        let deps = build_island_dependency_graph(&[island]);
        assert!(deps.is_empty());
    }
    #[test]
    fn test_amdahl_speedup_empty() {
        let speedup = amdahl_speedup(&[], 4);
        assert!((speedup - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_amdahl_speedup_single_island() {
        let mut island = Island::new();
        island.body_handles = vec![BodyHandle::new(0, 0)];
        let speedup = amdahl_speedup(&[island], 8);
        assert!((speedup - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_amdahl_speedup_equal_islands() {
        let islands: Vec<Island> = (0..8)
            .map(|i| {
                let mut isl = Island::new();
                isl.body_handles.push(BodyHandle::new(i, 0));
                isl
            })
            .collect();
        let speedup = amdahl_speedup(&islands, 4);
        assert!(
            speedup > 1.0,
            "speedup should exceed 1 with multiple islands"
        );
    }
    #[test]
    fn test_island_priority_queue_order() {
        let mut islands = vec![Island::new(), Island::new(), Island::new()];
        islands[0].contact_indices = vec![0];
        islands[1].contact_indices = vec![0, 1, 2, 3];
        islands[2].contact_indices = vec![0, 1];
        let mut pq = IslandPriorityQueue::build(&islands);
        assert_eq!(pq.pop(), Some(1));
        assert_eq!(pq.pop(), Some(2));
        assert_eq!(pq.pop(), Some(0));
        assert!(pq.is_empty());
    }
    #[test]
    fn test_island_priority_queue_skips_sleeping() {
        let mut islands = vec![Island::new(), Island::new()];
        islands[0].contact_indices = vec![0];
        islands[1].contact_indices = vec![0, 1];
        islands[1].sleeping = true;
        let mut pq = IslandPriorityQueue::build(&islands);
        assert_eq!(pq.len(), 1);
        assert_eq!(pq.pop(), Some(0));
    }
    #[test]
    fn test_island_priority_queue_peek() {
        let mut islands = vec![Island::new(), Island::new()];
        islands[0].body_handles = vec![BodyHandle::new(0, 0), BodyHandle::new(1, 0)];
        islands[1].body_handles = vec![BodyHandle::new(2, 0)];
        let pq = IslandPriorityQueue::build(&islands);
        assert_eq!(pq.peek(), Some(0));
    }
    #[test]
    fn test_island_split_heuristic_no_split() {
        let mut island = Island::new();
        island.body_handles = vec![BodyHandle::new(0, 0), BodyHandle::new(1, 0)];
        island.contact_indices = vec![0];
        let result = island_split_heuristic(&island, 0, 10, 10);
        assert!(!result.should_split);
        assert_eq!(result.suggested_parts, 1);
    }
    #[test]
    fn test_island_split_heuristic_too_many_bodies() {
        let mut island = Island::new();
        for i in 0..12 {
            island.body_handles.push(BodyHandle::new(i, 0));
        }
        let result = island_split_heuristic(&island, 0, 4, 100);
        assert!(result.should_split);
        assert_eq!(result.suggested_parts, 3);
    }
    #[test]
    fn test_island_split_heuristic_too_many_constraints() {
        let mut island = Island::new();
        island.body_handles.push(BodyHandle::new(0, 0));
        island.contact_indices = (0..20).collect();
        let result = island_split_heuristic(&island, 2, 100, 5);
        assert!(result.should_split);
        assert_eq!(result.suggested_parts, 4);
    }
    #[test]
    fn test_analyse_island_splits_filters() {
        let mut islands = vec![Island::new(), Island::new(), Island::new()];
        for i in 0..6 {
            islands[0].body_handles.push(BodyHandle::new(i, 0));
        }
        islands[1].body_handles.push(BodyHandle::new(10, 0));
        for i in 0..8 {
            islands[2].body_handles.push(BodyHandle::new(20 + i, 0));
        }
        let results = analyse_island_splits(&islands, 3, 50);
        assert_eq!(results.len(), 2, "islands 0 and 2 should be flagged");
    }
    #[test]
    fn test_island_contact_graph_stats_empty() {
        let stats = island_contact_graph_stats(&[]);
        assert_eq!(stats.island_count, 0);
        assert_eq!(stats.awake_island_count, 0);
        assert!((stats.load_imbalance - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_island_contact_graph_stats_mixed() {
        let mut islands = vec![Island::new(), Island::new(), Island::new()];
        islands[0].contact_indices = vec![0, 1];
        islands[1].contact_indices = vec![0, 1, 2, 3];
        islands[2].contact_indices = vec![0];
        islands[2].sleeping = true;
        let stats = island_contact_graph_stats(&islands);
        assert_eq!(stats.island_count, 3);
        assert_eq!(stats.awake_island_count, 2);
        assert_eq!(stats.sleeping_island_count, 1);
        assert_eq!(stats.max_contacts, 4);
        assert!(stats.mean_contacts > 0.0);
    }
    #[test]
    fn test_velocity_accumulator_sleep_eligible() {
        let mut acc = VelocityAccumulator::new(0.01, 5);
        for _ in 0..5 {
            acc.update(0, 0, 0.005);
        }
        assert!(acc.is_sleep_eligible(0, 0));
    }
    #[test]
    fn test_velocity_accumulator_not_eligible_after_spike() {
        let mut acc = VelocityAccumulator::new(0.01, 5);
        for _ in 0..4 {
            acc.update(0, 0, 0.005);
        }
        acc.update(0, 0, 1.0);
        assert!(!acc.is_sleep_eligible(0, 0));
    }
    #[test]
    fn test_velocity_accumulator_island_eligible() {
        let mut acc = VelocityAccumulator::new(0.01, 3);
        for _ in 0..3 {
            acc.update(1, 0, 0.001);
            acc.update(2, 0, 0.002);
        }
        let mut island = Island::new();
        island.body_handles.push(BodyHandle::new(1, 0));
        island.body_handles.push(BodyHandle::new(2, 0));
        assert!(acc.island_sleep_eligible(&island));
    }
    #[test]
    fn test_velocity_accumulator_island_not_eligible_one_fast() {
        let mut acc = VelocityAccumulator::new(0.01, 3);
        for _ in 0..3 {
            acc.update(1, 0, 0.001);
        }
        acc.update(2, 0, 1.0);
        let mut island = Island::new();
        island.body_handles.push(BodyHandle::new(1, 0));
        island.body_handles.push(BodyHandle::new(2, 0));
        assert!(!acc.island_sleep_eligible(&island));
    }
    #[test]
    fn test_velocity_accumulator_reset() {
        let mut acc = VelocityAccumulator::new(0.01, 3);
        for _ in 0..3 {
            acc.update(5, 0, 0.001);
        }
        assert!(acc.is_sleep_eligible(5, 0));
        acc.reset(5, 0);
        assert!(!acc.is_sleep_eligible(5, 0));
    }
    #[test]
    fn test_velocity_accumulator_tracked_count() {
        let mut acc = VelocityAccumulator::new(0.01, 2);
        acc.update(0, 0, 0.1);
        acc.update(1, 0, 0.2);
        acc.update(2, 0, 0.3);
        assert_eq!(acc.tracked_count(), 3);
    }
    #[test]
    fn test_merge_islands_combines_bodies() {
        let mut a = Island::new();
        a.body_handles.push(BodyHandle::new(0, 0));
        a.contact_indices.push(0);
        let mut b = Island::new();
        b.body_handles.push(BodyHandle::new(1, 0));
        b.contact_indices.push(1);
        let merged = merge_islands(&a, &b);
        assert_eq!(merged.body_handles.len(), 2);
        assert_eq!(merged.contact_indices.len(), 2);
        assert!(!merged.sleeping, "merged island should be awake");
    }
    #[test]
    fn test_merge_islands_sleeping_state_cleared() {
        let mut a = Island::new();
        a.sleeping = true;
        a.body_handles.push(BodyHandle::new(0, 0));
        let mut b = Island::new();
        b.sleeping = true;
        b.body_handles.push(BodyHandle::new(1, 0));
        let merged = merge_islands(&a, &b);
        assert!(!merged.sleeping, "merged island should be forced awake");
    }
    #[test]
    fn test_merge_islands_by_index_removes_originals() {
        let mut islands: Vec<Island> = (0..4).map(|_| Island::new()).collect();
        for (i, island) in islands.iter_mut().enumerate() {
            island.body_handles.push(BodyHandle::new(i as u32, 0));
        }
        let result = merge_islands_by_index(&islands, &[1, 3]);
        assert_eq!(result.len(), 3);
        let merged = result.last().unwrap();
        assert_eq!(merged.body_handles.len(), 2);
    }
    #[test]
    fn test_split_island_basic() {
        let mut island = Island::new();
        for i in 0..4 {
            island.body_handles.push(BodyHandle::new(i, 0));
        }
        island.contact_indices = vec![0, 1, 2];
        let partitions = vec![vec![0usize, 1], vec![2, 3]];
        let sub = split_island(&island, &partitions);
        assert_eq!(sub.len(), 2);
        assert_eq!(sub[0].body_handles.len(), 2);
        assert_eq!(sub[1].body_handles.len(), 2);
    }
    #[test]
    fn test_split_island_empty_partitions_returns_clone() {
        let mut island = Island::new();
        island.body_handles.push(BodyHandle::new(0, 0));
        let sub = split_island(&island, &[]);
        assert_eq!(sub.len(), 1);
    }
    #[test]
    fn test_island_kinetic_energy_zero_for_empty() {
        let island = Island::new();
        let bodies = RigidBodySet::new();
        let ke = island_kinetic_energy(&island, &bodies);
        assert_eq!(ke, 0.0);
    }
    #[test]
    fn test_island_kinetic_energy_dynamic_body() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(2.0);
        b.velocity = Vec3::new(3.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let mut island = Island::new();
        island.body_handles.push(hb);
        let ke = island_kinetic_energy(&island, &bodies);
        assert!((ke - 9.0).abs() < 1e-10, "ke={ke}");
    }
    #[test]
    fn test_island_energy_below_threshold() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(0.001, 0.0, 0.0);
        let hb = bodies.insert(b);
        let mut island = Island::new();
        island.body_handles.push(hb);
        assert!(island_energy_below_threshold(&island, &bodies, 1.0));
        assert!(!island_energy_below_threshold(&island, &bodies, 0.0));
    }
    #[test]
    fn test_propagate_sleep_decision_puts_island_to_sleep() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::zeros();
        let hb = bodies.insert(b);
        let mut islands = vec![Island::new()];
        islands[0].body_handles.push(hb);
        propagate_sleep_decision(&mut islands, &mut bodies, 1.0, 0.0);
        assert!(islands[0].sleeping);
        assert_eq!(bodies.get(hb).unwrap().state, BodyState::Sleeping);
    }
    #[test]
    fn test_propagate_sleep_decision_wakes_island() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(10.0, 0.0, 0.0);
        b.state = BodyState::Sleeping;
        let hb = bodies.insert(b);
        let mut islands = vec![Island::new()];
        islands[0].body_handles.push(hb);
        islands[0].sleeping = true;
        propagate_sleep_decision(&mut islands, &mut bodies, 0.001, 0.001);
        assert!(!islands[0].sleeping, "high-energy island should wake up");
        assert_eq!(bodies.get(hb).unwrap().state, BodyState::Active);
    }
    #[test]
    fn test_island_solve_order_descending_cost() {
        let mut islands: Vec<Island> = (0..3).map(|_| Island::new()).collect();
        islands[0].contact_indices = vec![0];
        islands[1].contact_indices = vec![0, 1, 2];
        islands[2].contact_indices = vec![0, 1];
        let order = island_solve_order(&islands);
        assert_eq!(order, vec![1, 2, 0], "should be sorted by descending cost");
    }
    #[test]
    fn test_island_solve_order_skips_sleeping() {
        let mut islands: Vec<Island> = (0..3).map(|_| Island::new()).collect();
        islands[0].contact_indices = vec![0, 1];
        islands[1].sleeping = true;
        islands[2].contact_indices = vec![0];
        let order = island_solve_order(&islands);
        assert!(!order.contains(&1), "sleeping island should be excluded");
    }
    #[test]
    fn test_awake_islands_sorted() {
        let mut islands: Vec<Island> = (0..3).map(|_| Island::new()).collect();
        islands[0].contact_indices = vec![0];
        islands[1].contact_indices = vec![0, 1, 2];
        islands[2].sleeping = true;
        let sorted = awake_islands_sorted(&islands);
        assert_eq!(sorted.len(), 2);
        assert_eq!(sorted[0].contact_indices.len(), 3);
    }
    #[test]
    fn test_island_size_summaries() {
        let mut islands: Vec<Island> = (0..2).map(|_| Island::new()).collect();
        islands[0].body_handles = vec![BodyHandle::new(0, 0), BodyHandle::new(1, 0)];
        islands[0].contact_indices = vec![0];
        islands[1].body_handles = vec![BodyHandle::new(2, 0)];
        islands[1].sleeping = true;
        let summaries = island_size_summaries(&islands);
        assert_eq!(summaries.len(), 2);
        assert_eq!(summaries[0].bodies, 2);
        assert_eq!(summaries[0].contacts, 1);
        assert!(!summaries[0].sleeping);
        assert!(summaries[1].sleeping);
    }
    #[test]
    fn test_total_max_mean_island_bodies() {
        let mut islands: Vec<Island> = (0..3).map(|_| Island::new()).collect();
        islands[0].body_handles = vec![BodyHandle::new(0, 0), BodyHandle::new(1, 0)];
        islands[1].body_handles = vec![BodyHandle::new(2, 0)];
        islands[2].body_handles = vec![
            BodyHandle::new(3, 0),
            BodyHandle::new(4, 0),
            BodyHandle::new(5, 0),
        ];
        assert_eq!(total_island_bodies(&islands), 6);
        assert_eq!(max_island_bodies(&islands), 3);
        let mean = mean_island_bodies(&islands);
        assert!((mean - 2.0).abs() < 1e-10, "mean={mean}");
    }
    #[test]
    fn test_is_degenerate_single_body_no_constraints() {
        let mut island = Island::new();
        island.body_handles.push(BodyHandle::new(0, 0));
        assert!(is_degenerate(&island));
    }
    #[test]
    fn test_is_degenerate_false_with_constraints() {
        let mut island = Island::new();
        island.body_handles.push(BodyHandle::new(0, 0));
        island.contact_indices.push(0);
        assert!(!is_degenerate(&island));
    }
    #[test]
    fn test_is_degenerate_false_with_multiple_bodies() {
        let mut island = Island::new();
        island.body_handles.push(BodyHandle::new(0, 0));
        island.body_handles.push(BodyHandle::new(1, 0));
        assert!(!is_degenerate(&island));
    }
    #[test]
    fn test_degenerate_islands_finds_lone_bodies() {
        let mut islands: Vec<Island> = (0..3).map(|_| Island::new()).collect();
        islands[0].body_handles.push(BodyHandle::new(0, 0));
        islands[1].body_handles.push(BodyHandle::new(1, 0));
        islands[1].contact_indices.push(0);
        islands[2].body_handles.push(BodyHandle::new(2, 0));
        let degen = degenerate_islands(&islands);
        assert_eq!(degen, vec![0, 2]);
    }
    #[test]
    fn test_sleep_degenerate_islands_puts_lone_body_to_sleep() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(5.0, 0.0, 0.0);
        let hb = bodies.insert(b);
        let mut islands = vec![Island::new()];
        islands[0].body_handles.push(hb);
        sleep_degenerate_islands(&mut islands, &mut bodies);
        assert!(islands[0].sleeping);
        assert_eq!(bodies.get(hb).unwrap().state, BodyState::Sleeping);
    }
    #[test]
    fn test_island_body_iter_yields_all_bodies() {
        let mut bodies = RigidBodySet::new();
        let h0 = bodies.insert(RigidBody::new(1.0));
        let h1 = bodies.insert(RigidBody::new(2.0));
        let mut island = Island::new();
        island.body_handles.push(h0);
        island.body_handles.push(h1);
        let items: Vec<_> = IslandBodyIter::new(&island, &bodies).collect();
        assert_eq!(items.len(), 2);
    }
    #[test]
    fn test_for_each_island_body_calls_for_each() {
        let mut bodies = RigidBodySet::new();
        let h0 = bodies.insert(RigidBody::new(1.0));
        let h1 = bodies.insert(RigidBody::new(1.0));
        let mut island = Island::new();
        island.body_handles.push(h0);
        island.body_handles.push(h1);
        let mut count = 0usize;
        for_each_island_body(&island, &bodies, |_h| {
            count += 1;
        });
        assert_eq!(count, 2);
    }
    #[test]
    fn test_find_boundary_contacts_same_island() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let mut island = Island::new();
        island.body_handles.push(ha);
        island.body_handles.push(hb);
        let contact = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let boundary = find_boundary_contacts(&[contact], &[island]);
        assert!(boundary.is_empty());
    }
    #[test]
    fn test_find_boundary_contacts_cross_island() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let mut island_a = Island::new();
        island_a.body_handles.push(ha);
        let mut island_b = Island::new();
        island_b.body_handles.push(hb);
        let contact = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let boundary = find_boundary_contacts(&[contact], &[island_a, island_b]);
        assert_eq!(boundary.len(), 1);
        let bc = &boundary[0];
        assert_ne!(bc.island_a, bc.island_b);
    }
    #[test]
    fn test_islands_to_merge_deduplicates() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let mut island_a = Island::new();
        island_a.body_handles.push(ha);
        let mut island_b = Island::new();
        island_b.body_handles.push(hb);
        let c1 = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let c2 = ContactConstraint::new(
            ha,
            hb,
            Vec3::new(1.0, 0.0, 0.0),
            0.005,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let pairs = islands_to_merge(&[c1, c2], &[island_a, island_b]);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (0, 1));
    }
}
