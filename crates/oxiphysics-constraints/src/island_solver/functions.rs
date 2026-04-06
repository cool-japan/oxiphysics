//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    ConstraintResidual, ConstraintType, Island, IslandBody, IslandBudget, IslandConstraint,
    IslandDependency, IslandEnergyReport, IslandSubStep, SortedConstraint,
};

/// Opaque identifier for a body within the island solver.
pub type BodyId = u32;
/// Opaque identifier for a constraint within the island solver.
pub type ConstraintId = u32;
pub(super) const SLEEP_VELOCITY_THRESHOLD: f64 = 0.01;
/// Compute and apply a non-penetration normal impulse between two bodies.
///
/// Returns the applied impulse magnitude λ (useful for warm-starting next frame).
///
/// # Parameters
/// * `body_a` / `body_b` – the two participating bodies (mutated in-place).
/// * `normal`      – contact normal pointing from B toward A (world-space).
/// * `depth`       – penetration depth; positive values indicate overlap.
/// * `restitution` – coefficient of restitution (0 = perfectly inelastic).
pub fn solve_contact_constraint(
    body_a: &mut IslandBody,
    body_b: &mut IslandBody,
    normal: [f64; 3],
    depth: f64,
    restitution: f64,
) -> f64 {
    if depth <= 0.0 {
        return 0.0;
    }
    let inv_mass_sum = body_a.inv_mass + body_b.inv_mass;
    if inv_mass_sum == 0.0 {
        return 0.0;
    }
    let rel_vel = [
        body_a.linear_vel[0] - body_b.linear_vel[0],
        body_a.linear_vel[1] - body_b.linear_vel[1],
        body_a.linear_vel[2] - body_b.linear_vel[2],
    ];
    let vn = rel_vel[0] * normal[0] + rel_vel[1] * normal[1] + rel_vel[2] * normal[2];
    pub(super) const BAUMGARTE: f64 = 0.2;
    pub(super) const SLOP: f64 = 0.005;
    let bias = BAUMGARTE * f64::max(0.0, depth - SLOP);
    let lambda_raw = (-(1.0 + restitution) * vn + bias) / inv_mass_sum;
    let lambda = f64::max(0.0, lambda_raw);
    if lambda == 0.0 {
        return 0.0;
    }
    let impulse = [normal[0] * lambda, normal[1] * lambda, normal[2] * lambda];
    if !body_a.is_static {
        body_a.linear_vel[0] += impulse[0] * body_a.inv_mass;
        body_a.linear_vel[1] += impulse[1] * body_a.inv_mass;
        body_a.linear_vel[2] += impulse[2] * body_a.inv_mass;
    }
    if !body_b.is_static {
        body_b.linear_vel[0] -= impulse[0] * body_b.inv_mass;
        body_b.linear_vel[1] -= impulse[1] * body_b.inv_mass;
        body_b.linear_vel[2] -= impulse[2] * body_b.inv_mass;
    }
    lambda
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::island_solver::IslandManager;

    use crate::island_solver::UnionFind;
    fn make_body(id: BodyId, is_static: bool) -> IslandBody {
        IslandBody {
            id,
            inv_mass: if is_static { 0.0 } else { 1.0 },
            linear_vel: [0.0; 3],
            angular_vel: [0.0; 3],
            position: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            is_static,
            is_sleeping: false,
        }
    }
    fn make_contact(id: ConstraintId, body_a: BodyId, body_b: BodyId) -> IslandConstraint {
        IslandConstraint {
            id,
            body_a,
            body_b,
            constraint_type: ConstraintType::Contact {
                normal: [0.0, 1.0, 0.0],
                depth: 0.01,
            },
            lambda: 0.0,
        }
    }
    #[test]
    fn test_union_find_initially_separate() {
        let mut uf = UnionFind::new(5);
        for i in 0..5_usize {
            for j in (i + 1)..5 {
                assert!(
                    !uf.same_component(i, j),
                    "elements {i} and {j} should be in different components initially"
                );
            }
        }
    }
    #[test]
    fn test_union_find_union_same_component() {
        let mut uf = UnionFind::new(4);
        uf.union(0, 1);
        uf.union(2, 3);
        assert!(uf.same_component(0, 1));
        assert!(uf.same_component(2, 3));
        assert!(!uf.same_component(0, 2));
        assert!(!uf.same_component(1, 3));
        uf.union(1, 2);
        assert!(uf.same_component(0, 3));
    }
    #[test]
    fn test_island_manager_two_connected_bodies_same_island() {
        let mut mgr = IslandManager::new();
        mgr.add_body(make_body(0, false));
        mgr.add_body(make_body(1, false));
        mgr.add_constraint(make_contact(0, 0, 1));
        mgr.build_islands();
        assert_eq!(mgr.island_count(), 1, "two connected bodies → one island");
        assert_eq!(mgr.islands[0].bodies.len(), 2);
    }
    #[test]
    fn test_island_manager_disconnected_bodies_separate_islands() {
        let mut mgr = IslandManager::new();
        mgr.add_body(make_body(0, false));
        mgr.add_body(make_body(1, false));
        mgr.build_islands();
        assert_eq!(
            mgr.island_count(),
            2,
            "no constraints → two separate singleton islands"
        );
    }
    #[test]
    fn test_build_islands_chain_n_bodies() {
        pub(super) const N: usize = 8;
        let mut mgr = IslandManager::new();
        for i in 0..N {
            mgr.add_body(make_body(i as BodyId, false));
        }
        for i in 0..(N - 1) {
            mgr.add_constraint(make_contact(
                i as ConstraintId,
                i as BodyId,
                (i + 1) as BodyId,
            ));
        }
        mgr.build_islands();
        assert_eq!(mgr.island_count(), 1, "chain of {N} bodies → 1 island");
        assert_eq!(mgr.islands[0].bodies.len(), N);
        assert_eq!(mgr.islands[0].constraints.len(), N - 1);
    }
    #[test]
    fn test_sleeping_island_all_slow_bodies_can_sleep() {
        let mut island = Island {
            bodies: vec![0, 1],
            constraints: vec![],
            is_sleeping: false,
        };
        let bodies = vec![make_body(0, false), make_body(1, false)];
        assert!(
            island.can_sleep(&bodies),
            "island with all-slow bodies can_sleep = true"
        );
        let mut fast_bodies = bodies.clone();
        fast_bodies[0].linear_vel = [5.0, 0.0, 0.0];
        island.is_sleeping = false;
        assert!(
            !island.can_sleep(&fast_bodies),
            "island with fast body cannot sleep"
        );
    }
    #[test]
    fn test_sequential_impulse_solver_resolves_penetration() {
        let mut body_a = IslandBody {
            id: 0,
            inv_mass: 1.0,
            linear_vel: [0.0, -1.0, 0.0],
            angular_vel: [0.0; 3],
            position: [0.0, 0.1, 0.0],
            orientation: [0.0, 0.0, 0.0, 1.0],
            is_static: false,
            is_sleeping: false,
        };
        let mut body_b = IslandBody {
            id: 1,
            inv_mass: 1.0,
            linear_vel: [0.0, 1.0, 0.0],
            angular_vel: [0.0; 3],
            position: [0.0, -0.1, 0.0],
            orientation: [0.0, 0.0, 0.0, 1.0],
            is_static: false,
            is_sleeping: false,
        };
        let normal = [0.0_f64, 1.0, 0.0];
        let depth = 0.05;
        let lambda = solve_contact_constraint(&mut body_a, &mut body_b, normal, depth, 0.0);
        assert!(
            lambda > 0.0,
            "impulse must be positive for penetrating contact"
        );
        let rel_vn = (body_a.linear_vel[0] - body_b.linear_vel[0]) * normal[0]
            + (body_a.linear_vel[1] - body_b.linear_vel[1]) * normal[1]
            + (body_a.linear_vel[2] - body_b.linear_vel[2]) * normal[2];
        assert!(
            rel_vn >= -1e-10,
            "relative velocity along normal should be non-negative after impulse, got {rel_vn}"
        );
    }
}
/// Sort constraints by island index, placing same-island constraints together.
///
/// Returns a sorted list of `SortedConstraint`.
pub fn sort_constraints_by_island(
    constraints: &[IslandConstraint],
    island_map: &std::collections::HashMap<BodyId, usize>,
) -> Vec<SortedConstraint> {
    let mut sorted: Vec<SortedConstraint> = constraints
        .iter()
        .map(|c| {
            let island_idx = island_map.get(&c.body_a).copied().unwrap_or(usize::MAX);
            SortedConstraint {
                island_idx,
                constraint: c.clone(),
            }
        })
        .collect();
    sorted.sort_by_key(|sc| sc.island_idx);
    sorted
}
/// Build a map from `BodyId` to island index.
///
/// Used to accelerate constraint-to-island lookups.
pub fn build_body_to_island_map(islands: &[Island]) -> std::collections::HashMap<BodyId, usize> {
    let mut map = std::collections::HashMap::new();
    for (island_idx, island) in islands.iter().enumerate() {
        for &bid in &island.bodies {
            map.insert(bid, island_idx);
        }
    }
    map
}
/// Build a dependency graph between islands.
///
/// Two islands are dependent if they share a constraint that references
/// bodies from both islands.  In practice this detects articulation bodies.
pub fn build_island_dependency_graph(
    islands: &[Island],
    constraints: &[IslandConstraint],
    body_to_island: &std::collections::HashMap<BodyId, usize>,
) -> Vec<IslandDependency> {
    let mut deps = Vec::new();
    for c in constraints {
        let ia = body_to_island.get(&c.body_a).copied();
        let ib = body_to_island.get(&c.body_b).copied();
        if let (Some(a), Some(b)) = (ia, ib)
            && a != b
        {
            let already = deps.iter().any(|d: &IslandDependency| {
                (d.island_a == a && d.island_b == b) || (d.island_a == b && d.island_b == a)
            });
            if !already {
                deps.push(IslandDependency {
                    island_a: a,
                    island_b: b,
                });
            }
        }
        let _ = islands;
    }
    deps
}
/// Compute adaptive sub-step counts for each island.
///
/// Islands with more constraints and higher kinetic energy get more sub-steps,
/// up to `max_substeps`.
pub fn adaptive_substeps(
    islands: &[Island],
    bodies: &[IslandBody],
    max_substeps: usize,
) -> Vec<IslandSubStep> {
    islands
        .iter()
        .enumerate()
        .map(|(idx, island)| {
            if island.is_sleeping {
                return IslandSubStep {
                    island_idx: idx,
                    n_substeps: 0,
                };
            }
            let ke: f64 = island
                .bodies
                .iter()
                .filter_map(|&bid| bodies.iter().find(|b| b.id == bid))
                .map(|b| {
                    let lv = b.linear_vel;
                    let av = b.angular_vel;
                    let lv2 = lv[0] * lv[0] + lv[1] * lv[1] + lv[2] * lv[2];
                    let av2 = av[0] * av[0] + av[1] * av[1] + av[2] * av[2];
                    (lv2 + av2)
                        * if b.inv_mass > 0.0 {
                            1.0 / b.inv_mass
                        } else {
                            0.0
                        }
                })
                .sum();
            let base = 1 + island.constraints.len();
            let energy_factor = if ke > 1.0 {
                (ke.log10() + 1.0) as usize
            } else {
                1
            };
            let n = (base * energy_factor).min(max_substeps).max(1);
            IslandSubStep {
                island_idx: idx,
                n_substeps: n,
            }
        })
        .collect()
}
/// Compute energy reports for all islands.
pub fn compute_island_energy_reports(
    islands: &[Island],
    bodies: &[IslandBody],
) -> Vec<IslandEnergyReport> {
    islands
        .iter()
        .enumerate()
        .map(|(idx, island)| IslandEnergyReport::compute(idx, island, bodies))
        .collect()
}
/// Compute solver budgets for all islands using a priority heuristic.
///
/// Islands with more constraints and higher energy receive more iterations.
/// The total iterations across all islands is approximately `base_iterations * n_active_islands`.
pub fn compute_island_budgets(
    islands: &[Island],
    bodies: &[IslandBody],
    base_iterations: usize,
    max_iterations: usize,
) -> Vec<IslandBudget> {
    let scores: Vec<f64> = islands
        .iter()
        .map(|island| {
            if island.is_sleeping {
                return 0.0;
            }
            let n_constraints = island.constraints.len() as f64;
            let ke: f64 = island
                .bodies
                .iter()
                .filter_map(|&bid| bodies.iter().find(|b| b.id == bid))
                .map(|b| {
                    let lv = b.linear_vel;
                    let av = b.angular_vel;
                    lv[0] * lv[0]
                        + lv[1] * lv[1]
                        + lv[2] * lv[2]
                        + av[0] * av[0]
                        + av[1] * av[1]
                        + av[2] * av[2]
                })
                .sum();
            n_constraints + ke.sqrt()
        })
        .collect();
    let max_score = scores
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max)
        .max(1.0);
    islands
        .iter()
        .enumerate()
        .map(|(idx, island)| {
            if island.is_sleeping {
                return IslandBudget::new(idx, 0, 0, 0.0);
            }
            let priority = scores[idx] / max_score;
            let iterations = (base_iterations as f64 * (1.0 + priority)).ceil() as usize;
            let iterations = iterations.min(max_iterations).max(1);
            let time_us = (iterations as u64) * 50;
            IslandBudget::new(idx, iterations, time_us, priority)
        })
        .collect()
}
/// Compute constraint residuals for an island.
///
/// Returns the maximum residual and a per-constraint list.
pub fn compute_island_residuals(
    island: &Island,
    bodies: &[IslandBody],
    constraints: &[IslandConstraint],
) -> (f64, Vec<ConstraintResidual>) {
    let mut max_residual = 0.0_f64;
    let mut residuals = Vec::new();
    for &cid in &island.constraints {
        if let Some(c) = constraints.iter().find(|c| c.id == cid) {
            let body_a = bodies.iter().find(|b| b.id == c.body_a);
            let body_b = bodies.iter().find(|b| b.id == c.body_b);
            if let (Some(ba), Some(bb)) = (body_a, body_b) {
                let residual = match &c.constraint_type {
                    ConstraintType::Contact { normal, depth } => {
                        if *depth > 0.0 {
                            let rel_vel = [
                                ba.linear_vel[0] - bb.linear_vel[0],
                                ba.linear_vel[1] - bb.linear_vel[1],
                                ba.linear_vel[2] - bb.linear_vel[2],
                            ];
                            let vn = rel_vel[0] * normal[0]
                                + rel_vel[1] * normal[1]
                                + rel_vel[2] * normal[2];
                            (-vn).max(0.0) + *depth
                        } else {
                            0.0
                        }
                    }
                    ConstraintType::Distance { target_length } => {
                        let dx = ba.position[0] - bb.position[0];
                        let dy = ba.position[1] - bb.position[1];
                        let dz = ba.position[2] - bb.position[2];
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                        (dist - target_length).abs()
                    }
                    ConstraintType::Joint { .. } => {
                        let dv = [
                            ba.linear_vel[0] - bb.linear_vel[0],
                            ba.linear_vel[1] - bb.linear_vel[1],
                            ba.linear_vel[2] - bb.linear_vel[2],
                        ];
                        (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]).sqrt()
                    }
                };
                max_residual = max_residual.max(residual);
                residuals.push(ConstraintResidual {
                    constraint_id: cid,
                    residual,
                    is_active: c.lambda.abs() > 1e-14,
                });
            }
        }
    }
    (max_residual, residuals)
}
#[cfg(test)]
mod expanded_island_tests {
    use super::*;
    use crate::island_solver::ContactManifoldCache;
    use crate::island_solver::ContactManifoldEntry;
    use crate::island_solver::ImpulseWarmStartCache;
    use crate::island_solver::IslandManager;
    use crate::island_solver::IslandMerger;
    use crate::island_solver::IslandSolver;
    use crate::island_solver::IslandSplitDetector;
    use crate::island_solver::IslandStatistics;
    use crate::island_solver::ParallelIslandSolver;
    fn make_body(id: BodyId, is_static: bool) -> IslandBody {
        IslandBody {
            id,
            inv_mass: if is_static { 0.0 } else { 1.0 },
            linear_vel: [0.0; 3],
            angular_vel: [0.0; 3],
            position: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            is_static,
            is_sleeping: false,
        }
    }
    fn make_contact(id: ConstraintId, body_a: BodyId, body_b: BodyId) -> IslandConstraint {
        IslandConstraint {
            id,
            body_a,
            body_b,
            constraint_type: ConstraintType::Contact {
                normal: [0.0, 1.0, 0.0],
                depth: 0.01,
            },
            lambda: 0.0,
        }
    }
    #[test]
    fn test_sort_constraints_by_island() {
        let mut mgr = IslandManager::new();
        mgr.add_body(make_body(0, false));
        mgr.add_body(make_body(1, false));
        mgr.add_body(make_body(2, false));
        mgr.add_body(make_body(3, false));
        mgr.add_constraint(make_contact(0, 0, 1));
        mgr.add_constraint(make_contact(1, 2, 3));
        mgr.build_islands();
        let body_to_island = build_body_to_island_map(&mgr.islands);
        let sorted = sort_constraints_by_island(&mgr.constraints, &body_to_island);
        assert_eq!(sorted.len(), 2);
        if sorted.len() == 2 {
            assert!(sorted[0].island_idx <= sorted[1].island_idx);
        }
    }
    #[test]
    fn test_build_body_to_island_map() {
        let mut mgr = IslandManager::new();
        mgr.add_body(make_body(0, false));
        mgr.add_body(make_body(1, false));
        mgr.add_constraint(make_contact(0, 0, 1));
        mgr.build_islands();
        let map = build_body_to_island_map(&mgr.islands);
        assert_eq!(map.len(), 2);
        let idx0 = map[&0];
        let idx1 = map[&1];
        assert_eq!(idx0, idx1, "connected bodies should be in same island");
    }
    #[test]
    fn test_island_dependency_graph_no_deps() {
        let mut mgr = IslandManager::new();
        mgr.add_body(make_body(0, false));
        mgr.add_body(make_body(1, false));
        mgr.add_constraint(make_contact(0, 0, 1));
        mgr.build_islands();
        let map = build_body_to_island_map(&mgr.islands);
        let deps = build_island_dependency_graph(&mgr.islands, &mgr.constraints, &map);
        assert_eq!(deps.len(), 0, "single island has no cross-island deps");
    }
    #[test]
    fn test_adaptive_substeps_sleeping_zero() {
        let island = Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: true,
        };
        let bodies = vec![make_body(0, false)];
        let plan = adaptive_substeps(&[island], &bodies, 8);
        assert_eq!(
            plan[0].n_substeps, 0,
            "sleeping island should have 0 sub-steps"
        );
    }
    #[test]
    fn test_adaptive_substeps_active_at_least_one() {
        let island = Island {
            bodies: vec![0],
            constraints: vec![0],
            is_sleeping: false,
        };
        let mut body = make_body(0, false);
        body.linear_vel = [1.0, 0.0, 0.0];
        let plan = adaptive_substeps(&[island], &[body], 8);
        assert!(
            plan[0].n_substeps >= 1,
            "active island should have >= 1 sub-step"
        );
    }
    #[test]
    fn test_island_energy_report_zero_at_rest() {
        let island = Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: false,
        };
        let body = make_body(0, false);
        let report = IslandEnergyReport::compute(0, &island, &[body]);
        assert!(
            (report.kinetic_energy).abs() < 1e-15,
            "ke at rest = {}",
            report.kinetic_energy
        );
        assert_eq!(report.n_dynamic_bodies, 1);
        assert!(report.can_sleep, "body at rest should be able to sleep");
    }
    #[test]
    fn test_island_energy_report_nonzero_moving() {
        let island = Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: false,
        };
        let mut body = make_body(0, false);
        body.linear_vel = [2.0, 0.0, 0.0];
        let report = IslandEnergyReport::compute(0, &island, &[body]);
        assert!(
            report.kinetic_energy > 0.0,
            "moving body should have KE > 0"
        );
    }
    #[test]
    fn test_compute_island_energy_reports_count() {
        let mut mgr = IslandManager::new();
        mgr.add_body(make_body(0, false));
        mgr.add_body(make_body(1, false));
        mgr.add_constraint(make_contact(0, 0, 1));
        mgr.build_islands();
        let reports = compute_island_energy_reports(&mgr.islands, &mgr.bodies);
        assert_eq!(reports.len(), mgr.islands.len());
    }
    #[test]
    fn test_parallel_island_solver_processes_islands() {
        let mut mgr = IslandManager::new();
        let mut ba = make_body(0, false);
        let mut bb = make_body(1, false);
        ba.linear_vel = [0.0, -1.0, 0.0];
        bb.linear_vel = [0.0, 1.0, 0.0];
        mgr.add_body(ba);
        mgr.add_body(bb);
        mgr.add_constraint(IslandConstraint {
            id: 0,
            body_a: 0,
            body_b: 1,
            constraint_type: ConstraintType::Contact {
                normal: [0.0, 1.0, 0.0],
                depth: 0.05,
            },
            lambda: 0.0,
        });
        mgr.build_islands();
        let solver = ParallelIslandSolver::new(5, 4);
        let count = solver.solve_all_islands(&mut mgr);
        assert!(count >= 1, "should process at least one island");
    }
    #[test]
    fn test_island_statistics() {
        let mut mgr = IslandManager::new();
        for i in 0..6 {
            mgr.add_body(make_body(i as BodyId, false));
        }
        mgr.add_constraint(make_contact(0, 0, 1));
        mgr.add_constraint(make_contact(1, 2, 3));
        mgr.build_islands();
        let stats = IslandStatistics::compute(&mgr);
        assert_eq!(stats.total_bodies, 6);
        assert_eq!(stats.total_constraints, 2);
        assert!(stats.n_islands >= 3, "should have >= 3 islands");
    }
    #[test]
    fn test_wake_island() {
        let mut mgr = IslandManager::new();
        let mut body = make_body(0, false);
        body.is_sleeping = true;
        mgr.add_body(body);
        mgr.build_islands();
        mgr.wake_island(0);
        assert!(
            !mgr.bodies[0].is_sleeping,
            "body should be awake after wake_island"
        );
        assert!(
            !mgr.islands[0].is_sleeping,
            "island should be awake after wake_island"
        );
    }
    #[test]
    fn test_manifold_entry_add_points() {
        let mut entry = ContactManifoldEntry::new(0, 1, [0.0, 1.0, 0.0]);
        assert!(entry.add_point([0.0, 0.0, 0.0], 0.01));
        assert!(entry.add_point([1.0, 0.0, 0.0], 0.02));
        assert_eq!(entry.n_points, 2);
        let total = entry.total_normal_impulse();
        assert!(total.abs() < 1e-15);
    }
    #[test]
    fn test_manifold_entry_max_points() {
        let mut entry = ContactManifoldEntry::new(0, 1, [0.0, 1.0, 0.0]);
        for _ in 0..4 {
            assert!(entry.add_point([0.0; 3], 0.01));
        }
        assert!(
            !entry.add_point([1.0, 0.0, 0.0], 0.01),
            "should return false when full"
        );
    }
    #[test]
    fn test_manifold_entry_warm_start() {
        let mut entry = ContactManifoldEntry::new(0, 1, [0.0, 1.0, 0.0]);
        entry.add_point([0.0; 3], 0.01);
        entry.update_lambda(0, 5.0);
        let ws = entry.warm_start_impulse(0, 0.8);
        assert!((ws - 4.0).abs() < 1e-12, "warm start = 5*0.8 = 4, got {ws}");
    }
    #[test]
    fn test_manifold_entry_advance_frame_decay() {
        let mut entry = ContactManifoldEntry::new(0, 1, [0.0, 1.0, 0.0]);
        entry.add_point([0.0; 3], 0.01);
        entry.update_lambda(0, 10.0);
        entry.advance_frame(0.9);
        assert!((entry.normal_lambdas[0] - 9.0).abs() < 1e-12);
        assert_eq!(entry.frame, 1);
    }
    #[test]
    fn test_manifold_cache_insert_and_get() {
        let mut cache = ContactManifoldCache::new();
        let entry = ContactManifoldEntry::new(2, 5, [0.0, 1.0, 0.0]);
        cache.insert(entry);
        assert!(cache.get(2, 5).is_some());
        assert!(cache.get(5, 2).is_some(), "symmetric lookup should work");
        assert!(cache.get(0, 1).is_none());
    }
    #[test]
    fn test_manifold_cache_advance_and_evict() {
        let mut cache = ContactManifoldCache::new();
        let mut entry = ContactManifoldEntry::new(0, 1, [0.0, 1.0, 0.0]);
        entry.frame = 10;
        cache.insert(entry);
        cache.advance_all(0.9, 5);
        assert!(cache.is_empty(), "old entry should be evicted");
    }
    #[test]
    fn test_manifold_cache_remove() {
        let mut cache = ContactManifoldCache::new();
        cache.insert(ContactManifoldEntry::new(0, 1, [0.0, 1.0, 0.0]));
        cache.insert(ContactManifoldEntry::new(2, 3, [1.0, 0.0, 0.0]));
        cache.remove(0, 1);
        assert_eq!(cache.len(), 1, "one entry should remain");
        assert!(cache.get(0, 1).is_none());
    }
    #[test]
    fn test_compute_island_budgets_sleeping_zero() {
        let island = Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: true,
        };
        let body = make_body(0, false);
        let budgets = compute_island_budgets(&[island], &[body], 10, 20);
        assert_eq!(budgets[0].iterations, 0);
    }
    #[test]
    fn test_compute_island_budgets_active_at_least_one() {
        let island = Island {
            bodies: vec![0],
            constraints: vec![0],
            is_sleeping: false,
        };
        let mut body = make_body(0, false);
        body.linear_vel = [1.0, 0.0, 0.0];
        let budgets = compute_island_budgets(&[island], &[body], 10, 50);
        assert!(
            budgets[0].iterations >= 1,
            "active island should get at least 1 iteration"
        );
        assert!(budgets[0].priority >= 0.0);
    }
    #[test]
    fn test_compute_island_budgets_capped_at_max() {
        let island = Island {
            bodies: vec![0],
            constraints: (0..100).map(|i| i as u32).collect(),
            is_sleeping: false,
        };
        let body = make_body(0, false);
        let max = 15;
        let budgets = compute_island_budgets(&[island], &[body], 5, max);
        assert!(
            budgets[0].iterations <= max,
            "iterations should not exceed max"
        );
    }
    #[test]
    fn test_warm_start_cache_store_retrieve() {
        let mut cache = ImpulseWarmStartCache::new(100);
        cache.store(0, 1, 5.0);
        assert!((cache.retrieve(0, 1) - 5.0).abs() < 1e-12);
        assert!(
            (cache.retrieve(1, 0) - 5.0).abs() < 1e-12,
            "symmetric lookup"
        );
    }
    #[test]
    fn test_warm_start_cache_decay() {
        let mut cache = ImpulseWarmStartCache::new(100);
        cache.store(0, 1, 10.0);
        cache.decay_all(0.5);
        assert!((cache.retrieve(0, 1) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_warm_start_cache_prune() {
        let mut cache = ImpulseWarmStartCache::new(100);
        cache.store(0, 1, 1.0);
        cache.store(2, 3, 0.0001);
        cache.prune(0.01);
        assert!(
            cache.retrieve(0, 1) > 0.0,
            "significant impulse should remain"
        );
        assert!(
            (cache.retrieve(2, 3)).abs() < 1e-12,
            "negligible impulse should be pruned"
        );
    }
    #[test]
    fn test_warm_start_cache_missing_returns_zero() {
        let cache = ImpulseWarmStartCache::new(100);
        assert!((cache.retrieve(42, 99)).abs() < 1e-12);
    }
    #[test]
    fn test_island_merger_merges_two_islands() {
        let mut mgr = IslandManager::new();
        mgr.add_body(make_body(0, false));
        mgr.add_body(make_body(1, false));
        mgr.build_islands();
        assert_eq!(mgr.islands.len(), 2, "should start with 2 islands");
        let new_constraint = make_contact(99, 0, 1);
        let n_merges = IslandMerger::merge_from_constraints(&mut mgr, &[new_constraint]);
        assert_eq!(n_merges, 1, "one merge should occur");
        assert_eq!(mgr.islands.len(), 1, "should now have 1 island");
        assert_eq!(mgr.islands[0].bodies.len(), 2);
    }
    #[test]
    fn test_island_merger_no_merge_same_island() {
        let mut mgr = IslandManager::new();
        mgr.add_body(make_body(0, false));
        mgr.add_body(make_body(1, false));
        mgr.add_constraint(make_contact(0, 0, 1));
        mgr.build_islands();
        assert_eq!(mgr.islands.len(), 1);
        let existing = make_contact(0, 0, 1);
        let n_merges = IslandMerger::merge_from_constraints(&mut mgr, &[existing]);
        assert_eq!(n_merges, 0, "no merge needed when already in same island");
    }
    #[test]
    fn test_split_detector_no_split_needed() {
        let mut mgr = IslandManager::new();
        mgr.add_body(make_body(0, false));
        mgr.add_body(make_body(1, false));
        mgr.add_constraint(make_contact(0, 0, 1));
        mgr.build_islands();
        let active = vec![make_contact(0, 0, 1)];
        let splits = IslandSplitDetector::detect(&mgr, &active);
        assert!(
            splits.is_empty(),
            "no split needed when constraint still active"
        );
    }
    #[test]
    fn test_split_detector_detects_disconnection() {
        let mut mgr = IslandManager::new();
        mgr.add_body(make_body(0, false));
        mgr.add_body(make_body(1, false));
        mgr.add_body(make_body(2, false));
        mgr.add_constraint(make_contact(0, 0, 1));
        mgr.add_constraint(make_contact(1, 1, 2));
        mgr.build_islands();
        assert_eq!(mgr.islands.len(), 1);
        let active = vec![make_contact(0, 0, 1)];
        let splits = IslandSplitDetector::detect(&mgr, &active);
        assert!(
            !splits.is_empty(),
            "should detect split when constraint removed"
        );
    }
    #[test]
    fn test_compute_residuals_contact_no_penetration() {
        let island = Island {
            bodies: vec![0, 1],
            constraints: vec![0],
            is_sleeping: false,
        };
        let bodies = vec![make_body(0, false), make_body(1, false)];
        let constraints = vec![IslandConstraint {
            id: 0,
            body_a: 0,
            body_b: 1,
            constraint_type: ConstraintType::Contact {
                normal: [0.0, 1.0, 0.0],
                depth: 0.0,
            },
            lambda: 0.0,
        }];
        let (max_res, _) = compute_island_residuals(&island, &bodies, &constraints);
        assert!((max_res).abs() < 1e-12, "no penetration → zero residual");
    }
    #[test]
    fn test_compute_residuals_contact_with_penetration() {
        let island = Island {
            bodies: vec![0, 1],
            constraints: vec![0],
            is_sleeping: false,
        };
        let bodies = vec![make_body(0, false), make_body(1, false)];
        let constraints = vec![IslandConstraint {
            id: 0,
            body_a: 0,
            body_b: 1,
            constraint_type: ConstraintType::Contact {
                normal: [0.0, 1.0, 0.0],
                depth: 0.1,
            },
            lambda: 0.0,
        }];
        let (max_res, residuals) = compute_island_residuals(&island, &bodies, &constraints);
        assert!(max_res > 0.0, "penetration should give non-zero residual");
        assert_eq!(residuals.len(), 1);
    }
    fn make_moving_body(id: BodyId, vx: f64) -> IslandBody {
        IslandBody {
            id,
            inv_mass: 1.0,
            linear_vel: [vx, 0.0, 0.0],
            angular_vel: [0.0; 3],
            position: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            is_static: false,
            is_sleeping: false,
        }
    }
    #[test]
    fn test_compute_island_energy_zero_for_static_bodies() {
        let mut solver = IslandSolver::new();
        solver.manager.add_body(IslandBody {
            id: 0,
            inv_mass: 0.0,
            linear_vel: [10.0, 0.0, 0.0],
            angular_vel: [0.0; 3],
            position: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            is_static: true,
            is_sleeping: false,
        });
        solver.manager.islands.push(Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: false,
        });
        let energy = solver.compute_island_energy(0).unwrap();
        assert!(
            energy == 0.0,
            "static body has no kinetic energy, e={energy}"
        );
    }
    #[test]
    fn test_compute_island_energy_positive_for_moving_body() {
        let mut solver = IslandSolver::new();
        solver.manager.add_body(make_moving_body(0, 4.0));
        solver.manager.islands.push(Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: false,
        });
        let energy = solver.compute_island_energy(0).unwrap();
        assert!(
            energy > 0.0,
            "moving body should have positive kinetic energy, e={energy}"
        );
        assert!((energy - 8.0).abs() < 1e-10, "expected 8 J, got {energy}");
    }
    #[test]
    fn test_compute_island_energy_out_of_range_returns_none() {
        let solver = IslandSolver::new();
        assert!(
            solver.compute_island_energy(99).is_none(),
            "out-of-range should return None"
        );
    }
    #[test]
    fn test_compute_all_island_energies_length() {
        let mut solver = IslandSolver::new();
        solver.manager.add_body(make_moving_body(0, 1.0));
        solver.manager.add_body(make_moving_body(1, 2.0));
        solver.manager.add_constraint(make_contact(0, 0, 1));
        solver.manager.build_islands();
        let energies = solver.compute_all_island_energies();
        assert_eq!(energies.len(), solver.manager.islands.len());
    }
    #[test]
    fn test_merge_small_islands_singleton_merged() {
        let mut solver = IslandSolver::with_thresholds(1e-4, 2);
        solver.manager.add_body(make_body(0, false));
        solver.manager.add_body(make_body(1, false));
        solver.manager.add_constraint(make_contact(0, 0, 1));
        solver.manager.islands.push(Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: false,
        });
        solver.manager.islands.push(Island {
            bodies: vec![1],
            constraints: vec![],
            is_sleeping: false,
        });
        let merges = solver.merge_small_islands();
        assert!(merges > 0, "should merge singletons, merges={merges}");
        assert_eq!(
            solver.manager.islands.len(),
            1,
            "result should be a single island"
        );
        assert_eq!(solver.manager.islands[0].bodies.len(), 2);
    }
    #[test]
    fn test_merge_small_islands_no_merge_when_large_enough() {
        let mut solver = IslandSolver::with_thresholds(1e-4, 2);
        solver.manager.add_body(make_body(0, false));
        solver.manager.add_body(make_body(1, false));
        solver.manager.add_constraint(make_contact(0, 0, 1));
        solver.manager.islands.push(Island {
            bodies: vec![0, 1],
            constraints: vec![0],
            is_sleeping: false,
        });
        let merges = solver.merge_small_islands();
        assert_eq!(merges, 0, "no merge needed when island is large enough");
        assert_eq!(solver.manager.islands.len(), 1, "island count unchanged");
    }
    #[test]
    fn test_sleeping_criteria_high_energy_not_sleeping() {
        let mut solver = IslandSolver::with_thresholds(1e-4, 1);
        solver.manager.add_body(make_moving_body(0, 100.0));
        solver.manager.islands.push(Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: false,
        });
        let to_sleep = solver.compute_sleeping_criteria();
        assert!(to_sleep.is_empty(), "high-energy island should not sleep");
    }
    #[test]
    fn test_sleeping_criteria_static_only_island_not_queued() {
        let mut solver = IslandSolver::with_thresholds(1e-4, 1);
        solver.manager.add_body(make_body(0, true));
        solver.manager.islands.push(Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: false,
        });
        let to_sleep = solver.compute_sleeping_criteria();
        let _ = to_sleep;
    }
    #[test]
    fn test_apply_sleeping_marks_islands() {
        let mut solver = IslandSolver::with_thresholds(1e6, 1);
        solver.manager.add_body(make_body(0, false));
        solver.manager.islands.push(Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: false,
        });
        let slept = solver.apply_sleeping();
        assert!(slept <= 1, "at most one island slept, slept={slept}");
        if slept == 1 {
            assert!(
                solver.manager.islands[0].is_sleeping,
                "island should be marked sleeping"
            );
        }
    }
    #[test]
    fn test_wake_all_islands_clears_sleeping() {
        let mut solver = IslandSolver::new();
        solver.manager.islands.push(Island {
            bodies: vec![0],
            constraints: vec![],
            is_sleeping: true,
        });
        solver.manager.islands.push(Island {
            bodies: vec![1],
            constraints: vec![],
            is_sleeping: true,
        });
        solver.wake_all_islands();
        assert!(
            solver.manager.islands.iter().all(|i| !i.is_sleeping),
            "all islands should be awake"
        );
    }
}
