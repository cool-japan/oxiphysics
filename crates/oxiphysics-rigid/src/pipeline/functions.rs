//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::BodySnapshot;

#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
#[inline]
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n > 1e-12 {
        scale3(a, 1.0 / n)
    } else {
        [0.0, 1.0, 0.0]
    }
}
/// Compute an AABB for a body represented as a sphere (bounding sphere).
///
/// The sphere radius is `1.0 / inv_mass.cbrt() * 0.3` (heuristic) clamped to
/// at least 0.1 m.  For static bodies (inv_mass == 0) the radius is 0.5 m.
pub(super) fn body_aabb(pos: [f64; 3], inv_mass: f64) -> ([f64; 3], [f64; 3]) {
    let r = if inv_mass > 0.0 {
        (1.0 / inv_mass).cbrt() * 0.3
    } else {
        0.5_f64
    }
    .max(0.1);
    (
        [pos[0] - r, pos[1] - r, pos[2] - r],
        [pos[0] + r, pos[1] + r, pos[2] + r],
    )
}
/// Returns `true` if two AABBs overlap.
pub(super) fn aabb_overlap(
    min_a: [f64; 3],
    max_a: [f64; 3],
    min_b: [f64; 3],
    max_b: [f64; 3],
) -> bool {
    min_a[0] <= max_b[0]
        && max_a[0] >= min_b[0]
        && min_a[1] <= max_b[1]
        && max_a[1] >= min_b[1]
        && min_a[2] <= max_b[2]
        && max_a[2] >= min_b[2]
}
/// Builds a union-find parent array for `n` nodes connected by `edges`.
pub fn union_find_build(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..n).collect();
    for &(a, b) in edges {
        union_find_union(&mut parent, a, b);
    }
    parent
}
/// Path-compressed find: returns the root of node `i`.
pub fn union_find_root(parent: &mut Vec<usize>, i: usize) -> usize {
    if parent[i] != i {
        parent[i] = union_find_root(parent, parent[i]);
    }
    parent[i]
}
/// Unions the sets containing `a` and `b`.
pub fn union_find_union(parent: &mut Vec<usize>, a: usize, b: usize) {
    let ra = union_find_root(parent, a);
    let rb = union_find_root(parent, b);
    if ra != rb {
        parent[rb] = ra;
    }
}
#[cfg(test)]
mod tests {

    use crate::pipeline::BodySnapshot;
    use crate::pipeline::BroadphaseDetector;

    use crate::pipeline::Contact;

    use crate::pipeline::IslandManager;
    use crate::pipeline::NarrowphaseDetector;

    use crate::pipeline::PhysicsPipeline;
    use crate::pipeline::PhysicsPipelineConfig;
    use crate::pipeline::PipelineStep;

    use crate::pipeline::SequentialImpulseSolver;
    use crate::pipeline::SleepManager;

    fn make_body(id: u64, pos: [f64; 3], vel: [f64; 3], inv_mass: f64) -> BodySnapshot {
        BodySnapshot {
            id,
            position: pos,
            velocity: vel,
            ang_velocity: [0.0; 3],
            inv_mass,
            active: true,
        }
    }
    fn default_config() -> PhysicsPipelineConfig {
        PhysicsPipelineConfig::new()
    }
    #[test]
    fn test_island_chain_3_bodies() {
        let bodies = vec![10u64, 20, 30];
        let contacts = vec![(10u64, 20u64), (20u64, 30u64)];
        let mgr = IslandManager::build_islands(&bodies, &contacts);
        assert_eq!(mgr.island_count(), 1);
        assert_eq!(mgr.bodies_in_island(0).len(), 3);
    }
    #[test]
    fn test_island_no_contacts_4_bodies() {
        let bodies = vec![1u64, 2, 3, 4];
        let contacts: Vec<(u64, u64)> = vec![];
        let mgr = IslandManager::build_islands(&bodies, &contacts);
        assert_eq!(mgr.island_count(), 4);
    }
    #[test]
    fn test_island_merge() {
        let bodies = vec![1u64, 2, 3, 4];
        let contacts: Vec<(u64, u64)> = vec![];
        let mut mgr = IslandManager::build_islands(&bodies, &contacts);
        assert_eq!(mgr.island_count(), 4);
        mgr.merge_islands(0, 1);
        assert_eq!(mgr.island_count(), 3);
    }
    #[test]
    fn test_island_two_separate_pairs() {
        let bodies = vec![0u64, 1, 2, 3];
        let contacts = vec![(0u64, 1u64), (2u64, 3u64)];
        let mgr = IslandManager::build_islands(&bodies, &contacts);
        assert_eq!(mgr.island_count(), 2);
    }
    #[test]
    fn test_apply_gravity_dynamic_body() {
        let config = default_config();
        let dt = config.dt;
        let g_y = config.gravity[1];
        let step = PipelineStep::new(config);
        let mut bodies = vec![make_body(1, [0.0; 3], [0.0; 3], 1.0)];
        step.apply_gravity(&mut bodies);
        let expected_vy = g_y * dt;
        assert!((bodies[0].velocity[1] - expected_vy).abs() < 1e-12);
        assert!((bodies[0].velocity[0]).abs() < 1e-12);
        assert!((bodies[0].velocity[2]).abs() < 1e-12);
    }
    #[test]
    fn test_static_body_unaffected_by_gravity() {
        let config = default_config();
        let step = PipelineStep::new(config);
        let mut bodies = vec![make_body(2, [0.0; 3], [0.0; 3], 0.0)];
        step.apply_gravity(&mut bodies);
        assert_eq!(bodies[0].velocity, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_integrate_positions_euler_step() {
        let config = default_config();
        let dt = config.dt;
        let step = PipelineStep::new(config);
        let vx = 2.0_f64;
        let vy = -3.0_f64;
        let vz = 1.0_f64;
        let mut bodies = vec![BodySnapshot {
            id: 3,
            position: [1.0, 2.0, 3.0],
            velocity: [vx, vy, vz],
            ang_velocity: [0.0; 3],
            inv_mass: 1.0,
            active: true,
        }];
        step.integrate_positions(&mut bodies);
        assert!((bodies[0].position[0] - (1.0 + vx * dt)).abs() < 1e-12);
        assert!((bodies[0].position[1] - (2.0 + vy * dt)).abs() < 1e-12);
        assert!((bodies[0].position[2] - (3.0 + vz * dt)).abs() < 1e-12);
    }
    #[test]
    fn test_broadphase_overlap_detected() {
        let bodies = vec![
            make_body(0, [0.0; 3], [0.0; 3], 1.0),
            make_body(1, [0.1, 0.0, 0.0], [0.0; 3], 1.0),
        ];
        let pairs = BroadphaseDetector::detect(&bodies);
        assert!(!pairs.is_empty(), "expected at least one pair");
    }
    #[test]
    fn test_broadphase_no_overlap() {
        let bodies = vec![
            make_body(0, [0.0; 3], [0.0; 3], 1.0),
            make_body(1, [1000.0, 0.0, 0.0], [0.0; 3], 1.0),
        ];
        let pairs = BroadphaseDetector::detect(&bodies);
        assert!(pairs.is_empty(), "expected no pairs");
    }
    #[test]
    fn test_broadphase_empty_bodies() {
        let bodies: Vec<BodySnapshot> = vec![];
        let pairs = BroadphaseDetector::detect(&bodies);
        assert!(pairs.is_empty());
    }
    #[test]
    fn test_narrowphase_generates_contact() {
        let bodies = vec![
            make_body(0, [0.0; 3], [0.0; 3], 1.0),
            make_body(1, [0.05, 0.0, 0.0], [0.0; 3], 1.0),
        ];
        let pairs = vec![(0usize, 1usize)];
        let contacts = NarrowphaseDetector::generate_contacts(&bodies, &pairs);
        assert!(!contacts.is_empty(), "should produce a contact");
        assert!(contacts[0].depth > 0.0);
    }
    #[test]
    fn test_narrowphase_no_contact_far_apart() {
        let bodies = vec![
            make_body(0, [0.0; 3], [0.0; 3], 1.0),
            make_body(1, [100.0, 0.0, 0.0], [0.0; 3], 1.0),
        ];
        let pairs = vec![(0usize, 1usize)];
        let contacts = NarrowphaseDetector::generate_contacts(&bodies, &pairs);
        assert!(contacts.is_empty(), "should not produce a contact");
    }
    #[test]
    fn test_solver_resolves_approach() {
        let mut bodies = vec![
            make_body(0, [-0.05, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0),
            make_body(1, [0.05, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0),
        ];
        let contacts = vec![Contact {
            idx_a: 0,
            idx_b: 1,
            normal: [-1.0, 0.0, 0.0],
            depth: 0.1,
            contact_point: [0.0; 3],
        }];
        let solver = SequentialImpulseSolver::new(10, 0.0);
        solver.solve(&mut bodies, &contacts);
        let rel_vx = bodies[0].velocity[0] - bodies[1].velocity[0];
        let n = [-1.0_f64, 0.0, 0.0];
        let vn = rel_vx * n[0];
        assert!(vn >= -1e-6, "still approaching after solve: vn={vn}");
    }
    #[test]
    fn test_solver_no_impulse_for_separating() {
        let mut bodies = vec![
            make_body(0, [0.0; 3], [-1.0, 0.0, 0.0], 1.0),
            make_body(1, [0.1, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0),
        ];
        let vx0_before = bodies[0].velocity[0];
        let vx1_before = bodies[1].velocity[0];
        let contacts = vec![Contact {
            idx_a: 0,
            idx_b: 1,
            normal: [-1.0, 0.0, 0.0],
            depth: 0.01,
            contact_point: [0.05; 3],
        }];
        let solver = SequentialImpulseSolver::new(5, 0.5);
        solver.solve(&mut bodies, &contacts);
        assert!((bodies[0].velocity[0] - vx0_before).abs() < 1e-12);
        assert!((bodies[1].velocity[0] - vx1_before).abs() < 1e-12);
    }
    #[test]
    fn test_sleep_manager_puts_body_to_sleep() {
        let mgr = SleepManager::new(0.1, 0.1, 0.5);
        let mut bodies = vec![make_body(0, [0.0; 3], [0.0; 3], 1.0)];
        let mut timers = vec![0.0_f64];
        let dt = 0.016;
        let steps = ((0.5 / dt) as usize) + 5;
        for _ in 0..steps {
            mgr.update(&mut bodies, &mut timers, dt);
        }
        assert!(!bodies[0].active, "body should be sleeping");
        assert_eq!(SleepManager::sleeping_count(&bodies), 1);
    }
    #[test]
    fn test_sleep_manager_wakes_body_on_motion() {
        let mgr = SleepManager::new(0.1, 0.1, 0.5);
        let mut bodies = vec![make_body(0, [0.0; 3], [0.0; 3], 1.0)];
        let mut timers = vec![0.49_f64];
        bodies[0].active = false;
        bodies[0].velocity = [5.0, 0.0, 0.0];
        mgr.update(&mut bodies, &mut timers, 0.016);
        assert!(
            bodies[0].active,
            "body should be awake after velocity spike"
        );
    }
    #[test]
    fn test_pipeline_step_runs_without_bodies() {
        let mut pipeline = PhysicsPipeline::new(default_config());
        let mut bodies: Vec<BodySnapshot> = vec![];
        let report = pipeline.step(&mut bodies);
        assert_eq!(report.num_contacts, 0);
        assert_eq!(report.num_sleeping, 0);
    }
    #[test]
    fn test_pipeline_step_advances_velocity() {
        let mut config = default_config();
        config.gravity = [0.0, -10.0, 0.0];
        let dt = config.dt;
        let mut pipeline = PhysicsPipeline::new(config);
        let mut bodies = vec![make_body(0, [0.0; 3], [0.0; 3], 1.0)];
        pipeline.step(&mut bodies);
        assert!(
            bodies[0].velocity[1] < 0.0,
            "vy should be negative after gravity: vy={}",
            bodies[0].velocity[1]
        );
        let expected_vy = -10.0 * dt;
        assert!(
            (bodies[0].velocity[1] - expected_vy).abs() < 0.1,
            "vy={} expected≈{}",
            bodies[0].velocity[1],
            expected_vy
        );
    }
    #[test]
    fn test_pipeline_step_report_contact_count() {
        let mut pipeline = PhysicsPipeline::new(default_config());
        let mut bodies = vec![
            make_body(0, [0.0; 3], [1.0, 0.0, 0.0], 1.0),
            make_body(1, [0.05, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0),
        ];
        let report = pipeline.step(&mut bodies);
        assert!(report.num_contacts > 0, "expected contacts but got 0");
    }
    #[test]
    fn test_pipeline_step_report_timing() {
        let mut pipeline = PhysicsPipeline::new(default_config());
        let mut bodies = vec![make_body(0, [0.0; 3], [0.0; 3], 1.0)];
        let report = pipeline.step(&mut bodies);
        let _ = report.step_time_us;
    }
    #[test]
    fn test_pipeline_static_body_not_moved() {
        let mut config = default_config();
        config.gravity = [0.0, -9.81, 0.0];
        let mut pipeline = PhysicsPipeline::new(config);
        let mut bodies = vec![make_body(0, [0.0; 3], [0.0; 3], 0.0)];
        for _ in 0..10 {
            pipeline.step(&mut bodies);
        }
        assert_eq!(bodies[0].position, [0.0; 3], "static body should not move");
    }
}
#[cfg(test)]
mod tests_pipeline_ext {

    use crate::pipeline::BodySnapshot;

    use crate::pipeline::ContactManifold;

    use crate::pipeline::SubStepInterpolator;
    use crate::pipeline::VelocityContactConstraint;
    use crate::pipeline::WarmStartCache;
    use crate::pipeline::WarmStartEntry;
    #[test]
    fn test_warm_start_cache_insert_and_find() {
        let mut cache = WarmStartCache::new();
        let e = WarmStartEntry::new(1, 2, [0.0, 1.0, 0.0]);
        cache.upsert(e);
        assert!(cache.find(1, 2).is_some());
        assert!(cache.find(2, 1).is_some(), "order-independent lookup");
        assert!(cache.find(1, 3).is_none(), "missing pair");
    }
    #[test]
    fn test_warm_start_cache_upsert_updates() {
        let mut cache = WarmStartCache::new();
        let mut e = WarmStartEntry::new(0, 1, [1.0, 0.0, 0.0]);
        e.lambda_n = 5.0;
        cache.upsert(e);
        let mut e2 = WarmStartEntry::new(0, 1, [1.0, 0.0, 0.0]);
        e2.lambda_n = 10.0;
        cache.upsert(e2);
        assert_eq!(cache.len(), 1, "should not duplicate");
        assert!((cache.find(0, 1).unwrap().lambda_n - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_warm_start_cache_scale_all() {
        let mut cache = WarmStartCache::new();
        let mut e = WarmStartEntry::new(0, 1, [0.0, 1.0, 0.0]);
        e.lambda_n = 8.0;
        e.lambda_t1 = 2.0;
        cache.upsert(e);
        cache.scale_all(0.5);
        let found = cache.find(0, 1).unwrap();
        assert!((found.lambda_n - 4.0).abs() < 1e-12);
        assert!((found.lambda_t1 - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_warm_start_cache_prune() {
        let mut cache = WarmStartCache::new();
        cache.upsert(WarmStartEntry::new(0, 1, [0.0, 1.0, 0.0]));
        cache.upsert(WarmStartEntry::new(2, 3, [1.0, 0.0, 0.0]));
        cache.prune(&[(0, 1)]);
        assert_eq!(cache.len(), 1);
        assert!(cache.find(0, 1).is_some());
        assert!(cache.find(2, 3).is_none());
    }
    #[test]
    fn test_warm_start_entry_clamp_normal() {
        let mut e = WarmStartEntry::new(0, 1, [0.0, 1.0, 0.0]);
        e.lambda_n = -3.0;
        e.clamp_normal();
        assert_eq!(e.lambda_n, 0.0);
    }
    fn make_snapshot(id: u64, pos: [f64; 3]) -> BodySnapshot {
        BodySnapshot {
            id,
            position: pos,
            velocity: [0.0; 3],
            ang_velocity: [0.0; 3],
            inv_mass: 1.0,
            active: true,
        }
    }
    #[test]
    fn test_interpolator_at_alpha_zero_gives_prev() {
        let mut interp = SubStepInterpolator::new();
        interp.advance(vec![make_snapshot(0, [0.0, 0.0, 0.0])]);
        interp.advance(vec![make_snapshot(0, [2.0, 0.0, 0.0])]);
        let pos = interp.interpolate_position(0, 0.0).unwrap();
        assert!(
            (pos[0] - 0.0).abs() < 1e-12,
            "alpha=0 should give prev: {}",
            pos[0]
        );
    }
    #[test]
    fn test_interpolator_at_alpha_one_gives_curr() {
        let mut interp = SubStepInterpolator::new();
        interp.advance(vec![make_snapshot(0, [0.0, 0.0, 0.0])]);
        interp.advance(vec![make_snapshot(0, [2.0, 0.0, 0.0])]);
        let pos = interp.interpolate_position(0, 1.0).unwrap();
        assert!(
            (pos[0] - 2.0).abs() < 1e-12,
            "alpha=1 should give curr: {}",
            pos[0]
        );
    }
    #[test]
    fn test_interpolator_midpoint() {
        let mut interp = SubStepInterpolator::new();
        interp.advance(vec![make_snapshot(0, [0.0, 0.0, 0.0])]);
        interp.advance(vec![make_snapshot(0, [4.0, 0.0, 0.0])]);
        let pos = interp.interpolate_position(0, 0.5).unwrap();
        assert!(
            (pos[0] - 2.0).abs() < 1e-12,
            "midpoint should be 2: {}",
            pos[0]
        );
    }
    #[test]
    fn test_interpolator_out_of_range_returns_none() {
        let interp = SubStepInterpolator::new();
        assert!(interp.interpolate_position(0, 0.5).is_none());
    }
    #[test]
    fn test_interpolator_body_count() {
        let mut interp = SubStepInterpolator::new();
        interp.advance(vec![make_snapshot(0, [0.0; 3]), make_snapshot(1, [1.0; 3])]);
        assert_eq!(interp.body_count(), 2);
    }
    #[test]
    fn test_velocity_contact_effective_mass() {
        let c =
            VelocityContactConstraint::new(0, 1, [0.0, 1.0, 0.0], -1.0, -0.01, 0.0, 0.3, 0.5, 0.5);
        assert!(
            (c.effective_mass - 1.0).abs() < 1e-10,
            "eff mass = {}",
            c.effective_mass
        );
    }
    #[test]
    fn test_velocity_contact_baumgarte_bias() {
        let c =
            VelocityContactConstraint::new(0, 1, [0.0, 1.0, 0.0], -1.0, -0.02, 0.0, 0.3, 1.0, 1.0);
        let bias = c.baumgarte_bias(1.0 / 60.0, 0.2, 0.005);
        assert!(bias > 0.0, "bias should be positive: {bias}");
    }
    #[test]
    fn test_velocity_contact_normal_impulse_clamp() {
        let mut c =
            VelocityContactConstraint::new(0, 1, [0.0, 1.0, 0.0], -2.0, -0.01, 0.5, 0.3, 1.0, 1.0);
        let delta = c.solve_normal(-2.0, 0.0);
        assert!(
            c.lambda_n >= 0.0,
            "normal impulse must be non-negative: {}",
            c.lambda_n
        );
        assert!(delta.is_finite());
    }
    #[test]
    fn test_contact_manifold_total_impulse() {
        let mut manifold = ContactManifold::new();
        let mut c1 =
            VelocityContactConstraint::new(0, 1, [0.0, 1.0, 0.0], -1.0, -0.01, 0.0, 0.3, 1.0, 1.0);
        c1.lambda_n = 3.0;
        let mut c2 =
            VelocityContactConstraint::new(0, 1, [0.0, 1.0, 0.0], -1.0, -0.01, 0.0, 0.3, 1.0, 1.0);
        c2.lambda_n = 7.0;
        manifold.add(c1);
        manifold.add(c2);
        assert!((manifold.total_normal_impulse() - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_contact_manifold_max_penetration() {
        let mut manifold = ContactManifold::new();
        let mut c1 =
            VelocityContactConstraint::new(0, 1, [0.0, 1.0, 0.0], -1.0, -0.005, 0.0, 0.3, 1.0, 1.0);
        c1.penetration = -0.005;
        let mut c2 =
            VelocityContactConstraint::new(0, 1, [0.0, 1.0, 0.0], -1.0, -0.02, 0.0, 0.3, 1.0, 1.0);
        c2.penetration = -0.02;
        manifold.add(c1);
        manifold.add(c2);
        let max_p = manifold.max_penetration();
        assert!((max_p - 0.02).abs() < 1e-12, "max_p = {max_p}");
    }
    #[test]
    fn test_contact_manifold_clear() {
        let mut manifold = ContactManifold::new();
        manifold.add(VelocityContactConstraint::new(
            0,
            1,
            [0.0, 1.0, 0.0],
            -1.0,
            -0.01,
            0.0,
            0.3,
            1.0,
            1.0,
        ));
        manifold.clear();
        assert!(manifold.is_empty());
    }
}
/// Trait for objects that want to apply external forces or torques before the
/// physics integration step.
///
/// Implementors receive a mutable slice of body snapshots and the current
/// simulation time, and may modify velocities / positions directly.
pub trait PreStepHook {
    /// Called once per sub-step *before* gravity integration.
    fn pre_step(&self, bodies: &mut [BodySnapshot], sim_time: f64, dt: f64);
}
/// Trait for objects that want to read back body state after the physics step.
pub trait PostStepHook {
    /// Called once per sub-step *after* position integration.
    fn post_step(&self, bodies: &[BodySnapshot], sim_time: f64, dt: f64);
}
#[cfg(test)]
mod tests_pipeline_new {
    use super::*;
    use crate::pipeline::BodySnapshot;

    use crate::pipeline::ConstantForceHook;

    use crate::pipeline::DampingHook;
    use crate::pipeline::DeterministicStepper;
    use crate::pipeline::EnergyMonitorHook;
    use crate::pipeline::EventAwarePipeline;
    use crate::pipeline::EventQueue;

    use crate::pipeline::PhaseTimings;
    use crate::pipeline::PhysicsEvent;
    use crate::pipeline::PhysicsEventKind;
    use crate::pipeline::PhysicsPipeline;
    use crate::pipeline::PhysicsPipelineConfig;

    use crate::pipeline::ProfiledPipeline;

    fn make_dyn_body(id: u64, pos: [f64; 3], vel: [f64; 3]) -> BodySnapshot {
        BodySnapshot {
            id,
            position: pos,
            velocity: vel,
            ang_velocity: [0.0; 3],
            inv_mass: 1.0,
            active: true,
        }
    }
    #[test]
    fn test_deterministic_stepper_no_step_small_frame() {
        let mut stepper = DeterministicStepper::new(1.0 / 60.0, 4);
        let mut pipeline = PhysicsPipeline::new(PhysicsPipelineConfig::new());
        let mut bodies = vec![make_dyn_body(0, [0.0; 3], [0.0; 3])];
        let (steps, _alpha) = stepper.advance(0.001, &mut bodies, &mut pipeline);
        assert_eq!(steps, 0);
    }
    #[test]
    fn test_deterministic_stepper_one_step_per_frame() {
        let fixed_dt = 1.0 / 60.0;
        let mut stepper = DeterministicStepper::new(fixed_dt, 4);
        let mut pipeline = PhysicsPipeline::new(PhysicsPipelineConfig::new());
        let mut bodies = vec![make_dyn_body(0, [0.0; 3], [0.0; 3])];
        let (steps, alpha) = stepper.advance(fixed_dt, &mut bodies, &mut pipeline);
        assert_eq!(steps, 1);
        assert!(alpha < 1.0, "alpha = {alpha}");
        assert!((stepper.sim_time - fixed_dt).abs() < 1e-12);
    }
    #[test]
    fn test_deterministic_stepper_accumulator_carries_over() {
        let fixed_dt = 1.0 / 60.0;
        let mut stepper = DeterministicStepper::new(fixed_dt, 8);
        let mut pipeline = PhysicsPipeline::new(PhysicsPipelineConfig::new());
        let mut bodies = vec![make_dyn_body(0, [0.0; 3], [0.0; 3])];
        let (steps, _) = stepper.advance(fixed_dt * 2.5, &mut bodies, &mut pipeline);
        assert_eq!(steps, 2, "should take 2 full steps");
        assert!(stepper.accumulator > 0.0, "should carry over 0.5 fixed_dt");
    }
    #[test]
    fn test_deterministic_stepper_max_steps_cap() {
        let fixed_dt = 1.0 / 60.0;
        let mut stepper = DeterministicStepper::new(fixed_dt, 2);
        let mut pipeline = PhysicsPipeline::new(PhysicsPipelineConfig::new());
        let mut bodies = vec![make_dyn_body(0, [0.0; 3], [0.0; 3])];
        let (steps, _) = stepper.advance(fixed_dt * 10.0, &mut bodies, &mut pipeline);
        assert!(steps <= 2, "should be capped at 2, got {steps}");
    }
    #[test]
    fn test_deterministic_stepper_reset() {
        let mut stepper = DeterministicStepper::new(1.0 / 60.0, 4);
        stepper.accumulator = 1.0;
        stepper.sim_time = 5.0;
        stepper.total_steps = 100;
        stepper.reset();
        assert_eq!(stepper.accumulator, 0.0);
        assert_eq!(stepper.sim_time, 0.0);
        assert_eq!(stepper.total_steps, 0);
    }
    #[test]
    fn test_deterministic_stepper_is_ready() {
        let mut stepper = DeterministicStepper::new(1.0 / 60.0, 4);
        assert!(!stepper.is_ready());
        stepper.accumulator = 1.0 / 60.0;
        assert!(stepper.is_ready());
    }
    #[test]
    fn test_profiled_pipeline_step_runs() {
        let mut pp = ProfiledPipeline::new(PhysicsPipelineConfig::new());
        let mut bodies = vec![make_dyn_body(0, [0.0; 3], [0.0; 3])];
        let _report = pp.step(&mut bodies);
        assert_eq!(pp.step_count, 1);
        assert!(
            pp.last_timings.total_us < 1_000_000,
            "step shouldn't take > 1 s"
        );
    }
    #[test]
    fn test_profiled_pipeline_avg_step_us() {
        let mut pp = ProfiledPipeline::new(PhysicsPipelineConfig::new());
        let mut bodies = vec![make_dyn_body(0, [0.0; 3], [0.0; 3])];
        pp.step(&mut bodies);
        pp.step(&mut bodies);
        let avg = pp.avg_step_us();
        assert!(avg >= 0.0, "avg step time should be non-negative");
    }
    #[test]
    fn test_profiled_pipeline_reset_stats() {
        let mut pp = ProfiledPipeline::new(PhysicsPipelineConfig::new());
        let mut bodies = vec![make_dyn_body(0, [0.0; 3], [0.0; 3])];
        pp.step(&mut bodies);
        pp.reset_stats();
        assert_eq!(pp.step_count, 0);
        assert_eq!(pp.cumulative_timings.total_us, 0);
    }
    #[test]
    fn test_phase_timings_fractions() {
        let t = PhaseTimings {
            gravity_us: 1,
            broadphase_us: 2,
            narrowphase_us: 3,
            solver_us: 4,
            integrate_us: 0,
            sleep_us: 0,
            total_us: 10,
        };
        assert!((t.solver_fraction() - 0.4).abs() < 1e-10);
        assert!((t.collision_fraction() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_event_pipeline_collision_begin() {
        let mut pipeline = EventAwarePipeline::new(PhysicsPipelineConfig::new());
        let mut queue = EventQueue::new();
        let mut bodies = vec![
            make_dyn_body(10, [0.0; 3], [0.0; 3]),
            make_dyn_body(20, [0.05, 0.0, 0.0], [0.0; 3]),
        ];
        pipeline.step(&mut bodies, &mut queue);
        let count = queue.count_kind(PhysicsEventKind::CollisionBegin);
        assert!(count >= 1, "expected CollisionBegin, got {count}");
    }
    #[test]
    fn test_event_queue_drain() {
        let mut q = EventQueue::new();
        q.push(PhysicsEvent::body(PhysicsEventKind::BodySlept, 1, 0.0));
        q.push(PhysicsEvent::body(PhysicsEventKind::BodyWoke, 2, 0.0));
        assert_eq!(q.len(), 2);
        let drained = q.drain();
        assert_eq!(drained.len(), 2);
        assert!(q.is_empty());
    }
    #[test]
    fn test_event_queue_count_kind() {
        let mut q = EventQueue::new();
        q.push(PhysicsEvent::collision(
            PhysicsEventKind::CollisionBegin,
            1,
            2,
            0.0,
        ));
        q.push(PhysicsEvent::collision(
            PhysicsEventKind::CollisionEnd,
            1,
            2,
            0.1,
        ));
        q.push(PhysicsEvent::collision(
            PhysicsEventKind::CollisionBegin,
            3,
            4,
            0.1,
        ));
        assert_eq!(q.count_kind(PhysicsEventKind::CollisionBegin), 2);
        assert_eq!(q.count_kind(PhysicsEventKind::CollisionEnd), 1);
        assert_eq!(q.count_kind(PhysicsEventKind::BodySlept), 0);
    }
    #[test]
    fn test_constant_force_hook_accelerates_bodies() {
        let hook = ConstantForceHook {
            force: [10.0, 0.0, 0.0],
        };
        let mut bodies = vec![make_dyn_body(0, [0.0; 3], [0.0; 3])];
        hook.pre_step(&mut bodies, 0.0, 0.1);
        assert!(
            (bodies[0].velocity[0] - 1.0).abs() < 1e-12,
            "vx = {}",
            bodies[0].velocity[0]
        );
    }
    #[test]
    fn test_constant_force_hook_static_body_unaffected() {
        let hook = ConstantForceHook {
            force: [100.0, 0.0, 0.0],
        };
        let mut bodies = vec![BodySnapshot {
            id: 0,
            position: [0.0; 3],
            velocity: [0.0; 3],
            ang_velocity: [0.0; 3],
            inv_mass: 0.0,
            active: true,
        }];
        hook.pre_step(&mut bodies, 0.0, 0.1);
        assert_eq!(
            bodies[0].velocity[0], 0.0,
            "static body should not be moved"
        );
    }
    #[test]
    fn test_damping_hook_reduces_velocity() {
        let hook = DampingHook { coeff: 1.0 };
        let mut bodies = vec![make_dyn_body(0, [0.0; 3], [10.0, 0.0, 0.0])];
        hook.pre_step(&mut bodies, 0.0, 0.5);
        assert!(
            (bodies[0].velocity[0] - 5.0).abs() < 1e-12,
            "vx = {}",
            bodies[0].velocity[0]
        );
    }
    #[test]
    fn test_energy_monitor_hook_records_samples() {
        let hook = EnergyMonitorHook::new();
        let bodies = vec![make_dyn_body(0, [0.0; 3], [2.0, 0.0, 0.0])];
        hook.post_step(&bodies, 0.0, 0.016);
        let samples = hook.snapshot();
        assert_eq!(samples.len(), 1);
        assert!((samples[0] - 2.0).abs() < 1e-12, "KE = {}", samples[0]);
    }
}
