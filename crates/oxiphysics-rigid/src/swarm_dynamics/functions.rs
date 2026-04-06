//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::SwarmAgent;

#[inline]
pub(super) fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}
#[inline]
pub(super) fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n > 1e-15 {
        vec3_scale(a, 1.0 / n)
    } else {
        [0.0, 0.0, 0.0]
    }
}
#[inline]
pub(super) fn vec3_dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    vec3_norm(vec3_sub(a, b))
}
#[inline]
pub(super) fn vec3_clamp_magnitude(v: [f64; 3], max_mag: f64) -> [f64; 3] {
    let n = vec3_norm(v);
    if n > max_mag && n > 1e-15 {
        vec3_scale(v, max_mag / n)
    } else {
        v
    }
}
/// Compute the separation steering force for `agent` given `neighbors`.
///
/// For each neighbor within `radius`, a repulsive force is accumulated that
/// is inversely proportional to the distance.  The result is a unit-mass
/// steering vector (not yet multiplied by any weight).
pub fn separation_force(agent: &SwarmAgent, neighbors: &[SwarmAgent], radius: f64) -> [f64; 3] {
    let mut steer = [0.0f64; 3];
    let mut count = 0usize;
    for n in neighbors {
        if n.id == agent.id {
            continue;
        }
        let d = vec3_dist(agent.position, n.position);
        if d < radius && d > 1e-9 {
            let diff = vec3_sub(agent.position, n.position);
            let diff_norm = vec3_scale(vec3_normalize(diff), 1.0 / d);
            steer = vec3_add(steer, diff_norm);
            count += 1;
        }
    }
    if count > 0 {
        steer = vec3_scale(steer, 1.0 / count as f64);
    }
    steer
}
/// Compute the average position (centroid) of `neighbors`.
///
/// Returns the zero vector when the slice is empty.
pub fn cohesion_center(neighbors: &[SwarmAgent]) -> [f64; 3] {
    if neighbors.is_empty() {
        return [0.0; 3];
    }
    let sum = neighbors
        .iter()
        .fold([0.0f64; 3], |acc, n| vec3_add(acc, n.position));
    vec3_scale(sum, 1.0 / neighbors.len() as f64)
}
/// Compute the alignment steering force: steer toward average neighbor heading.
pub(super) fn alignment_force(
    agent: &SwarmAgent,
    neighbors: &[SwarmAgent],
    radius: f64,
) -> [f64; 3] {
    let mut avg_vel = [0.0f64; 3];
    let mut count = 0usize;
    for n in neighbors {
        if n.id == agent.id {
            continue;
        }
        if vec3_dist(agent.position, n.position) < radius {
            avg_vel = vec3_add(avg_vel, n.velocity);
            count += 1;
        }
    }
    if count > 0 {
        avg_vel = vec3_scale(avg_vel, 1.0 / count as f64);
        avg_vel = vec3_normalize(avg_vel);
        vec3_sub(avg_vel, vec3_normalize(agent.velocity))
    } else {
        [0.0; 3]
    }
}
/// Cohesion force: steer toward average position of nearby neighbors.
pub(super) fn cohesion_force(
    agent: &SwarmAgent,
    neighbors: &[SwarmAgent],
    radius: f64,
) -> [f64; 3] {
    let nearby: Vec<&SwarmAgent> = neighbors
        .iter()
        .filter(|n| n.id != agent.id && vec3_dist(agent.position, n.position) < radius)
        .collect();
    if nearby.is_empty() {
        return [0.0; 3];
    }
    let center = nearby
        .iter()
        .fold([0.0f64; 3], |acc, n| vec3_add(acc, n.position));
    let center = vec3_scale(center, 1.0 / nearby.len() as f64);
    vec3_normalize(vec3_sub(center, agent.position))
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::swarm_dynamics::AgentState;
    use crate::swarm_dynamics::AntColonyOptimizer;
    use crate::swarm_dynamics::BoidFlocking;
    use crate::swarm_dynamics::BoidSimulation;
    use crate::swarm_dynamics::CollectiveTransport;
    use crate::swarm_dynamics::FlockingRules;
    use crate::swarm_dynamics::FormationControl;
    use crate::swarm_dynamics::FormationShape;
    use crate::swarm_dynamics::ParticleSwarmOptimizer;
    use crate::swarm_dynamics::PotentialFieldSwarm;
    use crate::swarm_dynamics::PsoParticle;
    use crate::swarm_dynamics::SwarmAgentFull;
    use crate::swarm_dynamics::SwarmFormationControl;
    use crate::swarm_dynamics::SwarmMetrics;
    fn agent_at(id: u32, pos: [f64; 3]) -> SwarmAgent {
        SwarmAgent::new(id, pos, 5.0, 10.0)
    }
    fn agent_with_vel(id: u32, pos: [f64; 3], vel: [f64; 3]) -> SwarmAgent {
        let mut a = agent_at(id, pos);
        a.velocity = vel;
        a
    }
    #[test]
    fn agent_initial_speed_is_zero() {
        let a = agent_at(0, [0.0; 3]);
        assert!(a.speed() < 1e-12);
    }
    #[test]
    fn agent_heading_zero_when_stationary() {
        let a = agent_at(0, [1.0, 2.0, 3.0]);
        let h = a.heading();
        assert!(vec3_norm(h) < 1e-12);
    }
    #[test]
    fn agent_heading_unit_when_moving() {
        let mut a = agent_at(0, [0.0; 3]);
        a.velocity = [3.0, 0.0, 4.0];
        let h = a.heading();
        assert!((vec3_norm(h) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn apply_steering_moves_agent() {
        let mut a = agent_at(0, [0.0; 3]);
        a.apply_steering([1.0, 0.0, 0.0], 1.0);
        assert!(a.position[0] > 0.0);
    }
    #[test]
    fn apply_steering_clamps_to_max_speed() {
        let mut a = agent_at(0, [0.0; 3]);
        for _ in 0..100 {
            a.apply_steering([1000.0, 0.0, 0.0], 0.1);
        }
        assert!(a.speed() <= a.max_speed + 1e-9);
    }
    #[test]
    fn apply_steering_clamps_force() {
        let mut a = agent_at(0, [0.0; 3]);
        a.apply_steering([1000.0, 0.0, 0.0], 0.01);
        assert!((a.velocity[0] - 0.1).abs() < 1e-9);
    }
    #[test]
    fn separation_force_zero_with_no_neighbors() {
        let a = agent_at(0, [0.0; 3]);
        let f = separation_force(&a, &[], 3.0);
        assert!(vec3_norm(f) < 1e-12);
    }
    #[test]
    fn separation_force_pushes_away_from_neighbor() {
        let a = agent_at(0, [0.0, 0.0, 0.0]);
        let n = agent_at(1, [1.0, 0.0, 0.0]);
        let f = separation_force(&a, &[n], 3.0);
        assert!(f[0] < 0.0, "separation should push in -x: {:?}", f);
    }
    #[test]
    fn separation_force_ignores_distant_neighbors() {
        let a = agent_at(0, [0.0; 3]);
        let far = agent_at(1, [100.0, 0.0, 0.0]);
        let f = separation_force(&a, &[far], 3.0);
        assert!(vec3_norm(f) < 1e-12);
    }
    #[test]
    fn separation_force_ignores_self() {
        let a = agent_at(0, [0.0; 3]);
        let same = agent_at(0, [1.0, 0.0, 0.0]);
        let f = separation_force(&a, &[same], 3.0);
        assert!(vec3_norm(f) < 1e-12);
    }
    #[test]
    fn separation_force_symmetric_two_agents() {
        let a = agent_at(0, [0.0, 0.0, 0.0]);
        let b = agent_at(1, [1.0, 0.0, 0.0]);
        let fa = separation_force(&a, &[a.clone(), b.clone()], 3.0);
        let fb = separation_force(&b, &[a.clone(), b.clone()], 3.0);
        assert!(
            fa[0] * fb[0] < 0.0,
            "forces should be in opposite x directions"
        );
    }
    #[test]
    fn separation_force_stronger_when_closer() {
        let a = agent_at(0, [0.0; 3]);
        let close = agent_at(1, [0.5, 0.0, 0.0]);
        let far = agent_at(2, [2.5, 0.0, 0.0]);
        let f_close = separation_force(&a, &[close], 5.0);
        let f_far = separation_force(&a, &[far], 5.0);
        assert!(vec3_norm(f_close) > vec3_norm(f_far));
    }
    #[test]
    fn cohesion_center_empty_returns_zero() {
        let c = cohesion_center(&[]);
        assert!(vec3_norm(c) < 1e-12);
    }
    #[test]
    fn cohesion_center_single_agent() {
        let a = agent_at(0, [3.0, 4.0, 5.0]);
        let c = cohesion_center(&[a]);
        assert!((c[0] - 3.0).abs() < 1e-9);
        assert!((c[1] - 4.0).abs() < 1e-9);
        assert!((c[2] - 5.0).abs() < 1e-9);
    }
    #[test]
    fn cohesion_center_two_agents_midpoint() {
        let a = agent_at(0, [0.0, 0.0, 0.0]);
        let b = agent_at(1, [4.0, 0.0, 0.0]);
        let c = cohesion_center(&[a, b]);
        assert!((c[0] - 2.0).abs() < 1e-9);
    }
    #[test]
    fn cohesion_center_three_agents() {
        let a = agent_at(0, [0.0, 0.0, 0.0]);
        let b = agent_at(1, [3.0, 0.0, 0.0]);
        let c = agent_at(2, [6.0, 0.0, 0.0]);
        let center = cohesion_center(&[a, b, c]);
        assert!((center[0] - 3.0).abs() < 1e-9);
    }
    #[test]
    fn order_parameter_zero_agents() {
        assert_eq!(SwarmMetrics::order_parameter(&[]), 0.0);
    }
    #[test]
    fn order_parameter_perfectly_aligned() {
        let agents = vec![
            agent_with_vel(0, [0.0; 3], [1.0, 0.0, 0.0]),
            agent_with_vel(1, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            agent_with_vel(2, [2.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        ];
        let op = SwarmMetrics::order_parameter(&agents);
        assert!(
            (op - 1.0).abs() < 1e-9,
            "perfect alignment => OP=1, got {op}"
        );
    }
    #[test]
    fn order_parameter_opposite_directions() {
        let agents = vec![
            agent_with_vel(0, [0.0; 3], [1.0, 0.0, 0.0]),
            agent_with_vel(1, [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]),
        ];
        let op = SwarmMetrics::order_parameter(&agents);
        assert!(op < 1e-9, "opposite headings => OP≈0, got {op}");
    }
    #[test]
    fn order_parameter_single_moving_agent() {
        let agents = vec![agent_with_vel(0, [0.0; 3], [2.0, 0.0, 0.0])];
        let op = SwarmMetrics::order_parameter(&agents);
        assert!((op - 1.0).abs() < 1e-9);
    }
    #[test]
    fn polarization_equals_order_parameter() {
        let agents = vec![
            agent_with_vel(0, [0.0; 3], [1.0, 0.0, 0.0]),
            agent_with_vel(1, [1.0, 0.0, 0.0], [0.5, 0.0, 0.0]),
        ];
        let op = SwarmMetrics::order_parameter(&agents);
        let pol = SwarmMetrics::polarization(&agents);
        assert!((op - pol).abs() < 1e-12);
    }
    #[test]
    fn dispersion_zero_all_at_same_point() {
        let agents = vec![agent_at(0, [1.0, 2.0, 3.0]), agent_at(1, [1.0, 2.0, 3.0])];
        assert!(SwarmMetrics::dispersion(&agents) < 1e-9);
    }
    #[test]
    fn dispersion_two_agents_symmetric() {
        let agents = vec![agent_at(0, [-1.0, 0.0, 0.0]), agent_at(1, [1.0, 0.0, 0.0])];
        let d = SwarmMetrics::dispersion(&agents);
        assert!((d - 1.0).abs() < 1e-9, "dispersion={d}");
    }
    #[test]
    fn dispersion_zero_no_agents() {
        assert_eq!(SwarmMetrics::dispersion(&[]), 0.0);
    }
    #[test]
    fn boid_simulation_step_does_not_panic() {
        let agents = vec![
            agent_at(0, [0.0, 0.0, 0.0]),
            agent_at(1, [1.0, 0.0, 0.0]),
            agent_at(2, [2.0, 0.0, 0.0]),
        ];
        let mut sim = BoidSimulation::new(agents, FlockingRules::default_rules());
        sim.step(0.016);
    }
    #[test]
    fn boid_simulation_len_and_is_empty() {
        let agents = vec![agent_at(0, [0.0; 3]), agent_at(1, [1.0, 0.0, 0.0])];
        let sim = BoidSimulation::new(agents, FlockingRules::default_rules());
        assert_eq!(sim.len(), 2);
        assert!(!sim.is_empty());
    }
    #[test]
    fn boid_simulation_empty_is_empty() {
        let sim = BoidSimulation::new(vec![], FlockingRules::default_rules());
        assert!(sim.is_empty());
    }
    #[test]
    fn boid_separation_keeps_agents_apart() {
        let agents = vec![agent_at(0, [0.0, 0.0, 0.0]), agent_at(1, [0.1, 0.0, 0.0])];
        let rules = FlockingRules {
            separation_radius: 2.0,
            alignment_radius: 0.0,
            cohesion_radius: 0.0,
            separation_weight: 5.0,
            alignment_weight: 0.0,
            cohesion_weight: 0.0,
        };
        let mut sim = BoidSimulation::new(agents, rules);
        let init_dist = vec3_dist(sim.agents[0].position, sim.agents[1].position);
        for _ in 0..20 {
            sim.step(0.016);
        }
        let final_dist = vec3_dist(sim.agents[0].position, sim.agents[1].position);
        assert!(
            final_dist >= init_dist,
            "separation should increase distance: {init_dist} -> {final_dist}"
        );
    }
    #[test]
    fn boid_cohesion_pulls_agents_together() {
        let agents = vec![agent_at(0, [0.0, 0.0, 0.0]), agent_at(1, [10.0, 0.0, 0.0])];
        let rules = FlockingRules {
            separation_radius: 0.0,
            alignment_radius: 0.0,
            cohesion_radius: 20.0,
            separation_weight: 0.0,
            alignment_weight: 0.0,
            cohesion_weight: 2.0,
        };
        let mut sim = BoidSimulation::new(agents, rules);
        let init_dist = vec3_dist(sim.agents[0].position, sim.agents[1].position);
        for _ in 0..50 {
            sim.step(0.016);
        }
        let final_dist = vec3_dist(sim.agents[0].position, sim.agents[1].position);
        assert!(
            final_dist <= init_dist,
            "cohesion should decrease distance: {init_dist} -> {final_dist}"
        );
    }
    #[test]
    fn boid_speeds_never_exceed_max() {
        let agents = (0..10).map(|i| agent_at(i, [i as f64, 0.0, 0.0])).collect();
        let mut sim = BoidSimulation::new(agents, FlockingRules::default_rules());
        for _ in 0..100 {
            sim.step(0.016);
        }
        for a in &sim.agents {
            assert!(
                a.speed() <= a.max_speed + 1e-9,
                "agent {} speed {} exceeds max {}",
                a.id,
                a.speed(),
                a.max_speed
            );
        }
    }
    #[test]
    fn potential_field_attracts_toward_goal() {
        let swarm = PotentialFieldSwarm::new([10.0, 0.0, 0.0], vec![], 1.0, 0.0, 1.0);
        let agent = agent_at(0, [0.0, 0.0, 0.0]);
        let f = swarm.compute_force(&agent);
        assert!(f[0] > 0.0, "force should point toward goal (+x): {:?}", f);
    }
    #[test]
    fn potential_field_repels_from_obstacle() {
        let swarm = PotentialFieldSwarm::new([0.0, 0.0, 0.0], vec![[1.0, 0.0, 0.0]], 0.0, 1.0, 5.0);
        let agent = agent_at(0, [0.5, 0.0, 0.0]);
        let f = swarm.compute_force(&agent);
        assert!(f[0] < 0.0, "force should repel in -x: {:?}", f);
    }
    #[test]
    fn potential_field_no_force_at_goal_no_obstacles() {
        let swarm = PotentialFieldSwarm::new([0.0, 0.0, 0.0], vec![], 1.0, 0.0, 1.0);
        let agent = agent_at(0, [0.0, 0.0, 0.0]);
        let f = swarm.compute_force(&agent);
        assert!(vec3_norm(f) < 1e-12);
    }
    #[test]
    fn potential_field_obstacle_outside_radius_no_repulsion() {
        let swarm =
            PotentialFieldSwarm::new([100.0, 0.0, 0.0], vec![[50.0, 0.0, 0.0]], 1.0, 5.0, 1.0);
        let agent = agent_at(0, [0.0, 0.0, 0.0]);
        let f = swarm.compute_force(&agent);
        assert!(f[0] > 0.0);
        assert!(f[1].abs() < 1e-9);
    }
    #[test]
    fn potential_field_multiple_obstacles() {
        let swarm = PotentialFieldSwarm::new(
            [0.0, 0.0, 0.0],
            vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            0.0,
            1.0,
            5.0,
        );
        let agent = agent_at(0, [0.0, 0.0, 0.0]);
        let f = swarm.compute_force(&agent);
        assert!(f[0].abs() < 1e-9, "symmetric obstacles: {:?}", f);
    }
    #[test]
    fn formation_control_reduces_error() {
        let targets = vec![[5.0, 0.0, 0.0], [5.0, 2.0, 0.0]];
        let agents = vec![agent_at(0, [0.0; 3]), agent_at(1, [0.0, 0.0, 0.0])];
        let mut fc = FormationControl::new(targets, agents, 2.0);
        let init_err = fc.mean_position_error();
        for _ in 0..100 {
            fc.consensus_step(0.016);
        }
        let final_err = fc.mean_position_error();
        assert!(
            final_err < init_err,
            "error should decrease: {init_err} -> {final_err}"
        );
    }
    #[test]
    fn formation_control_len_and_is_empty() {
        let fc = FormationControl::new(vec![[1.0, 0.0, 0.0]], vec![agent_at(0, [0.0; 3])], 1.0);
        assert_eq!(fc.len(), 1);
        assert!(!fc.is_empty());
    }
    #[test]
    fn formation_control_empty_step_no_panic() {
        let mut fc = FormationControl::new(vec![], vec![], 1.0);
        fc.consensus_step(0.016);
        assert!(fc.is_empty());
    }
    #[test]
    fn formation_control_mean_error_zero_when_at_target() {
        let pos = [3.0, 4.0, 5.0];
        let fc = FormationControl::new(vec![pos], vec![agent_at(0, pos)], 1.0);
        assert!(fc.mean_position_error() < 1e-9);
    }
    #[test]
    fn formation_control_mismatched_lengths_uses_shorter() {
        let targets = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let agents = vec![agent_at(0, [0.0; 3])];
        let mut fc = FormationControl::new(targets, agents, 1.0);
        fc.consensus_step(0.016);
        assert_eq!(fc.len(), 1);
    }
    #[test]
    fn agent_clone_is_independent() {
        let a = agent_at(0, [1.0, 2.0, 3.0]);
        let mut b = a.clone();
        b.position[0] = 99.0;
        assert!((a.position[0] - 1.0).abs() < 1e-12);
        let _ = b.position[0];
    }
    #[test]
    fn vec3_helpers_roundtrip() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let sum = vec3_add(a, b);
        let diff = vec3_sub(sum, b);
        assert!((diff[0] - a[0]).abs() < 1e-12);
        assert!((diff[1] - a[1]).abs() < 1e-12);
        assert!((diff[2] - a[2]).abs() < 1e-12);
    }
    #[test]
    fn vec3_normalize_unit_vector() {
        let v = [3.0f64, 0.0, 4.0];
        let n = vec3_normalize(v);
        assert!((vec3_norm(n) - 1.0).abs() < 1e-9);
    }
    #[test]
    fn vec3_normalize_zero_returns_zero() {
        let n = vec3_normalize([0.0; 3]);
        assert!(vec3_norm(n) < 1e-12);
    }
    #[test]
    fn order_parameter_between_zero_and_one() {
        let agents = vec![
            agent_with_vel(0, [0.0; 3], [1.0, 0.0, 0.0]),
            agent_with_vel(1, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            agent_with_vel(2, [2.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ];
        let op = SwarmMetrics::order_parameter(&agents);
        assert!((0.0..=1.0 + 1e-9).contains(&op), "OP={op}");
    }
    #[test]
    fn boid_simulation_many_agents_stable() {
        let agents: Vec<SwarmAgent> = (0..20)
            .map(|i| agent_at(i, [(i as f64) * 0.5, 0.0, 0.0]))
            .collect();
        let mut sim = BoidSimulation::new(agents, FlockingRules::default_rules());
        for _ in 0..50 {
            sim.step(0.016);
        }
        for a in &sim.agents {
            for &x in &a.position {
                assert!(x.is_finite(), "agent {} position not finite", a.id);
            }
        }
    }
    fn full_agent(id: u32, pos: [f64; 3]) -> SwarmAgentFull {
        SwarmAgentFull::new(id, pos, 5.0, 10.0)
    }
    #[test]
    fn agent_full_initial_energy_is_one() {
        let a = full_agent(0, [0.0; 3]);
        assert!((a.energy - 1.0).abs() < 1e-12);
    }
    #[test]
    fn agent_full_initial_state_is_idle() {
        let a = full_agent(0, [0.0; 3]);
        assert_eq!(a.state, AgentState::Idle);
    }
    #[test]
    fn agent_full_set_state_works() {
        let mut a = full_agent(0, [0.0; 3]);
        a.set_state(AgentState::Active);
        assert_eq!(a.state, AgentState::Active);
    }
    #[test]
    fn agent_full_depleted_when_energy_reaches_zero() {
        let mut a = full_agent(0, [0.0; 3]);
        a.energy_drain = 10.0;
        a.max_speed = 100.0;
        a.max_force = 100.0;
        for _ in 0..1000 {
            a.apply_steering([100.0, 0.0, 0.0], 0.01);
        }
        assert_eq!(
            a.state,
            AgentState::Depleted,
            "should be Depleted after energy runs out"
        );
    }
    #[test]
    fn agent_full_recovers_from_depleted() {
        let mut a = full_agent(0, [0.0; 3]);
        a.state = AgentState::Depleted;
        a.energy = 0.0;
        a.energy_recovery = 1.0;
        for _ in 0..100 {
            a.apply_steering([0.0; 3], 0.01);
        }
        assert!(a.energy > 0.0, "energy should recover, got {}", a.energy);
    }
    #[test]
    fn agent_full_heading_updated_from_velocity() {
        let mut a = full_agent(0, [0.0; 3]);
        a.apply_steering([10.0, 0.0, 0.0], 0.1);
        let h = a.heading;
        let n = (h[0] * h[0] + h[1] * h[1] + h[2] * h[2]).sqrt();
        assert!(
            (n - 1.0).abs() < 1e-9 || n < 1e-9,
            "heading should be unit or zero, got {n}"
        );
        if n > 1e-9 {
            assert!(h[0] > 0.0, "heading should point in +x after +x force");
        }
    }
    fn boid_flock(n: usize) -> BoidFlocking {
        let agents: Vec<SwarmAgentFull> = (0..n)
            .map(|i| full_agent(i as u32, [(i as f64) * 1.0, 0.0, 0.0]))
            .collect();
        BoidFlocking::new(agents)
    }
    #[test]
    fn boid_flocking_step_no_panic() {
        let mut flock = boid_flock(5);
        flock.step(0.016);
    }
    #[test]
    fn boid_flocking_len_is_empty() {
        let flock = boid_flock(3);
        assert_eq!(flock.len(), 3);
        assert!(!flock.is_empty());
        let empty = BoidFlocking::new(vec![]);
        assert!(empty.is_empty());
    }
    #[test]
    fn boid_flocking_obstacle_avoidance_pushes_away() {
        let mut flock = BoidFlocking::new(vec![full_agent(0, [0.5, 0.0, 0.0])]);
        flock.add_obstacle([1.0, 0.0, 0.0]);
        flock.w_avoid = 10.0;
        flock.w_sep = 0.0;
        flock.w_ali = 0.0;
        flock.w_coh = 0.0;
        let pos_before = flock.agents[0].position[0];
        for _ in 0..20 {
            flock.step(0.016);
        }
        assert!(
            flock.agents[0].position[0] <= pos_before + 1e-9,
            "agent should be pushed away from obstacle (in -x direction)"
        );
    }
    #[test]
    fn boid_flocking_positions_finite() {
        let mut flock = boid_flock(10);
        flock.w_sep = 1.5;
        flock.w_ali = 1.0;
        flock.w_coh = 0.8;
        for _ in 0..100 {
            flock.step(0.016);
        }
        for a in &flock.agents {
            for &x in &a.position {
                assert!(
                    x.is_finite(),
                    "position component should be finite, got {x}"
                );
            }
        }
    }
    #[test]
    fn boid_flocking_separation_increases_distance() {
        let agents = vec![
            full_agent(0, [0.0, 0.0, 0.0]),
            full_agent(1, [0.1, 0.0, 0.0]),
        ];
        let mut flock = BoidFlocking::new(agents);
        flock.w_sep = 5.0;
        flock.w_ali = 0.0;
        flock.w_coh = 0.0;
        flock.w_avoid = 0.0;
        let d0 = vec3_dist(flock.agents[0].position, flock.agents[1].position);
        for _ in 0..30 {
            flock.step(0.016);
        }
        let d1 = vec3_dist(flock.agents[0].position, flock.agents[1].position);
        assert!(
            d1 >= d0,
            "separation should increase distance: {d0} -> {d1}"
        );
    }
    #[test]
    fn pso_finds_minimum_sphere_function() {
        let mut pso = ParticleSwarmOptimizer::new(30, 3, -5.0, 5.0);
        pso.run(100, &|x: &[f64]| x.iter().map(|v| v * v).sum::<f64>());
        let best = pso.global_best_value;
        assert!(best < 1.0, "PSO should find minimum near 0, got {best}");
    }
    #[test]
    fn pso_global_best_value_non_negative() {
        let mut pso = ParticleSwarmOptimizer::new(10, 2, -2.0, 2.0);
        pso.run(20, &|x: &[f64]| x.iter().map(|v| v * v).sum::<f64>());
        assert!(pso.global_best_value >= 0.0);
    }
    #[test]
    fn pso_iteration_count_matches() {
        let mut pso = ParticleSwarmOptimizer::new(5, 2, -1.0, 1.0);
        let n = 10;
        pso.run(n, &|x: &[f64]| x.iter().map(|v| v * v).sum::<f64>());
        assert_eq!(pso.iteration, n, "iteration should equal n after run(n)");
    }
    #[test]
    fn pso_mean_distance_to_best_non_negative() {
        let mut pso = ParticleSwarmOptimizer::new(10, 2, -3.0, 3.0);
        pso.run(50, &|x: &[f64]| x.iter().map(|v| v * v).sum::<f64>());
        let dist = pso.mean_distance_to_best();
        assert!(dist >= 0.0, "mean distance should be >= 0, got {dist}");
    }
    #[test]
    fn pso_particle_personal_best_updates() {
        let mut p = PsoParticle::new(vec![10.0, 10.0]);
        p.update_personal_best(100.0);
        assert_eq!(p.best_value, 100.0);
        p.update_personal_best(50.0);
        assert_eq!(p.best_value, 50.0);
        p.update_personal_best(200.0);
        assert_eq!(
            p.best_value, 50.0,
            "best should not increase on worse value"
        );
    }
    #[test]
    fn pso_convergence_decreases_distance() {
        let mut pso = ParticleSwarmOptimizer::new(20, 2, -5.0, 5.0);
        let obj = |x: &[f64]| x.iter().map(|v| v * v).sum::<f64>();
        pso.run(10, &obj);
        let d_early = pso.mean_distance_to_best();
        pso.run(100, &obj);
        let d_late = pso.mean_distance_to_best();
        assert!(
            d_late <= d_early + 1e-6,
            "swarm should converge: {d_early} -> {d_late}"
        );
    }
    fn small_distance_matrix() -> Vec<Vec<f64>> {
        vec![
            vec![0.0, 1.0, 1.414, 1.0],
            vec![1.0, 0.0, 1.0, 1.414],
            vec![1.414, 1.0, 0.0, 1.0],
            vec![1.0, 1.414, 1.0, 0.0],
        ]
    }
    #[test]
    fn aco_tour_length_correct_known() {
        let dm = small_distance_matrix();
        let aco = AntColonyOptimizer::new(dm, 5, 1.0, 2.0, 0.5);
        let tour = vec![0, 1, 2, 3];
        let len = aco.tour_length(&tour);
        assert!(
            (len - 4.0).abs() < 1e-9,
            "expected tour length 4.0, got {len}"
        );
    }
    #[test]
    fn aco_best_length_non_negative() {
        let dm = small_distance_matrix();
        let mut aco = AntColonyOptimizer::new(dm, 5, 1.0, 2.0, 0.5);
        aco.run(20);
        assert!(
            aco.best_length > 0.0,
            "best length should be positive, got {}",
            aco.best_length
        );
    }
    #[test]
    fn aco_best_length_improves_or_stays_same() {
        let dm = small_distance_matrix();
        let mut aco = AntColonyOptimizer::new(dm, 10, 1.0, 2.0, 0.1);
        aco.run(5);
        let mid_best = aco.best_length;
        aco.run(20);
        let final_best = aco.best_length;
        assert!(
            final_best <= mid_best + 1e-9,
            "best length should not get worse: {mid_best} -> {final_best}"
        );
    }
    #[test]
    fn aco_best_tour_covers_all_nodes() {
        let dm = small_distance_matrix();
        let mut aco = AntColonyOptimizer::new(dm, 5, 1.0, 2.0, 0.5);
        aco.run(10);
        let n = aco.n_nodes;
        assert_eq!(aco.best_tour.len(), n, "best tour should have all nodes");
        let mut seen = vec![false; n];
        for &node in &aco.best_tour {
            seen[node] = true;
        }
        assert!(
            seen.iter().all(|&s| s),
            "all nodes should appear in best tour"
        );
    }
    #[test]
    fn aco_pheromone_non_negative_after_evaporation() {
        let dm = small_distance_matrix();
        let mut aco = AntColonyOptimizer::new(dm, 5, 1.0, 2.0, 0.9);
        aco.run(50);
        for row in &aco.pheromone {
            for &tau in row {
                assert!(tau >= 0.0, "pheromone should be >= 0, got {tau}");
            }
        }
    }
    #[test]
    fn aco_mean_pheromone_positive() {
        let dm = small_distance_matrix();
        let mut aco = AntColonyOptimizer::new(dm, 5, 1.0, 2.0, 0.5);
        aco.run(5);
        assert!(aco.mean_pheromone() > 0.0);
    }
    fn make_full_agents(n: usize) -> Vec<SwarmAgentFull> {
        (0..n)
            .map(|i| full_agent(i as u32, [0.0, 0.0, i as f64 * 0.5]))
            .collect()
    }
    #[test]
    fn formation_line_desired_positions_count() {
        let fc = SwarmFormationControl::new(
            make_full_agents(5),
            FormationShape::Line,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            2.0,
            0,
        );
        assert_eq!(fc.desired_positions().len(), 5);
    }
    #[test]
    fn formation_v_desired_positions_count() {
        let fc = SwarmFormationControl::new(
            make_full_agents(7),
            FormationShape::V,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            2.0,
            0,
        );
        assert_eq!(fc.desired_positions().len(), 7);
    }
    #[test]
    fn formation_grid_desired_positions_count() {
        let fc = SwarmFormationControl::new(
            make_full_agents(9),
            FormationShape::Grid,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            2.0,
            0,
        );
        assert_eq!(fc.desired_positions().len(), 9);
    }
    #[test]
    fn formation_error_zero_when_at_target() {
        let mut agents = make_full_agents(4);
        let tmp_fc = SwarmFormationControl::new(
            agents.clone(),
            FormationShape::Line,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            2.0,
            1.0,
            0,
        );
        let targets = tmp_fc.desired_positions();
        for (a, t) in agents.iter_mut().zip(targets.iter()) {
            a.position = *t;
        }
        let fc = SwarmFormationControl::new(
            agents,
            FormationShape::Line,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            2.0,
            1.0,
            0,
        );
        assert!(
            fc.formation_error() < 1e-9,
            "error should be zero when at targets"
        );
    }
    #[test]
    fn formation_step_reduces_error() {
        let agents = make_full_agents(4);
        let mut fc = SwarmFormationControl::new(
            agents,
            FormationShape::Line,
            [10.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            2.0,
            3.0,
            0,
        );
        let e0 = fc.formation_error();
        for _ in 0..100 {
            fc.step(0.016);
        }
        let e1 = fc.formation_error();
        assert!(
            e1 < e0 + 1e-6,
            "formation error should decrease: {e0} -> {e1}"
        );
    }
    #[test]
    fn formation_empty_no_panic() {
        let mut fc = SwarmFormationControl::new(
            vec![],
            FormationShape::Grid,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            1.0,
            1.0,
            0,
        );
        fc.step(0.016);
        assert!(fc.is_empty());
    }
    #[test]
    fn formation_line_centred_correctly() {
        let fc = SwarmFormationControl::new(
            make_full_agents(3),
            FormationShape::Line,
            [0.0; 3],
            [1.0, 0.0, 0.0],
            2.0,
            1.0,
            0,
        );
        let pos = fc.desired_positions();
        let zs: Vec<f64> = pos.iter().map(|p| p[2]).collect();
        let min_z = zs.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_z = zs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            (max_z - min_z - 4.0).abs() < 1e-9,
            "line span should be 4m, got {}",
            max_z - min_z
        );
    }
    #[test]
    fn collective_transport_initial_distance_to_goal() {
        let ct = CollectiveTransport::new(4, [0.0; 3], 1.0, [10.0, 0.0, 0.0], 1.0);
        let d = ct.distance_to_goal();
        assert!(
            (d - 10.0).abs() < 1e-9,
            "initial distance should be 10, got {d}"
        );
    }
    #[test]
    fn collective_transport_payload_moves_toward_goal() {
        let mut ct = CollectiveTransport::new(4, [0.0; 3], 1.0, [5.0, 0.0, 0.0], 0.5);
        let d0 = ct.distance_to_goal();
        for _ in 0..200 {
            ct.step(0.016);
        }
        let d1 = ct.distance_to_goal();
        assert!(d1 < d0, "payload should move toward goal: {d0} -> {d1}");
    }
    #[test]
    fn collective_transport_load_variance_zero_after_rebalance() {
        let mut ct = CollectiveTransport::new(4, [0.0; 3], 1.0, [5.0, 0.0, 0.0], 1.0);
        ct.slots[0].load_fraction = 0.5;
        ct.slots[1].load_fraction = 0.1;
        ct.slots[2].load_fraction = 0.3;
        ct.slots[3].load_fraction = 0.1;
        ct.rebalance_loads();
        assert!(
            ct.load_variance() < 1e-12,
            "variance should be 0 after rebalance"
        );
    }
    #[test]
    fn collective_transport_carrier_count() {
        let ct = CollectiveTransport::new(6, [0.0; 3], 2.0, [1.0, 0.0, 0.0], 1.0);
        assert_eq!(ct.agents.len(), 6);
        assert_eq!(ct.slots.len(), 6);
    }
    #[test]
    fn collective_transport_step_no_nan() {
        let mut ct = CollectiveTransport::new(4, [0.0; 3], 1.0, [3.0, 0.0, 0.0], 0.5);
        for _ in 0..50 {
            ct.step(0.016);
        }
        for &x in &ct.payload {
            assert!(x.is_finite(), "payload position should be finite, got {x}");
        }
    }
    #[test]
    fn collective_transport_zero_carriers_no_panic() {
        let mut ct = CollectiveTransport::new(0, [0.0; 3], 1.0, [5.0, 0.0, 0.0], 1.0);
        ct.step(0.016);
        assert_eq!(ct.agents.len(), 0);
    }
}
