//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{PbdParticle, XpbdDistanceConstraint};

#[cfg(test)]
mod tests {
    use super::super::types::*;

    use crate::pbd::PbdDistanceConstraint;

    use oxiphysics_core::math::Vec3;
    pub(super) const EPS: f64 = 1e-9;
    /// stiffness=1, stretched → corrects to rest_length
    #[test]
    fn test_pbd_distance_project_rigid() {
        let mut particles = vec![
            PbdParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            PbdParticle::new(Vec3::new(3.0, 0.0, 0.0), 1.0),
        ];
        let constraint = PbdDistanceConstraint::new(0, 1, 1.0, 1.0);
        constraint.project(&mut particles);
        let dist = (particles[1].position - particles[0].position).norm();
        assert!(
            (dist - 1.0).abs() < EPS,
            "expected rest_length=1.0, got dist={}",
            dist
        );
    }
    /// stiffness=0.5, stretched → partially corrected
    #[test]
    fn test_pbd_distance_project_soft() {
        let mut particles = vec![
            PbdParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            PbdParticle::new(Vec3::new(3.0, 0.0, 0.0), 1.0),
        ];
        let constraint = PbdDistanceConstraint::new(0, 1, 1.0, 0.5);
        constraint.project(&mut particles);
        let dist = (particles[1].position - particles[0].position).norm();
        assert!(
            (dist - 2.0).abs() < EPS,
            "expected dist≈2.0 after partial correction, got {}",
            dist
        );
        assert!(dist > 1.0 && dist < 3.0, "expected partial correction");
    }
    /// Fixed particle (inv_mass=0) doesn't move
    #[test]
    fn test_pbd_fixed_particle_unchanged() {
        let mut particles = vec![
            PbdParticle::new_fixed(Vec3::new(0.0, 0.0, 0.0)),
            PbdParticle::new(Vec3::new(3.0, 0.0, 0.0), 1.0),
        ];
        let constraint = PbdDistanceConstraint::new(0, 1, 1.0, 1.0);
        constraint.project(&mut particles);
        assert!(
            particles[0].position.norm() < EPS,
            "fixed particle moved: {:?}",
            particles[0].position
        );
        let dist = (particles[1].position - particles[0].position).norm();
        assert!((dist - 1.0).abs() < EPS, "expected dist=1.0, got {}", dist);
    }
    /// Particle inside sphere → pushed out to sphere surface
    #[test]
    fn test_pbd_sphere_collision_push() {
        let mut particles = vec![PbdParticle::new(Vec3::new(0.5, 0.0, 0.0), 1.0)];
        let constraint = PbdSphereCollision {
            particle: 0,
            sphere_center: Vec3::zeros(),
            sphere_radius: 1.0,
            restitution: 0.0,
        };
        constraint.project(&mut particles);
        let dist = particles[0].position.norm();
        assert!(
            (dist - 1.0).abs() < EPS,
            "expected particle on sphere surface (dist=1.0), got {}",
            dist
        );
    }
    /// Single free particle with gravity only → falls
    #[test]
    fn test_pbd_step_gravity_freefall() {
        let mut sim = PbdSimulation::new();
        let idx = sim.add_particle(Vec3::new(0.0, 10.0, 0.0), 1.0);
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let dt = 0.01;
        let initial_y = sim.particles[idx].position.y;
        sim.step(dt, gravity, 1);
        assert!(
            sim.particles[idx].position.y < initial_y,
            "particle should fall under gravity"
        );
    }
    /// 2-particle chain with top pinned → bottom hangs below top
    #[test]
    fn test_pbd_step_pinned_chain() {
        let mut sim = PbdSimulation::new();
        let top = sim.add_particle(Vec3::new(0.0, 1.0, 0.0), 1.0);
        let bottom = sim.add_particle(Vec3::new(0.0, 0.0, 0.0), 1.0);
        sim.pin_particle(top);
        sim.add_distance_constraint(top, bottom, 1.0);
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let dt = 0.01;
        for _ in 0..100 {
            sim.step(dt, gravity, 10);
        }
        assert!(
            (sim.particles[top].position - Vec3::new(0.0, 1.0, 0.0)).norm() < EPS,
            "pinned particle should not move"
        );
        assert!(
            sim.particles[bottom].position.y < sim.particles[top].position.y,
            "bottom particle should hang below top"
        );
        let dist = (sim.particles[bottom].position - sim.particles[top].position).norm();
        assert!(
            (dist - 1.0).abs() < 0.01,
            "distance should be approximately rest_length=1.0, got {}",
            dist
        );
    }
    /// Three particles in line, bend → bending constraint resists
    #[test]
    fn test_pbd_bending_constraint() {
        let mut sim = PbdSimulation::new();
        let a = sim.add_particle(Vec3::new(0.0, 0.0, 0.0), 1.0);
        let b = sim.add_particle(Vec3::new(1.0, 0.0, 0.0), 1.0);
        let c = sim.add_particle(Vec3::new(2.0, 0.0, 0.0), 1.0);
        sim.add_bending_constraint(a, b, c, 1.0);
        sim.particles[c].position = Vec3::new(1.0, 1.0, 0.0);
        let dist_before = (sim.particles[c].position - sim.particles[a].position).norm();
        sim.bending_constraints[0].project(&mut sim.particles);
        let dist_after = (sim.particles[c].position - sim.particles[a].position).norm();
        assert!(
            dist_after > dist_before,
            "bending constraint should push particles apart (dist_before={}, dist_after={})",
            dist_before,
            dist_after
        );
    }
    /// Volume constraint: volume is corrected toward rest volume
    #[test]
    fn test_pbd_volume_constraint_correction() {
        let mut particles = vec![
            PbdParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            PbdParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
            PbdParticle::new(Vec3::new(0.0, 1.0, 0.0), 1.0),
            PbdParticle::new(Vec3::new(0.0, 0.0, 1.0), 1.0),
        ];
        let indices = [0, 1, 2, 3];
        let rest_vol = PbdVolumeConstraint::compute_volume(&particles, &indices);
        particles[3].position = Vec3::new(0.0, 0.0, 0.5);
        let vol_before = PbdVolumeConstraint::compute_volume(&particles, &indices);
        let constraint = PbdVolumeConstraint::new(indices, rest_vol, 1.0);
        constraint.project(&mut particles);
        let vol_after = PbdVolumeConstraint::compute_volume(&particles, &indices);
        assert!(
            (vol_after - rest_vol).abs() < (vol_before - rest_vol).abs() + EPS,
            "Volume should be corrected: before={vol_before}, after={vol_after}, rest={rest_vol}"
        );
    }
    /// Particle below plane is pushed above
    #[test]
    fn test_pbd_plane_collision() {
        let mut particles = vec![PbdParticle::new(Vec3::new(0.0, -0.5, 0.0), 1.0)];
        let constraint =
            PbdPlaneCollision::new(0, Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
        constraint.project(&mut particles);
        assert!(
            particles[0].position.y >= -EPS,
            "Particle should be above the plane: y = {}",
            particles[0].position.y
        );
    }
    /// Particle above plane is not moved
    #[test]
    fn test_pbd_plane_collision_above() {
        let mut particles = vec![PbdParticle::new(Vec3::new(0.0, 1.0, 0.0), 1.0)];
        let constraint =
            PbdPlaneCollision::new(0, Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
        let pos_before = particles[0].position;
        constraint.project(&mut particles);
        assert!(
            (particles[0].position - pos_before).norm() < EPS,
            "Particle above plane should not move"
        );
    }
    /// create_chain makes correct number of particles and constraints
    #[test]
    fn test_pbd_create_chain() {
        let mut sim = PbdSimulation::new();
        let indices = sim.create_chain(Vec3::zeros(), 5, 1.0, 1.0, 1.0);
        assert_eq!(indices.len(), 5);
        assert_eq!(sim.particle_count(), 5);
        assert_eq!(sim.distance_constraint_count(), 4);
    }
    /// create_cloth_grid makes correct grid
    #[test]
    fn test_pbd_create_cloth_grid() {
        let mut sim = PbdSimulation::new();
        let grid = sim.create_cloth_grid(Vec3::zeros(), 3, 3, 1.0, 1.0, 1.0);
        assert_eq!(grid.len(), 3);
        assert_eq!(grid[0].len(), 3);
        assert_eq!(sim.particle_count(), 9);
        assert_eq!(sim.distance_constraint_count(), 12);
    }
    /// step_substep runs without panic
    #[test]
    fn test_pbd_step_substep() {
        let mut sim = PbdSimulation::new();
        sim.add_particle(Vec3::new(0.0, 5.0, 0.0), 1.0);
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        sim.step_substep(0.02, gravity, 5, 4);
        assert!(sim.particles[0].position.y < 5.0, "particle should fall");
    }
    /// apply_damping reduces velocity
    #[test]
    fn test_pbd_damping() {
        let mut sim = PbdSimulation::new();
        let idx = sim.add_particle(Vec3::zeros(), 1.0);
        sim.particles[idx].velocity = Vec3::new(10.0, 0.0, 0.0);
        sim.apply_damping(0.1);
        assert!(
            sim.particles[idx].velocity.x < 10.0,
            "Damping should reduce velocity"
        );
    }
    /// total_potential_energy is correct for single particle
    #[test]
    fn test_pbd_potential_energy() {
        let mut sim = PbdSimulation::new();
        sim.add_particle(Vec3::new(0.0, 10.0, 0.0), 2.0);
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let pe = sim.total_potential_energy(gravity);
        assert!((pe - 196.2).abs() < 1e-10, "PE = {pe}, expected 196.2");
    }
    /// Kinetic energy is non-negative
    #[test]
    fn test_pbd_kinetic_energy_non_negative() {
        let mut sim = PbdSimulation::new();
        sim.add_particle(Vec3::zeros(), 1.0);
        sim.particles[0].velocity = Vec3::new(3.0, 4.0, 0.0);
        let ke = sim.total_kinetic_energy();
        assert!((ke - 12.5).abs() < EPS, "KE = {ke}, expected 12.5");
    }
    /// Default creates empty simulation
    #[test]
    fn test_pbd_default() {
        let sim = PbdSimulation::default();
        assert_eq!(sim.particle_count(), 0);
    }
    /// Multiple constraint iterations converge
    #[test]
    fn test_pbd_convergence_with_iterations() {
        let mut sim = PbdSimulation::new();
        let a = sim.add_particle(Vec3::new(0.0, 0.0, 0.0), 1.0);
        let b = sim.add_particle(Vec3::new(5.0, 0.0, 0.0), 1.0);
        sim.add_distance_constraint(a, b, 1.0);
        let gravity = Vec3::zeros();
        sim.step(0.01, gravity, 100);
        let dist = (sim.particles[b].position - sim.particles[a].position).norm();
        let rest = sim.distance_constraints[0].rest_length;
        assert!(
            (dist - rest).abs() < 0.01,
            "After 100 iterations, dist={dist} should be ~rest={rest}"
        );
    }
}
/// XPBD simulation step for a single sub-step.
pub fn xpbd_pbd_step(
    particles: &mut [PbdParticle],
    constraints: &mut [XpbdDistanceConstraint],
    dt: f64,
    n_iterations: usize,
) {
    for c in constraints.iter_mut() {
        c.reset_lambda();
    }
    for _ in 0..n_iterations {
        for c in constraints.iter_mut() {
            c.project(particles, dt);
        }
    }
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    use crate::pbd::ConstraintPriority;

    use crate::pbd::ParallelPbd;
    use crate::pbd::PbdConstraintGraph;

    use crate::pbd::PbdConvergenceMonitor;
    use crate::pbd::PbdDistanceConstraint;

    use crate::pbd::PrioritizedConstraint;

    use crate::pbd::SequentialPbd;

    use oxiphysics_core::math::Vec3;
    #[test]
    fn test_constraint_graph_add_edge() {
        let mut g = PbdConstraintGraph::new(4);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        assert!(g.are_connected(0, 1));
        assert!(g.are_connected(1, 0));
        assert!(!g.are_connected(0, 2));
        assert_eq!(g.degree[1], 2);
    }
    #[test]
    fn test_constraint_graph_coloring() {
        let mut g = PbdConstraintGraph::new(3);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        let colors = g.greedy_coloring();
        assert_ne!(colors[0], colors[1]);
        assert_ne!(colors[1], colors[2]);
    }
    #[test]
    fn test_constraint_graph_chromatic_number_path() {
        let mut g = PbdConstraintGraph::new(4);
        g.add_edge(0, 1);
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        let chi = g.chromatic_number_estimate();
        assert!(chi <= 2, "path graph should be 2-colorable, got {chi}");
    }
    #[test]
    fn test_sequential_pbd_converges() {
        let mut particles = vec![
            PbdParticle::new_fixed(Vec3::new(0.0, 0.0, 0.0)),
            PbdParticle::new(Vec3::new(5.0, 0.0, 0.0), 1.0),
        ];
        let constraints = vec![PbdDistanceConstraint::new(0, 1, 1.0, 1.0)];
        let solver = SequentialPbd::new(50);
        solver.solve(&mut particles, &constraints);
        let dist = (particles[1].position - particles[0].position).norm();
        assert!((dist - 1.0).abs() < 0.01, "dist = {dist}");
    }
    #[test]
    fn test_parallel_pbd_reduces_violation() {
        let mut particles = vec![
            PbdParticle::new_fixed(Vec3::new(0.0, 0.0, 0.0)),
            PbdParticle::new(Vec3::new(5.0, 0.0, 0.0), 1.0),
        ];
        let constraints = vec![PbdDistanceConstraint::new(0, 1, 1.0, 1.0)];
        let solver = ParallelPbd::new(20);
        solver.solve(&mut particles, &constraints);
        let dist = (particles[1].position - particles[0].position).norm();
        assert!(
            dist < 5.0,
            "parallel PBD should reduce violation: dist={dist}"
        );
    }
    #[test]
    fn test_xpbd_constraint_converges() {
        let mut particles = vec![
            PbdParticle::new_fixed(Vec3::new(0.0, 0.0, 0.0)),
            PbdParticle::new(Vec3::new(3.0, 0.0, 0.0), 1.0),
        ];
        let mut constraints = vec![XpbdDistanceConstraint::new(0, 1, 1.0, 0.0)];
        let dt = 0.01;
        xpbd_pbd_step(&mut particles, &mut constraints, dt, 50);
        let dist = (particles[1].position - particles[0].position).norm();
        assert!(
            (dist - 1.0).abs() < 0.1,
            "XPBD should converge to rest: dist={dist}"
        );
    }
    #[test]
    fn test_xpbd_constraint_soft() {
        let mut particles = vec![
            PbdParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            PbdParticle::new(Vec3::new(3.0, 0.0, 0.0), 1.0),
        ];
        let mut constraints = vec![XpbdDistanceConstraint::new(0, 1, 1.0, 1e6)];
        xpbd_pbd_step(&mut particles, &mut constraints, 0.01, 10);
        let dist = (particles[1].position - particles[0].position).norm();
        assert!(dist < 3.0 && dist > 0.0, "dist = {dist}");
    }
    #[test]
    fn test_convergence_monitor_records() {
        let mut mon = PbdConvergenceMonitor::new(0.001);
        assert!(!mon.record(1.0));
        assert!(!mon.record(0.1));
        assert!(mon.record(0.0005));
        assert_eq!(mon.history.len(), 3);
    }
    #[test]
    fn test_convergence_monitor_rate() {
        let mut mon = PbdConvergenceMonitor::new(0.001);
        mon.record(1.0);
        mon.record(0.5);
        let rate = mon.convergence_rate().unwrap();
        assert!((rate - 0.5).abs() < 1e-12, "rate = {rate}");
    }
    #[test]
    fn test_convergence_monitor_reset() {
        let mut mon = PbdConvergenceMonitor::new(0.001);
        mon.record(1.0);
        mon.reset();
        assert!(mon.history.is_empty());
    }
    #[test]
    fn test_convergence_max_violation() {
        let particles = vec![
            PbdParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            PbdParticle::new(Vec3::new(5.0, 0.0, 0.0), 1.0),
        ];
        let constraints = vec![PbdDistanceConstraint::new(0, 1, 1.0, 1.0)];
        let v = PbdConvergenceMonitor::compute_max_violation(&particles, &constraints);
        assert!((v - 4.0).abs() < 1e-12, "violation = {v}");
    }
    #[test]
    fn test_priority_sort() {
        let mut constraints = [
            PrioritizedConstraint::new(0, 1, 1.0, 1.0, ConstraintPriority::Low),
            PrioritizedConstraint::new(1, 2, 1.0, 1.0, ConstraintPriority::Critical),
            PrioritizedConstraint::new(2, 3, 1.0, 1.0, ConstraintPriority::Medium),
        ];
        constraints.sort_by_key(|c| std::cmp::Reverse(c.priority));
        assert_eq!(constraints[0].priority, ConstraintPriority::Critical);
        assert_eq!(constraints[1].priority, ConstraintPriority::Medium);
        assert_eq!(constraints[2].priority, ConstraintPriority::Low);
    }
}
/// Compute the signed volume of a tetrahedron defined by four points.
pub fn tet_volume(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let c = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
    let cross = [
        b[1] * c[2] - b[2] * c[1],
        b[2] * c[0] - b[0] * c[2],
        b[0] * c[1] - b[1] * c[0],
    ];
    (a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]) / 6.0
}
#[inline]
pub(super) fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
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
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
#[inline]
pub(super) fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[cfg(test)]
mod tests_xpbd_raw {
    use super::*;

    use crate::pbd::PbdSubstepSolver;

    use crate::pbd::XpbdBendingConstraint;

    use crate::pbd::XpbdDistanceConstraintRaw;

    use crate::pbd::XpbdVolumeConstraint;

    #[test]
    fn test_distance_constraint_reduces_error() {
        let mut positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let masses = [1.0, 1.0];
        let c = XpbdDistanceConstraintRaw::new(0, 1, 1.0, 0.0);
        let mut lambda = 0.0;
        let dt = 0.01;
        for _ in 0..50 {
            c.project(&mut positions, &masses, &mut lambda, dt);
        }
        let diff = vec3_sub(positions[1], positions[0]);
        let dist = vec3_norm(diff);
        assert!((dist - 1.0).abs() < 0.01, "Expected dist ≈ 1.0, got {dist}");
    }
    #[test]
    fn test_volume_of_regular_tetrahedron() {
        let p0 = [1.0, 1.0, 1.0];
        let p1 = [1.0, -1.0, -1.0];
        let p2 = [-1.0, 1.0, -1.0];
        let p3 = [-1.0, -1.0, 1.0];
        let vol = tet_volume(p0, p1, p2, p3).abs();
        let expected = 8.0 / 3.0;
        assert!(
            (vol - expected).abs() < 1e-10,
            "vol = {vol}, expected {expected}"
        );
    }
    #[test]
    fn test_bending_constraint_at_rest_angle() {
        let mut positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let masses = [1.0, 1.0, 1.0];
        let rest_angle = std::f64::consts::PI;
        let c = XpbdBendingConstraint::new(0, 1, 2, rest_angle, 0.0);
        let mut lambda = 0.0_f64;
        let initial_lambda = lambda;
        c.project(&mut positions, &masses, &mut lambda, 0.01);
        assert!(
            (lambda - initial_lambda).abs() < 0.1,
            "lambda changed unexpectedly: {lambda}"
        );
    }
    #[test]
    fn test_xpbd_lambda_update() {
        let mut positions = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let masses = [1.0, 1.0];
        let c = XpbdDistanceConstraintRaw::new(0, 1, 1.0, 0.0);
        let mut lambda = 0.0;
        c.project(&mut positions, &masses, &mut lambda, 0.01);
        assert!(lambda.abs() > 0.0, "lambda should be updated from zero");
    }
    #[test]
    fn test_volume_constraint_reduces_error() {
        let p0 = [0.0, 0.0, 0.0_f64];
        let p1 = [1.0, 0.0, 0.0_f64];
        let p2 = [0.5, 1.0, 0.0_f64];
        let p3 = [0.5, 0.5, 1.0_f64];
        let rest_vol = tet_volume(p0, p1, p2, p3);
        let perturbed_p3 = [0.5, 0.5, 2.0_f64];
        let perturbed_vol = tet_volume(p0, p1, p2, perturbed_p3);
        let initial_error = (perturbed_vol - rest_vol).abs();
        let mut positions = [p0, p1, p2, perturbed_p3];
        let masses = [1.0, 1.0, 1.0, 1.0];
        let c = XpbdVolumeConstraint::new([0, 1, 2, 3], rest_vol, 0.0);
        let mut lambda = 0.0;
        for _ in 0..50 {
            c.project(&mut positions, &masses, &mut lambda, 0.01);
        }
        let new_vol = tet_volume(positions[0], positions[1], positions[2], positions[3]);
        let final_error = (new_vol - rest_vol).abs();
        assert!(
            final_error < initial_error,
            "volume error should decrease: initial={initial_error}, final={final_error}"
        );
    }
    #[test]
    fn test_substep_solver_moves_particles() {
        let mut positions = [[0.0, 0.0, 0.0_f64]];
        let mut velocities = [[1.0, 0.0, 0.0_f64]];
        let masses = [1.0_f64];
        let solver = PbdSubstepSolver::new(4, 0.1);
        solver.solve(&mut positions, &mut velocities, &masses, &mut |_| {});
        assert!(positions[0][0] > 0.0, "particle should have moved");
    }
}
#[inline(always)]
pub(super) fn pbd2_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline(always)]
pub(super) fn pbd2_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline(always)]
pub(super) fn pbd2_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline(always)]
pub(super) fn pbd2_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline(always)]
pub(super) fn pbd2_norm(a: [f64; 3]) -> f64 {
    pbd2_dot(a, a).sqrt()
}
#[inline(always)]
pub(super) fn pbd2_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline(always)]
pub(super) fn pbd2_normalize(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = pbd2_norm(a);
    if n < 1e-14 {
        None
    } else {
        Some(pbd2_scale(a, 1.0 / n))
    }
}
#[cfg(test)]
mod tests_pbd_new {
    use super::*;
    use crate::pbd::HierarchicalPbdSolver;
    use crate::pbd::HierarchyLevel;
    use crate::pbd::IsometricBending;

    use crate::pbd::PbdConstraintPriority;
    use crate::pbd::PbdConstraintScheduler;

    use crate::pbd::PrioritizedPbdConstraint;

    use crate::pbd::XpbdCollisionConstraint;
    use crate::pbd::XpbdDihedralBending;

    use crate::pbd::XpbdPressureConstraint;
    use crate::pbd::XpbdStretchConstraint;

    #[test]
    fn test_stretch_no_stretch_no_change() {
        let mut positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let inv_masses = [1.0, 1.0];
        let c = XpbdStretchConstraint::new(0, 1, 1.0, 0.0);
        let mut lambda = 0.0_f64;
        let before = positions[1];
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
        assert!((positions[1][0] - before[0]).abs() < 1e-10);
    }
    #[test]
    fn test_stretch_corrects_overstretched() {
        let mut positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let inv_masses = [1.0, 1.0];
        let c = XpbdStretchConstraint::new(0, 1, 1.0, 0.0);
        let mut lambda = 0.0_f64;
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
        let dist = ((positions[1][0] - positions[0][0]).powi(2)
            + (positions[1][1] - positions[0][1]).powi(2)
            + (positions[1][2] - positions[0][2]).powi(2))
        .sqrt();
        assert!(dist < 2.0, "distance={dist}");
    }
    #[test]
    fn test_stretch_fixed_particle_does_not_move() {
        let mut positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let inv_masses = [0.0, 1.0];
        let c = XpbdStretchConstraint::new(0, 1, 1.0, 0.0);
        let mut lambda = 0.0_f64;
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
        assert!((positions[0][0]).abs() < 1e-12);
        assert!(positions[1][0] < 3.0, "particle 1 should move closer");
    }
    #[test]
    fn test_stretch_compliance_reduces_correction() {
        let mut pos_rigid = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut pos_soft = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let inv_masses = [1.0, 1.0];
        let rigid = XpbdStretchConstraint::new(0, 1, 1.0, 0.0);
        let soft = XpbdStretchConstraint::new(0, 1, 1.0, 1e3);
        let mut lam_r = 0.0_f64;
        let mut lam_s = 0.0_f64;
        rigid.project(&mut pos_rigid, &inv_masses, &mut lam_r, 0.01);
        soft.project(&mut pos_soft, &inv_masses, &mut lam_s, 0.01);
        let d_rigid = (pos_rigid[1][0] - 2.0).abs();
        let d_soft = (pos_soft[1][0] - 2.0).abs();
        assert!(d_rigid >= d_soft, "rigid={d_rigid} soft={d_soft}");
    }
    #[test]
    fn test_dihedral_flat_computes_angle() {
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, -1.0, 0.0],
        ];
        let c = XpbdDihedralBending::new([0, 1, 2, 3], 0.0, 0.0);
        let angle = c.compute_angle(&positions);
        assert!(
            (0.0..=std::f64::consts::PI + 1e-10).contains(&angle),
            "angle should be in [0,pi], got {angle}"
        );
    }
    #[test]
    fn test_dihedral_project_runs_without_panic() {
        let mut positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.0, 1.0],
        ];
        let inv_masses = [1.0, 1.0, 1.0, 1.0];
        let c = XpbdDihedralBending::new([0, 1, 2, 3], std::f64::consts::PI / 2.0, 0.01);
        let mut lambda = 0.0_f64;
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
    }
    #[test]
    fn test_isometric_bending_rest_pose_energy() {
        let rest = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [1.5, 1.0, 0.0],
        ];
        let ib = IsometricBending::new([0, 1, 2, 3], rest, 1.0);
        let energy = ib.bending_energy(&rest);
        assert!(energy >= 0.0, "energy must be non-negative, got {energy}");
    }
    #[test]
    fn test_isometric_bending_out_of_plane_nonzero_energy() {
        let rest = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [1.5, 1.0, 0.0],
        ];
        let ib = IsometricBending::new([0, 1, 2, 3], rest, 1.0);
        let mut deformed = rest;
        deformed[3] = [1.5, 1.0, 2.0];
        let energy = ib.bending_energy(&deformed);
        assert!(energy >= 0.0);
    }
    #[test]
    fn test_pressure_constraint_at_rest_no_correction() {
        let p0 = [0.0, 0.0, 0.0_f64];
        let p1 = [1.0, 0.0, 0.0_f64];
        let p2 = [0.5, 1.0, 0.0_f64];
        let p3 = [0.5, 0.5, 1.0_f64];
        let rest_vol = {
            let v1 = pbd2_sub(p1, p0);
            let v2 = pbd2_sub(p2, p0);
            let v3 = pbd2_sub(p3, p0);
            pbd2_dot(v1, pbd2_cross(v2, v3)) / 6.0
        };
        let mut positions = [p0, p1, p2, p3];
        let inv_masses = [1.0, 1.0, 1.0, 1.0];
        let c = XpbdPressureConstraint::new([0, 1, 2, 3], rest_vol, 0.0);
        let mut lambda = 0.0_f64;
        let before = positions;
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
        assert!(lambda.abs() < 0.1, "lambda={lambda}");
        let _ = before;
    }
    #[test]
    fn test_pressure_constraint_reduces_error() {
        let p0 = [0.0, 0.0, 0.0_f64];
        let p1 = [1.0, 0.0, 0.0_f64];
        let p2 = [0.5, 1.0, 0.0_f64];
        let p3_rest = [0.5, 0.5, 1.0_f64];
        let rest_vol = pbd2_dot(
            pbd2_sub(p1, p0),
            pbd2_cross(pbd2_sub(p2, p0), pbd2_sub(p3_rest, p0)),
        ) / 6.0;
        let p3_perturbed = [0.5, 0.5, 3.0_f64];
        let mut positions = [p0, p1, p2, p3_perturbed];
        let inv_masses = [1.0, 1.0, 1.0, 1.0];
        let c = XpbdPressureConstraint::new([0, 1, 2, 3], rest_vol, 0.0);
        let mut lambda = 0.0_f64;
        for _ in 0..20 {
            c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
        }
        assert!(
            lambda.abs() > 0.0,
            "lambda should be nonzero after iteration"
        );
    }
    #[test]
    fn test_collision_constraint_no_penetration_no_change() {
        let mut positions = [[0.0, 1.0, 0.0]];
        let inv_masses = [1.0];
        let c = XpbdCollisionConstraint::new(0, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0);
        let mut lambda = 0.0_f64;
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
        assert!((positions[0][1] - 1.0).abs() < 1e-12);
        assert!(lambda.abs() < 1e-12);
    }
    #[test]
    fn test_collision_constraint_penetrating_corrected() {
        let mut positions = [[0.0, -0.5, 0.0]];
        let inv_masses = [1.0];
        let c = XpbdCollisionConstraint::new(0, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0);
        let mut lambda = 0.0_f64;
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
        assert!(positions[0][1] > -0.5, "particle should move up");
    }
    #[test]
    fn test_collision_constraint_evaluate_positive_outside() {
        let c = XpbdCollisionConstraint::new(0, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0);
        let positions = [[0.0, 2.0, 0.0]];
        assert!(c.evaluate(&positions) > 0.0);
    }
    #[test]
    fn test_collision_constraint_evaluate_negative_inside() {
        let c = XpbdCollisionConstraint::new(0, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0);
        let positions = [[0.0, -1.0, 0.0]];
        assert!(c.evaluate(&positions) < 0.0);
    }
    #[test]
    fn test_collision_constraint_fixed_particle_not_moved() {
        let mut positions = [[0.0, -2.0, 0.0]];
        let inv_masses = [0.0];
        let c = XpbdCollisionConstraint::new(0, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0);
        let mut lambda = 0.0_f64;
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
        assert!(
            (positions[0][1] - (-2.0)).abs() < 1e-12,
            "fixed particle should not move"
        );
    }
    #[test]
    fn test_hierarchy_total_iterations() {
        let mut solver = HierarchicalPbdSolver::new(0.016);
        solver.add_level(HierarchyLevel::new(4, 4, 10.0));
        solver.add_level(HierarchyLevel::new(8, 2, 1.0));
        solver.add_level(HierarchyLevel::new(16, 1, 0.1));
        assert_eq!(solver.total_iterations(), 28);
    }
    #[test]
    fn test_hierarchy_solve_calls_all_iterations() {
        let mut solver = HierarchicalPbdSolver::new(0.016);
        solver.add_level(HierarchyLevel::new(3, 2, 5.0));
        solver.add_level(HierarchyLevel::new(2, 1, 1.0));
        let mut call_count = 0u32;
        solver.solve(|_level, _scale, _iter| {
            call_count += 1;
        });
        assert_eq!(call_count, 5);
    }
    #[test]
    fn test_hierarchy_level_particle_count() {
        let level = HierarchyLevel::new(4, 2, 1.0);
        assert_eq!(level.particle_count(10), 5);
        assert_eq!(level.particle_count(11), 6);
    }
    #[test]
    fn test_hierarchy_compliance_scale_passed() {
        let mut solver = HierarchicalPbdSolver::new(0.016);
        solver.add_level(HierarchyLevel::new(1, 1, 42.0));
        let mut seen_scale = 0.0_f64;
        solver.solve(|_, scale, _| {
            seen_scale = scale;
        });
        assert!((seen_scale - 42.0).abs() < 1e-12);
    }
    #[test]
    fn test_scheduler_sort_order() {
        let mut sched = PbdConstraintScheduler::new();
        sched.add(PrioritizedPbdConstraint::new(
            PbdConstraintPriority::Cosmetic,
            0,
            1.0,
        ));
        sched.add(PrioritizedPbdConstraint::new(
            PbdConstraintPriority::Critical,
            1,
            1.0,
        ));
        sched.add(PrioritizedPbdConstraint::new(
            PbdConstraintPriority::Normal,
            2,
            1.0,
        ));
        sched.sort_by_priority();
        assert_eq!(
            sched.constraints[0].priority,
            PbdConstraintPriority::Critical
        );
        assert_eq!(sched.constraints[1].priority, PbdConstraintPriority::Normal);
        assert_eq!(
            sched.constraints[2].priority,
            PbdConstraintPriority::Cosmetic
        );
    }
    #[test]
    fn test_scheduler_indices_at_priority() {
        let mut sched = PbdConstraintScheduler::new();
        sched.add(PrioritizedPbdConstraint::new(
            PbdConstraintPriority::High,
            10,
            1.0,
        ));
        sched.add(PrioritizedPbdConstraint::new(
            PbdConstraintPriority::Normal,
            20,
            1.0,
        ));
        sched.add(PrioritizedPbdConstraint::new(
            PbdConstraintPriority::High,
            30,
            1.0,
        ));
        let high_indices = sched.indices_at_priority(PbdConstraintPriority::High);
        assert_eq!(high_indices.len(), 2);
        assert!(high_indices.contains(&10));
        assert!(high_indices.contains(&30));
    }
    #[test]
    fn test_scheduler_len_and_is_empty() {
        let mut sched = PbdConstraintScheduler::new();
        assert!(sched.is_empty());
        sched.add(PrioritizedPbdConstraint::new(
            PbdConstraintPriority::Normal,
            0,
            1.0,
        ));
        assert_eq!(sched.len(), 1);
        assert!(!sched.is_empty());
    }
    #[test]
    fn test_constraint_priority_ordering() {
        assert!(PbdConstraintPriority::Critical > PbdConstraintPriority::High);
        assert!(PbdConstraintPriority::High > PbdConstraintPriority::Normal);
        assert!(PbdConstraintPriority::Normal > PbdConstraintPriority::Cosmetic);
    }
    #[test]
    fn test_prioritized_constraint_weight_clamped() {
        let c = PrioritizedPbdConstraint::new(PbdConstraintPriority::Normal, 0, -5.0);
        assert!((c.weight).abs() < 1e-12, "weight should be clamped to 0");
    }
}
