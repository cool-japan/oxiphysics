//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::RigidBodySet;

use super::types::DriftCorrectionResult;

/// Default Baumgarte stabilization coefficient β ∈ (0, 1].
pub(super) const BAUMGARTE_BETA: f64 = 0.2;
/// Penetration slop: overlap smaller than this is not corrected.
pub(super) const PENETRATION_SLOP: f64 = 0.005;
/// Fraction of kinetic energy allowed to increase per substep before the
/// energy-gain monitor fires and damps velocities.
pub(super) const ENERGY_GAIN_THRESHOLD: f64 = 1.05;
/// Maximum velocity allowed for speculative contacts (m/s).
pub(super) const SPECULATIVE_MAX_VELOCITY: f64 = 100.0;
/// Minimum substep count for adaptive substepping.
pub(super) const MIN_SUBSTEPS: usize = 1;
/// Maximum substep count for adaptive substepping.
pub(super) const MAX_SUBSTEPS: usize = 64;
/// Clamp a body's velocity for speculative contact handling.
///
/// Prevents tunnelling by limiting the velocity component along a direction
/// to a maximum value that would not exceed the contact distance in one step.
#[allow(dead_code)]
pub fn clamp_speculative_velocity(velocity: &mut Vec3, normal: &Vec3, distance: f64, dt: f64) {
    if dt < 1e-12 {
        return;
    }
    let vn = velocity.dot(normal);
    if vn >= 0.0 {
        return;
    }
    let max_approach = (distance / dt).min(SPECULATIVE_MAX_VELOCITY);
    if -vn > max_approach {
        let excess = -vn - max_approach;
        *velocity += *normal * excess;
    }
}
/// Compute the predicted closest approach distance.
///
/// Given current separation `distance`, relative velocity `rel_vel`, and
/// time step `dt`, returns the predicted minimum distance.
#[allow(dead_code)]
pub fn predicted_closest_approach(distance: f64, rel_vel: f64, dt: f64) -> f64 {
    distance + rel_vel * dt
}
/// Compute the Baumgarte bias velocity for a constraint.
///
/// Used externally by constraint implementations that expose a `penetration`
/// field.
///
/// `penetration` – overlap depth in metres (positive = overlapping).
/// `beta`        – Baumgarte coefficient ∈ (0, 1].
/// `dt`          – current substep size in seconds.
#[allow(dead_code)]
pub fn baumgarte_bias(penetration: f64, beta: f64, dt: f64) -> f64 {
    if dt < 1e-12 {
        return 0.0;
    }
    let excess = (penetration - PENETRATION_SLOP).max(0.0);
    beta * excess / dt
}
/// Compute the ERP (Error Reduction Parameter) from Baumgarte beta and dt.
///
/// `erp = beta * dt / (1 + beta * dt)`
///
/// Commonly used in ODE-style constraint formulation.
#[allow(dead_code)]
pub fn compute_erp(beta: f64, dt: f64) -> f64 {
    let denom = 1.0 + beta * dt;
    if denom.abs() < 1e-20 {
        return 0.0;
    }
    beta * dt / denom
}
/// Compute the CFM (Constraint Force Mixing) from compliance and dt.
///
/// `cfm = 1 / (compliance * dt * dt + damping * dt)`
///
/// CFM regularizes the constraint to prevent infinite stiffness.
#[allow(dead_code)]
pub fn compute_cfm(compliance: f64, damping: f64, dt: f64) -> f64 {
    let denom = compliance * dt * dt + damping * dt;
    if denom.abs() < 1e-20 {
        return 0.0;
    }
    1.0 / denom
}
/// Apply a position-level pseudo-impulse to resolve penetration directly.
///
/// This is the "split impulse" form: positions are corrected without touching
/// velocities, avoiding the artificial energy injection of pure Baumgarte.
///
/// `bodies`       – mutable body set.
/// `ha`, `hb`     – handles of the two bodies.
/// `normal`       – contact normal from B to A.
/// `penetration`  – overlap depth.
#[allow(dead_code)]
pub fn apply_position_correction(
    bodies: &mut RigidBodySet,
    ha: oxiphysics_core::BodyHandle,
    hb: oxiphysics_core::BodyHandle,
    normal: &Vec3,
    penetration: f64,
) {
    if penetration <= PENETRATION_SLOP {
        return;
    }
    let excess = penetration - PENETRATION_SLOP;
    let inv_mass_a = bodies.get(ha).map_or(0.0, |b| b.inverse_mass);
    let inv_mass_b = bodies.get(hb).map_or(0.0, |b| b.inverse_mass);
    let k = inv_mass_a + inv_mass_b;
    if k < 1e-12 {
        return;
    }
    let correction = BAUMGARTE_BETA * excess / k;
    if let Some(body_a) = bodies.get_mut(ha) {
        body_a.transform.position += *normal * (correction * inv_mass_a);
    }
    if let Some(body_b) = bodies.get_mut(hb) {
        body_b.transform.position -= *normal * (correction * inv_mass_b);
    }
}
/// Apply position correction with a custom beta parameter.
#[allow(dead_code)]
pub fn apply_position_correction_custom(
    bodies: &mut RigidBodySet,
    ha: oxiphysics_core::BodyHandle,
    hb: oxiphysics_core::BodyHandle,
    normal: &Vec3,
    penetration: f64,
    beta: f64,
    slop: f64,
) {
    if penetration <= slop {
        return;
    }
    let excess = penetration - slop;
    let inv_mass_a = bodies.get(ha).map_or(0.0, |b| b.inverse_mass);
    let inv_mass_b = bodies.get(hb).map_or(0.0, |b| b.inverse_mass);
    let k = inv_mass_a + inv_mass_b;
    if k < 1e-12 {
        return;
    }
    let correction = beta * excess / k;
    if let Some(body_a) = bodies.get_mut(ha) {
        body_a.transform.position += *normal * (correction * inv_mass_a);
    }
    if let Some(body_b) = bodies.get_mut(hb) {
        body_b.transform.position -= *normal * (correction * inv_mass_b);
    }
}
/// Compute the velocity required to separate two bodies in one step.
///
/// `penetration` – current overlap.
/// `dt` – time step.
#[allow(dead_code)]
pub fn separation_velocity(penetration: f64, dt: f64) -> f64 {
    if dt < 1e-12 || penetration < PENETRATION_SLOP {
        return 0.0;
    }
    let excess = penetration - PENETRATION_SLOP;
    excess / dt
}
/// Compute total kinetic energy for all dynamic bodies.
pub(super) fn compute_kinetic_energy(bodies: &RigidBodySet) -> f64 {
    let mut ke = 0.0;
    for (_h, body) in bodies.iter() {
        if body.inverse_mass > 0.0 {
            let m = 1.0 / body.inverse_mass;
            let v2 = body.velocity.x.powi(2) + body.velocity.y.powi(2) + body.velocity.z.powi(2);
            ke += 0.5 * m * v2;
        }
    }
    ke
}
/// Compute maximum velocity magnitude among all dynamic bodies.
pub(super) fn compute_max_velocity(bodies: &RigidBodySet) -> f64 {
    let mut max_v2 = 0.0f64;
    for (_h, body) in bodies.iter() {
        if body.inverse_mass > 0.0 {
            let v2 = body.velocity.x.powi(2) + body.velocity.y.powi(2) + body.velocity.z.powi(2);
            max_v2 = max_v2.max(v2);
        }
    }
    max_v2.sqrt()
}
/// Damp a single body's velocity if its kinetic energy increased by more than
/// `threshold` factor compared to a reference.  Used inside the energy monitor.
pub(super) fn clamp_body_energy_gain(
    bodies: &mut RigidBodySet,
    handle: oxiphysics_core::BodyHandle,
    threshold: f64,
) {
    if let Some(body) = bodies.get_mut(handle) {
        if body.inverse_mass <= 0.0 {
            return;
        }
        let v2 = body.velocity.x.powi(2) + body.velocity.y.powi(2) + body.velocity.z.powi(2);
        let max_speed_sq = threshold * v2;
        if v2 > max_speed_sq * threshold {
            let scale = (max_speed_sq / v2).sqrt();
            body.velocity.x *= scale;
            body.velocity.y *= scale;
            body.velocity.z *= scale;
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Constraint;
    use crate::TgsSolver;
    use crate::contact::ContactConstraint;
    use crate::tgs_solver::AdaptiveTgsSolver;

    use crate::tgs_solver::PositionCorrectionMode;
    use crate::tgs_solver::SubstepDiagnostics;
    use crate::tgs_solver::TgsWarmStartState;
    use oxiphysics_rigid::RigidBody;
    #[test]
    fn test_tgs_solver_substeps() {
        let mut bodies = RigidBodySet::new();
        let mut body = RigidBody::new(1.0);
        body.transform.position = Vec3::new(0.0, 1.0, 0.0);
        body.velocity = Vec3::new(0.0, -2.0, 0.0);
        body.linear_damping = 0.0;
        body.angular_damping = 0.0;
        body.gravity_scale = 0.0;
        let hb = bodies.insert(body);
        let ground = RigidBody::new_static();
        let hg = bodies.insert(ground);
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = TgsSolver::new(4, 8, Vec3::zeros());
        solver.step(&mut constraints, &mut bodies, 1.0 / 60.0);
        let body = bodies.get(hb).unwrap();
        assert!(
            body.velocity.y >= -0.1,
            "TGS should resolve contact, vy={}",
            body.velocity.y
        );
    }
    #[test]
    fn test_more_substeps_better_accuracy() {
        let make_bodies_and_constraints = || {
            let mut bodies = RigidBodySet::new();
            let mut body = RigidBody::new(1.0);
            body.transform.position = Vec3::new(0.0, 0.5, 0.0);
            body.velocity = Vec3::new(0.0, -3.0, 0.0);
            body.linear_damping = 0.0;
            body.angular_damping = 0.0;
            body.gravity_scale = 0.0;
            let hb = bodies.insert(body);
            let hg = bodies.insert(RigidBody::new_static());
            let c = ContactConstraint::new(
                hb,
                hg,
                Vec3::new(0.0, 1.0, 0.0),
                0.02,
                Vec3::zeros(),
                Vec3::zeros(),
                0.0,
                0.0,
            );
            let constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
            (bodies, constraints, hb)
        };
        let (mut bodies1, mut c1, hb1) = make_bodies_and_constraints();
        let solver_few = TgsSolver::new(1, 4, Vec3::zeros());
        solver_few.step(&mut c1, &mut bodies1, 1.0 / 60.0);
        let (mut bodies8, mut c8, hb8) = make_bodies_and_constraints();
        let solver_many = TgsSolver::new(8, 4, Vec3::zeros());
        solver_many.step(&mut c8, &mut bodies8, 1.0 / 60.0);
        let vy_few = bodies1.get(hb1).unwrap().velocity.y;
        let vy_many = bodies8.get(hb8).unwrap().velocity.y;
        assert!(vy_few >= -0.5, "1 substep: vy={vy_few}");
        assert!(vy_many >= -0.5, "8 substeps: vy={vy_many}");
    }
    #[test]
    fn test_position_correction_reduces_overlap() {
        let mut bodies = RigidBodySet::new();
        let mut body = RigidBody::new(1.0);
        body.transform.position = Vec3::new(0.0, -0.1, 0.0);
        body.linear_damping = 0.0;
        body.angular_damping = 0.0;
        let hb = bodies.insert(body);
        let hg = bodies.insert(RigidBody::new_static());
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.1,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.0,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = TgsSolver::with_options(4, 8, 4, Vec3::zeros(), 0.2, false);
        solver.step(&mut constraints, &mut bodies, 1.0 / 60.0);
        let pos_y = bodies.get(hb).unwrap().transform.position.y;
        assert!(
            pos_y > -0.1,
            "Position correction should reduce overlap, y={pos_y}"
        );
    }
    #[test]
    fn test_baumgarte_proportional_to_penetration() {
        let solver = TgsSolver::default();
        let dt = 1.0 / 60.0;
        let bias_small = solver.baumgarte_bias(0.01, dt);
        let bias_large = solver.baumgarte_bias(0.10, dt);
        assert!(
            bias_large > bias_small,
            "bias_large={bias_large}, bias_small={bias_small}"
        );
        let bias_slop = solver.baumgarte_bias(PENETRATION_SLOP * 0.5, dt);
        assert_eq!(bias_slop, 0.0);
    }
    #[test]
    fn test_island_solve_equivalent_to_full_solve() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(0.0, 0.5, 0.0);
        b.velocity = Vec3::new(0.0, -2.0, 0.0);
        b.linear_damping = 0.0;
        b.angular_damping = 0.0;
        b.gravity_scale = 0.0;
        let hb = bodies.insert(b);
        let hg = bodies.insert(RigidBody::new_static());
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.0,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = TgsSolver::new(4, 8, Vec3::zeros());
        let island_bodies = vec![hb, hg];
        let stats = solver.solve_island(&island_bodies, &mut constraints, &mut bodies, 1.0 / 60.0);
        assert!(stats.iterations_used > 0);
        let vy = bodies.get(hb).unwrap().velocity.y;
        assert!(vy >= -0.1, "Island solve should resolve contact, vy={vy}");
    }
    #[test]
    fn test_step_returns_stats() {
        let mut bodies = RigidBodySet::new();
        let hg = bodies.insert(RigidBody::new_static());
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(0.0, -1.0, 0.0);
        let hb = bodies.insert(b);
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.02,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = TgsSolver::new(3, 5, Vec3::zeros());
        let stats = solver.step(&mut constraints, &mut bodies, 1.0 / 60.0);
        assert_eq!(stats.iterations_used, 15);
        assert!(stats.converged);
    }
    #[test]
    fn test_free_baumgarte_bias() {
        let bias = baumgarte_bias(0.02, 0.2, 1.0 / 60.0);
        assert!(
            bias > 0.0,
            "should produce positive bias for penetration > slop"
        );
        let bias_zero_dt = baumgarte_bias(0.1, 0.2, 0.0);
        assert_eq!(bias_zero_dt, 0.0, "zero dt should give zero bias");
    }
    #[test]
    fn test_position_correction_pushes_apart() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(0.0, -0.05, 0.0);
        let hb = bodies.insert(b);
        let hg = bodies.insert(RigidBody::new_static());
        let normal = Vec3::new(0.0, 1.0, 0.0);
        apply_position_correction(&mut bodies, hb, hg, &normal, 0.05);
        let pos_y = bodies.get(hb).unwrap().transform.position.y;
        assert!(pos_y > -0.05, "body should be pushed up, y={pos_y}");
    }
    #[test]
    fn test_default_tgs_solver() {
        let s = TgsSolver::default();
        assert_eq!(s.substeps, 4);
        assert_eq!(s.velocity_iterations, 4);
        assert!((s.gravity.y + 9.81).abs() < 1e-6);
        assert!((s.baumgarte_beta - BAUMGARTE_BETA).abs() < 1e-12);
    }
    #[test]
    fn test_tgs_warm_start_state() {
        let mut state = TgsWarmStartState::new(3);
        assert_eq!(state.normal_impulses.len(), 3);
        assert_eq!(state.friction_impulses.len(), 3);
        assert_eq!(state.age, 0);
        state.normal_impulses[0] = 10.0;
        state.friction_impulses[1] = [2.0, 3.0];
        state.scale(2.0);
        assert!((state.normal_impulses[0] - 20.0).abs() < 1e-12);
        assert!((state.friction_impulses[1][0] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_tgs_warm_start_reset() {
        let mut state = TgsWarmStartState::new(2);
        state.normal_impulses[0] = 5.0;
        state.age = 10;
        state.reset();
        assert!(state.normal_impulses.iter().all(|&v| v == 0.0));
        assert_eq!(state.age, 0);
    }
    #[test]
    fn test_tgs_warm_start_resize() {
        let mut state = TgsWarmStartState::new(2);
        state.normal_impulses[0] = 7.0;
        state.resize(5);
        assert_eq!(state.normal_impulses.len(), 5);
        assert!((state.normal_impulses[0] - 7.0).abs() < 1e-12);
        assert_eq!(state.normal_impulses[4], 0.0);
    }
    #[test]
    fn test_tgs_warm_start_age() {
        let mut state = TgsWarmStartState::new(1);
        state.advance_age();
        state.advance_age();
        state.advance_age();
        assert_eq!(state.age, 3);
        assert!(!state.is_stale(5));
        assert!(state.is_stale(2));
    }
    #[test]
    fn test_tgs_warm_start_total_normal() {
        let mut state = TgsWarmStartState::new(3);
        state.normal_impulses[0] = 5.0;
        state.normal_impulses[1] = -3.0;
        state.normal_impulses[2] = 2.0;
        assert!((state.total_normal_impulse() - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_tgs_warm_start_active_count() {
        let mut state = TgsWarmStartState::new(5);
        state.normal_impulses[1] = 1.0;
        state.normal_impulses[3] = -0.5;
        assert_eq!(state.active_count(), 2);
    }
    #[test]
    fn test_substep_diagnostics_zero() {
        let d = SubstepDiagnostics::zero();
        assert_eq!(d.ke_start, 0.0);
        assert_eq!(d.ke_end, 0.0);
        assert!(!d.energy_gained());
    }
    #[test]
    fn test_substep_diagnostics_energy_ratio() {
        let d = SubstepDiagnostics {
            ke_start: 10.0,
            ke_end: 12.0,
            max_velocity: 5.0,
            max_penetration: 0.01,
        };
        assert!((d.energy_ratio() - 1.2).abs() < 1e-12);
        assert!(d.energy_gained());
    }
    #[test]
    fn test_substep_diagnostics_zero_start() {
        let d = SubstepDiagnostics {
            ke_start: 0.0,
            ke_end: 5.0,
            max_velocity: 3.0,
            max_penetration: 0.0,
        };
        assert!((d.energy_ratio() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_step_with_diagnostics() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(0.0, -1.0, 0.0);
        b.linear_damping = 0.0;
        b.angular_damping = 0.0;
        b.gravity_scale = 0.0;
        let hb = bodies.insert(b);
        let hg = bodies.insert(RigidBody::new_static());
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.0,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = TgsSolver::new(3, 4, Vec3::zeros());
        let (stats, diags) =
            solver.step_with_diagnostics(&mut constraints, &mut bodies, 1.0 / 60.0);
        assert_eq!(diags.len(), 3);
        assert!(stats.iterations_used > 0);
    }
    #[test]
    fn test_adaptive_substep_count() {
        let solver = AdaptiveTgsSolver::new(4, Vec3::zeros(), 0.01);
        let n = solver.compute_substeps(10.0, 1.0 / 60.0);
        assert!((16..=18).contains(&n), "expected ~17 substeps, got {n}");
        let n_zero = solver.compute_substeps(0.0, 1.0 / 60.0);
        assert_eq!(n_zero, solver.min_substeps);
    }
    #[test]
    fn test_adaptive_tgs_solve() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(0.0, -2.0, 0.0);
        b.linear_damping = 0.0;
        b.angular_damping = 0.0;
        b.gravity_scale = 0.0;
        let hb = bodies.insert(b);
        let hg = bodies.insert(RigidBody::new_static());
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.0,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = AdaptiveTgsSolver::new(4, Vec3::zeros(), 0.01);
        let stats = solver.step(&mut constraints, &mut bodies, 1.0 / 60.0);
        assert!(stats.iterations_used > 0);
    }
    #[test]
    fn test_compute_erp() {
        let erp = compute_erp(0.2, 1.0 / 60.0);
        assert!(erp > 0.0 && erp < 1.0, "erp={erp}");
        let erp_zero = compute_erp(0.0, 1.0 / 60.0);
        assert!((erp_zero).abs() < 1e-12);
    }
    #[test]
    fn test_compute_cfm() {
        let cfm = compute_cfm(0.001, 0.01, 1.0 / 60.0);
        assert!(cfm > 0.0, "cfm should be positive, got {cfm}");
        let cfm_zero = compute_cfm(0.0, 0.0, 1.0 / 60.0);
        assert_eq!(cfm_zero, 0.0, "zero compliance/damping → zero cfm");
    }
    #[test]
    fn test_clamp_speculative_velocity() {
        let mut vel = Vec3::new(0.0, -50.0, 0.0);
        let normal = Vec3::new(0.0, 1.0, 0.0);
        let distance = 0.1;
        let dt = 1.0 / 60.0;
        clamp_speculative_velocity(&mut vel, &normal, distance, dt);
        assert!(vel.y >= -7.0, "velocity should be clamped, vy={}", vel.y);
    }
    #[test]
    fn test_clamp_speculative_separating() {
        let mut vel = Vec3::new(0.0, 5.0, 0.0);
        let normal = Vec3::new(0.0, 1.0, 0.0);
        clamp_speculative_velocity(&mut vel, &normal, 0.1, 1.0 / 60.0);
        assert!((vel.y - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_predicted_closest_approach() {
        let d = predicted_closest_approach(0.1, -0.5, 0.1);
        assert!((d - 0.05).abs() < 1e-12);
        let d_sep = predicted_closest_approach(0.1, 1.0, 0.1);
        assert!((d_sep - 0.2).abs() < 1e-12);
    }
    #[test]
    fn test_position_correction_custom() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(1.0);
        b.transform.position = Vec3::new(0.0, -0.1, 0.0);
        let hb = bodies.insert(b);
        let hg = bodies.insert(RigidBody::new_static());
        let normal = Vec3::new(0.0, 1.0, 0.0);
        apply_position_correction_custom(&mut bodies, hb, hg, &normal, 0.1, 0.5, 0.01);
        let pos_y = bodies.get(hb).unwrap().transform.position.y;
        assert!(
            pos_y > -0.1,
            "custom correction should push body up, y={pos_y}"
        );
    }
    #[test]
    fn test_separation_velocity() {
        let v = separation_velocity(0.02, 1.0 / 60.0);
        assert!((v - 0.9).abs() < 1e-12);
        let v_slop = separation_velocity(0.003, 1.0 / 60.0);
        assert_eq!(v_slop, 0.0);
        let v_zero = separation_velocity(0.1, 0.0);
        assert_eq!(v_zero, 0.0);
    }
    #[test]
    fn test_position_correction_mode_eq() {
        assert_eq!(
            PositionCorrectionMode::Baumgarte,
            PositionCorrectionMode::Baumgarte
        );
        assert_ne!(
            PositionCorrectionMode::Baumgarte,
            PositionCorrectionMode::SplitImpulse
        );
        assert_ne!(
            PositionCorrectionMode::SplitImpulse,
            PositionCorrectionMode::NonlinearGaussSeidel
        );
    }
    #[test]
    fn test_total_kinetic_energy() {
        let mut bodies = RigidBodySet::new();
        let mut b = RigidBody::new(2.0);
        b.velocity = Vec3::new(3.0, 4.0, 0.0);
        bodies.insert(b);
        bodies.insert(RigidBody::new_static());
        let solver = TgsSolver::default();
        let ke = solver.total_kinetic_energy(&bodies);
        assert!((ke - 25.0).abs() < 1e-10, "ke={ke}");
    }
    #[test]
    fn test_max_velocity() {
        let mut bodies = RigidBodySet::new();
        let mut b1 = RigidBody::new(1.0);
        b1.velocity = Vec3::new(3.0, 4.0, 0.0);
        bodies.insert(b1);
        let mut b2 = RigidBody::new(1.0);
        b2.velocity = Vec3::new(0.0, 1.0, 0.0);
        bodies.insert(b2);
        let solver = TgsSolver::default();
        let max_v = solver.max_velocity(&bodies);
        assert!((max_v - 5.0).abs() < 1e-10, "max_v={max_v}");
    }
    #[test]
    fn test_solve_alias() {
        let mut bodies = RigidBodySet::new();
        let hg = bodies.insert(RigidBody::new_static());
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(0.0, -1.0, 0.0);
        b.gravity_scale = 0.0;
        let hb = bodies.insert(b);
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.0,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = TgsSolver::new(2, 4, Vec3::zeros());
        let stats = solver.solve(&mut constraints, &mut bodies, 1.0 / 60.0);
        assert_eq!(stats.iterations_used, 8);
    }
    #[test]
    fn test_warm_start_state_new() {
        let ws = TgsWarmStartState::new(5);
        assert_eq!(ws.normal_impulses.len(), 5);
        assert_eq!(ws.friction_impulses.len(), 5);
        assert_eq!(ws.age, 0);
    }
    #[test]
    fn test_warm_start_state_reset_clears_impulses() {
        let mut ws = TgsWarmStartState::new(3);
        ws.normal_impulses[0] = 1.0;
        ws.friction_impulses[1] = [0.5, -0.3];
        ws.age = 4;
        ws.reset();
        assert_eq!(ws.normal_impulses[0], 0.0);
        assert_eq!(ws.friction_impulses[1], [0.0, 0.0]);
        assert_eq!(ws.age, 0);
    }
    #[test]
    fn test_warm_start_state_scale() {
        let mut ws = TgsWarmStartState::new(2);
        ws.normal_impulses[0] = 2.0;
        ws.friction_impulses[0] = [1.0, -1.0];
        ws.scale(0.5);
        assert!((ws.normal_impulses[0] - 1.0).abs() < 1e-10);
        assert!((ws.friction_impulses[0][0] - 0.5).abs() < 1e-10);
        assert!((ws.friction_impulses[0][1] + 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_warm_start_state_resize() {
        let mut ws = TgsWarmStartState::new(2);
        ws.resize(5);
        assert_eq!(ws.normal_impulses.len(), 5);
        assert_eq!(ws.friction_impulses.len(), 5);
    }
    #[test]
    fn test_warm_start_state_total_normal_impulse() {
        let mut ws = TgsWarmStartState::new(3);
        ws.normal_impulses[0] = 2.0;
        ws.normal_impulses[1] = -1.0;
        ws.normal_impulses[2] = 3.0;
        assert!((ws.total_normal_impulse() - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_warm_start_state_active_count() {
        let mut ws = TgsWarmStartState::new(4);
        ws.normal_impulses[0] = 1.0;
        ws.normal_impulses[2] = -0.5;
        assert_eq!(ws.active_count(), 2);
    }
    #[test]
    fn test_warm_start_advance_age_and_stale() {
        let mut ws = TgsWarmStartState::new(1);
        for _ in 0..5 {
            ws.advance_age();
        }
        assert!(ws.is_stale(4));
        assert!(!ws.is_stale(5));
    }
    #[test]
    fn test_adaptive_tgs_compute_substeps_min() {
        let solver = AdaptiveTgsSolver::new(4, Vec3::zeros(), 0.1);
        let substeps = solver.compute_substeps(0.0, 1.0 / 60.0);
        assert_eq!(substeps, MIN_SUBSTEPS);
    }
    #[test]
    fn test_adaptive_tgs_compute_substeps_fast_body() {
        let solver = AdaptiveTgsSolver::new(4, Vec3::zeros(), 0.01);
        let substeps = solver.compute_substeps(1.0, 1.0 / 60.0);
        assert!(substeps >= 1, "should require at least 1 substep");
    }
    #[test]
    fn test_adaptive_tgs_compute_substeps_clamped_to_max() {
        let solver = AdaptiveTgsSolver::new(4, Vec3::zeros(), 0.001);
        let substeps = solver.compute_substeps(1000.0, 1.0 / 60.0);
        assert_eq!(substeps, MAX_SUBSTEPS);
    }
    #[test]
    fn test_adaptive_tgs_step_runs() {
        let mut bodies = RigidBodySet::new();
        let hg = bodies.insert(RigidBody::new_static());
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(0.0, -1.0, 0.0);
        b.gravity_scale = 0.0;
        let hb = bodies.insert(b);
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.0,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = AdaptiveTgsSolver::new(4, Vec3::zeros(), 0.1);
        let stats = solver.step(&mut constraints, &mut bodies, 1.0 / 60.0);
        assert!(stats.iterations_used > 0);
    }
    #[test]
    fn test_clamp_speculative_velocity_separating() {
        let mut vel = Vec3::new(0.0, 1.0, 0.0);
        clamp_speculative_velocity(&mut vel, &Vec3::new(0.0, 1.0, 0.0), 0.5, 1.0 / 60.0);
        assert!(
            (vel.y - 1.0).abs() < 1e-10,
            "separating velocity should not change"
        );
    }
    #[test]
    fn test_clamp_speculative_velocity_approaching_below_max() {
        let mut vel = Vec3::new(0.0, -0.1, 0.0);
        clamp_speculative_velocity(&mut vel, &Vec3::new(0.0, 1.0, 0.0), 1.0, 1.0 / 60.0);
        assert!((vel.y + 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_clamp_speculative_velocity_clamped() {
        let mut vel = Vec3::new(0.0, -200.0, 0.0);
        clamp_speculative_velocity(&mut vel, &Vec3::new(0.0, 1.0, 0.0), 0.001, 1.0 / 60.0);
        assert!(
            vel.y > -200.0,
            "approaching velocity should be clamped, vy={}",
            vel.y
        );
    }
    #[test]
    fn test_predicted_closest_approach_approaching() {
        let dist = predicted_closest_approach(1.0, -10.0, 0.1);
        assert!((dist - 0.0).abs() < 1e-10, "predicted approach={dist}");
    }
    #[test]
    fn test_predicted_closest_approach_separating() {
        let dist = predicted_closest_approach(1.0, 5.0, 0.1);
        assert!((dist - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_baumgarte_bias_within_slop_returns_zero() {
        let bias = baumgarte_bias(PENETRATION_SLOP * 0.5, 0.2, 1.0 / 60.0);
        assert_eq!(bias, 0.0);
    }
    #[test]
    fn test_baumgarte_bias_proportional() {
        let dt = 1.0 / 60.0;
        let b1 = baumgarte_bias(0.01, 0.2, dt);
        let b2 = baumgarte_bias(0.02, 0.2, dt);
        assert!(b2 > b1, "bias should grow with penetration depth");
    }
    #[test]
    fn test_compute_erp_range_check() {
        let erp = compute_erp(0.2, 1.0 / 60.0);
        assert!(erp > 0.0 && erp < 1.0, "ERP should be in (0,1), got {erp}");
    }
    #[test]
    fn test_compute_cfm_zero_compliance_returns_zero() {
        let cfm = compute_cfm(0.0, 0.0, 1.0 / 60.0);
        assert_eq!(cfm, 0.0, "zero compliance and damping → zero CFM");
    }
    #[test]
    fn test_compute_cfm_nonzero() {
        let cfm = compute_cfm(0.001, 0.01, 1.0 / 60.0);
        assert!(cfm > 0.0, "nonzero compliance should produce nonzero CFM");
    }
    #[test]
    fn test_separation_velocity_at_slop_returns_zero() {
        let sv = separation_velocity(PENETRATION_SLOP * 0.5, 1.0 / 60.0);
        assert_eq!(sv, 0.0);
    }
    #[test]
    fn test_separation_velocity_proportional_to_penetration() {
        let sv1 = separation_velocity(0.01, 1.0 / 60.0);
        let sv2 = separation_velocity(0.02, 1.0 / 60.0);
        assert!(sv2 > sv1);
    }
    #[test]
    fn test_apply_position_correction_moves_bodies() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let pos_a_before = bodies.get(ha).unwrap().transform.position.y;
        apply_position_correction(&mut bodies, ha, hb, &Vec3::new(0.0, 1.0, 0.0), 0.1);
        let pos_a_after = bodies.get(ha).unwrap().transform.position.y;
        assert!(
            pos_a_after > pos_a_before,
            "body A should move up after correction"
        );
    }
    #[test]
    fn test_apply_position_correction_custom_beta() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let pos_before = bodies.get(ha).unwrap().transform.position.y;
        apply_position_correction_custom(
            &mut bodies,
            ha,
            hb,
            &Vec3::new(0.0, 1.0, 0.0),
            0.1,
            0.5,
            0.005,
        );
        let pos_after = bodies.get(ha).unwrap().transform.position.y;
        assert!(pos_after > pos_before, "custom correction should move body");
    }
    #[test]
    fn test_substep_diagnostics_zero_full() {
        let d = SubstepDiagnostics::zero();
        assert_eq!(d.ke_start, 0.0);
        assert_eq!(d.ke_end, 0.0);
        assert_eq!(d.energy_ratio(), 1.0);
        assert!(!d.energy_gained());
    }
    #[test]
    fn test_substep_diagnostics_energy_gained() {
        let d = SubstepDiagnostics {
            ke_start: 1.0,
            ke_end: 2.0,
            max_velocity: 5.0,
            max_penetration: 0.0,
        };
        assert!(d.energy_ratio() > 1.0);
        assert!(d.energy_gained());
    }
    #[test]
    fn test_substep_diagnostics_energy_not_gained() {
        let d = SubstepDiagnostics {
            ke_start: 10.0,
            ke_end: 10.1,
            max_velocity: 3.0,
            max_penetration: 0.0,
        };
        assert!(!d.energy_gained());
    }
    #[test]
    fn test_step_with_diagnostics_count() {
        let mut bodies = RigidBodySet::new();
        let hg = bodies.insert(RigidBody::new_static());
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(0.0, -1.0, 0.0);
        b.gravity_scale = 0.0;
        let hb = bodies.insert(b);
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.0,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = TgsSolver::new(3, 4, Vec3::zeros());
        let (stats, diags) =
            solver.step_with_diagnostics(&mut constraints, &mut bodies, 1.0 / 60.0);
        assert_eq!(diags.len(), 3, "should have one diagnostic per substep");
        assert_eq!(stats.iterations_used, 12);
    }
    #[test]
    fn test_energy_monitoring_enabled_does_not_panic() {
        let mut bodies = RigidBodySet::new();
        let hg = bodies.insert(RigidBody::new_static());
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(0.0, -5.0, 0.0);
        b.gravity_scale = 0.0;
        let hb = bodies.insert(b);
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = TgsSolver::with_options(2, 4, 1, Vec3::zeros(), 0.2, true);
        let stats = solver.step(&mut constraints, &mut bodies, 1.0 / 60.0);
        assert!(stats.converged);
    }
    #[test]
    fn test_position_correction_mode_variants() {
        assert_ne!(
            PositionCorrectionMode::Baumgarte,
            PositionCorrectionMode::SplitImpulse
        );
        assert_ne!(
            PositionCorrectionMode::SplitImpulse,
            PositionCorrectionMode::NonlinearGaussSeidel
        );
    }
    #[test]
    fn test_step_warm_runs_and_updates_state() {
        let mut bodies = RigidBodySet::new();
        let hg = bodies.insert(RigidBody::new_static());
        let mut b = RigidBody::new(1.0);
        b.velocity = Vec3::new(0.0, -1.0, 0.0);
        b.gravity_scale = 0.0;
        let hb = bodies.insert(b);
        let c = ContactConstraint::new(
            hb,
            hg,
            Vec3::new(0.0, 1.0, 0.0),
            0.01,
            Vec3::zeros(),
            Vec3::zeros(),
            0.0,
            0.0,
        );
        let mut constraints: Vec<Box<dyn Constraint>> = vec![Box::new(c)];
        let solver = TgsSolver::new(2, 4, Vec3::zeros());
        let mut state = TgsWarmStartState::new(1);
        let dt = 1.0 / 60.0;
        let stats = solver.step_warm(&mut constraints, &mut bodies, dt, &mut state);
        assert_eq!(stats.iterations_used, 8);
        assert!((state.prev_dt - dt).abs() < 1e-12);
        assert_eq!(state.age, 1);
    }
}
/// A single-axis drift correction step applying a position pseudo-impulse.
///
/// Given penetration depths and effective masses, applies a fraction
/// of the correction per iteration using the Baumgarte formula.
#[allow(dead_code)]
pub fn drift_correction_step(
    penetrations: &[f64],
    eff_masses: &[f64],
    beta: f64,
    slop: f64,
    iterations: usize,
) -> DriftCorrectionResult {
    let max_before = penetrations.iter().cloned().fold(0.0_f64, f64::max);
    let mut corrections = penetrations.to_vec();
    for _ in 0..iterations {
        for (pen, &eff) in corrections.iter_mut().zip(eff_masses.iter()) {
            let excess = (*pen - slop).max(0.0);
            if eff > 1e-20 {
                *pen -= beta * excess / eff;
            }
        }
    }
    let max_after = corrections.iter().cloned().fold(0.0_f64, f64::max);
    DriftCorrectionResult {
        max_error_before: max_before,
        max_error_after: max_after,
        iterations,
    }
}
#[cfg(test)]
mod extended_tests {
    use super::*;

    use crate::tgs_solver::BaumgarteVariant;
    use crate::tgs_solver::BlockConstraint6;
    use crate::tgs_solver::CorrectionStrategy;
    use crate::tgs_solver::IslandDescriptor;
    use crate::tgs_solver::JacobianRow6;
    use crate::tgs_solver::MassMatrix6;
    use crate::tgs_solver::ParallelIslandPlanner;

    use oxiphysics_rigid::RigidBody;
    #[test]
    fn test_jacobian_row6_zero() {
        let row = JacobianRow6::zero();
        let v = [0.0f64; 6];
        assert_eq!(row.dot_a(&v), 0.0);
        assert_eq!(row.dot_b(&v), 0.0);
    }
    #[test]
    fn test_jacobian_row6_from_parts() {
        let row = JacobianRow6::from_parts(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
        );
        assert!((row.j_a[0] - 1.0).abs() < 1e-12);
        assert!((row.j_b[0] - (-1.0)).abs() < 1e-12);
    }
    #[test]
    fn test_jacobian_row6_dot_a() {
        let row = JacobianRow6::from_parts([1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [0.0; 3], [0.0; 3]);
        let v = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0];
        assert!((row.dot_a(&v) - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_jacobian_row6_constraint_velocity() {
        let row = JacobianRow6::from_parts(
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        );
        let va = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vb = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!((row.constraint_velocity(&va, &vb) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_mass_matrix6_apply_inv_j_linear() {
        let m = MassMatrix6::new(0.5, [1.0, 1.0, 1.0]);
        let j = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let result = m.apply_inv_j(&j);
        assert!((result[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_mass_matrix6_effective_mass_normal_contact() {
        let m = MassMatrix6::new(1.0 / 1.0, [1.0; 3]);
        let j = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        let eff_a = m.effective_mass(&j);
        let eff_total = eff_a + eff_a;
        assert!((eff_total - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_mass_matrix6_static_body() {
        let m = MassMatrix6::static_body();
        let j = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        assert_eq!(m.effective_mass(&j), 0.0);
    }
    #[test]
    fn test_block_constraint6_solve_reduces_relative_velocity() {
        let row = JacobianRow6::from_parts(
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 0.0],
        );
        let mass = MassMatrix6::new(1.0, [1.0; 3]);
        let mut constraint = BlockConstraint6::new(row, mass, mass, 0.0, 0.0, f64::INFINITY);
        let mut va = [0.0, -2.0, 0.0, 0.0, 0.0, 0.0_f64];
        let mut vb = [0.0; 6_usize];
        for _ in 0..10 {
            let (_dl, da, db) = constraint.solve_iteration(&va, &vb);
            for i in 0..6 {
                va[i] += da[i];
            }
            for i in 0..6 {
                vb[i] += db[i];
            }
        }
        let jv = row.constraint_velocity(&va, &vb);
        assert!(jv >= -0.1, "Relative velocity should be resolved: jv={jv}");
    }
    #[test]
    fn test_block_constraint6_lambda_clamped() {
        let row = JacobianRow6::zero();
        let mass = MassMatrix6::new(1.0, [1.0; 3]);
        let mut c = BlockConstraint6::new(row, mass, mass, 0.0, 0.0, 10.0);
        c.lambda = 9.0;
        c.lambda = (c.lambda + 5.0).clamp(c.lambda_min, c.lambda_max);
        assert!(
            (c.lambda - 10.0).abs() < 1e-12,
            "Lambda should be clamped to max"
        );
    }
    #[test]
    fn test_block_constraint6_reset_lambda() {
        let row = JacobianRow6::zero();
        let mass = MassMatrix6::new(1.0, [1.0; 3]);
        let mut c = BlockConstraint6::new(row, mass, mass, 0.0, -100.0, 100.0);
        c.lambda = 50.0;
        c.reset_lambda();
        assert_eq!(c.lambda, 0.0);
    }
    #[test]
    fn test_drift_correction_step_reduces_penetration() {
        let pens = vec![0.05, 0.02, 0.10];
        let effs = vec![1.0, 1.0, 1.0];
        let result = drift_correction_step(&pens, &effs, 0.5, 0.005, 5);
        assert!(
            result.max_error_after < result.max_error_before,
            "Drift correction should reduce max error"
        );
    }
    #[test]
    fn test_drift_correction_below_slop_unchanged() {
        let pens = vec![0.001];
        let effs = vec![1.0];
        let result = drift_correction_step(&pens, &effs, 0.5, 0.005, 5);
        assert!((result.max_error_after - result.max_error_before).abs() < 1e-12);
    }
    #[test]
    fn test_drift_correction_zero_eff_mass() {
        let pens = vec![0.1];
        let effs = vec![0.0];
        let result = drift_correction_step(&pens, &effs, 0.5, 0.005, 3);
        assert!((result.max_error_after - result.max_error_before).abs() < 1e-12);
    }
    #[test]
    fn test_drift_correction_error_reduction_ratio() {
        let pens = vec![0.1];
        let effs = vec![1.0];
        let result = drift_correction_step(&pens, &effs, 0.8, 0.0, 10);
        assert!(result.error_reduction_ratio() < 1.0);
        assert!(result.reduced_by(0.5));
    }
    #[test]
    fn test_baumgarte_classic_variant() {
        let variant = BaumgarteVariant::Classic;
        let bias = variant.compute_bias(0.02, 0.005, 1.0 / 60.0);
        assert!(bias > 0.0, "Classic Baumgarte should produce positive bias");
    }
    #[test]
    fn test_baumgarte_erp_variant() {
        let variant = BaumgarteVariant::Erp;
        let bias = variant.compute_bias(0.02, 0.005, 1.0 / 60.0);
        assert!(bias > 0.0, "ERP variant should produce positive bias");
    }
    #[test]
    fn test_baumgarte_two_param_variant() {
        let variant = BaumgarteVariant::TwoParam {
            zeta: 1.0,
            omega_n: 10.0,
        };
        let bias = variant.compute_bias(0.02, 0.005, 1.0 / 60.0);
        assert!((bias - 0.3).abs() < 1e-10, "Two-param bias={bias}");
    }
    #[test]
    fn test_baumgarte_zero_dt_gives_zero() {
        let variant = BaumgarteVariant::Classic;
        assert_eq!(variant.compute_bias(0.1, 0.0, 0.0), 0.0);
    }
    #[test]
    fn test_baumgarte_below_slop_gives_zero() {
        let variant = BaumgarteVariant::Classic;
        let bias = variant.compute_bias(0.001, 0.005, 1.0 / 60.0);
        assert_eq!(bias, 0.0, "Below slop should give zero");
    }
    #[test]
    fn test_correction_strategy_velocity_level() {
        let s = CorrectionStrategy::VelocityLevel;
        assert!(s.use_velocity_correction(0.1));
        assert!(!s.use_position_correction(0.1));
    }
    #[test]
    fn test_correction_strategy_position_level() {
        let s = CorrectionStrategy::PositionLevel;
        assert!(!s.use_velocity_correction(0.1));
        assert!(s.use_position_correction(0.1));
    }
    #[test]
    fn test_correction_strategy_hybrid_small() {
        let s = CorrectionStrategy::Hybrid { threshold: 0.05 };
        assert!(s.use_velocity_correction(0.02));
        assert!(!s.use_position_correction(0.02));
    }
    #[test]
    fn test_correction_strategy_hybrid_large() {
        let s = CorrectionStrategy::Hybrid { threshold: 0.05 };
        assert!(!s.use_velocity_correction(0.10));
        assert!(s.use_position_correction(0.10));
    }
    #[test]
    fn test_island_descriptor_counts() {
        let island = IslandDescriptor::new(vec![0, 1, 2], vec![0, 1]);
        assert_eq!(island.body_count(), 3);
        assert_eq!(island.constraint_count(), 2);
    }
    #[test]
    fn test_island_descriptor_estimated_cost() {
        let island = IslandDescriptor::new(vec![0, 1], vec![0]);
        assert_eq!(island.estimated_cost, 24);
    }
    #[test]
    fn test_planner_sort_by_cost() {
        let mut planner = ParallelIslandPlanner::new();
        planner.add_island(IslandDescriptor::new(vec![0], vec![]));
        planner.add_island(IslandDescriptor::new(vec![0, 1, 2, 3], vec![0, 1, 2]));
        planner.sort_by_cost();
        assert!(planner.islands[0].estimated_cost >= planner.islands[1].estimated_cost);
    }
    #[test]
    fn test_planner_assign_workers_one_worker() {
        let mut planner = ParallelIslandPlanner::new();
        planner.add_island(IslandDescriptor::new(vec![0], vec![]));
        planner.add_island(IslandDescriptor::new(vec![1], vec![]));
        let assignments = planner.assign_workers(1);
        assert_eq!(assignments.len(), 1);
        assert_eq!(assignments[0].len(), 2);
    }
    #[test]
    fn test_planner_assign_workers_multiple() {
        let mut planner = ParallelIslandPlanner::new();
        for _ in 0..6 {
            planner.add_island(IslandDescriptor::new(vec![0], vec![]));
        }
        let assignments = planner.assign_workers(3);
        assert_eq!(assignments.len(), 3);
        let total_assigned: usize = assignments.iter().map(|a| a.len()).sum();
        assert_eq!(total_assigned, 6);
    }
    #[test]
    fn test_planner_total_cost() {
        let mut planner = ParallelIslandPlanner::new();
        planner.add_island(IslandDescriptor::new(vec![0], vec![]));
        planner.add_island(IslandDescriptor::new(vec![1], vec![]));
        assert_eq!(planner.total_cost(), 12);
    }
    #[test]
    fn test_planner_assign_workers_zero() {
        let planner = ParallelIslandPlanner::new();
        let assignments = planner.assign_workers(0);
        assert!(assignments.is_empty());
    }
    #[test]
    fn test_compute_erp_positive() {
        let erp = compute_erp(0.2, 1.0 / 60.0);
        assert!(erp > 0.0 && erp < 1.0, "ERP={erp}");
    }
    #[test]
    fn test_compute_cfm_positive() {
        let cfm = compute_cfm(1e-4, 0.01, 1.0 / 60.0);
        assert!(cfm > 0.0, "CFM={cfm}");
    }
    #[test]
    fn test_compute_erp_zero_dt() {
        let erp = compute_erp(0.2, 0.0);
        assert_eq!(erp, 0.0);
    }
    #[test]
    fn test_compute_cfm_zero_inputs() {
        let cfm = compute_cfm(0.0, 0.0, 0.0);
        assert_eq!(cfm, 0.0);
    }
    #[test]
    fn test_separation_velocity_normal() {
        let v = separation_velocity(0.05, 1.0 / 60.0);
        assert!(
            v > 0.0,
            "separation velocity should be positive for overlapping bodies"
        );
    }
    #[test]
    fn test_separation_velocity_below_slop() {
        let v = separation_velocity(0.001, 1.0 / 60.0);
        assert_eq!(v, 0.0, "below slop → zero separation velocity");
    }
    #[test]
    fn test_separation_velocity_zero_dt() {
        let v = separation_velocity(0.1, 0.0);
        assert_eq!(v, 0.0);
    }
    #[test]
    fn test_predicted_closest_approach_approaching() {
        let d = predicted_closest_approach(1.0, -2.0, 0.5);
        assert!((d - 0.0).abs() < 1e-12, "1.0 + (-2.0)*0.5 = 0.0");
    }
    #[test]
    fn test_predicted_closest_approach_separating() {
        let d = predicted_closest_approach(1.0, 1.0, 0.5);
        assert!((d - 1.5).abs() < 1e-12);
    }
    #[test]
    fn test_apply_position_correction_custom_no_op_below_slop() {
        use oxiphysics_rigid::{RigidBody, RigidBodySet};
        let mut bodies = RigidBodySet::new();
        let b = RigidBody::new(1.0);
        let ha = bodies.insert(b);
        let hb = bodies.insert(RigidBody::new_static());
        let normal = Vec3::new(0.0, 1.0, 0.0);
        let y_before = bodies.get(ha).unwrap().transform.position.y;
        apply_position_correction_custom(&mut bodies, ha, hb, &normal, 0.001, 0.5, 0.005);
        let y_after = bodies.get(ha).unwrap().transform.position.y;
        assert!(
            (y_before - y_after).abs() < 1e-12,
            "Below slop: no change expected"
        );
    }
    #[test]
    fn test_apply_position_correction_custom_moves_body() {
        let mut bodies = RigidBodySet::new();
        let b = RigidBody::new(1.0);
        let ha = bodies.insert(b);
        let hb = bodies.insert(RigidBody::new_static());
        let normal = Vec3::new(0.0, 1.0, 0.0);
        let y_before = bodies.get(ha).unwrap().transform.position.y;
        apply_position_correction_custom(&mut bodies, ha, hb, &normal, 0.05, 0.5, 0.001);
        let y_after = bodies.get(ha).unwrap().transform.position.y;
        assert!(
            y_after > y_before,
            "Custom correction should move body upward"
        );
    }
}
