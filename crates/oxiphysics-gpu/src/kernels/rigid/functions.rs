//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{AccumulatedImpulse, RigidBodyState};

/// Quaternion multiplication: `q_out = q_a * q_b`.
///
/// Both quaternions are `[qx, qy, qz, qw]`.
pub(super) fn quat_mul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    let [ax, ay, az, aw] = a;
    let [bx, by, bz, bw] = b;
    [
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
        aw * bw - ax * bx - ay * by - az * bz,
    ]
}
/// Normalise a quaternion to unit length.
pub(super) fn quat_normalise(q: [f64; 4]) -> [f64; 4] {
    let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if len < 1e-30 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    [q[0] / len, q[1] / len, q[2] / len, q[3] / len]
}
/// Compute dq/dt = 0.5 * \[omega, 0\] * q (pure quaternion derivative).
pub(super) fn quat_derivative(q: [f64; 4], omega: [f64; 3]) -> [f64; 4] {
    let omega_q = [omega[0], omega[1], omega[2], 0.0];
    let dq = quat_mul(omega_q, q);
    [dq[0] * 0.5, dq[1] * 0.5, dq[2] * 0.5, dq[3] * 0.5]
}
/// Rotate a vector by a quaternion: v' = q * v * q^{-1}.
pub(super) fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let [qx, qy, qz, qw] = q;
    let cx = qy * v[2] - qz * v[1];
    let cy = qz * v[0] - qx * v[2];
    let cz = qx * v[1] - qy * v[0];
    [
        v[0] + 2.0 * (qw * cx + qy * cz - qz * cy),
        v[1] + 2.0 * (qw * cy + qz * cx - qx * cz),
        v[2] + 2.0 * (qw * cz + qx * cy - qy * cx),
    ]
}
/// Euler-step all rigid bodies by `dt` given external forces and torques.
pub fn integrate_euler(
    states: &mut [RigidBodyState],
    forces: &[[f64; 3]],
    torques: &[[f64; 3]],
    dt: f64,
) {
    for (s, (f, tau)) in states.iter_mut().zip(forces.iter().zip(torques.iter())) {
        let acc = [
            f[0] * s.inverse_mass,
            f[1] * s.inverse_mass,
            f[2] * s.inverse_mass,
        ];
        s.velocity[0] += acc[0] * dt;
        s.velocity[1] += acc[1] * dt;
        s.velocity[2] += acc[2] * dt;
        s.position[0] += s.velocity[0] * dt;
        s.position[1] += s.velocity[1] * dt;
        s.position[2] += s.velocity[2] * dt;
        let alpha = [
            tau[0] * s.inverse_mass,
            tau[1] * s.inverse_mass,
            tau[2] * s.inverse_mass,
        ];
        s.angular_velocity[0] += alpha[0] * dt;
        s.angular_velocity[1] += alpha[1] * dt;
        s.angular_velocity[2] += alpha[2] * dt;
        let dq = quat_derivative(s.orientation, s.angular_velocity);
        s.orientation[0] += dq[0] * dt;
        s.orientation[1] += dq[1] * dt;
        s.orientation[2] += dq[2] * dt;
        s.orientation[3] += dq[3] * dt;
        s.orientation = quat_normalise(s.orientation);
    }
}
/// RK4 integration of a single rigid body for one time step.
pub fn integrate_rk4(
    state: &RigidBodyState,
    force: [f64; 3],
    torque: [f64; 3],
    dt: f64,
) -> RigidBodyState {
    let deriv = |s: &RigidBodyState| -> ([f64; 3], [f64; 3], [f64; 4], [f64; 3]) {
        let dpos = s.velocity;
        let dvel = [
            force[0] * s.inverse_mass,
            force[1] * s.inverse_mass,
            force[2] * s.inverse_mass,
        ];
        let dq = quat_derivative(s.orientation, s.angular_velocity);
        let domega = [
            torque[0] * s.inverse_mass,
            torque[1] * s.inverse_mass,
            torque[2] * s.inverse_mass,
        ];
        (dpos, dvel, dq, domega)
    };
    let step = |s: &RigidBodyState,
                dp: [f64; 3],
                dv: [f64; 3],
                dq: [f64; 4],
                dw: [f64; 3],
                h: f64|
     -> RigidBodyState {
        let mut ns = *s;
        for k in 0..3 {
            ns.position[k] = s.position[k] + dp[k] * h;
            ns.velocity[k] = s.velocity[k] + dv[k] * h;
            ns.angular_velocity[k] = s.angular_velocity[k] + dw[k] * h;
        }
        for k in 0..4 {
            ns.orientation[k] = s.orientation[k] + dq[k] * h;
        }
        ns.orientation = quat_normalise(ns.orientation);
        ns
    };
    let (dp1, dv1, dq1, dw1) = deriv(state);
    let s2 = step(state, dp1, dv1, dq1, dw1, dt * 0.5);
    let (dp2, dv2, dq2, dw2) = deriv(&s2);
    let s3 = step(state, dp2, dv2, dq2, dw2, dt * 0.5);
    let (dp3, dv3, dq3, dw3) = deriv(&s3);
    let s4 = step(state, dp3, dv3, dq3, dw3, dt);
    let (dp4, dv4, dq4, dw4) = deriv(&s4);
    let mut out = *state;
    for k in 0..3 {
        out.position[k] =
            state.position[k] + dt / 6.0 * (dp1[k] + 2.0 * dp2[k] + 2.0 * dp3[k] + dp4[k]);
        out.velocity[k] =
            state.velocity[k] + dt / 6.0 * (dv1[k] + 2.0 * dv2[k] + 2.0 * dv3[k] + dv4[k]);
        out.angular_velocity[k] =
            state.angular_velocity[k] + dt / 6.0 * (dw1[k] + 2.0 * dw2[k] + 2.0 * dw3[k] + dw4[k]);
    }
    for k in 0..4 {
        out.orientation[k] =
            state.orientation[k] + dt / 6.0 * (dq1[k] + 2.0 * dq2[k] + 2.0 * dq3[k] + dq4[k]);
    }
    out.orientation = quat_normalise(out.orientation);
    out
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::compute::ComputeKernel;
    use crate::kernels::rigid::Aabb;

    use crate::kernels::rigid::BroadphaseUpdateKernel;
    use crate::kernels::rigid::ConstraintSolverKernel;

    use crate::kernels::rigid::ContactGenerationKernel;
    use crate::kernels::rigid::ContactPoint;
    use crate::kernels::rigid::DistanceConstraint;
    use crate::kernels::rigid::IntegratePositionKernel;
    use crate::kernels::rigid::IntegrateVelocityKernel;
    use crate::kernels::rigid::IslandSolver;
    use crate::kernels::rigid::QuaternionNormKernel;
    use crate::kernels::rigid::RigidBodyState;
    use crate::kernels::rigid::SemiImplicitEulerKernel;

    use crate::kernels::rigid::SoaRigidBody;

    #[test]
    fn test_rigid_gravity_integration() {
        let n = 10_usize;
        let g = -9.81_f64;
        let dt = 0.1_f64;
        let vel = vec![0.0_f64; n * 3];
        let mut force = vec![0.0_f64; n * 3];
        for i in 0..n {
            force[i * 3 + 1] = g;
        }
        let inv_mass = vec![1.0_f64; n];
        let dt_slice = vec![dt];
        let mut outputs = vec![Vec::new()];
        IntegrateVelocityKernel.execute(&[&vel, &force, &inv_mass, &dt_slice], &mut outputs, n);
        assert_eq!(outputs[0].len(), n * 3, "output length should be 3n");
        let expected_vy = g * dt;
        for i in 0..n {
            let vx = outputs[0][i * 3];
            let vy = outputs[0][i * 3 + 1];
            let vz = outputs[0][i * 3 + 2];
            assert!(vx.abs() < 1e-15, "body {i}: vx should be 0, got {vx}");
            assert!(
                (vy - expected_vy).abs() < 1e-12,
                "body {i}: vy should be {expected_vy}, got {vy}"
            );
            assert!(vz.abs() < 1e-15, "body {i}: vz should be 0, got {vz}");
        }
    }
    #[test]
    fn integrate_velocity_updates_correctly() {
        let vel = vec![1.0, 0.0, 0.0];
        let force = vec![10.0, 0.0, 0.0];
        let inv_mass = vec![0.5];
        let dt = vec![0.1];
        let mut outputs = vec![Vec::new()];
        IntegrateVelocityKernel.execute(&[&vel, &force, &inv_mass, &dt], &mut outputs, 1);
        assert!((outputs[0][0] - 1.5).abs() < 1e-12);
        assert!((outputs[0][1]).abs() < 1e-12);
        assert!((outputs[0][2]).abs() < 1e-12);
    }
    #[test]
    fn integrate_position_updates_correctly() {
        let pos = vec![0.0, 0.0, 0.0];
        let vel = vec![3.0, 4.0, 0.0];
        let dt = vec![0.5];
        let mut outputs = vec![Vec::new()];
        IntegratePositionKernel.execute(&[&pos, &vel, &dt], &mut outputs, 1);
        assert!((outputs[0][0] - 1.5).abs() < 1e-12);
        assert!((outputs[0][1] - 2.0).abs() < 1e-12);
        assert!((outputs[0][2]).abs() < 1e-12);
    }
    #[test]
    fn integrate_euler_free_fall() {
        let g = -9.81_f64;
        let dt = 0.1_f64;
        let mut state = RigidBodyState::at_rest([0.0, 0.0, 0.0], 1.0);
        let force = [0.0, g, 0.0];
        let torque = [0.0; 3];
        integrate_euler(std::slice::from_mut(&mut state), &[force], &[torque], dt);
        let expected_vy = g * dt;
        assert!((state.velocity[1] - expected_vy).abs() < 1e-12);
        let expected_y = expected_vy * dt;
        assert!((state.position[1] - expected_y).abs() < 1e-12);
        assert!(state.position[0].abs() < 1e-15);
        assert!(state.position[2].abs() < 1e-15);
    }
    #[test]
    fn integrate_rk4_free_fall() {
        let g = -9.81_f64;
        let dt = 0.1_f64;
        let state = RigidBodyState::at_rest([0.0, 0.0, 0.0], 1.0);
        let force = [0.0, g, 0.0];
        let torque = [0.0; 3];
        let new_state = integrate_rk4(&state, force, torque, dt);
        let expected_vy = g * dt;
        let expected_y = 0.5 * g * dt * dt;
        assert!((new_state.velocity[1] - expected_vy).abs() < 1e-10);
        assert!((new_state.position[1] - expected_y).abs() < 1e-10);
    }
    #[test]
    fn integrate_euler_orientation_stays_normalised() {
        let mut state = RigidBodyState::at_rest([0.0, 0.0, 0.0], 1.0);
        let torque = [0.0, 0.0, 1.0];
        let force = [0.0; 3];
        let dt = 0.01;
        for _ in 0..100 {
            integrate_euler(std::slice::from_mut(&mut state), &[force], &[torque], dt);
        }
        let q = state.orientation;
        let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
        assert!((len - 1.0).abs() < 1e-10, "quaternion norm = {len}");
    }
    /// Distance constraint should pull bodies together.
    #[test]
    fn test_distance_constraint_converges() {
        let mut positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let inv_masses = [1.0, 1.0];
        let constraints = [DistanceConstraint {
            body_a: 0,
            body_b: 1,
            rest_length: 1.0,
            compliance: 0.0,
        }];
        for _ in 0..10 {
            ConstraintSolverKernel::solve_distance_constraints(
                &mut positions,
                &inv_masses,
                &constraints,
                0.01,
            );
        }
        let dx = positions[1][0] - positions[0][0];
        let dist = dx.abs();
        assert!(
            (dist - 1.0).abs() < 0.1,
            "distance should approach 1.0, got {dist}"
        );
    }
    /// Static body should not move under constraint.
    #[test]
    fn test_constraint_static_body() {
        let mut positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let inv_masses = [0.0, 1.0];
        let constraints = [DistanceConstraint {
            body_a: 0,
            body_b: 1,
            rest_length: 1.0,
            compliance: 0.0,
        }];
        ConstraintSolverKernel::solve_distance_constraints(
            &mut positions,
            &inv_masses,
            &constraints,
            0.01,
        );
        assert!((positions[0][0]).abs() < 1e-14, "static body moved!");
    }
    /// Overlapping AABBs should be detected.
    #[test]
    fn test_aabb_overlap() {
        let a = Aabb::from_center([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::from_center([1.5, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(a.overlaps(&b));
    }
    /// Separated AABBs should not overlap.
    #[test]
    fn test_aabb_no_overlap() {
        let a = Aabb::from_center([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::from_center([5.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(!a.overlaps(&b));
    }
    /// Broadphase should find overlapping pairs.
    #[test]
    fn test_broadphase_finds_pairs() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let half_extents = [[1.0, 1.0, 1.0]; 3];
        let pairs = BroadphaseUpdateKernel::find_overlapping_pairs(&positions, &half_extents, 0.0);
        assert!(pairs.contains(&(0, 1)), "bodies 0 and 1 should overlap");
        assert!(
            !pairs.contains(&(0, 2)),
            "bodies 0 and 2 should not overlap"
        );
    }
    /// AABB volume.
    #[test]
    fn test_aabb_volume() {
        let a = Aabb::from_center([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        let vol = a.volume();
        assert!((vol - 48.0).abs() < 1e-12, "volume = {vol}, expected 48");
    }
    /// Two separate pairs should form two islands.
    #[test]
    fn test_island_solver_two_islands() {
        let constraints = [(0, 1), (2, 3)];
        let islands = IslandSolver::build(4, &constraints);
        assert_eq!(islands.num_islands, 2);
        assert_eq!(islands.island_ids[0], islands.island_ids[1]);
        assert_eq!(islands.island_ids[2], islands.island_ids[3]);
        assert_ne!(islands.island_ids[0], islands.island_ids[2]);
    }
    /// Connected chain should form one island.
    #[test]
    fn test_island_solver_chain() {
        let constraints = [(0, 1), (1, 2), (2, 3)];
        let islands = IslandSolver::build(4, &constraints);
        assert_eq!(islands.num_islands, 1);
        let all_same = islands
            .island_ids
            .iter()
            .all(|&id| id == islands.island_ids[0]);
        assert!(all_same, "all bodies should be in the same island");
    }
    /// Bodies in island query.
    #[test]
    fn test_bodies_in_island() {
        let constraints = [(0, 1), (2, 3)];
        let islands = IslandSolver::build(4, &constraints);
        let island_0 = islands.bodies_in_island(islands.island_ids[0]);
        assert_eq!(island_0.len(), 2);
        assert!(island_0.contains(&0));
        assert!(island_0.contains(&1));
    }
    /// Overlapping spheres should generate a contact.
    #[test]
    fn test_contact_generation_overlap() {
        let positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let radii = [1.0, 1.0];
        let pairs = [(0, 1)];
        let contacts =
            ContactGenerationKernel::generate_sphere_contacts(&positions, &radii, &pairs);
        assert_eq!(contacts.len(), 1);
        assert!(contacts[0].depth > 0.0, "depth should be positive");
        assert!(
            (contacts[0].depth - 0.5).abs() < 1e-12,
            "depth = {}",
            contacts[0].depth
        );
    }
    /// Separated spheres should not generate contacts.
    #[test]
    fn test_contact_generation_no_overlap() {
        let positions = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let radii = [1.0, 1.0];
        let pairs = [(0, 1)];
        let contacts =
            ContactGenerationKernel::generate_sphere_contacts(&positions, &radii, &pairs);
        assert!(contacts.is_empty());
    }
    /// Contact resolution should separate overlapping bodies.
    #[test]
    fn test_contact_resolution() {
        let mut positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let mut velocities = [[0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let inv_masses = [1.0, 1.0];
        let contacts = [ContactPoint {
            position: [0.75, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            depth: 0.5,
            body_a: 0,
            body_b: 1,
        }];
        ContactGenerationKernel::resolve_contacts(
            &mut positions,
            &mut velocities,
            &inv_masses,
            &contacts,
            0.5,
        );
        let dx = positions[1][0] - positions[0][0];
        assert!(dx > 1.5, "bodies should be pushed apart, dx = {dx}");
    }
    /// Semi-implicit Euler should update both velocity and position.
    #[test]
    fn test_semi_implicit_euler() {
        let pos = vec![0.0, 0.0, 0.0];
        let vel = vec![0.0, 0.0, 0.0];
        let force = vec![10.0, 0.0, 0.0];
        let inv_mass = vec![1.0];
        let dt = vec![0.1];
        let mut outputs = vec![Vec::new(), Vec::new()];
        SemiImplicitEulerKernel.execute(&[&pos, &vel, &force, &inv_mass, &dt], &mut outputs, 1);
        assert!((outputs[0][0] - 1.0).abs() < 1e-12, "v = {}", outputs[0][0]);
        assert!((outputs[1][0] - 0.1).abs() < 1e-12, "p = {}", outputs[1][0]);
    }
    /// Quaternion rotation should preserve vector length.
    #[test]
    fn test_quat_rotate_preserves_length() {
        let angle = std::f64::consts::FRAC_PI_4;
        let q = [0.0, 0.0, angle.sin(), angle.cos()];
        let v = [1.0, 0.0, 0.0];
        let rotated = quat_rotate(q, v);
        let len =
            (rotated[0] * rotated[0] + rotated[1] * rotated[1] + rotated[2] * rotated[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10, "rotated length = {len}");
    }
    #[test]
    fn test_soa_rigid_body_from_vec() {
        let states = vec![
            RigidBodyState::at_rest([0.0, 0.0, 0.0], 1.0),
            RigidBodyState::at_rest([1.0, 0.0, 0.0], 0.5),
        ];
        let soa = SoaRigidBody::from_slice(&states);
        assert_eq!(soa.count, 2);
        assert!((soa.pos_x[0] - 0.0).abs() < 1e-14);
        assert!((soa.pos_x[1] - 1.0).abs() < 1e-14);
        assert!((soa.inv_mass[0] - 1.0).abs() < 1e-14);
        assert!((soa.inv_mass[1] - 0.5).abs() < 1e-14);
    }
    #[test]
    fn test_soa_rigid_body_to_vec() {
        let states = vec![
            RigidBodyState::at_rest([1.0, 2.0, 3.0], 0.25),
            RigidBodyState::at_rest([4.0, 5.0, 6.0], 0.1),
        ];
        let soa = SoaRigidBody::from_slice(&states);
        let back = soa.to_vec();
        assert!((back[0].position[0] - 1.0).abs() < 1e-14);
        assert!((back[1].position[1] - 5.0).abs() < 1e-14);
    }
    #[test]
    fn test_soa_integrate_euler_gravity() {
        let states = vec![RigidBodyState::at_rest([0.0, 10.0, 0.0], 1.0)];
        let mut soa = SoaRigidBody::from_slice(&states);
        let forces = vec![[0.0f64, -9.81, 0.0]];
        let torques = vec![[0.0f64; 3]];
        soa.integrate_euler(&forces, &torques, 0.1);
        let back = soa.to_vec();
        assert!(
            back[0].velocity[1] < 0.0,
            "velocity should be negative after gravity"
        );
        assert!(back[0].position[1] < 10.0, "y should decrease");
    }
    #[test]
    fn test_soa_quaternion_stays_normalised() {
        let states = vec![RigidBodyState::at_rest([0.0; 3], 1.0)];
        let mut soa = SoaRigidBody::from_slice(&states);
        let forces = vec![[0.0f64; 3]];
        let torques = vec![[0.0, 0.0, 2.0]];
        for _ in 0..50 {
            soa.integrate_euler(&forces, &torques, 0.01);
        }
        let q = [soa.quat_x[0], soa.quat_y[0], soa.quat_z[0], soa.quat_w[0]];
        let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
        assert!((len - 1.0).abs() < 1e-10, "quaternion norm = {len}");
    }
    #[test]
    fn test_soa_angular_velocity_accumulates() {
        let states = vec![RigidBodyState::at_rest([0.0; 3], 1.0)];
        let mut soa = SoaRigidBody::from_slice(&states);
        let forces = vec![[0.0f64; 3]];
        let torques = vec![[0.0, 0.0, 1.0]];
        soa.integrate_euler(&forces, &torques, 0.5);
        let back = soa.to_vec();
        assert!(
            back[0].angular_velocity[2] > 0.0,
            "angular velocity z should increase"
        );
    }
    #[test]
    fn test_quat_norm_kernel_normalizes() {
        let quats = vec![[2.0f64, 0.0, 0.0, 0.0], [0.0, 3.0, 0.0, 0.0]];
        let normed = QuaternionNormKernel::normalize_batch(&quats);
        for q in &normed {
            let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
            assert!((len - 1.0).abs() < 1e-12, "not unit: {len}");
        }
    }
    #[test]
    fn test_quat_norm_kernel_identity_unchanged() {
        let quats = vec![[0.0, 0.0, 0.0, 1.0f64]];
        let normed = QuaternionNormKernel::normalize_batch(&quats);
        assert!((normed[0][3] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_quat_norm_kernel_zero_fallback() {
        let quats = vec![[0.0f64; 4]];
        let normed = QuaternionNormKernel::normalize_batch(&quats);
        assert!((normed[0][3] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_angular_velocity_integration_increases_omega() {
        let mut state = RigidBodyState::at_rest([0.0; 3], 1.0);
        let force = [0.0; 3];
        let torque = [1.0, 0.0, 0.0];
        let dt = 0.1;
        integrate_euler(std::slice::from_mut(&mut state), &[force], &[torque], dt);
        assert!(state.angular_velocity[0] > 0.0, "omega_x should increase");
    }
    #[test]
    fn test_angular_velocity_rk4_matches_euler_roughly() {
        let state = RigidBodyState::at_rest([0.0; 3], 1.0);
        let force = [0.0; 3];
        let torque = [0.0, 1.0, 0.0];
        let dt = 0.001;
        let rk4 = integrate_rk4(&state, force, torque, dt);
        let mut euler = state;
        integrate_euler(std::slice::from_mut(&mut euler), &[force], &[torque], dt);
        let diff = (rk4.angular_velocity[1] - euler.angular_velocity[1]).abs();
        assert!(
            diff < 0.01,
            "RK4 and Euler should roughly agree for small dt, diff={diff}"
        );
    }
    #[test]
    fn test_mock_cuda_integrate_kernel() {
        let n = 100;
        let mut states: Vec<RigidBodyState> = (0..n)
            .map(|_| RigidBodyState::at_rest([0.0, 100.0, 0.0], 1.0))
            .collect();
        let forces: Vec<[f64; 3]> = vec![[0.0, -9.81, 0.0]; n];
        let torques: Vec<[f64; 3]> = vec![[0.0; 3]; n];
        let dt = 1.0 / 60.0;
        for _ in 0..3 {
            integrate_euler(&mut states, &forces, &torques, dt);
        }
        for s in &states {
            assert!(s.position[1] < 100.0, "y should have decreased");
        }
    }
    #[test]
    fn test_mock_cuda_kernel_simd_batch() {
        let n = 64;
        let states: Vec<RigidBodyState> = (0..n)
            .map(|i| RigidBodyState::at_rest([i as f64, 0.0, 0.0], 1.0))
            .collect();
        let soa = SoaRigidBody::from_slice(&states);
        assert_eq!(soa.count, n);
        for i in 0..n {
            assert!((soa.pos_x[i] - i as f64).abs() < 1e-12);
        }
    }
}
/// Integrate a single [`RigidBodyState`] using the semi-implicit Euler method.
///
/// The velocity is updated first (using forces/torques at time *t*), then the
/// position / orientation is updated using the new velocity at time *t + dt*.
/// This is symplectic, energy-preserving in the absence of dissipation, and is
/// the standard method used in real-time physics engines.
pub fn integrate_semi_implicit(
    state: &RigidBodyState,
    force: [f64; 3],
    torque: [f64; 3],
    dt: f64,
) -> RigidBodyState {
    let mut s = *state;
    if s.inverse_mass < 1e-30 {
        return s;
    }
    let a = [
        force[0] * s.inverse_mass,
        force[1] * s.inverse_mass,
        force[2] * s.inverse_mass,
    ];
    s.velocity[0] += a[0] * dt;
    s.velocity[1] += a[1] * dt;
    s.velocity[2] += a[2] * dt;
    s.position[0] += s.velocity[0] * dt;
    s.position[1] += s.velocity[1] * dt;
    s.position[2] += s.velocity[2] * dt;
    let alpha = [
        torque[0] * s.inverse_mass,
        torque[1] * s.inverse_mass,
        torque[2] * s.inverse_mass,
    ];
    s.angular_velocity[0] += alpha[0] * dt;
    s.angular_velocity[1] += alpha[1] * dt;
    s.angular_velocity[2] += alpha[2] * dt;
    let dq = quat_derivative(s.orientation, s.angular_velocity);
    for k in 0..4 {
        s.orientation[k] += dq[k] * dt;
    }
    s.orientation = quat_normalise(s.orientation);
    s
}
/// GPU-batch semi-implicit Euler: integrate all bodies in a slice.
#[allow(clippy::too_many_arguments)]
pub fn batch_integrate_semi_implicit(
    states: &mut [RigidBodyState],
    forces: &[[f64; 3]],
    torques: &[[f64; 3]],
    dt: f64,
) {
    for (s, (f, tau)) in states.iter_mut().zip(forces.iter().zip(torques.iter())) {
        *s = integrate_semi_implicit(s, *f, *tau, dt);
    }
}
/// Integrate only the angular velocity of each body (no linear dynamics).
///
/// This mimics a GPU kernel that handles only rotational degrees of freedom
/// — useful when linear and angular integration are split into separate passes.
pub fn integrate_angular_velocity_only(
    states: &mut [RigidBodyState],
    torques: &[[f64; 3]],
    dt: f64,
) {
    for (s, tau) in states.iter_mut().zip(torques.iter()) {
        if s.inverse_mass < 1e-30 {
            continue;
        }
        let alpha = [
            tau[0] * s.inverse_mass,
            tau[1] * s.inverse_mass,
            tau[2] * s.inverse_mass,
        ];
        s.angular_velocity[0] += alpha[0] * dt;
        s.angular_velocity[1] += alpha[1] * dt;
        s.angular_velocity[2] += alpha[2] * dt;
        let dq = quat_derivative(s.orientation, s.angular_velocity);
        for k in 0..4 {
            s.orientation[k] += dq[k] * dt;
        }
        s.orientation = quat_normalise(s.orientation);
    }
}
/// Compute the world-space AABB of a rigid body given its OBB half-extents.
///
/// The OBB is transformed by the body's orientation quaternion to produce a
/// world AABB that tightly (conservatively) bounds the oriented box.
pub fn compute_world_aabb(
    position: [f64; 3],
    orientation: [f64; 4],
    half_extents: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let signs = [
        [-1.0f64, -1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, -1.0],
        [1.0, 1.0, 1.0],
    ];
    let mut aabb_min = [f64::INFINITY; 3];
    let mut aabb_max = [f64::NEG_INFINITY; 3];
    for s in &signs {
        let local = [
            s[0] * half_extents[0],
            s[1] * half_extents[1],
            s[2] * half_extents[2],
        ];
        let world = quat_rotate(orientation, local);
        for k in 0..3 {
            let v = position[k] + world[k];
            aabb_min[k] = f64::min(aabb_min[k], v);
            aabb_max[k] = f64::max(aabb_max[k], v);
        }
    }
    (aabb_min, aabb_max)
}
/// Batch AABB update kernel.
///
/// Returns a `Vec<(min, max)>` of world AABBs for all bodies.
pub fn batch_update_world_aabbs(
    states: &[RigidBodyState],
    half_extents: &[[f64; 3]],
) -> Vec<([f64; 3], [f64; 3])> {
    assert_eq!(states.len(), half_extents.len());
    states
        .iter()
        .zip(half_extents.iter())
        .map(|(s, he)| compute_world_aabb(s.position, s.orientation, *he))
        .collect()
}
/// Batch impulse application: apply accumulated impulses to all bodies.
pub fn apply_impulses(states: &mut [RigidBodyState], impulses: &[AccumulatedImpulse]) {
    assert_eq!(states.len(), impulses.len());
    for (s, imp) in states.iter_mut().zip(impulses.iter()) {
        imp.apply(s);
    }
}
#[cfg(test)]
mod extended_rigid_tests {

    use crate::kernels::rigid::AccumulatedImpulse;

    use crate::kernels::rigid::ContactBatchProcessor;

    use crate::kernels::rigid::ContactPoint;

    use crate::kernels::rigid::RigidBodyState;

    use crate::kernels::rigid::SleepParams;
    use crate::kernels::rigid::SleepState;
    use crate::kernels::rigid::SleepTest;
    use crate::kernels::rigid::SoaRigidBody;
    use crate::kernels::rigid::apply_impulses;
    use crate::kernels::rigid::batch_integrate_semi_implicit;
    use crate::kernels::rigid::batch_update_world_aabbs;
    use crate::kernels::rigid::compute_world_aabb;
    use crate::kernels::rigid::integrate_angular_velocity_only;
    use crate::kernels::rigid::integrate_semi_implicit;
    #[test]
    fn semi_implicit_linear_motion() {
        let state = RigidBodyState::at_rest([0.0, 0.0, 0.0], 1.0);
        let force = [10.0, 0.0, 0.0];
        let torque = [0.0; 3];
        let dt = 0.1;
        let next = integrate_semi_implicit(&state, force, torque, dt);
        assert!((next.velocity[0] - 1.0).abs() < 1e-12);
        assert!((next.position[0] - 0.1).abs() < 1e-12);
    }
    #[test]
    fn semi_implicit_static_body_unchanged() {
        let state = RigidBodyState::at_rest([5.0, 5.0, 5.0], 0.0);
        let next = integrate_semi_implicit(&state, [100.0; 3], [100.0; 3], 1.0);
        assert_eq!(next.position, state.position);
        assert_eq!(next.velocity, state.velocity);
    }
    #[test]
    fn semi_implicit_gravity_free_fall() {
        let dt = 0.01;
        let mut s = RigidBodyState::at_rest([0.0, 100.0, 0.0], 1.0);
        let force = [0.0, -9.81, 0.0];
        let torque = [0.0; 3];
        for _ in 0..100 {
            s = integrate_semi_implicit(&s, force, torque, dt);
        }
        assert!(
            (s.velocity[1] + 9.81).abs() < 0.01,
            "v_y = {}",
            s.velocity[1]
        );
    }
    #[test]
    fn batch_semi_implicit_matches_single() {
        let n = 5;
        let mut states: Vec<RigidBodyState> = (0..n)
            .map(|i| RigidBodyState::at_rest([i as f64, 0.0, 0.0], 1.0))
            .collect();
        let forces = vec![[1.0, 0.0, 0.0]; n];
        let torques = vec![[0.0; 3]; n];
        let dt = 0.01;
        let singles: Vec<RigidBodyState> = states
            .iter()
            .map(|s| integrate_semi_implicit(s, [1.0, 0.0, 0.0], [0.0; 3], dt))
            .collect();
        batch_integrate_semi_implicit(&mut states, &forces, &torques, dt);
        for (s, e) in states.iter().zip(singles.iter()) {
            for k in 0..3 {
                assert!((s.position[k] - e.position[k]).abs() < 1e-12);
                assert!((s.velocity[k] - e.velocity[k]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn angular_only_changes_omega_not_position() {
        let mut states = vec![RigidBodyState::at_rest([1.0, 2.0, 3.0], 1.0)];
        let torques = vec![[0.0, 0.0, 1.0]];
        integrate_angular_velocity_only(&mut states, &torques, 0.1);
        assert_eq!(states[0].position, [1.0, 2.0, 3.0]);
        assert!(states[0].angular_velocity[2] > 0.0);
    }
    #[test]
    fn angular_only_static_body_unchanged() {
        let mut states = vec![RigidBodyState::at_rest([0.0; 3], 0.0)];
        let torques = vec![[10.0; 3]];
        integrate_angular_velocity_only(&mut states, &torques, 1.0);
        assert_eq!(states[0].angular_velocity, [0.0; 3]);
    }
    #[test]
    fn world_aabb_identity_orientation() {
        let pos = [1.0, 2.0, 3.0];
        let q = [0.0, 0.0, 0.0, 1.0];
        let he = [1.0, 2.0, 3.0];
        let (mn, mx) = compute_world_aabb(pos, q, he);
        assert!((mn[0] - 0.0).abs() < 1e-10);
        assert!((mx[0] - 2.0).abs() < 1e-10);
        assert!((mn[1] - 0.0).abs() < 1e-10);
        assert!((mx[1] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn world_aabb_bounds_contain_corners() {
        let pos = [0.0; 3];
        let q = [0.0, 0.0, 0.0, 1.0];
        let he = [1.0, 1.0, 1.0];
        let (mn, mx) = compute_world_aabb(pos, q, he);
        for k in 0..3 {
            assert!(mn[k] <= -0.99, "min[{k}] = {}", mn[k]);
            assert!(mx[k] >= 0.99, "max[{k}] = {}", mx[k]);
        }
    }
    #[test]
    fn batch_world_aabbs_correct_count() {
        let states = vec![
            RigidBodyState::at_rest([0.0; 3], 1.0),
            RigidBodyState::at_rest([5.0; 3], 1.0),
        ];
        let hes = vec![[1.0; 3], [0.5; 3]];
        let aabbs = batch_update_world_aabbs(&states, &hes);
        assert_eq!(aabbs.len(), 2);
    }
    #[test]
    fn sleep_test_body_falls_asleep() {
        let mut test = SleepTest::new(1);
        let params = SleepParams {
            linear_threshold: 0.1,
            angular_threshold: 0.1,
            sleep_frames: 3,
        };
        let state = RigidBodyState::at_rest([0.0; 3], 1.0);
        for _ in 0..3 {
            test.update(&[state], &params);
        }
        assert_eq!(test.sleeping_count(), 1);
    }
    #[test]
    fn sleep_test_moving_body_stays_awake() {
        let mut test = SleepTest::new(1);
        let params = SleepParams::default();
        let mut state = RigidBodyState::at_rest([0.0; 3], 1.0);
        state.velocity = [1.0, 0.0, 0.0];
        for _ in 0..20 {
            test.update(&[state], &params);
        }
        assert_eq!(test.sleeping_count(), 0);
    }
    #[test]
    fn sleep_test_wake_all_resets() {
        let mut test = SleepTest::new(2);
        let params = SleepParams {
            sleep_frames: 1,
            ..Default::default()
        };
        let state = RigidBodyState::at_rest([0.0; 3], 1.0);
        test.update(&[state, state], &params);
        test.wake_all();
        assert_eq!(test.sleeping_count(), 0);
        assert!(test.dormant_frames.iter().all(|&f| f == 0));
    }
    #[test]
    fn sleep_test_static_body_is_sleeping() {
        let mut test = SleepTest::new(1);
        let params = SleepParams::default();
        let state = RigidBodyState::at_rest([0.0; 3], 0.0);
        test.update(&[state], &params);
        assert_eq!(test.sleep_states[0], SleepState::Sleeping);
    }
    #[test]
    fn accumulated_impulse_apply_changes_velocity() {
        let mut state = RigidBodyState::at_rest([0.0; 3], 1.0);
        let mut imp = AccumulatedImpulse::default();
        imp.add_linear([5.0, 0.0, 0.0]);
        imp.apply(&mut state);
        assert!((state.velocity[0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn accumulated_impulse_magnitude() {
        let mut imp = AccumulatedImpulse::default();
        imp.add_linear([3.0, 4.0, 0.0]);
        assert!((imp.linear_magnitude() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn accumulated_impulse_add_angular() {
        let mut imp = AccumulatedImpulse::default();
        imp.add_angular([0.0, 0.0, 1.0]);
        assert!((imp.angular_magnitude() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn apply_impulses_batch() {
        let n = 3;
        let mut states: Vec<RigidBodyState> = (0..n)
            .map(|_| RigidBodyState::at_rest([0.0; 3], 1.0))
            .collect();
        let impulses: Vec<AccumulatedImpulse> = (0..n)
            .map(|i| {
                let mut imp = AccumulatedImpulse::default();
                imp.add_linear([i as f64, 0.0, 0.0]);
                imp
            })
            .collect();
        apply_impulses(&mut states, &impulses);
        for (i, s) in states.iter().enumerate() {
            assert!(
                (s.velocity[0] - i as f64).abs() < 1e-12,
                "body {i}: v_x = {}",
                s.velocity[0]
            );
        }
    }
    #[test]
    fn contact_batch_generates_nonzero_impulses() {
        let positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let velocities = [[0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let inv_masses = [1.0, 1.0];
        let contacts = [ContactPoint {
            position: [0.75, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            depth: 0.5,
            body_a: 0,
            body_b: 1,
        }];
        let proc = ContactBatchProcessor::new(0.5, 0.2);
        let impulses =
            proc.process_contacts(&contacts, &positions, &velocities, &inv_masses, 2, 0.01);
        assert!(impulses[0].linear_magnitude() > 0.0 || impulses[1].linear_magnitude() > 0.0);
    }
    #[test]
    fn contact_batch_no_contacts_zero_impulse() {
        let positions = [[0.0; 3], [5.0, 0.0, 0.0]];
        let velocities = [[0.0; 3]; 2];
        let inv_masses = [1.0, 1.0];
        let contacts: &[ContactPoint] = &[];
        let proc = ContactBatchProcessor::new(0.5, 0.2);
        let impulses =
            proc.process_contacts(contacts, &positions, &velocities, &inv_masses, 2, 0.01);
        for imp in &impulses {
            assert!((imp.linear_magnitude()).abs() < 1e-12);
        }
    }
    #[test]
    fn world_aabb_90_degree_rotation_x() {
        let angle = std::f64::consts::FRAC_PI_4;
        let q = [angle.sin(), 0.0, 0.0, angle.cos()];
        let he = [1.0, 2.0, 3.0];
        let (mn, mx) = compute_world_aabb([0.0; 3], q, he);
        for k in 0..3 {
            assert!(mn[k] < mx[k], "min[{k}] >= max[{k}]");
        }
    }
    #[test]
    fn soa_batch_integrate_all_fall_equally() {
        let n = 10;
        let states: Vec<RigidBodyState> = (0..n)
            .map(|_| RigidBodyState::at_rest([0.0, 100.0, 0.0], 1.0))
            .collect();
        let mut soa = SoaRigidBody::from_slice(&states);
        let forces = vec![[0.0f64, -9.81, 0.0]; n];
        let torques = vec![[0.0f64; 3]; n];
        soa.integrate_euler(&forces, &torques, 0.1);
        let back = soa.to_vec();
        let vy0 = back[0].velocity[1];
        for s in &back {
            assert!(
                (s.velocity[1] - vy0).abs() < 1e-12,
                "all bodies should fall equally"
            );
        }
    }
    #[test]
    fn sleep_test_not_sleeping_before_threshold() {
        let mut test = SleepTest::new(1);
        let params = SleepParams {
            sleep_frames: 5,
            ..Default::default()
        };
        let state = RigidBodyState::at_rest([0.0; 3], 1.0);
        for _ in 0..4 {
            test.update(&[state], &params);
        }
        assert_eq!(
            test.sleeping_count(),
            0,
            "body should not sleep before threshold"
        );
    }
    #[test]
    fn contact_batch_static_pair_no_impulse() {
        let positions = [[0.0; 3], [0.5, 0.0, 0.0]];
        let velocities = [[0.0; 3]; 2];
        let inv_masses = [0.0, 0.0];
        let contacts = [ContactPoint {
            position: [0.25, 0.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            depth: 0.5,
            body_a: 0,
            body_b: 1,
        }];
        let proc = ContactBatchProcessor::new(0.5, 0.2);
        let impulses =
            proc.process_contacts(&contacts, &positions, &velocities, &inv_masses, 2, 0.01);
        for imp in &impulses {
            assert!(
                (imp.linear_magnitude()).abs() < 1e-12,
                "static pair should produce zero impulse"
            );
        }
    }
}
