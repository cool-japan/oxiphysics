//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::types::*;
#[cfg(test)]
mod tests_dynamics_ext {
    use super::*;
    #[test]
    fn test_spatial_inertia6_thin_rod_mass() {
        let si = SpatialInertia6::thin_rod_x(2.0, 1.0);
        assert_eq!(si.mass, 2.0);
    }
    #[test]
    fn test_spatial_inertia6_thin_rod_com() {
        let si = SpatialInertia6::thin_rod_x(1.0, 2.0);
        assert!((si.com[0] - 1.0).abs() < 1e-12, "com[0] = {}", si.com[0]);
        assert!(si.com[1].abs() < 1e-12);
    }
    #[test]
    fn test_spatial_inertia6_thin_rod_inertia() {
        let si = SpatialInertia6::thin_rod_x(3.0, 2.0);
        let expected = 3.0 * 4.0 / 3.0;
        assert!(
            (si.i_rot[1][1] - expected).abs() < 1e-10,
            "I_yy={}",
            si.i_rot[1][1]
        );
    }
    #[test]
    fn test_abi_backward_single_link() {
        let abi = abi_backward_pass(1, &[1.0], &[1.0], &[1.0]);
        assert_eq!(abi.len(), 1);
        assert!((abi[0] - 1.25).abs() < 1e-10, "ABI = {}", abi[0]);
    }
    #[test]
    fn test_abi_backward_accumulates_children() {
        let abi = abi_backward_pass(2, &[1.0, 0.5], &[1.0, 1.0], &[1.0, 1.0]);
        assert_eq!(abi.len(), 2);
        assert!((abi[1] - 0.75).abs() < 1e-10, "ABI[1]={}", abi[1]);
        assert!((abi[0] - 2.0).abs() < 1e-10, "ABI[0]={}", abi[0]);
    }
    #[test]
    fn test_abi_backward_empty() {
        let abi = abi_backward_pass(0, &[], &[], &[]);
        assert!(abi.is_empty());
    }
    #[test]
    fn test_coriolis_bias_zero_velocity() {
        let bias = coriolis_bias_forces(2, &[1.0, 1.0], &[1.0, 1.0], &[0.0, 0.0], &[0.0, 0.0]);
        for &b in &bias {
            assert!(b.abs() < 1e-12, "bias should be zero, got {b}");
        }
    }
    #[test]
    fn test_coriolis_bias_nonzero_for_rotating_chain() {
        let bias = coriolis_bias_forces(2, &[1.0, 1.0], &[1.0, 1.0], &[0.0, 0.0], &[1.0, 0.0]);
        let total: f64 = bias.iter().map(|&b| b.abs()).sum();
        assert!(
            total > 0.0,
            "nonzero velocity should produce Coriolis forces"
        );
    }
    #[test]
    fn test_coriolis_bias_empty() {
        let bias = coriolis_bias_forces(0, &[], &[], &[], &[]);
        assert!(bias.is_empty());
    }
    #[test]
    fn test_forward_dynamics_full_zero_torque_zero_gravity() {
        let qdd = forward_dynamics_full(
            2,
            &[1.0, 0.5],
            &[1.0, 0.5],
            &[1.0, 0.8],
            &[0.0, 0.0],
            &[0.0, 0.0],
            &[0.0, 0.0],
            [0.0, 0.0, 0.0],
        );
        assert_eq!(qdd.len(), 2);
        for &a in &qdd {
            assert!(a.abs() < 1e-10, "qdd should be zero, got {a}");
        }
    }
    #[test]
    fn test_forward_dynamics_full_positive_torque() {
        let qdd = forward_dynamics_full(
            1,
            &[1.0],
            &[1.0],
            &[1.0],
            &[0.0],
            &[0.0],
            &[5.0],
            [0.0, 0.0, 0.0],
        );
        assert_eq!(qdd.len(), 1);
        assert!(
            qdd[0] > 0.0,
            "positive torque should give positive acceleration"
        );
    }
    #[test]
    fn test_forward_dynamics_full_empty() {
        let qdd = forward_dynamics_full(0, &[], &[], &[], &[], &[], &[], [0.0; 3]);
        assert!(qdd.is_empty());
    }
    #[test]
    fn test_kinetic_energy_zero_velocity() {
        let ke = kinetic_energy(2, &[1.0, 0.5], &[0.0, 0.0], &[0.0, 0.0]);
        assert!(ke.abs() < 1e-12, "KE with zero velocity = {ke}");
    }
    #[test]
    fn test_kinetic_energy_single_link_known() {
        let ke = kinetic_energy(1, &[1.0], &[0.0], &[2.0]);
        assert!((ke - 2.0).abs() < 1e-10, "KE = {ke}");
    }
    #[test]
    fn test_kinetic_energy_positive() {
        let ke = kinetic_energy(2, &[1.0, 1.0], &[0.0, 0.0], &[1.0, 1.0]);
        assert!(ke > 0.0, "nonzero velocity should give positive KE");
    }
    #[test]
    fn test_potential_energy_zero_gravity() {
        let pe = potential_energy_chain(2, &[1.0, 1.0], &[1.0, 1.0], &[0.0, 0.0], 0.0, [0.0, 0.0]);
        assert!(pe.abs() < 1e-12, "PE with zero gravity = {pe}");
    }
    #[test]
    fn test_potential_energy_vertical_link() {
        use std::f64::consts::FRAC_PI_2;
        let pe = potential_energy_chain(1, &[1.0], &[1.0], &[FRAC_PI_2], -9.81, [0.0, 0.0]);
        let expected = 1.0 * (-9.81) * 0.5;
        assert!(
            (pe - expected).abs() < 1e-8,
            "PE = {pe}, expected {expected}"
        );
    }
    #[test]
    fn test_potential_energy_empty_chain() {
        let pe = potential_energy_chain(0, &[], &[], &[], -9.81, [0.0, 0.0]);
        assert_eq!(pe, 0.0);
    }
    #[test]
    fn test_tip_velocities_zero_qd() {
        let _chain = {
            let mut c = ArticulatedChain2D::new();
            c.add_link(1.0, 1.0, None);
            c.add_link(1.0, 1.0, Some(0));
            c
        };
        let lengths = vec![1.0, 1.0];
        let q = vec![0.0, 0.0];
        let qd = vec![0.0, 0.0];
        let vels = link_tip_velocities(2, &lengths, &q, &qd, [0.0, 0.0]);
        assert_eq!(vels.len(), 2);
        for &v in vels.iter().flatten() {
            assert!(
                v.abs() < 1e-12,
                "zero velocity → zero tip velocity, got {v}"
            );
        }
    }
    #[test]
    fn test_tip_velocities_single_link_known() {
        let vels = link_tip_velocities(1, &[1.0], &[0.0], &[1.0], [0.0, 0.0]);
        assert_eq!(vels.len(), 1);
        assert!(vels[0][0].abs() < 1e-12, "vx = {}", vels[0][0]);
        assert!((vels[0][1] - 1.0).abs() < 1e-12, "vy = {}", vels[0][1]);
    }
    #[test]
    fn test_tip_velocities_empty() {
        let vels = link_tip_velocities(0, &[], &[], &[], [0.0, 0.0]);
        assert!(vels.is_empty());
    }
    #[test]
    fn test_enforce_joint_limits_no_violation() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        chain.links[0].joint_angle = 0.5;
        chain.links[0].joint_velocity = 1.0;
        let triggered = enforce_joint_limits(&mut chain.links, &[(-1.0, 1.0)]);
        assert!(!triggered, "no limit violation should not trigger");
        assert!((chain.links[0].joint_angle - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_enforce_joint_limits_lower_bound() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        chain.links[0].joint_angle = -2.0;
        chain.links[0].joint_velocity = -1.0;
        let triggered = enforce_joint_limits(&mut chain.links, &[(-1.0, 1.0)]);
        assert!(triggered);
        assert!(
            (chain.links[0].joint_angle - (-1.0)).abs() < 1e-12,
            "clamped to q_min"
        );
        assert!(
            chain.links[0].joint_velocity >= 0.0,
            "outward velocity zeroed"
        );
    }
    #[test]
    fn test_enforce_joint_limits_upper_bound() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        chain.links[0].joint_angle = 2.0;
        chain.links[0].joint_velocity = 1.0;
        let triggered = enforce_joint_limits(&mut chain.links, &[(-1.0, 1.0)]);
        assert!(triggered);
        assert!(
            (chain.links[0].joint_angle - 1.0).abs() < 1e-12,
            "clamped to q_max"
        );
        assert!(
            chain.links[0].joint_velocity <= 0.0,
            "outward velocity zeroed"
        );
    }
    #[test]
    fn test_workspace_sampling_count() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        let pts = sample_workspace(
            &chain,
            [0.0, 0.0],
            20,
            &[-std::f64::consts::PI],
            &[std::f64::consts::PI],
        );
        assert_eq!(pts.len(), 20);
    }
    #[test]
    fn test_workspace_sampling_all_finite() {
        let mut chain = ArticulatedChain2D::new();
        chain.add_link(1.0, 1.0, None);
        chain.add_link(1.0, 0.8, Some(0));
        let pts = sample_workspace(
            &chain,
            [0.0, 0.0],
            50,
            &[-std::f64::consts::PI, -std::f64::consts::PI],
            &[std::f64::consts::PI, std::f64::consts::PI],
        );
        for pt in &pts {
            assert!(
                pt[0].is_finite() && pt[1].is_finite(),
                "workspace point must be finite"
            );
        }
    }
    #[test]
    fn test_workspace_sampling_empty_chain() {
        let chain = ArticulatedChain2D::new();
        let pts = sample_workspace(&chain, [0.0, 0.0], 10, &[], &[]);
        assert!(pts.is_empty());
    }
    #[test]
    fn test_cubic_spline_start_equals_q0() {
        assert!((cubic_spline_position(1.0, 5.0, 0.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_cubic_spline_end_equals_q1() {
        assert!((cubic_spline_position(1.0, 5.0, 1.0) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_cubic_spline_midpoint() {
        let mid = cubic_spline_position(0.0, 2.0, 0.5);
        assert!((mid - 1.0).abs() < 1e-12, "midpoint = {mid}");
    }
    #[test]
    fn test_cubic_spline_velocity_zero_at_endpoints() {
        let v0 = cubic_spline_velocity(0.0, 1.0, 0.0, 1.0);
        let v1 = cubic_spline_velocity(0.0, 1.0, 1.0, 1.0);
        assert!(v0.abs() < 1e-12, "velocity at t=0 should be 0, got {v0}");
        assert!(v1.abs() < 1e-12, "velocity at t=1 should be 0, got {v1}");
    }
    #[test]
    fn test_cubic_spline_velocity_peak_at_midpoint() {
        let v_mid = cubic_spline_velocity(0.0, 1.0, 0.5, 1.0);
        let v_quarter = cubic_spline_velocity(0.0, 1.0, 0.25, 1.0);
        assert!(v_mid > v_quarter, "velocity should peak at midpoint");
    }
    #[test]
    fn test_sample_cubic_trajectory_shape() {
        let q_start = vec![0.0, 0.0];
        let q_end = vec![1.0, 0.5];
        let traj = sample_cubic_trajectory(&q_start, &q_end, 1.0, 10);
        assert_eq!(traj.len(), 10);
        for config in &traj {
            assert_eq!(config.len(), 2);
        }
    }
    #[test]
    fn test_sample_cubic_trajectory_endpoints() {
        let q_start = vec![0.0, 1.0];
        let q_end = vec![2.0, 3.0];
        let traj = sample_cubic_trajectory(&q_start, &q_end, 1.0, 5);
        assert!((traj[0][0] - 0.0).abs() < 1e-12, "start[0]={}", traj[0][0]);
        assert!((traj[0][1] - 1.0).abs() < 1e-12, "start[1]={}", traj[0][1]);
        assert!((traj[4][0] - 2.0).abs() < 1e-12, "end[0]={}", traj[4][0]);
        assert!((traj[4][1] - 3.0).abs() < 1e-12, "end[1]={}", traj[4][1]);
    }
    #[test]
    fn test_sample_cubic_trajectory_single_step() {
        let traj = sample_cubic_trajectory(&[1.0], &[5.0], 1.0, 1);
        assert_eq!(traj.len(), 1);
        assert!((traj[0][0] - 1.0).abs() < 1e-12);
    }
}
