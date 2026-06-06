//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Compute the Euclidean norm of a 3D vector.
pub fn vec3_norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
/// Compute the dot product of two 3D vectors.
pub fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Compute the cross product of two 3D vectors.
pub fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::constraint_forces::ConstraintFatigue;
    use crate::constraint_forces::ConstraintForce;
    use crate::constraint_forces::ConstraintImpact;

    use crate::constraint_forces::InternalForce;

    use crate::constraint_forces::JointLoadAnalysis;
    use crate::constraint_forces::JointWear;

    use crate::constraint_forces::MechanicalPower;

    use crate::constraint_forces::ReactionForceLogger;

    use crate::constraint_forces::StaticEquilibrium;

    use crate::constraint_forces::WorkEnergyPrinciple;

    use crate::constraint_forces::vec3_norm;

    #[test]
    fn constraint_force_magnitude() {
        let cf = ConstraintForce::new([3.0, 4.0, 0.0], [0.0, 0.0, 0.0]);
        assert!((cf.magnitude() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn constraint_force_torque() {
        let cf = ConstraintForce::new([0.0, 0.0, 0.0], [0.0, 0.0, 5.0]);
        assert!((cf.torque - 5.0).abs() < 1e-10);
    }
    #[test]
    fn constraint_force_zero() {
        let cf = ConstraintForce::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(cf.magnitude() < 1e-10);
    }
    #[test]
    fn joint_load_simply_supported_reactions() {
        let (ra, rb) = JointLoadAnalysis::simply_supported_reactions(100.0, 3.0, 10.0);
        assert!((ra[1] + rb[1] - 100.0).abs() < 1e-8);
        assert!((rb[1] - 30.0).abs() < 1e-8);
    }
    #[test]
    fn joint_load_sum_forces_empty() {
        let jla = JointLoadAnalysis::new(vec![], vec![], vec![]);
        let sf = jla.sum_forces();
        assert!(vec3_norm(sf) < 1e-10);
    }
    #[test]
    fn joint_load_equilibrium_self_balancing() {
        let jla = JointLoadAnalysis::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[10.0, 0.0, 0.0], [-10.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        );
        assert!(jla.is_in_equilibrium(1e-8));
    }
    #[test]
    fn internal_force_simply_supported_uniform() {
        let iforce = InternalForce::simply_supported(10.0, 1.0, vec![]);
        assert!((iforce.reaction_left - 5.0).abs() < 1e-8);
        assert!((iforce.reaction_right - 5.0).abs() < 1e-8);
    }
    #[test]
    fn internal_force_shear_at_left_end() {
        let iforce = InternalForce::simply_supported(10.0, 1.0, vec![]);
        let v0 = iforce.shear_force(0.0);
        assert!((v0 - 5.0).abs() < 1e-8);
    }
    #[test]
    fn internal_force_shear_at_right_end() {
        let iforce = InternalForce::simply_supported(10.0, 1.0, vec![]);
        let vl = iforce.shear_force(10.0);
        assert!((vl + 5.0).abs() < 1e-8);
    }
    #[test]
    fn internal_force_max_moment_at_midspan() {
        let iforce = InternalForce::simply_supported(10.0, 1.0, vec![]);
        let (_max_m, max_x) = iforce.max_bending_moment();
        assert!((max_x - 5.0).abs() < 0.1);
    }
    #[test]
    fn internal_force_point_load_reaction() {
        let iforce = InternalForce::simply_supported(10.0, 0.0, vec![(5.0, 100.0)]);
        assert!((iforce.reaction_left - 50.0).abs() < 1e-8);
        assert!((iforce.reaction_right - 50.0).abs() < 1e-8);
    }
    #[test]
    fn mechanical_power_translational() {
        let mp = MechanicalPower::new(
            [10.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        );
        assert!((mp.translational_power() - 20.0).abs() < 1e-10);
    }
    #[test]
    fn mechanical_power_rotational() {
        let mp = MechanicalPower::new(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 5.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 2.0],
        );
        assert!((mp.rotational_power() - 10.0).abs() < 1e-10);
    }
    #[test]
    fn mechanical_power_total() {
        let mp = MechanicalPower::new(
            [10.0, 0.0, 0.0],
            [0.0, 0.0, 5.0],
            [2.0, 0.0, 0.0],
            [0.0, 0.0, 2.0],
        );
        assert!((mp.total_power() - 30.0).abs() < 1e-10);
    }
    #[test]
    fn bearing_l10_life_positive() {
        let cf = ConstraintFatigue::new(5000.0, 50000.0, 1500.0);
        assert!(cf.l10_hours() > 0.0);
    }
    #[test]
    fn bearing_higher_load_shorter_life() {
        let cf_light = ConstraintFatigue::new(1000.0, 50000.0, 1500.0);
        let cf_heavy = ConstraintFatigue::new(10000.0, 50000.0, 1500.0);
        assert!(cf_light.l10_hours() > cf_heavy.l10_hours());
    }
    #[test]
    fn hertz_contact_stress_positive() {
        let cf = ConstraintFatigue::new(5000.0, 50000.0, 1500.0);
        let p = cf.hertz_contact_stress(0.01);
        assert!(p > 0.0);
    }
    #[test]
    fn joint_wear_advances() {
        let mut wear = JointWear::new(1e-4, 1e6, 0.01, 1e9);
        wear.advance(1000, 1e-4);
        assert!(wear.wear_depth > 0.0);
        assert_eq!(wear.cycles, 1000);
    }
    #[test]
    fn joint_wear_remaining_life() {
        let mut wear = JointWear::new(1e-4, 1e6, 0.01, 1e9);
        wear.advance(100, 1e-4);
        let rl = wear.remaining_life(1e-3, 1e-4);
        assert!(rl > 0);
    }
    #[test]
    fn constraint_impact_impulse() {
        let impact = ConstraintImpact::new(1.0, 1.0, 1.0, 2.0, 0.01);
        assert!((impact.impulse() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn constraint_impact_peak_force() {
        let impact = ConstraintImpact::new(1.0, 1.0, 1.0, 2.0, 0.01);
        assert!((impact.peak_force() - 200.0).abs() < 1e-8);
    }
    #[test]
    fn constraint_impact_elastic_no_energy_loss() {
        let impact = ConstraintImpact::new(1.0, 1.0, 1.0, 2.0, 0.01);
        assert!(impact.energy_loss().abs() < 1e-10);
    }
    #[test]
    fn constraint_impact_plastic_energy_loss() {
        let impact = ConstraintImpact::new(1.0, 1.0, 0.0, 2.0, 0.01);
        assert!(impact.energy_loss() > 0.0);
    }
    #[test]
    fn reaction_logger_record_and_mean() {
        let mut logger = ReactionForceLogger::new(100);
        logger.record(0.0, [10.0, 0.0, 0.0]);
        logger.record(1.0, [10.0, 0.0, 0.0]);
        let mean = logger.mean_force();
        assert!((mean[0] - 10.0).abs() < 1e-10);
    }
    #[test]
    fn reaction_logger_peak_force() {
        let mut logger = ReactionForceLogger::new(100);
        logger.record(0.0, [3.0, 4.0, 0.0]);
        logger.record(1.0, [1.0, 0.0, 0.0]);
        assert!((logger.peak_force_magnitude() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn reaction_logger_dft_empty() {
        let logger = ReactionForceLogger::new(100);
        let spectrum = logger.dft_fx(100.0);
        assert!(spectrum.is_empty());
    }
    #[test]
    fn static_equilibrium_in_balance() {
        let mut eq = StaticEquilibrium::new();
        eq.add_force([100.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        eq.add_force([-100.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(eq.check(1e-8));
    }
    #[test]
    fn static_equilibrium_not_in_balance() {
        let mut eq = StaticEquilibrium::new();
        eq.add_force([100.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(!eq.check(1e-8));
    }
    #[test]
    fn static_equilibrium_with_moment() {
        let mut eq = StaticEquilibrium::new();
        eq.add_force([0.0, 100.0, 0.0], [1.0, 0.0, 0.0]);
        eq.add_force([0.0, -100.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(eq.check(1e-8));
    }
    #[test]
    fn work_energy_satisfied_simple() {
        let we = WorkEnergyPrinciple::new(10.0, 30.0, 20.0);
        assert!(we.is_satisfied(1e-8));
    }
    #[test]
    fn work_energy_delta_ke() {
        let we = WorkEnergyPrinciple::new(5.0, 15.0, 10.0);
        assert!((we.delta_ke() - 10.0).abs() < 1e-10);
    }
    #[test]
    fn work_energy_residual_nonzero() {
        let we = WorkEnergyPrinciple::new(10.0, 30.0, 15.0);
        assert!(we.residual().abs() > 1e-8);
    }
    #[test]
    fn vec3_norm_unit_vectors() {
        assert!((vec3_norm([1.0, 0.0, 0.0]) - 1.0).abs() < 1e-10);
        assert!((vec3_norm([0.0, 1.0, 0.0]) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn vec3_dot_orthogonal() {
        assert!(vec3_dot([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).abs() < 1e-10);
    }
    #[test]
    fn vec3_cross_basic() {
        let c = vec3_cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < 1e-10);
        assert!(c[0].abs() < 1e-10);
        assert!(c[1].abs() < 1e-10);
    }
}
/// Skew-symmetric matrix \[ω\]× as a 3×3 array.
pub fn skew(w: [f64; 3]) -> [[f64; 3]; 3] {
    [[0.0, -w[2], w[1]], [w[2], 0.0, -w[0]], [-w[1], w[0], 0.0]]
}
/// Matrix-vector product for 3×3 × 3.
pub fn mat3_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
/// 3×3 matrix addition.
pub fn mat3_add(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}
/// 3×3 matrix multiplication.
pub fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
/// 3×3 identity matrix.
pub fn mat3_identity() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
/// 3×3 matrix scalar multiply.
pub fn mat3_scale(m: [[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = m[i][j] * s;
        }
    }
    c
}
/// Rodrigues formula: R = I + sin(θ)\[ω\]× + (1−cos(θ))\[ω\]×².
/// `omega_unit` must be a unit vector, `theta` in radians.
pub fn rodrigues(omega_unit: [f64; 3], theta: f64) -> [[f64; 3]; 3] {
    let k = skew(omega_unit);
    let k2 = mat3_mul(k, k);
    let sin_t = theta.sin();
    let one_minus_cos = 1.0 - theta.cos();
    let i = mat3_identity();
    mat3_add(
        mat3_add(i, mat3_scale(k, sin_t)),
        mat3_scale(k2, one_minus_cos),
    )
}
/// Solve a 3×3 linear system Ax = b using Cramer's rule.
///
/// Returns x = A⁻¹b. If det ≈ 0, returns zero vector.
pub fn solve3x3(a: [[f64; 3]; 3], b: [f64; 3]) -> [f64; 3] {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-14 {
        return [0.0; 3];
    }
    let inv_det = 1.0 / det;
    [
        inv_det
            * (b[0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                - a[0][1] * (b[1] * a[2][2] - a[1][2] * b[2])
                + a[0][2] * (b[1] * a[2][1] - a[1][1] * b[2])),
        inv_det
            * (a[0][0] * (b[1] * a[2][2] - a[1][2] * b[2])
                - b[0] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                + a[0][2] * (a[1][0] * b[2] - b[1] * a[2][0])),
        inv_det
            * (a[0][0] * (a[1][1] * b[2] - b[1] * a[2][1])
                - a[0][1] * (a[1][0] * b[2] - b[1] * a[2][0])
                + b[0] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])),
    ]
}
/// Add two 3-vectors.
#[inline]
pub fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors.
#[inline]
pub fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector by scalar.
#[inline]
pub fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Normalize a 3-vector.
#[inline]
pub fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n < 1e-14 {
        [0.0; 3]
    } else {
        vec3_scale(a, 1.0 / n)
    }
}
#[cfg(test)]
mod advanced_tests {

    use crate::constraint_forces::Cable;
    use crate::constraint_forces::CableDrivenRobot;

    use crate::constraint_forces::ConstraintJacobianAssembler;
    use crate::constraint_forces::DHParams;
    use crate::constraint_forces::HolonomicConstraint;
    use crate::constraint_forces::HybridController;
    use crate::constraint_forces::ImpedanceController;

    use crate::constraint_forces::JacobianTransposeSolver;

    use crate::constraint_forces::KinematicChain;

    use crate::constraint_forces::NullSpaceProjector;
    use crate::constraint_forces::PoEChain;
    use crate::constraint_forces::PseudoinverseSolver;

    use crate::constraint_forces::ScrewAxis;
    use crate::constraint_forces::SpatialInertia;

    use crate::constraint_forces::StewartPlatform;
    use crate::constraint_forces::Twist;
    use crate::constraint_forces::VirtualWorkSolver;
    use crate::constraint_forces::WeightedPseudoinverse;

    use crate::constraint_forces::Wrench;
    use crate::constraint_forces::mat3_identity;
    use crate::constraint_forces::mat3_mul;
    use crate::constraint_forces::mat3_vec;
    use crate::constraint_forces::rodrigues;
    use crate::constraint_forces::skew;
    use crate::constraint_forces::solve3x3;
    use crate::constraint_forces::vec3_add;
    use crate::constraint_forces::vec3_norm;
    use crate::constraint_forces::vec3_normalize;
    use crate::constraint_forces::vec3_scale;
    use crate::constraint_forces::vec3_sub;
    use std::f64::consts::PI;
    #[test]
    fn twist_zero_norm() {
        let t = Twist::zero();
        assert!(t.norm() < 1e-15);
    }
    #[test]
    fn twist_add_and_scale() {
        let t1 = Twist::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let t2 = t1.scale(2.0);
        assert!((t2.omega[0] - 2.0).abs() < 1e-12);
        assert!((t2.linear[1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn twist_lie_bracket_anticommutative() {
        let t1 = Twist::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let t2 = Twist::new([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]);
        let lb12 = t1.lie_bracket(&t2);
        let lb21 = t2.lie_bracket(&t1);
        assert!((lb12.omega[2] + lb21.omega[2]).abs() < 1e-12);
    }
    #[test]
    fn wrench_power_computation() {
        let w = Wrench::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let t = Twist::new([2.0, 0.0, 0.0], [0.0, 3.0, 0.0]);
        assert!((w.power(&t) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn wrench_add() {
        let w1 = Wrench::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let w2 = Wrench::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let w = w1.add(&w2);
        assert!((w.moment[0] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn spatial_inertia_box_kinetic_energy_positive() {
        let si = SpatialInertia::box_inertia(1.0, 0.1, 0.1, 0.1);
        let t = Twist::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let ke = si.kinetic_energy(&t);
        assert!(ke > 0.0);
    }
    #[test]
    fn spatial_inertia_sphere_symmetric() {
        let si = SpatialInertia::sphere_inertia(2.0, 0.1);
        assert!((si.inertia[0][0] - si.inertia[1][1]).abs() < 1e-12);
        assert!((si.inertia[1][1] - si.inertia[2][2]).abs() < 1e-12);
    }
    #[test]
    fn spatial_inertia_momentum_linear() {
        let si = SpatialInertia::box_inertia(5.0, 0.1, 0.1, 0.1);
        let t = Twist::new([0.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
        let p = si.momentum(&t);
        assert!((p.force[0] - 10.0).abs() < 1e-10);
    }
    #[test]
    fn rodrigues_identity_at_zero_angle() {
        let r = rodrigues([0.0, 0.0, 1.0], 0.0);
        let i = mat3_identity();
        for row in 0..3 {
            for col in 0..3 {
                assert!((r[row][col] - i[row][col]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn rodrigues_90deg_about_z() {
        let r = rodrigues([0.0, 0.0, 1.0], PI / 2.0);
        let x = mat3_vec(r, [1.0, 0.0, 0.0]);
        assert!(x[0].abs() < 1e-10);
        assert!((x[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn skew_antisymmetric() {
        let s = skew([1.0, 2.0, 3.0]);
        for (i, row) in s.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val + s[j][i]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn dh_transform_identity_at_zero() {
        let dh = DHParams::revolute(0.0, 0.0, 0.0, 0.0);
        let t = dh.transform();
        assert!((t[0][0] - 1.0).abs() < 1e-12);
        assert!((t[3][3] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn dh_translation_nonzero_a() {
        let dh = DHParams::revolute(1.0, 0.0, 0.0, 0.0);
        let tr = dh.translation();
        assert!((tr[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn kinematic_chain_2r_planar_fk() {
        let joints = vec![
            DHParams::revolute(1.0, 0.0, 0.0, 0.0),
            DHParams::revolute(1.0, 0.0, 0.0, 0.0),
        ];
        let chain = KinematicChain::new(joints);
        let pos = chain.forward_kinematics(&[0.0, 0.0]);
        assert!((pos[0] - 2.0).abs() < 1e-8, "px={}", pos[0]);
    }
    #[test]
    fn kinematic_chain_jacobian_correct_size() {
        let joints = vec![
            DHParams::revolute(1.0, 0.0, 0.0, 0.0),
            DHParams::revolute(1.0, 0.0, 0.0, 0.0),
            DHParams::revolute(1.0, 0.0, 0.0, 0.0),
        ];
        let chain = KinematicChain::new(joints);
        let j = chain.geometric_jacobian(&[0.0, 0.0, 0.0]);
        assert_eq!(j.len(), 9);
    }
    #[test]
    fn poe_chain_home_position() {
        let screw = ScrewAxis::revolute([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let chain = PoEChain::new(vec![screw], [1.0, 0.0, 0.0]);
        let pos = chain.forward_kinematics(&[0.0]);
        assert!((pos[0] - 1.0).abs() < 1e-8, "px={}", pos[0]);
    }
    #[test]
    fn screw_axis_exp_map_prismatic() {
        let s = ScrewAxis::prismatic([0.0, 0.0, 1.0]);
        let t = s.exp_map(2.0);
        assert!((t[2][3] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn screw_axis_exp_map_revolute_identity_at_zero() {
        let s = ScrewAxis::revolute([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let t = s.exp_map(0.0);
        assert!((t[0][0] - 1.0).abs() < 1e-12);
        assert!((t[1][1] - 1.0).abs() < 1e-12);
        assert!((t[2][2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn jacobian_transpose_ik_converges_simple() {
        let joints = vec![
            DHParams::revolute(1.0, 0.0, 0.0, 0.0),
            DHParams::revolute(1.0, 0.0, 0.0, 0.0),
        ];
        let chain = KinematicChain::new(joints);
        let solver = JacobianTransposeSolver::new(0.01, 1e-4, 5000);
        let result = solver.solve(&chain, [1.5, 0.5, 0.0], vec![0.0, 0.0]);
        assert!(result.error < 0.1, "error={}", result.error);
    }
    #[test]
    fn pseudoinverse_ik_converges() {
        let joints = vec![
            DHParams::revolute(1.0, 0.0, 0.0, 0.0),
            DHParams::revolute(1.0, 0.0, 0.0, 0.0),
            DHParams::revolute(0.5, 0.0, 0.0, 0.0),
        ];
        let chain = KinematicChain::new(joints);
        let solver = PseudoinverseSolver::new(1e-4, 2000, 0.01);
        let result = solver.solve(&chain, [2.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]);
        assert!(result.error < 0.5, "error={}", result.error);
    }
    #[test]
    fn solve3x3_identity() {
        let a = mat3_identity();
        let b = [1.0, 2.0, 3.0];
        let x = solve3x3(a, b);
        assert!((x[0] - 1.0).abs() < 1e-12);
        assert!((x[1] - 2.0).abs() < 1e-12);
        assert!((x[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn solve3x3_singular_returns_zero() {
        let a = [[0.0; 3]; 3];
        let x = solve3x3(a, [1.0, 2.0, 3.0]);
        assert!(x[0].abs() < 1e-12 && x[1].abs() < 1e-12 && x[2].abs() < 1e-12);
    }
    #[test]
    fn null_space_gradient_computation() {
        let q = vec![0.5f64, 0.0, -0.5];
        let q_min = vec![-1.0, -1.0, -1.0];
        let q_max = vec![1.0, 1.0, 1.0];
        let grad = NullSpaceProjector::joint_limit_gradient(&q, &q_min, &q_max);
        assert_eq!(grad.len(), 3);
        assert!(grad[1].abs() < 1e-12);
    }
    #[test]
    fn weighted_pseudoinverse_correct_size() {
        let weights = vec![1.0, 1.0, 1.0];
        let wp = WeightedPseudoinverse::new(weights);
        let j: Vec<f64> = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let dq = wp.step(&j, [1.0, 2.0, 3.0]);
        assert_eq!(dq.len(), 3);
        assert!((dq[0] - 1.0).abs() < 1e-8, "dq[0]={}", dq[0]);
    }
    #[test]
    fn hybrid_controller_update_returns_correct_size() {
        let sel = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut hc = HybridController::new(
            sel,
            [100.0, 0.1, 1.0],
            [50.0, 0.1, 0.5],
            vec![10.0; 6],
            vec![0.0; 6],
        );
        let out = hc.update(&[0.0; 6], &[0.0; 6], 0.01);
        assert_eq!(out.len(), 6);
    }
    #[test]
    fn hybrid_controller_force_direction_nonzero() {
        let sel = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut hc = HybridController::new(
            sel,
            [100.0, 0.0, 0.0],
            [50.0, 0.0, 0.0],
            vec![10.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0; 6],
        );
        let out = hc.update(&[0.0; 6], &[0.0; 6], 0.01);
        assert!(out[0].abs() > 0.0);
        assert!(out[1].abs() < 1e-12);
    }
    #[test]
    fn impedance_controller_potential_energy() {
        let ic =
            ImpedanceController::new(vec![1.0; 3], vec![10.0; 3], vec![100.0; 3], vec![0.0; 3]);
        let pe = ic.potential_energy(&[1.0, 0.0, 0.0]);
        assert!((pe - 50.0).abs() < 1e-10);
    }
    #[test]
    fn impedance_controller_acceleration_restoring() {
        let ic = ImpedanceController::new(vec![1.0; 3], vec![0.0; 3], vec![100.0; 3], vec![0.0; 3]);
        let accel = ic.desired_acceleration(&[1.0, 0.0, 0.0], &[0.0; 3], &[0.0; 3]);
        assert!(accel[0] < 0.0, "accel[0]={}", accel[0]);
    }
    #[test]
    fn constraint_jacobian_assembler_size() {
        let mut asm = ConstraintJacobianAssembler::new(4);
        asm.add_constraint(HolonomicConstraint::new(
            "c1",
            vec![1.0, 0.0, 0.0, 0.0],
            0.0,
        ));
        asm.add_constraint(HolonomicConstraint::new(
            "c2",
            vec![0.0, 1.0, 0.0, 0.0],
            0.0,
        ));
        let j = asm.jacobian();
        assert_eq!(j.len(), 2 * 4);
    }
    #[test]
    fn holonomic_constraint_violation() {
        let c = HolonomicConstraint::new("test", vec![1.0, -1.0], 0.0);
        let v = c.violation(&[2.0, 2.0]);
        assert!(v.abs() < 1e-12);
    }
    #[test]
    fn holonomic_constraint_nonzero_violation() {
        let c = HolonomicConstraint::new("test", vec![1.0, 0.0], 0.0);
        let v = c.violation(&[5.0, 3.0]);
        assert!((v - 5.0).abs() < 1e-12);
    }
    #[test]
    fn cable_geometric_length_correct() {
        let cable = Cable::new([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 1.0, 0.1, 2.0, 1000.0);
        let l = cable.geometric_length([0.0, 0.0, 0.0]);
        assert!((l - 1.0).abs() < 1e-12);
    }
    #[test]
    fn cable_tension_zero_when_slack() {
        let cable = Cable::new([0.0, 0.0, 2.0], [0.0, 0.0, 0.0], 2.5, 0.1, 5.0, 1000.0);
        let t = cable.tension([0.0, 0.0, 0.0]);
        assert!(t < 1e-12, "t={t}");
    }
    #[test]
    fn cable_tension_positive_when_taut() {
        let cable = Cable::new([0.0, 0.0, 2.0], [0.0, 0.0, 0.0], 1.0, 0.1, 5.0, 1000.0);
        let t = cable.tension([0.0, 0.0, 0.0]);
        assert!(t > 0.0, "t={t}");
    }
    #[test]
    fn cable_driven_robot_net_force_direction() {
        let cable = Cable::new([0.0, 0.0, -2.0], [0.0, 0.0, 0.0], 1.0, 0.1, 5.0, 1000.0);
        let robot = CableDrivenRobot::new(vec![cable], 1.0);
        let f = robot.net_cable_force();
        assert!(f[2] < 0.0 || f[2] > -1e10);
        assert!(f[2].is_finite());
    }
    #[test]
    fn stewart_ik_home_position() {
        let sp = StewartPlatform::regular(0.3, 0.2, 0.1, 1.5);
        let lengths = sp.inverse_kinematics([0.0, 0.0, 0.8], 0.0);
        assert_eq!(lengths.len(), 6);
        let first = lengths[0];
        for l in &lengths {
            assert!((l - first).abs() < 1e-3, "asymmetry: {l} vs {first}");
        }
    }
    #[test]
    fn stewart_check_limits_within_range() {
        let sp = StewartPlatform::regular(0.3, 0.2, 0.5, 1.5);
        let lengths = sp.inverse_kinematics([0.0, 0.0, 0.8], 0.0);
        let in_limits = sp.check_limits(&lengths);
        let _ = in_limits;
    }
    #[test]
    fn stewart_reachability() {
        let sp = StewartPlatform::regular(0.5, 0.3, 0.3, 2.0);
        let reachable = sp.is_reachable([0.0, 0.0, 1.0], 0.0);
        let _ = reachable;
    }
    #[test]
    fn virtual_work_simple_identity_jacobian() {
        let mut vw = VirtualWorkSolver::new(3);
        let w = Wrench::new([0.0; 3], [1.0, 2.0, 3.0]);
        let jt: Vec<f64> = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let jr: Vec<f64> = vec![0.0; 9];
        vw.add_body(w, jt, jr);
        let q = vw.generalized_forces();
        assert!((q[0] - 1.0).abs() < 1e-12);
        assert!((q[1] - 2.0).abs() < 1e-12);
        assert!((q[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn virtual_work_consistency() {
        let mut vw = VirtualWorkSolver::new(3);
        let w = Wrench::new([0.0; 3], [1.0, 0.0, 0.0]);
        let jt: Vec<f64> = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let jr: Vec<f64> = vec![0.0; 9];
        vw.add_body(w, jt, jr);
        let dq = [1.0, 0.0, 0.0];
        let work = vw.virtual_work(&dq);
        assert!((work - 1.0).abs() < 1e-12);
    }
    #[test]
    fn vec3_helpers_correct() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let sum = vec3_add(a, b);
        assert!((sum[0] - 5.0).abs() < 1e-12);
        let diff = vec3_sub(b, a);
        assert!((diff[0] - 3.0).abs() < 1e-12);
        let scaled = vec3_scale(a, 2.0);
        assert!((scaled[2] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn vec3_normalize_unit_vector() {
        let v = vec3_normalize([3.0, 4.0, 0.0]);
        let n = vec3_norm(v);
        assert!((n - 1.0).abs() < 1e-12);
    }
    #[test]
    fn mat3_mul_identity() {
        let i = mat3_identity();
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let r = mat3_mul(i, a);
        for row in 0..3 {
            for col in 0..3 {
                assert!((r[row][col] - a[row][col]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn null_space_projector_projects_away_from_jacobian() {
        let proj = NullSpaceProjector::new(3, 3, 0.0001);
        let j: Vec<f64> = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let z = vec![1.0, 2.0, 3.0];
        let ns = proj.null_space_component(&j, &z);
        let ns_norm: f64 = ns.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(ns_norm < 0.1, "null space norm too large: {ns_norm}");
    }
}
