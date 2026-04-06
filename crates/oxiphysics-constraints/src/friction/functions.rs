//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
pub(super) fn normalize(a: [f64; 3]) -> [f64; 3] {
    let len = length(a);
    if len < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        scale(a, 1.0 / len)
    }
}
#[allow(dead_code)]
pub(super) fn negate(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}
/// Build an arbitrary tangent vector perpendicular to `n` (unit normal).
pub(super) fn build_tangent(n: [f64; 3]) -> [f64; 3] {
    if n[0].abs() < 0.9 {
        normalize(cross(n, [1.0, 0.0, 0.0]))
    } else {
        normalize(cross(n, [0.0, 1.0, 0.0]))
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::AnisotropicFriction;

    use crate::CoulombFriction;
    use crate::FrictionCacheEntry;
    use crate::FrictionCone;
    use crate::FrictionConeLinearization;
    use crate::FrictionContactCache;
    use crate::FrictionDissipation;
    use crate::FrictionForceAccumulator;
    use crate::FrictionMaterial;
    use crate::FrictionRegularizer;
    use crate::FrictionWarmstart;
    use crate::PairwiseFrictionSolver;

    use crate::RollingResistance;

    use crate::StaticFrictionConstraint;

    use crate::TangentBasisCache;

    use crate::VelocityDependentFriction;
    use crate::ViscousFriction;
    #[test]
    fn coulomb_static_no_slide() {
        let cf = CoulombFriction::new(0.5, 0.4);
        let tangential_force = 0.1;
        let normal_force = 1.0;
        assert!(!cf.is_sliding(tangential_force, normal_force));
    }
    #[test]
    fn coulomb_kinetic_opposes_velocity() {
        let cf = CoulombFriction::new(0.5, 0.4);
        let normal_force = 10.0;
        let tang_vel = [1.0, 0.0, 0.0];
        let f = cf.friction_force(tang_vel, normal_force);
        assert!(
            f[0] < 0.0,
            "friction force should oppose velocity, got {:?}",
            f
        );
        let expected_mag = cf.mu_kinetic * normal_force;
        let actual_mag = length(f);
        assert!(
            (actual_mag - expected_mag).abs() < 1e-9,
            "expected magnitude {}, got {}",
            expected_mag,
            actual_mag
        );
    }
    #[test]
    fn friction_cone_zero_tangential_inside() {
        let cone = FrictionCone::new(0.5, [0.0, 1.0, 0.0]);
        let impulse = [0.0, 5.0, 0.0];
        assert!(cone.inside_cone(impulse, 5.0));
    }
    #[test]
    fn friction_cone_project_clamps() {
        let cone = FrictionCone::new(0.3, [0.0, 1.0, 0.0]);
        let normal_impulse = 10.0;
        let impulse = [50.0, 0.0, 0.0];
        let projected = cone.project_impulse(impulse, normal_impulse);
        assert!(
            cone.inside_cone(projected, normal_impulse),
            "projected impulse should be inside cone, got {:?}",
            projected
        );
    }
    #[test]
    fn friction_material_combine_geometric_mean() {
        let a = FrictionMaterial::new(0.4, 0.3, 0.01, 0.8);
        let b = FrictionMaterial::new(0.9, 0.8, 0.02, 0.6);
        let c = FrictionMaterial::combine(&a, &b);
        let expected_static = (0.4_f64 * 0.9).sqrt();
        let expected_kinetic = (0.3_f64 * 0.8).sqrt();
        let expected_restitution = 0.6_f64;
        assert!((c.static_friction - expected_static).abs() < 1e-9);
        assert!((c.kinetic_friction - expected_kinetic).abs() < 1e-9);
        assert!((c.restitution - expected_restitution).abs() < 1e-9);
    }
    #[test]
    fn viscous_force_opposes_velocity() {
        let vf = ViscousFriction::new(2.0);
        let vel = [3.0, 0.0, 0.0];
        let f = vf.force(vel);
        assert!(f[0] < 0.0, "viscous force must oppose velocity");
        assert!((f[0] + 6.0).abs() < 1e-9);
    }
    #[test]
    fn anisotropic_friction_asymmetric() {
        let af = AnisotropicFriction::new(1.0, 0.1);
        let tangent1 = [1.0, 0.0, 0.0];
        let tangent2 = [0.0, 0.0, 1.0];
        let normal_impulse = 10.0;
        let imp_x = [9.0, 0.0, 0.0];
        let out_x = af.anisotropic_impulse(imp_x, tangent1, tangent2, normal_impulse);
        assert!(
            (out_x[0] - 9.0).abs() < 1e-6,
            "x-impulse should not be clamped"
        );
        let imp_z = [0.0, 0.0, 5.0];
        let out_z = af.anisotropic_impulse(imp_z, tangent1, tangent2, normal_impulse);
        assert!(
            out_z[2].abs() < 5.0,
            "z-impulse should be clamped, got {}",
            out_z[2]
        );
    }
    #[test]
    fn velocity_dependent_at_zero_speed() {
        let vdf = VelocityDependentFriction::new(0.6, 0.4, 0.01, 0.001);
        let mu = vdf.effective_mu(0.0);
        assert!(
            (mu - 0.6).abs() < 1e-12,
            "At zero speed, μ should be μ_s, got {mu}"
        );
    }
    #[test]
    fn velocity_dependent_stribeck_dip() {
        let vdf = VelocityDependentFriction::new(0.6, 0.4, 0.01, 0.0);
        let mu_fast = vdf.effective_mu(0.1);
        assert!(mu_fast < 0.6, "At moderate speed, μ should drop below μ_s");
        assert!(mu_fast > 0.39, "But should remain near μ_k");
    }
    #[test]
    fn velocity_dependent_viscous_increases() {
        let vdf = VelocityDependentFriction::new(0.6, 0.4, 0.01, 0.1);
        let mu_slow = vdf.effective_mu(0.01);
        let mu_fast = vdf.effective_mu(10.0);
        assert!(
            mu_fast > mu_slow,
            "With viscous term, μ should increase at high speed"
        );
    }
    #[test]
    fn velocity_dependent_friction_force() {
        let vdf = VelocityDependentFriction::new(0.6, 0.4, 0.01, 0.0);
        let vel = [1.0, 0.0, 0.0];
        let f = vdf.friction_force(vel, 10.0);
        assert!(f[0] < 0.0, "Force should oppose velocity");
    }
    #[test]
    fn warmstart_zero_factor() {
        let ws = FrictionWarmstart::new(0.0);
        assert_eq!(ws.warmstart_factor, 0.0);
    }
    #[test]
    fn warmstart_clamped_factor() {
        let ws = FrictionWarmstart::new(1.5);
        assert_eq!(ws.warmstart_factor, 1.0);
    }
    #[test]
    fn warmstart_get_impulse() {
        let mut ws = FrictionWarmstart::new(0.8);
        ws.resize(2);
        ws.update(0, [10.0, 0.0, 0.0]);
        let imp = ws.get_impulse(0);
        assert!((imp[0] - 8.0).abs() < 1e-12);
    }
    #[test]
    fn warmstart_clear() {
        let mut ws = FrictionWarmstart::new(0.8);
        ws.resize(2);
        ws.update(0, [10.0, 5.0, 3.0]);
        ws.clear();
        let imp = ws.get_impulse(0);
        assert!((length(imp)).abs() < 1e-12);
    }
    #[test]
    fn warmstart_out_of_bounds() {
        let ws = FrictionWarmstart::new(0.8);
        let imp = ws.get_impulse(100);
        assert_eq!(imp, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn cone_linearization_facets_count() {
        let cl = FrictionConeLinearization::new([0.0, 1.0, 0.0], 8);
        assert_eq!(cl.n_facets, 8);
        assert_eq!(cl.tangent_dirs.len(), 8);
    }
    #[test]
    fn cone_linearization_dirs_orthogonal_to_normal() {
        let normal = [0.0, 1.0, 0.0];
        let cl = FrictionConeLinearization::new(normal, 6);
        for dir in &cl.tangent_dirs {
            let d = dot(*dir, normal);
            assert!(
                d.abs() < 1e-10,
                "Tangent dir should be orthogonal to normal, dot={d}"
            );
        }
    }
    #[test]
    fn cone_linearization_dirs_unit_length() {
        let cl = FrictionConeLinearization::new([0.0, 0.0, 1.0], 4);
        for dir in &cl.tangent_dirs {
            let len = length(*dir);
            assert!(
                (len - 1.0).abs() < 1e-10,
                "Tangent dirs should be unit vectors"
            );
        }
    }
    #[test]
    fn cone_linearization_project_inside() {
        let cl = FrictionConeLinearization::new([0.0, 1.0, 0.0], 4);
        let small_impulse = [0.1, 0.0, 0.0];
        let (projected, _) = cl.project(small_impulse, 10.0);
        assert!((length(sub(projected, small_impulse))).abs() < 1e-10);
    }
    #[test]
    fn anisotropic_effective_mu_along_x() {
        let af = AnisotropicFriction::new(0.5, 0.3);
        let mu = af.effective_mu(0.0);
        assert!((mu - 0.5).abs() < 1e-10, "μ along x should be mu_x");
    }
    #[test]
    fn anisotropic_effective_mu_along_y() {
        let af = AnisotropicFriction::new(0.5, 0.3);
        let mu = af.effective_mu(std::f64::consts::FRAC_PI_2);
        assert!((mu - 0.3).abs() < 1e-10, "μ along y should be mu_y");
    }
    #[test]
    fn accumulator_single_impulse() {
        let mut acc = FrictionForceAccumulator::new(100.0);
        acc.accumulate([3.0, 4.0, 0.0]);
        let imp = acc.clamped_impulse();
        assert!((imp[0] - 3.0).abs() < 1e-12);
        assert!((imp[1] - 4.0).abs() < 1e-12);
        assert_eq!(acc.count, 1);
    }
    #[test]
    fn accumulator_multiple_impulses() {
        let mut acc = FrictionForceAccumulator::new(100.0);
        acc.accumulate([3.0, 0.0, 0.0]);
        acc.accumulate([0.0, 4.0, 0.0]);
        let imp = acc.clamped_impulse();
        assert!((imp[0] - 3.0).abs() < 1e-12);
        assert!((imp[1] - 4.0).abs() < 1e-12);
        assert_eq!(acc.count, 2);
    }
    #[test]
    fn accumulator_clamping() {
        let mut acc = FrictionForceAccumulator::new(5.0);
        acc.accumulate([30.0, 40.0, 0.0]);
        let imp = acc.clamped_impulse();
        let mag = length(imp);
        assert!(
            (mag - 5.0).abs() < 1e-10,
            "Should clamp to max_impulse, got {mag}"
        );
    }
    #[test]
    fn accumulator_reset() {
        let mut acc = FrictionForceAccumulator::new(100.0);
        acc.accumulate([10.0, 20.0, 30.0]);
        acc.reset();
        assert_eq!(acc.count, 0);
        assert!(length(acc.total_impulse) < 1e-12);
    }
    #[test]
    fn accumulator_saturation() {
        let mut acc = FrictionForceAccumulator::new(5.0);
        acc.accumulate([3.0, 4.0, 0.0]);
        assert!(acc.is_saturated());
        let mut acc2 = FrictionForceAccumulator::new(10.0);
        acc2.accumulate([1.0, 0.0, 0.0]);
        assert!(!acc2.is_saturated());
    }
    #[test]
    fn material_combine_arithmetic() {
        let a = FrictionMaterial::new(0.4, 0.3, 0.01, 0.8);
        let b = FrictionMaterial::new(0.6, 0.5, 0.03, 0.6);
        let c = FrictionMaterial::combine_arithmetic(&a, &b);
        assert!((c.static_friction - 0.5).abs() < 1e-12);
        assert!((c.kinetic_friction - 0.4).abs() < 1e-12);
        assert!((c.restitution - 0.8).abs() < 1e-12);
    }
    #[test]
    fn material_combine_max() {
        let a = FrictionMaterial::new(0.4, 0.3, 0.01, 0.8);
        let b = FrictionMaterial::new(0.6, 0.5, 0.03, 0.6);
        let c = FrictionMaterial::combine_max(&a, &b);
        assert!((c.static_friction - 0.6).abs() < 1e-12);
        assert!((c.kinetic_friction - 0.5).abs() < 1e-12);
        assert!((c.restitution - 0.6).abs() < 1e-12);
    }
    #[test]
    fn viscous_impulse_test() {
        let vf = ViscousFriction::new(2.0);
        let vel = [1.0, 0.0, 0.0];
        let imp = vf.impulse(vel, 0.01);
        assert!((imp[0] - (-0.02)).abs() < 1e-12);
    }
    #[test]
    fn rolling_resistance_torque_opposes_rotation() {
        let rr = RollingResistance::new(0.01);
        let omega = [0.0, 10.0, 0.0];
        let t = rr.torque(100.0, 0.5, omega);
        assert!(t[1] < 0.0, "Torque should oppose angular velocity");
        let expected_mag = 0.01 * 100.0 * 0.5;
        assert!((length(t) - expected_mag).abs() < 1e-9);
    }
    #[test]
    fn test_negate() {
        let v = negate([1.0, -2.0, 3.0]);
        assert_eq!(v, [-1.0, 2.0, -3.0]);
    }
    #[test]
    fn test_cross_product() {
        let c = cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[0]).abs() < 1e-12);
        assert!((c[1]).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn tangent_basis_cache_orthogonal() {
        let cache = TangentBasisCache::new([0.0, 1.0, 0.0]);
        assert!(dot(cache.tangent1, cache.normal).abs() < 1e-10);
        assert!(dot(cache.tangent2, cache.normal).abs() < 1e-10);
        assert!(dot(cache.tangent1, cache.tangent2).abs() < 1e-10);
    }
    #[test]
    fn tangent_basis_cache_unit_vectors() {
        let cache = TangentBasisCache::new([0.0, 0.0, 1.0]);
        assert!((length(cache.tangent1) - 1.0).abs() < 1e-10);
        assert!((length(cache.tangent2) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn tangent_basis_cache_freshness() {
        let cache = TangentBasisCache::new([0.0, 1.0, 0.0]);
        assert!(cache.is_fresh([0.0, 1.0, 0.0]));
        assert!(!cache.is_fresh([1.0, 0.0, 0.0]));
    }
    #[test]
    fn tangent_basis_cache_decompose_reconstruct() {
        let cache = TangentBasisCache::new([0.0, 1.0, 0.0]);
        let jt1 = 3.0;
        let jt2 = -2.0;
        let reconstructed = cache.reconstruct(jt1, jt2);
        let (decomp_jt1, decomp_jt2) = cache.decompose(reconstructed);
        assert!((decomp_jt1 - jt1).abs() < 1e-10, "jt1={decomp_jt1}");
        assert!((decomp_jt2 - jt2).abs() < 1e-10, "jt2={decomp_jt2}");
    }
    #[test]
    fn tangent_basis_cache_update_stale() {
        let mut cache = TangentBasisCache::new([0.0, 1.0, 0.0]);
        let old_t1 = cache.tangent1;
        cache.update([1.0, 0.0, 0.0]);
        assert!(dot(cache.tangent1, cache.normal).abs() < 1e-10);
        let diff = length(sub(cache.tangent1, old_t1));
        assert!(
            diff > 0.1,
            "Tangent basis should change when normal changes"
        );
    }
    #[test]
    fn friction_regularizer_zero_velocity_returns_zero() {
        let reg = FrictionRegularizer::new(0.5, 0.0);
        let result = reg.solve([0.0, 0.0, 0.0], 1.0, 10.0);
        assert!(length(result) < 1e-12);
    }
    #[test]
    fn friction_regularizer_opposes_velocity() {
        let reg = FrictionRegularizer::new(0.5, 0.0);
        let vel = [1.0, 0.0, 0.0];
        let result = reg.solve(vel, 1.0, 10.0);
        assert!(result[0] < 0.0, "regularizer should oppose velocity");
    }
    #[test]
    fn friction_regularizer_clamped_by_normal_impulse() {
        let reg = FrictionRegularizer::new(0.5, 0.0);
        let vel = [100.0, 0.0, 0.0];
        let normal_impulse = 1.0;
        let result = reg.solve(vel, 1.0, normal_impulse);
        assert!(
            length(result) <= reg.mu * normal_impulse + 1e-10,
            "impulse should be clamped by μ*Fn"
        );
    }
    #[test]
    fn friction_regularizer_compliance_reduces_impulse() {
        let reg_hard = FrictionRegularizer::new(0.5, 0.0);
        let reg_soft = FrictionRegularizer::new(0.5, 10.0);
        let vel = [1.0, 0.0, 0.0];
        let imp_hard = length(reg_hard.solve(vel, 1.0, 10.0));
        let imp_soft = length(reg_soft.solve(vel, 1.0, 10.0));
        assert!(
            imp_soft <= imp_hard + 1e-10,
            "compliance should reduce or equal the impulse"
        );
    }
    #[test]
    fn static_friction_constraint_zero_vel_no_impulse() {
        let mut sfc = StaticFrictionConstraint::new(
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            2.0,
            0.5,
            10.0,
        );
        let impulse = sfc.solve_velocity([0.0, 0.0, 0.0]);
        assert!(length(impulse) < 1e-12, "zero velocity → zero impulse");
    }
    #[test]
    fn static_friction_constraint_opposes_tangential_velocity() {
        let mut sfc = StaticFrictionConstraint::new(
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            2.0,
            0.5,
            10.0,
        );
        let rel_vel = [2.0, 0.0, 0.0];
        let impulse = sfc.solve_velocity(rel_vel);
        assert!(
            impulse[0] < 0.0,
            "impulse should oppose tangential velocity"
        );
    }
    #[test]
    fn static_friction_constraint_clamped_to_cone() {
        let mu = 0.3;
        let normal_impulse = 1.0;
        let mut sfc = StaticFrictionConstraint::new(
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            0.1,
            mu,
            normal_impulse,
        );
        sfc.solve_velocity([100.0, 0.0, 100.0]);
        let total = (sfc.lambda_t1 * sfc.lambda_t1 + sfc.lambda_t2 * sfc.lambda_t2).sqrt();
        assert!(
            total <= mu * normal_impulse + 1e-10,
            "accumulated impulse should stay within cone: |λ|={total}"
        );
    }
    #[test]
    fn static_friction_constraint_at_limit_detection() {
        let mu = 0.3;
        let mut sfc = StaticFrictionConstraint::new(
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            0.01,
            mu,
            1.0,
        );
        sfc.solve_velocity([1000.0, 0.0, 0.0]);
        assert!(
            sfc.is_at_limit(),
            "should be at friction limit after large velocity"
        );
    }
    #[test]
    fn static_friction_constraint_reset() {
        let mut sfc = StaticFrictionConstraint::new(
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            2.0,
            0.5,
            10.0,
        );
        sfc.solve_velocity([5.0, 0.0, 3.0]);
        sfc.reset();
        assert!(sfc.lambda_t1.abs() < 1e-15);
        assert!(sfc.lambda_t2.abs() < 1e-15);
    }
    #[test]
    fn pairwise_friction_solver_single_contact() {
        let mut solver = PairwiseFrictionSolver::new(0.4);
        solver.set_contact_count(1);
        solver.set_total_normal_impulse(10.0);
        let imp = solver.apply_contact(0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0);
        assert!(imp[0] < 0.0, "friction should oppose tangential velocity");
    }
    #[test]
    fn pairwise_friction_solver_budget_per_contact() {
        let mut solver = PairwiseFrictionSolver::new(0.5);
        solver.set_contact_count(4);
        solver.set_total_normal_impulse(8.0);
        let imp = solver.apply_contact(0, [100.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0);
        assert!(
            length(imp) <= 1.0 + 1e-10,
            "impulse should be clamped to per-contact budget"
        );
    }
    #[test]
    fn pairwise_friction_solver_total_impulse() {
        let mut solver = PairwiseFrictionSolver::new(0.3);
        solver.set_contact_count(2);
        solver.set_total_normal_impulse(10.0);
        solver.apply_contact(0, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0);
        solver.apply_contact(1, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0);
        let total = solver.total_tangential_impulse();
        assert!(total[0] < 0.0, "total friction should oppose motion");
    }
    #[test]
    fn friction_dissipation_records_positive_work() {
        let mut diss = FrictionDissipation::new();
        diss.record([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(diss.current_step_dissipation > 0.0);
    }
    #[test]
    fn friction_dissipation_advance_step() {
        let mut diss = FrictionDissipation::new();
        diss.record([-2.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        diss.advance_step();
        assert_eq!(diss.step_count, 1);
        assert!(diss.total_dissipation > 0.0);
        assert!((diss.current_step_dissipation).abs() < 1e-15);
    }
    #[test]
    fn friction_dissipation_peak_tracking() {
        let mut diss = FrictionDissipation::new();
        diss.record([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        diss.advance_step();
        diss.record([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        diss.advance_step();
        diss.record([-0.1, 0.0, 0.0], [1.0, 0.0, 0.0]);
        diss.advance_step();
        assert!(
            (diss.peak_dissipation - 5.0).abs() < 1e-10,
            "peak={}",
            diss.peak_dissipation
        );
    }
    #[test]
    fn friction_dissipation_average() {
        let mut diss = FrictionDissipation::new();
        diss.record([-2.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        diss.advance_step();
        diss.record([-4.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        diss.advance_step();
        let avg = diss.average_dissipation();
        assert!((avg - 3.0).abs() < 1e-10, "avg={avg}");
    }
    #[test]
    fn friction_cache_entry_warm_start() {
        let mut entry = FrictionCacheEntry::new(0, 1, 42);
        entry.accumulated_impulse = [4.0, 0.0, 0.0];
        let ws = entry.warm_start_impulse(0.8);
        assert!((ws[0] - 3.2).abs() < 1e-12, "warm_start scaled={}", ws[0]);
    }
    #[test]
    fn friction_cache_entry_decay() {
        let mut entry = FrictionCacheEntry::new(0, 1, 99);
        entry.accumulated_impulse = [10.0, 0.0, 0.0];
        entry.decay(0.9);
        assert!((entry.accumulated_impulse[0] - 9.0).abs() < 1e-12);
    }
    #[test]
    fn friction_cache_entry_update_age() {
        let mut entry = FrictionCacheEntry::new(0, 1, 7);
        entry.update([1.0, 2.0, 3.0], true);
        assert_eq!(entry.age, 1);
        assert!(entry.was_static);
    }
    #[test]
    fn friction_contact_cache_insert_and_get() {
        let mut cache = FrictionContactCache::new(10);
        let entry = FrictionCacheEntry::new(0, 1, 100);
        cache.insert(entry);
        assert!(cache.get(0, 1, 100).is_some());
        assert!(cache.get(0, 1, 999).is_none());
    }
    #[test]
    fn friction_contact_cache_symmetric_lookup() {
        let mut cache = FrictionContactCache::new(10);
        let entry = FrictionCacheEntry::new(5, 3, 7);
        cache.insert(entry);
        assert!(
            cache.get(3, 5, 7).is_some(),
            "reversed body order should also find entry"
        );
    }
    #[test]
    fn friction_contact_cache_update_existing() {
        let mut cache = FrictionContactCache::new(10);
        let mut entry = FrictionCacheEntry::new(0, 1, 1);
        entry.accumulated_impulse = [1.0, 0.0, 0.0];
        cache.insert(entry);
        let mut updated = FrictionCacheEntry::new(0, 1, 1);
        updated.accumulated_impulse = [5.0, 0.0, 0.0];
        cache.insert(updated);
        assert_eq!(cache.len(), 1, "should not duplicate");
        assert!((cache.get(0, 1, 1).unwrap().accumulated_impulse[0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn friction_contact_cache_age_and_evict() {
        let mut cache = FrictionContactCache::new(100);
        let mut e = FrictionCacheEntry::new(0, 1, 1);
        e.age = 10;
        cache.insert(e);
        cache.age_and_evict(5);
        assert!(cache.is_empty(), "old entries should be evicted");
    }
    #[test]
    fn friction_contact_cache_clear() {
        let mut cache = FrictionContactCache::new(10);
        cache.insert(FrictionCacheEntry::new(0, 1, 1));
        cache.insert(FrictionCacheEntry::new(2, 3, 2));
        cache.clear();
        assert!(cache.is_empty());
    }
}
#[cfg(test)]
mod tests_friction_new {

    use crate::AnisotropicFrictionCoeffs;
    use crate::AnisotropicFrictionSolver;
    use crate::ConeFrictionProjector;

    use crate::PatchFriction;
    use crate::SpinFriction;
    use crate::StribeckParams;
    use crate::ThermalFriction;
    #[test]
    fn test_anisotropic_isotropic_constructor() {
        let c = AnisotropicFrictionCoeffs::isotropic(0.8, 0.6);
        assert!((c.mu_s_primary - 0.8).abs() < 1e-12);
        assert!((c.mu_k_primary - 0.6).abs() < 1e-12);
        assert!((c.mu_s_secondary - 0.8).abs() < 1e-12);
        assert!((c.mu_k_secondary - 0.6).abs() < 1e-12);
    }
    #[test]
    fn test_anisotropic_solver_inside_ellipse_no_slip() {
        let coeffs = AnisotropicFrictionCoeffs::new(0.8, 0.6, 0.4, 0.3);
        let solver = AnisotropicFrictionSolver::new(coeffs, 10.0);
        let (clamped, slipping) = solver.clamp_impulse(1.0, 0.5);
        assert!(!slipping, "should be static");
        assert!((clamped[0] - 1.0).abs() < 1e-12);
        assert!((clamped[1] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_anisotropic_solver_outside_ellipse_slips() {
        let coeffs = AnisotropicFrictionCoeffs::new(0.8, 0.6, 0.4, 0.3);
        let solver = AnisotropicFrictionSolver::new(coeffs, 10.0);
        let (_clamped, slipping) = solver.clamp_impulse(100.0, 100.0);
        assert!(slipping, "should slip");
    }
    #[test]
    fn test_anisotropic_solver_kinetic_clamped_values() {
        let coeffs = AnisotropicFrictionCoeffs::new(0.8, 0.6, 0.4, 0.3);
        let solver = AnisotropicFrictionSolver::new(coeffs, 10.0);
        let (clamped, slipping) = solver.clamp_impulse(100.0, 100.0);
        assert!(slipping);
        assert!(clamped[0].abs() <= 6.0 + 1e-12);
        assert!(clamped[1].abs() <= 3.0 + 1e-12);
    }
    #[test]
    fn test_anisotropic_friction_force_3d_direction() {
        let coeffs = AnisotropicFrictionCoeffs::isotropic(0.8, 0.6);
        let solver = AnisotropicFrictionSolver::new(coeffs, 5.0);
        let f = solver.friction_force_3d([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.01, 0.0, 0.016);
        assert!(f[0] < 0.0 || f[0].abs() < 1e-12);
    }
    #[test]
    fn test_stribeck_at_zero_speed_equals_static() {
        let p = StribeckParams::new(0.4, 0.8, 0.1, 0.0);
        assert!((p.mu_at(0.0) - 0.8).abs() < 1e-10, "mu={}", p.mu_at(0.0));
    }
    #[test]
    fn test_stribeck_high_speed_approaches_kinetic() {
        let p = StribeckParams::new(0.4, 0.8, 0.1, 0.0);
        let mu_high = p.mu_at(100.0);
        assert!((mu_high - 0.4).abs() < 1e-4, "mu={mu_high}");
    }
    #[test]
    fn test_stribeck_friction_force_scales_with_fn() {
        let p = StribeckParams::new(0.4, 0.6, 0.5, 0.0);
        let f1 = p.friction_force(10.0, 0.0);
        let f2 = p.friction_force(20.0, 0.0);
        assert!((f2 - 2.0 * f1).abs() < 1e-10, "f1={f1} f2={f2}");
    }
    #[test]
    fn test_stribeck_signed_force_opposes_motion() {
        let p = StribeckParams::new(0.4, 0.8, 0.1, 0.0);
        let f_pos = p.friction_force_signed(10.0, 2.0);
        let f_neg = p.friction_force_signed(10.0, -2.0);
        assert!(f_pos < 0.0, "friction should oppose positive velocity");
        assert!(f_neg > 0.0, "friction should oppose negative velocity");
    }
    #[test]
    fn test_stribeck_zero_velocity_returns_zero() {
        let p = StribeckParams::new(0.4, 0.8, 0.1, 0.0);
        assert!((p.friction_force_signed(10.0, 0.0)).abs() < 1e-14);
    }
    #[test]
    fn test_cone_projector_inside_no_clamp() {
        let proj = ConeFrictionProjector::new(0.5);
        let (imp, slipping) = proj.project([0.5, 0.5], 10.0);
        assert!(!slipping);
        assert!((imp[0] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_cone_projector_outside_clamped() {
        let proj = ConeFrictionProjector::new(0.5);
        let (imp, slipping) = proj.project([10.0, 0.0], 5.0);
        assert!(slipping);
        assert!((imp[0] - 2.5).abs() < 1e-10, "clamped={}", imp[0]);
        assert!(imp[1].abs() < 1e-12);
    }
    #[test]
    fn test_cone_projector_3d_inside() {
        let proj = ConeFrictionProjector::new(0.8);
        let t = [0.1, 0.1, 0.1];
        let (_out, slipping) = proj.project_3d(t, 10.0);
        assert!(!slipping);
    }
    #[test]
    fn test_cone_projector_3d_clamped_magnitude() {
        let proj = ConeFrictionProjector::new(0.5);
        let t = [0.0, 0.0, 20.0];
        let (out, slipping) = proj.project_3d(t, 5.0);
        assert!(slipping);
        let mag = (out[0] * out[0] + out[1] * out[1] + out[2] * out[2]).sqrt();
        assert!((mag - 2.5).abs() < 1e-10, "mag={mag}");
    }
    #[test]
    fn test_patch_friction_translational_two_thirds() {
        let pf = PatchFriction::new(0.6, 0.05);
        let f = pf.translational_force(10.0);
        assert!((f - 4.0).abs() < 1e-10, "f={f}");
    }
    #[test]
    fn test_patch_friction_spin_torque() {
        let pf = PatchFriction::new(0.6, 0.05);
        let t = pf.spin_torque(10.0);
        assert!((t - 0.2).abs() < 1e-10, "t={t}");
    }
    #[test]
    fn test_patch_friction_force_3d_opposes_velocity() {
        let pf = PatchFriction::new(0.8, 0.1);
        let v = [1.0, 0.0, 0.0];
        let f = pf.friction_force_3d(v, 5.0);
        assert!(f[0] < 0.0, "force should oppose velocity in x");
        assert!(f[1].abs() < 1e-12);
    }
    #[test]
    fn test_patch_friction_force_zero_velocity() {
        let pf = PatchFriction::new(0.8, 0.1);
        let f = pf.friction_force_3d([0.0, 0.0, 0.0], 10.0);
        assert!(f[0].abs() < 1e-14);
        assert!(f[1].abs() < 1e-14);
        assert!(f[2].abs() < 1e-14);
    }
    #[test]
    fn test_spin_friction_zero_omega_returns_zero() {
        let sf = SpinFriction::new(1.0, 0.01);
        assert!(sf.torque_signed(0.0).abs() < 1e-14);
    }
    #[test]
    fn test_spin_friction_high_omega_approaches_max() {
        let sf = SpinFriction::new(5.0, 0.001);
        let t = sf.torque_magnitude(100.0);
        assert!((t - 5.0).abs() < 1e-6, "t={t}");
    }
    #[test]
    fn test_spin_friction_signed_opposes_omega() {
        let sf = SpinFriction::new(2.0, 0.1);
        assert!(sf.torque_signed(1.0) < 0.0);
        assert!(sf.torque_signed(-1.0) > 0.0);
    }
    #[test]
    fn test_thermal_friction_at_reference_temp() {
        let tf = ThermalFriction::new(0.7, 300.0, -0.001, 0.1, 1.0);
        assert!((tf.mu_at(300.0) - 0.7).abs() < 1e-12);
    }
    #[test]
    fn test_thermal_friction_decreases_with_temperature() {
        let tf = ThermalFriction::new(0.7, 300.0, -0.001, 0.1, 1.0);
        assert!(tf.mu_at(400.0) < tf.mu_at(300.0));
    }
    #[test]
    fn test_thermal_friction_clamped_to_min() {
        let tf = ThermalFriction::new(0.7, 300.0, -1.0, 0.1, 1.0);
        assert!((tf.mu_at(1000.0) - 0.1).abs() < 1e-12);
    }
    #[test]
    fn test_thermal_friction_clamped_to_max() {
        let tf = ThermalFriction::new(0.7, 300.0, 1.0, 0.1, 1.0);
        assert!((tf.mu_at(1000.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_thermal_friction_force_scales_with_normal_force() {
        let tf = ThermalFriction::new(0.5, 20.0, 0.0, 0.0, 2.0);
        let f1 = tf.friction_force(10.0, 20.0);
        let f2 = tf.friction_force(20.0, 20.0);
        assert!((f2 - 2.0 * f1).abs() < 1e-10, "f1={f1} f2={f2}");
    }
}
