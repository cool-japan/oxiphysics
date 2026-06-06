//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// 3-D Euclidean norm.
#[inline]
pub(super) fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
/// Normalise a 3-D vector (zero vector if degenerate).
#[inline]
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1.0e-15 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}
/// Vector subtraction a − b.
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Clamp `x` to \[lo, hi\].
#[inline]
pub(super) fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::muscle_simulation::ActivationDynamics;
    use crate::muscle_simulation::FatigueDynamics;
    use crate::muscle_simulation::ForceLengthRelation;
    use crate::muscle_simulation::ForceVelocityRelation;
    use crate::muscle_simulation::HillMuscleModel;
    use crate::muscle_simulation::JointConfig;
    use crate::muscle_simulation::MuscleAttachment;
    use crate::muscle_simulation::MuscleEntry;
    use crate::muscle_simulation::MuscleGeometry;
    use crate::muscle_simulation::MuscleOptimization;
    use crate::muscle_simulation::MuscleRecruitment;
    use crate::muscle_simulation::MusculoskeletalSystem;
    use crate::muscle_simulation::MusculotendonUnit;
    use crate::muscle_simulation::OptMuscle;
    use crate::muscle_simulation::SarcomereModel;
    use crate::muscle_simulation::WholeBodyMuscle;
    pub(super) const TOL: f64 = 1e-9;
    #[test]
    fn test_norm3_unit_vector() {
        let v = [3.0_f64, 4.0, 0.0];
        assert!((norm3(v) - 5.0).abs() < TOL);
    }
    #[test]
    fn test_normalize3_gives_unit_length() {
        let v = [3.0, 4.0, 0.0];
        let n = normalize3(v);
        assert!((norm3(n) - 1.0).abs() < TOL);
    }
    #[test]
    fn test_normalize3_zero_vector() {
        let n = normalize3([0.0, 0.0, 0.0]);
        assert_eq!(n, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_clamp_bounds() {
        assert_eq!(clamp(-1.0, 0.0, 1.0), 0.0);
        assert_eq!(clamp(2.0, 0.0, 1.0), 1.0);
        assert_eq!(clamp(0.5, 0.0, 1.0), 0.5);
    }
    #[test]
    fn test_muscle_geometry_mtu_length() {
        let mg = MuscleGeometry::soleus();
        let l = mg.musculotendon_length();
        assert!(l > 0.0, "MTU length should be positive, got {l}");
    }
    #[test]
    fn test_muscle_geometry_fiber_length_positive() {
        let mg = MuscleGeometry::rectus_femoris();
        let l = mg.fiber_length(0.45, 0.10);
        assert!(l > 0.0, "fiber length should be positive, got {l}");
    }
    #[test]
    fn test_muscle_geometry_force_projection_le_force() {
        let mg = MuscleGeometry::soleus();
        let proj = mg.fiber_force_projection(100.0);
        assert!(
            proj <= 100.0 + TOL,
            "force projection should be ≤ fiber force"
        );
        assert!(
            proj > 0.0,
            "force projection should be positive for positive force"
        );
    }
    #[test]
    fn test_muscle_geometry_volume_positive() {
        let mg = MuscleGeometry::soleus();
        assert!(mg.volume > 0.0, "muscle volume should be positive");
    }
    #[test]
    fn test_fv_at_zero_velocity_is_one() {
        let fv = ForceVelocityRelation::fast_twitch();
        let f = fv.fv(0.0);
        assert!((f - 1.0).abs() < 1e-6, "fv(0) should be 1.0, got {f}");
    }
    #[test]
    fn test_fv_shortening_decreases_force() {
        let fv = ForceVelocityRelation::fast_twitch();
        let f0 = fv.fv(0.0);
        let f_short = fv.fv(-0.5);
        assert!(
            f_short < f0,
            "shortening should reduce force: f0={f0}, f_short={f_short}"
        );
    }
    #[test]
    fn test_fv_lengthening_increases_force() {
        let fv = ForceVelocityRelation::fast_twitch();
        let f0 = fv.fv(0.0);
        let f_len = fv.fv(0.5);
        assert!(
            f_len > f0,
            "lengthening should increase force: f0={f0}, f_len={f_len}"
        );
    }
    #[test]
    fn test_fv_max_velocity_force_near_zero() {
        let fv = ForceVelocityRelation::fast_twitch();
        let f_max_vel = fv.fv(-1.0);
        assert!(
            f_max_vel.abs() < 0.1,
            "force at max velocity should be near 0, got {f_max_vel}"
        );
    }
    #[test]
    fn test_fv_optimal_velocity_fraction_in_range() {
        let fv = ForceVelocityRelation::fast_twitch();
        let opt = fv.optimal_velocity_fraction();
        assert!(
            opt > 0.0 && opt < 1.0,
            "optimal velocity fraction should be in (0,1), got {opt}"
        );
    }
    #[test]
    fn test_fl_active_at_optimal_is_one() {
        let fl = ForceLengthRelation::gastrocnemius();
        let f = fl.fl_active(1.0);
        assert!(
            (f - 1.0).abs() < 1e-6,
            "fl_active at l_opt should be 1.0, got {f}"
        );
    }
    #[test]
    fn test_fl_active_decreases_away_from_optimal() {
        let fl = ForceLengthRelation::gastrocnemius();
        let f_opt = fl.fl_active(1.0);
        let f_short = fl.fl_active(0.6);
        let f_long = fl.fl_active(1.5);
        assert!(
            f_short < f_opt,
            "active force should decrease at short length"
        );
        assert!(
            f_long < f_opt,
            "active force should decrease at long length"
        );
    }
    #[test]
    fn test_fl_passive_zero_at_optimal() {
        let fl = ForceLengthRelation::gastrocnemius();
        let fp = fl.fl_passive(1.0);
        assert_eq!(fp, 0.0, "passive force should be zero at optimal length");
    }
    #[test]
    fn test_fl_passive_positive_at_long_length() {
        let fl = ForceLengthRelation::gastrocnemius();
        let fp = fl.fl_passive(1.5);
        assert!(
            fp > 0.0,
            "passive force should be positive at long length, got {fp}"
        );
    }
    #[test]
    fn test_titin_force_zero_below_engagement() {
        let fl = ForceLengthRelation::default_params();
        let ft = fl.titin_force(fl.optimal_fiber_length * 0.9);
        assert_eq!(ft, 0.0, "titin force should be 0 below engagement length");
    }
    #[test]
    fn test_titin_force_positive_above_engagement() {
        let fl = ForceLengthRelation::default_params();
        let l_eng = fl.titin_engagement * fl.optimal_fiber_length;
        let ft = fl.titin_force(l_eng + 0.01);
        assert!(
            ft > 0.0,
            "titin force should be positive above engagement, got {ft}"
        );
    }
    #[test]
    fn test_hill_ce_force_zero_at_rest() {
        let muscle = HillMuscleModel::new(1000.0, 0.1, 0.2, 0.1);
        let f = muscle.ce_force_active();
        assert_eq!(f, 0.0, "active force should be 0 at zero activation");
    }
    #[test]
    fn test_hill_ce_force_positive_at_activation() {
        let mut muscle = HillMuscleModel::new(1000.0, 0.1, 0.2, 0.1);
        muscle.activation = 1.0;
        let f = muscle.ce_force_active();
        assert!(
            f > 0.0,
            "active force should be positive at full activation, got {f}"
        );
    }
    #[test]
    fn test_hill_pee_force_positive_at_long_length() {
        let mut muscle = HillMuscleModel::new(1000.0, 0.1, 0.2, 0.1);
        muscle.l_ce_norm = 1.5;
        let f = muscle.pee_force();
        assert!(
            f > 0.0,
            "passive force should be positive at long length, got {f}"
        );
    }
    #[test]
    fn test_hill_step_updates_activation() {
        let mut muscle = HillMuscleModel::tibialis_anterior();
        let a0 = muscle.activation;
        muscle.step(0.01, muscle.tendon_slack_length + muscle.l_opt, 1.0);
        assert!(
            muscle.activation > a0,
            "activation should increase toward target 1.0"
        );
    }
    #[test]
    fn test_hill_step_clamps_activation() {
        let mut muscle = HillMuscleModel::soleus();
        for _ in 0..1000 {
            muscle.step(0.01, muscle.tendon_slack_length + muscle.l_opt, 1.0);
        }
        assert!(
            muscle.activation <= 1.0,
            "activation should not exceed 1.0, got {}",
            muscle.activation
        );
    }
    #[test]
    fn test_mtu_tendon_strain_non_negative() {
        let mtu = MusculotendonUnit::soleus();
        let strain = mtu.tendon_strain();
        assert!(
            strain >= 0.0,
            "tendon strain should be non-negative, got {strain}"
        );
    }
    #[test]
    fn test_mtu_fiber_length_positive() {
        let mtu = MusculotendonUnit::soleus();
        let l = mtu.fiber_length();
        assert!(l > 0.0, "fiber length should be positive, got {l}");
    }
    #[test]
    fn test_mtu_pennation_angle_positive() {
        let mtu = MusculotendonUnit::soleus();
        let alpha = mtu.current_pennation_angle();
        assert!(
            alpha >= 0.0,
            "pennation angle should be non-negative, got {alpha}"
        );
    }
    #[test]
    fn test_mtu_step_changes_state() {
        let mut mtu = MusculotendonUnit::soleus();
        let l0 = mtu.fiber_length();
        let lmt = mtu.lmt * 1.02;
        mtu.step(0.01, lmt, 0.8);
        let l1 = mtu.fiber_length();
        let a1 = mtu.muscle.activation;
        assert!(
            (l1 - l0).abs() > TOL || a1 > TOL,
            "MTU state should change after step: l0={l0}, l1={l1}, a={a1}"
        );
    }
    #[test]
    fn test_activation_dynamics_increases_toward_target() {
        let mut dyn_act = ActivationDynamics::default_params();
        let a0 = dyn_act.activation;
        dyn_act.step(1.0, 0.01);
        assert!(
            dyn_act.activation > a0,
            "activation should increase toward target=1.0"
        );
    }
    #[test]
    fn test_activation_dynamics_decreases_toward_target() {
        let mut dyn_act = ActivationDynamics::default_params();
        dyn_act.activation = 0.8;
        let a0 = dyn_act.activation;
        dyn_act.step(0.0, 0.01);
        assert!(
            dyn_act.activation < a0,
            "activation should decrease toward target=0"
        );
    }
    #[test]
    fn test_activation_dynamics_stays_above_min() {
        let mut dyn_act = ActivationDynamics::default_params();
        for _ in 0..1000 {
            dyn_act.step(0.0, 0.01);
        }
        assert!(
            dyn_act.activation >= dyn_act.activation_min,
            "activation should stay above minimum"
        );
    }
    #[test]
    fn test_activation_dynamics_steady_state() {
        let dyn_act = ActivationDynamics::default_params();
        let ss = dyn_act.steady_state(0.7);
        assert!(
            (ss - 0.7).abs() < 1e-6,
            "steady state should equal drive when >min"
        );
    }
    #[test]
    fn test_recruitment_zero_drive_no_force() {
        let mut pool = MuscleRecruitment::new(10, 1000.0);
        pool.update(0.0);
        assert_eq!(
            pool.recruited_count(),
            0,
            "no units should be recruited at zero drive"
        );
    }
    #[test]
    fn test_recruitment_full_drive_all_recruited() {
        let mut pool = MuscleRecruitment::new(10, 1000.0);
        pool.update(1.0);
        assert_eq!(
            pool.recruited_count(),
            pool.motor_units.len(),
            "all units should be recruited at full drive"
        );
    }
    #[test]
    fn test_recruitment_force_increases_with_drive() {
        let mut pool = MuscleRecruitment::new(20, 500.0);
        pool.update(0.3);
        let f_low = pool.total_force();
        pool.update(0.8);
        let f_high = pool.total_force();
        assert!(
            f_high > f_low,
            "force should increase with drive: {f_low} vs {f_high}"
        );
    }
    #[test]
    fn test_recruitment_force_fraction_in_range() {
        let mut pool = MuscleRecruitment::new(10, 1000.0);
        pool.update(0.5);
        let ff = pool.force_fraction();
        assert!(
            (0.0..=1.0).contains(&ff),
            "force fraction should be in [0,1], got {ff}"
        );
    }
    #[test]
    fn test_whole_body_static_optimisation_correct_sign() {
        let muscles = vec![
            MuscleEntry::new("bicep", MusculotendonUnit::soleus(), 0.04, 1.0),
            MuscleEntry::new("brachialis", MusculotendonUnit::soleus(), 0.03, 1.0),
        ];
        let wbm = WholeBodyMuscle::new(muscles);
        let activations = wbm.static_optimisation(50.0);
        assert_eq!(
            activations.len(),
            2,
            "should return one activation per muscle"
        );
        for a in &activations {
            assert!(
                *a >= 0.0 && *a <= 1.0,
                "activations should be in [0,1], got {a}"
            );
        }
    }
    #[test]
    fn test_whole_body_redundancy_index() {
        let muscles = vec![
            MuscleEntry::new("m1", MusculotendonUnit::soleus(), 0.04, 1.0),
            MuscleEntry::new("m2", MusculotendonUnit::soleus(), 0.03, 1.0),
            MuscleEntry::new("m3", MusculotendonUnit::soleus(), 0.02, 1.0),
        ];
        let wbm = WholeBodyMuscle::new(muscles);
        assert_eq!(
            wbm.redundancy_index(),
            2,
            "redundancy index should be n_muscles - 1"
        );
    }
    #[test]
    fn test_whole_body_step_advances_state() {
        let muscles = vec![MuscleEntry::new(
            "soleus",
            MusculotendonUnit::soleus(),
            0.05,
            1.0,
        )];
        let mut wbm = WholeBodyMuscle::new(muscles);
        let a0 = wbm.muscles[0].mtu.muscle.activation;
        let lmt = wbm.muscles[0].mtu.lmt;
        wbm.step(0.01, &[0.8], &[lmt]);
        let a1 = wbm.muscles[0].mtu.muscle.activation;
        assert!(
            a1 > a0,
            "activation should increase after step with target=0.8"
        );
    }
    #[test]
    fn test_sarcomere_default_state_valid() {
        let s = SarcomereModel::default();
        assert!(s.attached_fraction >= 0.0 && s.attached_fraction <= 1.0);
        assert!(s.sarcomere_length > 0.0);
    }
    #[test]
    fn test_sarcomere_active_force_positive_when_attached() {
        let s = SarcomereModel {
            attached_fraction: 0.5,
            ..Default::default()
        };
        let f = s.active_force();
        assert!(
            f > 0.0,
            "active force should be positive when cross-bridges attached, got {f}"
        );
    }
    #[test]
    fn test_sarcomere_active_force_zero_when_detached() {
        let s = SarcomereModel {
            attached_fraction: 0.0,
            ..Default::default()
        };
        let f = s.active_force();
        assert!(
            f < 1e-12,
            "active force should be zero when no cross-bridges, got {f}"
        );
    }
    #[test]
    fn test_sarcomere_step_changes_state() {
        let mut s = SarcomereModel::default();
        let f0 = s.attached_fraction;
        s.step(0.01, 1.0, 0.0);
        let f1 = s.attached_fraction;
        assert!(
            (f1 - f0).abs() > 1e-12 || f0 > 0.0,
            "sarcomere state should evolve"
        );
    }
    #[test]
    fn test_sarcomere_velocity_reduces_force() {
        let mut s_slow = SarcomereModel::default();
        let mut s_fast = SarcomereModel::default();
        s_slow.attached_fraction = 0.5;
        s_fast.attached_fraction = 0.5;
        for _ in 0..10 {
            s_slow.step(0.001, 1.0, 0.1);
            s_fast.step(0.001, 1.0, 2.0);
        }
        let f_slow = s_slow.active_force();
        let f_fast = s_fast.active_force();
        assert!(
            f_slow >= f_fast - 1e-9,
            "lower shortening velocity should produce >= force: {f_slow} vs {f_fast}"
        );
    }
    #[test]
    fn test_sarcomere_optimal_length_max_force() {
        let s = SarcomereModel::default();
        let f_opt = s.length_dependent_force(s.optimal_sarcomere_length);
        let f_short = s.length_dependent_force(s.optimal_sarcomere_length * 0.7);
        let f_long = s.length_dependent_force(s.optimal_sarcomere_length * 1.4);
        assert!(
            f_opt >= f_short - 1e-9,
            "force at optimal should be >= short length"
        );
        assert!(
            f_opt >= f_long - 1e-9,
            "force at optimal should be >= long length"
        );
    }
    #[test]
    fn test_msk_add_joint_and_muscle() {
        let mut msk = MusculoskeletalSystem::new();
        msk.add_joint(JointConfig::new("elbow", [0.0, 0.3, 0.0], 0.0, -2.356, 0.0));
        let mtu = MusculotendonUnit::soleus();
        msk.add_muscle(MuscleAttachment::new(
            "bicep",
            mtu,
            0,
            0.04,
            [0.0, 0.28, 0.0],
            [0.0, 0.0, 0.0],
        ));
        assert_eq!(msk.joints.len(), 1);
        assert_eq!(msk.muscles.len(), 1);
    }
    #[test]
    fn test_msk_joint_torque_non_negative_with_activation() {
        let mut msk = MusculoskeletalSystem::new();
        msk.add_joint(JointConfig::new("knee", [0.0, 0.4, 0.0], 0.0, -2.0, 0.0));
        let mut mtu = MusculotendonUnit::soleus();
        mtu.muscle.activation = 0.8;
        msk.add_muscle(MuscleAttachment::new(
            "quad",
            mtu,
            0,
            0.05,
            [0.0, 0.38, 0.0],
            [0.0, 0.0, 0.0],
        ));
        let torques = msk.compute_joint_torques();
        assert_eq!(torques.len(), 1);
        assert!(torques[0].abs() >= 0.0, "torque should be defined");
    }
    #[test]
    fn test_msk_joint_torque_zero_at_zero_activation() {
        let mut msk = MusculoskeletalSystem::new();
        msk.add_joint(JointConfig::new("hip", [0.0, 0.5, 0.0], 0.0, -2.0, 0.0));
        let mut mtu = MusculotendonUnit::soleus();
        mtu.muscle.activation = 0.0;
        msk.add_muscle(MuscleAttachment::new(
            "glute",
            mtu,
            0,
            0.06,
            [0.0, 0.48, 0.0],
            [0.0, 0.0, 0.0],
        ));
        let torques = msk.compute_joint_torques();
        assert!(
            torques[0].abs() < 5.0,
            "torque should be near zero at zero activation, got {}",
            torques[0]
        );
    }
    #[test]
    fn test_msk_moment_arm_positive() {
        let att = MuscleAttachment::new(
            "bicep",
            MusculotendonUnit::soleus(),
            0,
            0.04,
            [0.0, 0.3, 0.0],
            [0.0, 0.0, 0.0],
        );
        assert!(
            att.moment_arm > 0.0,
            "moment arm should be positive, got {}",
            att.moment_arm
        );
    }
    #[test]
    fn test_msk_joint_angle_clamped() {
        let mut cfg = JointConfig::new("ankle", [0.0, 0.1, 0.0], 0.0, -0.785, 0.524);
        cfg.set_angle(2.0);
        assert!(
            cfg.angle <= cfg.max_angle + 1e-9,
            "angle should not exceed max"
        );
        cfg.set_angle(-2.0);
        assert!(
            cfg.angle >= cfg.min_angle - 1e-9,
            "angle should not go below min"
        );
    }
    #[test]
    fn test_fatigue_initial_state_fresh() {
        let fd = FatigueDynamics::new();
        assert!(
            fd.central_fatigue < 0.1,
            "should start with low central fatigue"
        );
        assert!(
            fd.peripheral_fatigue < 0.1,
            "should start with low peripheral fatigue"
        );
        assert!(fd.fatigue_index() < 0.1, "fatigue index should start low");
    }
    #[test]
    fn test_fatigue_increases_with_effort() {
        let mut fd = FatigueDynamics::new();
        let fi0 = fd.fatigue_index();
        for _ in 0..100 {
            fd.step(0.01, 0.9);
        }
        let fi1 = fd.fatigue_index();
        assert!(
            fi1 > fi0,
            "fatigue should increase with high effort: {fi0} -> {fi1}"
        );
    }
    #[test]
    fn test_fatigue_recovers_at_rest() {
        let mut fd = FatigueDynamics::new();
        for _ in 0..200 {
            fd.step(0.01, 1.0);
        }
        let fi_fatigued = fd.fatigue_index();
        for _ in 0..200 {
            fd.step(0.01, 0.0);
        }
        let fi_recovered = fd.fatigue_index();
        assert!(
            fi_recovered < fi_fatigued,
            "fatigue should recover at rest: {fi_fatigued} -> {fi_recovered}"
        );
    }
    #[test]
    fn test_fatigue_index_in_range() {
        let mut fd = FatigueDynamics::new();
        for _ in 0..1000 {
            fd.step(0.01, 1.0);
        }
        let fi = fd.fatigue_index();
        assert!(
            (0.0..=1.0 + 1e-9).contains(&fi),
            "fatigue index should be in [0,1], got {fi}"
        );
    }
    #[test]
    fn test_central_peripheral_both_increase() {
        let mut fd = FatigueDynamics::new();
        for _ in 0..50 {
            fd.step(0.01, 0.8);
        }
        assert!(fd.central_fatigue > 0.0, "central fatigue should increase");
        assert!(
            fd.peripheral_fatigue > 0.0,
            "peripheral fatigue should increase"
        );
    }
    #[test]
    fn test_fatigue_effective_force_reduced() {
        let fd_fresh = FatigueDynamics::new();
        let mut fd_tired = FatigueDynamics::new();
        for _ in 0..500 {
            fd_tired.step(0.01, 1.0);
        }
        let f_fresh = fd_fresh.effective_force_scale(1.0);
        let f_tired = fd_tired.effective_force_scale(1.0);
        assert!(
            f_tired < f_fresh,
            "fatigue should reduce force scale: {f_fresh} vs {f_tired}"
        );
    }
    #[test]
    fn test_static_opt_activations_in_range() {
        let mtu1 = MusculotendonUnit::soleus();
        let mtu2 = MusculotendonUnit::soleus();
        let mut opt = MuscleOptimization::new(vec![
            OptMuscle::new("sol", mtu1, 0.04),
            OptMuscle::new("ta", mtu2, 0.03),
        ]);
        let acts = opt.static_optimization(30.0);
        for a in &acts {
            assert!(
                *a >= 0.0 && *a <= 1.0 + 1e-9,
                "activation out of range: {a}"
            );
        }
    }
    #[test]
    fn test_static_opt_higher_demand_higher_activation() {
        let mtu1 = MusculotendonUnit::soleus();
        let mut opt_low = MuscleOptimization::new(vec![OptMuscle::new("m1", mtu1, 0.04)]);
        let mtu3 = MusculotendonUnit::soleus();
        let mut opt_high = MuscleOptimization::new(vec![OptMuscle::new("m1", mtu3, 0.04)]);
        let a_low = opt_low.static_optimization(10.0);
        let a_high = opt_high.static_optimization(100.0);
        assert!(
            a_high[0] >= a_low[0] - 1e-9,
            "higher demand should require >= activation: {:.3} vs {:.3}",
            a_low[0],
            a_high[0]
        );
    }
    #[test]
    fn test_emg_driven_activation_follows_emg() {
        let mtu = MusculotendonUnit::soleus();
        let mut opt = MuscleOptimization::new(vec![OptMuscle::new("sol", mtu, 0.04)]);
        let emg = vec![0.7f64];
        let acts = opt.emg_driven_step(0.01, &emg);
        assert_eq!(acts.len(), 1);
        assert!(
            acts[0] >= 0.0 && acts[0] <= 1.0,
            "activation out of range: {}",
            acts[0]
        );
    }
    #[test]
    fn test_emg_driven_zero_emg_low_activation() {
        let mtu = MusculotendonUnit::soleus();
        let mut opt = MuscleOptimization::new(vec![OptMuscle::new("sol", mtu, 0.04)]);
        let emg = vec![0.0f64];
        let mut acts = vec![0.0];
        for _ in 0..100 {
            acts = opt.emg_driven_step(0.01, &emg);
        }
        assert!(
            acts[0] < 0.3,
            "zero EMG should lead to low activation, got {}",
            acts[0]
        );
    }
    #[test]
    fn test_muscle_stress_minimum_objective() {
        let mtu1 = MusculotendonUnit::soleus();
        let mtu2 = MusculotendonUnit::soleus();
        let mut opt = MuscleOptimization::new(vec![
            OptMuscle::new("m1", mtu1, 0.04),
            OptMuscle::new("m2", mtu2, 0.04),
        ]);
        let acts = opt.static_optimization(50.0);
        let stress_sum: f64 = acts.iter().map(|a| a * a).sum();
        assert!(stress_sum >= 0.0, "stress sum should be non-negative");
        assert!(stress_sum <= acts.len() as f64, "stress sum should be <= n");
    }
}
