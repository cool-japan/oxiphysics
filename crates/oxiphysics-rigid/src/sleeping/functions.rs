//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Kinetic energy of a rigid body.
///
/// KE = 0.5 * mass * |v|² + 0.5 * dot(inertia, omega²)
pub fn kinetic_energy_rigid(mass: f64, vel: [f64; 3], inertia: [f64; 3], omega: [f64; 3]) -> f64 {
    let v2 = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
    let rot = inertia[0] * omega[0] * omega[0]
        + inertia[1] * omega[1] * omega[1]
        + inertia[2] * omega[2] * omega[2];
    0.5 * mass * v2 + 0.5 * rot
}
/// Gravitational potential energy of a rigid body.
///
/// PE = mass * g * height (pos\[1\] is taken as height, Y-up convention).
pub fn potential_energy_rigid(mass: f64, pos: [f64; 3], g: f64) -> f64 {
    mass * g * pos[1]
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::AdaptiveSleepThreshold;

    use crate::EnergyTracker;

    use crate::IslandSleepManager;

    use crate::SleepConfig;

    use crate::SleepHistory;

    use crate::SleepManager;
    use crate::SleepPredictor;
    use crate::SleepState;
    use crate::SleepStats;

    use crate::SleepWakePattern;
    use crate::SleepingSystem;

    fn fast_config() -> SleepConfig {
        SleepConfig {
            linear_threshold: 0.1,
            angular_threshold: 0.1,
            drowsy_frames: 3,
            energy_threshold: 0.0,
            wake_damping: 1.0,
            sleep_time_seconds: 0.0,
        }
    }
    #[test]
    fn test_body_goes_to_sleep() {
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body(1, [0.0; 3]);
        mgr.tick();
        assert!(matches!(mgr.bodies[idx].state, SleepState::Drowsy { .. }));
        mgr.tick();
        assert!(matches!(mgr.bodies[idx].state, SleepState::Drowsy { .. }));
        mgr.tick();
        assert_eq!(mgr.bodies[idx].state, SleepState::Sleeping);
    }
    #[test]
    fn test_body_above_threshold_stays_active() {
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body(2, [0.0; 3]);
        mgr.update_body_velocity(idx, [5.0, 0.0, 0.0], [0.0; 3]);
        for _ in 0..10 {
            mgr.tick();
        }
        assert_eq!(mgr.bodies[idx].state, SleepState::Active);
    }
    #[test]
    fn test_wake_sleeping_body() {
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body(3, [0.0; 3]);
        mgr.bodies[idx].state = SleepState::Sleeping;
        mgr.wake(idx);
        assert_eq!(mgr.bodies[idx].state, SleepState::Active);
    }
    #[test]
    fn test_wake_nearby() {
        let mut mgr = SleepManager::new(fast_config());
        let near = mgr.add_body(4, [1.0, 0.0, 0.0]);
        let far = mgr.add_body(5, [100.0, 0.0, 0.0]);
        mgr.bodies[near].state = SleepState::Sleeping;
        mgr.bodies[far].state = SleepState::Sleeping;
        mgr.wake_nearby([0.0; 3], 5.0);
        assert_eq!(mgr.bodies[near].state, SleepState::Active);
        assert_eq!(mgr.bodies[far].state, SleepState::Sleeping);
    }
    #[test]
    fn test_group_into_islands() {
        let mut mgr = SleepManager::new(fast_config());
        mgr.add_body(10, [0.0; 3]);
        mgr.add_body(11, [1.0, 0.0, 0.0]);
        mgr.add_body(12, [50.0, 0.0, 0.0]);
        let pairs = vec![(0_usize, 1_usize)];
        let mut islands = IslandSleepManager::group_into_islands(&mgr.bodies, &pairs);
        for island in &mut islands {
            island.sort_unstable();
        }
        islands.sort_by_key(|v| v[0]);
        assert_eq!(islands.len(), 2);
        assert_eq!(islands[0], vec![0, 1]);
        assert_eq!(islands[1], vec![2]);
    }
    #[test]
    fn test_sleep_stats_compute() {
        let mut mgr = SleepManager::new(fast_config());
        let a = mgr.add_body(20, [0.0; 3]);
        let b = mgr.add_body(21, [1.0, 0.0, 0.0]);
        let c = mgr.add_body(22, [2.0, 0.0, 0.0]);
        mgr.bodies[a].state = SleepState::Active;
        mgr.bodies[b].state = SleepState::Drowsy { frames: 2 };
        mgr.bodies[c].state = SleepState::Sleeping;
        let stats = SleepStats::compute(&mgr);
        assert_eq!(stats.total, 3);
        assert_eq!(stats.active, 1);
        assert_eq!(stats.drowsy, 1);
        assert_eq!(stats.sleeping, 1);
        assert!((stats.sleep_ratio() - 1.0 / 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_island_cannot_sleep_fast_body() {
        let cfg = fast_config();
        let mut mgr = SleepManager::new(cfg.clone());
        let slow = mgr.add_body(30, [0.0; 3]);
        let fast = mgr.add_body(31, [1.0, 0.0, 0.0]);
        mgr.update_body_velocity(fast, [5.0, 0.0, 0.0], [0.0; 3]);
        let island = vec![slow, fast];
        assert!(!IslandSleepManager::island_can_sleep(
            &island,
            &mgr.bodies,
            &cfg
        ));
    }
    #[test]
    fn test_sleep_config_fast() {
        let cfg = SleepConfig::fast();
        assert!((cfg.linear_threshold - 0.05).abs() < 1e-10);
        assert_eq!(cfg.drowsy_frames, 10);
    }
    #[test]
    fn test_energy_based_sleep() {
        let cfg = SleepConfig::energy_based(0.01);
        let mut mgr = SleepManager::new(cfg);
        let idx = mgr.add_body_with_mass(1, [0.0; 3], 1.0, 1.0);
        for _ in 0..65 {
            mgr.tick();
        }
        assert_eq!(mgr.bodies[idx].state, SleepState::Sleeping);
    }
    #[test]
    fn test_energy_based_sleep_stays_active_with_velocity() {
        let cfg = SleepConfig::energy_based(0.01);
        let mut mgr = SleepManager::new(cfg);
        let idx = mgr.add_body_with_mass(2, [0.0; 3], 1.0, 1.0);
        mgr.update_body_velocity(idx, [5.0, 0.0, 0.0], [0.0; 3]);
        for _ in 0..100 {
            mgr.tick();
        }
        assert_eq!(mgr.bodies[idx].state, SleepState::Active);
    }
    #[test]
    fn test_timer_based_sleep() {
        let cfg = SleepConfig::timer_based(0.5);
        let mut mgr = SleepManager::new(cfg);
        let idx = mgr.add_body(1, [0.0; 3]);
        for _ in 0..6 {
            mgr.tick_timer(0.1);
        }
        assert_eq!(mgr.bodies[idx].state, SleepState::Sleeping);
    }
    #[test]
    fn test_timer_based_sleep_resets_on_motion() {
        let cfg = SleepConfig::timer_based(0.5);
        let mut mgr = SleepManager::new(cfg);
        let idx = mgr.add_body(1, [0.0; 3]);
        for _ in 0..3 {
            mgr.tick_timer(0.1);
        }
        mgr.update_body_velocity(idx, [5.0, 0.0, 0.0], [0.0; 3]);
        mgr.tick_timer(0.1);
        assert_eq!(mgr.bodies[idx].state, SleepState::Active);
        assert!((mgr.bodies[idx].sleep_timer - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_gradual_wake() {
        let mut cfg = fast_config();
        cfg.wake_damping = 0.5;
        let mut mgr = SleepManager::new(cfg);
        let idx = mgr.add_body(1, [0.0; 3]);
        mgr.update_body_velocity(idx, [10.0, 0.0, 0.0], [0.0; 3]);
        mgr.bodies[idx].state = SleepState::Sleeping;
        mgr.wake_gradual(idx);
        assert_eq!(mgr.bodies[idx].state, SleepState::Active);
        assert!((mgr.bodies[idx].linear_velocity[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_wake_nearby_gradual() {
        let mut cfg = fast_config();
        cfg.wake_damping = 0.25;
        let mut mgr = SleepManager::new(cfg);
        let idx = mgr.add_body(1, [1.0, 0.0, 0.0]);
        mgr.update_body_velocity(idx, [8.0, 0.0, 0.0], [4.0, 0.0, 0.0]);
        mgr.bodies[idx].state = SleepState::Sleeping;
        mgr.wake_nearby_gradual([0.0; 3], 5.0);
        assert_eq!(mgr.bodies[idx].state, SleepState::Active);
        assert!((mgr.bodies[idx].linear_velocity[0] - 2.0).abs() < 1e-10);
        assert!((mgr.bodies[idx].angular_velocity[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_sleeping_count() {
        let mut mgr = SleepManager::new(fast_config());
        mgr.add_body(1, [0.0; 3]);
        mgr.add_body(2, [1.0, 0.0, 0.0]);
        mgr.add_body(3, [2.0, 0.0, 0.0]);
        mgr.bodies[0].state = SleepState::Sleeping;
        mgr.bodies[2].state = SleepState::Sleeping;
        assert_eq!(mgr.sleeping_count(), 2);
        assert_eq!(mgr.active_count(), 1);
    }
    #[test]
    fn test_total_kinetic_energy() {
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body_with_mass(1, [0.0; 3], 2.0, 1.0);
        mgr.update_body_velocity(idx, [3.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let ke = mgr.total_kinetic_energy();
        assert!((ke - 9.0).abs() < 1e-10);
    }
    #[test]
    fn test_set_thresholds_at_runtime() {
        let mut mgr = SleepManager::new(fast_config());
        mgr.set_thresholds(0.5, 0.5);
        assert!((mgr.config.linear_threshold - 0.5).abs() < 1e-10);
        assert!((mgr.config.angular_threshold - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_sleep_stats_energy() {
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body_with_mass(1, [0.0; 3], 1.0, 1.0);
        mgr.update_body_velocity(idx, [2.0, 0.0, 0.0], [0.0; 3]);
        let stats = SleepStats::compute(&mgr);
        assert!((stats.total_energy - 2.0).abs() < 1e-10);
        assert!((stats.avg_energy - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_sleep_stats_active_ratio() {
        let mut mgr = SleepManager::new(fast_config());
        mgr.add_body(1, [0.0; 3]);
        mgr.add_body(2, [0.0; 3]);
        mgr.add_body(3, [0.0; 3]);
        mgr.add_body(4, [0.0; 3]);
        mgr.bodies[0].state = SleepState::Active;
        mgr.bodies[1].state = SleepState::Active;
        mgr.bodies[2].state = SleepState::Sleeping;
        mgr.bodies[3].state = SleepState::Sleeping;
        let stats = SleepStats::compute(&mgr);
        assert!((stats.active_ratio() - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_island_sleep_propagation_wake() {
        let mut island_mgr = IslandSleepManager::new(fast_config());
        island_mgr.manager.add_body(10, [0.0; 3]);
        island_mgr.manager.add_body(11, [1.0, 0.0, 0.0]);
        island_mgr.manager.add_body(12, [2.0, 0.0, 0.0]);
        for b in &mut island_mgr.manager.bodies {
            b.state = SleepState::Sleeping;
        }
        island_mgr
            .manager
            .update_body_velocity(1, [5.0, 0.0, 0.0], [0.0; 3]);
        let contacts = vec![(0, 1), (1, 2)];
        island_mgr.process_islands(&contacts);
        assert_eq!(island_mgr.manager.bodies[0].state, SleepState::Active);
        assert_eq!(island_mgr.manager.bodies[1].state, SleepState::Active);
        assert_eq!(island_mgr.manager.bodies[2].state, SleepState::Active);
    }
    #[test]
    fn test_island_sleep_propagation_sleep() {
        let mut island_mgr = IslandSleepManager::new(fast_config());
        island_mgr.manager.add_body(10, [0.0; 3]);
        island_mgr.manager.add_body(11, [1.0, 0.0, 0.0]);
        let contacts = vec![(0, 1)];
        island_mgr.process_islands(&contacts);
        assert_eq!(island_mgr.manager.bodies[0].state, SleepState::Sleeping);
        assert_eq!(island_mgr.manager.bodies[1].state, SleepState::Sleeping);
    }
    #[test]
    fn test_wake_island() {
        let mut island_mgr = IslandSleepManager::new(fast_config());
        island_mgr.manager.add_body(10, [0.0; 3]);
        island_mgr.manager.add_body(11, [1.0, 0.0, 0.0]);
        for b in &mut island_mgr.manager.bodies {
            b.state = SleepState::Sleeping;
        }
        island_mgr.wake_island(&[0, 1]);
        assert_eq!(island_mgr.manager.bodies[0].state, SleepState::Active);
        assert_eq!(island_mgr.manager.bodies[1].state, SleepState::Active);
    }
    #[test]
    fn test_body_kinetic_energy_self() {
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body_with_mass(1, [0.0; 3], 4.0, 2.0);
        mgr.update_body_velocity(idx, [1.0, 0.0, 0.0], [0.0, 0.0, 3.0]);
        let ke = mgr.bodies[idx].kinetic_energy_self();
        assert!((ke - 11.0).abs() < 1e-10);
    }
    #[test]
    fn test_is_below_threshold() {
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body(1, [0.0; 3]);
        assert!(mgr.bodies[idx].is_below_threshold(0.1, 0.1));
        mgr.update_body_velocity(idx, [5.0, 0.0, 0.0], [0.0; 3]);
        assert!(!mgr.bodies[idx].is_below_threshold(0.1, 0.1));
    }
    #[test]
    fn test_sleep_timer_reset_on_wake() {
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body(1, [0.0; 3]);
        mgr.bodies[idx].sleep_timer = 5.0;
        mgr.bodies[idx].state = SleepState::Sleeping;
        mgr.wake(idx);
        assert!((mgr.bodies[idx].sleep_timer - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_energy_island_can_sleep() {
        let cfg = SleepConfig::energy_based(0.01);
        let mut mgr = SleepManager::new(cfg.clone());
        let a = mgr.add_body_with_mass(1, [0.0; 3], 1.0, 1.0);
        let b = mgr.add_body_with_mass(2, [1.0, 0.0, 0.0], 1.0, 1.0);
        let island = vec![a, b];
        assert!(IslandSleepManager::island_can_sleep(
            &island,
            &mgr.bodies,
            &cfg
        ));
    }
    #[test]
    fn test_energy_island_cannot_sleep_with_velocity() {
        let cfg = SleepConfig::energy_based(0.01);
        let mut mgr = SleepManager::new(cfg.clone());
        let a = mgr.add_body_with_mass(1, [0.0; 3], 1.0, 1.0);
        let b = mgr.add_body_with_mass(2, [1.0, 0.0, 0.0], 1.0, 1.0);
        mgr.update_body_velocity(b, [5.0, 0.0, 0.0], [0.0; 3]);
        let island = vec![a, b];
        assert!(!IslandSleepManager::island_can_sleep(
            &island,
            &mgr.bodies,
            &cfg
        ));
    }
    #[test]
    fn test_predictor_already_below_threshold() {
        let p = SleepPredictor::new(1.0, 1.0, 0.1, 0.1);
        assert_eq!(p.time_to_linear_sleep(0.05), 0.0);
    }
    #[test]
    fn test_predictor_no_damping_infinite() {
        let p = SleepPredictor::new(0.0, 0.0, 0.1, 0.1);
        assert_eq!(p.time_to_linear_sleep(5.0), f64::INFINITY);
    }
    #[test]
    fn test_predictor_time_to_linear_sleep() {
        let p = SleepPredictor::new(1.0, 1.0, 1.0, 1.0);
        let t = p.time_to_linear_sleep(10.0);
        assert!((t - 10_f64.ln()).abs() < 1e-10, "t = {t}");
    }
    #[test]
    fn test_predictor_time_to_angular_sleep() {
        let p = SleepPredictor::new(1.0, 2.0, 1.0, 1.0);
        let t = p.time_to_angular_sleep(4.0);
        let expected = 4_f64.ln() / 2.0;
        assert!((t - expected).abs() < 1e-10, "t_ang = {t}");
    }
    #[test]
    fn test_predictor_time_to_sleep_takes_max() {
        let p = SleepPredictor::new(0.5, 2.0, 0.1, 0.1);
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body(1, [0.0; 3]);
        mgr.update_body_velocity(idx, [5.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let t = p.time_to_sleep(&mgr.bodies[idx]);
        let t_lin = p.time_to_linear_sleep(5.0);
        let t_ang = p.time_to_angular_sleep(1.0);
        assert!((t - t_lin.max(t_ang)).abs() < 1e-10, "max time = {t}");
    }
    #[test]
    fn test_predictor_predicted_speed_decays() {
        let p = SleepPredictor::new(1.0, 1.0, 0.1, 0.1);
        let v0 = 10.0;
        let v1 = p.predicted_linear_speed(v0, 1.0);
        assert!((v1 - v0 * (-1.0_f64).exp()).abs() < 1e-10, "v(1) = {v1}");
    }
    #[test]
    fn test_adaptive_threshold_initial_midpoint() {
        let a = AdaptiveSleepThreshold::new(0.01, 1.0, 0.1, 0.5);
        let expected = 0.505;
        assert!((a.current_threshold - expected).abs() < 0.001);
    }
    #[test]
    fn test_adaptive_threshold_raises_when_too_few_sleeping() {
        let mut a = AdaptiveSleepThreshold::new(0.01, 1.0, 0.1, 0.8);
        let before = a.current_threshold;
        a.update(0.2);
        assert!(
            a.current_threshold > before,
            "Should raise threshold when too few sleeping"
        );
    }
    #[test]
    fn test_adaptive_threshold_lowers_when_too_many_sleeping() {
        let mut a = AdaptiveSleepThreshold::new(0.01, 1.0, 0.1, 0.2);
        let before = a.current_threshold;
        a.update(0.9);
        assert!(
            a.current_threshold < before,
            "Should lower threshold when too many sleeping"
        );
    }
    #[test]
    fn test_adaptive_threshold_clamped() {
        let mut a = AdaptiveSleepThreshold::new(0.01, 1.0, 10.0, 0.0);
        for _ in 0..100 {
            a.update(0.0);
        }
        assert!(
            a.current_threshold <= 1.0,
            "Threshold should be clamped: {}",
            a.current_threshold
        );
    }
    #[test]
    fn test_adaptive_threshold_apply_to_config() {
        let a = AdaptiveSleepThreshold::new(0.1, 0.5, 0.01, 0.5);
        let mut cfg = SleepConfig::default();
        a.apply_to_config(&mut cfg);
        assert!((cfg.linear_threshold - a.current_threshold).abs() < 1e-14);
        assert!((cfg.angular_threshold - a.current_threshold).abs() < 1e-14);
    }
    #[test]
    fn test_sleep_wake_pattern_starts_active() {
        let p = SleepWakePattern::new(1.0, 0.5);
        assert!(p.is_active, "Should start active");
    }
    #[test]
    fn test_sleep_wake_pattern_transitions_to_rest() {
        let mut p = SleepWakePattern::new(1.0, 0.5);
        p.step(1.0);
        assert!(
            !p.is_active,
            "Should transition to rest after active_duration"
        );
    }
    #[test]
    fn test_sleep_wake_pattern_transitions_back_to_active() {
        let mut p = SleepWakePattern::new(1.0, 0.5);
        p.step(1.0);
        p.step(0.5);
        assert!(p.is_active, "Should return to active after rest_duration");
    }
    #[test]
    fn test_sleep_wake_pattern_duty_cycle() {
        let p = SleepWakePattern::new(3.0, 1.0);
        assert!(
            (p.duty_cycle() - 0.75).abs() < 1e-10,
            "Duty cycle = {}",
            p.duty_cycle()
        );
    }
    #[test]
    fn test_sleep_wake_pattern_apply_wakes_sleeping() {
        let p = SleepWakePattern::new(1.0, 0.5);
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body(1, [0.0; 3]);
        mgr.bodies[idx].state = SleepState::Sleeping;
        p.apply_to_manager(&mut mgr, &[idx]);
        assert_eq!(mgr.bodies[idx].state, SleepState::Active);
    }
    #[test]
    fn test_sleep_wake_pattern_rest_does_not_force_wake() {
        let mut p = SleepWakePattern::new(1.0, 0.5);
        p.step(1.0);
        let mut mgr = SleepManager::new(fast_config());
        let idx = mgr.add_body(1, [0.0; 3]);
        mgr.bodies[idx].state = SleepState::Sleeping;
        p.apply_to_manager(&mut mgr, &[idx]);
        assert_eq!(mgr.bodies[idx].state, SleepState::Sleeping);
    }
    #[test]
    fn test_propagate_wake_wakes_island_members() {
        let mut island_mgr = IslandSleepManager::new(fast_config());
        island_mgr.manager.add_body(1, [0.0; 3]);
        island_mgr.manager.add_body(2, [1.0, 0.0, 0.0]);
        island_mgr.manager.add_body(3, [2.0, 0.0, 0.0]);
        island_mgr.manager.bodies[1].state = SleepState::Sleeping;
        island_mgr.manager.bodies[2].state = SleepState::Sleeping;
        let contacts = vec![(0_usize, 1_usize), (1_usize, 2_usize)];
        island_mgr.propagate_wake(&contacts);
        assert_eq!(island_mgr.manager.bodies[1].state, SleepState::Active);
        assert_eq!(island_mgr.manager.bodies[2].state, SleepState::Active);
    }
    #[test]
    fn test_propagate_wake_isolated_bodies_stay_sleeping() {
        let mut island_mgr = IslandSleepManager::new(fast_config());
        island_mgr.manager.add_body(1, [0.0; 3]);
        island_mgr.manager.add_body(2, [100.0, 0.0, 0.0]);
        island_mgr.manager.bodies[1].state = SleepState::Sleeping;
        island_mgr.manager.bodies[0].state = SleepState::Active;
        island_mgr.propagate_wake(&[]);
        assert_eq!(island_mgr.manager.bodies[1].state, SleepState::Sleeping);
    }
    #[test]
    fn test_island_count() {
        let mut island_mgr = IslandSleepManager::new(fast_config());
        island_mgr.manager.add_body(1, [0.0; 3]);
        island_mgr.manager.add_body(2, [1.0, 0.0, 0.0]);
        island_mgr.manager.add_body(3, [10.0, 0.0, 0.0]);
        let contacts = vec![(0_usize, 1_usize)];
        let count = island_mgr.island_count(&contacts);
        assert_eq!(count, 2, "Should have 2 islands");
    }
    #[test]
    fn test_sleeping_and_active_indices() {
        let mut island_mgr = IslandSleepManager::new(fast_config());
        island_mgr.manager.add_body(1, [0.0; 3]);
        island_mgr.manager.add_body(2, [1.0, 0.0, 0.0]);
        island_mgr.manager.bodies[0].state = SleepState::Sleeping;
        island_mgr.manager.bodies[1].state = SleepState::Active;
        let sleeping = island_mgr.sleeping_indices();
        let active = island_mgr.active_indices();
        assert_eq!(sleeping, vec![0]);
        assert_eq!(active, vec![1]);
    }
    #[test]
    fn test_force_all_to_sleep() {
        let mut island_mgr = IslandSleepManager::new(fast_config());
        island_mgr.manager.add_body(1, [0.0; 3]);
        island_mgr.manager.add_body(2, [1.0, 0.0, 0.0]);
        island_mgr.force_all_to_sleep();
        for b in &island_mgr.manager.bodies {
            assert_eq!(b.state, SleepState::Sleeping);
        }
    }
    #[test]
    fn test_tick_islands_propagates_wake() {
        let mut island_mgr = IslandSleepManager::new(fast_config());
        island_mgr.manager.add_body(1, [0.0; 3]);
        island_mgr.manager.add_body(2, [1.0, 0.0, 0.0]);
        island_mgr.manager.bodies[0].state = SleepState::Sleeping;
        island_mgr.manager.bodies[1].state = SleepState::Sleeping;
        island_mgr.manager.bodies[0].state = SleepState::Active;
        island_mgr.tick_islands(&[(0, 1)]);
        assert_eq!(island_mgr.manager.bodies[1].state, SleepState::Active);
    }
    #[test]
    fn test_sleep_history_record_and_count() {
        let mut h = SleepHistory::new();
        h.record_sleep(1.0, 1);
        h.record_sleep(2.0, 2);
        h.record_wake(3.0, 1);
        assert_eq!(h.sleep_count(), 2);
        assert_eq!(h.wake_count(), 1);
    }
    #[test]
    fn test_sleep_history_events_for_body() {
        let mut h = SleepHistory::new();
        h.record_sleep(1.0, 1);
        h.record_wake(2.0, 1);
        h.record_sleep(3.0, 2);
        let evs = h.events_for(1);
        assert_eq!(evs.len(), 2);
    }
    #[test]
    fn test_sleep_history_clear() {
        let mut h = SleepHistory::new();
        h.record_sleep(1.0, 1);
        h.clear();
        assert_eq!(h.events.len(), 0);
    }
    #[test]
    fn test_sleep_history_empty() {
        let h = SleepHistory::new();
        assert_eq!(h.sleep_count(), 0);
        assert_eq!(h.wake_count(), 0);
        assert_eq!(h.events_for(1).len(), 0);
    }
    #[test]
    fn test_energy_tracker_below_threshold_after_low_values() {
        let mut tracker = EnergyTracker::new(5, 1.0);
        for _ in 0..5 {
            tracker.push(0.001);
        }
        assert!(
            tracker.is_below_threshold(),
            "all values below threshold, should return true"
        );
    }
    #[test]
    fn test_energy_tracker_not_below_threshold_with_high_value() {
        let mut tracker = EnergyTracker::new(3, 1.0);
        tracker.push(0.1);
        tracker.push(5.0);
        tracker.push(0.1);
        assert!(!tracker.is_below_threshold());
    }
    #[test]
    fn test_energy_tracker_mean_energy() {
        let mut tracker = EnergyTracker::new(4, 100.0);
        tracker.push(2.0);
        tracker.push(4.0);
        let mean = tracker.mean_energy();
        assert!((mean - 3.0).abs() < 1e-10, "mean: {mean}");
    }
    #[test]
    fn test_sleeping_system_wake_resets_tracker() {
        let mut sys = SleepingSystem::new(3, 10.0);
        for _ in 0..3 {
            sys.update_body(1, 0.001, 0.0);
        }
        sys.wake_body(1);
        let tracker = sys.trackers.get(&1).unwrap();
        assert_eq!(tracker.count, 0, "count should be reset to 0 after wake");
    }
    #[test]
    fn test_sleeping_system_sleep_delay() {
        let mut sys = SleepingSystem::new(4, 1.0);
        for _ in 0..3 {
            sys.update_body(42, 0.0, 0.0);
        }
        assert!(!sys.should_sleep(42), "not enough history to sleep yet");
        sys.update_body(42, 0.0, 0.0);
        assert!(
            sys.should_sleep(42),
            "should sleep after full history of low energy"
        );
    }
    #[test]
    fn test_kinetic_energy_rigid_formula() {
        let ke = kinetic_energy_rigid(2.0, [3.0, 0.0, 0.0], [1.0, 2.0, 3.0], [0.0, 1.0, 0.0]);
        assert!((ke - 10.0).abs() < 1e-10, "KE: {ke}");
    }
    #[test]
    fn test_potential_energy_rigid_formula() {
        let pe = potential_energy_rigid(5.0, [0.0, 3.0, 0.0], 9.81);
        let expected = 5.0 * 9.81 * 3.0;
        assert!((pe - expected).abs() < 1e-10, "PE: {pe}");
    }
}
#[cfg(test)]
mod tests_sleeping_extra {

    use crate::ActivityLevel;

    use crate::SleepConfig;
    use crate::SleepDurationTracker;

    use crate::SleepHysteresis;
    use crate::SleepManager;

    use crate::SleepState;

    use crate::VelocityDamper;
    use crate::WakePropagator;
    fn fast_cfg() -> SleepConfig {
        SleepConfig {
            linear_threshold: 0.1,
            angular_threshold: 0.1,
            drowsy_frames: 3,
            energy_threshold: 0.0,
            wake_damping: 1.0,
            sleep_time_seconds: 0.0,
        }
    }
    #[test]
    fn test_activity_level_resting() {
        assert_eq!(ActivityLevel::from_speed(0.01), ActivityLevel::Resting);
    }
    #[test]
    fn test_activity_level_minimal() {
        assert_eq!(ActivityLevel::from_speed(0.1), ActivityLevel::Minimal);
    }
    #[test]
    fn test_activity_level_moderate() {
        assert_eq!(ActivityLevel::from_speed(1.0), ActivityLevel::Moderate);
    }
    #[test]
    fn test_activity_level_high() {
        assert_eq!(ActivityLevel::from_speed(5.0), ActivityLevel::High);
    }
    #[test]
    fn test_activity_level_needs_simulation_resting_false() {
        assert!(!ActivityLevel::Resting.needs_simulation());
    }
    #[test]
    fn test_activity_level_needs_simulation_active_true() {
        assert!(ActivityLevel::High.needs_simulation());
        assert!(ActivityLevel::Moderate.needs_simulation());
        assert!(ActivityLevel::Minimal.needs_simulation());
    }
    #[test]
    fn test_activity_level_from_kinetic_energy() {
        let level = ActivityLevel::from_kinetic_energy(0.001, 0.01, 1.0, 10.0);
        assert_eq!(level, ActivityLevel::Resting);
        let level2 = ActivityLevel::from_kinetic_energy(0.5, 0.01, 1.0, 10.0);
        assert_eq!(level2, ActivityLevel::Minimal);
        let level3 = ActivityLevel::from_kinetic_energy(5.0, 0.01, 1.0, 10.0);
        assert_eq!(level3, ActivityLevel::Moderate);
        let level4 = ActivityLevel::from_kinetic_energy(50.0, 0.01, 1.0, 10.0);
        assert_eq!(level4, ActivityLevel::High);
    }
    #[test]
    fn test_velocity_damper_apply_reduces_speed() {
        let damper = VelocityDamper::new(0.9);
        let (lin, ang) = damper.apply([10.0, 0.0, 0.0], [0.0, 5.0, 0.0]);
        assert!((lin[0] - 9.0).abs() < 1e-10, "lin[0] = {}", lin[0]);
        assert!((ang[1] - 4.5).abs() < 1e-10, "ang[1] = {}", ang[1]);
    }
    #[test]
    fn test_velocity_damper_identity_no_change() {
        let damper = VelocityDamper::identity();
        let (lin, ang) = damper.apply([10.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
        assert!((lin[0] - 10.0).abs() < 1e-10);
        assert!((ang[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_velocity_damper_only_affects_drowsy() {
        let damper = VelocityDamper::new(0.5);
        let mut mgr = SleepManager::new(fast_cfg());
        let idx_drowsy = mgr.add_body(1, [0.0; 3]);
        let idx_active = mgr.add_body(2, [1.0, 0.0, 0.0]);
        mgr.update_body_velocity(idx_drowsy, [10.0, 0.0, 0.0], [0.0; 3]);
        mgr.update_body_velocity(idx_active, [10.0, 0.0, 0.0], [0.0; 3]);
        mgr.bodies[idx_drowsy].state = SleepState::Drowsy { frames: 1 };
        damper.apply_to_manager(&mut mgr);
        assert!(
            (mgr.bodies[idx_drowsy].linear_velocity[0] - 5.0).abs() < 1e-10,
            "drowsy body should be damped: {}",
            mgr.bodies[idx_drowsy].linear_velocity[0]
        );
        assert!(
            (mgr.bodies[idx_active].linear_velocity[0] - 10.0).abs() < 1e-10,
            "active body should not be damped: {}",
            mgr.bodies[idx_active].linear_velocity[0]
        );
    }
    #[test]
    fn test_velocity_damper_clamps_factor() {
        let damper = VelocityDamper::new(1.5);
        assert!((damper.frame_damping - 1.0).abs() < 1e-10);
        let damper2 = VelocityDamper::new(-0.5);
        assert!((damper2.frame_damping).abs() < 1e-10);
    }
    #[test]
    fn test_hysteresis_should_sleep_below_threshold() {
        let h = SleepHysteresis::new(0.1, 0.3);
        assert!(h.should_sleep(0.05), "below sleep_threshold should sleep");
        assert!(
            !h.should_sleep(0.2),
            "above sleep_threshold should not sleep"
        );
    }
    #[test]
    fn test_hysteresis_should_wake_above_threshold() {
        let h = SleepHysteresis::new(0.1, 0.3);
        assert!(h.should_wake(0.5), "above wake_threshold should wake");
        assert!(!h.should_wake(0.2), "below wake_threshold should not wake");
    }
    #[test]
    fn test_hysteresis_band_prevents_oscillation() {
        let h = SleepHysteresis::new(0.1, 0.3);
        assert!(!h.should_sleep(0.2), "in band: should not sleep");
        assert!(!h.should_wake(0.2), "in band: should not wake");
    }
    #[test]
    fn test_hysteresis_update_active_to_drowsy() {
        let h = SleepHysteresis::new(0.1, 0.3);
        let new_state = h.update_state(&SleepState::Active, 0.05);
        assert!(
            matches!(new_state, SleepState::Drowsy { .. }),
            "should become drowsy: {:?}",
            new_state
        );
    }
    #[test]
    fn test_hysteresis_update_sleeping_to_active() {
        let h = SleepHysteresis::new(0.1, 0.3);
        let new_state = h.update_state(&SleepState::Sleeping, 0.5);
        assert_eq!(
            new_state,
            SleepState::Active,
            "should wake: {:?}",
            new_state
        );
    }
    #[test]
    fn test_hysteresis_update_sleeping_stays_sleeping_in_band() {
        let h = SleepHysteresis::new(0.1, 0.3);
        let new_state = h.update_state(&SleepState::Sleeping, 0.2);
        assert_eq!(
            new_state,
            SleepState::Sleeping,
            "should stay sleeping in band: {:?}",
            new_state
        );
    }
    #[test]
    fn test_hysteresis_band_width() {
        let h = SleepHysteresis::new(0.1, 0.4);
        assert!((h.band_width() - 0.3).abs() < 1e-12);
    }
    #[test]
    fn test_wake_propagator_wakes_adjacent_sleeping() {
        let prop = WakePropagator::new(5);
        let mut mgr = SleepManager::new(fast_cfg());
        mgr.add_body(1, [0.0; 3]);
        mgr.add_body(2, [1.0, 0.0, 0.0]);
        mgr.bodies[0].state = SleepState::Active;
        mgr.bodies[1].state = SleepState::Sleeping;
        let woken = prop.propagate(&mut mgr, &[(0, 1)]);
        assert_eq!(woken, 1, "one body should be woken");
        assert_eq!(mgr.bodies[1].state, SleepState::Active);
    }
    #[test]
    fn test_wake_propagator_chain_propagation() {
        let prop = WakePropagator::new(5);
        let mut mgr = SleepManager::new(fast_cfg());
        mgr.add_body(1, [0.0; 3]);
        mgr.add_body(2, [1.0, 0.0, 0.0]);
        mgr.add_body(3, [2.0, 0.0, 0.0]);
        mgr.bodies[0].state = SleepState::Active;
        mgr.bodies[1].state = SleepState::Sleeping;
        mgr.bodies[2].state = SleepState::Sleeping;
        let contacts = vec![(0usize, 1usize), (1usize, 2usize)];
        prop.propagate(&mut mgr, &contacts);
        assert_eq!(mgr.bodies[1].state, SleepState::Active);
        assert_eq!(mgr.bodies[2].state, SleepState::Active);
    }
    #[test]
    fn test_wake_propagator_does_not_wake_disconnected() {
        let prop = WakePropagator::new(5);
        let mut mgr = SleepManager::new(fast_cfg());
        mgr.add_body(1, [0.0; 3]);
        mgr.add_body(2, [100.0, 0.0, 0.0]);
        mgr.bodies[0].state = SleepState::Active;
        mgr.bodies[1].state = SleepState::Sleeping;
        prop.propagate(&mut mgr, &[]);
        assert_eq!(mgr.bodies[1].state, SleepState::Sleeping);
    }
    #[test]
    fn test_wake_propagator_depth_limit() {
        let prop = WakePropagator::new(1);
        let mut mgr = SleepManager::new(fast_cfg());
        mgr.add_body(1, [0.0; 3]);
        mgr.add_body(2, [1.0, 0.0, 0.0]);
        mgr.add_body(3, [2.0, 0.0, 0.0]);
        mgr.bodies[0].state = SleepState::Active;
        mgr.bodies[1].state = SleepState::Sleeping;
        mgr.bodies[2].state = SleepState::Sleeping;
        let contacts = vec![(0usize, 1usize), (1usize, 2usize)];
        prop.propagate(&mut mgr, &contacts);
        assert_eq!(mgr.bodies[1].state, SleepState::Active);
    }
    #[test]
    fn test_duration_tracker_basic_accumulation() {
        let mut tracker = SleepDurationTracker::new();
        tracker.record_sleep_start(1, 0.0);
        tracker.record_wake(1, 5.0);
        assert!((tracker.total_sleep_time(1) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_duration_tracker_multiple_sleep_cycles() {
        let mut tracker = SleepDurationTracker::new();
        tracker.record_sleep_start(1, 0.0);
        tracker.record_wake(1, 3.0);
        tracker.record_sleep_start(1, 4.0);
        tracker.record_wake(1, 7.0);
        assert!((tracker.total_sleep_time(1) - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_duration_tracker_currently_sleeping() {
        let mut tracker = SleepDurationTracker::new();
        tracker.record_sleep_start(2, 1.0);
        assert!(tracker.is_currently_sleeping(2));
        tracker.record_wake(2, 4.0);
        assert!(!tracker.is_currently_sleeping(2));
    }
    #[test]
    fn test_duration_tracker_unknown_body_returns_zero() {
        let tracker = SleepDurationTracker::new();
        assert!((tracker.total_sleep_time(99) - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_duration_tracker_current_sleep_duration() {
        let mut tracker = SleepDurationTracker::new();
        tracker.record_sleep_start(3, 2.0);
        let dur = tracker.current_sleep_duration(3, 5.0);
        assert!((dur - 3.0).abs() < 1e-10, "current duration: {dur}");
    }
    #[test]
    fn test_duration_tracker_reset() {
        let mut tracker = SleepDurationTracker::new();
        tracker.record_sleep_start(1, 0.0);
        tracker.record_wake(1, 5.0);
        tracker.reset();
        assert!((tracker.total_sleep_time(1) - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_duration_tracker_wake_without_start_is_noop() {
        let mut tracker = SleepDurationTracker::new();
        tracker.record_wake(5, 10.0);
        assert!((tracker.total_sleep_time(5) - 0.0).abs() < 1e-12);
    }
}
#[cfg(test)]
mod new_sleep_tests {

    use crate::EnergyDecayTracker;

    use crate::FrozenBody;
    use crate::IslandSleepDecision;
    use crate::IslandSleepEvaluator;

    use crate::PseudoSleepBody;

    use crate::SleepThresholdAdapter;

    #[test]
    fn test_frozen_body_new_is_frozen() {
        let b = FrozenBody::new(1, [0.0, 0.0, 0.0], 5.0);
        assert!(b.frozen);
    }
    #[test]
    fn test_frozen_body_inv_mass_zero_when_frozen() {
        let b = FrozenBody::new(1, [0.0; 3], 5.0);
        assert_eq!(b.inv_mass(), 0.0);
    }
    #[test]
    fn test_frozen_body_inv_mass_positive_when_unfrozen() {
        let mut b = FrozenBody::new(1, [0.0; 3], 5.0);
        b.unfreeze();
        assert!(b.inv_mass() > 0.0, "inv_mass={}", b.inv_mass());
    }
    #[test]
    fn test_frozen_body_freeze_unfreeze_cycle() {
        let mut b = FrozenBody::new(2, [1.0, 2.0, 3.0], 10.0);
        b.unfreeze();
        assert!(!b.frozen);
        b.freeze();
        assert!(b.frozen);
    }
    #[test]
    fn test_frozen_body_set_position() {
        let mut b = FrozenBody::new(1, [0.0; 3], 5.0);
        b.set_position([1.0, 2.0, 3.0]);
        assert_eq!(b.position, [1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_frozen_body_set_orientation() {
        let mut b = FrozenBody::new(1, [0.0; 3], 5.0);
        let q = [0.707, 0.0, 0.707, 0.0];
        b.set_orientation(q);
        assert_eq!(b.orientation, q);
    }
    #[test]
    fn test_pseudo_sleep_disabled_always_updates() {
        let mut b = PseudoSleepBody::new(1, 0.016, 0.25);
        b.disable();
        for _ in 0..10 {
            let result = b.tick();
            assert!(result.is_some(), "disabled pseudo-sleep should always tick");
        }
    }
    #[test]
    fn test_pseudo_sleep_rate_fraction_limits_updates() {
        let mut b = PseudoSleepBody::new(1, 0.016, 0.25);
        let mut update_count = 0;
        for _ in 0..8 {
            if b.tick().is_some() {
                update_count += 1;
            }
        }
        assert_eq!(update_count, 2, "update_count={update_count}");
    }
    #[test]
    fn test_pseudo_sleep_effective_dt_larger_than_full() {
        let b = PseudoSleepBody::new(1, 0.016, 0.25);
        assert!(
            (b.effective_dt() - 0.064).abs() < 1e-10,
            "eff_dt={}",
            b.effective_dt()
        );
    }
    #[test]
    fn test_pseudo_sleep_enable_disable() {
        let mut b = PseudoSleepBody::new(1, 0.016, 0.25);
        b.disable();
        assert!(!b.active);
        b.enable();
        assert!(b.active);
    }
    #[test]
    fn test_threshold_adapter_is_still_below_thresholds() {
        let adapter = SleepThresholdAdapter::new(0.1, 0.1, 0.001, 1.0);
        assert!(adapter.is_still(0.05, 0.05));
    }
    #[test]
    fn test_threshold_adapter_not_still_above_threshold() {
        let adapter = SleepThresholdAdapter::new(0.1, 0.1, 0.001, 1.0);
        assert!(!adapter.is_still(0.5, 0.05));
    }
    #[test]
    fn test_threshold_adapter_update_stays_in_bounds() {
        let mut adapter = SleepThresholdAdapter::new(0.1, 0.1, 0.001, 1.0);
        for _ in 0..1000 {
            adapter.update(1000.0);
        }
        assert!(adapter.linear_threshold <= 1.0 + 1e-9);
        assert!(adapter.linear_threshold >= 0.001 - 1e-9);
    }
    #[test]
    fn test_threshold_adapter_update_quiet_scene() {
        let mut adapter = SleepThresholdAdapter::new(0.5, 0.5, 0.001, 1.0);
        for _ in 0..1000 {
            adapter.update(0.0);
        }
        assert!(
            adapter.linear_threshold <= 0.5 + 1e-9,
            "threshold should not exceed initial for quiet scene"
        );
    }
    #[test]
    fn test_island_evaluator_all_sleeping() {
        let eval = IslandSleepEvaluator::new(0.01);
        let result = eval.evaluate(&[0.001, 0.002, 0.003]);
        assert_eq!(result, IslandSleepDecision::SleepAll);
    }
    #[test]
    fn test_island_evaluator_all_active() {
        let eval = IslandSleepEvaluator::new(0.01);
        let result = eval.evaluate(&[1.0, 2.0, 3.0]);
        assert_eq!(result, IslandSleepDecision::WakeAll);
    }
    #[test]
    fn test_island_evaluator_mixed() {
        let eval = IslandSleepEvaluator::new(0.01);
        let result = eval.evaluate(&[0.001, 2.0]);
        assert_eq!(result, IslandSleepDecision::Mixed);
    }
    #[test]
    fn test_island_evaluator_empty_island_sleeps() {
        let eval = IslandSleepEvaluator::new(0.01);
        assert_eq!(eval.evaluate(&[]), IslandSleepDecision::SleepAll);
    }
    #[test]
    fn test_island_evaluator_count_candidates() {
        let eval = IslandSleepEvaluator::new(0.01);
        let count = eval.count_sleeping_candidates(&[0.001, 0.005, 1.0, 2.0]);
        assert_eq!(count, 2);
    }
    #[test]
    fn test_decay_tracker_is_decaying() {
        let mut t = EnergyDecayTracker::new(1);
        t.update(1.0, 0.01);
        t.update(0.5, 0.01);
        assert!(t.is_decaying());
    }
    #[test]
    fn test_decay_tracker_not_decaying_when_increasing() {
        let mut t = EnergyDecayTracker::new(1);
        t.update(0.5, 0.01);
        t.update(1.0, 0.01);
        assert!(!t.is_decaying());
    }
    #[test]
    fn test_decay_tracker_time_to_threshold_none_if_already_below() {
        let mut t = EnergyDecayTracker::new(1);
        t.update(0.001, 0.01);
        assert!(t.time_to_threshold(0.01).is_none());
    }
    #[test]
    fn test_decay_tracker_time_to_threshold_positive() {
        let mut t = EnergyDecayTracker::new(1);
        t.update(1.0, 0.1);
        t.update(0.9, 0.1);
        let tte = t.time_to_threshold(0.01);
        if let Some(v) = tte {
            assert!(v > 0.0, "time to threshold should be positive: {v}");
        }
    }
    #[test]
    fn test_decay_tracker_tau_initially_infinite() {
        let t = EnergyDecayTracker::new(1);
        assert!(t.tau.is_infinite());
    }
    #[test]
    fn test_decay_tracker_update_sets_current_energy() {
        let mut t = EnergyDecayTracker::new(1);
        t.update(3.125, 0.01);
        assert!((t.current_energy - 3.125).abs() < 1e-12);
    }
}
