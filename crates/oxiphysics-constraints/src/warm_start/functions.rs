//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{BlendMode, CacheSelectionPolicy, WarmStartCache, WarmStartCandidate};

/// Apply a linear and angular impulse to a body.
///
/// Updates `vel` and `ang_vel` in-place:
/// - `vel += impulse * inv_mass`
/// - `ang_vel += inv_inertia * (r × impulse)`
#[allow(dead_code)]
pub fn apply_impulse(
    vel: &mut [f64; 3],
    ang_vel: &mut [f64; 3],
    r: [f64; 3],
    impulse: [f64; 3],
    inv_mass: f64,
    inv_inertia: f64,
) {
    vel[0] += impulse[0] * inv_mass;
    vel[1] += impulse[1] * inv_mass;
    vel[2] += impulse[2] * inv_mass;
    let torque = cross3(r, impulse);
    ang_vel[0] += torque[0] * inv_inertia;
    ang_vel[1] += torque[1] * inv_inertia;
    ang_vel[2] += torque[2] * inv_inertia;
}
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[allow(dead_code)]
pub(super) fn dist_sq_3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::warm_start::AdaptiveWarmStartScaler;

    use crate::warm_start::ConstraintPairMatcher;
    use crate::warm_start::ContactFingerprint;

    use crate::warm_start::ContactVelocitySolver;
    use crate::warm_start::ImpulseAging;

    use crate::warm_start::WarmStartMap;

    use crate::warm_start::WarmStartQualityMetrics;

    use crate::warm_start::WarmStartStrategy;
    #[test]
    fn test_warm_start_cache_magnitude() {
        let c = WarmStartCache::with_impulses(3.0, 4.0, 0.0);
        assert!((c.magnitude() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_warm_start_cache_scale() {
        let c = WarmStartCache::with_impulses(10.0, 20.0, 30.0);
        let s = c.scale(0.5);
        assert!((s.lambda_n - 5.0).abs() < 1e-12);
        assert!((s.lambda_t1 - 10.0).abs() < 1e-12);
        assert!((s.lambda_t2 - 15.0).abs() < 1e-12);
    }
    #[test]
    fn test_warm_start_cache_is_negligible() {
        let c = WarmStartCache::with_impulses(1e-15, 1e-15, 1e-15);
        assert!(c.is_negligible(1e-10));
        let c2 = WarmStartCache::with_impulses(1.0, 0.0, 0.0);
        assert!(!c2.is_negligible(1e-10));
    }
    #[test]
    fn test_warm_start_cache_lerp() {
        let a = WarmStartCache::with_impulses(0.0, 0.0, 0.0);
        let b = WarmStartCache::with_impulses(10.0, 20.0, 30.0);
        let mid = a.lerp(&b, 0.5);
        assert!((mid.lambda_n - 5.0).abs() < 1e-12);
        assert!((mid.lambda_t1 - 10.0).abs() < 1e-12);
        assert!((mid.lambda_t2 - 15.0).abs() < 1e-12);
    }
    #[test]
    fn test_warm_start_cache_lerp_clamped() {
        let a = WarmStartCache::with_impulses(0.0, 0.0, 0.0);
        let b = WarmStartCache::with_impulses(10.0, 10.0, 10.0);
        let result = a.lerp(&b, 2.0);
        assert!((result.lambda_n - 10.0).abs() < 1e-12);
    }
    /// WarmStartMap key ordering: (5,3) and (3,5) must map to the same key.
    #[test]
    fn test_warm_start_map_key_ordering() {
        let k1 = WarmStartMap::key(5, 3);
        let k2 = WarmStartMap::key(3, 5);
        assert_eq!(
            k1, k2,
            "Keys must be identical regardless of argument order"
        );
        let mut map = WarmStartMap::new();
        map.store(
            5,
            3,
            WarmStartCache {
                lambda_n: 42.0,
                lambda_t1: 1.0,
                lambda_t2: 2.0,
            },
        );
        let retrieved = map.get_or_default(3, 5);
        assert!(
            (retrieved.lambda_n - 42.0).abs() < 1e-12,
            "Should retrieve the same entry via reversed IDs"
        );
    }
    #[test]
    fn test_warm_start_map_remove() {
        let mut map = WarmStartMap::new();
        map.store(1, 2, WarmStartCache::with_impulses(5.0, 0.0, 0.0));
        assert_eq!(map.len(), 1);
        let removed = map.remove(2, 1);
        assert!(removed.is_some());
        assert!((removed.unwrap().lambda_n - 5.0).abs() < 1e-12);
        assert!(map.is_empty());
    }
    #[test]
    fn test_warm_start_map_scale_all() {
        let mut map = WarmStartMap::new();
        map.store(1, 2, WarmStartCache::with_impulses(10.0, 20.0, 30.0));
        map.store(3, 4, WarmStartCache::with_impulses(5.0, 10.0, 15.0));
        map.scale_all(0.5);
        assert!((map.get_or_default(1, 2).lambda_n - 5.0).abs() < 1e-12);
        assert!((map.get_or_default(3, 4).lambda_t2 - 7.5).abs() < 1e-12);
    }
    #[test]
    fn test_warm_start_map_prune() {
        let mut map = WarmStartMap::new();
        map.store(1, 2, WarmStartCache::with_impulses(10.0, 0.0, 0.0));
        map.store(3, 4, WarmStartCache::with_impulses(0.001, 0.0, 0.0));
        map.prune(0.01);
        assert_eq!(map.len(), 1);
        assert!(map.get(1, 2).is_some());
        assert!(map.get(3, 4).is_none());
    }
    #[test]
    fn test_warm_start_map_stats() {
        let mut map = WarmStartMap::new();
        map.store(1, 2, WarmStartCache::with_impulses(3.0, 4.0, 0.0));
        map.store(3, 4, WarmStartCache::with_impulses(0.0, 0.0, 10.0));
        let (total, max_mag) = map.stats();
        assert!((total - 15.0).abs() < 1e-12);
        assert!((max_mag - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_warm_start_map_apply_strategy_none() {
        let mut map = WarmStartMap::new();
        map.store(1, 2, WarmStartCache::with_impulses(10.0, 0.0, 0.0));
        map.apply_strategy(WarmStartStrategy::None);
        assert!(map.is_empty());
    }
    #[test]
    fn test_warm_start_map_apply_strategy_decayed() {
        let mut map = WarmStartMap::new();
        map.store(1, 2, WarmStartCache::with_impulses(10.0, 0.0, 0.0));
        map.apply_strategy(WarmStartStrategy::Decayed(0.8));
        assert!((map.get_or_default(1, 2).lambda_n - 8.0).abs() < 1e-12);
    }
    #[test]
    fn test_impulse_aging_basic() {
        let mut aging = ImpulseAging::new(3, 0.9);
        aging.touch(1, 2);
        assert_eq!(aging.age_of(1, 2), Some(0));
        assert_eq!(aging.len(), 1);
        aging.tick();
        assert_eq!(aging.age_of(1, 2), Some(1));
        aging.tick();
        assert_eq!(aging.age_of(1, 2), Some(2));
        aging.tick();
        assert_eq!(aging.age_of(1, 2), Some(3));
        let evicted = aging.tick();
        assert_eq!(evicted.len(), 1);
        assert!(aging.is_empty());
    }
    #[test]
    fn test_impulse_aging_decay_multiplier() {
        let mut aging = ImpulseAging::new(10, 0.8);
        aging.touch(1, 2);
        assert!((aging.decay_multiplier(1, 2) - 1.0).abs() < 1e-12);
        aging.tick();
        assert!((aging.decay_multiplier(1, 2) - 0.8).abs() < 1e-12);
        aging.tick();
        assert!((aging.decay_multiplier(1, 2) - 0.64).abs() < 1e-12);
    }
    #[test]
    fn test_impulse_aging_touch_resets() {
        let mut aging = ImpulseAging::new(3, 0.9);
        aging.touch(1, 2);
        aging.tick();
        aging.tick();
        aging.touch(1, 2);
        assert_eq!(aging.age_of(1, 2), Some(0));
    }
    #[test]
    fn test_impulse_aging_apply_to_map() {
        let mut aging = ImpulseAging::new(2, 0.5);
        let mut map = WarmStartMap::new();
        map.store(1, 2, WarmStartCache::with_impulses(100.0, 0.0, 0.0));
        map.store(3, 4, WarmStartCache::with_impulses(50.0, 0.0, 0.0));
        aging.touch(1, 2);
        aging.touch(3, 4);
        aging.apply_to_map(&mut map);
        assert!((map.get_or_default(1, 2).lambda_n - 50.0).abs() < 1e-12);
        assert!((map.get_or_default(3, 4).lambda_n - 25.0).abs() < 1e-12);
    }
    #[test]
    fn test_quality_metrics_hit_rate() {
        let mut m = WarmStartQualityMetrics::new();
        m.record_hit(10.0, 10.0);
        m.record_hit(5.0, 5.0);
        m.record_miss();
        assert!((m.hit_rate() - 2.0 / 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_quality_metrics_perfect() {
        let mut m = WarmStartQualityMetrics::new();
        m.record_hit(10.0, 10.0);
        m.record_hit(5.0, 5.0);
        assert!((m.mean_relative_error()).abs() < 1e-12);
        assert!((m.quality_score() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_quality_metrics_empty() {
        let m = WarmStartQualityMetrics::new();
        assert!((m.hit_rate()).abs() < 1e-12);
        assert!((m.quality_score()).abs() < 1e-12);
    }
    #[test]
    fn test_quality_metrics_with_error() {
        let mut m = WarmStartQualityMetrics::new();
        m.record_hit(8.0, 10.0);
        assert!((m.mean_relative_error() - 0.2).abs() < 1e-12);
        assert!((m.quality_score() - 0.8).abs() < 1e-12);
    }
    #[test]
    fn test_adaptive_scaler_increases_on_good_quality() {
        let mut scaler = AdaptiveWarmStartScaler::with_params(0.5, 0.3, 1.0, 0.1, 0.5);
        let mut metrics = WarmStartQualityMetrics::new();
        metrics.record_hit(10.0, 10.0);
        let old_scale = scaler.scale;
        scaler.update(&metrics);
        assert!(
            scaler.scale > old_scale,
            "Scale should increase on good quality"
        );
    }
    #[test]
    fn test_adaptive_scaler_decreases_on_poor_quality() {
        let mut scaler = AdaptiveWarmStartScaler::with_params(0.8, 0.3, 1.0, 0.1, 0.5);
        let metrics = WarmStartQualityMetrics::new();
        let old_scale = scaler.scale;
        scaler.update(&metrics);
        assert!(
            scaler.scale < old_scale,
            "Scale should decrease on poor quality"
        );
    }
    #[test]
    fn test_adaptive_scaler_stays_in_bounds() {
        let mut scaler = AdaptiveWarmStartScaler::with_params(0.3, 0.3, 1.0, 1.0, 0.5);
        let empty_metrics = WarmStartQualityMetrics::new();
        for _ in 0..100 {
            scaler.update(&empty_metrics);
        }
        assert!(scaler.scale >= scaler.min_scale - 1e-12);
        let mut good_metrics = WarmStartQualityMetrics::new();
        good_metrics.record_hit(10.0, 10.0);
        for _ in 0..100 {
            scaler.update(&good_metrics);
        }
        assert!(scaler.scale <= scaler.max_scale + 1e-12);
    }
    #[test]
    fn test_adaptive_scaler_apply() {
        let scaler = AdaptiveWarmStartScaler::with_params(0.5, 0.3, 1.0, 0.1, 0.5);
        let cache = WarmStartCache::with_impulses(10.0, 20.0, 30.0);
        let scaled = scaler.apply(&cache);
        assert!((scaled.lambda_n - 5.0).abs() < 1e-12);
        assert!((scaled.lambda_t1 - 10.0).abs() < 1e-12);
        assert!((scaled.lambda_t2 - 15.0).abs() < 1e-12);
    }
    #[test]
    fn test_pair_matcher_exact_match() {
        let mut matcher = ConstraintPairMatcher::new(0.01, 0.99);
        let fp = ContactFingerprint {
            body_pair: WarmStartMap::key(1, 2),
            local_point_a: [1.0, 0.0, 0.0],
            local_point_b: [0.0, 1.0, 0.0],
            normal: [0.0, 1.0, 0.0],
        };
        matcher.push_frame(std::slice::from_ref(&fp));
        let result = matcher.find_match(&fp);
        assert_eq!(result, Some(0));
    }
    #[test]
    fn test_pair_matcher_no_match_distance() {
        let mut matcher = ConstraintPairMatcher::new(0.01, 0.99);
        let fp_prev = ContactFingerprint {
            body_pair: WarmStartMap::key(1, 2),
            local_point_a: [1.0, 0.0, 0.0],
            local_point_b: [0.0, 1.0, 0.0],
            normal: [0.0, 1.0, 0.0],
        };
        matcher.push_frame(&[fp_prev]);
        let fp_curr = ContactFingerprint {
            body_pair: WarmStartMap::key(1, 2),
            local_point_a: [2.0, 0.0, 0.0],
            local_point_b: [0.0, 1.0, 0.0],
            normal: [0.0, 1.0, 0.0],
        };
        assert_eq!(matcher.find_match(&fp_curr), None);
    }
    #[test]
    fn test_pair_matcher_no_match_normal() {
        let mut matcher = ConstraintPairMatcher::new(0.01, 0.99);
        let fp_prev = ContactFingerprint {
            body_pair: WarmStartMap::key(1, 2),
            local_point_a: [1.0, 0.0, 0.0],
            local_point_b: [0.0, 1.0, 0.0],
            normal: [0.0, 1.0, 0.0],
        };
        matcher.push_frame(&[fp_prev]);
        let fp_curr = ContactFingerprint {
            body_pair: WarmStartMap::key(1, 2),
            local_point_a: [1.0, 0.0, 0.0],
            local_point_b: [0.0, 1.0, 0.0],
            normal: [0.0, -1.0, 0.0],
        };
        assert_eq!(matcher.find_match(&fp_curr), None);
    }
    #[test]
    fn test_pair_matcher_clear() {
        let mut matcher = ConstraintPairMatcher::new(0.01, 0.99);
        let fp = ContactFingerprint {
            body_pair: WarmStartMap::key(1, 2),
            local_point_a: [0.0; 3],
            local_point_b: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
        };
        matcher.push_frame(std::slice::from_ref(&fp));
        assert_eq!(matcher.prev_pair_count(), 1);
        matcher.clear();
        assert_eq!(matcher.prev_pair_count(), 0);
    }
    /// solve_normal_impulse: zero penetration + separating velocity → zero impulse.
    #[test]
    fn test_solve_normal_impulse_separating() {
        let warm = WarmStartCache::default();
        let rel_vel_n = 2.0;
        let (delta_lambda, updated) = ContactVelocitySolver::solve_normal_impulse(
            rel_vel_n,
            1.0,
            1.0,
            1.0,
            1.0,
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.0,
            1.0 / 60.0,
            &warm,
        );
        assert!(
            delta_lambda <= 0.0 || delta_lambda.abs() < 1e-12,
            "Separating contact should produce no (or zero) impulse, got {delta_lambda}"
        );
        assert!(
            updated.lambda_n >= 0.0,
            "Accumulated normal impulse must stay non-negative"
        );
    }
    /// Penetrating contact should produce a positive normal impulse.
    #[test]
    fn test_solve_normal_impulse_penetrating() {
        let warm = WarmStartCache::default();
        let (delta_lambda, updated) = ContactVelocitySolver::solve_normal_impulse(
            -1.0,
            1.0,
            1.0,
            0.0,
            0.0,
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.01,
            1.0 / 60.0,
            &warm,
        );
        assert!(
            delta_lambda > 0.0,
            "Approaching contact should produce positive impulse, got {delta_lambda}"
        );
        assert!(updated.lambda_n > 0.0);
    }
    /// Warm-started normal solve should accumulate correctly.
    #[test]
    fn test_solve_normal_impulse_accumulation() {
        let warm = WarmStartCache::with_impulses(5.0, 0.0, 0.0);
        let (delta, updated) = ContactVelocitySolver::solve_normal_impulse(
            -1.0,
            1.0,
            1.0,
            0.0,
            0.0,
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.01,
            1.0 / 60.0,
            &warm,
        );
        assert!(
            updated.lambda_n >= 5.0 - 1e-12,
            "Accumulated impulse should grow from warm start"
        );
        assert!(
            (updated.lambda_n - warm.lambda_n - delta).abs() < 1e-12,
            "Delta should equal difference in accumulated impulse"
        );
    }
    /// Friction cone clamping: μ=0.5, λ_n=10 → |λ_t| ≤ 5.
    #[test]
    fn test_friction_cone_clamping() {
        let mu = 0.5_f64;
        let lambda_n = 10.0_f64;
        let friction_limit = mu * lambda_n;
        let (lt1, lt2) =
            ContactVelocitySolver::solve_friction_impulse(100.0, 100.0, friction_limit, 1.0, 1.0);
        let mag = (lt1 * lt1 + lt2 * lt2).sqrt();
        assert!(
            mag <= friction_limit + 1e-10,
            "Friction magnitude {mag} must be <= μ*λ_n = {friction_limit}"
        );
    }
    /// Friction within cone should not be clamped.
    #[test]
    fn test_friction_within_cone() {
        let (lt1, lt2) = ContactVelocitySolver::solve_friction_impulse(0.1, 0.0, 100.0, 1.0, 1.0);
        assert!((lt1 + 0.05).abs() < 1e-12);
        assert!(lt2.abs() < 1e-12);
    }
    /// apply_impulse: linear velocity change = impulse / mass.
    #[test]
    fn test_apply_impulse_linear() {
        let mut vel = [0.0_f64; 3];
        let mut ang_vel = [0.0_f64; 3];
        let r = [0.0_f64; 3];
        let impulse = [3.0, 0.0, 0.0];
        let inv_mass = 0.5;
        let inv_inertia = 1.0;
        apply_impulse(&mut vel, &mut ang_vel, r, impulse, inv_mass, inv_inertia);
        assert!(
            (vel[0] - 1.5).abs() < 1e-12,
            "Expected vx = 1.5, got {}",
            vel[0]
        );
        assert!(vel[1].abs() < 1e-12, "vy should be zero");
        assert!(vel[2].abs() < 1e-12, "vz should be zero");
        assert!(
            ang_vel[0].abs() < 1e-12 && ang_vel[1].abs() < 1e-12 && ang_vel[2].abs() < 1e-12,
            "No angular change expected"
        );
    }
    /// apply_impulse: angular velocity change from off-center impulse.
    #[test]
    fn test_apply_impulse_angular() {
        let mut vel = [0.0_f64; 3];
        let mut ang_vel = [0.0_f64; 3];
        let r = [1.0, 0.0, 0.0];
        let impulse = [0.0, 1.0, 0.0];
        let inv_mass = 1.0;
        let inv_inertia = 1.0;
        apply_impulse(&mut vel, &mut ang_vel, r, impulse, inv_mass, inv_inertia);
        assert!(
            (ang_vel[2] - 1.0).abs() < 1e-12,
            "Expected ang_vel_z = 1.0, got {}",
            ang_vel[2]
        );
    }
    #[test]
    fn test_dist_sq_3() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 6.0, 3.0];
        assert!((dist_sq_3(a, b) - 25.0).abs() < 1e-12);
    }
    #[test]
    fn test_dist_sq_3_same_point() {
        let a = [1.0, 2.0, 3.0];
        assert!(dist_sq_3(a, a).abs() < 1e-12);
    }
}
#[cfg(test)]
mod tests_warm_start_extended {

    use crate::warm_start::AdaptiveWarmStart;

    use crate::warm_start::WarmStartImpulseCache;

    use crate::warm_start::WarmStartQuality;
    use crate::warm_start::WarmStartRecord;
    #[test]
    fn test_insert_and_lookup() {
        let mut cache = WarmStartImpulseCache::new(5);
        cache.insert((1, 2), [1.0, 2.0, 3.0]);
        let result = cache.lookup((1, 2));
        assert!(result.is_some());
        let imp = result.unwrap();
        assert!((imp[0] - 1.0).abs() < 1e-12);
        assert!((imp[1] - 2.0).abs() < 1e-12);
        assert!((imp[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_lookup_missing_returns_none() {
        let cache = WarmStartImpulseCache::new(5);
        assert!(cache.lookup((9, 9)).is_none());
    }
    #[test]
    fn test_age_all_increments() {
        let mut cache = WarmStartImpulseCache::new(10);
        cache.insert((0, 1), [1.0, 0.0, 0.0]);
        assert_eq!(cache.records[&(0, 1)].age, 0);
        cache.age_all();
        assert_eq!(cache.records[&(0, 1)].age, 1);
        cache.age_all();
        assert_eq!(cache.records[&(0, 1)].age, 2);
    }
    #[test]
    fn test_evict_removes_old() {
        let mut cache = WarmStartImpulseCache::new(2);
        cache.insert((1, 2), [0.0, 0.0, 0.0]);
        cache.age_all();
        cache.age_all();
        cache.age_all();
        cache.evict_aged();
        assert!(
            cache.lookup((1, 2)).is_none(),
            "old record should be evicted"
        );
    }
    #[test]
    fn test_evict_keeps_young() {
        let mut cache = WarmStartImpulseCache::new(5);
        cache.insert((1, 2), [0.0, 0.0, 0.0]);
        cache.age_all();
        cache.evict_aged();
        assert!(cache.lookup((1, 2)).is_some(), "young record should remain");
    }
    #[test]
    fn test_decay_reduces_impulse() {
        let mut record = WarmStartRecord::new((0, 1), [4.0, 2.0, 1.0]);
        record.decay(0.5);
        assert!((record.impulse[0] - 2.0).abs() < 1e-12);
        assert!((record.impulse[1] - 1.0).abs() < 1e-12);
        assert!((record.impulse[2] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_quality_fresh() {
        let mut cache = WarmStartImpulseCache::new(5);
        cache.insert((3, 4), [0.0, 0.0, 0.0]);
        assert_eq!(cache.quality((3, 4)), WarmStartQuality::Fresh);
    }
    #[test]
    fn test_quality_aged() {
        let mut cache = WarmStartImpulseCache::new(10);
        cache.insert((3, 4), [0.0, 0.0, 0.0]);
        cache.age_all();
        match cache.quality((3, 4)) {
            WarmStartQuality::Aged(score) => {
                assert!(score > 0.0 && score < 1.0, "score = {score}");
            }
            other => panic!("expected Aged, got {:?}", other),
        }
    }
    #[test]
    fn test_quality_stale_missing() {
        let cache = WarmStartImpulseCache::new(5);
        assert_eq!(cache.quality((99, 99)), WarmStartQuality::Stale);
    }
    #[test]
    fn test_quality_stale_over_max_age() {
        let mut cache = WarmStartImpulseCache::new(2);
        cache.insert((1, 2), [1.0, 0.0, 0.0]);
        for _ in 0..3 {
            cache.age_all();
        }
        assert_eq!(cache.quality((1, 2)), WarmStartQuality::Stale);
    }
    #[test]
    fn test_adaptive_warm_start_uses_cache_when_fresh() {
        let mut aws = AdaptiveWarmStart::new(5, 0.5);
        aws.cache.insert((1, 2), [9.0, 8.0, 7.0]);
        let result = aws.apply_if_valid((1, 2), [0.0, 0.0, 0.0]);
        assert!((result[0] - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_adaptive_warm_start_falls_back_when_stale() {
        let aws = AdaptiveWarmStart::new(5, 0.5);
        let fallback = [1.0, 2.0, 3.0];
        let result = aws.apply_if_valid((99, 100), fallback);
        assert_eq!(result, fallback);
    }
}
#[cfg(test)]
mod tests_warm_start_new {
    use super::*;
    use crate::warm_start::AdaptiveVelocityGate;

    use crate::warm_start::AlphaImpulseAging;
    use crate::warm_start::BodyVelocitySnapshot;
    use crate::warm_start::ContactPoint3D;
    use crate::warm_start::IslandId;
    use crate::warm_start::IslandWarmStartData;
    use crate::warm_start::IslandWarmStartManager;
    use crate::warm_start::ManifoldQuality;
    use crate::warm_start::PersistentManifoldMatcher;
    use crate::warm_start::QualityBasedAging;
    use crate::warm_start::WarmStartContribStats;

    use crate::warm_start::WarmStartMap;

    #[test]
    fn test_alpha_aging_zero_alpha_preserves_impulse() {
        let aging = AlphaImpulseAging::new(0.0);
        let imp = [3.0, 2.0, 1.0];
        let aged = aging.age(imp);
        assert!((aged[0] - 3.0).abs() < 1e-12);
        assert!((aged[1] - 2.0).abs() < 1e-12);
        assert!((aged[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_alpha_aging_full_alpha_zeroes_impulse() {
        let aging = AlphaImpulseAging::new(1.0);
        let aged = aging.age([5.0, 5.0, 5.0]);
        assert!(aged[0].abs() < 1e-12);
        assert!(aged[1].abs() < 1e-12);
        assert!(aged[2].abs() < 1e-12);
    }
    #[test]
    fn test_alpha_aging_half_alpha_halves_impulse() {
        let aging = AlphaImpulseAging::new(0.5);
        let aged = aging.age([10.0, 4.0, 2.0]);
        assert!((aged[0] - 5.0).abs() < 1e-12);
        assert!((aged[1] - 2.0).abs() < 1e-12);
        assert!((aged[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_alpha_aging_clamps_above_one() {
        let aging = AlphaImpulseAging::new(2.0);
        let aged = aging.age([7.0, 0.0, 0.0]);
        assert!(aged[0].abs() < 1e-12);
    }
    #[test]
    fn test_alpha_aging_map() {
        let mut map = WarmStartMap::default();
        map.cache.insert(
            (1, 2),
            WarmStartCache {
                lambda_n: 10.0,
                lambda_t1: 4.0,
                lambda_t2: 2.0,
            },
        );
        let aging = AlphaImpulseAging::new(0.5);
        aging.age_map(&mut map);
        let c = &map.cache[&(1, 2)];
        assert!((c.lambda_n - 5.0).abs() < 1e-12);
        assert!((c.lambda_t1 - 2.0).abs() < 1e-12);
        assert!((c.lambda_t2 - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_manifold_quality_scale_full() {
        let q = ManifoldQuality::new(1.0);
        let scaled = q.scale_impulse([4.0, 2.0, 1.0]);
        assert!((scaled[0] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_manifold_quality_scale_half() {
        let q = ManifoldQuality::new(0.5);
        let scaled = q.scale_impulse([8.0, 0.0, 0.0]);
        assert!((scaled[0] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_manifold_quality_clamp_negative() {
        let q = ManifoldQuality::new(-1.0);
        let scaled = q.scale_impulse([9.0, 9.0, 9.0]);
        assert!(scaled[0].abs() < 1e-12);
    }
    #[test]
    fn test_quality_based_aging_above_floor() {
        let policy = QualityBasedAging::new(0.2, 0.3);
        let scale = policy.effective_scale(0.8);
        assert!((scale - 0.64).abs() < 1e-12, "scale={scale}");
    }
    #[test]
    fn test_quality_based_aging_below_floor_zeroes() {
        let policy = QualityBasedAging::new(0.1, 0.5);
        let scale = policy.effective_scale(0.2);
        assert!(scale.abs() < 1e-12, "should be zero below floor");
    }
    #[test]
    fn test_quality_based_aging_age_impulse() {
        let policy = QualityBasedAging::new(0.0, 0.0);
        let aged = policy.age_impulse([6.0, 0.0, 0.0], 1.0);
        assert!((aged[0] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_velocity_gate_allows_small_change() {
        let gate = AdaptiveVelocityGate::new(1.0);
        let snap = BodyVelocitySnapshot::new([0.0, 0.0, 0.0], [0.0; 3]);
        assert!(gate.should_apply(&snap, [0.5, 0.0, 0.0]));
    }
    #[test]
    fn test_velocity_gate_blocks_large_change() {
        let gate = AdaptiveVelocityGate::new(0.1);
        let snap = BodyVelocitySnapshot::new([0.0, 0.0, 0.0], [0.0; 3]);
        assert!(!gate.should_apply(&snap, [5.0, 0.0, 0.0]));
    }
    #[test]
    fn test_velocity_gate_apply_or_zero_returns_impulse_when_valid() {
        let gate = AdaptiveVelocityGate::new(2.0);
        let snap = BodyVelocitySnapshot::new([1.0, 0.0, 0.0], [0.0; 3]);
        let imp = [3.0, 0.0, 0.0];
        let result = gate.apply_or_zero(&snap, [1.5, 0.0, 0.0], imp);
        assert!((result[0] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_velocity_gate_apply_or_zero_returns_zero_when_blocked() {
        let gate = AdaptiveVelocityGate::new(0.01);
        let snap = BodyVelocitySnapshot::new([0.0, 0.0, 0.0], [0.0; 3]);
        let imp = [3.0, 0.0, 0.0];
        let result = gate.apply_or_zero(&snap, [5.0, 0.0, 0.0], imp);
        assert!(result[0].abs() < 1e-12);
    }
    #[test]
    fn test_body_velocity_snapshot_linear_delta_sq() {
        let snap = BodyVelocitySnapshot::new([1.0, 2.0, 3.0], [0.0; 3]);
        let delta_sq = snap.linear_delta_sq([4.0, 2.0, 3.0]);
        assert!((delta_sq - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_manifold_matcher_exact_match() {
        let matcher = PersistentManifoldMatcher::new(0.1);
        let pt = ContactPoint3D::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0);
        let prev = vec![pt];
        let curr = vec![pt];
        let matches = matcher.match_points(&curr, &prev);
        assert_eq!(matches[0], Some(0));
    }
    #[test]
    fn test_manifold_matcher_no_match_far_apart() {
        let matcher = PersistentManifoldMatcher::new(0.1);
        let prev = vec![ContactPoint3D::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0)];
        let curr = vec![ContactPoint3D::new([5.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0)];
        let matches = matcher.match_points(&curr, &prev);
        assert_eq!(matches[0], None);
    }
    #[test]
    fn test_manifold_matcher_picks_closest() {
        let matcher = PersistentManifoldMatcher::new(2.0);
        let prev = vec![
            ContactPoint3D::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0),
            ContactPoint3D::new([1.5, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0),
        ];
        let curr = vec![ContactPoint3D::new([1.4, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0)];
        let matches = matcher.match_points(&curr, &prev);
        assert_eq!(matches[0], Some(1));
    }
    #[test]
    fn test_manifold_matcher_match_count() {
        let matcher = PersistentManifoldMatcher::new(0.5);
        let prev = vec![
            ContactPoint3D::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0),
            ContactPoint3D::new([10.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0),
        ];
        let curr = vec![
            ContactPoint3D::new([0.1, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0),
            ContactPoint3D::new([99.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0),
        ];
        assert_eq!(matcher.match_count(&curr, &prev), 1);
    }
    #[test]
    fn test_contrib_stats_initial_hit_rate_zero() {
        let stats = WarmStartContribStats::new();
        assert!((stats.hit_rate()).abs() < 1e-12);
    }
    #[test]
    fn test_contrib_stats_hit_rate() {
        let mut stats = WarmStartContribStats::new();
        stats.record_hit([1.0, 0.0, 0.0]);
        stats.record_hit([2.0, 0.0, 0.0]);
        stats.record_miss();
        assert!(
            (stats.hit_rate() - 2.0 / 3.0).abs() < 1e-10,
            "rate={}",
            stats.hit_rate()
        );
    }
    #[test]
    fn test_contrib_stats_peak_impulse() {
        let mut stats = WarmStartContribStats::new();
        stats.record_hit([3.0, 4.0, 0.0]);
        stats.record_hit([1.0, 0.0, 0.0]);
        assert!(
            (stats.peak_impulse - 5.0).abs() < 1e-10,
            "peak={}",
            stats.peak_impulse
        );
    }
    #[test]
    fn test_contrib_stats_rms_impulse() {
        let mut stats = WarmStartContribStats::new();
        stats.record_hit([3.0, 0.0, 0.0]);
        stats.record_hit([4.0, 0.0, 0.0]);
        let expected_rms = 12.5_f64.sqrt();
        assert!((stats.rms_impulse() - expected_rms).abs() < 1e-10);
    }
    #[test]
    fn test_contrib_stats_reset() {
        let mut stats = WarmStartContribStats::new();
        stats.record_hit([1.0, 1.0, 1.0]);
        stats.record_miss();
        stats.reset();
        assert_eq!(stats.hits, 0);
        assert_eq!(stats.misses, 0);
        assert!((stats.total_impulse_sq).abs() < 1e-12);
    }
    #[test]
    fn test_island_data_insert_and_get() {
        let mut data = IslandWarmStartData::new(IslandId(1));
        data.insert_pair((3, 5), [1.0, 2.0, 3.0]);
        let imp = data.get_pair((3, 5)).unwrap();
        assert!((imp[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_island_data_symmetric_key() {
        let mut data = IslandWarmStartData::new(IslandId(0));
        data.insert_pair((7, 2), [5.0, 0.0, 0.0]);
        assert!(data.get_pair((7, 2)).is_some());
        assert!(data.get_pair((2, 7)).is_some());
    }
    #[test]
    fn test_island_data_decay() {
        let mut data = IslandWarmStartData::new(IslandId(0));
        data.insert_pair((1, 2), [8.0, 0.0, 0.0]);
        data.decay_all(0.5);
        let imp = data.get_pair((1, 2)).unwrap();
        assert!((imp[0] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_island_data_age_increment() {
        let mut data = IslandWarmStartData::new(IslandId(42));
        assert_eq!(data.age, 0);
        data.increment_age();
        data.increment_age();
        assert_eq!(data.age, 2);
    }
    #[test]
    fn test_island_manager_store_and_retrieve() {
        let mut mgr = IslandWarmStartManager::new(5);
        let mut d = IslandWarmStartData::new(IslandId(10));
        d.insert_pair((0, 1), [1.0, 0.0, 0.0]);
        mgr.store_island(d);
        assert!(mgr.get_island(IslandId(10)).is_some());
        assert!(mgr.get_island(IslandId(99)).is_none());
    }
    #[test]
    fn test_island_manager_advance_frame_evicts_old() {
        let mut mgr = IslandWarmStartManager::new(2);
        mgr.store_island(IslandWarmStartData::new(IslandId(1)));
        mgr.advance_frame();
        mgr.advance_frame();
        mgr.advance_frame();
        assert!(
            mgr.get_island(IslandId(1)).is_none(),
            "old island should be evicted"
        );
    }
    #[test]
    fn test_island_manager_advance_frame_keeps_young() {
        let mut mgr = IslandWarmStartManager::new(5);
        mgr.store_island(IslandWarmStartData::new(IslandId(2)));
        mgr.advance_frame();
        assert!(mgr.get_island(IslandId(2)).is_some());
    }
    #[test]
    fn test_island_manager_decay_all() {
        let mut mgr = IslandWarmStartManager::new(10);
        let mut d = IslandWarmStartData::new(IslandId(3));
        d.insert_pair((0, 1), [10.0, 0.0, 0.0]);
        mgr.store_island(d);
        mgr.decay_all(0.5);
        let imp = mgr
            .get_island(IslandId(3))
            .unwrap()
            .get_pair((0, 1))
            .unwrap();
        assert!((imp[0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_island_manager_island_count() {
        let mut mgr = IslandWarmStartManager::new(10);
        assert_eq!(mgr.island_count(), 0);
        mgr.store_island(IslandWarmStartData::new(IslandId(1)));
        mgr.store_island(IslandWarmStartData::new(IslandId(2)));
        assert_eq!(mgr.island_count(), 2);
    }
}
/// Select the best warm-start candidate from a slice according to a policy.
///
/// Returns `None` if the slice is empty.
#[allow(dead_code)]
pub fn select_best_candidate(
    candidates: &[WarmStartCandidate],
    policy: CacheSelectionPolicy,
) -> Option<&WarmStartCandidate> {
    if candidates.is_empty() {
        return None;
    }
    match policy {
        CacheSelectionPolicy::MostRecent => candidates.iter().min_by_key(|c| c.age),
        CacheSelectionPolicy::LargestMagnitude => candidates.iter().max_by(|a, b| {
            a.cache
                .magnitude()
                .partial_cmp(&b.cache.magnitude())
                .expect("operation should succeed")
        }),
        CacheSelectionPolicy::SmallestMagnitude => candidates.iter().min_by(|a, b| {
            a.cache
                .magnitude()
                .partial_cmp(&b.cache.magnitude())
                .expect("operation should succeed")
        }),
        CacheSelectionPolicy::ClosestQuality(target) => candidates.iter().min_by(|a, b| {
            let da = (a.quality - target).abs();
            let db = (b.quality - target).abs();
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        }),
    }
}
/// Blend two `WarmStartCache` entries according to `mode`.
#[allow(dead_code)]
pub fn blend_caches(
    prev: &WarmStartCache,
    current: &WarmStartCache,
    mode: BlendMode,
) -> WarmStartCache {
    match mode {
        BlendMode::Lerp(t) => prev.lerp(current, t),
        BlendMode::Average => prev.lerp(current, 0.5),
        BlendMode::Max => WarmStartCache {
            lambda_n: if prev.lambda_n.abs() >= current.lambda_n.abs() {
                prev.lambda_n
            } else {
                current.lambda_n
            },
            lambda_t1: if prev.lambda_t1.abs() >= current.lambda_t1.abs() {
                prev.lambda_t1
            } else {
                current.lambda_t1
            },
            lambda_t2: if prev.lambda_t2.abs() >= current.lambda_t2.abs() {
                prev.lambda_t2
            } else {
                current.lambda_t2
            },
        },
        BlendMode::Min => WarmStartCache {
            lambda_n: if prev.lambda_n.abs() <= current.lambda_n.abs() {
                prev.lambda_n
            } else {
                current.lambda_n
            },
            lambda_t1: if prev.lambda_t1.abs() <= current.lambda_t1.abs() {
                prev.lambda_t1
            } else {
                current.lambda_t1
            },
            lambda_t2: if prev.lambda_t2.abs() <= current.lambda_t2.abs() {
                prev.lambda_t2
            } else {
                current.lambda_t2
            },
        },
    }
}
#[cfg(test)]
mod tests_warm_start_expanded {
    use super::*;
    use crate::warm_start::AdaptiveImpulseThreshold;

    use crate::warm_start::CacheStatistics;

    use crate::warm_start::ExponentialImpulseAging;
    use crate::warm_start::ImpulseHistory;

    use crate::warm_start::WarmStartMap;

    #[test]
    fn exp_aging_zero_alpha_no_decay() {
        let aging = ExponentialImpulseAging::new(0.0, 0.0);
        assert!((aging.age_scalar(5.0) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn exp_aging_full_alpha_zeroes() {
        let aging = ExponentialImpulseAging::new(1.0, 0.0);
        assert!(aging.age_scalar(5.0).abs() < 1e-12);
    }
    #[test]
    fn exp_aging_half_alpha_halves() {
        let aging = ExponentialImpulseAging::new(0.5, 0.0);
        assert!((aging.age_scalar(4.0) - 2.0).abs() < 1e-12);
    }
    #[test]
    fn exp_aging_cache_returns_none_below_threshold() {
        let aging = ExponentialImpulseAging::new(0.9, 0.5);
        let cache = WarmStartCache::with_impulses(0.1, 0.0, 0.0);
        assert!(aging.age_cache(&cache).is_none());
    }
    #[test]
    fn exp_aging_cache_returns_some_above_threshold() {
        let aging = ExponentialImpulseAging::new(0.1, 0.001);
        let cache = WarmStartCache::with_impulses(10.0, 0.0, 0.0);
        assert!(aging.age_cache(&cache).is_some());
    }
    #[test]
    fn exp_aging_apply_to_map_prunes_small() {
        let aging = ExponentialImpulseAging::new(0.5, 1.0);
        let mut map = WarmStartMap::new();
        map.store(0, 1, WarmStartCache::with_impulses(0.5, 0.0, 0.0));
        aging.apply_to_map(&mut map);
        assert_eq!(map.len(), 0, "small entry should be pruned");
    }
    #[test]
    fn exp_aging_apply_to_map_keeps_large() {
        let aging = ExponentialImpulseAging::new(0.1, 0.001);
        let mut map = WarmStartMap::new();
        map.store(0, 1, WarmStartCache::with_impulses(100.0, 0.0, 0.0));
        aging.apply_to_map(&mut map);
        assert_eq!(map.len(), 1, "large entry should remain");
    }
    #[test]
    fn adaptive_threshold_starts_at_floor() {
        let t = AdaptiveImpulseThreshold::new(0.01, 10.0);
        assert!((t.threshold - 0.01).abs() < 1e-12);
    }
    #[test]
    fn adaptive_threshold_observe_increases() {
        let mut t = AdaptiveImpulseThreshold::new(0.01, 10.0);
        t.observe(5.0);
        assert!(
            t.threshold > 0.01,
            "threshold should increase: {}",
            t.threshold
        );
    }
    #[test]
    fn adaptive_threshold_tick_decays() {
        let mut t = AdaptiveImpulseThreshold::new(0.01, 10.0);
        t.observe(5.0);
        let after_observe = t.threshold;
        t.tick();
        assert!(t.threshold < after_observe, "tick should decay threshold");
    }
    #[test]
    fn adaptive_threshold_stays_above_floor() {
        let mut t = AdaptiveImpulseThreshold::new(0.01, 10.0);
        for _ in 0..1000 {
            t.tick();
        }
        assert!(t.threshold >= 0.01 - 1e-12, "should not go below floor");
    }
    #[test]
    fn adaptive_threshold_stays_below_ceiling() {
        let mut t = AdaptiveImpulseThreshold::new(0.01, 10.0);
        for _ in 0..1000 {
            t.observe(100.0);
        }
        assert!(t.threshold <= 10.0 + 1e-12, "should not exceed ceiling");
    }
    #[test]
    fn adaptive_threshold_exceeds() {
        let t = AdaptiveImpulseThreshold::new(1.0, 10.0);
        assert!(t.exceeds(2.0));
        assert!(!t.exceeds(0.5));
    }
    #[test]
    fn adaptive_threshold_reset() {
        let mut t = AdaptiveImpulseThreshold::new(0.01, 10.0);
        t.observe(5.0);
        t.reset();
        assert!((t.threshold - 0.01).abs() < 1e-12);
    }
    #[test]
    fn select_most_recent_picks_smallest_age() {
        let candidates = vec![
            WarmStartCandidate::new(WarmStartCache::with_impulses(1.0, 0.0, 0.0), 5, 0.5),
            WarmStartCandidate::new(WarmStartCache::with_impulses(2.0, 0.0, 0.0), 1, 0.8),
        ];
        let best = select_best_candidate(&candidates, CacheSelectionPolicy::MostRecent).unwrap();
        assert_eq!(best.age, 1);
    }
    #[test]
    fn select_largest_magnitude() {
        let candidates = vec![
            WarmStartCandidate::new(WarmStartCache::with_impulses(3.0, 0.0, 0.0), 0, 1.0),
            WarmStartCandidate::new(WarmStartCache::with_impulses(10.0, 0.0, 0.0), 0, 0.5),
            WarmStartCandidate::new(WarmStartCache::with_impulses(5.0, 0.0, 0.0), 0, 0.8),
        ];
        let best =
            select_best_candidate(&candidates, CacheSelectionPolicy::LargestMagnitude).unwrap();
        assert!((best.cache.lambda_n - 10.0).abs() < 1e-12);
    }
    #[test]
    fn select_smallest_magnitude() {
        let candidates = vec![
            WarmStartCandidate::new(WarmStartCache::with_impulses(3.0, 0.0, 0.0), 0, 1.0),
            WarmStartCandidate::new(WarmStartCache::with_impulses(10.0, 0.0, 0.0), 0, 0.5),
        ];
        let best =
            select_best_candidate(&candidates, CacheSelectionPolicy::SmallestMagnitude).unwrap();
        assert!((best.cache.lambda_n - 3.0).abs() < 1e-12);
    }
    #[test]
    fn select_closest_quality() {
        let candidates = vec![
            WarmStartCandidate::new(WarmStartCache::default(), 0, 0.2),
            WarmStartCandidate::new(WarmStartCache::default(), 0, 0.7),
            WarmStartCandidate::new(WarmStartCache::default(), 0, 0.9),
        ];
        let best =
            select_best_candidate(&candidates, CacheSelectionPolicy::ClosestQuality(0.75)).unwrap();
        assert!((best.quality - 0.7).abs() < 1e-10 || (best.quality - 0.9).abs() < 1e-10);
    }
    #[test]
    fn select_returns_none_for_empty() {
        let candidates: Vec<WarmStartCandidate> = vec![];
        assert!(select_best_candidate(&candidates, CacheSelectionPolicy::MostRecent).is_none());
    }
    #[test]
    fn candidate_composite_score_decreases_with_age() {
        let c0 = WarmStartCandidate::new(WarmStartCache::default(), 0, 1.0);
        let c5 = WarmStartCandidate::new(WarmStartCache::default(), 5, 1.0);
        assert!(c0.composite_score() > c5.composite_score());
    }
    #[test]
    fn cache_stats_initial_all_zero() {
        let s = CacheStatistics::new();
        assert_eq!(s.total_lookups, 0);
        assert_eq!(s.hits, 0);
        assert_eq!(s.misses, 0);
        assert_eq!(s.hit_rate(), 0.0);
    }
    #[test]
    fn cache_stats_hit_rate() {
        let mut s = CacheStatistics::new();
        s.record_hit(1.0);
        s.record_hit(2.0);
        s.record_miss();
        assert!((s.hit_rate() - 2.0 / 3.0).abs() < 1e-10);
    }
    #[test]
    fn cache_stats_mean_impulse() {
        let mut s = CacheStatistics::new();
        s.record_hit(2.0);
        s.record_hit(4.0);
        assert!((s.mean_impulse() - 3.0).abs() < 1e-10);
    }
    #[test]
    fn cache_stats_rms_impulse() {
        let mut s = CacheStatistics::new();
        s.record_hit(3.0);
        s.record_hit(4.0);
        let expected = (12.5_f64).sqrt();
        assert!(
            (s.rms_impulse() - expected).abs() < 1e-10,
            "rms: {}",
            s.rms_impulse()
        );
    }
    #[test]
    fn cache_stats_peak_impulse() {
        let mut s = CacheStatistics::new();
        s.record_hit(5.0);
        s.record_hit(3.0);
        s.record_hit(8.0);
        assert!((s.peak_impulse - 8.0).abs() < 1e-12);
    }
    #[test]
    fn cache_stats_reset() {
        let mut s = CacheStatistics::new();
        s.record_hit(5.0);
        s.record_miss();
        s.reset();
        assert_eq!(s.total_lookups, 0);
        assert_eq!(s.hits, 0);
        assert_eq!(s.peak_impulse, 0.0);
    }
    #[test]
    fn cache_stats_prune_count() {
        let mut s = CacheStatistics::new();
        s.record_prune();
        s.record_prune();
        assert_eq!(s.prunes, 2);
    }
    #[test]
    fn impulse_history_push_and_mean() {
        let mut h = ImpulseHistory::new(4);
        h.push(1.0);
        h.push(3.0);
        assert!((h.mean() - 2.0).abs() < 1e-12, "mean: {}", h.mean());
    }
    #[test]
    fn impulse_history_circular_overwrite() {
        let mut h = ImpulseHistory::new(3);
        h.push(1.0);
        h.push(2.0);
        h.push(3.0);
        h.push(10.0);
        assert_eq!(h.len(), 3);
    }
    #[test]
    fn impulse_history_max_val() {
        let mut h = ImpulseHistory::new(5);
        h.push(1.0);
        h.push(7.0);
        h.push(3.0);
        assert!((h.max_val() - 7.0).abs() < 1e-12);
    }
    #[test]
    fn impulse_history_all_below() {
        let mut h = ImpulseHistory::new(4);
        h.push(0.001);
        h.push(0.002);
        assert!(h.all_below(0.01));
        assert!(!h.all_below(0.001));
    }
    #[test]
    fn impulse_history_empty() {
        let h = ImpulseHistory::new(4);
        assert!(h.is_empty());
        assert_eq!(h.len(), 0);
        assert_eq!(h.mean(), 0.0);
        assert_eq!(h.max_val(), 0.0);
    }
    #[test]
    fn impulse_history_clear() {
        let mut h = ImpulseHistory::new(4);
        h.push(5.0);
        h.push(10.0);
        h.clear();
        assert!(h.is_empty());
        assert_eq!(h.mean(), 0.0);
    }
    #[test]
    fn impulse_history_capacity() {
        let h = ImpulseHistory::new(8);
        assert_eq!(h.capacity(), 8);
    }
    #[test]
    fn blend_caches_lerp_midpoint() {
        let a = WarmStartCache::with_impulses(0.0, 0.0, 0.0);
        let b = WarmStartCache::with_impulses(10.0, 4.0, 2.0);
        let blended = blend_caches(&a, &b, BlendMode::Lerp(0.5));
        assert!((blended.lambda_n - 5.0).abs() < 1e-12);
        assert!((blended.lambda_t1 - 2.0).abs() < 1e-12);
    }
    #[test]
    fn blend_caches_average() {
        let a = WarmStartCache::with_impulses(2.0, 0.0, 0.0);
        let b = WarmStartCache::with_impulses(8.0, 0.0, 0.0);
        let blended = blend_caches(&a, &b, BlendMode::Average);
        assert!((blended.lambda_n - 5.0).abs() < 1e-12);
    }
    #[test]
    fn blend_caches_max() {
        let a = WarmStartCache::with_impulses(3.0, 0.0, 0.0);
        let b = WarmStartCache::with_impulses(7.0, 0.0, 0.0);
        let blended = blend_caches(&a, &b, BlendMode::Max);
        assert!((blended.lambda_n - 7.0).abs() < 1e-12);
    }
    #[test]
    fn blend_caches_min() {
        let a = WarmStartCache::with_impulses(3.0, 0.0, 0.0);
        let b = WarmStartCache::with_impulses(7.0, 0.0, 0.0);
        let blended = blend_caches(&a, &b, BlendMode::Min);
        assert!((blended.lambda_n - 3.0).abs() < 1e-12);
    }
    #[test]
    fn blend_caches_lerp_zero_returns_prev() {
        let a = WarmStartCache::with_impulses(5.0, 1.0, 2.0);
        let b = WarmStartCache::with_impulses(0.0, 0.0, 0.0);
        let blended = blend_caches(&a, &b, BlendMode::Lerp(0.0));
        assert!((blended.lambda_n - 5.0).abs() < 1e-12);
    }
    #[test]
    fn blend_caches_lerp_one_returns_current() {
        let a = WarmStartCache::with_impulses(5.0, 0.0, 0.0);
        let b = WarmStartCache::with_impulses(9.0, 0.0, 0.0);
        let blended = blend_caches(&a, &b, BlendMode::Lerp(1.0));
        assert!((blended.lambda_n - 9.0).abs() < 1e-12);
    }
}
