//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod extra_tests {
    use super::super::*;
    use crate::ccd_constraints::CcdBroadphasePair;
    use crate::ccd_constraints::CcdStats;
    use crate::ccd_constraints::SpeculativeContactConstraint;
    use crate::ccd_constraints::SubStepCcd;
    use crate::ccd_constraints::ToiPair;
    use crate::ccd_constraints::ToiPositionCorrector;
    use crate::ccd_constraints::ToiQueue;
    use oxiphysics_core::math::Vec3;
    use oxiphysics_rigid::RigidBody;
    use oxiphysics_rigid::RigidBodySet;
    #[test]
    fn speculative_contact_target_velocity_penetrating() {
        let mut bodies = RigidBodySet::new();
        let h = bodies.insert(RigidBody::new(1.0));
        let sc = SpeculativeContactConstraint::new(
            h,
            Vec3::new(0.0, 1.0, 0.0),
            -0.01,
            Vec3::zeros(),
            0.0,
            0.3,
        );
        let v_target = sc.target_normal_velocity(1.0 / 60.0);
        assert!(
            v_target < 0.0,
            "target velocity should be negative (constraint active): {v_target}"
        );
    }
    #[test]
    fn speculative_contact_not_active_when_no_penetration() {
        let mut bodies = RigidBodySet::new();
        let h = bodies.insert(RigidBody::new(1.0));
        let sc = SpeculativeContactConstraint::new(
            h,
            Vec3::new(0.0, 1.0, 0.0),
            0.05,
            Vec3::zeros(),
            0.0,
            0.3,
        );
        assert!(!sc.is_active(), "should not be active when not penetrating");
    }
    #[test]
    fn speculative_contact_active_when_penetrating() {
        let mut bodies = RigidBodySet::new();
        let h = bodies.insert(RigidBody::new(1.0));
        let sc = SpeculativeContactConstraint::new(
            h,
            Vec3::new(0.0, 1.0, 0.0),
            -0.02,
            Vec3::zeros(),
            0.0,
            0.3,
        );
        assert!(sc.is_active());
    }
    #[test]
    fn speculative_contact_zero_dt_returns_zero_velocity() {
        let mut bodies = RigidBodySet::new();
        let h = bodies.insert(RigidBody::new(1.0));
        let sc = SpeculativeContactConstraint::new(
            h,
            Vec3::new(0.0, 1.0, 0.0),
            -0.01,
            Vec3::zeros(),
            0.0,
            0.0,
        );
        assert_eq!(sc.target_normal_velocity(0.0), 0.0);
    }
    #[test]
    fn toi_corrector_zero_penetration() {
        let c = ToiPositionCorrector::default();
        assert_eq!(c.correction_magnitude(0.0, 2.0, 1.0 / 60.0), 0.0);
    }
    #[test]
    fn toi_corrector_nonzero_penetration() {
        let c = ToiPositionCorrector::default();
        let mag = c.correction_magnitude(0.05, 2.0, 1.0 / 60.0);
        assert!(
            mag > 0.0,
            "deep penetration should produce positive correction: {mag}"
        );
    }
    #[test]
    fn toi_corrector_zero_inv_mass() {
        let c = ToiPositionCorrector::default();
        assert_eq!(c.correction_magnitude(0.1, 0.0, 0.016), 0.0);
    }
    #[test]
    fn toi_corrector_zero_dt() {
        let c = ToiPositionCorrector::default();
        assert_eq!(c.correction_magnitude(0.1, 2.0, 0.0), 0.0);
    }
    #[test]
    fn toi_corrector_custom_baumgarte() {
        let c = ToiPositionCorrector::new(0.5, 0.0, 1.0);
        let mag = c.correction_magnitude(0.1, 1.0, 1.0);
        assert!((mag - 0.05).abs() < 1e-10, "mag={mag}");
    }
    #[test]
    fn toi_queue_empty_initially() {
        let q = ToiQueue::new();
        assert_eq!(q.total_count(), 0);
        assert_eq!(q.pending_count(), 0);
        assert!(q.is_done());
    }
    #[test]
    fn toi_queue_push_and_pop() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let mut q = ToiQueue::new();
        q.push(ToiPair::new(ha, hb, 0.3));
        assert_eq!(q.pending_count(), 1);
        let pair = q.pop_earliest().unwrap();
        assert!((pair.toi - 0.3).abs() < 1e-12);
        assert_eq!(q.pending_count(), 0);
    }
    #[test]
    fn toi_queue_pop_earliest_returns_smallest_toi() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let hc = bodies.insert(RigidBody::new(1.0));
        let mut q = ToiQueue::new();
        q.push(ToiPair::new(ha, hb, 0.8));
        q.push(ToiPair::new(hb, hc, 0.2));
        q.push(ToiPair::new(ha, hc, 0.5));
        let first = q.pop_earliest().unwrap();
        assert!((first.toi - 0.2).abs() < 1e-12, "first={}", first.toi);
    }
    #[test]
    fn toi_queue_sorted_on_insert() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let mut q = ToiQueue::new();
        q.push(ToiPair::new(ha, hb, 0.9));
        q.push(ToiPair::new(ha, hb, 0.1));
        q.push(ToiPair::new(ha, hb, 0.5));
        let p1 = q.pop_earliest().unwrap();
        let p2 = q.pop_earliest().unwrap();
        let p3 = q.pop_earliest().unwrap();
        assert!(
            p1.toi <= p2.toi && p2.toi <= p3.toi,
            "not sorted: {},{},{}",
            p1.toi,
            p2.toi,
            p3.toi
        );
    }
    #[test]
    fn toi_queue_clear() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let mut q = ToiQueue::new();
        q.push(ToiPair::new(ha, hb, 0.5));
        q.clear();
        assert!(q.is_done());
        assert_eq!(q.total_count(), 0);
    }
    #[test]
    fn substep_ccd_no_events_single_step() {
        let ccd = SubStepCcd::default();
        let steps = ccd.build_substeps(&[]);
        assert_eq!(steps.len(), 1);
        assert!((steps[0].0).abs() < 1e-12);
        assert!((steps[0].1 - 1.0).abs() < 1e-12);
    }
    #[test]
    fn substep_ccd_one_event_two_steps() {
        let ccd = SubStepCcd::default();
        let steps = ccd.build_substeps(&[0.3]);
        assert_eq!(steps.len(), 2, "steps: {steps:?}");
        assert!((steps[0].0).abs() < 1e-12);
        assert!((steps[0].1 - 0.3).abs() < 1e-12);
        assert!((steps[1].0 - 0.3).abs() < 1e-12);
        assert!((steps[1].1 - 1.0).abs() < 1e-12);
    }
    #[test]
    fn substep_ccd_multiple_events() {
        let ccd = SubStepCcd::default();
        let steps = ccd.build_substeps(&[0.2, 0.5, 0.8]);
        assert_eq!(steps.len(), 4, "steps: {steps:?}");
    }
    #[test]
    fn substep_ccd_duplicate_events_deduped() {
        let ccd = SubStepCcd::default();
        let steps = ccd.build_substeps(&[0.5, 0.5, 0.5]);
        assert_eq!(steps.len(), 2, "duplicates should be deduped: {steps:?}");
    }
    #[test]
    fn substep_ccd_max_substeps_respected() {
        let ccd = SubStepCcd::new(2, 4);
        let steps = ccd.build_substeps(&[0.1, 0.2, 0.3, 0.4, 0.5]);
        assert!(
            steps.len() <= 2,
            "steps={} should be <= max_substeps=2",
            steps.len()
        );
    }
    #[test]
    fn toi_restitution_elastic() {
        let v_a = Vec3::new(-1.0, 0.0, 0.0);
        let v_b = Vec3::new(1.0, 0.0, 0.0);
        let n = Vec3::new(1.0, 0.0, 0.0);
        let (j, va_post, vb_post) = toi_restitution_impulse(v_a, 1.0, v_b, 1.0, n, 1.0);
        assert!(j > 0.0, "impulse should be positive: {j}");
        assert!((va_post.x - 1.0).abs() < 1e-10, "va_post.x={}", va_post.x);
        assert!(
            (vb_post.x - (-1.0)).abs() < 1e-10,
            "vb_post.x={}",
            vb_post.x
        );
    }
    #[test]
    fn toi_restitution_inelastic() {
        let v_a = Vec3::new(-1.0, 0.0, 0.0);
        let v_b = Vec3::new(1.0, 0.0, 0.0);
        let n = Vec3::new(1.0, 0.0, 0.0);
        let (_, va_post, vb_post) = toi_restitution_impulse(v_a, 1.0, v_b, 1.0, n, 0.0);
        assert!(
            (va_post.x - vb_post.x).abs() < 1e-10,
            "bodies should have same velocity: va={} vb={}",
            va_post.x,
            vb_post.x
        );
    }
    #[test]
    fn toi_restitution_separating_no_impulse() {
        let v_a = Vec3::new(2.0, 0.0, 0.0);
        let v_b = Vec3::new(0.0, 0.0, 0.0);
        let n = Vec3::new(1.0, 0.0, 0.0);
        let (j, va_post, _) = toi_restitution_impulse(v_a, 1.0, v_b, 1.0, n, 0.5);
        assert_eq!(j, 0.0, "separating bodies should get zero impulse: {j}");
        assert!((va_post.x - 2.0).abs() < 1e-12);
    }
    #[test]
    fn toi_restitution_static_body() {
        let v_a = Vec3::new(0.0, -5.0, 0.0);
        let v_b = Vec3::zeros();
        let n = Vec3::new(0.0, 1.0, 0.0);
        let (j, va_post, _) = toi_restitution_impulse(v_a, 1.0, v_b, 0.0, n, 0.5);
        assert!(j > 0.0, "impulse should be positive: {j}");
        assert!(
            va_post.y > 0.0,
            "body A should bounce up: va_post.y={}",
            va_post.y
        );
    }
    #[test]
    fn ke_elastic_collision_conserves_energy() {
        let v_a = Vec3::new(1.0, 0.0, 0.0);
        let v_b = Vec3::new(-1.0, 0.0, 0.0);
        let n = Vec3::new(1.0, 0.0, 0.0);
        let (_, va_post, vb_post) = toi_restitution_impulse(v_a, 1.0, v_b, 1.0, n, 1.0);
        let (ke_before, ke_after) = toi_kinetic_energy_change(v_a, 1.0, v_b, 1.0, va_post, vb_post);
        assert!(
            (ke_before - ke_after).abs() < 1e-9,
            "elastic: ke_before={ke_before}, ke_after={ke_after}"
        );
    }
    #[test]
    fn ke_inelastic_collision_loses_energy() {
        let v_a = Vec3::new(-2.0, 0.0, 0.0);
        let v_b = Vec3::new(2.0, 0.0, 0.0);
        let n = Vec3::new(1.0, 0.0, 0.0);
        let (_, va_post, vb_post) = toi_restitution_impulse(v_a, 1.0, v_b, 1.0, n, 0.0);
        let (ke_before, ke_after) = toi_kinetic_energy_change(v_a, 1.0, v_b, 1.0, va_post, vb_post);
        assert!(
            ke_before > ke_after,
            "inelastic: should lose energy: ke_before={ke_before}, ke_after={ke_after}"
        );
    }
    #[test]
    fn ccd_stats_initially_zero() {
        let s = CcdStats::new();
        assert_eq!(s.constraint_count, 0);
        assert_eq!(s.active_impacts, 0);
        assert!(!s.had_active_impacts());
    }
    #[test]
    fn ccd_stats_record_event() {
        let mut s = CcdStats::new();
        s.record_event(0.3, 5.0);
        assert_eq!(s.constraint_count, 1);
        assert_eq!(s.active_impacts, 1);
        assert!((s.earliest_toi - 0.3).abs() < 1e-12);
        assert!((s.total_impulse - 5.0).abs() < 1e-12);
    }
    #[test]
    fn ccd_stats_record_zero_impulse_not_active() {
        let mut s = CcdStats::new();
        s.record_event(0.5, 0.0);
        assert_eq!(s.constraint_count, 1);
        assert_eq!(s.active_impacts, 0);
    }
    #[test]
    fn ccd_stats_merge() {
        let mut a = CcdStats::new();
        a.record_event(0.3, 2.0);
        let mut b = CcdStats::new();
        b.record_event(0.1, 3.0);
        a.merge(&b);
        assert_eq!(a.constraint_count, 2);
        assert!((a.earliest_toi - 0.1).abs() < 1e-12);
        assert!((a.total_impulse - 5.0).abs() < 1e-12);
    }
    #[test]
    fn ccd_stats_earliest_toi_tracks_minimum() {
        let mut s = CcdStats::new();
        s.record_event(0.8, 1.0);
        s.record_event(0.2, 1.0);
        s.record_event(0.5, 1.0);
        assert!(
            (s.earliest_toi - 0.2).abs() < 1e-12,
            "earliest_toi={}",
            s.earliest_toi
        );
    }
    #[test]
    fn filter_ccd_pairs_threshold() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let pairs = vec![
            CcdBroadphasePair {
                body_a: ha,
                body_b: hb,
                toi_upper_bound: 0.2,
            },
            CcdBroadphasePair {
                body_a: ha,
                body_b: hb,
                toi_upper_bound: 0.5,
            },
            CcdBroadphasePair {
                body_a: ha,
                body_b: hb,
                toi_upper_bound: 0.8,
            },
        ];
        let filtered = filter_ccd_pairs(&pairs, 0.5);
        assert_eq!(filtered.len(), 2, "filtered={}", filtered.len());
    }
    #[test]
    fn sort_ccd_pairs_ascending() {
        let mut bodies = RigidBodySet::new();
        let ha = bodies.insert(RigidBody::new(1.0));
        let hb = bodies.insert(RigidBody::new(1.0));
        let mut pairs = vec![
            CcdBroadphasePair {
                body_a: ha,
                body_b: hb,
                toi_upper_bound: 0.8,
            },
            CcdBroadphasePair {
                body_a: ha,
                body_b: hb,
                toi_upper_bound: 0.1,
            },
            CcdBroadphasePair {
                body_a: ha,
                body_b: hb,
                toi_upper_bound: 0.5,
            },
        ];
        sort_ccd_pairs_by_toi(&mut pairs);
        assert!(pairs[0].toi_upper_bound <= pairs[1].toi_upper_bound);
        assert!(pairs[1].toi_upper_bound <= pairs[2].toi_upper_bound);
    }
}
