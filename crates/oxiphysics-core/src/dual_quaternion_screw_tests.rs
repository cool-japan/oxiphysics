// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for ScLERP, the screw power `DualQuaternion::pow` and the
//! `ScrewMotion` decomposition.
//!
//! Split out of `dual_quaternion.rs` to keep that file under the 2000-line
//! limit.  Included via `#[path]` as the body of `dual_quaternion::screw_tests`.

use super::*;
use std::f64::consts::PI;

fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

fn vec3_approx_eq(a: [f64; 3], b: [f64; 3], eps: f64) -> bool {
    approx_eq(a[0], b[0], eps) && approx_eq(a[1], b[1], eps) && approx_eq(a[2], b[2], eps)
}

// ── ScLERP: pure translations, shortest path, exactness ───────────────────

/// Probe points that are not collinear, so two rigid motions agree on all
/// of them only if they are the same motion.  (`dq_approx_eq` negates the
/// real and dual parts independently and can accept different poses.)
const PROBES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

fn same_pose(a: &DualQuaternion, b: &DualQuaternion, eps: f64) -> bool {
    PROBES
        .iter()
        .all(|&p| vec3_approx_eq(a.transform_point(p), b.transform_point(p), eps))
}

fn negated(dq: &DualQuaternion) -> DualQuaternion {
    DualQuaternion {
        real: dq.real.scale(-1.0),
        dual: dq.dual.scale(-1.0),
    }
}

fn translation_dq(t: [f64; 3]) -> DualQuaternion {
    DualQuaternion::from_rotation_translation(Quaternion::identity(), t)
}

#[test]
fn test_sclerp_pure_translation_interpolates_linearly() {
    // The ScrewMotion route used to return the identity for every t because
    // the pure-translation pitch was multiplied by a zero angle.
    let target = translation_dq([1.0, 2.0, 3.0]);
    let identity = DualQuaternion::identity();
    for (t, expected) in [
        (0.0, [0.0, 0.0, 0.0]),
        (0.25, [0.25, 0.5, 0.75]),
        (0.5, [0.5, 1.0, 1.5]),
        (1.0, [1.0, 2.0, 3.0]),
        (2.0, [2.0, 4.0, 6.0]),
        (-1.0, [-1.0, -2.0, -3.0]),
    ] {
        let mid = identity.sclerp(&target, t);
        let tr = mid.to_translation();
        assert!(vec3_approx_eq(tr, expected, 1e-12), "t = {t}: {tr:?}");
        assert!(
            vec3_approx_eq(mid.transform_point([0.0, 0.0, 0.0]), expected, 1e-12),
            "t = {t}"
        );
        assert!(approx_eq(mid.rotation_angle(), 0.0, 1e-7), "t = {t}");
    }
}

#[test]
fn test_sclerp_between_two_translations_from_nonzero_start() {
    let a = translation_dq([1.0, 0.0, -1.0]);
    let b = translation_dq([3.0, 2.0, 1.0]);
    let mid = a.sclerp(&b, 0.5);
    assert!(vec3_approx_eq(mid.to_translation(), [2.0, 1.0, 0.0], 1e-12));
    assert!(same_pose(&a.sclerp(&b, 1.0), &b, 1e-12));
    assert!(same_pose(&a.sclerp(&b, 0.0), &a, 1e-12));
}

#[test]
fn test_sclerp_takes_the_shortest_rotation() {
    // A rotation by 300 degrees about z (scalar part cos 150 < 0) is the
    // rotation by -60 degrees: the midpoint of the shortest path is -30.
    let long = DualQuaternion::from_rotation_translation(
        Quaternion::from_axis_angle([0.0, 0.0, 1.0], 300.0_f64.to_radians()),
        [0.0, 0.0, 0.0],
    );
    assert!(long.real.w < 0.0);
    let mid = DualQuaternion::identity().sclerp(&long, 0.5);
    let p = mid.transform_point([1.0, 0.0, 0.0]);
    let half = 30.0_f64.to_radians();
    assert!(
        vec3_approx_eq(p, [half.cos(), -half.sin(), 0.0], 1e-12),
        "midpoint of the shortest arc expected, got {p:?}"
    );
    // The long way round (+150 degrees) would give a rotation angle of 150.
    assert!(approx_eq(mid.rotation_angle(), 30.0_f64.to_radians(), 1e-9));
}

#[test]
fn test_sclerp_is_invariant_under_the_double_cover() {
    let a = DualQuaternion::from_rotation_translation(
        Quaternion::from_axis_angle([1.0, 0.0, 0.0], 0.4),
        [0.3, -0.2, 0.1],
    );
    let b = DualQuaternion::from_rotation_translation(
        Quaternion::from_axis_angle([0.0, 1.0, 1.0], 2.9),
        [1.0, 2.0, -0.5],
    );
    for t in [0.0, 0.2, 0.5, 0.8, 1.0, 1.5] {
        let direct = a.sclerp(&b, t);
        let via_negative = a.sclerp(&negated(&b), t);
        assert!(same_pose(&direct, &via_negative, 1e-10), "t = {t}");
    }
}

#[test]
fn test_sclerp_matches_closed_form_screw_half_step() {
    // Screw: 90 degrees about z through (1, 0, 0) with pitch 1, i.e. an axial
    // advance of pi/2.  Halving the motion halves both rotation and advance.
    let axis = [0.0, 0.0, 1.0];
    let point = [1.0, 0.0, 0.0];
    let full = DualQuaternion::from_twist(axis, PI / 2.0, 1.0, point);
    let half = DualQuaternion::from_twist(axis, PI / 4.0, 1.0, point);
    let mid = DualQuaternion::identity().sclerp(&full, 0.5);
    assert!(same_pose(&mid, &half, 1e-12));
    let quarter = DualQuaternion::from_twist(axis, PI / 8.0, 1.0, point);
    assert!(same_pose(
        &DualQuaternion::identity().sclerp(&full, 0.25),
        &quarter,
        1e-12
    ));
}

#[test]
fn test_sclerp_general_poses_hit_endpoints_and_square_to_the_full_motion() {
    let p0 = DualQuaternion::from_rotation_translation(
        Quaternion::from_axis_angle([0.0, 1.0, 0.0], 30.0_f64.to_radians()),
        [1.0, -2.0, 0.5],
    );
    let p1 = DualQuaternion::from_rotation_translation(
        Quaternion::from_axis_angle([1.0, 1.0, 0.0], 100.0_f64.to_radians()),
        [-1.0, 0.3, 2.0],
    );
    assert!(same_pose(&p0.sclerp(&p1, 0.0), &p0, 1e-12));
    assert!(same_pose(&p0.sclerp(&p1, 1.0), &p1, 1e-12));
    // M = sclerp(0.5): the step p0 -> M applied twice must reproduce p0 -> p1.
    let mid = p0.sclerp(&p1, 0.5);
    let step = p0.conjugate().mul(&mid);
    let two_steps = p0.mul(&step.mul(&step));
    assert!(same_pose(&two_steps, &p1, 1e-12));
    // Geodesic property: the quarter point is the midpoint of p0 and the midpoint.
    let quarter = p0.sclerp(&p1, 0.25);
    assert!(same_pose(&quarter, &p0.sclerp(&mid, 0.5), 1e-12));
    // The motion is uniform: equal steps along the screw.
    let three_quarter = p0.sclerp(&p1, 0.75);
    assert!(same_pose(&three_quarter, &mid.sclerp(&p1, 0.5), 1e-12));
}

#[test]
fn test_sclerp_matches_se3_matrix_logarithm_reference() {
    // Reference poses computed independently as T0 expm(t logm(T0^-1 T1))
    // with 4x4 homogeneous matrices (scipy.linalg), probed at the origin
    // and the three unit points.
    let p0 = DualQuaternion::from_rotation_translation(
        Quaternion::from_axis_angle([0.0, 1.0, 0.0], 30.0_f64.to_radians()),
        [1.0, -2.0, 0.5],
    );
    let p1 = DualQuaternion::from_rotation_translation(
        Quaternion::from_axis_angle([1.0, 1.0, 0.0], 100.0_f64.to_radians()),
        [-1.0, 0.3, 2.0],
    );
    let at_0_3 = [
        [
            2.183260016987586e-01,
            -1.229879927449269e+00,
            4.759198170806395e-01,
        ],
        [
            9.522846456072748e-01,
            -1.101860023633737e+00,
            -1.911001325903084e-01,
        ],
        [
            3.463459055142905e-01,
            -2.914834916985145e-01,
            7.968918308097832e-01,
        ],
        [
            8.853459513697064e-01,
            -1.550851941178413e+00,
            1.148274896739911e+00,
        ],
    ];
    let at_0_7 = [
        [
            -6.246662874309541e-01,
            -2.311490211121097e-01,
            1.115870856403882e+00,
        ],
        [
            -8.520194088976996e-02,
            1.489059201638511e-01,
            3.645144426942375e-01,
        ],
        [
            -2.446113461549934e-01,
            4.552123771782162e-01,
            1.735924299704277e+00,
        ],
        [
            1.266901262786906e-01,
            -8.512024644125032e-01,
            1.341696601235392e+00,
        ],
    ];
    for (t, expected) in [(0.3, at_0_3), (0.7, at_0_7)] {
        let pose = p0.sclerp(&p1, t);
        for (probe, want) in PROBES.iter().zip(expected.iter()) {
            let got = pose.transform_point(*probe);
            assert!(
                vec3_approx_eq(got, *want, 1e-12),
                "t = {t}, probe {probe:?}: {got:?} vs {want:?}"
            );
        }
    }
}

#[test]
fn test_sclerp_long_rotation_matches_reference_on_the_short_arc() {
    // Rotation by 300 degrees about z (= -60 degrees) with a translation:
    // the reference is the screw motion for -60 degrees, halved.
    let target = DualQuaternion::from_rotation_translation(
        Quaternion::from_axis_angle([0.0, 0.0, 1.0], 300.0_f64.to_radians()),
        [0.5, -1.0, 2.0],
    );
    let expected = [
        [3.839745962155613e-01, -4.330127018922192e-01, 1.0],
        [1.25, -9.330127018922191e-01, 1.0],
        [8.839745962155612e-01, 4.330127018922195e-01, 1.0],
        [3.839745962155613e-01, -4.330127018922192e-01, 2.0],
    ];
    let mid = DualQuaternion::identity().sclerp(&target, 0.5);
    for (probe, want) in PROBES.iter().zip(expected.iter()) {
        let got = mid.transform_point(*probe);
        assert!(
            vec3_approx_eq(got, *want, 1e-12),
            "probe {probe:?}: {got:?} vs {want:?}"
        );
    }
}

#[test]
fn test_sclerp_keeps_perpendicular_translation_for_tiny_rotations() {
    // 1e-7 rad about x with a translation that is mostly perpendicular to
    // the axis: the old screw-axis reconstruction fell back to a zero axis
    // point below ~1e-5 rad and dropped the perpendicular component.
    let target = DualQuaternion::from_rotation_translation(
        Quaternion::from_axis_angle([1.0, 0.0, 0.0], 1e-7),
        [1.0, 2.0, 3.0],
    );
    let mid = DualQuaternion::identity().sclerp(&target, 0.5);
    assert!(
        vec3_approx_eq(mid.to_translation(), [0.5, 1.0, 1.5], 1e-6),
        "got {:?}",
        mid.to_translation()
    );
    // Half of the 1e-7 rad rotation: the quaternion vector part is sin(2.5e-8).
    assert!(approx_eq(mid.real.x, 2.5e-8, 1e-15), "got {}", mid.real.x);
    let end = DualQuaternion::identity().sclerp(&target, 1.0);
    assert!(same_pose(&end, &target, 1e-12));
}

/// Minimal deterministic xorshift generator for reproducible random poses.
struct XorShift(u64);

impl XorShift {
    fn next_f64(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }

    fn pose(&mut self) -> DualQuaternion {
        let axis = [
            self.uniform(-1.0, 1.0),
            self.uniform(-1.0, 1.0),
            self.uniform(-1.0, 1.0) + 0.1,
        ];
        let angle = self.uniform(-3.1, 3.1);
        let t = [
            self.uniform(-5.0, 5.0),
            self.uniform(-5.0, 5.0),
            self.uniform(-5.0, 5.0),
        ];
        DualQuaternion::from_rotation_translation(Quaternion::from_axis_angle(axis, angle), t)
    }
}

#[test]
fn test_sclerp_random_poses_satisfy_the_geodesic_identities() {
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    for case in 0..200 {
        let p0 = rng.pose();
        let p1 = rng.pose();
        // Half of the cases use the opposite sign of the same target motion.
        let target = if case % 2 == 0 { p1 } else { negated(&p1) };
        assert!(
            same_pose(&p0.sclerp(&target, 0.0), &p0, 1e-9),
            "case {case}"
        );
        assert!(
            same_pose(&p0.sclerp(&target, 1.0), &p1, 1e-9),
            "case {case}"
        );
        // Midpoint M: stepping p0 -> M twice must reproduce p0 -> p1.
        let mid = p0.sclerp(&target, 0.5);
        let step = p0.conjugate().mul(&mid);
        assert!(
            same_pose(&p0.mul(&step.mul(&step)), &p1, 1e-9),
            "case {case}: midpoint does not halve the motion"
        );
        // Result must not depend on the sign of the target.
        let t = rng.uniform(-0.5, 1.5);
        assert!(
            same_pose(&p0.sclerp(&p1, t), &p0.sclerp(&negated(&p1), t), 1e-9),
            "case {case}: t = {t}"
        );
        // Geodesic additivity: p0 -> (t1 + t2) equals p0 -> t1 followed by the
        // fraction t2 / (1 - t1) of the rest of the way.
        let (t1, t2) = (0.3, 0.45);
        let first = p0.sclerp(&p1, t1);
        let composed = first.sclerp(&p1, t2 / (1.0 - t1));
        assert!(
            same_pose(&composed, &p0.sclerp(&p1, t1 + t2), 1e-9),
            "case {case}: additivity"
        );
    }
}

#[test]
fn test_screw_motion_round_trips_random_poses() {
    let mut rng = XorShift(0xD1B5_4A32_D192_ED03);
    for case in 0..200 {
        let pose = rng.pose();
        let screw = ScrewMotion::from_dual_quaternion(&pose);
        assert!(
            same_pose(&screw.to_dual_quaternion(), &pose, 1e-8),
            "case {case}: angle {}",
            screw.angle
        );
    }
}

#[test]
fn test_pow_is_a_one_parameter_group() {
    let axis = [2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0];
    let dq = DualQuaternion::from_twist(axis, 1.1, 0.4, [0.2, 0.7, -0.3]);
    assert!(same_pose(&dq.pow(0.0), &DualQuaternion::identity(), 1e-12));
    assert!(same_pose(&dq.pow(1.0), &dq, 1e-12));
    assert!(same_pose(&dq.pow(2.0), &dq.mul(&dq), 1e-12));
    assert!(same_pose(&dq.pow(-1.0), &dq.conjugate(), 1e-12));
    assert!(same_pose(
        &dq.pow(0.3).mul(&dq.pow(0.45)),
        &dq.pow(0.75),
        1e-12
    ));
    // pow keeps the dual quaternion a unit rigid motion.
    let unit = dq.pow(0.37);
    assert!(approx_eq(unit.real.norm(), 1.0, 1e-12));
    let orthogonality = unit.real.dot(&unit.dual);
    assert!(orthogonality.abs() < 1e-12, "real . dual = {orthogonality}");
}

#[test]
fn test_pow_pure_translation_and_pure_rotation() {
    let shifted = translation_dq([2.0, -4.0, 6.0]);
    assert!(vec3_approx_eq(
        shifted.pow(0.5).to_translation(),
        [1.0, -2.0, 3.0],
        1e-12
    ));
    let spin = DualQuaternion::from_twist([0.0, 0.0, 1.0], PI / 2.0, 0.0, [0.0; 3]);
    let half = spin.pow(0.5);
    let p = half.transform_point([1.0, 0.0, 0.0]);
    assert!(vec3_approx_eq(
        p,
        [(PI / 4.0).cos(), (PI / 4.0).sin(), 0.0],
        1e-12
    ));
    assert!(vec3_approx_eq(half.to_translation(), [0.0; 3], 1e-12));
}

// ── ScrewMotion round trips ──────────────────────────────────────────────

#[test]
fn test_screw_motion_round_trips_pure_translation() {
    let dq = translation_dq([1.0, 2.0, 3.0]);
    let sm = ScrewMotion::from_dual_quaternion(&dq);
    assert_eq!(sm.angle, 0.0);
    assert!(approx_eq(sm.pitch, 14.0_f64.sqrt(), 1e-12));
    let back = sm.to_dual_quaternion();
    assert!(
        same_pose(&back, &dq, 1e-12),
        "got {:?}",
        back.to_translation()
    );
}

#[test]
fn test_screw_motion_round_trips_pure_translation_double_cover() {
    // -(1 + eps d) encodes the same translation.
    let dq = negated(&translation_dq([-0.5, 0.25, 4.0]));
    let sm = ScrewMotion::from_dual_quaternion(&dq);
    assert!(sm.axis_dir.iter().all(|c| c.is_finite()));
    assert!(same_pose(
        &sm.to_dual_quaternion(),
        &translation_dq([-0.5, 0.25, 4.0]),
        1e-12
    ));
}

#[test]
fn test_screw_motion_round_trips_small_rotations_with_perpendicular_translation() {
    // For rotations between ~1e-8 and ~1e-5 rad the old decomposition
    // (acos-based angle, axis point forced to zero below a threshold)
    // dropped the translation perpendicular to the axis: the recomposed
    // pose was off by |τ⊥| = 3.6 here.
    for angle in [1e-3, 1e-5, 1e-6, 1e-7, 1e-9] {
        let dq = DualQuaternion::from_rotation_translation(
            Quaternion::from_axis_angle([1.0, 0.0, 0.0], angle),
            [1.0, 2.0, 3.0],
        );
        let sm = ScrewMotion::from_dual_quaternion(&dq);
        assert!(
            ((sm.angle - angle) / angle).abs() < 1e-9,
            "angle {angle}: recovered {}",
            sm.angle
        );
        assert!(
            same_pose(&sm.to_dual_quaternion(), &dq, 1e-7),
            "angle {angle}: round trip lost the translation"
        );
    }
}

#[test]
fn test_screw_motion_round_trips_general_screw() {
    let dq = DualQuaternion::from_twist([0.0, 0.0, 1.0], 1.2, 0.5, [1.0, -1.0, 0.0]);
    let sm = ScrewMotion::from_dual_quaternion(&dq);
    assert!(approx_eq(sm.angle, 1.2, 1e-12));
    assert!(approx_eq(sm.pitch, 0.5, 1e-12));
    assert!(same_pose(&sm.to_dual_quaternion(), &dq, 1e-10));
}
