//! Regression tests for the numerical 6-DOF Jacobian.
//!
//! The angular rows are the world-frame directions of the joint axes. An earlier
//! version extracted them from `R₀ᵀR₁`, which is the body-frame rotation, and
//! read the skew-symmetric entries with the wrong sign.

use super::*;

const TOLERANCE: f64 = 1e-7;

/// World-frame joint axes of `arm` at `q` as a flat row-major `3 × dof` array
/// (zero columns for prismatic joints).
fn analytic_angular_rows(arm: &RobotArm, q: &[f64]) -> Vec<f64> {
    let frames = arm.forward_frames(q);
    let n = arm.dof();
    let mut rows = vec![0.0; 3 * n];
    let mut column = 0;
    for (i, link) in arm.links.iter().enumerate() {
        let axis = match link.joint_type {
            JointType::Revolute => [frames[i][2], frames[i][6], frames[i][10]],
            JointType::Prismatic => [0.0; 3],
            JointType::Fixed => continue,
        };
        for (row, component) in axis.iter().enumerate() {
            rows[row * n + column] = *component;
        }
        column += 1;
    }
    rows
}

fn assert_rows_close(actual: &[f64], expected: &[f64], what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: length");
    for (index, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (a - e).abs() < TOLERANCE,
            "{what}: entry {index} is {a}, expected {e}"
        );
    }
}

/// Checks both blocks of the numerical Jacobian against the analytic ones.
fn assert_matches_analytic(arm: &RobotArm, q: &[f64]) {
    let n = arm.dof();
    let numerical = arm
        .jacobian_6dof_numerical(q)
        .expect("configuration has one value per joint");
    assert_eq!(numerical.len(), 6 * n);
    assert_rows_close(
        &numerical[..3 * n],
        &arm.jacobian_analytical(q),
        "linear rows",
    );
    assert_rows_close(
        &numerical[3 * n..],
        &analytic_angular_rows(arm, q),
        "angular rows",
    );
}

#[test]
fn planar_arm_angular_rows_are_the_positive_z_axis() {
    let arm = build_planar_3dof([0.5, 0.4, 0.3]);
    let q = [0.3, -0.7, 0.5];
    let j = arm.jacobian_6dof_numerical(&q).expect("three joint values");
    let n = 3;
    for joint in 0..n {
        assert!(j[3 * n + joint].abs() < TOLERANCE, "wx of joint {joint}");
        assert!(j[4 * n + joint].abs() < TOLERANCE, "wy of joint {joint}");
        assert!(
            (j[5 * n + joint] - 1.0).abs() < TOLERANCE,
            "wz of joint {joint} is {}",
            j[5 * n + joint]
        );
    }
}

#[test]
fn planar_arm_linear_rows_match_the_analytic_jacobian() {
    let arm = build_planar_3dof([0.5, 0.4, 0.3]);
    assert_matches_analytic(&arm, &[0.3, -0.7, 0.5]);
    assert_matches_analytic(&arm, &[0.0, 0.0, 0.0]);
    assert_matches_analytic(&arm, &[-1.2, 2.1, 0.9]);
}

#[test]
fn second_joint_of_a_twisted_arm_reports_its_world_frame_axis() {
    // Joint 1 turns about the world z axis; its link is twisted by 90 degrees so
    // that the axis of joint 2 is z1 = Rz(q1) (0, -1, 0) = (sin q1, -cos q1, 0),
    // which differs from the body-frame axis (0, 0, 1) of the end effector.
    let lim = JointLimits::new(-PI, PI, 3.0, 50.0);
    let arm = RobotArm::new(
        "twisted",
        vec![
            RobotLink::new(
                "a",
                DhParam::new(0.0, 0.0, 0.5, PI / 2.0),
                JointType::Revolute,
                Some(lim),
                0.03,
                0.25,
            ),
            RobotLink::new(
                "b",
                DhParam::new(0.0, 0.0, 0.4, 0.0),
                JointType::Revolute,
                Some(lim),
                0.03,
                0.2,
            ),
        ],
    );
    let q1 = 0.6;
    let j = arm
        .jacobian_6dof_numerical(&[q1, 0.35])
        .expect("two joint values");
    let n = 2;
    let expected = [[0.0, 0.0, 1.0], [q1.sin(), -q1.cos(), 0.0]];
    for (joint, axis) in expected.iter().enumerate() {
        for (row, component) in axis.iter().enumerate() {
            let value = j[(3 + row) * n + joint];
            assert!(
                (value - component).abs() < TOLERANCE,
                "joint {joint}, angular row {row}: {value} vs {component}"
            );
        }
    }
    assert_matches_analytic(&arm, &[q1, 0.35]);
}

#[test]
fn puma_like_arm_matches_the_analytic_jacobian() {
    let arm = build_6dof_puma_like();
    assert_matches_analytic(&arm, &[0.1, -0.4, 0.8, 0.3, -0.6, 0.9]);
    assert_matches_analytic(&arm, &[1.3, 0.2, -0.5, -1.1, 0.4, -0.2]);
    assert_matches_analytic(&arm, &[0.0; 6]);
}

#[test]
fn prismatic_fixed_and_base_transform_are_handled_in_the_world_frame() {
    let lim = JointLimits::new(-PI, PI, 3.0, 50.0);
    let base = mat4_mul(&rot_z(0.9), &mat4_mul(&rot_x(0.4), &trans(0.2, -0.1, 0.3)));
    let arm = RobotArm::new(
        "mixed",
        vec![
            RobotLink::new(
                "shoulder",
                DhParam::new(0.2, 0.1, 0.0, PI / 2.0),
                JointType::Revolute,
                Some(lim),
                0.03,
                0.1,
            ),
            RobotLink::new(
                "mount",
                DhParam::new(0.0, 0.05, 0.1, 0.0),
                JointType::Fixed,
                None,
                0.03,
                0.1,
            ),
            RobotLink::new(
                "slide",
                DhParam::new(0.0, 0.2, 0.0, -PI / 2.0),
                JointType::Prismatic,
                Some(lim),
                0.03,
                0.1,
            ),
            RobotLink::new(
                "wrist",
                DhParam::new(0.1, 0.0, 0.3, 0.0),
                JointType::Revolute,
                Some(lim),
                0.03,
                0.1,
            ),
        ],
    )
    .with_base_transform(base);
    assert_eq!(arm.dof(), 3);
    let q = [0.5, 0.15, -0.8];
    assert_matches_analytic(&arm, &q);

    // The prismatic column must not rotate the end effector.
    let j = arm.jacobian_6dof_numerical(&q).expect("three joint values");
    for row in 3..6 {
        assert!(j[row * 3 + 1].abs() < TOLERANCE, "angular row {row}");
    }
}

#[test]
fn wrong_configuration_length_is_rejected() {
    let arm = build_planar_3dof([0.5, 0.4, 0.3]);
    for q in [&[0.1, 0.2][..], &[0.1, 0.2, 0.3, 0.4][..], &[][..]] {
        match arm.jacobian_6dof_numerical(q) {
            Err(Error::LengthMismatch { expected, actual }) => {
                assert_eq!(expected, 3);
                assert_eq!(actual, q.len());
            }
            other => panic!("expected a length mismatch, got {other:?}"),
        }
    }
}

#[test]
fn jacobian_size_follows_the_degrees_of_freedom_not_the_link_count() {
    let lim = JointLimits::new(-PI, PI, 3.0, 50.0);
    let arm = RobotArm::new(
        "one_joint_one_fixed",
        vec![
            RobotLink::new(
                "joint",
                DhParam::new(0.0, 0.0, 0.5, 0.0),
                JointType::Revolute,
                Some(lim),
                0.03,
                0.25,
            ),
            RobotLink::new(
                "tool",
                DhParam::new(0.0, 0.0, 0.2, 0.0),
                JointType::Fixed,
                None,
                0.03,
                0.1,
            ),
        ],
    );
    assert_eq!(arm.links.len(), 2);
    let j = arm
        .jacobian_6dof_numerical(&[0.4])
        .expect("one joint value");
    assert_eq!(j.len(), 6);
    // Rotating a 0.7 m arm about z moves the tip at 0.7 m/rad.
    let tip_speed = (j[0] * j[0] + j[1] * j[1]).sqrt();
    assert!((tip_speed - 0.7).abs() < TOLERANCE, "tip speed {tip_speed}");
    assert!((j[5] - 1.0).abs() < TOLERANCE);
}
