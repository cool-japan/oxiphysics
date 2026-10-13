//! Regression tests for the XPBD dihedral / pressure constraints, the
//! tetrahedral volume constraints of [`PbdSimulation`] and the PBD sub-step
//! solver.

use super::functions::{
    pbd2_add, pbd2_cross, pbd2_dot, pbd2_norm, pbd2_scale, pbd2_sub, tet_volume,
};
use super::types::*;
use oxiphysics_core::math::Vec3;

const TIGHT: f64 = 1e-12;

fn weighted_centroid(positions: &[[f64; 3]], inv_masses: &[f64]) -> [f64; 3] {
    let mut total_mass = 0.0;
    let mut acc = [0.0; 3];
    for (p, &w) in positions.iter().zip(inv_masses.iter()) {
        let m = 1.0 / w;
        total_mass += m;
        acc = pbd2_add(acc, pbd2_scale(*p, m));
    }
    pbd2_scale(acc, 1.0 / total_mass)
}

fn max_component_diff(a: [f64; 3], b: [f64; 3]) -> f64 {
    pbd2_sub(a, b).iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

/// Central-difference gradient of the dihedral angle.
fn numeric_angle_gradient(c: &XpbdDihedralBending, positions: &[[f64; 3]; 4]) -> [[f64; 3]; 4] {
    let h = 1e-6;
    let mut grad = [[0.0; 3]; 4];
    for (k, &idx) in c.indices.iter().enumerate() {
        for axis in 0..3 {
            let mut plus = *positions;
            let mut minus = *positions;
            plus[idx][axis] += h;
            minus[idx][axis] -= h;
            grad[k][axis] = (c.compute_angle(&plus) - c.compute_angle(&minus)) / (2.0 * h);
        }
    }
    grad
}

/// Central-difference gradient of the tetrahedron volume.
fn numeric_volume_gradient(c: &XpbdPressureConstraint, positions: &[[f64; 3]; 4]) -> [[f64; 3]; 4] {
    let h = 1e-6;
    let mut grad = [[0.0; 3]; 4];
    for (k, &idx) in c.tet_indices.iter().enumerate() {
        for axis in 0..3 {
            let mut plus = *positions;
            let mut minus = *positions;
            plus[idx][axis] += h;
            minus[idx][axis] -= h;
            grad[k][axis] = (c.current_volume(&plus) - c.current_volume(&minus)) / (2.0 * h);
        }
    }
    grad
}

/// Poses of the hinge covering both fold directions and generic skew poses.
fn dihedral_poses() -> Vec<[[f64; 3]; 4]> {
    vec![
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.0, 1.0],
        ],
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.0, -1.0],
        ],
        [
            [0.1, -0.2, 0.3],
            [1.3, 0.2, -0.1],
            [0.4, 1.1, 0.5],
            [0.7, -0.8, 0.9],
        ],
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.5, 0.0],
            [1.0, 1.5, 0.3],
            [0.8, -0.4, 1.2],
        ],
        [
            [-0.3, 0.2, 0.1],
            [0.9, -0.1, 0.4],
            [0.2, 0.9, -0.6],
            [0.5, 0.3, -1.1],
        ],
    ]
}

/// Hinge folded to exactly `angle` radians (`p2` fixed, `p3` rotated).
fn folded_hinge(angle: f64) -> [[f64; 3]; 4] {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, -angle.cos(), angle.sin()],
    ]
}

// ── XpbdDihedralBending ──────────────────────────────────────────────────────

#[test]
fn dihedral_gradients_match_central_differences() {
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 0.0, 0.0);
    for pose in dihedral_poses() {
        let angle = c.compute_angle(&pose);
        assert!(
            angle > 0.05 && angle < std::f64::consts::PI - 0.05,
            "test pose must stay away from the flat / folded singularities, angle={angle}"
        );
        let analytic = c
            .angle_gradients(&pose)
            .expect("non-degenerate pose has gradients");
        let numeric = numeric_angle_gradient(&c, &pose);
        for k in 0..4 {
            for axis in 0..3 {
                assert!(
                    (analytic[k][axis] - numeric[k][axis]).abs() < 1e-6,
                    "pose {pose:?}: d(angle)/d(p{k}[{axis}]) analytic={} numeric={}",
                    analytic[k][axis],
                    numeric[k][axis]
                );
            }
        }
    }
}

#[test]
fn dihedral_gradients_are_translation_and_rotation_invariant() {
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 0.0, 0.0);
    for pose in dihedral_poses() {
        let g = c.angle_gradients(&pose).expect("gradients");
        let mut force_sum = [0.0; 3];
        let mut torque_sum = [0.0; 3];
        for k in 0..4 {
            force_sum = pbd2_add(force_sum, g[k]);
            torque_sum = pbd2_add(torque_sum, pbd2_cross(pose[k], g[k]));
        }
        assert!(
            pbd2_norm(force_sum) < 1e-10,
            "gradients must sum to zero, got {force_sum:?}"
        );
        assert!(
            pbd2_norm(torque_sum) < 1e-10,
            "gradients must be torque free, got {torque_sum:?}"
        );
    }
}

#[test]
fn dihedral_gradients_are_none_for_degenerate_poses() {
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 0.0, 0.0);
    // Zero-length shared edge.
    let collapsed_edge = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, 0.0, 1.0],
    ];
    assert!(c.angle_gradients(&collapsed_edge).is_none());
    // Apex of triangle A lies on the shared edge line (zero area).
    let collinear_apex = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.5, 0.0, 1.0],
    ];
    assert!(c.angle_gradients(&collinear_apex).is_none());
}

#[test]
fn dihedral_project_folds_to_rest_angle_and_preserves_centroid() {
    let mut positions = folded_hinge(1.0);
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 0.5, 0.0);
    assert!((c.compute_angle(&positions) - 1.0).abs() < 1e-12);
    let inv_masses = [1.0; 4];
    let centroid_before = weighted_centroid(&positions, &inv_masses);
    let mut lambda = 0.0;
    for _ in 0..50 {
        c.project(&mut positions, &inv_masses, &mut lambda, 1.0 / 60.0);
    }
    let angle = c.compute_angle(&positions);
    assert!(
        (angle - 0.5).abs() < 1e-3,
        "angle should converge to the rest angle 0.5, got {angle}"
    );
    let centroid_after = weighted_centroid(&positions, &inv_masses);
    assert!(
        max_component_diff(centroid_before, centroid_after) < TIGHT,
        "the hinge must fold about its centroid, drift {:?}",
        pbd2_sub(centroid_after, centroid_before)
    );
}

#[test]
fn dihedral_project_converges_for_both_fold_directions() {
    let inv_masses = [1.0; 4];
    for pose in dihedral_poses() {
        for rest in [0.3_f64, 1.2, 2.4] {
            let mut positions = pose;
            let c = XpbdDihedralBending::new([0, 1, 2, 3], rest, 0.0);
            let mut lambda = 0.0;
            for _ in 0..200 {
                c.project(&mut positions, &inv_masses, &mut lambda, 1.0 / 60.0);
            }
            let angle = c.compute_angle(&positions);
            assert!(
                (angle - rest).abs() < 1e-6,
                "pose {pose:?} rest {rest}: converged angle {angle}"
            );
        }
    }
}

#[test]
fn dihedral_project_flattens_to_zero_rest_angle() {
    let mut positions = folded_hinge(0.8);
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 0.0, 0.0);
    let inv_masses = [1.0; 4];
    let mut lambda = 0.0;
    for _ in 0..200 {
        c.project(&mut positions, &inv_masses, &mut lambda, 1.0 / 60.0);
    }
    assert!(
        c.compute_angle(&positions) < 1e-3,
        "folded hinge should flatten, angle={}",
        c.compute_angle(&positions)
    );
}

#[test]
fn dihedral_project_keeps_pinned_particles_fixed() {
    let mut positions = folded_hinge(1.0);
    let before = positions;
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 0.4, 0.0);
    let inv_masses = [0.0, 0.0, 1.0, 1.0];
    let mut lambda = 0.0;
    for _ in 0..100 {
        c.project(&mut positions, &inv_masses, &mut lambda, 1.0 / 60.0);
    }
    assert_eq!(positions[0], before[0]);
    assert_eq!(positions[1], before[1]);
    assert!(
        (c.compute_angle(&positions) - 0.4).abs() < 1e-6,
        "free apexes must still reach the rest angle"
    );
}

#[test]
fn dihedral_project_conserves_mass_weighted_centroid() {
    let mut positions = dihedral_poses()[2];
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 1.4, 0.0);
    let inv_masses = [1.0, 2.0, 0.5, 4.0];
    let before = weighted_centroid(&positions, &inv_masses);
    let mut lambda = 0.0;
    for _ in 0..30 {
        c.project(&mut positions, &inv_masses, &mut lambda, 1.0 / 60.0);
    }
    let after = weighted_centroid(&positions, &inv_masses);
    assert!(max_component_diff(before, after) < 1e-11);
    assert!((c.compute_angle(&positions) - 1.4).abs() < 1e-6);
}

#[test]
fn dihedral_project_compliance_softens_correction() {
    let inv_masses = [1.0; 4];
    let dt = 1.0 / 60.0;
    let rigid = XpbdDihedralBending::new([0, 1, 2, 3], 0.5, 0.0);
    let soft = XpbdDihedralBending::new([0, 1, 2, 3], 0.5, 1e-2);
    let mut pos_rigid = folded_hinge(1.0);
    let mut pos_soft = folded_hinge(1.0);
    let (mut lam_rigid, mut lam_soft) = (0.0, 0.0);
    rigid.project(&mut pos_rigid, &inv_masses, &mut lam_rigid, dt);
    soft.project(&mut pos_soft, &inv_masses, &mut lam_soft, dt);
    let err_rigid = (rigid.compute_angle(&pos_rigid) - 0.5).abs();
    let err_soft = (soft.compute_angle(&pos_soft) - 0.5).abs();
    assert!(
        err_rigid < err_soft,
        "rigid error {err_rigid} must be smaller than compliant error {err_soft}"
    );
    assert!(err_soft < 0.5, "compliant constraint must still correct");
    assert!(
        lam_rigid < 0.0 && lam_soft < 0.0,
        "folding more than rest gives a negative multiplier"
    );
    assert!(lam_soft.abs() < lam_rigid.abs());
}

#[test]
fn dihedral_project_ignores_degenerate_pose() {
    let mut positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.5, 0.0, 1.0],
    ];
    let before = positions;
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 0.7, 0.0);
    let mut lambda = 0.0;
    c.project(&mut positions, &[1.0; 4], &mut lambda, 1.0 / 60.0);
    assert_eq!(positions, before);
    assert_eq!(lambda, 0.0);
}

#[test]
fn dihedral_project_with_zero_dt_is_finite() {
    let mut positions = folded_hinge(1.0);
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 0.5, 0.0);
    let mut lambda = 0.0;
    c.project(&mut positions, &[1.0; 4], &mut lambda, 0.0);
    assert!(positions.iter().flatten().all(|v| v.is_finite()));
    assert!(lambda.is_finite());
}

#[test]
fn dihedral_zero_time_step_treats_a_compliant_constraint_as_rigid() {
    // dt = 0 carries no stiffness information: the compliance must not turn
    // into an infinite alpha-tilde (inf * lambda = NaN once lambda != 0).
    let inv_masses = [1.0; 4];
    let soft = XpbdDihedralBending::new([0, 1, 2, 3], 0.5, 50.0);
    let rigid = XpbdDihedralBending::new([0, 1, 2, 3], 0.5, 0.0);
    let mut pos_soft = folded_hinge(1.0);
    let mut pos_rigid = folded_hinge(1.0);
    let (mut lam_soft, mut lam_rigid) = (3.0, 3.0);
    soft.project(&mut pos_soft, &inv_masses, &mut lam_soft, 0.0);
    rigid.project(&mut pos_rigid, &inv_masses, &mut lam_rigid, 1.0 / 60.0);
    assert!(pos_soft.iter().flatten().all(|v| v.is_finite()));
    assert_eq!(pos_soft, pos_rigid);
    assert_eq!(lam_soft, lam_rigid);
}

#[test]
fn dihedral_project_unfolds_a_flat_hinge() {
    // theta = 0 is the kink of the unsigned angle: either side is a valid
    // direction, and the constraint must pick one and converge to the rest angle.
    let mut positions = folded_hinge(0.0);
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 0.6, 0.0);
    assert!(c.compute_angle(&positions) < TIGHT);
    let inv_masses = [1.0; 4];
    let centroid_before = weighted_centroid(&positions, &inv_masses);
    let mut lambda = 0.0;
    for _ in 0..100 {
        c.project(&mut positions, &inv_masses, &mut lambda, 1.0 / 60.0);
    }
    let angle = c.compute_angle(&positions);
    assert!(
        (angle - 0.6).abs() < 1e-6,
        "flat hinge must open to the rest angle, got {angle}"
    );
    let centroid_after = weighted_centroid(&positions, &inv_masses);
    assert!(max_component_diff(centroid_before, centroid_after) < 1e-11);
}

#[test]
fn dihedral_project_opens_a_hinge_folded_flat_onto_itself() {
    // theta = pi: both triangles coincide.  Both sides are equivalent, the
    // constraint opens the hinge towards the rest angle.
    let mut positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, 1.0, 0.0],
    ];
    let c = XpbdDihedralBending::new([0, 1, 2, 3], 2.5, 0.0);
    assert!((c.compute_angle(&positions) - std::f64::consts::PI).abs() < 1e-12);
    let inv_masses = [1.0; 4];
    let mut lambda = 0.0;
    for _ in 0..100 {
        c.project(&mut positions, &inv_masses, &mut lambda, 1.0 / 60.0);
    }
    let angle = c.compute_angle(&positions);
    assert!(
        (angle - 2.5).abs() < 1e-6,
        "folded hinge must open to the rest angle, got {angle}"
    );
}

// ── XpbdPressureConstraint ───────────────────────────────────────────────────

fn pressure_tet(apex_z: f64) -> [[f64; 3]; 4] {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, 0.5, apex_z],
    ]
}

#[test]
fn pressure_current_volume_matches_tet_volume() {
    let p = pressure_tet(1.0);
    let c = XpbdPressureConstraint::new([0, 1, 2, 3], 0.0, 0.0);
    let expected = tet_volume(p[0], p[1], p[2], p[3]);
    assert!((c.current_volume(&p) - expected).abs() < TIGHT);
    assert!((expected - 1.0 / 6.0).abs() < TIGHT);
}

#[test]
fn pressure_gradients_match_central_differences() {
    let c = XpbdPressureConstraint::new([0, 1, 2, 3], 0.0, 0.0);
    for pose in dihedral_poses() {
        let analytic = c.volume_gradients(&pose);
        let numeric = numeric_volume_gradient(&c, &pose);
        for k in 0..4 {
            for axis in 0..3 {
                assert!(
                    (analytic[k][axis] - numeric[k][axis]).abs() < 1e-8,
                    "dV/dp{k}[{axis}] analytic={} numeric={}",
                    analytic[k][axis],
                    numeric[k][axis]
                );
            }
        }
        let sum = analytic.iter().fold([0.0; 3], |acc, g| pbd2_add(acc, *g));
        assert!(pbd2_norm(sum) < TIGHT, "volume gradients must sum to zero");
    }
}

#[test]
fn pressure_constraint_restores_rest_volume() {
    let rest_volume = tet_volume(
        pressure_tet(1.0)[0],
        pressure_tet(1.0)[1],
        pressure_tet(1.0)[2],
        pressure_tet(1.0)[3],
    );
    let mut positions = pressure_tet(3.0);
    let c = XpbdPressureConstraint::new([0, 1, 2, 3], rest_volume, 0.0);
    let initial_error = (c.current_volume(&positions) - rest_volume).abs();
    assert!(
        initial_error > 0.3,
        "lifted apex must violate the constraint"
    );
    let inv_masses = [1.0; 4];
    let mut lambda = 0.0;
    for _ in 0..20 {
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
    }
    let volume = c.current_volume(&positions);
    assert!(
        (volume - rest_volume).abs() < 1e-6 * rest_volume,
        "volume {volume} must reach the rest volume {rest_volume}"
    );
}

#[test]
fn pressure_constraint_inflates_a_compressed_tetrahedron() {
    let rest = pressure_tet(1.0);
    let rest_volume = tet_volume(rest[0], rest[1], rest[2], rest[3]);
    let mut positions = rest.map(|p| pbd2_scale(p, 0.5));
    let c = XpbdPressureConstraint::new([0, 1, 2, 3], rest_volume, 0.0);
    assert!(c.current_volume(&positions) < 0.2 * rest_volume);
    let inv_masses = [1.0; 4];
    let mut lambda = 0.0;
    for _ in 0..40 {
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
    }
    let volume = c.current_volume(&positions);
    assert!(
        (volume - rest_volume).abs() < 1e-6 * rest_volume,
        "compressed tet must re-inflate, volume {volume} vs {rest_volume}"
    );
}

#[test]
fn pressure_constraint_at_rest_changes_nothing() {
    let mut positions = pressure_tet(1.0);
    let before = positions;
    let rest_volume = tet_volume(positions[0], positions[1], positions[2], positions[3]);
    let c = XpbdPressureConstraint::new([0, 1, 2, 3], rest_volume, 0.0);
    let mut lambda = 0.0;
    c.project(&mut positions, &[1.0; 4], &mut lambda, 0.01);
    assert!(lambda.abs() < TIGHT, "lambda={lambda}");
    for k in 0..4 {
        assert!(max_component_diff(positions[k], before[k]) < TIGHT);
    }
}

#[test]
fn pressure_constraint_conserves_mass_weighted_centroid() {
    let rest = pressure_tet(1.0);
    let rest_volume = tet_volume(rest[0], rest[1], rest[2], rest[3]);
    let mut positions = pressure_tet(2.5);
    let inv_masses = [1.0, 0.25, 2.0, 0.5];
    let before = weighted_centroid(&positions, &inv_masses);
    let c = XpbdPressureConstraint::new([0, 1, 2, 3], rest_volume, 0.0);
    let mut lambda = 0.0;
    for _ in 0..20 {
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
    }
    let after = weighted_centroid(&positions, &inv_masses);
    assert!(max_component_diff(before, after) < 1e-11);
    let volume = c.current_volume(&positions);
    assert!((volume - rest_volume).abs() < 1e-6 * rest_volume);
}

#[test]
fn pressure_constraint_respects_pinned_particles_and_masses() {
    let rest = pressure_tet(1.0);
    let rest_volume = tet_volume(rest[0], rest[1], rest[2], rest[3]);
    let mut positions = pressure_tet(2.0);
    let before = positions;
    // Particles 0..=2 are pinned; only the apex may move.
    let inv_masses = [0.0, 0.0, 0.0, 1.0];
    let c = XpbdPressureConstraint::new([0, 1, 2, 3], rest_volume, 0.0);
    let mut lambda = 0.0;
    for _ in 0..10 {
        c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
    }
    for k in 0..3 {
        assert_eq!(positions[k], before[k]);
    }
    assert!((c.current_volume(&positions) - rest_volume).abs() < 1e-9);
    // The apex moves along the gradient direction: only its z coordinate changes.
    assert!((positions[3][0] - before[3][0]).abs() < 1e-12);
    assert!((positions[3][1] - before[3][1]).abs() < 1e-12);
    assert!(positions[3][2] < before[3][2]);
}

#[test]
fn pressure_constraint_compliance_softens_correction() {
    let rest = pressure_tet(1.0);
    let rest_volume = tet_volume(rest[0], rest[1], rest[2], rest[3]);
    let dt = 0.01;
    let mut rigid_pos = pressure_tet(3.0);
    let mut soft_pos = pressure_tet(3.0);
    let rigid = XpbdPressureConstraint::new([0, 1, 2, 3], rest_volume, 0.0);
    let soft = XpbdPressureConstraint::new([0, 1, 2, 3], rest_volume, 1e-3);
    let (mut lam_rigid, mut lam_soft) = (0.0, 0.0);
    rigid.project(&mut rigid_pos, &[1.0; 4], &mut lam_rigid, dt);
    soft.project(&mut soft_pos, &[1.0; 4], &mut lam_soft, dt);
    let err_rigid = (rigid.current_volume(&rigid_pos) - rest_volume).abs();
    let err_soft = (soft.current_volume(&soft_pos) - rest_volume).abs();
    assert!(err_rigid < err_soft, "rigid {err_rigid} vs soft {err_soft}");
    assert!(lam_soft.abs() < lam_rigid.abs());
}

#[test]
fn pressure_zero_time_step_treats_a_compliant_constraint_as_rigid() {
    let rest = pressure_tet(1.0);
    let rest_volume = tet_volume(rest[0], rest[1], rest[2], rest[3]);
    let soft = XpbdPressureConstraint::new([0, 1, 2, 3], rest_volume, 50.0);
    let rigid = XpbdPressureConstraint::new([0, 1, 2, 3], rest_volume, 0.0);
    let mut pos_soft = pressure_tet(3.0);
    let mut pos_rigid = pressure_tet(3.0);
    let (mut lam_soft, mut lam_rigid) = (3.0, 3.0);
    soft.project(&mut pos_soft, &[1.0; 4], &mut lam_soft, 0.0);
    rigid.project(&mut pos_rigid, &[1.0; 4], &mut lam_rigid, 0.01);
    assert!(pos_soft.iter().flatten().all(|v| v.is_finite()));
    assert_eq!(pos_soft, pos_rigid);
    assert_eq!(lam_soft, lam_rigid);
}

#[test]
fn pressure_constraint_is_independent_of_the_length_scale() {
    // The volume gradient is of the order edge^2, so an absolute threshold on
    // the denominator would silently disable millimetre sized tetrahedra.
    for scale in [1.0e-4, 1.0e-3, 1.0, 1.0e3] {
        let rest = pressure_tet(1.0).map(|p| pbd2_scale(p, scale));
        let rest_volume = tet_volume(rest[0], rest[1], rest[2], rest[3]);
        let c = XpbdPressureConstraint::new([0, 1, 2, 3], rest_volume, 0.0);
        let mut positions = rest.map(|p| pbd2_scale(p, 0.5));
        assert!(
            c.current_volume(&positions) < 0.2 * rest_volume,
            "scale {scale}: the tetrahedron starts compressed"
        );
        let inv_masses = [1.0; 4];
        let mut lambda = 0.0;
        for _ in 0..40 {
            c.project(&mut positions, &inv_masses, &mut lambda, 0.01);
        }
        let volume = c.current_volume(&positions);
        assert!(
            (volume - rest_volume).abs() < 1e-6 * rest_volume,
            "scale {scale}: volume {volume} must reach {rest_volume}"
        );
    }
}

#[test]
fn pressure_constraint_ignores_a_collapsed_tetrahedron() {
    // All four vertices coincide: the volume has no gradient to follow.
    let mut positions = [[0.5, 0.5, 0.5]; 4];
    let c = XpbdPressureConstraint::new([0, 1, 2, 3], 0.1, 0.0);
    let mut lambda = 0.0;
    c.project(&mut positions, &[1.0; 4], &mut lambda, 0.01);
    assert_eq!(positions, [[0.5, 0.5, 0.5]; 4]);
    assert_eq!(lambda, 0.0);
}

// ── PbdSimulation volume constraints ─────────────────────────────────────────

fn unit_tet_simulation() -> PbdSimulation {
    let mut sim = PbdSimulation::new();
    sim.add_particle(Vec3::new(0.0, 0.0, 0.0), 1.0);
    sim.add_particle(Vec3::new(1.0, 0.0, 0.0), 1.0);
    sim.add_particle(Vec3::new(0.0, 1.0, 0.0), 1.0);
    sim.add_particle(Vec3::new(0.0, 0.0, 1.0), 1.0);
    sim
}

fn sim_tet_volume(sim: &PbdSimulation) -> f64 {
    PbdVolumeConstraint::compute_volume(&sim.particles, &[0, 1, 2, 3])
}

#[test]
fn simulation_volume_constraint_is_registered() {
    let mut sim = unit_tet_simulation();
    assert_eq!(sim.volume_constraint_count(), 0);
    sim.add_volume_constraint(0, 1, 2, 3, 0.8)
        .expect("valid tetrahedron");
    assert_eq!(sim.volume_constraint_count(), 1);
    let c = &sim.volume_constraints[0];
    assert_eq!(c.indices, [0, 1, 2, 3]);
    assert!((c.rest_volume - 1.0 / 6.0).abs() < TIGHT);
    assert!((c.stiffness - 0.8).abs() < TIGHT);
}

#[test]
fn simulation_step_restores_compressed_tetrahedron_volume() {
    let mut sim = unit_tet_simulation();
    sim.add_volume_constraint(0, 1, 2, 3, 1.0)
        .expect("valid tetrahedron");
    let rest_volume = sim.volume_constraints[0].rest_volume;
    for p in &mut sim.particles {
        p.position *= 0.5;
        p.prev_position = p.position;
    }
    let compressed = sim_tet_volume(&sim);
    assert!(
        (compressed - rest_volume).abs() > 0.8 * rest_volume,
        "half-size tetrahedron has one eighth of the rest volume"
    );
    sim.step(1.0 / 60.0, Vec3::zeros(), 30);
    let volume = sim_tet_volume(&sim);
    assert!(
        (volume - rest_volume).abs() / rest_volume < 1e-2,
        "volume {volume} must be restored to {rest_volume}"
    );
}

#[test]
fn simulation_without_volume_constraint_keeps_compressed_volume() {
    let mut sim = unit_tet_simulation();
    let rest_volume = sim_tet_volume(&sim);
    for p in &mut sim.particles {
        p.position *= 0.5;
        p.prev_position = p.position;
    }
    sim.step(1.0 / 60.0, Vec3::zeros(), 30);
    let volume = sim_tet_volume(&sim);
    assert!(
        (volume - rest_volume / 8.0).abs() < 1e-12,
        "without a volume constraint nothing may change the volume"
    );
}

#[test]
fn simulation_volume_constraint_conserves_centroid() {
    let mut sim = unit_tet_simulation();
    sim.add_volume_constraint(0, 1, 2, 3, 1.0)
        .expect("valid tetrahedron");
    for p in &mut sim.particles {
        p.position *= 0.5;
        p.prev_position = p.position;
    }
    let centroid_before: Vec3 = sim.particles.iter().map(|p| p.position).sum::<Vec3>() / 4.0;
    sim.step(1.0 / 60.0, Vec3::zeros(), 30);
    let centroid_after: Vec3 = sim.particles.iter().map(|p| p.position).sum::<Vec3>() / 4.0;
    assert!((centroid_after - centroid_before).norm() < 1e-12);
}

#[test]
fn simulation_volume_constraint_rejects_out_of_range_index() {
    let mut sim = unit_tet_simulation();
    let err = sim
        .add_volume_constraint(0, 1, 2, 4, 1.0)
        .expect_err("particle 4 does not exist");
    assert!(
        matches!(err, crate::Error::ConfigurationError { value, .. } if (value - 4.0).abs() < TIGHT),
        "unexpected error {err:?}"
    );
    assert_eq!(sim.volume_constraint_count(), 0);
}

#[test]
fn simulation_volume_constraint_rejects_repeated_index() {
    let mut sim = unit_tet_simulation();
    let err = sim
        .add_volume_constraint(0, 1, 2, 1, 1.0)
        .expect_err("particle 1 is repeated");
    assert!(
        matches!(err, crate::Error::DegenerateConstraint { .. }),
        "unexpected error {err:?}"
    );
    assert_eq!(sim.volume_constraint_count(), 0);
}

// ── PbdSubstepSolver ─────────────────────────────────────────────────────────

#[test]
fn substep_solver_integrates_gravity_into_velocity_and_position() {
    let g = -9.81;
    let n = 10_usize;
    let solver = PbdSubstepSolver::new(n, 1.0).with_gravity([0.0, g, 0.0]);
    let mut positions = [[0.0, 10.0, 0.0]];
    let mut velocities = [[0.0, 0.0, 0.0]];
    solver.solve(&mut positions, &mut velocities, &[2.0], &mut |_| {});
    // Semi-implicit Euler: y_N = y0 + g h^2 N (N + 1) / 2 with h = 1 / N.
    let h = 1.0 / n as f64;
    let expected_y = 10.0 + g * h * h * (n * (n + 1)) as f64 / 2.0;
    assert!(
        (positions[0][1] - expected_y).abs() < 1e-12,
        "y={}",
        positions[0][1]
    );
    assert!(
        (velocities[0][1] - g).abs() < 1e-12,
        "vy={}",
        velocities[0][1]
    );
    assert_eq!(positions[0][0], 0.0);
    assert_eq!(velocities[0][0], 0.0);
}

#[test]
fn substep_solver_feeds_constraint_corrections_back_into_velocity() {
    // A particle falls onto the ground plane y = 0: the projection removes the
    // penetration, and the velocity update must remove the impact velocity.
    let solver = PbdSubstepSolver::new(8, 1.0 / 60.0).with_gravity([0.0, -9.81, 0.0]);
    let mut positions = [[0.0, 1.0, 0.0]];
    let mut velocities = [[0.0, 0.0, 0.0]];
    let mut ground = |p: &mut [[f64; 3]]| {
        for q in p.iter_mut() {
            q[1] = q[1].max(0.0);
        }
    };
    for _ in 0..120 {
        solver.solve(&mut positions, &mut velocities, &[1.0], &mut ground);
    }
    assert!(positions[0][1].abs() < 1e-12, "y={}", positions[0][1]);
    assert!(
        velocities[0][1].abs() < 1e-9,
        "resting particle keeps vy={}",
        velocities[0][1]
    );
}

#[test]
fn substep_solver_keeps_pendulum_on_its_circle() {
    // Anchor (mass 0) at the origin, bob released from rest at (1, 0, 0).
    let solver = PbdSubstepSolver::new(32, 1.0 / 60.0).with_gravity([0.0, -9.81, 0.0]);
    let masses = [0.0, 1.0];
    let mut positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let mut velocities = [[0.0; 3]; 2];
    let mut rod = |p: &mut [[f64; 3]]| {
        let anchor = p[0];
        let d = pbd2_sub(p[1], anchor);
        let len = pbd2_norm(d);
        if len > 1e-12 {
            p[1] = pbd2_add(anchor, pbd2_scale(d, 1.0 / len));
        }
    };
    let mut max_speed = 0.0_f64;
    for frame in 0..100 {
        solver.solve(&mut positions, &mut velocities, &masses, &mut rod);
        let radial = pbd2_sub(positions[1], positions[0]);
        let radius = pbd2_norm(radial);
        assert!(
            (radius - 1.0).abs() < 1e-3,
            "frame {frame}: radius drifted to {radius}"
        );
        let radial_velocity = pbd2_dot(velocities[1], pbd2_scale(radial, 1.0 / radius));
        assert!(
            radial_velocity.abs() < 1e-2,
            "frame {frame}: velocity must be tangential, radial component {radial_velocity}"
        );
        max_speed = max_speed.max(pbd2_norm(velocities[1]));
        assert_eq!(positions[0], [0.0, 0.0, 0.0], "anchor must not move");
        assert_eq!(velocities[0], [0.0; 3], "anchor velocity must stay zero");
    }
    // A bob released at rest from the horizontal reaches sqrt(2 g L) at the bottom.
    let bottom_speed = (2.0_f64 * 9.81).sqrt();
    assert!(
        (max_speed - bottom_speed).abs() < 0.05 * bottom_speed,
        "max speed {max_speed} should be close to {bottom_speed}"
    );
}

#[test]
fn substep_solver_never_moves_massless_particles() {
    let solver = PbdSubstepSolver::new(4, 0.1).with_gravity([0.0, -9.81, 0.0]);
    let mut positions = [[1.0, 2.0, 3.0]];
    let mut velocities = [[4.0, 5.0, 6.0]];
    solver.solve(&mut positions, &mut velocities, &[0.0], &mut |_| {});
    assert_eq!(positions[0], [1.0, 2.0, 3.0]);
    assert_eq!(velocities[0], [4.0, 5.0, 6.0]);
}

#[test]
fn substep_solver_tolerates_zero_substeps_and_short_slices() {
    let solver = PbdSubstepSolver::new(0, 0.1);
    let mut positions = [[0.0; 3], [1.0, 0.0, 0.0]];
    let mut velocities = [[1.0, 0.0, 0.0]];
    solver.solve(&mut positions, &mut velocities, &[1.0, 1.0], &mut |_| {});
    // Zero sub-steps behave as one: only the common prefix (particle 0) moves.
    assert!((positions[0][0] - 0.1).abs() < 1e-12);
    assert_eq!(positions[1], [1.0, 0.0, 0.0]);
    assert!((velocities[0][0] - 1.0).abs() < 1e-12);
}
