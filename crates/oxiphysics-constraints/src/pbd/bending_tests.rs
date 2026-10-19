//! Regression tests for [`XpbdBendingConstraint`]: exact angle gradient,
//! momentum conservation, convergence to the rest angle and the degenerate
//! (straight / folded) configurations.

use super::functions::{pbd2_add, pbd2_cross, pbd2_dot, pbd2_norm, pbd2_scale, pbd2_sub};
use super::types::XpbdBendingConstraint;
use std::f64::consts::{FRAC_PI_2, PI};

fn constraint(rest: f64, compliance: f64) -> XpbdBendingConstraint {
    XpbdBendingConstraint::new(0, 1, 2, rest, compliance)
}

/// Central-difference gradient of the bending angle.
fn numeric_gradient(c: &XpbdBendingConstraint, positions: &[[f64; 3]]) -> [[f64; 3]; 3] {
    let h = 1e-6;
    let mut grad = [[0.0; 3]; 3];
    for (k, &particle) in [c.p1, c.p2, c.p3].iter().enumerate() {
        for axis in 0..3 {
            let mut plus = positions.to_vec();
            let mut minus = positions.to_vec();
            plus[particle][axis] += h;
            minus[particle][axis] -= h;
            grad[k][axis] = (c.compute_angle(&plus) - c.compute_angle(&minus)) / (2.0 * h);
        }
    }
    grad
}

/// Hinge poses: right angle, generic 3-D skew, nearly straight, nearly folded.
fn poses() -> Vec<Vec<[f64; 3]>> {
    vec![
        vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[1.3, 0.2, -0.4], [0.1, -0.3, 0.2], [-0.5, 0.9, 0.7]],
        vec![[2.0, 0.0, 0.0], [0.0, 0.1, 0.0], [-1.5, 0.0, 0.3]],
        vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.9, 0.2, 0.1]],
        vec![[-0.7, 1.1, 0.4], [0.2, 0.2, 0.2], [0.9, -0.8, 1.3]],
    ]
}

/// Linear momentum and angular momentum (about the origin) of the displacement.
fn displacement_momenta(
    before: &[[f64; 3]],
    after: &[[f64; 3]],
    masses: &[f64],
) -> ([f64; 3], [f64; 3]) {
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for i in 0..before.len() {
        let dp = pbd2_sub(after[i], before[i]);
        linear = pbd2_add(linear, pbd2_scale(dp, masses[i]));
        angular = pbd2_add(angular, pbd2_scale(pbd2_cross(before[i], dp), masses[i]));
    }
    (linear, angular)
}

fn max_abs(v: [f64; 3]) -> f64 {
    v.iter().fold(0.0_f64, |m, x| m.max(x.abs()))
}

#[test]
fn angle_gradients_match_central_differences() {
    let c = constraint(PI, 0.0);
    for (n, pose) in poses().iter().enumerate() {
        let analytic = c.angle_gradients(pose).expect("well conditioned pose");
        let numeric = numeric_gradient(&c, pose);
        for k in 0..3 {
            let err = max_abs(pbd2_sub(analytic[k], numeric[k]));
            assert!(
                err < 1e-7,
                "pose {n}, particle {k}: analytic {:?} vs numeric {:?}",
                analytic[k],
                numeric[k]
            );
        }
    }
}

#[test]
fn angle_gradients_are_invariant_and_have_edge_length_norms() {
    let c = constraint(PI, 0.0);
    for (n, pose) in poses().iter().enumerate() {
        let g = c.angle_gradients(pose).expect("well conditioned pose");
        // Translation invariance: the gradients sum to zero.
        let sum = pbd2_add(pbd2_add(g[0], g[1]), g[2]);
        assert!(max_abs(sum) < 1e-12, "pose {n}: sum {sum:?}");
        // Rotation invariance: the gradients are torque free.
        let torque = (0..3).fold([0.0; 3], |acc, k| pbd2_add(acc, pbd2_cross(pose[k], g[k])));
        assert!(max_abs(torque) < 1e-12, "pose {n}: torque {torque:?}");
        // |dtheta/dp1| = 1/|p1 - p2| and |dtheta/dp3| = 1/|p3 - p2|.
        let len1 = pbd2_norm(pbd2_sub(pose[0], pose[1]));
        let len3 = pbd2_norm(pbd2_sub(pose[2], pose[1]));
        assert!((pbd2_norm(g[0]) * len1 - 1.0).abs() < 1e-12, "pose {n}");
        assert!((pbd2_norm(g[2]) * len3 - 1.0).abs() < 1e-12, "pose {n}");
        // Each end gradient is perpendicular to its own edge.
        assert!(pbd2_dot(g[0], pbd2_sub(pose[0], pose[1])).abs() < 1e-12);
        assert!(pbd2_dot(g[2], pbd2_sub(pose[2], pose[1])).abs() < 1e-12);
    }
}

#[test]
fn angle_is_resolved_next_to_pi_and_zero() {
    let c = constraint(PI, 0.0);
    let delta = 1e-9;
    // Hinge lifted by `delta` above the line between the end particles:
    // theta = pi - 2 atan(delta).  acos would round this to exactly pi.
    let near_straight = [[-1.0, 0.0, 0.0], [0.0, delta, 0.0], [1.0, 0.0, 0.0]];
    let deficit = PI - c.compute_angle(&near_straight);
    assert!(
        (deficit - 2.0 * delta).abs() < 1e-12,
        "deficit {deficit} should be {}",
        2.0 * delta
    );
    // Same next to zero: both end particles on one side of the hinge.
    let near_folded = [[1.0, delta, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let opening = c.compute_angle(&near_folded);
    assert!((opening - delta).abs() < 1e-12, "opening {opening}");
    // Exactly straight and exactly folded.
    assert_eq!(
        c.compute_angle(&[[-1.0, 0.0, 0.0], [0.0; 3], [1.0, 0.0, 0.0]]),
        PI
    );
    assert_eq!(
        c.compute_angle(&[[1.0, 0.0, 0.0], [0.0; 3], [2.0, 0.0, 0.0]]),
        0.0
    );
    // Right angle.
    let right = c.compute_angle(&poses()[0]);
    assert!((right - FRAC_PI_2).abs() < 1e-15, "right angle {right}");
}

#[test]
fn ninety_degree_bend_straightens_within_twenty_iterations() {
    let c = constraint(PI, 0.0);
    let masses = [1.0, 1.0, 1.0];
    let mut positions = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let centroid_before = pbd2_scale(
        positions.iter().fold([0.0; 3], |acc, p| pbd2_add(acc, *p)),
        1.0 / 3.0,
    );
    let mut lambda = 0.0;
    let mut straightened_at = None;
    for iteration in 0..20 {
        let before = positions.clone();
        c.project(&mut positions, &masses, &mut lambda, 0.01);
        let (linear, angular) = displacement_momenta(&before, &positions, &masses);
        assert!(max_abs(linear) < 1e-14, "linear momentum {linear:?}");
        assert!(max_abs(angular) < 1e-14, "angular momentum {angular:?}");
        if (PI - c.compute_angle(&positions)).abs() < 1e-9 && straightened_at.is_none() {
            straightened_at = Some(iteration);
        }
    }
    let n = straightened_at.expect("the chain must straighten within 20 iterations");
    assert!(n < 10, "straightened only after {n} iterations");
    // The three particles are now collinear ...
    let e1 = pbd2_sub(positions[0], positions[1]);
    let e2 = pbd2_sub(positions[2], positions[1]);
    assert!(pbd2_norm(pbd2_cross(e1, e2)) < 1e-8 * pbd2_norm(e1) * pbd2_norm(e2));
    // ... and the centre of mass has not moved (momentum conservation).
    let centroid_after = pbd2_scale(
        positions.iter().fold([0.0; 3], |acc, p| pbd2_add(acc, *p)),
        1.0 / 3.0,
    );
    assert!(max_abs(pbd2_sub(centroid_after, centroid_before)) < 1e-13);
}

#[test]
fn unequal_masses_conserve_linear_and_angular_momentum() {
    let masses = [1.0, 3.0, 2.0];
    for rest in [PI, 2.0, 0.7] {
        let c = constraint(rest, 0.0);
        let mut positions = vec![[1.0, 0.2, -0.3], [0.0, 0.0, 0.0], [-0.4, 1.0, 0.5]];
        let mut lambda = 0.0;
        for _ in 0..8 {
            let before = positions.clone();
            c.project(&mut positions, &masses, &mut lambda, 0.01);
            let (linear, angular) = displacement_momenta(&before, &positions, &masses);
            assert!(max_abs(linear) < 1e-13, "rest {rest}: linear {linear:?}");
            assert!(max_abs(angular) < 1e-13, "rest {rest}: angular {angular:?}");
        }
        let theta = c.compute_angle(&positions);
        assert!((theta - rest).abs() < 1e-9, "rest {rest}: theta {theta}");
    }
}

#[test]
fn pinned_end_particles_stay_put_while_the_hinge_moves() {
    let masses = [0.0, 1.0, 0.0];
    let ends = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    // Straighten: the hinge is pulled onto the line through the pinned ends.
    let c = constraint(PI, 0.0);
    let mut positions = vec![ends[0], [0.0, 0.5, 0.0], ends[1]];
    let mut lambda = 0.0;
    for _ in 0..30 {
        c.project(&mut positions, &masses, &mut lambda, 0.01);
    }
    assert_eq!(positions[0], ends[0]);
    assert_eq!(positions[2], ends[1]);
    assert!(positions[1][1].abs() < 1e-9, "hinge at {:?}", positions[1]);
    // A right angle between pinned ends lies on the Thales circle of radius 1.
    let c = constraint(FRAC_PI_2, 0.0);
    let mut positions = vec![ends[0], [0.0, 0.5, 0.0], ends[1]];
    let mut lambda = 0.0;
    for _ in 0..30 {
        c.project(&mut positions, &masses, &mut lambda, 0.01);
    }
    assert!(
        (positions[1][1] - 1.0).abs() < 1e-9,
        "hinge at {:?}",
        positions[1]
    );
    assert!(positions[1][0].abs() < 1e-9 && positions[1][2].abs() < 1e-9);
}

#[test]
fn fully_pinned_or_massless_constraints_do_nothing() {
    let c = constraint(PI, 0.0);
    let start = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    for masses in [[0.0; 3], [-1.0, 0.0, -2.0], [f64::INFINITY; 3]] {
        let mut positions = start.clone();
        let mut lambda = 0.0;
        c.project(&mut positions, &masses, &mut lambda, 0.01);
        assert_eq!(positions, start);
        assert_eq!(lambda, 0.0);
    }
}

#[test]
fn compliance_follows_the_xpbd_update() {
    let masses = [1.0, 2.0, 1.0];
    let start = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let (compliance, dt) = (0.05, 0.1);
    let c = constraint(PI, compliance);
    // Expected multiplier from the closed-form XPBD update at lambda = 0.
    let g = c.angle_gradients(&start).expect("well conditioned");
    let denom: f64 = (0..3).map(|k| pbd2_dot(g[k], g[k]) / masses[k]).sum();
    let alpha_tilde = compliance / (dt * dt);
    let theta_c = c.compute_angle(&start) - PI;
    let expected_lambda = -theta_c / (denom + alpha_tilde);
    let mut soft = start.clone();
    let mut lambda = 0.0;
    c.project(&mut soft, &masses, &mut lambda, dt);
    assert!(
        (lambda - expected_lambda).abs() < 1e-14,
        "{lambda} vs {expected_lambda}"
    );
    // The compliant correction is smaller than the rigid one.
    let mut rigid = start.clone();
    let mut rigid_lambda = 0.0;
    constraint(PI, 0.0).project(&mut rigid, &masses, &mut rigid_lambda, dt);
    assert!(lambda.abs() < rigid_lambda.abs());
    let moved = |p: &[[f64; 3]]| pbd2_norm(pbd2_sub(p[1], start[1]));
    assert!(moved(&soft) < moved(&rigid));
}

#[test]
fn compliant_constraint_converges_to_the_xpbd_fixed_point() {
    // At the fixed point of the accumulated update, C + alpha_tilde * lambda = 0.
    let masses = [1.0, 1.0, 1.0];
    let (compliance, dt) = (1e-3, 0.1);
    let c = constraint(PI, compliance);
    let mut positions = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let mut lambda = 0.0;
    for _ in 0..200 {
        c.project(&mut positions, &masses, &mut lambda, dt);
    }
    let residual = c.compute_angle(&positions) - PI + compliance / (dt * dt) * lambda;
    assert!(residual.abs() < 1e-12, "fixed-point residual {residual}");
    // Soft: the chain stops short of fully straight.
    assert!(PI - c.compute_angle(&positions) > 1e-3);
}

#[test]
fn straight_chain_with_a_bent_rest_angle_bends_about_a_perpendicular_axis() {
    let masses = [1.0, 2.0, 1.0];
    let c = constraint(FRAC_PI_2, 0.0);
    let start = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    assert!(
        c.angle_gradients(&start).is_none(),
        "gradient is undefined here"
    );
    let mut positions = start.clone();
    let mut lambda = 0.0;
    for _ in 0..20 {
        c.project(&mut positions, &masses, &mut lambda, 0.01);
    }
    let theta = c.compute_angle(&positions);
    assert!((theta - FRAC_PI_2).abs() < 1e-9, "theta {theta}");
    let (linear, angular) = displacement_momenta(&start, &positions, &masses);
    assert!(max_abs(linear) < 1e-12 && max_abs(angular) < 1e-12);
    assert!(positions.iter().flatten().all(|v| v.is_finite()));
}

#[test]
fn folded_chain_with_an_open_rest_angle_unfolds() {
    let masses = [1.0, 1.0, 1.0];
    let c = constraint(FRAC_PI_2, 0.0);
    let mut positions = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let mut lambda = 0.0;
    for _ in 0..20 {
        c.project(&mut positions, &masses, &mut lambda, 0.01);
    }
    let theta = c.compute_angle(&positions);
    assert!((theta - FRAC_PI_2).abs() < 1e-9, "theta {theta}");
}

#[test]
fn slightly_open_hinge_closes_to_a_zero_rest_angle() {
    let masses = [1.0, 1.0, 1.0];
    let c = constraint(0.0, 0.0);
    let mut positions = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.5, 0.0]];
    let mut lambda = 0.0;
    for _ in 0..30 {
        c.project(&mut positions, &masses, &mut lambda, 0.01);
    }
    assert!(c.compute_angle(&positions) < 1e-9);
}

#[test]
fn satisfied_degenerate_configurations_are_left_alone() {
    let masses = [1.0, 1.0, 1.0];
    // Straight chain at rest angle pi: C = 0 and the bend plane is undefined.
    let c = constraint(PI, 0.0);
    let straight = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    let mut positions = straight.clone();
    let mut lambda = 0.0;
    c.project(&mut positions, &masses, &mut lambda, 0.01);
    assert_eq!(positions, straight);
    assert_eq!(lambda, 0.0);
    // Coincident particles have no edge direction at all.
    let coincident = vec![[0.0; 3], [0.0; 3], [0.0, 1.0, 0.0]];
    let mut positions = coincident.clone();
    c.project(&mut positions, &masses, &mut lambda, 0.01);
    assert_eq!(positions, coincident);
    assert_eq!(c.compute_angle(&coincident), PI);
    assert!(c.angle_gradients(&coincident).is_none());
}
