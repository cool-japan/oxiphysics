//! Regression tests for the four-bar linkage kinematics.
//!
//! The ground link `l1` lies on the +x axis, the crank pivot is the origin and
//! the follower pivot is `(l1, 0)`, so a solved configuration must satisfy the
//! loop-closure condition `|(l2·cosθ2, l2·sinθ2) − (l1 + l4·cosθ4, l4·sinθ4)| = l3`.

use std::f64::consts::PI;

use super::functions::four_bar_angles;
use super::types::{FourBarAnalysis, FourBarLinkage};

/// Distance between the crank tip and the follower tip minus the coupler length.
fn loop_closure_error(links: [f64; 4], theta2: f64, theta4: f64) -> f64 {
    let [l1, l2, l3, l4] = links;
    let crank_tip = (l2 * theta2.cos(), l2 * theta2.sin());
    let follower_tip = (l1 + l4 * theta4.cos(), l4 * theta4.sin());
    let gap = (crank_tip.0 - follower_tip.0).hypot(crank_tip.1 - follower_tip.1);
    gap - l3
}

fn solve(links: [f64; 4], theta2: f64) -> Option<f64> {
    FourBarLinkage::new(links[0], links[1], links[2], links[3]).follower_angle(theta2)
}

#[test]
fn follower_angle_closes_the_loop_when_crank_and_follower_differ() {
    // l2 != l4 is exactly the case in which swapping K1 and K2 matters.
    let links = [4.0, 1.5, 4.5, 3.0];
    for theta2 in [0.3, 1.0, 2.0] {
        let theta4 = solve(links, theta2).expect("linkage assembles at this crank angle");
        let error = loop_closure_error(links, theta2, theta4);
        assert!(
            error.abs() < 1e-9,
            "loop closure violated at theta2={theta2}: error={error}"
        );
    }
}

#[test]
fn follower_angle_closes_the_loop_for_second_geometry() {
    // With K1/K2 swapped this geometry returned None at 0.3 and 2.0 and a
    // coupler length of 5.07 at 1.0 instead of 3.5.
    let links = [4.0, 1.5, 3.5, 3.0];
    for theta2 in [0.3, 1.0, 2.0] {
        let theta4 = solve(links, theta2).expect("linkage assembles at this crank angle");
        let error = loop_closure_error(links, theta2, theta4);
        assert!(
            error.abs() < 1e-9,
            "loop closure violated at theta2={theta2}: error={error}"
        );
    }
}

#[test]
fn follower_angle_closes_the_loop_over_a_full_crank_revolution() {
    // Crank-rocker: the crank can turn through every angle.
    let links = [4.0, 2.0, 4.0, 3.0];
    for step in 0..72 {
        let theta2 = -PI + 2.0 * PI * f64::from(step) / 72.0;
        let theta4 = solve(links, theta2).expect("crank-rocker assembles at every crank angle");
        let error = loop_closure_error(links, theta2, theta4);
        assert!(
            error.abs() < 1e-9,
            "loop closure violated at theta2={theta2}: error={error}"
        );
        assert!(
            theta4 > -PI - 1e-12 && theta4 <= PI + 1e-12,
            "theta4={theta4} is outside the principal range"
        );
    }
}

#[test]
fn follower_angle_stays_accurate_when_the_quadratic_coefficient_vanishes() {
    // l1 = l4 and l2 = l3 make Freudenstein's A coefficient identically zero,
    // so the quadratic in tan(theta4 / 2) degenerates to a linear equation.
    let links = [3.0, 2.0, 2.0, 3.0];
    for theta2 in [0.3, 1.0, 2.0, PI] {
        let theta4 = solve(links, theta2).expect("kite assembles at this crank angle");
        assert!(theta4.is_finite());
        let error = loop_closure_error(links, theta2, theta4);
        assert!(
            error.abs() < 1e-9,
            "loop closure violated at theta2={theta2} with A = 0: error={error}"
        );
    }
}

#[test]
fn follower_angle_is_none_when_the_linkage_cannot_assemble() {
    // The crank tip is always farther than l3 + l4 = 2 from the follower pivot.
    for step in 0..36 {
        let theta2 = 2.0 * PI * f64::from(step) / 36.0;
        assert!(solve([4.0, 1.5, 1.0, 1.0], theta2).is_none());
    }
    // A zero-length crank or follower has no Freudenstein constants.
    assert!(solve([4.0, 0.0, 3.0, 3.0], 0.5).is_none());
    assert!(solve([4.0, 1.5, 3.0, 0.0], 0.5).is_none());
}

#[test]
fn four_bar_angles_return_a_consistent_coupler_direction() {
    let links = [4.0, 2.0, 4.0, 3.0];
    for theta2 in [0.3, 1.0, 2.0, 3.0] {
        let (theta3, theta4) = four_bar_angles(links[0], links[1], links[2], links[3], theta2)
            .expect("crank-rocker assembles");
        let crank_tip = (links[1] * theta2.cos(), links[1] * theta2.sin());
        let follower_tip = (links[0] + links[3] * theta4.cos(), links[3] * theta4.sin());
        // crank tip + coupler = follower tip
        let coupler_end = (
            crank_tip.0 + links[2] * theta3.cos(),
            crank_tip.1 + links[2] * theta3.sin(),
        );
        let miss = (coupler_end.0 - follower_tip.0).hypot(coupler_end.1 - follower_tip.1);
        assert!(
            miss < 1e-9,
            "coupler direction does not reach the follower tip at theta2={theta2}: miss={miss}"
        );
    }
}

/// Smallest signed difference between two angles.
fn angle_difference(a: f64, b: f64) -> f64 {
    (a - b + PI).rem_euclid(2.0 * PI) - PI
}

#[test]
fn analysis_angular_velocities_match_finite_differences_of_the_angles() {
    let links = [4.0, 2.0, 4.0, 3.0];
    let omega2 = 1.0;
    let step = 1e-6;
    for theta2 in [0.3, 1.0, 2.0] {
        let analysis =
            FourBarAnalysis::analyse(links[0], links[1], links[2], links[3], theta2, omega2)
                .expect("crank-rocker assembles");
        let (theta3_plus, theta4_plus) =
            four_bar_angles(links[0], links[1], links[2], links[3], theta2 + step)
                .expect("assembles");
        let (theta3_minus, theta4_minus) =
            four_bar_angles(links[0], links[1], links[2], links[3], theta2 - step)
                .expect("assembles");
        let omega3_reference = angle_difference(theta3_plus, theta3_minus) / (2.0 * step) * omega2;
        let omega4_reference = angle_difference(theta4_plus, theta4_minus) / (2.0 * step) * omega2;
        assert!(
            (analysis.omega3 - omega3_reference).abs() < 1e-5,
            "omega3 at theta2={theta2}: analytic {} vs finite difference {omega3_reference}",
            analysis.omega3
        );
        assert!(
            (analysis.omega4 - omega4_reference).abs() < 1e-5,
            "omega4 at theta2={theta2}: analytic {} vs finite difference {omega4_reference}",
            analysis.omega4
        );
    }
}
