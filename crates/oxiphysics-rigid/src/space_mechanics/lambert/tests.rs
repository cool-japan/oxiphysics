// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the Izzo Lambert solver.
//!
//! Solutions are checked against numerical propagation with the module's own
//! RK4 integrator (two-body gravity only), against a textbook example and
//! against the analytic properties of the time-of-flight function.

use super::*;
use crate::rng_util::uniform;
use crate::space_mechanics::{
    MU_EARTH, OrbitState, OrbitalElements, R_EARTH, elements_to_cartesian, vec3_dot,
};
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::f64::consts::PI;

/// Propagate a two-body state with RK4 (no J2, no drag) for `tof` seconds.
fn rk4_propagate(r: [f64; 3], v: [f64; 3], tof: f64, mu: f64) -> ([f64; 3], [f64; 3]) {
    let steps = (tof.ceil() as usize).max(4_000);
    let dt = tof / steps as f64;
    let mut state = OrbitState::new(r, v);
    for _ in 0..steps {
        state.rk4_step(dt, mu, 0.0, R_EARTH, 0.0, 0.0);
    }
    (state.pos, state.vel)
}

/// Solve, propagate the departure state and require arrival at `r2`.
fn assert_reaches_target(
    label: &str,
    r1: [f64; 3],
    r2: [f64; 3],
    tof: f64,
    prograde: bool,
) -> LambertSolution {
    let solution = lambert_solver(r1, r2, tof, MU_EARTH, prograde)
        .unwrap_or_else(|| panic!("{label}: the solver found no transfer"));
    let (r_end, v_end) = rk4_propagate(r1, solution.v1, tof, MU_EARTH);
    let miss = vec3_norm(vec3_sub(r_end, r2));
    assert!(miss < 1.0, "{label}: arrival misses r2 by {miss} m");
    let velocity_error = vec3_norm(vec3_sub(v_end, solution.v2));
    assert!(
        velocity_error < 1e-3,
        "{label}: arrival velocity is off by {velocity_error} m/s"
    );
    // The conic conserves angular momentum and energy between the two ends.
    let h1 = vec3_cross(r1, solution.v1);
    let h2 = vec3_cross(r2, solution.v2);
    assert!(
        vec3_norm(vec3_sub(h1, h2)) <= 1e-8 * vec3_norm(h1),
        "{label}: angular momentum differs between the ends"
    );
    let energy = |r: [f64; 3], v: [f64; 3]| 0.5 * vec3_dot(v, v) - MU_EARTH / vec3_norm(r);
    let (e1, e2) = (energy(r1, solution.v1), energy(r2, solution.v2));
    // The specific energy of a parabolic transfer is zero, so the tolerance is
    // scaled by the potential energy at the departure point instead of by |e1|.
    let energy_scale = MU_EARTH / vec3_norm(r1);
    assert!(
        (e1 - e2).abs() <= 1e-8 * energy_scale,
        "{label}: energy differs between the ends ({e1} vs {e2})"
    );
    solution
}

// ── textbook example ──────────────────────────────────────────────────────

#[test]
fn reproduces_vallado_example_7_5() {
    // Vallado, Fundamentals of Astrodynamics and Applications, Example 7-5
    // (km, km/s): r1 = 15945.34 I, r2 = 12214.83899 I + 10249.46731 J,
    // 76 minutes, short way.
    let mu = 398_600.441_8;
    let solution = lambert_solver(
        [15_945.34, 0.0, 0.0],
        [12_214.838_99, 10_249.467_31, 0.0],
        76.0 * 60.0,
        mu,
        true,
    )
    .expect("the example is a regular elliptic transfer");
    let expected_v1 = [2.058_913, 2.915_965, 0.0];
    let expected_v2 = [-3.451_565, 0.910_315, 0.0];
    for axis in 0..3 {
        assert!(
            (solution.v1[axis] - expected_v1[axis]).abs() < 2e-6,
            "v1 = {:?}",
            solution.v1
        );
        assert!(
            (solution.v2[axis] - expected_v2[axis]).abs() < 2e-6,
            "v2 = {:?}",
            solution.v2
        );
    }
    assert!(solution.semi_major_axis > 0.0 && (solution.semi_major_axis - 10_699.57).abs() < 0.05);
    assert!(
        solution.iterations <= 6,
        "iterations: {}",
        solution.iterations
    );
}

// ── propagation oracle ────────────────────────────────────────────────────

#[test]
fn solutions_propagate_to_the_target_for_a_range_of_geometries() {
    let leo = R_EARTH + 400_000.0;
    let angle = 250.0_f64.to_radians();
    let leo_far = [7.0e6 * angle.cos(), 7.0e6 * angle.sin(), 0.0];
    let arc = 0.01_f64;

    // Parabolic transfer time from Euler's formula for the chord/semi-perimeter.
    let (p1, p2) = ([7.0e6, 0.0, 0.0], [0.0, 1.2e7, 5.0e6]);
    let chord = vec3_norm(vec3_sub(p2, p1));
    let semi = 0.5 * (vec3_norm(p1) + vec3_norm(p2) + chord);
    let parabolic_tof =
        (2.0_f64.sqrt() / 3.0 / MU_EARTH.sqrt()) * (semi.powf(1.5) - (semi - chord).powf(1.5));

    struct Case {
        label: &'static str,
        r1: [f64; 3],
        r2: [f64; 3],
        tof: f64,
        prograde: bool,
    }
    let cases = [
        Case {
            label: "leo quarter orbit",
            r1: [leo, 0.0, 0.0],
            r2: [0.0, leo, 0.0],
            tof: 1_500.0,
            prograde: true,
        },
        Case {
            label: "high transfer",
            r1: [7.0e6, 0.0, 0.0],
            r2: [-3.0e7, 2.0e7, 0.0],
            tof: 18_000.0,
            prograde: true,
        },
        Case {
            label: "250 degree long way",
            r1: [7.0e6, 0.0, 0.0],
            r2: leo_far,
            tof: 4_000.0,
            prograde: true,
        },
        Case {
            label: "retrograde 270 degrees",
            r1: [7.0e6, 0.0, 0.0],
            r2: [0.0, 7.0e6, 0.0],
            tof: 4_000.0,
            prograde: false,
        },
        Case {
            label: "retrograde 270 degrees, unequal radii",
            r1: [7.0e6, 0.0, 0.0],
            r2: [0.0, 7.5e6, 0.0],
            tof: 3_000.0,
            prograde: false,
        },
        Case {
            label: "prograde 90 degrees, unequal radii",
            r1: [7.0e6, 0.0, 0.0],
            r2: [0.0, 7.5e6, 0.0],
            tof: 3_000.0,
            prograde: true,
        },
        Case {
            label: "inclined plane",
            r1: [6.8e6, 5.0e5, 3.0e5],
            r2: [-2.0e6, 7.0e6, 2.5e6],
            tof: 2_500.0,
            prograde: true,
        },
        Case {
            label: "hyperbolic",
            r1: [7.0e6, 0.0, 0.0],
            r2: [2.0e7, 2.5e7, 0.0],
            tof: 1_500.0,
            prograde: true,
        },
        Case {
            label: "parabolic",
            r1: p1,
            r2: p2,
            tof: parabolic_tof,
            prograde: true,
        },
        Case {
            label: "short arc",
            r1: [7.0e6, 0.0, 0.0],
            r2: [7.0e6 * arc.cos(), 7.0e6 * arc.sin(), 0.0],
            tof: 9.3,
            prograde: true,
        },
        Case {
            label: "polar plane",
            r1: [7.0e6, 0.0, 0.0],
            r2: [0.0, 0.0, 7.0e6],
            tof: 1_300.0,
            prograde: true,
        },
        Case {
            label: "fast hyperbolic flyby",
            r1: [7.0e6, 0.0, 0.0],
            r2: [4.0e7, 3.0e7, 1.0e7],
            tof: 900.0,
            prograde: true,
        },
    ];
    for case in &cases {
        let solution = assert_reaches_target(case.label, case.r1, case.r2, case.tof, case.prograde);
        match case.label {
            "hyperbolic" | "fast hyperbolic flyby" => {
                assert!(
                    solution.semi_major_axis < 0.0,
                    "{}: a = {}",
                    case.label,
                    solution.semi_major_axis
                )
            }
            "parabolic" => assert!(
                solution.semi_major_axis.abs() > 1e12,
                "parabolic transfer has a = {}",
                solution.semi_major_axis
            ),
            _ => assert!(
                solution.semi_major_axis > 0.0,
                "{}: a = {}",
                case.label,
                solution.semi_major_axis
            ),
        }
    }
}

#[test]
fn recovers_the_true_velocities_of_random_orbits() {
    let mut rng = StdRng::seed_from_u64(20_150_601);
    let mut checked = 0;
    for case in 0..40 {
        let a = uniform(&mut rng, 7.0e6, 4.2e7);
        let e_max = (1.0 - 6.6e6 / a).clamp(0.0, 0.7);
        let elements = OrbitalElements::new(
            a,
            uniform(&mut rng, 0.0, e_max),
            uniform(&mut rng, 0.0, PI),
            uniform(&mut rng, 0.0, 2.0 * PI),
            uniform(&mut rng, 0.0, 2.0 * PI),
            uniform(&mut rng, 0.0, 2.0 * PI),
        );
        let (r1, v1) = elements_to_cartesian(&elements, MU_EARTH);
        let period = 2.0 * PI * (a.powi(3) / MU_EARTH).sqrt();
        let tof = period * uniform(&mut rng, 0.02, 0.98);
        let (r2, v2) = rk4_propagate(r1, v1, tof, MU_EARTH);
        // Skip transfers whose plane is nearly undefined (angle close to 0 or 180 degrees).
        let sin_angle = vec3_norm(vec3_cross(r1, r2)) / (vec3_norm(r1) * vec3_norm(r2));
        if sin_angle < 1e-3 {
            continue;
        }
        let prograde = vec3_cross(r1, v1)[2] >= 0.0;
        let solution = lambert_solver(r1, r2, tof, MU_EARTH, prograde)
            .unwrap_or_else(|| panic!("case {case}: no transfer found"));
        let e1 = vec3_norm(vec3_sub(solution.v1, v1)) / vec3_norm(v1);
        let e2 = vec3_norm(vec3_sub(solution.v2, v2)) / vec3_norm(v2);
        assert!(
            e1 < 1e-6 && e2 < 1e-6,
            "case {case}: relative velocity errors {e1:e}, {e2:e}"
        );
        assert!(
            solution.iterations <= 12,
            "case {case}: {} iterations",
            solution.iterations
        );
        checked += 1;
    }
    assert!(checked >= 35, "only {checked} of 40 orbits were usable");
}

#[test]
fn reversing_the_transfer_negates_and_swaps_the_velocities() {
    let cases = [
        ([7.0e6, 0.0, 0.0], [-1.0e6, 8.0e6, 1.0e6], 3_000.0),
        ([6.8e6, 5.0e5, 3.0e5], [-2.0e6, 7.0e6, 2.5e6], 2_500.0),
        ([1.0e7, 2.0e6, -1.0e6], [-9.0e6, -3.0e6, 2.0e6], 9_000.0),
    ];
    for (r1, r2, tof) in cases {
        for prograde in [true, false] {
            let forward =
                lambert_solver(r1, r2, tof, MU_EARTH, prograde).expect("forward transfer");
            let backward =
                lambert_solver(r2, r1, tof, MU_EARTH, !prograde).expect("reverse transfer");
            let scale = vec3_norm(forward.v1);
            assert!(vec3_norm(vec3_add(forward.v1, backward.v2)) <= 1e-9 * scale);
            assert!(vec3_norm(vec3_add(forward.v2, backward.v1)) <= 1e-9 * scale);
            assert!(
                (forward.semi_major_axis - backward.semi_major_axis).abs()
                    <= 1e-9 * forward.semi_major_axis.abs()
            );
        }
    }
}

#[test]
fn prograde_and_retrograde_transfers_differ_and_have_opposite_angular_momentum() {
    // Unequal radii (rho != 0) exercise the radial/tangential coupling of the
    // long-way branch, which equal-radius transfers cannot reveal.
    let (r1, r2, tof) = ([7.0e6, 0.0, 0.0], [0.0, 7.5e6, 0.0], 3_000.0);
    let prograde = assert_reaches_target("prograde, 90 degrees", r1, r2, tof, true);
    let retrograde = assert_reaches_target("retrograde, 270 degrees", r1, r2, tof, false);
    assert!(vec3_cross(r1, prograde.v1)[2] > 0.0);
    assert!(vec3_cross(r1, retrograde.v1)[2] < 0.0);
    assert!(vec3_cross(r2, prograde.v2)[2] > 0.0);
    assert!(vec3_cross(r2, retrograde.v2)[2] < 0.0);
    // Two different conics connect the same points in the same time.
    assert!(vec3_norm(vec3_sub(prograde.v1, retrograde.v1)) > 100.0);
    assert!((prograde.semi_major_axis - retrograde.semi_major_axis).abs() > 1.0e3);
}

// ── failure modes ─────────────────────────────────────────────────────────

#[test]
fn invalid_or_degenerate_requests_return_none() {
    let r1 = [7.0e6, 0.0, 0.0];
    let r2 = [0.0, 7.0e6, 0.0];
    let ok = |r1, r2, tof, mu| lambert_solver(r1, r2, tof, mu, true).is_some();
    assert!(ok(r1, r2, 1_500.0, MU_EARTH));
    assert!(!ok(r1, r2, 0.0, MU_EARTH));
    assert!(!ok(r1, r2, -10.0, MU_EARTH));
    assert!(!ok(r1, r2, f64::NAN, MU_EARTH));
    assert!(!ok(r1, r2, f64::INFINITY, MU_EARTH));
    assert!(!ok(r1, r2, 1_500.0, 0.0));
    assert!(!ok(r1, r2, 1_500.0, -MU_EARTH));
    assert!(!ok([0.0; 3], r2, 1_500.0, MU_EARTH));
    assert!(!ok(r1, [0.0; 3], 1_500.0, MU_EARTH));
    assert!(!ok([f64::NAN, 0.0, 0.0], r2, 1_500.0, MU_EARTH));
    // 0 and 180 degree transfers have no defined plane.
    assert!(!ok(r1, [1.4e7, 0.0, 0.0], 1_500.0, MU_EARTH));
    assert!(!ok(r1, [-9.0e6, 0.0, 0.0], 1_500.0, MU_EARTH));
    assert!(!ok(r1, r1, 1_500.0, MU_EARTH));
}

#[test]
fn unreachably_long_flight_times_do_not_report_a_fake_solution() {
    // The previous implementation clamped x to (-1, 1) and always answered.
    let r1 = [7.0e6, 0.0, 0.0];
    let r2 = [0.0, 7.0e6, 0.0];
    assert!(lambert_solver(r1, r2, 1.0e300, MU_EARTH, true).is_none());
    assert!(solve_for_x(0.5, f64::NAN).is_none());
    assert!(solve_for_x(0.5, 1.0e300).is_none());
    assert!(solve_for_x(1.0, 1.0).is_none());
    assert!(solve_for_x(0.5, -1.0).is_none());
}

// ── time-of-flight function ───────────────────────────────────────────────

const LAMBDAS: [f64; 9] = [-0.995, -0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9, 0.995];

#[test]
fn root_finder_converges_over_the_whole_lambda_time_plane() {
    for &lambda in &LAMBDAS {
        for k in 0..60 {
            let t_target = 10.0_f64.powf(-4.0 + 7.0 * f64::from(k) / 59.0);
            let (x, iterations) = solve_for_x(lambda, t_target)
                .unwrap_or_else(|| panic!("no root for lambda {lambda}, T {t_target}"));
            let t = tof_and_derivatives(x, lambda).time;
            assert!(
                (t - t_target).abs() <= 1e-11 * t_target.max(1.0),
                "lambda {lambda}, T {t_target}: residual {}",
                t - t_target
            );
            assert!(
                iterations <= 40,
                "lambda {lambda}, T {t_target}: {iterations} iterations"
            );
        }
    }
}

#[test]
fn time_of_flight_decreases_strictly_and_covers_all_conics() {
    for &lambda in &LAMBDAS {
        let mut previous = f64::INFINITY;
        for k in 0..400 {
            let x = -0.99 + 0.0125 * f64::from(k);
            let t = tof_and_derivatives(x, lambda).time;
            assert!(
                t > 0.0 && t < previous,
                "lambda {lambda}: T not decreasing at x = {x}"
            );
            previous = t;
        }
    }
}

#[test]
fn parabolic_time_of_flight_matches_euler_formula() {
    for &lambda in &LAMBDAS {
        let eval = tof_and_derivatives(1.0, lambda);
        let expected = 2.0 * (1.0 - lambda.powi(3)) / 3.0;
        assert!(
            (eval.time - expected).abs() < 1e-14,
            "lambda {lambda}: {} vs {expected}",
            eval.time
        );
        assert!(eval.d1.is_finite() && eval.d2.is_finite() && eval.d3.is_finite());
        assert!(eval.d1 < 0.0, "T decreases through the parabolic point");
    }
}

#[test]
fn derivatives_match_finite_differences_including_the_parabolic_point() {
    let xs = [
        -0.6, 0.3, 0.6, 0.75, 0.8, 0.9, 0.99, 1.0, 1.01, 1.1, 1.17, 1.19, 1.6, 3.0,
    ];
    for &lambda in &[-0.7, 0.2, 0.9] {
        for &x in &xs {
            let h = 1e-4;
            let at = |x: f64| tof_and_derivatives(x, lambda);
            let centre = at(x);
            let (plus, minus) = (at(x + h), at(x - h));
            let fd1 = (plus.time - minus.time) / (2.0 * h);
            let fd2 = (plus.d1 - minus.d1) / (2.0 * h);
            let fd3 = (plus.d2 - minus.d2) / (2.0 * h);
            let close = |got: f64, want: f64| (got - want).abs() <= 1e-5 * want.abs().max(1e-2);
            assert!(
                close(centre.d1, fd1),
                "lambda {lambda}, x {x}: T' {} vs {fd1}",
                centre.d1
            );
            assert!(
                close(centre.d2, fd2),
                "lambda {lambda}, x {x}: T'' {} vs {fd2}",
                centre.d2
            );
            assert!(
                close(centre.d3, fd3),
                "lambda {lambda}, x {x}: T''' {} vs {fd3}",
                centre.d3
            );
        }
    }
}

#[test]
fn battin_series_agrees_with_the_lagrange_form_where_both_are_valid() {
    for &lambda in &LAMBDAS {
        for &x in &[0.78, 0.85, 0.95, 1.05, 1.12, 1.18] {
            let series = tof_battin(x, lambda);
            let lagrange = tof_lagrange(x, lambda);
            let close = |a: f64, b: f64, tol: f64| (a - b).abs() <= tol * b.abs().max(1e-3);
            assert!(
                close(series.time, lagrange.time, 1e-11),
                "lambda {lambda}, x {x}: T"
            );
            assert!(
                close(series.d1, lagrange.d1, 1e-8),
                "lambda {lambda}, x {x}: T'"
            );
            assert!(
                close(series.d2, lagrange.d2, 1e-6),
                "lambda {lambda}, x {x}: T''"
            );
            assert!(
                close(series.d3, lagrange.d3, 1e-4),
                "lambda {lambda}, x {x}: T'''"
            );
        }
    }
}

#[test]
fn time_of_flight_is_continuous_across_the_switch_to_the_series() {
    for &lambda in &LAMBDAS {
        for edge in [BATTIN_LOWER, BATTIN_UPPER] {
            let inside = tof_and_derivatives(edge + if edge < 1.0 { 1e-9 } else { -1e-9 }, lambda);
            let outside = tof_and_derivatives(edge + if edge < 1.0 { -1e-9 } else { 1e-9 }, lambda);
            assert!(
                (inside.time - outside.time).abs() < 1e-8,
                "lambda {lambda}: T jumps at {edge}"
            );
            assert!(
                (inside.d1 - outside.d1).abs() < 1e-7,
                "lambda {lambda}: T' jumps at {edge}"
            );
        }
    }
}

#[test]
fn hypergeometric_series_matches_its_closed_form() {
    // 2F1(3, 1; 5/2; z) at z = 0 is 1 and its first Taylor coefficient is 3/(5/2).
    let at_zero = hypergeometric_3_1_five_halves(Taylor3::variable(0.0));
    assert!((at_zero.0[0] - 1.0).abs() < 1e-15);
    assert!((at_zero.0[1] - 1.2).abs() < 1e-14);
    // Second coefficient: (3*4 * 1*2) / ((5/2)(7/2) * 2!) = 24 / 17.5.
    assert!((at_zero.0[2] - 24.0 / 17.5).abs() < 1e-14);
    // Identity 2F1(1, 1; 3/2; z) relation is not needed; check against direct summation.
    let z = -0.35;
    let direct: f64 = {
        let mut term = 1.0;
        let mut sum = 0.0;
        for k in 0..200 {
            sum += term;
            term *= (3.0 + f64::from(k)) / (2.5 + f64::from(k)) * z;
        }
        sum
    };
    let composed = hypergeometric_3_1_five_halves(Taylor3::constant(z));
    assert!((composed.0[0] - direct).abs() < 1e-14);
}

#[test]
fn taylor_arithmetic_reproduces_polynomial_derivatives() {
    // f(h) = sqrt(1 + 2 x) at x = 0.5 + h: f = sqrt(2), f' = 1/sqrt(2), f'' = -1/(2 sqrt(2)), f''' = 3/(4 sqrt(2)).
    let x = Taylor3::variable(0.5);
    let f = (Taylor3::constant(1.0) + x.scale(2.0)).sqrt();
    let root2 = 2.0_f64.sqrt();
    assert!((f.0[0] - root2).abs() < 1e-15);
    assert!((f.0[1] - 1.0 / root2).abs() < 1e-15);
    assert!((2.0 * f.0[2] + 1.0 / (2.0 * root2)).abs() < 1e-15);
    assert!((6.0 * f.0[3] - 3.0 / (4.0 * root2)).abs() < 1e-15);
    // Product rule: (x^2)'' = 2.
    let sq = x * x;
    assert!(
        (sq.0[0] - 0.25).abs() < 1e-15
            && (sq.0[1] - 1.0).abs() < 1e-15
            && (sq.0[2] - 1.0).abs() < 1e-15
    );
}
