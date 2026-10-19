// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

// Tests of the exact Riemann solver of the shock tube problem.
//
// The reference values are those of Toro, "Riemann Solvers and Numerical
// Methods for Fluid Dynamics", Table 4.3 (γ = 1.4), cross-checked with an
// independent root finder. Structural properties (mirror symmetry, Galilean
// invariance, Riemann invariants in the fans, Rankine-Hugoniot jumps across the
// shocks) are checked on top so every branch of the wave pattern is exercised:
// left/right shock, left/right rarefaction and vacuum.

use super::*;

const GAMMA: f64 = 1.4;

/// A shock tube with `(ρ, u, p)` left and right states, γ = 1.4 and the
/// diaphragm at `x = 0.5`.
fn problem(left: (f64, f64, f64), right: (f64, f64, f64)) -> SodShockTube {
    SodShockTube {
        left: ShockTubeState {
            rho: left.0,
            velocity: left.1,
            pressure: left.2,
        },
        right: ShockTubeState {
            rho: right.0,
            velocity: right.1,
            pressure: right.2,
        },
        gamma: GAMMA,
        length: 1.0,
        diaphragm: 0.5,
    }
}

/// The five test problems of Toro, Table 4.1, with their reference star
/// pressure and velocity `(p*, u*)` (computed with an independent root finder).
fn toro_problems() -> Vec<(SodShockTube, f64, f64)> {
    vec![
        (
            problem((1.0, 0.0, 1.0), (0.125, 0.0, 0.1)),
            0.30313017805065,
            0.92745262004895,
        ),
        (
            problem((1.0, -2.0, 0.4), (1.0, 2.0, 0.4)),
            0.0018938734200548,
            0.0,
        ),
        (
            problem((1.0, 0.0, 1000.0), (1.0, 0.0, 0.01)),
            460.89378749138,
            19.597451388723,
        ),
        (
            problem((5.99924, 19.5975, 460.894), (5.99242, -6.19633, 46.0950)),
            1691.6469553991,
            8.6897744116324,
        ),
        (
            problem((1.0, -19.59745, 1000.0), (1.0, -19.59745, 0.01)),
            460.89378749138,
            1.3887230601028e-06,
        ),
    ]
}

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

fn sound_speed(rho: f64, p: f64) -> f64 {
    (GAMMA * p / rho).sqrt()
}

#[test]
fn test_star_state_matches_the_textbook_values() {
    for (k, (tube, p_star, u_star)) in toro_problems().into_iter().enumerate() {
        let (p, u) = tube
            .star_state()
            .unwrap_or_else(|| panic!("problem {}: no star state", k + 1));
        assert!(
            (p - p_star).abs() < 1e-9 * p_star,
            "problem {}: p* = {p}, expected {p_star}",
            k + 1
        );
        assert!(
            (u - u_star).abs() < 1e-9 * (1.0 + u_star.abs()),
            "problem {}: u* = {u}, expected {u_star}",
            k + 1
        );
    }
}

/// Sod at `t = 0.2`: the contact moves to the right (`u* ≈ +0.92745`; the old
/// `½ (f_L − f_R)` gave the negated value), and the left rarefaction, contact
/// and right shock are at their textbook positions.
#[test]
fn test_sod_solution_has_the_textbook_wave_pattern() {
    let sod = SodShockTube::standard();
    let t = 0.2;
    let (rho_star_l, rho_star_r, u_star, p_star) =
        (0.426_319_43, 0.265_573_71, 0.927_452_62, 0.303_130_18);
    // Undisturbed states with their velocities.
    assert_eq!(sod.exact_solution(0.1, t), (1.0, 0.0, 1.0));
    assert_eq!(sod.exact_solution(0.95, t), (0.125, 0.0, 0.1));
    // Star states on both sides of the contact discontinuity at x = 0.5 + u* t.
    let contact = 0.5 + u_star * t;
    for (x, rho) in [(0.6, rho_star_l), (contact - 1e-6, rho_star_l)] {
        let (r, u, p) = sod.exact_solution(x, t);
        assert!((r - rho).abs() < 1e-7, "x = {x}: rho = {r}");
        assert!((u - u_star).abs() < 1e-7 && (p - p_star).abs() < 1e-7);
    }
    for (x, rho) in [(contact + 1e-6, rho_star_r), (0.75, rho_star_r)] {
        let (r, u, p) = sod.exact_solution(x, t);
        assert!((r - rho).abs() < 1e-7, "x = {x}: rho = {r}");
        assert!((u - u_star).abs() < 1e-7 && (p - p_star).abs() < 1e-7);
    }
    // Rarefaction head x = 0.5 − c_L t = 0.26336, tail x = 0.5 + (u* − c*_L) t.
    let head = 0.5 - sound_speed(1.0, 1.0) * t;
    assert_eq!(sod.exact_solution(head - 1e-6, t), (1.0, 0.0, 1.0));
    let (r_in, u_in, p_in) = sod.exact_solution(head + 1e-3, t);
    assert!(r_in < 1.0 && r_in > 0.99 && u_in > 0.0 && p_in < 1.0);
    let tail = 0.5 + (-0.070_272_813) * t;
    let (r_tail, u_tail, p_tail) = sod.exact_solution(tail + 1e-7, t);
    assert!((r_tail - rho_star_l).abs() < 1e-5 && (u_tail - u_star).abs() < 1e-5);
    assert!((p_tail - p_star).abs() < 1e-5);
    // Right shock at x = 0.5 + 1.7521557 t.
    let shock = 0.5 + 1.752_155_7 * t;
    assert_eq!(sod.exact_solution(shock + 1e-6, t), (0.125, 0.0, 0.1));
    let (r_behind, _, p_behind) = sod.exact_solution(shock - 1e-6, t);
    assert!((r_behind - rho_star_r).abs() < 1e-6 && (p_behind - p_star).abs() < 1e-6);
}

/// Speeds `(head, tail)` of the left rarefaction fan of `tube`, if it has one.
fn left_fan_speeds(tube: &SodShockTube) -> Option<(f64, f64)> {
    let g = tube.gamma;
    let (p_star, u_star) = tube.star_state()?;
    if p_star >= tube.left.pressure {
        return None;
    }
    let c_l = tube.left.sound_speed(g);
    let c_star = c_l * (p_star / tube.left.pressure).powf((g - 1.0) / (2.0 * g));
    Some((tube.left.velocity - c_l, u_star - c_star))
}

/// Speeds `(head, tail)` of the right rarefaction fan of `tube`, if it has one.
fn right_fan_speeds(tube: &SodShockTube) -> Option<(f64, f64)> {
    let g = tube.gamma;
    let (p_star, u_star) = tube.star_state()?;
    if p_star >= tube.right.pressure {
        return None;
    }
    let c_r = tube.right.sound_speed(g);
    let c_star = c_r * (p_star / tube.right.pressure).powf((g - 1.0) / (2.0 * g));
    Some((tube.right.velocity + c_r, u_star + c_star))
}

/// Inside a rarefaction fan the solution is an isentropic simple wave: it
/// carries the Riemann invariant of the outgoing family, has the entropy of
/// the undisturbed state, and sits on the characteristic `x − x₀ = (u ∓ c) t`.
#[test]
fn test_rarefaction_fans_are_isentropic_simple_waves() {
    let t = 0.3;
    let mut fans_checked = 0;
    for (tube, _, _) in toro_problems() {
        let g = tube.gamma;
        let (left, right) = (tube.left, tube.right);
        if let Some((head, tail)) = left_fan_speeds(&tube) {
            let c_l = left.sound_speed(g);
            let invariant = left.velocity + 2.0 * c_l / (g - 1.0);
            let entropy = left.pressure / left.rho.powf(g);
            for k in 1..20 {
                let speed = head + (tail - head) * k as f64 / 20.0;
                let (rho, u, p) = tube.exact_solution(tube.diaphragm + speed * t, t);
                let c = sound_speed(rho, p);
                assert!(
                    (speed - (u - c)).abs() < 1e-9 * (1.0 + speed.abs()),
                    "left fan"
                );
                assert!(
                    (u + 2.0 * c / (g - 1.0) - invariant).abs() < 1e-9 * (1.0 + invariant.abs())
                );
                assert!((p / rho.powf(g) - entropy).abs() < 1e-9 * entropy);
                // Between the undisturbed and the star state.
                assert!(p < left.pressure && rho < left.rho);
            }
            fans_checked += 1;
        }
        if let Some((head, tail)) = right_fan_speeds(&tube) {
            let c_r = right.sound_speed(g);
            let invariant = right.velocity - 2.0 * c_r / (g - 1.0);
            let entropy = right.pressure / right.rho.powf(g);
            for k in 1..20 {
                let speed = head + (tail - head) * k as f64 / 20.0;
                let (rho, u, p) = tube.exact_solution(tube.diaphragm + speed * t, t);
                let c = sound_speed(rho, p);
                assert!(
                    (speed - (u + c)).abs() < 1e-9 * (1.0 + speed.abs()),
                    "right fan"
                );
                assert!(
                    (u - 2.0 * c / (g - 1.0) - invariant).abs() < 1e-9 * (1.0 + invariant.abs())
                );
                assert!((p / rho.powf(g) - entropy).abs() < 1e-9 * entropy);
                assert!(p < right.pressure && rho < right.rho);
            }
            fans_checked += 1;
        }
    }
    // Sod (left), 123 (both), Toro 3 (left), Toro 5 (left).
    assert_eq!(fans_checked, 5);
}

/// The symmetric double rarefaction ("123 problem") has `u* = 0` and a
/// solution that is mirror symmetric about the diaphragm; a sign error in the
/// right fan breaks the symmetry.
#[test]
fn test_double_rarefaction_is_mirror_symmetric() {
    let tube = problem((1.0, -2.0, 0.4), (1.0, 2.0, 0.4));
    let (p_star, u_star) = tube.star_state().expect("no vacuum");
    assert!(u_star.abs() < 1e-12, "u* = {u_star}");
    assert!((p_star - 0.001_893_873_4).abs() < 1e-8);
    let t = 0.15;
    for k in 1..400 {
        let xi = k as f64 * 0.005;
        let (r_l, u_l, p_l) = tube.exact_solution(0.5 - xi, t);
        let (r_r, u_r, p_r) = tube.exact_solution(0.5 + xi, t);
        assert!((r_l - r_r).abs() < 1e-10 * r_l, "density at xi = {xi}");
        assert!((u_l + u_r).abs() < 1e-10, "velocity at xi = {xi}");
        assert!((p_l - p_r).abs() < 1e-10 * p_l, "pressure at xi = {xi}");
    }
    // The centre is the evacuated low-pressure star region.
    let (r0, u0, p0) = tube.exact_solution(0.5, t);
    assert!(u0.abs() < 1e-12);
    assert!((p0 - p_star).abs() < 1e-12);
    assert!(r0 < 0.1, "central density {r0}");
}

/// Swapping the two states and reflecting `x` reflects the solution
/// (`ρ` and `p` even, `u` odd). The mirrored Sod problem has a left shock and a
/// right rarefaction, the branches the standard problem never reaches.
#[test]
fn test_mirrored_problem_gives_the_mirrored_solution() {
    for (tube, _, _) in toro_problems() {
        let mirrored = problem(
            (tube.right.rho, -tube.right.velocity, tube.right.pressure),
            (tube.left.rho, -tube.left.velocity, tube.left.pressure),
        );
        let t = 0.04;
        for k in 0..600 {
            let xi = -1.5 + 3.0 * k as f64 / 599.0;
            let (rho, u, p) = tube.exact_solution(0.5 + xi, t);
            let (rho_m, u_m, p_m) = mirrored.exact_solution(0.5 - xi, t);
            assert!(
                (rho - rho_m).abs() < 1e-9 * rho
                    && (u + u_m).abs() < 1e-9 * u.abs().max(1.0)
                    && (p - p_m).abs() < 1e-9 * p,
                "xi = {xi}: ({rho}, {u}, {p}) vs mirrored ({rho_m}, {u_m}, {p_m})"
            );
        }
    }
}

/// Galilean invariance: adding a uniform velocity `U` to both initial states
/// moves the whole solution with `U` (this exercises the `u_L`, `u_R` terms of
/// the shock and rarefaction speeds).
#[test]
fn test_uniform_velocity_boost_translates_the_solution() {
    let t = 0.05;
    for boost in [-2.5_f64, 3.0] {
        for (tube, _, _) in toro_problems() {
            let boosted = problem(
                (
                    tube.left.rho,
                    tube.left.velocity + boost,
                    tube.left.pressure,
                ),
                (
                    tube.right.rho,
                    tube.right.velocity + boost,
                    tube.right.pressure,
                ),
            );
            for k in 0..500 {
                let x = -1.0 + 3.0 * k as f64 / 499.0;
                let (rho, u, p) = tube.exact_solution(x - boost * t, t);
                let (rho_b, u_b, p_b) = boosted.exact_solution(x, t);
                assert!(
                    (rho - rho_b).abs() < 1e-9 * rho
                        && (u + boost - u_b).abs() < 1e-9 * u_b.abs().max(1.0)
                        && (p - p_b).abs() < 1e-9 * p,
                    "boost {boost}, x = {x}: ({rho}, {u}, {p}) vs ({rho_b}, {u_b}, {p_b})"
                );
            }
        }
    }
}

/// Mass, momentum and energy are conserved across every shock of the exact
/// solution (the jump location is found from the solution itself).
#[test]
fn test_rankine_hugoniot_conditions_hold_across_every_shock() {
    let problems = [
        (SodShockTube::standard(), 0.2),
        // Mirrored Sod: left shock and right rarefaction.
        (problem((0.125, 0.0, 0.1), (1.0, 0.0, 1.0)), 0.2),
        // Toro 3: left rarefaction, contact, right shock.
        (problem((1.0, 0.0, 1000.0), (1.0, 0.0, 0.01)), 0.012),
        // Toro 4: two shocks.
        (
            problem((5.99924, 19.5975, 460.894), (5.99242, -6.19633, 46.0950)),
            0.035,
        ),
    ];
    let mut shocks_checked = 0;
    for (tube, t) in problems {
        let g = tube.gamma;
        let grid = 4000;
        let (x_min, x_max) = (-1.0, 2.0);
        let sample = |x: f64| tube.exact_solution(x, t);
        let p_max = (0..=grid)
            .map(|k| sample(x_min + (x_max - x_min) * k as f64 / grid as f64).2)
            .fold(0.0f64, f64::max);
        for k in 0..grid {
            let xa = x_min + (x_max - x_min) * k as f64 / grid as f64;
            let xb = x_min + (x_max - x_min) * (k + 1) as f64 / grid as f64;
            let (pa, pb) = (sample(xa).2, sample(xb).2);
            if (pa - pb).abs() < 0.1 * p_max {
                continue;
            }
            // A jump between xa and xb: bisect for its location.
            let (mut lo, mut hi) = (xa, xb);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if (sample(mid).2 - pa).abs() < (sample(mid).2 - pb).abs() {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let location = 0.5 * (lo + hi);
            let speed = (location - tube.diaphragm) / t;
            let a = sample(location - 1e-7);
            let b = sample(location + 1e-7);
            let total_energy = |s: (f64, f64, f64)| s.2 / (g - 1.0) + 0.5 * s.0 * s.1 * s.1;
            let flux_mass = |s: (f64, f64, f64)| s.0 * (s.1 - speed);
            let flux_momentum = |s: (f64, f64, f64)| s.0 * (s.1 - speed) * s.1 + s.2;
            let flux_energy = |s: (f64, f64, f64)| total_energy(s) * (s.1 - speed) + s.2 * s.1;
            for (name, fa, fb) in [
                ("mass", flux_mass(a), flux_mass(b)),
                ("momentum", flux_momentum(a), flux_momentum(b)),
                ("energy", flux_energy(a), flux_energy(b)),
            ] {
                let scale = fa.abs().max(fb.abs()).max(1e-12);
                assert!(
                    (fa - fb).abs() < 1e-6 * scale,
                    "{name} not conserved across the shock at x = {location}: {fa} vs {fb}"
                );
            }
            shocks_checked += 1;
        }
    }
    assert_eq!(shocks_checked, 1 + 1 + 1 + 2, "unexpected number of shocks");
}

/// Two strong diverging rarefactions leave a vacuum: no star state, a
/// vacuum core, and fans whose density and pressure vanish at their tails.
#[test]
fn test_diverging_rarefactions_generate_a_vacuum() {
    let tube = problem((1.0, -10.0, 0.4), (1.0, 10.0, 0.4));
    assert_eq!(tube.star_state(), None);
    let t = 0.1;
    let g = GAMMA;
    let c = sound_speed(1.0, 0.4);
    let (head_left, tail_left) = (-10.0 - c, -10.0 + 2.0 * c / (g - 1.0));
    let (tail_right, head_right) = (10.0 - 2.0 * c / (g - 1.0), 10.0 + c);
    assert!(
        tail_left < tail_right,
        "the pair must be in the vacuum regime"
    );
    let at = |speed: f64| tube.exact_solution(0.5 + speed * t, t);
    // Undisturbed states outside the heads.
    assert_eq!(at(head_left - 0.5), (1.0, -10.0, 0.4));
    assert_eq!(at(head_right + 0.5), (1.0, 10.0, 0.4));
    // Empty core between the tails.
    assert_eq!(at(0.0), (0.0, 0.0, 0.0));
    assert_eq!(at(0.5 * (tail_left + tail_right)), (0.0, 0.0, 0.0));
    // Fans are simple waves whose density and pressure vanish at the tail.
    for k in 1..20 {
        let speed = head_left + (tail_left - head_left) * k as f64 / 20.0;
        let (rho, u, p) = at(speed);
        let c_local = sound_speed(rho, p);
        assert!((speed - (u - c_local)).abs() < 1e-9, "left fan at {speed}");
        let invariant = -10.0 + 2.0 * c / (g - 1.0);
        assert!((u + 2.0 * c_local / (g - 1.0) - invariant).abs() < 1e-9);
        let speed = head_right + (tail_right - head_right) * k as f64 / 20.0;
        let (rho, u, p) = at(speed);
        let c_local = sound_speed(rho, p);
        assert!((speed - (u + c_local)).abs() < 1e-9, "right fan at {speed}");
        let invariant = 10.0 - 2.0 * c / (g - 1.0);
        assert!((u - 2.0 * c_local / (g - 1.0) - invariant).abs() < 1e-9);
    }
    for edge in [tail_left - 1e-9, tail_right + 1e-9] {
        let (rho, _, p) = at(edge);
        assert!(rho < 1e-6 && p < 1e-6, "edge {edge}: rho = {rho}, p = {p}");
    }
}

/// For `t ≤ 0` the initial states are returned, with their velocities.
#[test]
fn test_initial_states_are_returned_at_time_zero() {
    let tube = problem((1.0, 0.7, 2.0), (0.5, -0.3, 0.4));
    assert_eq!(tube.exact_solution(0.4, 0.0), (1.0, 0.7, 2.0));
    assert_eq!(tube.exact_solution(0.6, 0.0), (0.5, -0.3, 0.4));
    assert_eq!(tube.exact_solution(0.4, -1.0), (1.0, 0.7, 2.0));
}

/// Unphysical states have no solution; every component is `NaN`.
#[test]
fn test_unphysical_states_give_nan() {
    let bad = [
        problem((-1.0, 0.0, 1.0), (0.125, 0.0, 0.1)),
        problem((1.0, 0.0, 0.0), (0.125, 0.0, 0.1)),
        problem((1.0, 0.0, 1.0), (0.0, 0.0, 0.1)),
        problem((1.0, f64::NAN, 1.0), (0.125, 0.0, 0.1)),
    ];
    for tube in bad {
        let (rho, u, p) = tube.exact_solution(0.7, 0.2);
        assert!(rho.is_nan() && u.is_nan() && p.is_nan());
        assert_eq!(tube.star_state(), None);
    }
    let mut cold = SodShockTube::standard();
    cold.gamma = 1.0;
    assert_eq!(cold.star_state(), None);
}

/// The Newton iteration converges to the root of the pressure function and
/// the sampled solution stays physical for arbitrary (random) states.
#[test]
fn test_star_pressure_is_the_root_for_random_problems() {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut vacuum_cases = 0;
    let mut solved = 0;
    for _ in 0..400 {
        let tube = problem(
            (
                0.1 + 5.0 * lcg_unit(&mut state),
                6.0 * (lcg_unit(&mut state) - 0.5),
                0.01 + 50.0 * lcg_unit(&mut state).powi(3),
            ),
            (
                0.1 + 5.0 * lcg_unit(&mut state),
                6.0 * (lcg_unit(&mut state) - 0.5),
                0.01 + 50.0 * lcg_unit(&mut state).powi(3),
            ),
        );
        let du = tube.right.velocity - tube.left.velocity;
        let c_l = tube.left.sound_speed(GAMMA);
        let c_r = tube.right.sound_speed(GAMMA);
        let Some((p, u)) = tube.star_state() else {
            assert!(
                2.0 * (c_l + c_r) / (GAMMA - 1.0) <= du,
                "no star state without a vacuum"
            );
            vacuum_cases += 1;
            continue;
        };
        solved += 1;
        let (f_l, _) = pressure_function(p, &tube.left, GAMMA);
        let (f_r, _) = pressure_function(p, &tube.right, GAMMA);
        let scale = f_l.abs() + f_r.abs() + du.abs();
        assert!(
            (f_l + f_r + du).abs() < 1e-9 * scale,
            "residual {} (scale {scale})",
            f_l + f_r + du
        );
        assert!(p > 0.0);
        assert!((u - (tube.left.velocity - f_l)).abs() < 1e-9 * u.abs().max(1.0));
        assert!((u - (tube.right.velocity + f_r)).abs() < 1e-9 * u.abs().max(1.0));
        for k in 0..60 {
            let x = 0.5 + (k as f64 / 59.0 * 2.0 - 1.0) * 5.0;
            let (rho, vel, pres) = tube.exact_solution(x, 0.5);
            assert!(
                rho > 0.0 && pres > 0.0 && vel.is_finite(),
                "x = {x}: {rho} {vel} {pres}"
            );
        }
    }
    assert!(
        solved > 300,
        "only {solved} problems solved ({vacuum_cases} vacuum)"
    );
}
