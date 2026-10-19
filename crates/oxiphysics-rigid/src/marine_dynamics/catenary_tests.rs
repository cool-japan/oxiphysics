// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the elastic catenary solver and of `MooringLine`.
//!
//! The solver is checked against closed-form inextensible catenary relations
//! and against an independent numerical integration of the catenary ODE.

use super::*;

const INEXTENSIBLE: f64 = f64::INFINITY;

fn solve(length: f64, weight: f64, stiffness: f64, span: f64, height: f64) -> CatenarySolution {
    match solve_static_catenary(length, weight, stiffness, span, height) {
        Ok(solution) => solution,
        Err(error) => {
            panic!("solver failed for ({length}, {weight}, {stiffness}, {span}, {height}): {error}")
        }
    }
}

/// Simpson integration of `f` over `[a, b]` with `n` (even) intervals.
fn simpson<F: Fn(f64) -> (f64, f64)>(f: F, a: f64, b: f64, n: usize) -> (f64, f64) {
    if b <= a {
        return (0.0, 0.0);
    }
    let h = (b - a) / n as f64;
    let (mut sx, mut sz) = (0.0, 0.0);
    for i in 0..=n {
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let (dx, dz) = f(a + i as f64 * h);
        sx += weight * dx;
        sz += weight * dz;
    }
    (sx * h / 3.0, sz * h / 3.0)
}

/// Fairlead position obtained by integrating the elastic catenary ODE
/// `dx/ds = (1 + T/EA)·H/T`, `dz/ds = (1 + T/EA)·V(s)/T` along the unstretched
/// arc length, with `V(s) = V_A + w·max(0, s - L_B)`.
///
/// The suspended part is integrated on panels whose width doubles away from
/// the touchdown point, starting with the catenary parameter `a = H/w`: the
/// integrand varies on the scale `a` next to the touchdown knee, so a uniform
/// step would under-resolve lines that are only just taut.
fn integrate_fairlead(solution: &CatenarySolution) -> (f64, f64) {
    let inv_ea = if solution.axial_stiffness.is_finite() {
        1.0 / solution.axial_stiffness
    } else {
        0.0
    };
    let h = solution.horizontal_tension;
    let w = solution.weight_per_length;
    let l_b = solution.grounded_length;
    let derivatives = |s: f64| {
        let v = solution.anchor_vertical_tension + w * (s - l_b).max(0.0);
        let t = h.hypot(v);
        let stretch = 1.0 + t * inv_ea;
        (stretch * h / t, stretch * v / t)
    };
    let (mut x, mut z) = simpson(derivatives, 0.0, l_b, 2000);
    let mut start = l_b;
    let mut width = (h / w).max(1e-9 * solution.length);
    while start < solution.length {
        let end = (start + width).min(solution.length);
        let (dx, dz) = simpson(derivatives, start, end, 400);
        x += dx;
        z += dz;
        start = end;
        width *= 2.0;
    }
    (x, z)
}

// ── inextensible closed forms ───────────────────────────────────────────

#[test]
fn inextensible_line_on_the_seabed_matches_the_closed_form() {
    let (length, w, z) = (150.0, 200.0, 50.0);
    for &x in &[101.0, 110.0, 125.0, 135.0, 138.0] {
        let solution = solve(length, w, INEXTENSIBLE, x, z);
        let h = solution.horizontal_tension;
        let a = h / w;
        let suspended = (z * (z + 2.0 * a)).sqrt();
        assert!(h > 0.0 && h.is_finite(), "X = {x}: H = {h}");
        // X = L - L_s + a asinh(L_s / a), V = w L_s.
        let span = length - suspended + a * (suspended / a).asinh();
        assert!((span - x).abs() < 1e-8, "X = {x}: closed-form span {span}");
        assert!(
            (solution.vertical_tension - w * suspended).abs() < 1e-6 * w * suspended,
            "X = {x}: V = {} vs w·L_s = {}",
            solution.vertical_tension,
            w * suspended
        );
        // z = a (cosh(x_s / a) - 1) on the suspended part.
        let x_s = a * (suspended / a).asinh();
        assert!((a * ((x_s / a).cosh() - 1.0) - z).abs() < 1e-6, "X = {x}");
        // Tension at the fairlead: T = w (a + z).
        assert!(
            (solution.fairlead_tension() - w * (a + z)).abs() < 1e-6 * w * (a + z),
            "X = {x}"
        );
        assert!((solution.grounded_length - (length - suspended)).abs() < 1e-8);
        assert_eq!(solution.anchor_vertical_tension, 0.0);
        assert!((solution.suspended_length() - suspended).abs() < 1e-8);
        // Horizontal tension is uniform along the grounded part.
        assert!((solution.anchor_tension() - h).abs() < 1e-9 * h);
        assert!((solution.catenary_parameter() - a).abs() < 1e-12);
    }
}

#[test]
fn inextensible_fully_suspended_line_matches_the_catenary_chord_relation() {
    // Between the touchdown limit (138.63 m) and the taut limit (141.42 m).
    let (length, w, z) = (150.0, 200.0, 50.0);
    for &x in &[139.0, 140.0, 141.0] {
        let solution = solve(length, w, INEXTENSIBLE, x, z);
        let a = solution.horizontal_tension / w;
        assert_eq!(solution.grounded_length, 0.0, "X = {x}");
        assert!(solution.anchor_vertical_tension > 0.0, "X = {x}");
        // sqrt(L² - Z²) = 2 a sinh(X / 2a) for a catenary segment of length L.
        let chord_relation = 2.0 * a * (x / (2.0 * a)).sinh();
        let expected = (length * length - z * z).sqrt();
        assert!(
            (chord_relation - expected).abs() < 1e-8 * expected,
            "X = {x}: {chord_relation} vs {expected}"
        );
        // V - V_A = w L and T_F - T_A = w Z.
        assert!(
            (solution.vertical_tension - solution.anchor_vertical_tension - w * length).abs()
                < 1e-8 * w * length
        );
        assert!(
            (solution.fairlead_tension() - solution.anchor_tension() - w * z).abs() < 1e-8 * w * z
        );
    }
}

#[test]
fn horizontal_tension_grows_with_the_span_and_is_continuous_at_touchdown() {
    let (length, w, z) = (150.0, 200.0, 50.0);
    // Touchdown limit: L_s = L, a* = (L² - Z²) / (2 Z) = 200.
    let a_star: f64 = 200.0;
    let touchdown_span = a_star * (length / a_star).asinh();
    assert!((touchdown_span - 138.6294).abs() < 1e-3);

    let mut previous = 0.0;
    let mut previous_grounded = f64::INFINITY;
    let mut x = 100.5;
    while x < 141.0 {
        let solution = solve(length, w, INEXTENSIBLE, x, z);
        let h = solution.horizontal_tension;
        assert!(
            h > previous,
            "H must increase with X: X = {x}, H = {h}, previous = {previous}"
        );
        assert!(solution.grounded_length <= previous_grounded + 1e-9);
        if x < touchdown_span - 1e-6 {
            assert!(solution.grounded_length > 0.0, "X = {x}");
            assert_eq!(solution.anchor_vertical_tension, 0.0);
        }
        previous = h;
        previous_grounded = solution.grounded_length;
        x += 0.25;
    }

    let below = solve(length, w, INEXTENSIBLE, touchdown_span - 1e-3, z);
    let above = solve(length, w, INEXTENSIBLE, touchdown_span + 1e-3, z);
    assert!(
        (above.horizontal_tension - below.horizontal_tension).abs()
            < 0.01 * below.horizontal_tension,
        "H jumps across the touchdown transition: {} -> {}",
        below.horizontal_tension,
        above.horizontal_tension
    );
    assert!(below.grounded_length < 0.1);
    assert!(above.anchor_vertical_tension < 0.1 * w * length);
    assert!(
        (above.vertical_tension - below.vertical_tension).abs() < 0.01 * below.vertical_tension
    );
}

// ── elastic catenary ────────────────────────────────────────────────────

#[test]
fn elastic_solutions_satisfy_the_catenary_ode() {
    // (length, weight N/m, EA N, X, Z): several partially grounded lines and
    // two lines that are stretched beyond their unstretched chord (fully
    // suspended, no seabed contact).
    let cases = [
        (150.0, 200.0, 5.0e6, 125.0, 50.0),
        (300.0, 1000.0, 1.0e8, 250.0, 100.0),
        (150.0, 200.0, 2.0e6, 145.0, 40.0),
        (150.0, 200.0, 5.0e6, 142.0, 50.0),
        (100.0, 50.0, 1.0e6, 80.0, 40.0),
        (500.0, 150.0, 8.0e7, 470.0, 120.0),
    ];
    for &(length, w, ea, x, z) in &cases {
        let solution = solve(length, w, ea, x, z);
        let (xi, zi) = integrate_fairlead(&solution);
        assert!(
            (xi - x).abs() < 1e-6 * x && (zi - z).abs() < 1e-6 * z.max(1.0),
            "case ({length}, {w}, {ea}, {x}, {z}): integrated fairlead ({xi}, {zi}), \
             solution {solution:?}"
        );
        // The fairlead carries the weight of the suspended length.
        let carried = w * solution.suspended_length();
        let balance = solution.vertical_tension - solution.anchor_vertical_tension;
        assert!(
            (balance - carried).abs() < 1e-6 * carried.max(1.0),
            "vertical force balance: {balance} vs {carried}"
        );
    }
}

#[test]
fn elasticity_lowers_the_tension_and_approaches_the_inextensible_limit() {
    let (length, w, x, z) = (150.0, 200.0, 135.0, 50.0);
    let rigid = solve(length, w, INEXTENSIBLE, x, z).horizontal_tension;
    let soft = solve(length, w, 2.0e5, x, z).horizontal_tension;
    let medium = solve(length, w, 2.0e7, x, z).horizontal_tension;
    let stiff = solve(length, w, 2.0e12, x, z).horizontal_tension;
    assert!(soft < medium && medium < stiff && stiff <= rigid * (1.0 + 1e-9));
    assert!(
        (stiff - rigid).abs() < 1e-4 * rigid,
        "EA = 2e12 should match the rigid line: {stiff} vs {rigid}"
    );
}

#[test]
fn stretched_taut_line_reaches_spans_beyond_the_unstretched_length() {
    // Chord 150.4 m > L = 150 m: only an elastic line has an equilibrium.
    let (length, w, ea, x, z) = (150.0_f64, 200.0_f64, 2.0e6_f64, 145.0_f64, 40.0_f64);
    assert!(x.hypot(z) > length);
    assert!(matches!(
        solve_static_catenary(length, w, INEXTENSIBLE, x, z),
        Err(CatenaryError::LineTooShort { .. })
    ));
    let solution = solve(length, w, ea, x, z);
    assert!(solution.horizontal_tension > 0.0 && solution.horizontal_tension.is_finite());
    // The line stretches under its own tension, so part of it can still rest
    // on the seabed even though the chord exceeds the unstretched length.
    let (xi, zi) = integrate_fairlead(&solution);
    assert!(
        (xi - x).abs() < 1e-6 * x && (zi - z).abs() < 1e-6 * z,
        "({xi}, {zi})"
    );
    // The mean tension must at least stretch the line to the chord, and the
    // fairlead carries the largest tension of the line.
    let strain = x.hypot(z) / length - 1.0;
    assert!(strain > 0.0);
    assert!(
        solution.fairlead_tension() > strain * ea,
        "fairlead tension {} must exceed EA·strain = {}",
        solution.fairlead_tension(),
        strain * ea
    );
    assert!(solution.fairlead_tension() > solution.anchor_tension());
}

#[test]
fn weightless_line_is_a_taut_elastic_string() {
    let (length, ea, x, z) = (100.0_f64, 1.0e6_f64, 60.0_f64, 90.0_f64);
    let chord = x.hypot(z);
    let solution = solve(length, 0.0, ea, x, z);
    let tension = ea * (chord / length - 1.0);
    assert!((solution.fairlead_tension() - tension).abs() < 1e-9 * tension);
    assert!((solution.horizontal_tension - tension * x / chord).abs() < 1e-9 * tension);
    assert!((solution.vertical_tension - tension * z / chord).abs() < 1e-9 * tension);
    // Slack weightless line: no tension at all.
    let slack = solve(length, 0.0, ea, 30.0, 40.0);
    assert_eq!(slack.horizontal_tension, 0.0);
    assert_eq!(slack.vertical_tension, 0.0);
    let slack_rigid = solve(length, 0.0, INEXTENSIBLE, 30.0, 40.0);
    assert_eq!(slack_rigid.fairlead_tension(), 0.0);
}

// ── slack lines and degenerate geometry ─────────────────────────────────

#[test]
fn slack_line_hangs_vertically_without_horizontal_tension() {
    let (length, w, z) = (150.0, 200.0, 50.0);
    for &stiffness in &[INEXTENSIBLE, 1.0e9, 1.0e6] {
        for &x in &[0.0, 40.0, 90.0] {
            let solution = solve(length, w, stiffness, x, z);
            assert_eq!(
                solution.horizontal_tension, 0.0,
                "EA = {stiffness}, X = {x}"
            );
            assert!(solution.vertical_tension.is_finite());
            // The hanging part carries its own weight: about w·Z.
            assert!(
                (solution.vertical_tension - w * z).abs() < 0.01 * w * z,
                "EA = {stiffness}, X = {x}: V = {}",
                solution.vertical_tension
            );
            assert!(solution.grounded_length > 0.0);
        }
    }
}

#[test]
fn line_on_the_seabed_at_zero_height_carries_no_vertical_load() {
    let solution = solve(100.0, 150.0, INEXTENSIBLE, 60.0, 0.0);
    assert_eq!(solution.horizontal_tension, 0.0);
    assert_eq!(solution.vertical_tension, 0.0);
    assert_eq!(solution.grounded_length, 100.0);
    // Taut along the seabed with an elastic line: T = EA (X / L - 1).
    let taut = solve(100.0, 150.0, 1.0e6, 100.5, 0.0);
    assert!((taut.horizontal_tension - 1.0e6 * 0.005).abs() < 1e-6 * 5000.0);
    assert_eq!(taut.vertical_tension, 0.0);
}

#[test]
fn invalid_parameters_are_reported() {
    let err = |l: f64, w: f64, ea: f64, x: f64, z: f64| solve_static_catenary(l, w, ea, x, z);
    for (result, name) in [
        (err(0.0, 1.0, 1e6, 1.0, 1.0), "length"),
        (err(-5.0, 1.0, 1e6, 1.0, 1.0), "length"),
        (err(f64::NAN, 1.0, 1e6, 1.0, 1.0), "length"),
        (err(f64::INFINITY, 1.0, 1e6, 1.0, 1.0), "length"),
        (err(10.0, -1.0, 1e6, 1.0, 1.0), "weight_per_length"),
        (err(10.0, f64::NAN, 1e6, 1.0, 1.0), "weight_per_length"),
        (err(10.0, 1.0, 0.0, 1.0, 1.0), "axial_stiffness"),
        (err(10.0, 1.0, -1.0, 1.0, 1.0), "axial_stiffness"),
        (err(10.0, 1.0, f64::NAN, 1.0, 1.0), "axial_stiffness"),
        (err(10.0, 1.0, 1e6, -1.0, 1.0), "horizontal_span"),
        (err(10.0, 1.0, 1e6, f64::NAN, 1.0), "horizontal_span"),
        (err(10.0, 1.0, 1e6, 1.0, -1.0), "vertical_span"),
        (err(10.0, 1.0, 1e6, 1.0, f64::INFINITY), "vertical_span"),
    ] {
        match result {
            Err(CatenaryError::InvalidParameter { name: reported, .. }) => {
                assert_eq!(reported, name);
            }
            other => panic!("expected an invalid `{name}`, got {other:?}"),
        }
    }
    // Inextensible lines cannot span a chord that is not shorter than the line.
    assert!(matches!(
        solve_static_catenary(100.0, 100.0, INEXTENSIBLE, 100.0, 0.0),
        Err(CatenaryError::LineTooShort { .. })
    ));
    assert!(matches!(
        solve_static_catenary(100.0, 100.0, INEXTENSIBLE, 80.0, 70.0),
        Err(CatenaryError::LineTooShort { .. })
    ));
    assert!(solve_static_catenary(100.0, 100.0, 1.0e7, 80.0, 70.0).is_ok());
    let message = solve_static_catenary(100.0, 100.0, INEXTENSIBLE, 100.0, 0.0)
        .map(|_| ())
        .map_err(|e| e.to_string());
    assert!(message.is_err_and(|m| m.contains("cannot span")));
}

// ── line profile ────────────────────────────────────────────────────────

#[test]
fn profile_starts_at_the_anchor_and_ends_at_the_fairlead() {
    let cases = [
        (150.0, 200.0, INEXTENSIBLE, 125.0, 50.0),
        (150.0, 200.0, 5.0e6, 125.0, 50.0),
        (150.0, 200.0, INEXTENSIBLE, 140.0, 50.0),
        (150.0, 200.0, 2.0e6, 145.0, 40.0),
        (150.0, 200.0, 5.0e6, 60.0, 50.0),
    ];
    for &(length, w, ea, x, z) in &cases {
        let solution = solve(length, w, ea, x, z);
        let start = solution.position_at(0.0);
        let end = solution.position_at(length);
        assert!(
            start[0].abs() < 1e-12 && start[1].abs() < 1e-12,
            "{start:?}"
        );
        assert!(
            (end[0] - x).abs() < 1e-6 * x && (end[1] - z).abs() < 1e-6 * z,
            "case ({length}, {w}, {ea}, {x}, {z}): end = {end:?}"
        );
        // Clamped outside the line.
        assert_eq!(solution.position_at(-3.0), start);
        assert_eq!(solution.position_at(length + 3.0), end);
        // The profile never dips below the seabed and x never moves backwards.
        let mut previous = start;
        for k in 1..=200 {
            let point = solution.position_at(length * k as f64 / 200.0);
            assert!(point[1] >= -1e-9, "case {x}: z = {}", point[1]);
            assert!(point[0] >= previous[0] - 1e-9, "x decreased");
            assert!(point[1] >= previous[1] - 1e-9, "z decreased");
            previous = point;
        }
    }
}

// ── randomized sweep ────────────────────────────────────────────────────

/// Deterministic generator of uniform numbers in `[0, 1)`.
fn next_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

/// A random line geometry `(length, weight, stiffness, span, height)`.
///
/// The span runs from 0.3 to 1.1 times the largest span an inextensible line
/// of that length can reach at the given height, so the sweep covers slack,
/// partly grounded, fully suspended and over-stretched configurations; a
/// quarter of the lines are inextensible and the others have `EA/(w·L)`
/// between 10 and 1e5.
fn random_case(state: &mut u64) -> (f64, f64, f64, f64, f64) {
    let length = 20.0 + 980.0 * next_unit(state);
    let weight = 10.0 + 2000.0 * next_unit(state);
    let height = 0.9 * length * next_unit(state);
    let reach = (length * length - height * height).sqrt();
    let span = reach * (0.3 + 0.8 * next_unit(state));
    let stiffness = if next_unit(state) < 0.25 {
        INEXTENSIBLE
    } else {
        weight * length * 10.0_f64.powf(1.0 + 4.0 * next_unit(state))
    };
    (length, weight, stiffness, span, height)
}

#[test]
fn randomized_geometries_close_the_catenary_equations() {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let (mut slack, mut grounded, mut suspended, mut rejected) = (0, 0, 0, 0);
    for case in 0..600 {
        let (length, w, ea, x, z) = random_case(&mut state);
        let solution = match solve_static_catenary(length, w, ea, x, z) {
            Ok(solution) => solution,
            Err(CatenaryError::LineTooShort { .. }) => {
                assert!(
                    ea.is_infinite() && x.hypot(z) >= length,
                    "case {case} ({length}, {w}, {ea}, {x}, {z}) wrongly reported as too short"
                );
                rejected += 1;
                continue;
            }
            Err(error) => panic!("case {case} ({length}, {w}, {ea}, {x}, {z}): {error}"),
        };
        let describe = format!("case {case} ({length}, {w}, {ea}, {x}, {z}) -> {solution:?}");
        // Vertical force balance: the fairlead carries the suspended weight.
        let carried = w * solution.suspended_length();
        let balance = solution.vertical_tension - solution.anchor_vertical_tension;
        assert!((balance - carried).abs() <= 1e-9 * w * length, "{describe}");
        if solution.horizontal_tension > 0.0 {
            let (xi, zi) = integrate_fairlead(&solution);
            assert!(
                (xi - x).abs() <= 1e-6 * x.max(1.0) && (zi - z).abs() <= 1e-6 * z.max(1.0),
                "integrated fairlead ({xi}, {zi}): {describe}"
            );
            if solution.grounded_length > 0.0 {
                grounded += 1;
            } else {
                suspended += 1;
            }
        } else {
            // Slack: the line hangs vertically from the fairlead, and the
            // surplus length is bunched up on the seabed within the span.
            slack += 1;
            let inv_ea = if ea.is_finite() { 1.0 / ea } else { 0.0 };
            let v = solution.vertical_tension;
            let hanging_height = v / w + 0.5 * v * v * inv_ea / w;
            assert!(
                (hanging_height - z).abs() <= 1e-9 * z.max(1.0),
                "{describe}"
            );
            assert!(solution.grounded_length >= x - 1e-9 * length, "{describe}");
            assert_eq!(solution.anchor_vertical_tension, 0.0, "{describe}");
        }
    }
    assert!(
        slack > 10 && grounded > 50 && suspended > 5 && rejected > 5,
        "the sweep must exercise every regime: slack {slack}, grounded {grounded}, \
         suspended {suspended}, rejected {rejected}"
    );
}

#[test]
fn tensions_never_decrease_when_the_fairlead_moves_away() {
    let mut state = 0xD1B5_4A32_D192_ED03_u64;
    for case in 0..120 {
        let (length, w, ea, _, z) = random_case(&mut state);
        let reach = (length * length - z * z).sqrt();
        let scale = w * length;
        let (mut horizontal, mut vertical, mut fairlead) = (0.0, 0.0, 0.0);
        for k in 0..=12 {
            let x = reach * (0.3 + 0.69 * k as f64 / 12.0);
            let solution = solve(length, w, ea, x, z);
            let describe = format!("case {case} k = {k} ({length}, {w}, {ea}, {x}, {z})");
            assert!(
                solution.horizontal_tension >= horizontal - 1e-9 * scale,
                "H decreased: {} -> {}: {describe}",
                horizontal,
                solution.horizontal_tension
            );
            assert!(
                solution.vertical_tension >= vertical - 1e-9 * scale,
                "V decreased: {} -> {}: {describe}",
                vertical,
                solution.vertical_tension
            );
            assert!(
                solution.fairlead_tension() >= fairlead - 1e-9 * scale,
                "T decreased: {describe}"
            );
            horizontal = solution.horizontal_tension;
            vertical = solution.vertical_tension;
            fairlead = solution.fairlead_tension();
        }
    }
}

// ── MooringLine ─────────────────────────────────────────────────────────

#[test]
fn mooring_line_tension_matches_the_closed_form_for_a_stiff_line() {
    let line = MooringLine::new(150.0, 200.0, 1.0e12, [0.0, 0.0, -50.0]);
    let (t_h, t_v) = line.catenary_tension(125.0, 50.0);
    let a = t_h / 200.0;
    let suspended = (50.0 * (50.0 + 2.0 * a)).sqrt();
    let span = 150.0 - suspended + a * (suspended / a).asinh();
    assert!((span - 125.0).abs() < 1e-5, "closure of X: {span}");
    assert!((t_v - 200.0 * suspended).abs() < 1e-3 * t_v);
}

#[test]
fn mooring_restoring_force_points_to_the_anchor_and_down() {
    let line = MooringLine::new(150.0, 200.0, 1.0e9, [10.0, 20.0, -50.0]);
    // Fairlead at 45 degrees from the anchor, 125 m away horizontally.
    let offset = 125.0 / 2.0_f64.sqrt();
    let fairlead = [10.0 + offset, 20.0 + offset, 0.0];
    let f = line.restoring_force(fairlead);
    let (t_h, t_v) = line.catenary_tension(125.0, 50.0);
    assert!(t_h > 1.0e3 && t_v > 0.0);
    assert!((f[0] + t_h / 2.0_f64.sqrt()).abs() < 1e-6 * t_h, "{f:?}");
    assert!((f[1] + t_h / 2.0_f64.sqrt()).abs() < 1e-6 * t_h, "{f:?}");
    assert!((f[2] + t_v).abs() < 1e-6 * t_v, "{f:?}");
}

#[test]
fn mooring_tension_handles_negative_height_and_invalid_lines() {
    let line = MooringLine::new(150.0, 200.0, 1.0e9, [0.0, 0.0, 0.0]);
    assert_eq!(
        line.catenary_tension(100.0, -5.0),
        line.catenary_tension(100.0, 0.0)
    );

    let broken = MooringLine::new(0.0, 200.0, 1.0e9, [0.0; 3]);
    assert_eq!(broken.catenary_tension(100.0, 10.0), (0.0, 0.0));
    assert_eq!(broken.restoring_force([100.0, 0.0, 10.0]), [0.0; 3]);

    // An inextensible line that is overstretched has infinite tension but the
    // restoring force must not turn into NaN.
    let taut = MooringLine::new(100.0, 200.0, f64::INFINITY, [0.0; 3]);
    let (t_h, t_v) = taut.catenary_tension(100.0, 0.0);
    assert!(t_h.is_infinite() && t_v.is_infinite());
    let f = taut.restoring_force([0.0, 100.0, 0.0]);
    assert!(f.iter().all(|c| !c.is_nan()), "{f:?}");
    assert_eq!(f[0], 0.0);
    assert!(f[1].is_infinite() && f[1] < 0.0);
}

#[test]
fn mooring_solve_catenary_exposes_the_full_solution() {
    let line = MooringLine::new(300.0, 1000.0, 1.0e8, [0.0, 0.0, 0.0]);
    let solution = line
        .solve_catenary(250.0, 100.0)
        .expect("feasible geometry");
    let (t_h, t_v) = line.catenary_tension(250.0, 100.0);
    assert_eq!(solution.horizontal_tension, t_h);
    assert_eq!(solution.vertical_tension, t_v);
    assert!(solution.grounded_length > 0.0 && solution.grounded_length < 300.0);
    assert!(line.solve_catenary(f64::NAN, 1.0).is_err());
}
