// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for the CFD-coupling structural models: the
//! vortex-induced-vibration oscillator and the catenary mooring line.
//!
//! Note: `oxiphysics_rigid::VortexInducedVibration` (crate root) is the
//! fluid-coupling type; the wake-oscillator model under test lives in
//! `oxiphysics_rigid::cfd_coupling`.

use oxiphysics_rigid::cfd_coupling::{MooringLine, VortexInducedVibration};
use std::f64::consts::PI;

/// Free-decay model: no fluid, no flow, so the only forces are the structural
/// damping and stiffness terms.
fn free_decay_model(mass_per_length: f64, wn: f64, zeta: f64, y0: f64) -> VortexInducedVibration {
    let mut viv = VortexInducedVibration::new(0.5, 10.0, mass_per_length, wn, zeta, 0.0, 0.0);
    viv.displacement = y0;
    viv
}

/// Times of the downward zero crossings (`y` going from positive to negative)
/// of the free decay, linearly interpolated inside the step.
fn downward_crossings(viv: &mut VortexInducedVibration, dt: f64, t_end: f64) -> Vec<f64> {
    let mut crossings = Vec::new();
    let mut previous_y = viv.displacement;
    let mut t = 0.0;
    while t < t_end {
        let (y, _q) = viv.step(dt, 0.0, 0.0);
        if previous_y > 0.0 && y <= 0.0 {
            let fraction = previous_y / (previous_y - y);
            crossings.push(t + fraction * dt);
        }
        previous_y = y;
        t += dt;
    }
    crossings
}

#[test]
fn free_decay_period_matches_the_damped_natural_period() {
    let (wn, zeta) = (10.0_f64, 0.05_f64);
    let mut viv = free_decay_model(50.0, wn, zeta, 0.01);
    let crossings = downward_crossings(&mut viv, 1.0e-4, 4.0);
    assert!(
        crossings.len() >= 4,
        "expected several oscillation cycles, got {crossings:?}"
    );
    // Same-sign crossings are one damped period apart; skip the first
    // half-cycle (the initial condition is not a crossing).
    let measured = crossings[2] - crossings[1];
    let expected = 2.0 * PI / (wn * (1.0 - zeta * zeta).sqrt());
    assert!(
        (measured - expected).abs() < 0.01 * expected,
        "period {measured} s vs expected {expected} s"
    );
}

#[test]
fn free_decay_period_does_not_depend_on_the_mass_per_length() {
    let (wn, zeta) = (6.0_f64, 0.02_f64);
    let expected = 2.0 * PI / (wn * (1.0 - zeta * zeta).sqrt());
    for mass in [0.5, 50.0, 5000.0] {
        let mut viv = free_decay_model(mass, wn, zeta, 0.02);
        let crossings = downward_crossings(&mut viv, 1.0e-4, 4.0);
        assert!(crossings.len() >= 3, "m={mass}: crossings {crossings:?}");
        let measured = crossings[2] - crossings[1];
        assert!(
            (measured - expected).abs() < 0.01 * expected,
            "m={mass}: period {measured} s vs expected {expected} s"
        );
    }
}

#[test]
fn free_decay_amplitude_follows_the_damping_ratio() {
    // Logarithmic decrement over one period: ln(y_n / y_{n+1}) = 2 pi zeta / sqrt(1 - zeta^2).
    let (wn, zeta) = (10.0_f64, 0.1_f64);
    let mut viv = free_decay_model(100.0, wn, zeta, 0.05);
    let dt = 1.0e-4;
    let mut peaks = Vec::new();
    let (mut y_prev2, mut y_prev) = (viv.displacement, viv.displacement);
    for _ in 0..40_000 {
        let (y, _q) = viv.step(dt, 0.0, 0.0);
        if y_prev > y_prev2 && y_prev >= y && y_prev > 0.0 {
            peaks.push(y_prev);
        }
        y_prev2 = y_prev;
        y_prev = y;
    }
    assert!(peaks.len() >= 3, "peaks: {peaks:?}");
    let decrement = (peaks[0] / peaks[1]).ln();
    let expected = 2.0 * PI * zeta / (1.0 - zeta * zeta).sqrt();
    assert!(
        (decrement - expected).abs() < 0.02 * expected,
        "log decrement {decrement} vs expected {expected}"
    );
}

#[test]
fn constant_load_settles_at_force_over_mass_times_omega_squared() {
    // Constant fluid acceleration, inertia force only: the static deflection is
    // F / (m wn^2) - the stiffness term must not be divided by the mass again.
    let (mass, wn, zeta) = (80.0_f64, 5.0_f64, 0.5_f64);
    let fluid_density = 1025.0;
    let mut viv = VortexInducedVibration::new(0.4, 10.0, mass, wn, zeta, 0.0, fluid_density);
    viv.cd_morison = 0.0;
    let acceleration = 0.3;
    let force_per_length = fluid_density * viv.cm * PI * 0.4 * 0.4 / 4.0 * acceleration;
    let expected = force_per_length / (mass * wn * wn);
    for _ in 0..200_000 {
        viv.step(1.0e-4, 0.0, acceleration);
    }
    assert!(
        (viv.displacement - expected).abs() < 1e-4 * expected,
        "static deflection {} vs expected {expected}",
        viv.displacement
    );
    assert!(viv.disp_velocity.abs() < 1e-6);
}

#[test]
fn massless_structure_keeps_the_restoring_terms() {
    // With no mass the force term is dropped, but damping and stiffness still
    // act on the displacement (no division by the vanishing mass).
    let mut viv = free_decay_model(0.0, 4.0, 0.1, 0.01);
    let (y, q) = viv.step(1.0e-3, 0.0, 1.0);
    assert!(y.is_finite() && q.is_finite());
    assert!(
        viv.disp_velocity < 0.0,
        "restoring acceleration points back"
    );
}

// ── MooringLine ─────────────────────────────────────────────────────────────

const WEIGHT: f64 = 1000.0;

fn line(length: f64, anchor: [f64; 3], fairlead: [f64; 3]) -> MooringLine {
    MooringLine::new(length, WEIGHT, 1.0e9, anchor, fairlead, 40)
}

/// Analytic arc length of the solved catenary between anchor and fairlead.
fn arc_length(ml: &MooringLine) -> f64 {
    let a = ml.catenary_parameter();
    let x0 = ml.vertex_offset;
    let h = ml.horizontal_span();
    a * (((h - x0) / a).sinh() + (x0 / a).sinh())
}

/// Assert every defining property of a solved catenary.
fn assert_consistent_catenary(ml: &MooringLine) {
    let l = ml.unstretched_length;
    let h = ml.horizontal_span();
    let v = ml.vertical_rise();
    let a = ml.catenary_parameter();
    let x0 = ml.vertex_offset;
    // The profile starts at the anchor and reaches the fairlead height.
    assert!(
        ml.profile_z(0.0).abs() < 1e-9,
        "z(0) = {}",
        ml.profile_z(0.0)
    );
    assert!(
        (ml.profile_z(h) - v).abs() < 1e-9 * (1.0 + v.abs()),
        "z(h) = {} vs rise {v}",
        ml.profile_z(h)
    );
    // The arc length between the end points is the line length.
    assert!(
        (arc_length(ml) - l).abs() < 1e-9 * l,
        "arc length {} vs line length {l}",
        arc_length(ml)
    );
    // 2a sinh(h/2a) = sqrt(L^2 - v^2).
    let lhs = 2.0 * a * (h / (2.0 * a)).sinh();
    let rhs = (l * l - v * v).sqrt();
    assert!((lhs - rhs).abs() < 1e-9 * rhs, "{lhs} vs {rhs}");
    // Vertical force balance: V_fairlead - V_anchor = w L, with
    // V_fairlead = w a sinh((h - x0)/a) and V_anchor = -w a sinh(x0/a).
    let w = ml.weight_per_length;
    let v_fair = w * a * ((h - x0) / a).sinh();
    let v_anchor = -w * a * (x0 / a).sinh();
    assert!(
        (v_fair - v_anchor - w * l).abs() < 1e-9 * w * l,
        "weight balance: {} vs {}",
        v_fair - v_anchor,
        w * l
    );
    assert!(
        (ml.vertical_tension - v_fair).abs() < 1e-8 * v_fair.abs().max(1.0),
        "T_v = {} vs w a sinh((h-x0)/a) = {v_fair}",
        ml.vertical_tension
    );
    assert!(ml.fairlead_tension() >= ml.horizontal_tension);
}

#[test]
fn solved_catenary_reaches_the_fairlead_with_the_line_length() {
    // L = 120 m between a fairlead 100 m away and 20 m above the anchor.
    let mut ml = line(120.0, [0.0; 3], [100.0, 20.0, 0.0]);
    let t_h = ml.solve_static().expect("slack line has a catenary");
    assert!(t_h > 0.0 && (t_h - ml.horizontal_tension).abs() < f64::EPSILON);
    assert_consistent_catenary(&ml);
}

#[test]
fn symmetric_catenary_has_its_vertex_at_midspan() {
    let mut ml = line(130.0, [0.0; 3], [100.0, 0.0, 0.0]);
    ml.solve_static().expect("slack line has a catenary");
    assert!((ml.vertex_offset - 50.0).abs() < 1e-9);
    assert!((ml.vertical_tension - 0.5 * WEIGHT * 130.0).abs() < 1e-6);
    assert_consistent_catenary(&ml);
    // The lowest point of the symmetric curve sags below the end points.
    assert!(ml.profile_z(50.0) < 0.0);
}

#[test]
fn catenary_with_fairlead_below_the_anchor() {
    let mut ml = line(140.0, [5.0, 10.0, -3.0], [105.0, -20.0, -3.0]);
    ml.solve_static().expect("slack line has a catenary");
    assert!((ml.vertical_rise() + 30.0).abs() < 1e-12);
    assert_consistent_catenary(&ml);
}

#[test]
fn deep_sag_and_nearly_taut_lines_are_solved_accurately() {
    // Very long line over a short span: tiny catenary parameter.
    let mut slack = line(1000.0, [0.0; 3], [10.0, 0.0, 0.0]);
    slack.solve_static().expect("slack line has a catenary");
    assert!(slack.catenary_parameter() < 1.0);
    assert_consistent_catenary(&slack);
    // A hair longer than the chord: huge catenary parameter, tiny sag.
    let mut tight = line(100.001, [0.0; 3], [100.0, 0.0, 0.0]);
    tight
        .solve_static()
        .expect("slightly slack line has a catenary");
    assert!(tight.catenary_parameter() > 1.0e3);
    assert_consistent_catenary(&tight);
    assert!(tight.horizontal_tension > slack.horizontal_tension);
}

#[test]
fn catenary_shape_does_not_depend_on_the_weight() {
    let mut light = MooringLine::new(120.0, 50.0, 1.0e9, [0.0; 3], [100.0, 20.0, 0.0], 10);
    let mut heavy = MooringLine::new(120.0, 5000.0, 1.0e9, [0.0; 3], [100.0, 20.0, 0.0], 10);
    light.solve_static().expect("slack line");
    heavy.solve_static().expect("slack line");
    assert!((light.catenary_parameter() - heavy.catenary_parameter()).abs() < 1e-9);
    assert!((light.vertex_offset - heavy.vertex_offset).abs() < 1e-9);
    assert!((heavy.horizontal_tension / light.horizontal_tension - 100.0).abs() < 1e-9);
}

#[test]
fn taut_line_is_an_error_and_leaves_the_state_untouched() {
    // Chord = sqrt(100^2 + 20^2) = 101.98 m > 100 m of line.
    let mut ml = line(100.0, [0.0; 3], [100.0, 20.0, 0.0]);
    ml.horizontal_tension = 123.0;
    ml.vertical_tension = 45.0;
    ml.vertex_offset = 6.0;
    assert!(ml.solve_static().is_err());
    // Line exactly as long as the chord: still taut.
    let mut exact = line(5.0, [0.0; 3], [3.0, 4.0, 0.0]);
    assert!(exact.solve_static().is_err());
    assert_eq!(ml.horizontal_tension, 123.0);
    assert_eq!(ml.vertical_tension, 45.0);
    assert_eq!(ml.vertex_offset, 6.0);
}

#[test]
fn degenerate_mooring_inputs_are_errors() {
    let mut weightless = MooringLine::new(120.0, 0.0, 1.0e9, [0.0; 3], [100.0, 20.0, 0.0], 8);
    assert!(weightless.solve_static().is_err());
    let mut negative = MooringLine::new(120.0, -5.0, 1.0e9, [0.0; 3], [100.0, 20.0, 0.0], 8);
    assert!(negative.solve_static().is_err());
    let mut no_length = line(0.0, [0.0; 3], [100.0, 20.0, 0.0]);
    assert!(no_length.solve_static().is_err());
    let mut nan_length = line(f64::NAN, [0.0; 3], [100.0, 20.0, 0.0]);
    assert!(nan_length.solve_static().is_err());
    // Fairlead straight above the anchor: no horizontal span.
    let mut vertical = line(120.0, [0.0; 3], [0.0, 50.0, 0.0]);
    assert!(vertical.solve_static().is_err());
}

#[test]
fn discretised_line_runs_from_anchor_to_fairlead_along_the_plane() {
    // Horizontal span 100 m in direction (0.6, 0.8) of the x-z plane.
    let anchor = [10.0, -50.0, 20.0];
    let fairlead = [70.0, -30.0, 100.0];
    let mut ml = MooringLine::new(120.0, WEIGHT, 1.0e9, anchor, fairlead, 2000);
    ml.solve_static().expect("slack line has a catenary");
    let pts = ml.discretised_positions();
    assert_eq!(pts.len(), 2001);
    for k in 0..3 {
        assert!((pts[0][k] - anchor[k]).abs() < 1e-12);
        assert!(
            (pts[2000][k] - fairlead[k]).abs() < 1e-9,
            "last point {:?} vs fairlead {fairlead:?}",
            pts[2000]
        );
    }
    for p in &pts {
        // Every point lies in the vertical plane through anchor and fairlead.
        let s = ((p[0] - anchor[0]) * 0.6 + (p[2] - anchor[2]) * 0.8).abs();
        let off_plane = ((p[0] - anchor[0]) * 0.8 - (p[2] - anchor[2]) * 0.6).abs();
        assert!(
            off_plane < 1e-9,
            "point {p:?} is off the plane by {off_plane}"
        );
        assert!((p[1] - anchor[1] - ml.profile_z(s)).abs() < 1e-9);
    }
    // The polyline approximates the arc length L.
    let polyline: f64 = pts
        .windows(2)
        .map(|w| {
            let d = [w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        })
        .sum();
    assert!(
        (polyline - 120.0).abs() < 1e-3,
        "polyline length {polyline}"
    );
}

#[test]
fn discretisation_of_a_line_without_segments_is_finite() {
    let mut ml = MooringLine::new(120.0, WEIGHT, 1.0e9, [0.0; 3], [100.0, 20.0, 0.0], 0);
    ml.solve_static().expect("slack line has a catenary");
    let pts = ml.discretised_positions();
    assert_eq!(pts.len(), 2);
    assert!(pts.iter().all(|p| p.iter().all(|c| c.is_finite())));
}

#[test]
fn linearised_stiffness_is_the_derivative_of_the_horizontal_tension() {
    let solve_at = |span: f64| {
        let mut ml = line(150.0, [0.0; 3], [span, 25.0, 0.0]);
        ml.solve_static().expect("slack line has a catenary")
    };
    for span in [60.0, 100.0, 140.0] {
        let mut ml = line(150.0, [0.0; 3], [span, 25.0, 0.0]);
        ml.solve_static().expect("slack line has a catenary");
        let delta = 1.0e-2;
        let numeric = (solve_at(span + delta) - solve_at(span - delta)) / (2.0 * delta);
        let analytic = ml.linearised_stiffness();
        assert!(analytic > 0.0);
        assert!(
            (analytic - numeric).abs() < 1e-5 * numeric.abs(),
            "span {span}: analytic {analytic} vs finite difference {numeric}"
        );
    }
    // Stiffening: pulling the fairlead further away raises the stiffness.
    let mut near = line(150.0, [0.0; 3], [80.0, 25.0, 0.0]);
    let mut far = line(150.0, [0.0; 3], [140.0, 25.0, 0.0]);
    near.solve_static().expect("slack line");
    far.solve_static().expect("slack line");
    assert!(far.linearised_stiffness() > near.linearised_stiffness());
    // Nothing solved yet: no stiffness.
    assert_eq!(
        line(150.0, [0.0; 3], [80.0, 25.0, 0.0]).linearised_stiffness(),
        0.0
    );
}

#[test]
fn profile_with_a_prescribed_tension_keeps_the_vertex_at_the_anchor() {
    let mut ml = line(120.0, [0.0; 3], [80.0, 50.0, 0.0]);
    ml.horizontal_tension = 1.0e5;
    let a = ml.catenary_parameter();
    for x in [0.0, 10.0, 25.0, 60.0] {
        let expected = a * ((x / a).cosh() - 1.0);
        assert!((ml.profile_z(x) - expected).abs() < 1e-12);
    }
}
