// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for the 2-D constant-strength source panel method.
//!
//! `Panel::source_influence` used a made-up kernel built on the source panel's
//! own normal.  The real influence of a constant-strength source sheet is
//! `u = ln(r1 / r2) / 2π` along the panel and `w = (θ2 − θ1) / 2π` across it.

use oxiphysics_rigid::Error;
use oxiphysics_rigid::fluid_dynamics::{Panel, SourcePanelMethod};
use std::f64::consts::PI;

const U_INF: f64 = 10.0;

fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

/// Panels through the points of a closed curve, ordered as given.
fn closed_polygon(points: &[[f64; 2]]) -> Vec<Panel> {
    (0..points.len())
        .map(|i| Panel::new(points[i], points[(i + 1) % points.len()]))
        .collect()
}

fn ellipse(a: f64, b: f64, n: usize) -> Vec<Panel> {
    let points: Vec<[f64; 2]> = (0..n)
        .map(|i| {
            let t = 2.0 * PI * i as f64 / n as f64;
            [a * t.cos(), b * t.sin()]
        })
        .collect();
    closed_polygon(&points)
}

// ── single panel ────────────────────────────────────────────────────────────

#[test]
fn far_field_of_a_panel_is_that_of_a_point_source() {
    // The panel runs along +x, so its normal is -y.
    let panel = Panel::new([0.0, 0.0], [1.0, 0.0]);
    assert_eq!(panel.normal(), [0.0, -1.0]);
    let expected = 1.0 / (2.0 * PI * 100.0);

    let below = panel.source_velocity([0.5, -100.0]);
    assert!(below[0].abs() < 1.0e-12);
    assert!((below[1] + expected).abs() < 1.0e-4 * expected, "{below:?}");

    let above = panel.source_velocity([0.5, 100.0]);
    assert!((above[1] - expected).abs() < 1.0e-4 * expected, "{above:?}");

    // Radial flow away from the panel in an oblique direction as well.
    let far = panel.source_velocity([70.0, 70.0]);
    // Distance from the panel centre (0.5, 0).
    let r = (69.5_f64 * 69.5 + 70.0 * 70.0).sqrt();
    let speed = (far[0] * far[0] + far[1] * far[1]).sqrt();
    assert!((speed - 1.0 / (2.0 * PI * r)).abs() < 1.0e-3 * speed);
    assert!(far[0] > 0.0 && far[1] > 0.0);
}

#[test]
fn normal_velocity_jumps_by_the_source_strength_across_the_sheet() {
    let panel = Panel::new([0.0, 0.0], [1.0, 0.0]);
    let n = panel.normal();
    let centre = panel.center();
    for eps in [1.0e-3, 1.0e-6, 1.0e-9] {
        let outside = [centre[0] + eps * n[0], centre[1] + eps * n[1]];
        let inside = [centre[0] - eps * n[0], centre[1] - eps * n[1]];
        let vn_out = dot(panel.source_velocity(outside), n);
        let vn_in = dot(panel.source_velocity(inside), n);
        assert!((vn_out - 0.5).abs() < 1.0e-3, "eps {eps}: outside {vn_out}");
        assert!((vn_in + 0.5).abs() < 1.0e-3, "eps {eps}: inside {vn_in}");
        assert!((vn_out + vn_in).abs() < 1.0e-9, "the jump is symmetric");
    }
    // On the panel itself the limit from the normal side is returned.
    let on_panel = panel.source_velocity(centre);
    assert!((dot(on_panel, n) - 0.5).abs() < 1.0e-12, "{on_panel:?}");
    assert!(dot(on_panel, panel.tangent()).abs() < 1.0e-12);
}

#[test]
fn normal_jump_follows_the_panel_orientation() {
    // A panel rotated by 135 degrees: the same jump across its own normal.
    let panel = Panel::new([1.0, 1.0], [0.0, 2.0]);
    let n = panel.normal();
    let centre = panel.center();
    let outside = [centre[0] + 1.0e-8 * n[0], centre[1] + 1.0e-8 * n[1]];
    let vn = dot(panel.source_velocity(outside), n);
    assert!((vn - 0.5).abs() < 1.0e-6, "{vn}");
}

#[test]
fn velocity_along_the_panel_line_beyond_its_ends_points_away() {
    let panel = Panel::new([0.0, 0.0], [1.0, 0.0]);
    let beyond_end = panel.source_velocity([2.0, 0.0]);
    let expected = 2.0_f64.ln() / (2.0 * PI);
    assert!((beyond_end[0] - expected).abs() < 1.0e-12, "{beyond_end:?}");
    assert!(beyond_end[1].abs() < 1.0e-12);
    let before_start = panel.source_velocity([-1.0, 0.0]);
    assert!(
        (before_start[0] + expected).abs() < 1.0e-12,
        "{before_start:?}"
    );
    assert!(before_start[1].abs() < 1.0e-12);
}

#[test]
fn panel_velocity_matches_quadrature_of_point_sources() {
    let panel = Panel::new([1.0, 0.5], [-0.5, 1.7]);
    let length = panel.length();
    let tangent = panel.tangent();
    let points = [
        [0.0, 0.0],
        [2.0, 3.0],
        [-1.0, -1.0],
        [0.3, 1.6],
        [-2.0, 2.5],
        [1.5, -0.7],
    ];
    let samples = 20_000;
    for p in points {
        let mut v = [0.0_f64; 2];
        for k in 0..samples {
            let s = (k as f64 + 0.5) * length / samples as f64;
            let source = [panel.p1[0] + s * tangent[0], panel.p1[1] + s * tangent[1]];
            let d = [p[0] - source[0], p[1] - source[1]];
            let r2 = dot(d, d);
            v[0] += d[0] / r2;
            v[1] += d[1] / r2;
        }
        let scale = length / samples as f64 / (2.0 * PI);
        v = [v[0] * scale, v[1] * scale];
        let analytic = panel.source_velocity(p);
        assert!(
            (analytic[0] - v[0]).abs() < 1.0e-7 && (analytic[1] - v[1]).abs() < 1.0e-7,
            "point {p:?}: analytic {analytic:?}, quadrature {v:?}"
        );
    }
}

#[test]
fn degenerate_panel_induces_no_velocity() {
    let panel = Panel::new([1.0, 1.0], [1.0, 1.0]);
    assert_eq!(panel.source_velocity([3.0, -2.0]), [0.0, 0.0]);
}

#[test]
fn field_point_at_a_panel_end_stays_finite() {
    let panel = Panel::new([0.0, 0.0], [1.0, 0.0]);
    for p in [[0.0, 0.0], [1.0, 0.0]] {
        let v = panel.source_velocity(p);
        assert!(v[0].is_finite() && v[1].is_finite(), "{p:?} -> {v:?}");
    }
}

// ── influence matrices ──────────────────────────────────────────────────────

#[test]
fn influence_matrices_have_the_expected_diagonals_and_symmetry() {
    let method = SourcePanelMethod::circle(1.0, 16, U_INF, 0.0);
    let a = method.influence_matrix();
    let b = method.tangential_influence_matrix();
    assert_eq!(a.len(), 16);
    for i in 0..16 {
        assert_eq!(a[i].len(), 16);
        assert_eq!(a[i][i], 0.5, "normal self-influence");
        assert_eq!(b[i][i], 0.0, "tangential self-influence");
        for j in 0..16 {
            // Rotating a regular polygon onto itself leaves the matrix alone.
            let (i2, j2) = ((i + 1) % 16, (j + 1) % 16);
            assert!((a[i][j] - a[i2][j2]).abs() < 1.0e-12);
            assert!((b[i][j] - b[i2][j2]).abs() < 1.0e-12);
        }
    }
    // A neighbouring panel pushes fluid towards the normal side of the next.
    assert!(a[0][1] != 0.0 && a[0][1].abs() < 0.5);
}

#[test]
fn influence_entries_come_from_the_source_velocity_of_each_panel() {
    let panels = ellipse(2.0, 1.0, 12);
    let method = SourcePanelMethod::new(panels.clone(), U_INF, 0.3);
    let a = method.influence_matrix();
    let b = method.tangential_influence_matrix();
    for (i, pi) in panels.iter().enumerate() {
        for (j, pj) in panels.iter().enumerate() {
            if i == j {
                continue;
            }
            let v = pj.source_velocity(pi.center());
            assert!((a[i][j] - dot(v, pi.normal())).abs() < 1.0e-12);
            assert!((b[i][j] - dot(v, pi.tangent())).abs() < 1.0e-12);
        }
    }
}

// ── solved flows ────────────────────────────────────────────────────────────

#[test]
fn circle_pressure_matches_potential_flow() {
    let n = 64;
    let method = SourcePanelMethod::circle(1.0, n, U_INF, 0.0);
    let solution = method.solve().expect("solvable");
    assert_eq!(solution.pressure_coefficient.len(), n);
    for (i, cp) in solution.pressure_coefficient.iter().enumerate() {
        let theta = 2.0 * PI * (i as f64 + 0.5) / n as f64;
        let exact = 1.0 - 4.0 * theta.sin().powi(2);
        assert!(
            (cp - exact).abs() < 0.02,
            "panel {i}: Cp {cp}, exact {exact}"
        );
        // Counter-clockwise tangent: the surface speed is -2 U sin(theta).
        let vt = solution.tangential_velocity[i];
        assert!(
            (vt + 2.0 * U_INF * theta.sin()).abs() < 0.02 * U_INF,
            "panel {i}: V_t {vt}"
        );
    }
    // A closed body neither creates nor destroys fluid.
    assert!(
        solution.net_source.abs() < 1.0e-9 * U_INF,
        "net source {}",
        solution.net_source
    );
}

#[test]
fn stagnation_points_sit_on_the_axis_of_the_stream() {
    let n = 64;
    let solution = SourcePanelMethod::circle(1.0, n, U_INF, 0.0)
        .solve()
        .expect("solvable");
    let cp = &solution.pressure_coefficient;
    let front = cp[n - 1].max(cp[0]);
    let back = cp[n / 2 - 1].max(cp[n / 2]);
    assert!(front > 0.97 && front <= 1.0 + 1.0e-6, "front Cp {front}");
    assert!(back > 0.97 && back <= 1.0 + 1.0e-6, "rear Cp {back}");
    let top = cp[n / 4 - 1].min(cp[n / 4]);
    assert!(top < -2.9 && top > -3.1, "shoulder Cp {top}");
}

#[test]
fn rotating_the_stream_rotates_the_pressure_pattern() {
    let n = 64;
    let alpha = 0.7;
    let solution = SourcePanelMethod::circle(1.0, n, U_INF, alpha)
        .solve()
        .expect("solvable");
    for (i, cp) in solution.pressure_coefficient.iter().enumerate() {
        let theta = 2.0 * PI * (i as f64 + 0.5) / n as f64;
        let exact = 1.0 - 4.0 * (theta - alpha).sin().powi(2);
        assert!(
            (cp - exact).abs() < 0.02,
            "panel {i}: Cp {cp}, exact {exact}"
        );
    }
}

#[test]
fn ellipse_surface_speed_matches_the_exact_solution() {
    // Flow along the major axis of an ellipse with semi-axes (a, b): the surface
    // speed at parameter t is U (a + b) |sin t| / sqrt(a^2 sin^2 t + b^2 cos^2 t).
    let (a, b, n) = (2.0, 1.0, 96);
    let method = SourcePanelMethod::new(ellipse(a, b, n), U_INF, 0.0);
    let solution = method.solve().expect("solvable");
    for (i, cp) in solution.pressure_coefficient.iter().enumerate() {
        let t = 2.0 * PI * (i as f64 + 0.5) / n as f64;
        let speed = U_INF * (a + b) * t.sin().abs()
            / (a * a * t.sin().powi(2) + b * b * t.cos().powi(2)).sqrt();
        let exact = 1.0 - (speed / U_INF).powi(2);
        assert!(
            (cp - exact).abs() < 0.01,
            "panel {i}: Cp {cp}, exact {exact}"
        );
    }
    assert!(solution.net_source.abs() < 1.0e-9 * U_INF);
}

#[test]
fn field_around_a_solved_circle_matches_the_exact_flow() {
    let n = 128;
    let method = SourcePanelMethod::circle(1.0, n, U_INF, 0.0);
    let solution = method.solve().expect("solvable");
    for p in [
        [2.0, 0.0],
        [0.0, 1.5],
        [-1.7, 1.2],
        [3.0, -2.0],
        [-4.0, -0.5],
    ] {
        let r2 = dot(p, p);
        let phi = p[1].atan2(p[0]);
        let exact = [
            U_INF * (1.0 - (2.0 * phi).cos() / r2),
            -U_INF * (2.0 * phi).sin() / r2,
        ];
        let v = method
            .velocity_at(&solution.strengths, p)
            .expect("one strength per panel");
        assert!(
            (v[0] - exact[0]).abs() < 0.01 * U_INF && (v[1] - exact[1]).abs() < 0.01 * U_INF,
            "point {p:?}: panels {v:?}, exact {exact:?}"
        );
    }
}

#[test]
fn no_fluid_crosses_the_body_surface() {
    let panels = ellipse(2.0, 1.0, 48);
    let method = SourcePanelMethod::new(panels.clone(), U_INF, 0.4);
    let solution = method.solve().expect("solvable");
    for (i, panel) in panels.iter().enumerate() {
        let n = panel.normal();
        let c = panel.center();
        let just_outside = [c[0] + 1.0e-9 * n[0], c[1] + 1.0e-9 * n[1]];
        let v = method
            .velocity_at(&solution.strengths, just_outside)
            .expect("one strength per panel");
        assert!(
            dot(v, n).abs() < 1.0e-6 * U_INF,
            "panel {i}: normal velocity {}",
            dot(v, n)
        );
    }
}

#[test]
fn still_fluid_gives_zero_strengths() {
    let solution = SourcePanelMethod::circle(1.0, 16, 0.0, 0.0)
        .solve()
        .expect("solvable");
    assert!(solution.strengths.iter().all(|s| s.abs() < 1.0e-15));
    assert!(solution.pressure_coefficient.iter().all(|c| *c == 0.0));
}

// ── error handling ──────────────────────────────────────────────────────────

#[test]
fn degenerate_inputs_are_reported() {
    let empty = SourcePanelMethod::new(Vec::new(), U_INF, 0.0);
    assert!(matches!(empty.solve(), Err(Error::General(_))));

    let mut panels = ellipse(2.0, 1.0, 8);
    panels[3] = Panel::new([0.5, 0.5], [0.5, 0.5]);
    let zero_length = SourcePanelMethod::new(panels, U_INF, 0.0);
    assert!(matches!(zero_length.solve(), Err(Error::General(_))));

    let mut panels = ellipse(2.0, 1.0, 8);
    panels[2] = Panel::new([f64::NAN, 0.0], [1.0, 0.0]);
    let not_finite = SourcePanelMethod::new(panels, U_INF, 0.0);
    assert!(matches!(not_finite.solve(), Err(Error::General(_))));
}

#[test]
fn coincident_panels_make_the_system_singular() {
    let mut panels = ellipse(2.0, 1.0, 8);
    panels.push(panels[0].clone());
    let method = SourcePanelMethod::new(panels, U_INF, 0.0);
    assert!(matches!(method.solve(), Err(Error::General(_))));
}

#[test]
fn velocity_field_needs_one_strength_per_panel() {
    let method = SourcePanelMethod::circle(1.0, 8, U_INF, 0.0);
    assert!(matches!(
        method.velocity_at(&[0.0; 7], [2.0, 0.0]),
        Err(Error::General(_))
    ));
    let v = method
        .velocity_at(&[0.0; 8], [2.0, 0.0])
        .expect("matching length");
    assert!((v[0] - U_INF).abs() < 1.0e-12 && v[1].abs() < 1.0e-12);
}
