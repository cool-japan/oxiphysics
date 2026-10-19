// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`DoublePendulum`]: Hamilton's equations, energy, and the
//! small-oscillation spectrum, checked independently of each other.
//!
//! Split out of `dynamical_systems.rs` to keep that file under the 2000-line
//! limit.  Included via `#[path]` as `dynamical_systems::double_pendulum_tests`.

use super::*;

const G: f64 = 9.81;

/// Deterministic states spread over the whole configuration space, including
/// momenta of both signs and angle differences near and beyond ±π/2.
fn sample_states() -> Vec<[f64; 4]> {
    (0..64)
        .map(|k| {
            let i = k as f64;
            [
                3.0 * (0.71 * i + 0.3).sin(),
                3.0 * (1.37 * i + 1.1).cos(),
                4.0 * (2.19 * i + 0.5).sin(),
                4.0 * (0.53 * i + 2.0).cos(),
            ]
        })
        .collect()
}

#[test]
fn rhs_is_hamiltons_equations_of_the_energy_function() {
    // q̇ = ∂H/∂p and ṗ = −∂H/∂q, with H differentiated numerically: the old
    // right-hand side (rod-pendulum factors, +g sin θ2) matched no energy function.
    let dp = DoublePendulum::new(G);
    let h = 1e-6;
    for s in sample_states() {
        let rhs = dp.rhs_arr(s);
        let mut gradient = [0.0; 4];
        for (k, slot) in gradient.iter_mut().enumerate() {
            let mut plus = s;
            let mut minus = s;
            plus[k] += h;
            minus[k] -= h;
            *slot = (dp.energy(plus) - dp.energy(minus)) / (2.0 * h);
        }
        let hamilton = [gradient[2], gradient[3], -gradient[0], -gradient[1]];
        for (got, want) in rhs.iter().zip(hamilton) {
            assert!(
                (got - want).abs() < 1e-6 * (1.0 + want.abs()),
                "state {s:?}: rhs {rhs:?} vs Hamilton {hamilton:?}"
            );
        }
    }
}

#[test]
fn energy_matches_the_lagrangian_kinetic_and_potential_energy() {
    // Velocities from the momenta, then T = θ̇1² + θ̇2²/2 + θ̇1 θ̇2 cos Δ for
    // point masses m = l = 1.
    let dp = DoublePendulum::new(G);
    for s in sample_states() {
        let (t1dot, t2dot, _) = DoublePendulum::angular_velocities(s);
        let delta = s[0] - s[1];
        let kinetic = t1dot * t1dot + 0.5 * t2dot * t2dot + t1dot * t2dot * delta.cos();
        let potential = -G * (2.0 * s[0].cos() + s[1].cos());
        assert!((dp.energy(s) - (kinetic + potential)).abs() < 1e-9 * (1.0 + kinetic.abs()));
        // And the momenta are those of the Lagrangian.
        let p1 = 2.0 * t1dot + t2dot * delta.cos();
        let p2 = t2dot + t1dot * delta.cos();
        assert!((p1 - s[2]).abs() < 1e-9 * (1.0 + s[2].abs()));
        assert!((p2 - s[3]).abs() < 1e-9 * (1.0 + s[3].abs()));
    }
    // Hanging at rest: H = V = -3 g.
    assert!((dp.energy([0.0; 4]) + 3.0 * G).abs() < 1e-12);
    // Inverted at rest: H = V = +3 g.
    let inverted = [std::f64::consts::PI, std::f64::consts::PI, 0.0, 0.0];
    assert!((dp.energy(inverted) - 3.0 * G).abs() < 1e-12);
}

#[test]
fn relative_energy_drift_stays_below_1e_8_over_1000_rk4_steps() {
    // The old right-hand side drifted by ~1e-2 (the old test tolerance hid it).
    let dp = DoublePendulum::new(G);
    for s0 in [
        [0.1, 0.05, 0.0, 0.0],
        [1.5, -0.5, 0.0, 0.0],
        [2.5, 1.0, 0.3, -0.2],
        [-1.0, 0.8, 1.0, 0.5],
        [0.5, 2.0, -1.0, 0.7],
    ] {
        let e0 = dp.energy(s0);
        // Relative to the energy scale, which keeps the measure well defined for a
        // start whose total energy happens to be close to zero.
        let scale = e0.abs().max(G);
        let trajectory = dp.integrate(s0, 1e-3, 1000);
        assert_eq!(trajectory.len(), 1001);
        let drift = trajectory
            .iter()
            .map(|s| (dp.energy(*s) - e0).abs() / scale)
            .fold(0.0_f64, f64::max);
        assert!(
            drift < 1e-8,
            "start {s0:?}: relative energy drift {drift:e}"
        );
    }
}

#[test]
fn hanging_equilibrium_is_stable() {
    // A small displacement must stay small (the wrong +g sin θ2 destabilised it).
    let dp = DoublePendulum::new(G);
    let trajectory = dp.integrate([0.01, -0.005, 0.0, 0.0], 1e-3, 5000);
    let largest = trajectory
        .iter()
        .map(|s| s[0].abs().max(s[1].abs()))
        .fold(0.0_f64, f64::max);
    assert!(largest < 0.05, "amplitude grew to {largest}");
}

#[test]
fn small_oscillations_follow_the_normal_modes_of_the_point_mass_pendulum() {
    // Linearised: M = [[2, 1], [1, 1]], K = g diag(2, 1), so ω² = g (2 ∓ √2) with
    // mode shapes θ2/θ1 = ±√2.  Start in a pure mode with amplitude 1e-3 and
    // integrate one period T = 2π/ω.  The residual phase error comes from the
    // amplitude dependence of the frequency (~1e-6 relative), so the state is
    // compared to the linear solution with a tolerance of 1e-4 of its amplitude;
    // a wrong ω (or the old rod-pendulum equations) misses by orders more.
    let dp = DoublePendulum::new(G);
    let sqrt2 = 2.0_f64.sqrt();
    let amplitude = 1e-3;
    for (omega_sq, ratio) in [(G * (2.0 - sqrt2), sqrt2), (G * (2.0 + sqrt2), -sqrt2)] {
        let omega = omega_sq.sqrt();
        let period = 2.0 * std::f64::consts::PI / omega;
        let steps = 4000;
        let start = [amplitude, ratio * amplitude, 0.0, 0.0];
        let trajectory = dp.integrate(start, period / steps as f64, steps);
        // Momentum scale of the mode: the largest |p| reached.
        let p_scale = trajectory
            .iter()
            .map(|s| s[2].abs().max(s[3].abs()))
            .fold(0.0_f64, f64::max);
        assert!(p_scale > 2.0 * amplitude, "p_scale = {p_scale}");
        let end = trajectory[steps];
        for (k, (got, want)) in end.iter().zip(start).enumerate() {
            let scale = if k < 2 { amplitude } else { p_scale };
            assert!(
                (got - want).abs() < 1e-4 * scale,
                "ω² = {omega_sq}: end {end:?} vs start {start:?}"
            );
        }
        // Half a period later the mode has flipped sign; a quarter period later
        // the angles pass through zero.
        let half = trajectory[steps / 2];
        assert!((half[0] + start[0]).abs() < 1e-4 * amplitude);
        assert!((half[1] + start[1]).abs() < 1e-4 * amplitude);
        let quarter = trajectory[steps / 4];
        assert!(quarter[0].abs() < 1e-4 * amplitude && quarter[1].abs() < 1e-4 * amplitude);
    }
}

#[test]
fn velocities_invert_the_momenta_exactly() {
    for s in sample_states() {
        let (t1dot, t2dot, sin_delta) = DoublePendulum::angular_velocities(s);
        let cos_delta = (s[0] - s[1]).cos();
        assert!((2.0 * t1dot + cos_delta * t2dot - s[2]).abs() < 1e-9 * (1.0 + s[2].abs()));
        assert!((cos_delta * t1dot + t2dot - s[3]).abs() < 1e-9 * (1.0 + s[3].abs()));
        assert!((sin_delta - (s[0] - s[1]).sin()).abs() < 1e-15);
    }
}
