// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Wheel-spin integration helpers of the raycast vehicle.
//!
//! The spin `ω` of a wheel with inertia `I` and radius `r` obeys
//!
//! ```text
//! I · dω/dt = T_drive − r · Fx(κ(ω)) − T_brake,    κ = (ω·r − v) / max(|v|, v_min)
//! ```
//!
//! where the longitudinal tire force `Fx` grows with the slip ratio `κ`.  Near
//! zero slip a tire is stiff (`dFx/dκ` of the order of 10⁵ N) while a wheel is
//! light, so the spin relaxes towards the rolling condition `ω·r = v` with a
//! rate `r²·(dFx/dκ) / (I·max(|v|, v_min))` of several hundred to tens of
//! thousands per second.  A plain explicit Euler step is only stable when that
//! rate times the step stays below 2; beyond it the spin flips sign every frame
//! and the tire force alternates between ±μ·Fz.
//!
//! [`advance_wheel_spin`] therefore takes a *linearly implicit* (Rosenbrock–
//! Euler) step that treats the spin dependence of the tire torque implicitly:
//! it is stable for any step size, converges monotonically to the rolling
//! condition and reduces to explicit Euler as the stiffness vanishes.  The
//! brake is applied afterwards as dry friction with the same effective inertia,
//! so it can stop a wheel but never reverse it, and the equilibria of the step
//! are the physical ones.

use crate::tire::TireModel;
use oxiphysics_core::math::Real;

/// Slip-ratio increment of the forward difference that estimates `dFx/dκ`.
const SLIP_STIFFNESS_PROBE: Real = 1.0e-4;

/// Longitudinal slip stiffness `∂Fx/∂κ` of a tire at its operating point, in N
/// per unit slip ratio.
///
/// It is the forward difference of [`TireModel::compute_forces`] over a small
/// slip-ratio increment with slip angle, load and friction held fixed; `fx` is
/// the longitudinal force already evaluated at `slip_ratio`.  A negative slope
/// (a tire past its force peak) or a non-finite estimate is reported as `0.0`,
/// so the result is always a valid, purely dissipative stiffness.
pub(crate) fn longitudinal_slip_stiffness<T: TireModel>(
    tire_model: &T,
    slip_angle: Real,
    slip_ratio: Real,
    normal_force: Real,
    friction: Real,
    fx: Real,
) -> Real {
    let (_, fx_probe) = tire_model.compute_forces(
        slip_angle,
        slip_ratio + SLIP_STIFFNESS_PROBE,
        normal_force,
        friction,
    );
    let slope = (fx_probe - fx) / SLIP_STIFFNESS_PROBE;
    if slope.is_finite() {
        slope.max(0.0)
    } else {
        0.0
    }
}

/// Apply a brake torque to a wheel spin speed as dry (Coulomb) friction.
///
/// A brake dissipates energy: it removes up to `brake_torque / inertia * dt`
/// from `|spin|` and can bring the wheel to rest, but it never reverses it.
/// A wheel whose speed is already below that amount (including a stationary
/// wheel) is held at exactly zero, which models static brake friction.
/// Non-positive or NaN brake torque leaves the spin untouched.
pub(crate) fn apply_brake_torque(spin: Real, brake_torque: Real, inertia: Real, dt: Real) -> Real {
    if brake_torque > 0.0 {
        let max_change = brake_torque / inertia * dt;
        if spin.abs() <= max_change {
            0.0
        } else {
            spin - max_change * spin.signum()
        }
    } else {
        spin
    }
}

/// Advance the spin of a wheel by one step of length `dt`.
///
/// * `spin` — current spin speed (rad/s, signed).
/// * `net_torque` — drive torque minus the tire reaction `r · Fx`, evaluated at
///   the start of the step (N·m); the brake is *not* part of it.
/// * `torque_stiffness` — `−∂(net torque)/∂ω` (N·m·s), i.e. how strongly the
///   tire reaction torque grows with the spin; `r² · (dFx/dκ) / max(|v|, v_min)`
///   for the slip-ratio definition of the vehicle.  Negative and NaN values are
///   treated as zero (no implicit damping).
/// * `brake_torque` — brake torque magnitude on this wheel (N·m, `>= 0`).
/// * `inertia` — wheel inertia about its axle (kg·m², `> 0`).
///
/// With the effective inertia `I' = I + dt · S` (`S` = `torque_stiffness`) the
/// free wheel moves to `ω_free = ω + dt · T_net / I'`, and the brake then
/// removes up to `dt · T_brake / I'` from `|ω_free|` without reversing it (see
/// [`apply_brake_torque`]).  Using the same `I'` for both keeps the fixed
/// points of the step identical to the physical ones for every `dt` and `S`: a
/// wheel in torque equilibrium (`T_net = T_brake · sign(ω)`) keeps its spin, and
/// a stationary wheel stays at rest while `|T_net| <= T_brake` (static brake
/// friction).  For `S = 0` the step is explicit Euler; for a large `S` the spin
/// moves to the rolling condition without overshooting it.
pub(crate) fn advance_wheel_spin(
    spin: Real,
    net_torque: Real,
    torque_stiffness: Real,
    brake_torque: Real,
    inertia: Real,
    dt: Real,
) -> Real {
    let effective_inertia = inertia + torque_stiffness.max(0.0) * dt;
    let free_spin = spin + net_torque * dt / effective_inertia;
    apply_brake_torque(free_spin, brake_torque, effective_inertia, dt)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tire::{LinearTire, PacejkaTire};

    // --- apply_brake_torque ---

    #[test]
    fn test_apply_brake_torque_removes_spin_magnitude_in_both_directions() {
        // 600 N m on I = 1 kg m^2 for 0.01 s removes 6 rad/s from |spin|.
        assert!((apply_brake_torque(10.0, 600.0, 1.0, 0.01) - 4.0).abs() < 1e-12);
        assert!((apply_brake_torque(-10.0, 600.0, 1.0, 0.01) + 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_apply_brake_torque_never_reverses_a_wheel() {
        assert_eq!(apply_brake_torque(3.0, 600.0, 1.0, 0.01), 0.0);
        assert_eq!(apply_brake_torque(-3.0, 600.0, 1.0, 0.01), 0.0);
        assert!((apply_brake_torque(6.1, 600.0, 1.0, 0.01) - 0.1).abs() < 1e-12);
        assert!((apply_brake_torque(-6.1, 600.0, 1.0, 0.01) + 0.1).abs() < 1e-12);
        // A wheel at rest stays at rest, however large the brake torque.
        assert_eq!(apply_brake_torque(0.0, 1.0e6, 1.0, 0.01), 0.0);
        assert_eq!(apply_brake_torque(-0.0, 1.0e6, 1.0, 0.01), 0.0);
    }

    #[test]
    fn test_apply_brake_torque_without_brake_leaves_spin_untouched() {
        assert_eq!(apply_brake_torque(7.5, 0.0, 1.0, 0.01), 7.5);
        assert_eq!(apply_brake_torque(-7.5, -100.0, 1.0, 0.01), -7.5);
        assert_eq!(apply_brake_torque(7.5, f64::NAN, 1.0, 0.01), 7.5);
        assert_eq!(apply_brake_torque(0.0, 0.0, 1.0, 0.01), 0.0);
    }

    // --- longitudinal_slip_stiffness ---

    #[test]
    fn test_slip_stiffness_of_linear_tire_is_its_longitudinal_stiffness() {
        let tire = LinearTire::new(40_000.0, 60_000.0);
        let (normal, mu, slip) = (5000.0, 1.0, 0.01);
        let (_, fx) = tire.compute_forces(0.0, slip, normal, mu);
        let stiffness = longitudinal_slip_stiffness(&tire, 0.0, slip, normal, mu, fx);
        assert!(
            (stiffness - 60_000.0).abs() < 1e-6 * 60_000.0,
            "stiffness = {stiffness}"
        );
    }

    #[test]
    fn test_slip_stiffness_vanishes_in_the_saturated_range() {
        let tire = LinearTire::new(40_000.0, 60_000.0);
        let (normal, mu, slip) = (5000.0, 1.0, 0.5);
        let (_, fx) = tire.compute_forces(0.0, slip, normal, mu);
        assert!((fx - mu * normal).abs() < 1e-9, "tire must be saturated");
        assert_eq!(
            longitudinal_slip_stiffness(&tire, 0.0, slip, normal, mu, fx),
            0.0
        );
    }

    #[test]
    fn test_slip_stiffness_of_pacejka_tire_matches_the_initial_slope() {
        // Magic Formula: Fx ~ B * C * D * mu * Fz * kappa for small slip, and a
        // zero slip angle leaves the longitudinal force unweighted.
        let tire = PacejkaTire::new();
        let (normal, mu) = (4000.0, 1.2);
        let slip = 1.0e-6;
        let (_, fx) = tire.compute_forces(0.0, slip, normal, mu);
        let stiffness = longitudinal_slip_stiffness(&tire, 0.0, slip, normal, mu, fx);
        let expected =
            tire.longitudinal.b * tire.longitudinal.c * tire.longitudinal.d * mu * normal;
        assert!(
            (stiffness - expected).abs() < 1e-3 * expected,
            "stiffness = {stiffness}, expected {expected}"
        );
    }

    #[test]
    fn test_slip_stiffness_is_never_negative_past_the_force_peak() {
        let tire = PacejkaTire::new();
        let (normal, mu, slip) = (4000.0, 1.0, 0.9);
        let (_, fx) = tire.compute_forces(0.0, slip, normal, mu);
        let (_, fx_probe) = tire.compute_forces(0.0, slip + SLIP_STIFFNESS_PROBE, normal, mu);
        assert!(
            fx_probe < fx,
            "the Magic Formula falls beyond its peak: {fx_probe} vs {fx}"
        );
        assert_eq!(
            longitudinal_slip_stiffness(&tire, 0.0, slip, normal, mu, fx),
            0.0
        );
    }

    #[test]
    fn test_slip_stiffness_without_grip_is_zero() {
        let tire = LinearTire::default();
        assert_eq!(
            longitudinal_slip_stiffness(&tire, 0.0, 0.01, 0.0, 1.0, 0.0),
            0.0
        );
        assert_eq!(
            longitudinal_slip_stiffness(&tire, 0.0, 0.01, 1000.0, f64::NAN, 0.0),
            0.0
        );
    }

    // --- advance_wheel_spin ---

    #[test]
    fn test_advance_wheel_spin_without_stiffness_is_explicit_euler() {
        // 100 N m on I = 2 kg m^2 for 0.01 s adds 0.5 rad/s.
        let spin = advance_wheel_spin(10.0, 100.0, 0.0, 0.0, 2.0, 0.01);
        assert!((spin - 10.5).abs() < 1e-12);
        // The same step with a brake of 600 N m (removes 3 rad/s).
        let braked = advance_wheel_spin(10.0, 100.0, 0.0, 600.0, 2.0, 0.01);
        assert!((braked - 7.5).abs() < 1e-12);
    }

    #[test]
    fn test_advance_wheel_spin_keeps_a_wheel_in_torque_equilibrium() {
        for stiffness in [0.0, 1.0, 1.0e4, 1.0e9] {
            assert_eq!(
                advance_wheel_spin(12.5, 0.0, stiffness, 0.0, 0.7, 0.02),
                12.5
            );
        }
    }

    #[test]
    fn test_advance_wheel_spin_brake_equilibrium_does_not_depend_on_the_stiffness() {
        // The tire reaction balances the brake: T_net = T_brake * sign(spin).  The
        // wheel must keep its spin for every stiffness and step, in both
        // directions, because the brake uses the same effective inertia as the
        // implicit step.
        for stiffness in [0.0, 1.0e3, 1.0e5] {
            for dt in [1.0e-4, 1.0 / 240.0, 1.0 / 60.0] {
                let forward = advance_wheel_spin(20.0, 600.0, stiffness, 600.0, 0.7, dt);
                assert!((forward - 20.0).abs() < 1e-12, "forward {forward}");
                let backward = advance_wheel_spin(-20.0, -600.0, stiffness, 600.0, 0.7, dt);
                assert!((backward + 20.0).abs() < 1e-12, "backward {backward}");
            }
        }
    }

    #[test]
    fn test_advance_wheel_spin_holds_a_stationary_wheel_within_the_static_brake_torque() {
        for stiffness in [0.0, 300.0, 1.0e5] {
            for net in [-600.0, -250.0, 0.0, 250.0, 600.0] {
                assert_eq!(
                    advance_wheel_spin(0.0, net, stiffness, 600.0, 1.0, 0.01),
                    0.0
                );
            }
        }
        // Beyond the static brake torque the wheel breaks away in the direction of
        // the net torque; dt * (|T_net| - T_brake) / (I + dt * S) = 0.25 here.
        let forward = advance_wheel_spin(0.0, 700.0, 300.0, 600.0, 1.0, 0.01);
        assert!((forward - 0.25).abs() < 1e-12, "forward {forward}");
        let backward = advance_wheel_spin(0.0, -700.0, 300.0, 600.0, 1.0, 0.01);
        assert!((backward + 0.25).abs() < 1e-12, "backward {backward}");
    }

    #[test]
    fn test_advance_wheel_spin_brake_never_reverses_a_stiffly_coupled_wheel() {
        // Rolling forward under a huge brake torque, with the tire pushing the
        // wheel forward: the wheel stops at zero instead of turning backwards.
        let spin = advance_wheel_spin(2.0, 50.0, 2.0e4, 1.0e5, 0.7, 1.0 / 60.0);
        assert_eq!(spin, 0.0);
        let reverse = advance_wheel_spin(-2.0, -50.0, 2.0e4, 1.0e5, 0.7, 1.0 / 60.0);
        assert_eq!(reverse, 0.0);
    }

    #[test]
    fn test_advance_wheel_spin_treats_invalid_stiffness_as_zero() {
        let explicit = advance_wheel_spin(3.0, 40.0, 0.0, 0.0, 1.0, 0.01);
        assert_eq!(
            advance_wheel_spin(3.0, 40.0, -5.0e3, 0.0, 1.0, 0.01),
            explicit
        );
        assert_eq!(
            advance_wheel_spin(3.0, 40.0, f64::NAN, 0.0, 1.0, 0.01),
            explicit
        );
    }

    #[test]
    fn test_advance_wheel_spin_damps_the_step_with_stiffness() {
        let free = advance_wheel_spin(0.0, 100.0, 0.0, 0.0, 1.0, 0.01);
        let damped = advance_wheel_spin(0.0, 100.0, 300.0, 0.0, 1.0, 0.01);
        // dt * S / I = 3: the step shrinks by a factor of 1 + 3.
        assert!((damped - free / 4.0).abs() < 1e-12, "{damped} vs {free}");
    }

    /// Relaxation of a wheel towards the rolling condition `ω·r = v` under a
    /// linear tire, iterated with the given step; returns the spin history.
    fn relax_towards_rolling(stiff_torque: f64, inertia: f64, dt: f64, steps: usize) -> Vec<f64> {
        // r = 0.3 m, v = 10 m/s: rolling spin 100/3 rad/s.  The tire torque is
        // -S * (ω - ω_roll) with S = r² * k / den.
        let omega_roll = 10.0 / 0.3;
        let mut spin = 0.0;
        let mut history = Vec::with_capacity(steps);
        for _ in 0..steps {
            let net_torque = -stiff_torque * (spin - omega_roll);
            spin = advance_wheel_spin(spin, net_torque, stiff_torque, 0.0, inertia, dt);
            history.push(spin);
        }
        history
    }

    #[test]
    fn test_advance_wheel_spin_converges_monotonically_for_any_step() {
        // S / I = 16 000 / s; steps from 1e-5 s (rate * dt = 0.16) to 0.1 s
        // (rate * dt = 1600) must all approach the rolling spin from below.
        let inertia = 0.675;
        let stiff_torque = 16_000.0 * inertia;
        let omega_roll = 10.0 / 0.3;
        for dt in [1.0e-5, 1.0e-4, 1.0 / 240.0, 1.0 / 120.0, 1.0 / 60.0, 0.1] {
            let history = relax_towards_rolling(stiff_torque, inertia, dt, 400);
            let mut previous = 0.0;
            for spin in &history {
                assert!(
                    *spin >= previous && *spin <= omega_roll * (1.0 + 1e-12),
                    "dt = {dt}: spin {spin} left the monotone approach (previous {previous})"
                );
                previous = *spin;
            }
            let last = history[history.len() - 1];
            assert!(
                (last - omega_roll).abs() < 1e-6 * omega_roll,
                "dt = {dt}: final spin {last}"
            );
        }
    }

    #[test]
    fn test_explicit_euler_would_diverge_where_the_implicit_step_converges() {
        // Same relaxation with the stiffness switched off in the update: the
        // amplification factor 1 - rate * dt = -265.7 at dt = 1/60 blows up.
        let inertia = 0.675;
        let stiff_torque = 16_000.0 * inertia;
        let omega_roll = 10.0 / 0.3;
        let dt = 1.0 / 60.0;
        let mut spin = 0.0;
        for _ in 0..4 {
            let net_torque = -stiff_torque * (spin - omega_roll);
            spin = advance_wheel_spin(spin, net_torque, 0.0, 0.0, inertia, dt);
        }
        assert!(
            (spin - omega_roll).abs() > 1.0e6,
            "explicit Euler must be unstable here, got {spin}"
        );
        let implicit = relax_towards_rolling(stiff_torque, inertia, dt, 4);
        assert!((implicit[3] - omega_roll).abs() < 0.01 * omega_roll);
    }

    #[test]
    fn test_advance_wheel_spin_reaches_the_traction_equilibrium_under_drive_torque() {
        // Linear tire: T_drive = S * (ω - ω_roll) at equilibrium, i.e. the wheel
        // settles above the rolling spin by T_drive / S.
        let inertia = 0.675;
        let stiff_torque = 16_000.0 * inertia;
        let omega_roll = 10.0 / 0.3;
        let drive = 300.0;
        let dt = 1.0 / 60.0;
        let mut spin = omega_roll;
        for _ in 0..200 {
            let net_torque = drive - stiff_torque * (spin - omega_roll);
            spin = advance_wheel_spin(spin, net_torque, stiff_torque, 0.0, inertia, dt);
        }
        let expected = omega_roll + drive / stiff_torque;
        assert!((spin - expected).abs() < 1e-9, "{spin} vs {expected}");
    }
}
