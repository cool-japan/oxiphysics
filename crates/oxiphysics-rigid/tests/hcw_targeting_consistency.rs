// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Cross-checks of the two public Hill-Clohessy-Wiltshire (HCW) two-impulse
//! entry points, `HillClohessyWiltshire::two_impulse_rendezvous` and
//! `SpacecraftFormation::two_impulse_transfer`.
//!
//! Both work in the frame `x` radial, `y` along-track, `z` cross-track and are
//! built on `hcw_two_impulse`.  Each burn pair is verified against a
//! propagator the solver does not use: the other module's closed-form
//! solution, and a direct numerical integration of the HCW differential
//! equations.

use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::planetary_mechanics::{HillClohessyWiltshire, hcw_two_impulse};
use oxiphysics_rigid::satellite_dynamics::{HillState, SpacecraftFormation};
use std::f64::consts::PI;

/// Radius of the chief's circular orbit \[m\].
const R_REF: f64 = 6.9e6;

/// One transfer problem: relative states of the deputy and of the target, and
/// the transfer time as a fraction of the orbital period.
struct Case {
    deputy: ([f64; 3], [f64; 3]),
    target: ([f64; 3], [f64; 3]),
    fraction_of_period: f64,
}

fn cases() -> [Case; 4] {
    [
        Case {
            deputy: ([500.0, -1_500.0, 200.0], [0.3, -0.5, 0.1]),
            target: ([0.0; 3], [0.0; 3]),
            fraction_of_period: 0.3,
        },
        Case {
            deputy: ([-250.0, 900.0, -40.0], [-0.2, 0.4, 0.05]),
            target: ([100.0, 50.0, 10.0], [0.01, 0.02, -0.01]),
            fraction_of_period: 0.45,
        },
        Case {
            deputy: ([0.0, -3_000.0, 0.0], [0.0; 3]),
            target: ([0.0; 3], [0.0; 3]),
            fraction_of_period: 0.25,
        },
        // More than one revolution, away from the singular transfer times.
        Case {
            deputy: ([800.0, 2_000.0, -300.0], [0.1, 0.0, 0.2]),
            target: ([-100.0, -400.0, 60.0], [0.0, 0.05, 0.0]),
            fraction_of_period: 1.2,
        },
    ]
}

fn state6(position: [f64; 3], velocity: [f64; 3]) -> [f64; 6] {
    [
        position[0],
        position[1],
        position[2],
        velocity[0],
        velocity[1],
        velocity[2],
    ]
}

fn vec3(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}

fn formation_for(case: &Case) -> SpacecraftFormation {
    SpacecraftFormation::new(
        R_REF,
        HillState::new(vec3(case.deputy.0), vec3(case.deputy.1)),
        HillState::new(vec3(case.target.0), vec3(case.target.1)),
    )
}

#[test]
fn formation_and_cw_model_return_the_same_burns() {
    let cw = HillClohessyWiltshire::new(R_REF);
    for (index, case) in cases().iter().enumerate() {
        let formation = formation_for(case);
        assert_eq!(formation.mean_motion(), cw.n, "case {index}: mean motion");
        let tf = case.fraction_of_period * 2.0 * PI / cw.n;
        let (dv1, dv2) = formation
            .two_impulse_transfer(tf)
            .unwrap_or_else(|| panic!("case {index}: formation transfer is singular"));
        let (cw_dv1, cw_dv2) = cw
            .two_impulse_rendezvous(
                state6(case.deputy.0, case.deputy.1),
                state6(case.target.0, case.target.1),
                tf,
            )
            .unwrap_or_else(|| panic!("case {index}: CW transfer is singular"));
        assert!(
            (dv1 - vec3(cw_dv1)).norm() <= 1e-12 * dv1.norm().max(1.0),
            "case {index}: first burn {dv1:?} vs {cw_dv1:?}"
        );
        assert!(
            (dv2 - vec3(cw_dv2)).norm() <= 1e-12 * dv2.norm().max(1.0),
            "case {index}: second burn {dv2:?} vs {cw_dv2:?}"
        );
    }
}

#[test]
fn formation_burns_reach_the_target_under_the_cw_propagator() {
    let cw = HillClohessyWiltshire::new(R_REF);
    for (index, case) in cases().iter().enumerate() {
        let tf = case.fraction_of_period * 2.0 * PI / cw.n;
        let (dv1, dv2) = formation_for(case)
            .two_impulse_transfer(tf)
            .unwrap_or_else(|| panic!("case {index}: formation transfer is singular"));
        let departure = state6(
            case.deputy.0,
            [
                case.deputy.1[0] + dv1.x,
                case.deputy.1[1] + dv1.y,
                case.deputy.1[2] + dv1.z,
            ],
        );
        let arrival = cw.propagate(departure, tf);
        let arrival_velocity = [arrival[3] + dv2.x, arrival[4] + dv2.y, arrival[5] + dv2.z];
        for axis in 0..3 {
            assert!(
                (arrival[axis] - case.target.0[axis]).abs() < 1e-6,
                "case {index}: position axis {axis} misses by {}",
                arrival[axis] - case.target.0[axis]
            );
            assert!(
                (arrival_velocity[axis] - case.target.1[axis]).abs() < 1e-8,
                "case {index}: velocity axis {axis} is off by {}",
                arrival_velocity[axis] - case.target.1[axis]
            );
        }
    }
}

#[test]
fn cw_burns_reach_the_target_under_the_formation_propagator() {
    let cw = HillClohessyWiltshire::new(R_REF);
    for (index, case) in cases().iter().enumerate() {
        let tf = case.fraction_of_period * 2.0 * PI / cw.n;
        let (dv1, dv2) = cw
            .two_impulse_rendezvous(
                state6(case.deputy.0, case.deputy.1),
                state6(case.target.0, case.target.1),
                tf,
            )
            .unwrap_or_else(|| panic!("case {index}: CW transfer is singular"));
        let mut formation = formation_for(case);
        formation.apply_impulse(vec3(dv1));
        let arrival = formation.propagate_exact(tf);
        let miss = (arrival.position - vec3(case.target.0)).norm();
        assert!(miss < 1e-6, "case {index}: arrival misses by {miss} m");
        let velocity_error = (arrival.velocity + vec3(dv2) - vec3(case.target.1)).norm();
        assert!(
            velocity_error < 1e-8,
            "case {index}: velocity after the second burn is off by {velocity_error} m/s"
        );
    }
}

#[test]
fn burns_agree_with_a_direct_integration_of_the_hcw_equations() {
    // The numerical stepper integrates x'' = 2 n y' + 3 n² x, y'' = -2 n x',
    // z'' = -n² z with a first-order scheme, so with a 5 ms step it tracks the
    // closed-form transfer to well under a metre over a fraction of an orbit.
    let dt = 0.005;
    for (index, case) in cases().iter().take(3).enumerate() {
        let mut formation = formation_for(case);
        let tf = case.fraction_of_period * 2.0 * PI / formation.mean_motion();
        let (dv1, dv2) = formation
            .two_impulse_transfer(tf)
            .unwrap_or_else(|| panic!("case {index}: formation transfer is singular"));
        formation.apply_impulse(dv1);
        let steps = (tf / dt).round() as usize;
        let step = tf / steps as f64;
        for _ in 0..steps {
            formation.step(step, Vec3::zeros());
        }
        let miss = (formation.deputy_state.position - vec3(case.target.0)).norm();
        assert!(
            miss < 1.0,
            "case {index}: numerically integrated arrival misses by {miss} m"
        );
        let velocity_error = (formation.deputy_state.velocity + dv2 - vec3(case.target.1)).norm();
        assert!(
            velocity_error < 1e-2,
            "case {index}: integrated arrival velocity is off by {velocity_error} m/s"
        );
    }
}

#[test]
fn solver_rejects_a_non_physical_mean_motion() {
    let (r0, v0) = ([100.0, 50.0, 0.0], [0.0; 3]);
    let (rf, vf) = ([0.0; 3], [0.0; 3]);
    assert!(hcw_two_impulse(1.1e-3, r0, v0, rf, vf, 1_000.0).is_some());
    assert!(hcw_two_impulse(0.0, r0, v0, rf, vf, 1_000.0).is_none());
    assert!(hcw_two_impulse(-1.1e-3, r0, v0, rf, vf, 1_000.0).is_none());
    assert!(hcw_two_impulse(f64::NAN, r0, v0, rf, vf, 1_000.0).is_none());
    assert!(hcw_two_impulse(f64::INFINITY, r0, v0, rf, vf, 1_000.0).is_none());
}
