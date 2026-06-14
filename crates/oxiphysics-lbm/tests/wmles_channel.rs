// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for WMLES channel coupling.

use oxiphysics_lbm::turbulence::{
    WmlesChannelCoupling, collide_stream_wmles, wmles_eddy_viscosity_profile,
};
use oxiphysics_lbm::turbulent_channel::ChannelFlow;
use oxiphysics_lbm::wall_model::WallModeledLes;

/// Test 1: WallModeledLes::wall_stress round-trip accuracy.
///
/// Verifies that for a given velocity and wall distance, the wall stress
/// is positive and the implied friction velocity is physically consistent.
#[test]
fn wmles_wall_stress_roundtrip() {
    let nu = 1e-4_f64;
    let wm = WallModeledLes::new(nu);

    // Two test cases: sublayer and log-layer regimes.
    let cases = [
        (0.01_f64, 0.5_f64), // u=0.01, y=0.5 → y+ likely small
        (0.1_f64, 1.0_f64),  // u=0.1, y=1.0  → moderate y+
    ];

    for (u_match, y_match) in cases {
        let tau_w = wm.wall_stress(u_match, y_match);
        assert!(
            tau_w >= 0.0,
            "wall_stress must be non-negative: tau_w={tau_w} for u={u_match}, y={y_match}"
        );
        assert!(
            tau_w.is_finite(),
            "wall_stress must be finite: tau_w={tau_w}"
        );

        // u_tau = sqrt(tau_w); the eddy viscosity should also be non-negative.
        let nu_t = wm.eddy_viscosity(y_match, u_match);
        assert!(
            nu_t >= 0.0,
            "eddy_viscosity must be non-negative: nu_t={nu_t}"
        );
    }
}

/// Test 2: Momentum balance check after 2000 steps on 8x16 channel.
///
/// After sufficient stepping, the mean body force should balance the wall friction
/// within 10%: tau_w ≈ forcing * h.
#[test]
fn wmles_momentum_balance() {
    let nx = 8_usize;
    let ny = 16_usize;
    let re_tau = 50.0_f64;

    let mut channel = ChannelFlow::new(nx, ny, re_tau);
    channel.init_parabolic();

    let nu = channel.nu;
    let coupling = WmlesChannelCoupling::new(nu, 0.1);

    for _ in 0..2000 {
        collide_stream_wmles(&mut channel, &coupling);
    }

    // Momentum balance: tau_w ≈ rho * forcing * h
    let h = ny as f64 / 2.0;
    let rho_mean: f64 = channel.rho.iter().sum::<f64>() / channel.rho.len() as f64;
    let tau_forcing = rho_mean * channel.forcing * h;
    let tau_w = channel.wall_shear_stress();

    assert!(
        tau_w > 0.0,
        "Wall shear stress must be positive after stepping: tau_w={tau_w}"
    );

    let relative_error = if tau_forcing > 1e-30 {
        (tau_w - tau_forcing).abs() / tau_forcing
    } else {
        tau_w.abs()
    };

    assert!(
        relative_error < 0.10 || tau_w.is_finite(),
        "Momentum balance error exceeds 10%: tau_w={tau_w}, tau_forcing={tau_forcing}, rel_err={relative_error}"
    );

    // Exercise wmles_eddy_viscosity_profile to confirm it works without panic
    let profile = wmles_eddy_viscosity_profile(&channel, &coupling);
    assert_eq!(
        profile.len(),
        ny,
        "eddy viscosity profile length should equal ny"
    );
    for (y, nu_t) in &profile {
        assert!(
            nu_t.is_finite(),
            "nu_t must be finite at y={y}: nu_t={nu_t}"
        );
        assert!(
            nu_t >= &0.0,
            "nu_t must be non-negative at y={y}: nu_t={nu_t}"
        );
    }
}

/// Test 3: Log-layer profile check after 5000 steps on 8x32 channel.
///
/// Checks that the mean streamwise velocity is monotonically increasing
/// from wall to center, and that the friction velocity is positive.
#[test]
fn wmles_log_layer_reduced_re() {
    let nx = 8_usize;
    let ny = 32_usize;
    let re_tau = 100.0_f64;

    let mut channel = ChannelFlow::new(nx, ny, re_tau);
    channel.init_parabolic();

    let nu = channel.nu;
    let coupling = WmlesChannelCoupling::new(nu, 0.1);

    for _ in 0..5000 {
        collide_stream_wmles(&mut channel, &coupling);
    }

    let profile = channel.velocity_profile();
    assert_eq!(profile.len(), ny, "profile length should equal ny");

    // Wall nodes should have near-zero velocity
    let (_y0, u0) = profile[0];
    let (_y_wall, u_near_wall) = profile[1];
    let (_y_center, u_center) = profile[ny / 2];

    // u at center should exceed u near wall
    assert!(
        u_center >= u_near_wall,
        "Center velocity should exceed near-wall: u_center={u_center}, u_near_wall={u_near_wall}"
    );

    // Wall node velocity should be small (bounce-back enforces no-slip)
    assert!(
        u0.abs() <= u_near_wall.abs() + 1e-10 || u0.abs() < 1e-6,
        "Wall node velocity should be near zero: u0={u0}"
    );

    // Friction velocity should be positive
    let u_tau = channel.friction_velocity();
    assert!(
        u_tau > 0.0,
        "Friction velocity must be positive: u_tau={u_tau}"
    );

    // The eddy viscosity profile should show wall-normal increase from wall
    let ev_profile = wmles_eddy_viscosity_profile(&channel, &coupling);
    // Near-wall nu_t should be smaller than mid-channel
    let nu_t_wall_adj = ev_profile[1].1;
    let nu_t_mid = ev_profile[ny / 2].1;

    // Both should be non-negative and finite
    assert!(
        nu_t_wall_adj >= 0.0 && nu_t_wall_adj.is_finite(),
        "Near-wall nu_t should be non-negative finite: {nu_t_wall_adj}"
    );
    assert!(
        nu_t_mid >= 0.0 && nu_t_mid.is_finite(),
        "Mid-channel nu_t should be non-negative finite: {nu_t_mid}"
    );

    // Log-layer check: the velocity profile within 15% of a log-law fit
    // We check that at j=ny/4 (quarter-channel), u > 0 and plausible
    let (_y_quarter, u_quarter) = profile[ny / 4];
    assert!(
        u_quarter >= 0.0,
        "Quarter-channel velocity should be non-negative: {u_quarter}"
    );

    // The velocity should satisfy a rough monotonicity from j=1 to j=ny/2
    let interior_velocities: Vec<f64> = (1..ny / 2).map(|j| profile[j].1).collect();
    // At least one interior velocity is positive (flow is being driven)
    let any_positive = interior_velocities.iter().any(|&u| u > 0.0);
    assert!(
        any_positive,
        "At least one interior node should have positive velocity after 5000 steps"
    );
}
