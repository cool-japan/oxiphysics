//! Regression tests for the Smagorinsky relaxation-time closure shared by the
//! 2-D LES collision ([`smagorinsky_omega`]) and the wall-modelled channel
//! solver ([`collide_stream_wmles`]).
//!
//! The closure `tau_eff = tau_0 + C_s^2 |S| / c_s^2` with
//! `S = -Pi_neq / (2 rho c_s^2 tau_eff)` and `|S| = sqrt(2 S:S)` is solved here
//! by fixed-point iteration, independently of the closed-form root used by the
//! solvers, so that a mis-scaled coefficient shows up as a large mismatch.

use super::functions::{smagorinsky_omega, smagorinsky_tau};
use super::wmles::{WmlesChannelCoupling, collide_stream_wmles};
use crate::grid::LbmGrid2D;
use crate::lattice::{CS2, D2Q9_VELOCITIES, D2Q9_WEIGHTS, LatticeType};
use crate::turbulent_channel::ChannelFlow;

/// D2Q9 populations with density `rho`, velocity `u` and the non-equilibrium
/// momentum flux `pi = [Pxx, Pxy, Pyy]` and nothing else:
/// `f_i = feq_i + w_i / (2 c_s^4) Pi_ab (c_ia c_ib - c_s^2 delta_ab)`.
///
/// The stress term has no mass and no momentum, so the macroscopic state of the
/// populations stays `(rho, u)` and their second moment of `f - feq` is `pi`.
fn stressed_populations(rho: f64, u: [f64; 2], pi: [f64; 3]) -> [f64; 9] {
    let mut f = [0.0_f64; 9];
    for (i, f_i) in f.iter_mut().enumerate() {
        let w = D2Q9_WEIGHTS[i];
        let cx = f64::from(D2Q9_VELOCITIES[i][0]);
        let cy = f64::from(D2Q9_VELOCITIES[i][1]);
        let cu = cx * u[0] + cy * u[1];
        let u_sq = u[0] * u[0] + u[1] * u[1];
        let feq = w * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
        let contraction = pi[0] * (cx * cx - CS2) + 2.0 * pi[1] * cx * cy + pi[2] * (cy * cy - CS2);
        *f_i = feq + w / (2.0 * CS2 * CS2) * contraction;
    }
    f
}

/// Frobenius norm `sqrt(Pxx^2 + 2 Pxy^2 + Pyy^2)` of a symmetric 2x2 tensor.
fn frobenius(pi: [f64; 3]) -> f64 {
    (pi[0] * pi[0] + 2.0 * pi[1] * pi[1] + pi[2] * pi[2]).sqrt()
}

/// Solve `tau = tau_0 + 3 C_s^2 |S(tau)|`, `S = -Pi / (2 rho c_s^2 tau)`, by
/// fixed-point iteration (a contraction for the stresses used here).
fn closure_tau_by_iteration(tau0: f64, cs_smag: f64, rho: f64, pi: [f64; 3]) -> f64 {
    let pi_sq = frobenius(pi).powi(2);
    let mut tau = tau0;
    for _ in 0..200 {
        let s_sq = pi_sq / (2.0 * rho * CS2 * tau).powi(2);
        let strain_magnitude = (2.0 * s_sq).sqrt();
        tau = tau0 + cs_smag * cs_smag * strain_magnitude / CS2;
    }
    tau
}

#[test]
fn smagorinsky_tau_matches_the_closed_form_coefficient() {
    // tau_eff = (tau_0 + sqrt(tau_0^2 + 18 sqrt(2) C_s^2 ||Pi|| / rho)) / 2.
    let (tau0, cs_smag, pi_norm, rho) = (0.8, 0.15, 3.0e-3, 1.1);
    let expected = 0.5
        * (tau0
            + (tau0 * tau0 + 18.0 * std::f64::consts::SQRT_2 * cs_smag * cs_smag * pi_norm / rho)
                .sqrt());
    let tau_eff = smagorinsky_tau(tau0, cs_smag, pi_norm, rho);
    assert!(
        (tau_eff - expected).abs() < 1e-13,
        "tau_eff = {tau_eff}, closed form = {expected}"
    );
    // The previous implementation divided by c_s^4 a second time, which made the
    // added relaxation time about 6.4 times too large.
    let doubly_scaled = 0.5
        * (tau0 + (tau0 * tau0 + 18.0 * cs_smag * cs_smag * pi_norm / (rho * CS2 * CS2)).sqrt());
    let ratio = (doubly_scaled - tau0) / (tau_eff - tau0);
    assert!(
        ratio > 6.0 && ratio < 6.6,
        "the old coefficient overshoots by {ratio}"
    );
}

#[test]
fn smagorinsky_tau_limits() {
    assert!((smagorinsky_tau(0.7, 0.17, 0.0, 1.0) - 0.7).abs() < 1e-15);
    assert!((smagorinsky_tau(0.7, 0.0, 1.0e-3, 1.0) - 0.7).abs() < 1e-15);
    // A cell without mass carries no strain information.
    assert!((smagorinsky_tau(0.7, 0.17, 1.0e-3, 0.0) - 0.7).abs() < 1e-15);
    let weak = smagorinsky_tau(0.7, 0.17, 1.0e-4, 1.0);
    let strong = smagorinsky_tau(0.7, 0.17, 1.0e-2, 1.0);
    assert!(0.7 < weak && weak < strong);
    // For a weak stress the added relaxation time is quadratic in C_s.
    let small = smagorinsky_tau(1.0, 0.05, 1.0e-6, 1.0) - 1.0;
    let large = smagorinsky_tau(1.0, 0.10, 1.0e-6, 1.0) - 1.0;
    assert!(
        (large / small - 4.0).abs() < 1e-6,
        "ratio {}",
        large / small
    );
}

#[test]
fn smagorinsky_omega_solves_the_eddy_viscosity_closure() {
    let rho = 1.02;
    let u = [0.03, -0.015];
    let pi = [4.0e-3, -6.0e-3, -2.5e-3];
    let cs_smag = 0.17;
    let tau0 = 0.7;

    let mut grid = LbmGrid2D::new(5, 5, LatticeType::D2Q9);
    let k = grid.idx(2, 3);
    let populations = stressed_populations(rho, u, pi);
    for (plane, &value) in grid.f.iter_mut().zip(&populations) {
        plane[k] = value;
    }
    grid.compute_macroscopic();
    assert!((grid.rho[k] - rho).abs() < 1e-14);
    assert!((grid.ux[k] - u[0]).abs() < 1e-14);
    assert!((grid.uy[k] - u[1]).abs() < 1e-14);

    let omega = smagorinsky_omega(&grid, cs_smag, 1.0 / tau0, 2, 3);
    let tau_eff = 1.0 / omega;
    let tau_reference = closure_tau_by_iteration(tau0, cs_smag, rho, pi);
    assert!(
        (tau_eff - tau_reference).abs() < 1e-12,
        "tau_eff = {tau_eff}, closure = {tau_reference}"
    );

    // The eddy viscosity equals (C_s)^2 |S| with |S| recovered from the stress.
    let strain_magnitude = frobenius(pi) / (std::f64::consts::SQRT_2 * rho * CS2 * tau_eff);
    let eddy_viscosity = CS2 * (tau_eff - tau0);
    assert!(
        (eddy_viscosity - cs_smag * cs_smag * strain_magnitude).abs() < 1e-12,
        "nu_t = {eddy_viscosity:e}, Cs^2 |S| = {:e}",
        cs_smag * cs_smag * strain_magnitude
    );
    assert!(eddy_viscosity > 0.0);
}

#[test]
fn smagorinsky_omega_without_stress_keeps_the_base_frequency() {
    let mut grid = LbmGrid2D::new(4, 4, LatticeType::D2Q9);
    let k = grid.idx(1, 1);
    let populations = stressed_populations(1.0, [0.04, 0.02], [0.0; 3]);
    for (plane, &value) in grid.f.iter_mut().zip(&populations) {
        plane[k] = value;
    }
    grid.compute_macroscopic();
    let omega = smagorinsky_omega(&grid, 0.2, 1.4, 1, 1);
    assert!((omega - 1.4).abs() < 1e-13, "omega = {omega}");
}

#[test]
fn bgk_collide_smagorinsky_relaxes_the_stress_with_the_local_frequency() {
    let rho = 1.0;
    let u = [0.02, 0.01];
    let pi = [3.0e-3, 2.0e-3, -1.0e-3];
    let (cs_smag, tau0) = (0.18, 0.65);

    let mut grid = LbmGrid2D::new(4, 4, LatticeType::D2Q9);
    let k = grid.idx(1, 2);
    let populations = stressed_populations(rho, u, pi);
    for (plane, &value) in grid.f.iter_mut().zip(&populations) {
        plane[k] = value;
    }
    super::functions::bgk_collide_smagorinsky_2d(&mut grid, 1.0 / tau0, cs_smag);

    // The stress after the collision is (1 - omega_local) times the stress before.
    let equilibrium = stressed_populations(rho, u, [0.0; 3]);
    let (mut p_xx, mut p_xy, mut p_yy) = (0.0, 0.0, 0.0);
    for i in 0..9 {
        let cx = f64::from(D2Q9_VELOCITIES[i][0]);
        let cy = f64::from(D2Q9_VELOCITIES[i][1]);
        let f_neq = grid.f[i][k] - equilibrium[i];
        p_xx += f_neq * cx * cx;
        p_xy += f_neq * cx * cy;
        p_yy += f_neq * cy * cy;
    }
    let tau_reference = closure_tau_by_iteration(tau0, cs_smag, rho, pi);
    let retained = 1.0 - 1.0 / tau_reference;
    assert!((p_xx - retained * pi[0]).abs() < 1e-12, "Pxx = {p_xx}");
    assert!((p_xy - retained * pi[1]).abs() < 1e-12, "Pxy = {p_xy}");
    assert!((p_yy - retained * pi[2]).abs() < 1e-12, "Pyy = {p_yy}");
}

#[test]
fn wmles_interior_collision_uses_the_smagorinsky_closure() {
    // A channel that is uniform in both directions.  Row 4 of 9 only receives
    // populations from interior rows (3, 4, 5), so after one call it holds the
    // post-collision state of the Smagorinsky cells.
    let (nx, ny) = (6, 9);
    let nu = 0.2;
    let cs_smag = 0.17;
    let rho = 1.03;
    let u = [0.025, -0.01];
    let pi = [5.0e-3, -4.0e-3, 2.0e-3];

    let mut channel = ChannelFlow::new(nx, ny, 180.0);
    channel.forcing = 0.0;
    let coupling = WmlesChannelCoupling::new(nu, cs_smag);
    let before = stressed_populations(rho, u, pi);
    for node in 0..nx * ny {
        channel.f_dist[node * 9..node * 9 + 9].copy_from_slice(&before);
    }

    collide_stream_wmles(&mut channel, &coupling);

    let tau0 = 3.0 * nu + 0.5;
    let tau_eff = closure_tau_by_iteration(tau0, cs_smag, rho, pi);
    assert!(tau_eff > tau0);
    let omega = 1.0 / tau_eff;
    let equilibrium = stressed_populations(rho, u, [0.0; 3]);
    for i in [0, 2, nx - 1] {
        for (q, (&f_q, &feq_q)) in before.iter().zip(&equilibrium).enumerate() {
            let expected = f_q + omega * (feq_q - f_q);
            let actual = channel.f_dist[channel.fi(i, 4, q)];
            assert!(
                (actual - expected).abs() < 1e-13,
                "column {i}, direction {q}: {actual} vs {expected}"
            );
        }
    }
}
