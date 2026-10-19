// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the Karlin–Bösch–Chikatamarla collision of
//! [`KbcD2Q9`]: natural-moment split, entropic stabiliser, viscosity and
//! stability.

use super::*;
use std::f64::consts::PI;

/// Natural moments `M[p][q] = Σ_i f_i c_ix^p c_iy^q` for `p, q ∈ {0, 1, 2}`.
fn natural_moments(f: &[f64; 9]) -> [[f64; 3]; 3] {
    let mut m = [[0.0; 3]; 3];
    for i in 0..9 {
        let cx = C[i].0 as f64;
        let cy = C[i].1 as f64;
        for (p, row) in m.iter_mut().enumerate() {
            for (q, value) in row.iter_mut().enumerate() {
                *value += f[i] * cx.powi(p as i32) * cy.powi(q as i32);
            }
        }
    }
    m
}

/// Deviatoric stress `N = M_20 − M_02` of a moment table.
fn n_moment(m: &[[f64; 3]; 3]) -> f64 {
    m[2][0] - m[0][2]
}

/// Trace `T = M_20 + M_02` of a moment table.
fn t_moment(m: &[[f64; 3]; 3]) -> f64 {
    m[2][0] + m[0][2]
}

/// Entropic equilibrium with the density and momentum of `f`.
fn local_equilibrium(f: &[f64; 9]) -> [f64; 9] {
    let m = natural_moments(f);
    let rho = m[0][0];
    EntropicD2Q9::equilibrium(rho, [m[1][0] / rho, m[0][1] / rho])
}

/// Equilibrium at `(rho, u)` plus a mass- and momentum-free perturbation of
/// relative size `scale` exciting the shear moments (`shear = true`) and the
/// higher-order moments `T`, `Q_xyy`, `Q_yxx`, `A` (`higher = true`).
fn excited_state(rho: f64, u: [f64; 2], scale: f64, shear: bool, higher: bool) -> [f64; 9] {
    let mut f = EntropicD2Q9::equilibrium(rho, u);
    for i in 0..9 {
        let cx = C[i].0 as f64;
        let cy = C[i].1 as f64;
        let mut p = 0.0;
        if shear {
            p += 0.9 * cx * cy + 0.6 * (cx * cx - cy * cy);
        }
        if higher {
            // Trace, third-order and fourth-order modes; each is mass- and
            // momentum-free and carries no shear stress (N = Pi_xy = 0).
            p += 0.3 * (1.5 * (cx * cx + cy * cy) - 1.0) + 0.5 * cx * (3.0 * cy * cy - 1.0)
                - 0.4 * cy * (3.0 * cx * cx - 1.0)
                + 0.2 * (3.0 * cx * cx - 1.0) * (3.0 * cy * cy - 1.0);
        }
        f[i] += W[i] * rho * scale * p;
    }
    f
}

/// Quadratic entropy estimate `⟨Δ'|Δ'⟩` of the post-collision deviation for
/// a given stabiliser `gamma`.
fn quadratic_entropy(f: &[f64; 9], f_eq: &[f64; 9], beta: f64, gamma: f64) -> f64 {
    let mut delta = [0.0; 9];
    for i in 0..9 {
        delta[i] = f[i] - f_eq[i];
    }
    let (ds, dh) = kbc_split(&delta);
    let mut post = [0.0; 9];
    for i in 0..9 {
        post[i] = delta[i] - beta * (2.0 * ds[i] + gamma * dh[i]);
    }
    entropic_scalar_product(&post, &post, f_eq).unwrap_or(f64::INFINITY)
}

// ---------------------------------------------------------------------------
// Natural-moment split
// ---------------------------------------------------------------------------

#[test]
fn shear_part_carries_exactly_the_deviatoric_stress() {
    let delta = [
        0.013, -0.027, 0.041, 0.006, -0.019, 0.022, -0.008, 0.017, -0.031,
    ];
    let (ds, dh) = kbc_split(&delta);
    let md = natural_moments(&delta);
    let ms = natural_moments(&ds);
    let mh = natural_moments(&dh);
    assert!((n_moment(&ms) - n_moment(&md)).abs() < 1e-15);
    assert!((ms[1][1] - md[1][1]).abs() < 1e-15);
    // Every other natural moment of the shear part vanishes.
    for (p, q) in [(0, 0), (1, 0), (0, 1), (2, 1), (1, 2), (2, 2)] {
        assert!(ms[p][q].abs() < 1e-15, "M_{p}{q} of s = {}", ms[p][q]);
    }
    assert!(
        t_moment(&ms).abs() < 1e-15,
        "trace of s = {}",
        t_moment(&ms)
    );
    // The higher-order part carries no shear stress ...
    assert!(n_moment(&mh).abs() < 1e-15);
    assert!(mh[1][1].abs() < 1e-15);
    // ... and everything else of delta.
    for (p, q) in [(0, 0), (1, 0), (0, 1), (2, 1), (1, 2), (2, 2)] {
        assert!((mh[p][q] - md[p][q]).abs() < 1e-15, "M_{p}{q}");
    }
    assert!((t_moment(&mh) - t_moment(&md)).abs() < 1e-15);
    for i in 0..9 {
        assert!((ds[i] + dh[i] - delta[i]).abs() < 1e-17);
    }
}

#[test]
fn entropic_scalar_product_is_weighted_by_the_equilibrium() {
    let f_eq = EntropicD2Q9::equilibrium(1.1, [0.04, -0.03]);
    let x = [0.1, -0.2, 0.3, 0.05, -0.04, 0.02, 0.01, -0.03, 0.07];
    let y = [0.02, 0.01, -0.05, 0.03, 0.06, -0.01, 0.02, 0.04, -0.02];
    let expected: f64 = (0..9).map(|i| x[i] * y[i] / f_eq[i]).sum();
    let got = entropic_scalar_product(&x, &y, &f_eq);
    assert!(matches!(got, Some(v) if (v - expected).abs() < 1e-14));
    let mut bad = f_eq;
    bad[4] = 0.0;
    assert_eq!(entropic_scalar_product(&x, &y, &bad), None);
    // Undefined scalar product: the stabiliser falls back to 1/beta.
    assert_eq!(kbc_stabiliser(0.5, &x, &y, &bad), 2.0);
}

// ---------------------------------------------------------------------------
// Stabiliser
// ---------------------------------------------------------------------------

/// At rest the shear and fourth-order modes are orthogonal in the entropic
/// scalar product, so gamma = 1/beta and the higher-order moments land
/// exactly on their equilibrium values.
#[test]
fn orthogonal_parts_give_gamma_one_over_beta_and_equilibrate_the_ghosts() {
    let mut f = EntropicD2Q9::equilibrium(1.0, [0.0, 0.0]);
    for i in 0..9 {
        let cx = C[i].0 as f64;
        let cy = C[i].1 as f64;
        f[i] += W[i] * (0.05 * cx * cy + 0.02 * (3.0 * cx * cx - 1.0) * (3.0 * cy * cy - 1.0));
    }
    let f_eq = local_equilibrium(&f);
    for &nu in &[1e-4, 0.02, 0.1, 0.3] {
        let beta = relaxation_beta(nu);
        let tau = 3.0 * nu + 0.5;
        let (post, gamma) = kbc_collide_populations(&f, &f_eq, beta);
        assert!(
            (gamma - 1.0 / beta).abs() < 1e-12,
            "nu = {nu}: gamma {gamma}"
        );
        let m_pre = natural_moments(&f);
        let m_eq = natural_moments(&f_eq);
        let m_post = natural_moments(&post);
        for (p, q) in [(2, 1), (1, 2), (2, 2)] {
            assert!(
                (m_post[p][q] - m_eq[p][q]).abs() < 1e-14,
                "nu = {nu}: M_{p}{q} not equilibrated"
            );
        }
        assert!((t_moment(&m_post) - t_moment(&m_eq)).abs() < 1e-14);
        let expected_pi = m_eq[1][1] + (1.0 - 1.0 / tau) * (m_pre[1][1] - m_eq[1][1]);
        assert!((m_post[1][1] - expected_pi).abs() < 1e-14, "nu = {nu}");
    }
}

/// A pure shear perturbation leaves only rounding noise in the higher-order
/// part; the noise floor must keep it from steering gamma.
#[test]
fn rounding_level_higher_order_part_does_not_steer_gamma() {
    for &scale in &[1e-7, 1e-3, 0.1] {
        let f = excited_state(1.0, [0.05, 0.02], scale, true, false);
        let f_eq = local_equilibrium(&f);
        let beta = relaxation_beta(0.05);
        let (_, gamma) = kbc_collide_populations(&f, &f_eq, beta);
        assert_eq!(gamma, 1.0 / beta, "scale {scale}");
    }
}

/// Away from rest the shear and higher-order parts are correlated, gamma
/// departs from 1/beta, yet the shear stress still relaxes with
/// tau = 3 nu + 1/2 and density and momentum are conserved.
#[test]
fn shear_relaxes_with_the_requested_viscosity_and_moments_are_conserved() {
    let f = excited_state(1.05, [0.08, -0.05], 0.1, true, true);
    let f_eq = local_equilibrium(&f);
    let m_pre = natural_moments(&f);
    let m_eq = natural_moments(&f_eq);
    // The perturbation itself is mass- and momentum-free.
    assert!((m_pre[0][0] - 1.05).abs() < 1e-14);
    assert!((m_pre[1][0] - 1.05 * 0.08).abs() < 1e-14);
    assert!((m_pre[0][1] + 1.05 * 0.05).abs() < 1e-14);
    for &nu in &[1e-5, 0.01, 0.05, 1.0 / 6.0, 0.4] {
        let beta = relaxation_beta(nu);
        let tau = 3.0 * nu + 0.5;
        let (post, gamma) = kbc_collide_populations(&f, &f_eq, beta);
        // At tau = 1 (beta = 1/2) the correlation term carries the factor
        // 2 - 1/beta = 0 and KBC coincides with BGK; elsewhere it is active.
        if (2.0 - 1.0 / beta).abs() > 0.1 {
            assert!(
                (gamma - 1.0 / beta).abs() > 1e-3,
                "nu = {nu}: gamma {gamma} is not stabilising"
            );
        } else {
            assert!(
                (gamma - 1.0 / beta).abs() < 1e-12,
                "nu = {nu}: gamma {gamma}"
            );
        }
        let m_post = natural_moments(&post);
        for (p, q) in [(0, 0), (1, 0), (0, 1)] {
            assert!(
                (m_post[p][q] - m_pre[p][q]).abs() < 1e-14,
                "nu = {nu}: conserved moment M_{p}{q} changed"
            );
        }
        let factor = 1.0 - 1.0 / tau;
        let expected_n = n_moment(&m_eq) + factor * (n_moment(&m_pre) - n_moment(&m_eq));
        let expected_pi = m_eq[1][1] + factor * (m_pre[1][1] - m_eq[1][1]);
        assert!(
            (n_moment(&m_post) - expected_n).abs() < 1e-14,
            "nu = {nu}: N"
        );
        assert!(
            (m_post[1][1] - expected_pi).abs() < 1e-14,
            "nu = {nu}: Pi_xy"
        );
    }
}

#[test]
fn gamma_minimises_the_quadratic_entropy_estimate() {
    let f = excited_state(1.0, [0.08, -0.05], 0.1, true, true);
    let f_eq = local_equilibrium(&f);
    for &nu in &[1e-4, 0.01, 0.05, 0.2] {
        let beta = relaxation_beta(nu);
        let (_, gamma) = kbc_collide_populations(&f, &f_eq, beta);
        let q_min = quadratic_entropy(&f, &f_eq, beta, gamma);
        for other in [gamma - 1e-3, gamma + 1e-3, 1.0 / beta, 2.0] {
            assert!(
                q_min < quadratic_entropy(&f, &f_eq, beta, other),
                "nu = {nu}: gamma {gamma} is not the minimiser (vs {other})"
            );
        }
    }
}

/// With excited higher-order moments the KBC collision dissipates at least
/// as much entropy as BGK (gamma = 2) and as relaxing the ghosts to
/// equilibrium (gamma = 1/beta), and it never increases H.
#[test]
fn kbc_dissipates_ghost_entropy_at_least_as_well_as_bgk() {
    for &scale in &[0.01, 0.05, 0.1, 0.3] {
        let f = excited_state(1.0, [0.08, -0.05], scale, true, true);
        let f_eq = local_equilibrium(&f);
        let h_pre = h_function(&f);
        for &nu in &[1e-4, 0.01, 0.05, 0.2] {
            let beta = relaxation_beta(nu);
            let (post, _) = kbc_collide_populations(&f, &f_eq, beta);
            let mut bgk = f;
            let mut ghosts_to_eq = f;
            let mut delta = [0.0; 9];
            for i in 0..9 {
                delta[i] = f[i] - f_eq[i];
                bgk[i] += 2.0 * beta * (f_eq[i] - f[i]);
            }
            let (ds, dh) = kbc_split(&delta);
            for i in 0..9 {
                ghosts_to_eq[i] = f[i] - 2.0 * beta * ds[i] - dh[i];
            }
            let h_kbc = h_function(&post);
            assert!(h_kbc < h_pre, "scale {scale}, nu {nu}: H grew");
            assert!(h_kbc <= h_function(&bgk), "scale {scale}, nu {nu}: vs BGK");
            assert!(
                h_kbc <= h_function(&ghosts_to_eq),
                "scale {scale}, nu {nu}: vs gamma = 1/beta"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Solver-level behaviour
// ---------------------------------------------------------------------------

const WAVE_NX: usize = 2;
const WAVE_NY: usize = 64;
/// Steps before the measurement starts: the wave is initialised at local
/// equilibrium (no viscous stress yet), and the start-up transient would
/// otherwise bias the decay rate by a few percent at small viscosities.
const WAVE_WARMUP: usize = 200;
const WAVE_STEPS: usize = 1000;

/// First Fourier amplitude of u_x(y) = A sin(2 pi y / ny).
fn shear_wave_amplitude(cells: &[[f64; 9]]) -> f64 {
    let mut acc = 0.0;
    for y in 0..WAVE_NY {
        let m = natural_moments(&cells[y * WAVE_NX]);
        acc += m[1][0] / m[0][0] * (2.0 * PI * y as f64 / WAVE_NY as f64).sin();
    }
    2.0 * acc / WAVE_NY as f64
}

/// Shear viscosity measured from the decay of a sinusoidal shear wave.
fn kbc_wave_viscosity(nu: f64) -> f64 {
    let mut solver = KbcD2Q9::new(WAVE_NX, WAVE_NY, nu);
    for y in 0..WAVE_NY {
        let ux = 0.01 * (2.0 * PI * y as f64 / WAVE_NY as f64).sin();
        let f = EntropicD2Q9::equilibrium(1.0, [ux, 0.0]);
        for x in 0..WAVE_NX {
            solver.f[y * WAVE_NX + x] = f;
        }
    }
    for _ in 0..WAVE_WARMUP {
        solver.step();
    }
    let a0 = shear_wave_amplitude(&solver.f);
    for _ in 0..WAVE_STEPS {
        solver.step();
    }
    let k = 2.0 * PI / WAVE_NY as f64;
    -(shear_wave_amplitude(&solver.f) / a0).ln() / (k * k * WAVE_STEPS as f64)
}

#[test]
fn kbc_shear_wave_viscosity_follows_the_requested_value() {
    for &nu in &[0.01, 0.05, 0.1] {
        let measured = kbc_wave_viscosity(nu);
        assert!(
            (measured - nu).abs() < 0.02 * nu,
            "nu = {nu}: measured {measured}"
        );
    }
}

/// Doubly periodic double shear layer (thin layers, Re = U L / nu large).
fn double_shear_layer(n: usize, u0: f64) -> Vec<[f64; 9]> {
    let l = n as f64;
    let kappa = 80.0;
    let mut cells = Vec::with_capacity(n * n);
    for y in 0..n {
        for x in 0..n {
            let yf = y as f64 / l;
            let xf = x as f64 / l;
            let ux = if yf <= 0.5 {
                u0 * (kappa * (yf - 0.25)).tanh()
            } else {
                u0 * (kappa * (0.75 - yf)).tanh()
            };
            let uy = 0.05 * u0 * (2.0 * PI * (xf + 0.25)).sin();
            cells.push(EntropicD2Q9::equilibrium(1.0, [ux, uy]));
        }
    }
    cells
}

/// At nu = 1e-4 (tau - 1/2 = 3e-4) on a coarse grid BGK diverges within a
/// few hundred steps while the KBC stabiliser keeps the flow bounded and
/// positive.
#[test]
fn kbc_stays_stable_where_bgk_diverges() {
    let n = 32;
    let nu = 1e-4;
    let steps = 500;
    let omega = 1.0 / (3.0 * nu + 0.5);

    let mut bgk = EntropicD2Q9::new(n, n, nu);
    bgk.f = double_shear_layer(n, 0.1);
    let mut bgk_diverged = false;
    for _ in 0..steps {
        for cell in bgk.f.iter_mut() {
            let f_eq = local_equilibrium(cell);
            for i in 0..9 {
                cell[i] += omega * (f_eq[i] - cell[i]);
            }
        }
        bgk.stream();
        let u_max = bgk.max_velocity();
        if !u_max.is_finite() || u_max > 1.0 {
            bgk_diverged = true;
            break;
        }
    }
    assert!(bgk_diverged, "reference BGK run was expected to diverge");

    let mut kbc = KbcD2Q9::new(n, n, nu);
    kbc.f = double_shear_layer(n, 0.1);
    let mass0 = kbc.total_density();
    for _ in 0..steps {
        kbc.step();
    }
    let u_max = kbc.max_velocity();
    assert!(u_max.is_finite() && u_max < 0.2, "KBC u_max = {u_max}");
    assert!(
        kbc.f.iter().flat_map(|c| c.iter()).all(|&v| v > 0.0),
        "KBC produced a non-positive population"
    );
    assert!((kbc.total_density() - mass0).abs() < 1e-9 * mass0);
}
