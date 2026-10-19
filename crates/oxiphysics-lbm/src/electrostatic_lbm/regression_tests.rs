// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Regression tests for the physics of the electrostatic LBM solvers
//! (Poisson scaling and Dirichlet cells, signed Einstein-consistent carrier
//! transport, the coupled p-n junction, Nernst-Planck initialisation/drift and
//! the two-fluid plasma coupling).

use super::*;
use crate::error::Error;

/// Centroid (in cells) of a 1-D density profile along x of an `nx × ny` grid.
fn centroid_x(n: &[f64], nx: usize, ny: usize) -> f64 {
    let mut mass = 0.0;
    let mut moment = 0.0;
    for y in 0..ny {
        for x in 0..nx {
            let v = n[y * nx + x];
            mass += v;
            moment += v * x as f64;
        }
    }
    moment / mass
}

/// Centroid (in cells) along y.
fn centroid_y(n: &[f64], nx: usize, ny: usize) -> f64 {
    let mut mass = 0.0;
    let mut moment = 0.0;
    for y in 0..ny {
        for x in 0..nx {
            let v = n[y * nx + x];
            mass += v;
            moment += v * y as f64;
        }
    }
    moment / mass
}

/// Gaussian bump `amp·exp(−(x−x0)²/(2σ²))` plus a small floor.
fn bump(x: f64, x0: f64, sigma: f64, amp: f64) -> f64 {
    let d = x - x0;
    amp * (-(d * d) / (2.0 * sigma * sigma)).exp() + 1e-12 * amp
}

// ── PoissonLbm (F2522) ──────────────────────────────────────────────────────

/// Uniform charge between grounded electrodes: φ = ρ·x·(L − x)/(2ε).
fn uniform_charge_error(nx: usize, ny: usize, omega: f64) -> f64 {
    let dx = 1.0e-3;
    let eps_r = 2.0;
    let eps = eps_r * constants::EPSILON_0;
    let len = (nx - 1) as f64 * dx;
    // choose ρ so that the peak potential is 1 V
    let rho = 8.0 * eps / (len * len);
    let mut solver = PoissonLbm::new(nx, ny, eps_r, omega).with_spacing(dx);
    for y in 0..ny {
        for x in 0..nx {
            solver.set_charge(x, y, rho);
        }
        solver.set_phi(0, y, 0.0);
        solver.set_phi(nx - 1, y, 0.0);
    }
    let (_iters, res) = solver
        .solve(400_000, 1e-13)
        .expect("grounded problem is solvable");
    assert!(res <= 1e-13, "not converged: residual {res:e}");
    let mut worst: f64 = 0.0;
    for y in 0..ny {
        for x in 0..nx {
            let xp = x as f64 * dx;
            let exact = rho * xp * (len - xp) / (2.0 * eps);
            worst = worst.max((solver.get_phi(x, y) - exact).abs());
        }
    }
    worst // relative to the 1 V peak
}

#[test]
fn poisson_uniform_charge_matches_parabola() {
    let err = uniform_charge_error(41, 1, 1.0);
    assert!(err < 0.01, "max error {err:e} V of a 1 V peak");
}

#[test]
fn poisson_uniform_charge_matches_parabola_in_2d_strip() {
    let err = uniform_charge_error(21, 4, 1.0);
    assert!(err < 0.01, "max error {err:e} V of a 1 V peak");
}

#[test]
fn poisson_uniform_charge_other_relaxation_rates() {
    for omega in [0.8, 1.4] {
        let err = uniform_charge_error(41, 1, omega);
        assert!(
            err < 0.002,
            "ω = {omega}: max error {err:e} V of a 1 V peak"
        );
        // second-order convergence of the Dirichlet treatment
        let err_fine = uniform_charge_error(81, 1, omega);
        assert!(err_fine < 0.35 * err, "ω = {omega}: {err:e} → {err_fine:e}");
    }
}

#[test]
fn poisson_potential_scales_with_spacing_squared() {
    let solve_peak = |dx: f64| {
        let nx = 21;
        let mut s = PoissonLbm::new(nx, 1, 1.0, 1.0).with_spacing(dx);
        for x in 0..nx {
            s.set_charge(x, 0, 1.0e-6);
        }
        s.set_phi(0, 0, 0.0);
        s.set_phi(nx - 1, 0, 0.0);
        s.solve(100_000, 1e-16).expect("solvable");
        s.get_phi(nx / 2, 0)
    };
    let p1 = solve_peak(1.0e-3);
    let p2 = solve_peak(2.0e-3);
    assert!(((p2 / p1) - 4.0).abs() < 1e-6, "ratio {}", p2 / p1);
}

#[test]
fn poisson_dirichlet_cells_are_held_every_step() {
    let nx = 12;
    let mut s = PoissonLbm::new(nx, 1, 1.0, 1.0).with_spacing(1e-3);
    for x in 0..nx {
        s.set_charge(x, 0, 1e-8);
    }
    s.set_phi(0, 0, 0.25);
    s.set_phi(nx - 1, 0, -0.5);
    assert!(s.has_dirichlet());
    for _ in 0..50 {
        s.step();
        assert_eq!(s.get_phi(0, 0), 0.25);
        assert_eq!(s.get_phi(nx - 1, 0), -0.5);
    }
    // releasing the cells removes the constraint
    s.release_phi(0, 0);
    s.release_phi(nx - 1, 0);
    assert!(!s.has_dirichlet());
}

#[test]
fn poisson_periodic_net_charge_is_rejected() {
    let mut s = PoissonLbm::new(8, 8, 1.0, 1.0);
    s.set_charge(3, 3, 1.0e-9);
    match s.solve(100, 1e-9) {
        Err(Error::ConvergenceFailure {
            residual,
            iterations,
            ..
        }) => {
            assert_eq!(iterations, 0);
            assert!(residual > 0.0);
        }
        other => panic!("expected ConvergenceFailure, got {other:?}"),
    }
}

#[test]
fn poisson_periodic_neutral_dipole_solves() {
    let n = 16;
    let mut s = PoissonLbm::new(n, n, 1.0, 1.0).with_spacing(1e-6);
    s.set_charge(4, 8, 1.0e3);
    s.set_charge(12, 8, -1.0e3);
    let (_it, res) = s
        .solve(100_000, 1e-13)
        .expect("neutral periodic problem is solvable");
    assert!(res <= 1e-13);
    let mean: f64 = s.phi.iter().sum::<f64>() / (n * n) as f64;
    assert!(mean.abs() < 1e-9, "mean φ {mean:e}");
    assert!(s.get_phi(4, 8) > 0.0 && s.get_phi(12, 8) < 0.0);
    assert!((s.get_phi(4, 8) + s.get_phi(12, 8)).abs() < 1e-9 * s.get_phi(4, 8));
}

#[test]
fn poisson_rejects_invalid_relaxation_rate() {
    let mut s = PoissonLbm::new(4, 1, 1.0, 2.5);
    s.set_phi(0, 0, 0.0);
    assert!(matches!(
        s.solve(10, 1e-6),
        Err(Error::ParameterOutOfBounds { .. })
    ));
}

// ── ChargeTransportLbm (F2528) ──────────────────────────────────────────────

#[test]
fn electrons_and_holes_drift_in_opposite_directions() {
    let nx = 128;
    let dx = 1e-7;
    let e_field = 2.0e4;
    let make = |sign: f64| {
        let mut ct = ChargeTransportLbm::new(nx, 1, 0.1, 300.0, 1.0, dx).with_charge_sign(sign);
        for x in 0..nx {
            ct.n[x] = bump(x as f64, 64.0, 4.0, 1e20);
            ct.ex[x] = e_field;
        }
        ct.init();
        ct
    };
    let mut electrons = make(-1.0);
    let mut holes = make(1.0);
    let c0 = centroid_x(&electrons.n, nx, 1);
    for _ in 0..200 {
        electrons.step();
        holes.step();
    }
    let de = centroid_x(&electrons.n, nx, 1) - c0;
    let dh = centroid_x(&holes.n, nx, 1) - c0;
    assert!(de < -1.0, "electrons must move against E: {de}");
    assert!(dh > 1.0, "holes must move along E: {dh}");
    assert!((de + dh).abs() < 1e-9 * dh.abs());
}

#[test]
fn drift_velocity_matches_mobility_times_field() {
    let nx = 200;
    let dx = 1e-7;
    let mu = 0.05;
    let e_field = 1.0e4;
    let mut ct = ChargeTransportLbm::new(nx, 1, mu, 300.0, 1.0, dx);
    for x in 0..nx {
        ct.n[x] = bump(x as f64, 60.0, 5.0, 1.0);
        ct.ex[x] = e_field;
    }
    ct.init();
    let c0 = centroid_x(&ct.n, nx, 1);
    let steps = 400;
    for _ in 0..steps {
        ct.step();
    }
    let dist = (centroid_x(&ct.n, nx, 1) - c0) * dx;
    let v = dist / (steps as f64 * ct.time_step());
    let expected = mu * e_field;
    assert!(
        (v - expected).abs() < 0.02 * expected,
        "v = {v}, μE = {expected}"
    );
}

#[test]
fn diffusion_spreading_obeys_einstein_relation() {
    let nx = 256;
    let dx = 1e-7;
    let mu = 0.03;
    let mut ct = ChargeTransportLbm::new(nx, 1, mu, 300.0, 1.2, dx);
    for x in 0..nx {
        ct.n[x] = bump(x as f64, 128.0, 4.0, 1.0);
    }
    ct.init();
    let variance = |n: &[f64]| {
        let c = centroid_x(n, nx, 1);
        let mass: f64 = n.iter().sum();
        n.iter()
            .enumerate()
            .map(|(x, v)| v * (x as f64 - c).powi(2))
            .sum::<f64>()
            / mass
            * dx
            * dx
    };
    let var0 = variance(&ct.n);
    let steps = 500;
    for _ in 0..steps {
        ct.step();
    }
    let t = steps as f64 * ct.time_step();
    let d_measured = (variance(&ct.n) - var0) / (2.0 * t);
    let d_einstein = ct.diffusion_coeff();
    assert!(
        (d_measured - d_einstein).abs() < 0.02 * d_einstein,
        "D = {d_measured:e}, μV_T = {d_einstein:e}"
    );
}

#[test]
fn set_time_step_gives_common_clock() {
    let mut fast = ChargeTransportLbm::new(4, 1, 0.14, 300.0, 1.0, 1e-8);
    let mut slow = ChargeTransportLbm::new(4, 1, 0.045, 300.0, 1.0, 1e-8);
    let dt = fast.time_step();
    slow.set_time_step(dt).expect("valid dt");
    fast.set_time_step(dt).expect("valid dt");
    assert!((slow.time_step() - dt).abs() < 1e-12 * dt);
    // the relaxation rate is kept; the moving weight scales with the mobility
    assert_eq!(fast.omega, 1.0);
    assert_eq!(slow.omega, 1.0);
    assert!((fast.lattice_weight - 1.0 / 6.0).abs() < 1e-12);
    assert!((slow.lattice_weight - (0.045 / 0.14) / 6.0).abs() < 1e-12);
    assert!(slow.set_time_step(-1.0).is_err());
    // a step needing a moving weight above 1/4 is rejected
    assert!(matches!(
        slow.set_time_step(10.0 * dt),
        Err(Error::ParameterOutOfBounds { .. })
    ));
}

#[test]
fn fixed_density_cells_are_held() {
    let nx = 16;
    let mut ct = ChargeTransportLbm::new(nx, 1, 0.1, 300.0, 1.0, 1e-7);
    ct.set_fixed_density(0, 0, 2.0);
    ct.set_fixed_density(nx - 1, 0, 1.0);
    for x in 1..nx - 1 {
        ct.n[x] = 1.0;
    }
    ct.init();
    for _ in 0..5000 {
        ct.step();
    }
    assert_eq!(ct.n[0], 2.0);
    assert_eq!(ct.n[nx - 1], 1.0);
    // steady pure diffusion between two reservoirs → linear profile
    for x in 0..nx {
        let exact = 2.0 - x as f64 / (nx - 1) as f64;
        assert!((ct.n[x] - exact).abs() < 1e-6, "x={x}: {}", ct.n[x]);
    }
}

// ── PlasmaLbm (F2528) ───────────────────────────────────────────────────────

#[test]
fn plasma_clouds_drift_oppositely_in_applied_field() {
    let nx = 128;
    let mut plasma = PlasmaLbm::new(nx, 1, 1e-2, 5e-3, 300.0, 1e-6, 1.0).expect("valid plasma");
    plasma.external_field = [1.0e4, 0.0];
    // dilute, well separated clouds: the self-field is negligible
    for x in 0..nx {
        plasma.electrons.n[x] = bump(x as f64, 84.0, 4.0, 1.0e6);
        plasma.ions.n[x] = bump(x as f64, 34.0, 4.0, 1.0e6);
    }
    let ce0 = centroid_x(&plasma.electrons.n, nx, 1);
    let ci0 = centroid_x(&plasma.ions.n, nx, 1);
    plasma.run(150);
    let de = centroid_x(&plasma.electrons.n, nx, 1) - ce0;
    let di = centroid_x(&plasma.ions.n, nx, 1) - ci0;
    assert!(de < -1.0, "electrons must move against E: {de}");
    assert!(di > 0.5, "ions must move along E: {di}");
    // equal |q|·V_T: displacement ratio equals the mobility ratio
    assert!(((de / di).abs() - 2.0).abs() < 0.05, "ratio {}", de / di);
}

#[test]
fn plasma_feeds_ey_and_separated_clouds_attract() {
    let (nx, ny) = (4, 48);
    let mut plasma = PlasmaLbm::new(nx, ny, 1e-2, 1e-2, 300.0, 1e-6, 1.0).expect("valid plasma");
    for y in 0..ny {
        for x in 0..nx {
            let i = y * nx + x;
            plasma.electrons.n[i] = bump(y as f64, 16.0, 3.0, 1.0e15);
            plasma.ions.n[i] = bump(y as f64, 32.0, 3.0, 1.0e15);
        }
    }
    let sep0 = centroid_y(&plasma.ions.n, nx, ny) - centroid_y(&plasma.electrons.n, nx, ny);
    plasma.run(400);
    let sep1 = centroid_y(&plasma.ions.n, nx, ny) - centroid_y(&plasma.electrons.n, nx, ny);
    let ey_max = plasma
        .electrons
        .ey
        .iter()
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    assert!(ey_max > 0.0, "E_y must be fed to the species");
    assert_eq!(plasma.electrons.ey, plasma.ions.ey);
    assert!(plasma.electrons.ex.iter().all(|v| v.abs() < 1e-9 * ey_max));
    assert!(sep1 < sep0 - 0.5, "clouds must attract: {sep0} → {sep1}");
}

#[test]
fn plasma_rejects_invalid_parameters() {
    assert!(matches!(
        PlasmaLbm::new(4, 4, 0.0, 0.0, 300.0, 1e-6, 1.0),
        Err(Error::ParameterOutOfBounds { .. })
    ));
    assert!(matches!(
        PlasmaLbm::new(0, 4, 1e-2, 1e-4, 300.0, 1e-6, 1.0),
        Err(Error::InvalidGridDimensions { .. })
    ));
}

// ── SpaceChargeLbm (F2526) ──────────────────────────────────────────────────

#[test]
fn space_charge_first_step_preserves_uniform_concentration() {
    let mut sc =
        SpaceChargeLbm::new(6, 5, vec![1.0, -2.0], vec![7e-8, 3e-8], 300.0, 1.0, 1e-6).expect("ok");
    sc.step();
    for s in 0..2 {
        for c in &sc.conc[s] {
            assert!((c - 1.0).abs() < 1e-12, "c = {c}");
        }
    }
}

#[test]
fn space_charge_component_mismatch_is_an_error() {
    let err = SpaceChargeLbm::new(4, 4, vec![1.0, -1.0], vec![1e-8], 300.0, 1.0, 1e-6);
    assert!(matches!(
        err,
        Err(Error::ComponentCountMismatch {
            expected: 2,
            found: 1,
            ..
        })
    ));
    assert!(matches!(
        SpaceChargeLbm::new(4, 4, vec![], vec![], 300.0, 1.0, 1e-6),
        Err(Error::InvalidSpeciesCount { .. })
    ));
}

#[test]
fn space_charge_init_after_editing_concentration() {
    let mut sc = SpaceChargeLbm::new(4, 1, vec![1.0], vec![1e-8], 300.0, 1.0, 1e-6).expect("ok");
    sc.conc[0] = vec![3.0; 4];
    sc.init();
    sc.step();
    for c in &sc.conc[0] {
        assert!((c - 3.0).abs() < 1e-12);
    }
}

#[test]
fn space_charge_drift_follows_valence_mobility_and_field() {
    let (nx, ny) = (96, 96);
    let dx = 1e-6;
    let mu = 5e-8;
    let mut sc =
        SpaceChargeLbm::new(nx, ny, vec![1.0, -1.0, 2.0], vec![mu; 3], 300.0, 1.0, dx).expect("ok");
    // φ decreasing along x and along y: E = (gx, gy) uniform
    let (gx, gy) = (500.0, 250.0);
    for y in 0..ny {
        for x in 0..nx {
            let i = y * nx + x;
            sc.phi[i] = -(gx * x as f64 + gy * y as f64) * dx;
            for s in 0..3 {
                sc.conc[s][i] = bump(x as f64, 48.0, 4.0, 1.0) * bump(y as f64, 48.0, 4.0, 1.0);
            }
        }
    }
    sc.init();
    let cx0 = centroid_x(&sc.conc[0], nx, ny);
    let cy0 = centroid_y(&sc.conc[0], nx, ny);
    let steps = 120;
    for _ in 0..steps {
        sc.step();
    }
    let t = steps as f64 * sc.time_step();
    for (s, z) in [(0, 1.0), (1, -1.0), (2, 2.0)] {
        let vx = (centroid_x(&sc.conc[s], nx, ny) - cx0) * dx / t;
        let vy = (centroid_y(&sc.conc[s], nx, ny) - cy0) * dx / t;
        let (ex_v, ey_v) = (z * mu * gx, z * mu * gy);
        assert!(
            (vx - ex_v).abs() < 0.03 * ex_v.abs(),
            "s={s}: vx {vx} vs {ex_v}"
        );
        assert!(
            (vy - ey_v).abs() < 0.03 * ey_v.abs(),
            "s={s}: vy {vy} vs {ey_v}"
        );
    }
}

// ── SemiconductorLbm (F2525) ────────────────────────────────────────────────

#[test]
fn pn_junction_reaches_self_consistent_equilibrium() {
    let (na, nd, ni) = (1.0e21, 1.0e21, 1.0e16);
    let temperature = 300.0;
    let vt = constants::K_BOLTZMANN * temperature / constants::ELEM_CHARGE;
    let eps = constants::SILICON_EPSILON_R * constants::EPSILON_0;
    let debye = (eps * vt / (constants::ELEM_CHARGE * na)).sqrt();
    let dx = debye / 4.0;
    let nx = 96;
    let mut semi =
        SemiconductorLbm::new(nx, na, nd, ni, temperature, 0.14, 0.045, dx).expect("valid device");
    semi.init_junction();
    let vbi = semi.builtin_voltage();
    let rho0: Vec<f64> = semi.poisson.rho_charge.clone();
    let phi0: Vec<f64> = semi.poisson.phi.clone();
    let (_steps, residual) = semi.solve_equilibrium(100_000, 1e-8);
    assert!(
        residual <= 1e-8,
        "coupled iteration did not settle: {residual:e}"
    );

    // carriers evolved: a depletion region with net space charge formed
    let junction = nx / 2;
    assert!(rho0[junction - 2].abs() < 1e-6 * constants::ELEM_CHARGE * na);
    assert!(semi.poisson.rho_charge[junction - 2] < -0.5 * constants::ELEM_CHARGE * na);
    assert!(semi.poisson.rho_charge[junction + 1] > 0.5 * constants::ELEM_CHARGE * nd);

    // depletion approximation
    let w = (2.0 * eps * (vbi - 2.0 * vt) * (na + nd) / (constants::ELEM_CHARGE * na * nd)).sqrt();
    let e_max_theory = constants::ELEM_CHARGE * na * nd * w / (eps * (na + nd));
    let half_w_cells = (0.5 * w / dx).ceil() as usize;
    let margin = 4 * 4; // four Debye lengths
    let left = junction - half_w_cells - margin;
    let right = junction + half_w_cells + margin;
    assert!(left > 0 && right < nx - 1, "domain too small for the test");
    // the initial (Laplace) potential does not yet carry the built-in drop
    // across the interior; the coupled solve has to build it
    let drop0 = phi0[right] - phi0[left];
    assert!((drop0 - vbi).abs() > 0.2 * vbi, "initial drop {drop0} V");
    let drop = semi.poisson.phi[right] - semi.poisson.phi[left];
    assert!(
        (drop - vbi).abs() < 0.1 * vbi,
        "interior drop {drop} V vs V_bi {vbi} V"
    );
    let e_max = semi
        .electrons
        .ex
        .iter()
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    assert!(
        (e_max - e_max_theory).abs() < 0.1 * e_max_theory,
        "E_max {e_max:e} vs depletion approximation {e_max_theory:e}"
    );
    // quasi-neutral regions: majority carriers back at the doping level
    assert!((semi.holes.n[left] - na).abs() < 0.05 * na);
    assert!((semi.electrons.n[right] - nd).abs() < 0.05 * nd);
    // mass-action law in thermal equilibrium
    for x in [left, junction - 2, junction, junction + 2, right] {
        let np = semi.electrons.n[x] * semi.holes.n[x] / (ni * ni);
        assert!((np - 1.0).abs() < 1e-3, "n·p/n_i² = {np} at x={x}");
    }
}

#[test]
fn semiconductor_rejects_invalid_parameters() {
    assert!(matches!(
        SemiconductorLbm::new(16, -1.0, 1e21, 1e16, 300.0, 0.1, 0.04, 1e-8),
        Err(Error::ParameterOutOfBounds { .. })
    ));
    assert!(matches!(
        SemiconductorLbm::new(2, 1e21, 1e21, 1e16, 300.0, 0.1, 0.04, 1e-8),
        Err(Error::InvalidGridDimensions { .. })
    ));
}
