// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Unit tests of the electrostatic LBM module.

use super::*;
use crate::error::Error;

// ── Poisson ──────────────────────────────────────────────────────────────

/// Uniform charge → quadratic potential profile (1-D slice).
#[test]
fn test_poisson_uniform_charge_quadratic() {
    let nx = 16;
    let ny = 1;
    let mut solver = PoissonLbm::new(nx, ny, 1.0, 1.0);
    // uniform charge density
    for x in 1..nx - 1 {
        solver.set_charge(x, 0, 1.0);
    }
    // Dirichlet BCs
    solver.set_phi(0, 0, 0.0);
    solver.set_phi(nx - 1, 0, 0.0);
    let (iters, _res) = solver
        .solve(2000, 1e-6)
        .expect("grounded problem is solvable");
    assert!(iters <= 2000);
    // positive charge between grounded electrodes → positive, symmetric bump
    let phi_mid = solver.get_phi(nx / 2, 0);
    assert!(phi_mid > 0.0, "φ_mid = {phi_mid}");
    for x in 1..nx - 1 {
        let mirror = solver.get_phi(nx - 1 - x, 0);
        let here = solver.get_phi(x, 0);
        assert!(
            (here - mirror).abs() <= 1e-9 * phi_mid,
            "asymmetric at x={x}"
        );
    }
}

/// Poisson solver: zero charge → potential stays at BCs.
#[test]
fn test_poisson_zero_charge() {
    let nx = 8;
    let mut solver = PoissonLbm::new(nx, 1, 1.0, 1.0);
    solver.set_phi(0, 0, 1.0);
    solver.set_phi(nx - 1, 0, 0.0);
    let (_i, _r) = solver
        .solve(5000, 1e-12)
        .expect("Dirichlet problem is solvable");
    // Dirichlet values are held exactly and the profile is linear.
    assert_eq!(solver.get_phi(0, 0), 1.0);
    assert_eq!(solver.get_phi(nx - 1, 0), 0.0);
    for x in 0..nx {
        let exact = 1.0 - x as f64 / (nx - 1) as f64;
        assert!((solver.get_phi(x, 0) - exact).abs() < 1e-6, "x={x}");
    }
}

/// Convergence: residual decreases with more iterations.
#[test]
fn test_poisson_bgk_convergence() {
    let nx = 12;
    let mut solver = PoissonLbm::new(nx, 1, 1.0, 1.0);
    for x in 2..nx - 2 {
        solver.set_charge(x, 0, 0.5);
    }
    solver.set_phi(0, 0, 0.0);
    solver.set_phi(nx - 1, 0, 0.0);
    let (_i, res1) = solver.solve(100, 0.0).expect("solvable");
    let (_i2, res2) = solver.solve(500, 0.0).expect("solvable");
    // more iterations → smaller residual (or already converged)
    assert!(res2 <= res1 + 1e-10);
}

/// Poisson: residual after full solve should be very small.
#[test]
fn test_poisson_residual_small() {
    let nx = 10;
    let mut solver = PoissonLbm::new(nx, 1, 1.0, 1.0);
    solver.set_charge(5, 0, 1.0);
    solver.set_phi(0, 0, 0.0);
    solver.set_phi(nx - 1, 0, 0.0);
    let (_i, res) = solver.solve(5000, 1e-9).expect("solvable");
    assert!(res < 1e-3, "residual {:.6} not small", res);
}

// ── Electric field ────────────────────────────────────────────────────────

/// E = −∇φ: linear potential → uniform field.
#[test]
fn test_electric_field_linear_phi() {
    let nx = 5;
    let ny = 1;
    let dx = 1.0;
    let mut egrid = ElectricFieldGrid::new(nx, ny, dx, dx);
    // φ(x) = x  →  Ex = -1, Ey = 0
    let phi: Vec<f64> = (0..nx).map(|x| x as f64).collect();
    egrid.compute_from_phi(&phi);
    // interior cells
    for x in 1..nx - 1 {
        let (ex, ey) = egrid.get(x, 0);
        assert!((ex - (-1.0)).abs() < 1e-10, "Ex wrong at x={}", x);
        assert!(ey.abs() < 1e-10);
    }
}

/// E = −∇φ: constant potential → zero field.
#[test]
fn test_electric_field_constant_phi() {
    let nx = 6;
    let ny = 2;
    let mut egrid = ElectricFieldGrid::new(nx, ny, 1.0, 1.0);
    let phi = vec![5.0_f64; nx * ny];
    egrid.compute_from_phi(&phi);
    for y in 0..ny {
        for x in 0..nx {
            let (ex, ey) = egrid.get(x, y);
            assert!(ex.abs() < 1e-12);
            assert!(ey.abs() < 1e-12);
        }
    }
}

/// Magnitude is sqrt(Ex²+Ey²).
#[test]
fn test_electric_field_magnitude() {
    let nx = 5;
    let ny = 5;
    let mut egrid = ElectricFieldGrid::new(nx, ny, 1.0, 1.0);
    let phi: Vec<f64> = (0..nx * ny).map(|i| i as f64).collect();
    egrid.compute_from_phi(&phi);
    let m = egrid.magnitude(2, 2);
    let (ex, ey) = egrid.get(2, 2);
    assert!((m - (ex * ex + ey * ey).sqrt()).abs() < 1e-12);
}

// ── Dielectric ────────────────────────────────────────────────────────────

/// Isotropic dielectric: energy density = ½ ε₀ ε_r E².
#[test]
fn test_dielectric_energy_density() {
    let eps_r = 4.0;
    let mat = DielectricMaterial::isotropic(eps_r);
    let ex = 1000.0_f64; // V/m
    let ey = 0.0;
    let u = mat.energy_density(ex, ey);
    let expected = 0.5 * constants::EPSILON_0 * eps_r * ex * ex;
    assert!((u - expected).abs() < 1e-25, "energy density mismatch");
}

/// Polarisation: P = ε₀ (ε_r − 1) E.
#[test]
fn test_dielectric_polarisation() {
    let eps_r = 3.0;
    let mat = DielectricMaterial::isotropic(eps_r);
    let ex = 500.0;
    let (px, _py) = mat.polarisation(ex, 0.0);
    let expected = constants::EPSILON_0 * (eps_r - 1.0) * ex;
    assert!((px - expected).abs() < 1e-25);
}

/// Displacement field D = ε₀ ε_r E.
#[test]
fn test_dielectric_displacement() {
    let mat = DielectricMaterial::isotropic(2.5);
    let (dx, _dy) = mat.displacement(200.0, 0.0);
    let expected = constants::EPSILON_0 * 2.5 * 200.0;
    assert!((dx - expected).abs() < 1e-25);
}

/// Anisotropic dielectric energy density.
#[test]
fn test_dielectric_anisotropic() {
    let mat = DielectricMaterial::anisotropic(2.0, 5.0);
    let ex = 100.0;
    let ey = 100.0;
    let u = mat.energy_density(ex, ey);
    let expected =
        0.5 * (constants::EPSILON_0 * 2.0 * ex * ex + constants::EPSILON_0 * 5.0 * ey * ey);
    assert!((u - expected).abs() < 1e-24);
}

// ── Einstein relation ─────────────────────────────────────────────────────

/// D / μ = V_T = kT/q.
#[test]
fn test_einstein_relation() {
    let temperature = 300.0;
    let mobility = 0.1; // m²/V·s
    let transport = ChargeTransportLbm::new(4, 1, mobility, temperature, 1.0, 1e-6);
    let vt_expected = constants::K_BOLTZMANN * temperature / constants::ELEM_CHARGE;
    let d = transport.diffusion_coeff();
    assert!((d / mobility - vt_expected).abs() < 1e-12);
}

/// Einstein relation holds at multiple temperatures.
#[test]
fn test_einstein_relation_temperatures() {
    for t in [200.0, 300.0, 400.0, 600.0] {
        let mu = 0.05;
        let transport = ChargeTransportLbm::new(4, 1, mu, t, 1.0, 1e-6);
        let vt = constants::K_BOLTZMANN * t / constants::ELEM_CHARGE;
        assert!((transport.diffusion_coeff() / mu - vt).abs() < 1e-12);
    }
}

// ── Semiconductor / p-n junction ──────────────────────────────────────────

/// Built-in voltage > 0 for a p-n junction.
#[test]
fn test_pn_junction_builtin_positive() {
    let ni = 1.5e10_f64; // cm⁻³ in SI: 1.5e16 m⁻³; using relative units
    let na = 1e17_f64;
    let nd = 1e17_f64;
    let semi =
        SemiconductorLbm::new(16, na, nd, ni, 300.0, 0.15, 0.045, 1e-7).expect("valid device");
    let vbi = semi.builtin_voltage();
    assert!(vbi > 0.0, "V_bi = {:.6} V should be positive", vbi);
}

/// Built-in voltage increases with doping.
#[test]
fn test_pn_junction_vbi_doping_dependence() {
    let ni = 1.5e10_f64;
    let s1 = SemiconductorLbm::new(8, 1e15, 1e15, ni, 300.0, 0.1, 0.04, 1e-7).expect("valid");
    let s2 = SemiconductorLbm::new(8, 1e17, 1e17, ni, 300.0, 0.1, 0.04, 1e-7).expect("valid");
    assert!(s2.builtin_voltage() > s1.builtin_voltage());
}

/// Junction steps keep the ohmic contacts pinned and all fields physical.
#[test]
fn test_semiconductor_step_runs() {
    let nx = 32;
    let (na, nd, ni) = (1e21, 1e21, 1e16);
    let mut semi =
        SemiconductorLbm::new(nx, na, nd, ni, 300.0, 0.14, 0.045, 2.5e-8).expect("valid device");
    semi.init_junction();
    let phi_p = semi.poisson.get_phi(0, 0);
    let phi_n = semi.poisson.get_phi(nx - 1, 0);
    assert!((phi_n - phi_p - semi.builtin_voltage()).abs() < 1e-9);
    for _ in 0..10 {
        semi.step();
    }
    assert_eq!(semi.poisson.get_phi(0, 0), phi_p);
    assert_eq!(semi.poisson.get_phi(nx - 1, 0), phi_n);
    assert!((semi.holes.n[0] - na).abs() < 1e-9 * na);
    assert!((semi.electrons.n[nx - 1] - nd).abs() < 1e-9 * nd);
    for x in 0..nx {
        assert!(semi.poisson.phi[x].is_finite());
        assert!(semi.electrons.n[x] > 0.0 && semi.holes.n[x] > 0.0, "x={x}");
    }
}

// ── Space charge ──────────────────────────────────────────────────────────

/// Two equal-and-opposite species → near-zero net charge density.
#[test]
fn test_space_charge_neutrality() {
    let nx = 8;
    let ny = 1;
    let valence = vec![1.0, -1.0];
    let mobility = vec![1e-8, 1e-8];
    let mut sc =
        SpaceChargeLbm::new(nx, ny, valence, mobility, 300.0, 1.0, 1e-6).expect("valid solver");
    // Run a few steps with zero potential (no drift); the constructor
    // initialises both species at equilibrium.
    for _ in 0..5 {
        sc.step();
    }
    let rho = sc.mean_charge_density();
    assert!(
        rho < 1e-9 * constants::FARADAY,
        "mean |ρ| = {:.6} should be near-zero with equal species",
        rho
    );
}

/// charge_density returns correct sign for positive species only.
#[test]
fn test_space_charge_density_sign() {
    let valence = vec![1.0];
    let mobility = vec![1e-8];
    let mut sc = SpaceChargeLbm::new(4, 1, valence, mobility, 300.0, 1.0, 1e-6).expect("valid");
    sc.conc[0][0] = 2.0;
    let rho = sc.charge_density(0);
    // 2 mol/m³ of a monovalent cation carries 2·F coulombs per m³.
    assert!((rho - 2.0 * constants::FARADAY).abs() < 1e-9 * rho);
}

// ── Coulomb explosion ─────────────────────────────────────────────────────

/// Two like charges repel: kinetic energy increases over time.
#[test]
fn test_coulomb_explosion_repulsion() {
    let x = vec![0.0, 1e-9];
    let y = vec![0.0, 0.0];
    let charge = vec![constants::ELEM_CHARGE; 2];
    let mass = vec![1.67e-27_f64; 2]; // proton mass
    let mut ce = CoulombExplosion::new(x, y, charge, mass, 1e-12);
    let ke0 = ce.kinetic_energy();
    for _ in 0..100 {
        ce.step(1e-18); // 1 as step
    }
    let ke1 = ce.kinetic_energy();
    assert!(ke1 > ke0, "KE should increase: {:.6} → {:.6}", ke0, ke1);
}

/// Energy conservation: KE + PE ≈ const (no dissipation).
#[test]
fn test_coulomb_explosion_energy_conservation() {
    let x = vec![0.0, 2e-9];
    let y = vec![0.0, 0.0];
    let q = constants::ELEM_CHARGE;
    let charge = vec![q, q];
    let mass = vec![1.67e-27_f64; 2];
    let mut ce = CoulombExplosion::new(x, y, charge, mass, 1e-13);
    let e0 = ce.kinetic_energy() + ce.potential_energy();
    for _ in 0..20 {
        ce.step(1e-19);
    }
    let e1 = ce.kinetic_energy() + ce.potential_energy();
    let rel_err = ((e1 - e0) / e0.abs()).abs();
    assert!(
        rel_err < 0.01,
        "energy conservation violated: rel_err = {:.6}",
        rel_err
    );
}

/// Opposite charges attract: particles move closer.
#[test]
fn test_coulomb_attraction() {
    let x0 = 0.0_f64;
    let x1 = 2e-9_f64;
    let x = vec![x0, x1];
    let y = vec![0.0, 0.0];
    let charge = vec![constants::ELEM_CHARGE, -constants::ELEM_CHARGE];
    let mass = vec![1.67e-27_f64; 2];
    let mut ce = CoulombExplosion::new(x, y, charge, mass, 1e-12);
    for _ in 0..50 {
        ce.step(1e-18);
    }
    let dist = (ce.x[1] - ce.x[0]).abs();
    assert!(dist < x1 - x0, "particles should attract");
}

/// Softening prevents force blow-up.
#[test]
fn test_coulomb_softening() {
    let x = vec![0.0, 1e-15]; // very close
    let y = vec![0.0, 0.0];
    let charge = vec![constants::ELEM_CHARGE; 2];
    let mass = vec![1.67e-27_f64; 2];
    let mut ce = CoulombExplosion::new(x, y, charge, mass, 1e-12);
    ce.step(1e-20);
    assert!(ce.vx[0].is_finite());
    assert!(ce.vx[1].is_finite());
}

// ── Plasma LBM ───────────────────────────────────────────────────────────

/// Plasma LBM: density stays positive after several steps.
#[test]
fn test_plasma_lbm_density_positive() {
    let nx = 8;
    let ny = 1;
    let mut plasma = PlasmaLbm::new(nx, ny, 1e-2, 1e-4, 300.0, 1e-6, 1.0).expect("valid plasma");
    // Set initial uniform density
    let n0 = 1e20_f64;
    for i in 0..nx {
        plasma.electrons.n[i] = n0;
        plasma.ions.n[i] = n0;
    }
    plasma.run(20);
    for i in 0..nx {
        assert!(plasma.electrons.n[i] >= 0.0);
        assert!(plasma.ions.n[i] >= 0.0);
    }
}

/// Plasma LBM: Poisson field updates without NaN.
#[test]
fn test_plasma_lbm_no_nan() {
    let mut plasma = PlasmaLbm::new(6, 1, 1e-2, 1e-4, 300.0, 1e-6, 1.0).expect("valid plasma");
    for i in 0..6 {
        plasma.electrons.n[i] = 1e18;
        plasma.ions.n[i] = 1e18;
    }
    plasma.run(10);
    for phi in &plasma.poisson.phi {
        assert!(phi.is_finite());
    }
    // a uniform neutral plasma stays uniform and neutral
    for i in 0..6 {
        assert!((plasma.electrons.n[i] - 1e18).abs() < 1e-6 * 1e18);
        assert!((plasma.ions.n[i] - 1e18).abs() < 1e-6 * 1e18);
    }
}

// ── Charge transport ──────────────────────────────────────────────────────

/// Carrier density remains finite after many drift-diffusion steps.
#[test]
fn test_charge_transport_finite() {
    let nx = 16;
    let ny = 1;
    let mut ct = ChargeTransportLbm::new(nx, ny, 0.1, 300.0, 1.0, 1e-7);
    for i in 0..nx {
        ct.n[i] = 1e20;
    }
    // Uniform field
    for i in 0..nx {
        ct.ex[i] = 1e4;
    }
    ct.run(50);
    for n in &ct.n {
        assert!(n.is_finite());
    }
}

/// No drift (zero field) → carrier profile stays uniform.
#[test]
fn test_charge_transport_no_drift_uniform() {
    let nx = 8;
    let mut ct = ChargeTransportLbm::new(nx, 1, 0.1, 300.0, 1.0, 1e-7);
    let n0 = 1.0_f64;
    for i in 0..nx {
        ct.n[i] = n0;
    }
    ct.run(20);
    for n in &ct.n {
        assert!((n - n0).abs() < 1e-12, "n = {:.6}", n);
    }
}

// ── General / constants ───────────────────────────────────────────────────

/// Physical constants have correct orders of magnitude.
#[test]
fn test_constants_magnitude() {
    let eps0 = constants::EPSILON_0;
    let elem = constants::ELEM_CHARGE;
    let kb = constants::K_BOLTZMANN;
    let vt = constants::VT_300K;
    assert!(eps0 > 8e-12 && eps0 < 9e-12);
    assert!(elem > 1.6e-19 && elem < 1.7e-19);
    assert!(kb > 1.3e-23 && kb < 1.4e-23);
    assert!((vt - 0.02585).abs() < 0.001);
}

/// V_T at 300 K matches kT/q.
#[test]
fn test_thermal_voltage_300k() {
    let vt = constants::K_BOLTZMANN * 300.0 / constants::ELEM_CHARGE;
    assert!((vt - constants::VT_300K).abs() < 1e-5);
}

/// Poisson: a charged periodic 2×2 grid without Dirichlet cells has no
/// solution and is reported as such.
#[test]
fn test_poisson_2x2() {
    let mut solver = PoissonLbm::new(2, 2, 1.0, 1.0);
    solver.set_charge(0, 0, 1.0);
    let result = solver.solve(10, 1e-6);
    assert!(
        matches!(result, Err(Error::ConvergenceFailure { .. })),
        "{result:?}"
    );
}

/// ElectricFieldGrid: 1×1 grid does not panic.
#[test]
fn test_electric_field_1x1() {
    let mut eg = ElectricFieldGrid::new(1, 1, 1.0, 1.0);
    let phi = vec![3.0_f64];
    eg.compute_from_phi(&phi);
    let (ex, ey) = eg.get(0, 0);
    assert!(ex.is_finite());
    assert!(ey.is_finite());
}
