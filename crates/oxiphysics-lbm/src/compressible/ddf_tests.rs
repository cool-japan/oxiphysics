// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the compressible equilibria and the
//! double-distribution solvers.

use super::*;

fn moments(feq: &[f64; 9]) -> (f64, f64, f64) {
    let mut rho = 0.0;
    let mut mx = 0.0;
    let mut my = 0.0;
    for q in 0..9 {
        rho += feq[q];
        mx += f64::from(EX9[q]) * feq[q];
        my += f64::from(EY9[q]) * feq[q];
    }
    (rho, mx, my)
}

fn third_moment_xyy(feq: &[f64; 9]) -> f64 {
    (0..9)
        .map(|q| {
            let cx = f64::from(EX9[q]);
            let cy = f64::from(EY9[q]);
            feq[q] * cx * cy * cy
        })
        .sum()
}

#[test]
fn equilibrium_compressible_conserves_mass_and_momentum() {
    let (rho, ux, uy) = (1.3, 0.1, 0.05);
    let feq = equilibrium_compressible(rho, ux, uy, LATTICE_CS2);
    let (m0, mx, my) = moments(&feq);
    assert!((m0 - rho).abs() < 1e-13);
    assert!((mx - rho * ux).abs() < 1e-13);
    assert!((my - rho * uy).abs() < 1e-13);
}

#[test]
fn equilibrium_compressible_pressure_tensor() {
    let (rho, ux, uy) = (0.9, 0.12, -0.07);
    let feq = equilibrium_compressible(rho, ux, uy, LATTICE_CS2);
    let mut pxx = 0.0;
    let mut pxy = 0.0;
    let mut pyy = 0.0;
    for q in 0..9 {
        let cx = f64::from(EX9[q]);
        let cy = f64::from(EY9[q]);
        pxx += feq[q] * cx * cx;
        pxy += feq[q] * cx * cy;
        pyy += feq[q] * cy * cy;
    }
    assert!((pxx - rho * (LATTICE_CS2 + ux * ux)).abs() < 1e-13);
    assert!((pxy - rho * ux * uy).abs() < 1e-13);
    assert!((pyy - rho * (LATTICE_CS2 + uy * uy)).abs() < 1e-13);
}

#[test]
fn equilibrium_compressible_matches_third_order_correction() {
    let corr = MaCorrection::new(3, LATTICE_CS2.sqrt());
    let a = equilibrium_compressible(1.1, 0.15, -0.08, LATTICE_CS2);
    let b = corr.corrected_feq(1.1, 0.15, -0.08);
    for q in 0..9 {
        assert!((a[q] - b[q]).abs() < 1e-14, "q={q}: {} vs {}", a[q], b[q]);
    }
}

#[test]
fn third_order_moment_discriminates_ma_order() {
    let (rho, ux, uy) = (1.0, 0.1, 0.05);
    let cs = LATTICE_CS2.sqrt();
    let exact = rho * (LATTICE_CS2 * ux + ux * uy * uy);
    let f3 = MaCorrection::new(3, cs).corrected_feq(rho, ux, uy);
    assert!((third_moment_xyy(&f3) - exact).abs() < 1e-13);
    let f2 = MaCorrection::new(2, cs).corrected_feq(rho, ux, uy);
    let second = third_moment_xyy(&f2);
    assert!((second - rho * LATTICE_CS2 * ux).abs() < 1e-13);
    assert!((second - exact).abs() > 1e-5);
}

#[test]
fn energy_equilibrium_moments() {
    let (rho, ux, uy, temp, gamma) = (1.2, 0.05, -0.03, 0.4, 1.4);
    let feq = equilibrium_compressible(rho, ux, uy, LATTICE_CS2);
    let geq = energy_equilibrium(&feq, rho, ux, uy, temp, gamma);
    let e = 0.5 * (ux * ux + uy * uy) + temp / (gamma - 1.0);
    let (g0, gx, gy) = moments(&geq);
    let p = rho * temp;
    assert!((g0 - rho * e).abs() < 1e-13);
    assert!((gx - (rho * e + p) * ux).abs() < 1e-13);
    assert!((gy - (rho * e + p) * uy).abs() < 1e-13);
}

#[test]
fn compressible_lbm_uses_lattice_sound_speed() {
    let sim = CompressibleLbm::new(4, 4, 1.5, 1.4, 3);
    assert!((sim.cs * sim.cs - LATTICE_CS2).abs() < 1e-15);
    // Rest state: p = ρT = ρ cs², adiabatic sound speed sqrt(γ T).
    assert!((sim.pressure_at(0, 0) - LATTICE_CS2).abs() < 1e-14);
    let feq_sum: f64 = sim.f[0][0].iter().sum();
    assert!((feq_sum - 1.0).abs() < 1e-14);
}

/// Smooth periodic perturbation of density, velocity and temperature.
fn perturb<F>(nx: usize, ny: usize, mut set: F)
where
    F: FnMut(usize, usize, f64, [f64; 2], f64),
{
    let two_pi = 2.0 * std::f64::consts::PI;
    for y in 0..ny {
        for x in 0..nx {
            let sx = (two_pi * x as f64 / nx as f64).sin();
            let cy = (two_pi * y as f64 / ny as f64).cos();
            let rho = 1.0 + 0.05 * sx * cy;
            let u = [0.1 + 0.02 * cy, 0.05 - 0.02 * sx];
            let temp = LATTICE_CS2 * (1.0 + 0.1 * sx);
            set(x, y, rho, u, temp);
        }
    }
}

#[test]
fn compressible_lbm_conserves_mass_momentum_energy() {
    for ma_order in [2usize, 3] {
        let mut sim = CompressibleLbm::new(16, 12, 1.4, 1.4, ma_order);
        let (nx, ny) = (sim.nx, sim.ny);
        perturb(nx, ny, |x, y, rho, u, temp| {
            sim.rho[y][x] = rho;
            sim.u[y][x] = u;
            sim.temp[y][x] = temp;
        });
        sim.initialize();
        sim.update_macros();
        let (m0, p0, e0) = sim.conserved_totals();
        assert!(p0[0].abs() > 1.0, "test must start with net momentum");
        for _ in 0..100 {
            sim.step();
        }
        let (m1, p1, e1) = sim.conserved_totals();
        assert!((m1 - m0).abs() < 1e-11 * m0, "mass drift {}", m1 - m0);
        assert!((p1[0] - p0[0]).abs() < 1e-11, "px drift {}", p1[0] - p0[0]);
        assert!((p1[1] - p0[1]).abs() < 1e-11, "py drift {}", p1[1] - p0[1]);
        assert!((e1 - e0).abs() < 1e-11 * e0, "energy drift {}", e1 - e0);
        assert!(sim.rho.iter().flatten().all(|r| r.is_finite() && *r > 0.0));
    }
}

#[test]
fn compressible_lbm_gamma_changes_dynamics() {
    let run = |gamma: f64| {
        let mut sim = CompressibleLbm::new(16, 4, 1.2, gamma, 3);
        for x in 0..16 {
            for y in 0..4 {
                sim.temp[y][x] = LATTICE_CS2 * if (4..8).contains(&x) { 1.2 } else { 1.0 };
            }
        }
        sim.initialize();
        sim.update_macros();
        for _ in 0..20 {
            sim.step();
        }
        sim.temp[0][6]
    };
    let t_air = run(1.4);
    let t_mono = run(5.0 / 3.0);
    assert!((t_air - t_mono).abs() > 1e-5, "γ has no effect");
}

#[test]
fn ns_uniform_rest_state_keeps_temperature() {
    let mut sim = NsCompressibleLbm::new(8, 8, 1.2, 1.4, 0.71, 1.4);
    for row in sim.temp.iter_mut() {
        for t in row.iter_mut() {
            *t = 1.0;
        }
    }
    sim.initialize();
    for _ in 0..100 {
        sim.step();
    }
    for row in &sim.temp {
        for &t in row {
            assert!((t - 1.0).abs() < 1e-12, "T drifted to {t}");
        }
    }
    for row in &sim.u {
        for v in row {
            assert!(v[0].abs() < 1e-13 && v[1].abs() < 1e-13);
        }
    }
}

#[test]
fn ns_conserves_totals_with_nonuniform_temperature() {
    let mut sim = NsCompressibleLbm::new(16, 16, 1.3, 1.1, 0.71, 1.4);
    let (nx, ny) = (sim.nx, sim.ny);
    perturb(nx, ny, |x, y, rho, u, temp| {
        sim.rho[y][x] = rho;
        sim.u[y][x] = u;
        sim.temp[y][x] = temp;
    });
    sim.initialize();
    sim.update_macros();
    let (m0, p0, e0) = sim.conserved_totals();
    for _ in 0..100 {
        sim.step();
    }
    let (m1, p1, e1) = sim.conserved_totals();
    assert!((m1 - m0).abs() < 1e-11 * m0);
    assert!((p1[0] - p0[0]).abs() < 1e-11);
    assert!((p1[1] - p0[1]).abs() < 1e-11);
    assert!((e1 - e0).abs() < 1e-11 * e0);
}

#[test]
fn ns_temperature_gradient_drives_flow() {
    let mut sim = NsCompressibleLbm::new(32, 4, 1.2, 1.2, 1.0, 1.4);
    for y in 0..4 {
        for x in 0..32 {
            let s = (2.0 * std::f64::consts::PI * x as f64 / 32.0).sin();
            sim.temp[y][x] = LATTICE_CS2 * (1.0 + 0.05 * s);
        }
    }
    sim.initialize();
    sim.update_macros();
    for _ in 0..5 {
        sim.step();
    }
    // Pressure is highest at x = 8 (sin peak); fluid accelerates away from it:
    // u_x > 0 on its right (x ≈ 16) and u_x < 0 on its left (x ≈ 0).
    assert!(sim.u[1][16][0] > 1e-5, "u_x(16) = {}", sim.u[1][16][0]);
    assert!(sim.u[1][0][0] < -1e-5, "u_x(0) = {}", sim.u[1][0][0]);
}

/// Shear-wave amplitude: projection of u_x onto sin(k y).
fn shear_amplitude(sim: &EulerLbm) -> f64 {
    let k = 2.0 * std::f64::consts::PI / sim.ny as f64;
    let mut acc = 0.0;
    for y in 0..sim.ny {
        for x in 0..sim.nx {
            acc += sim.u[y][x][0] * (k * y as f64).sin();
        }
    }
    2.0 * acc / (sim.nx * sim.ny) as f64
}

fn shear_wave(omega: f64) -> EulerLbm {
    let mut sim = EulerLbm::new(4, 32, omega, 1.4);
    let k = 2.0 * std::f64::consts::PI / sim.ny as f64;
    for y in 0..sim.ny {
        for x in 0..sim.nx {
            sim.u[y][x] = [1e-3 * (k * y as f64).sin(), 0.0];
        }
    }
    sim.initialize();
    sim.update_macros();
    sim
}

#[test]
fn euler_shear_wave_decay_matches_viscosity() {
    let mut sim = shear_wave(1.2);
    let k = 2.0 * std::f64::consts::PI / sim.ny as f64;
    for _ in 0..20 {
        sim.step();
    }
    let a1 = shear_amplitude(&sim);
    let steps = 200;
    for _ in 0..steps {
        sim.step();
    }
    let a2 = shear_amplitude(&sim);
    let nu_measured = (a1 / a2).ln() / (k * k * steps as f64);
    let nu = sim.viscosity();
    assert!(
        (nu_measured - nu).abs() / nu < 0.02,
        "ν measured {nu_measured} vs {nu}"
    );
}

#[test]
fn euler_default_omega_is_nearly_inviscid() {
    let mut sim = shear_wave(EulerLbm::DEFAULT_OMEGA);
    let a0 = shear_amplitude(&sim);
    for _ in 0..100 {
        sim.step();
    }
    let a1 = shear_amplitude(&sim);
    assert!(a1 / a0 > 0.99, "amplitude ratio {}", a1 / a0);
    assert!(a1 / a0 < 1.0 + 1e-9);
}

#[test]
fn euler_energy_population_is_evolved_and_conserved() {
    let mut sim = EulerLbm::with_default_omega(16, 16, 1.4);
    let (nx, ny) = (sim.nx, sim.ny);
    perturb(nx, ny, |x, y, rho, u, temp| {
        sim.rho[y][x] = rho;
        sim.u[y][x] = u;
        sim.temp[y][x] = temp;
    });
    sim.initialize();
    sim.update_macros();
    let g_init = sim.g.clone();
    let (m0, p0, e0) = sim.conserved_totals();
    for _ in 0..50 {
        sim.step();
    }
    let (m1, p1, e1) = sim.conserved_totals();
    assert!((m1 - m0).abs() < 1e-11 * m0);
    assert!((p1[0] - p0[0]).abs() < 1e-11 && (p1[1] - p0[1]).abs() < 1e-11);
    assert!((e1 - e0).abs() < 1e-11 * e0);
    // g actually streams: the field changes.
    assert!((sim.g[3][5][1] - g_init[3][5][1]).abs() > 1e-8);
    // energy field mirrors Σg.
    let sum: f64 = sim.g[3][5].iter().sum();
    assert!((sim.energy[3][5] - sum).abs() < 1e-14);
}

#[test]
fn supersonic_inlet_feq_carries_inflow_moments() {
    let inlet = SupersonicInlet::new(0.5, 0.3, 1.2, 1.2 * LATTICE_CS2, 1.4);
    let mut f = vec![vec![vec![0.0f64; 9]; 3]; 2];
    inlet.apply(&mut f, 2);
    let cell: [f64; 9] = std::array::from_fn(|q| f[1][0][q]);
    let (m0, mx, my) = moments(&cell);
    let [ux, uy] = inlet.velocity();
    assert!((m0 - 1.2).abs() < 1e-13);
    assert!((mx - 1.2 * ux).abs() < 1e-13);
    assert!((my - 1.2 * uy).abs() < 1e-13);
}
