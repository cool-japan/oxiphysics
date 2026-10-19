//! Regression tests for the forced / turbulent / time-step-aware stepping of
//! [`LbmSimulation`](super::LbmSimulation).
//!
//! They pin down the physics of the Guo forcing scheme (mass conservation and
//! exact momentum injection), the 3D forcing / Smagorinsky / wall support, the
//! meaning of the `dt` argument of `step`, and the ensemble perturbation.

use super::step_kernels::{
    CollisionSetup, LatticeTable, MAX_Q, build_wall_links_3d, collide_cell, smagorinsky_tau,
};
use super::types::{
    EnsembleRunner, FullStepLoop, GuoForceParams, LbmConfig, LbmSimulation, SimGrid,
    channel_walls_3d,
};
use crate::error::Error;
use crate::lattice::{CS2, D2Q9_VELOCITIES, Lattice, LatticeType};
use std::f64::consts::PI;

const ALL_LATTICES: [LatticeType; 3] = [LatticeType::D2Q9, LatticeType::D3Q19, LatticeType::D3Q27];

fn table_for(lattice_type: LatticeType) -> LatticeTable {
    LatticeTable::from_lattice(&Lattice::new(lattice_type))
}

fn config_for(lattice_type: LatticeType, n: usize, viscosity: f64) -> LbmConfig {
    match lattice_type {
        LatticeType::D2Q9 => LbmConfig::d2q9(n, n, viscosity),
        LatticeType::D3Q19 => LbmConfig::d3q19(n, n, n, viscosity),
        LatticeType::D3Q27 => LbmConfig::d3q27(n, n, n, viscosity),
    }
}

fn total_mass(sim: &LbmSimulation) -> f64 {
    sim.get_density_field().iter().sum()
}

// ---------------------------------------------------------------------------
// Guo force term
// ---------------------------------------------------------------------------

#[test]
fn guo_force_term_conserves_mass_and_injects_prescribed_momentum() {
    let omega = 1.3;
    let prefactor = 1.0 - 0.5 * omega;
    for lattice_type in ALL_LATTICES {
        let table = table_for(lattice_type);
        // The 2D lattice has no z direction: its force and velocity live in
        // the xy-plane.
        let (force, u) = if lattice_type == LatticeType::D2Q9 {
            ([3.0e-4, -2.0e-4, 0.0], [0.04, -0.03, 0.0])
        } else {
            ([3.0e-4, -2.0e-4, 5.0e-4], [0.04, -0.03, 0.02])
        };
        let guo = GuoForceParams::new(force, omega);
        let mut sum = 0.0;
        let mut first = [0.0f64; 3];
        let mut second = [[0.0f64; 3]; 3];
        for (&w, &c) in table.weights.iter().zip(&table.velocities) {
            let term = guo.force_term(w, c, u);
            sum += term;
            for a in 0..3 {
                first[a] += term * c[a];
                for b in 0..3 {
                    second[a][b] += term * c[a] * c[b];
                }
            }
        }
        assert!(
            sum.abs() < 1e-16,
            "{lattice_type:?}: sum of the Guo source must vanish, got {sum:e}"
        );
        for a in 0..3 {
            assert!(
                (first[a] - prefactor * force[a]).abs() < 1e-16,
                "{lattice_type:?}: first moment {a}: {} vs {}",
                first[a],
                prefactor * force[a]
            );
            for b in 0..3 {
                let expected = prefactor * (u[a] * force[b] + u[b] * force[a]);
                assert!(
                    (second[a][b] - expected).abs() < 1e-16,
                    "{lattice_type:?}: second moment ({a},{b}): {} vs {expected}",
                    second[a][b]
                );
            }
        }
    }
}

#[test]
fn d2q9_force_term_sums_to_zero_for_moving_fluid() {
    // The earlier implementation omitted the -u/cs^2 part, so the sum was
    // (1 - omega/2) u.F / cs^2 != 0 as soon as the fluid moved.
    let guo = GuoForceParams::new([1.0e-4, 5.0e-5, 0.0], 1.2);
    let (ux, uy) = (0.05, 0.02);
    let sum: f64 = (0..9).map(|i| guo.d2q9_force_term(i, ux, uy)).sum();
    assert!(sum.abs() < 1e-17, "sum = {sum:e}");
    let momentum_x: f64 = (0..9)
        .map(|i| guo.d2q9_force_term(i, ux, uy) * f64::from(D2Q9_VELOCITIES[i][0]))
        .sum();
    assert!(
        (momentum_x - 0.4 * 1.0e-4).abs() < 1e-17,
        "x-momentum = {momentum_x:e}"
    );
}

#[test]
fn generic_force_term_matches_d2q9_force_term() {
    let guo = GuoForceParams::new([2.0e-4, -1.0e-4, 0.0], 0.9);
    let table = table_for(LatticeType::D2Q9);
    for i in 0..9 {
        let generic = guo.force_term(table.weights[i], table.velocities[i], [0.03, 0.01, 0.0]);
        let specific = guo.d2q9_force_term(i, 0.03, 0.01);
        assert!((generic - specific).abs() < 1e-18, "direction {i}");
    }
}

#[test]
fn opposite_table_is_consistent_with_velocities() {
    for lattice_type in ALL_LATTICES {
        let table = table_for(lattice_type);
        for i in 0..table.q {
            let j = table.opposite[i];
            for d in 0..3 {
                assert_eq!(
                    table.velocities[j][d], -table.velocities[i][d],
                    "{lattice_type:?}: direction {i} vs opposite {j}"
                );
            }
            assert_eq!(table.opposite[j], i);
            assert!((table.weights[j] - table.weights[i]).abs() < 1e-18);
        }
    }
}

// ---------------------------------------------------------------------------
// Forced periodic flow: exact momentum bookkeeping
// ---------------------------------------------------------------------------

#[test]
fn periodic_2d_forced_flow_conserves_mass() {
    let mut cfg = LbmConfig::d2q9(32, 32, 1.0 / 6.0);
    cfg.body_force = Some([1e-5, 0.0, 0.0]);
    let mut sim = LbmSimulation::new(cfg);
    let mass_before = total_mass(&sim);
    sim.step_n(500);
    let mass_after = total_mass(&sim);
    assert!(
        ((mass_after - mass_before) / mass_before).abs() < 1e-10,
        "forcing must not create mass: {mass_before} -> {mass_after}"
    );
}

#[test]
fn forced_uniform_flow_accelerates_with_the_prescribed_force() {
    // In a periodic box a uniform force F accelerates the fluid uniformly:
    // every collision adds exactly F to the momentum, independent of omega.
    // The reported (force-shifted) velocity after n steps from rest is
    // (n + 1/2) F / rho.
    let force = [1.0e-5, -2.0e-5, 3.0e-5];
    for lattice_type in ALL_LATTICES {
        for viscosity in [1.0 / 6.0, 0.05, 0.4] {
            let mut cfg = config_for(lattice_type, 6, viscosity);
            cfg.body_force = Some(force);
            let is_3d = lattice_type != LatticeType::D2Q9;
            let mut expected_force = force;
            if !is_3d {
                expected_force[2] = 0.0;
            }
            let mut sim = LbmSimulation::new(cfg);
            let mass_before = total_mass(&sim);
            sim.step_n(100);
            assert!(((total_mass(&sim) - mass_before) / mass_before).abs() < 1e-12);
            for v in sim.get_velocity_field() {
                for d in 0..3 {
                    let expected = 100.5 * expected_force[d];
                    assert!(
                        (v[d] - expected).abs() < 1e-12,
                        "{lattice_type:?} nu={viscosity}: u[{d}] = {} expected {expected}",
                        v[d]
                    );
                }
            }
        }
    }
}

#[test]
fn two_dimensional_simulations_ignore_the_z_component_of_the_force() {
    // D2Q9 has no z direction.  A z component of the force used to enter the
    // force-shifted velocity (u_z = F_z / 2 rho) and spoil the moments of the
    // equilibrium and of the Guo source, so the mass leaked at a rate of about
    // 1e-9 per cell and step.
    let run = |force: [f64; 3]| {
        let mut cfg = LbmConfig::d2q9(8, 6, 0.1);
        cfg.body_force = Some(force);
        let mut sim = LbmSimulation::new(cfg);
        let mass_before = total_mass(&sim);
        sim.step_n(100);
        (sim, mass_before)
    };
    let (planar, _) = run([1.0e-5, -2.0e-5, 0.0]);
    let (with_fz, mass_before) = run([1.0e-5, -2.0e-5, 3.0e-5]);

    assert!(
        ((total_mass(&with_fz) - mass_before) / mass_before).abs() < 1e-13,
        "F_z must not change the mass: {mass_before} -> {}",
        total_mass(&with_fz)
    );
    assert_eq!(planar.get_density_field(), with_fz.get_density_field());
    assert_eq!(planar.get_velocity_field(), with_fz.get_velocity_field());
    for v in with_fz.get_velocity_field() {
        assert!((v[0] - 100.5e-5).abs() < 1e-12, "u_x = {}", v[0]);
        assert!((v[1] + 201.0e-5).abs() < 1e-12, "u_y = {}", v[1]);
        assert_eq!(v[2], 0.0);
    }
}

#[test]
fn two_dimensional_collision_kernel_drops_the_z_force() {
    let table = table_for(LatticeType::D2Q9);
    let (rho, u) = (1.02, [0.03, -0.02, 0.0]);
    let stress = [
        [1.0e-4, 5.0e-5, 0.0],
        [5.0e-5, -1.0e-4, 0.0],
        [0.0, 0.0, 0.0],
    ];
    let before = populations_with_stress(&table, rho, u, stress);
    let collide = |force: [f64; 3]| {
        let mut pop = before.clone();
        let setup = CollisionSetup {
            omega: 1.4,
            smagorinsky_cs: None,
            force,
        };
        collide_cell(&mut pop, &table, &setup);
        pop
    };
    let planar = collide([2.0e-4, -1.0e-4, 0.0]);
    let with_fz = collide([2.0e-4, -1.0e-4, 7.0e-4]);
    assert_eq!(planar, with_fz);
    let mass: f64 = with_fz.iter().sum();
    assert!((mass - rho).abs() < 1e-15, "mass {mass} vs {rho}");

    // The same kernel keeps the z force on a three-dimensional lattice.
    let table_3d = table_for(LatticeType::D3Q19);
    let mut pop = table_3d.weights.clone();
    let setup = CollisionSetup {
        omega: 1.4,
        smagorinsky_cs: None,
        force: [0.0, 0.0, 7.0e-4],
    };
    collide_cell(&mut pop, &table_3d, &setup);
    let momentum_z: f64 = pop
        .iter()
        .zip(&table_3d.velocities)
        .map(|(f, c)| f * c[2])
        .sum();
    assert!((momentum_z - 7.0e-4).abs() < 1e-15, "p_z = {momentum_z}");
}

#[test]
fn z_force_of_a_two_dimensional_run_survives_a_time_step_change() {
    // The rescaling to a new dt re-expresses the previous half-force velocity
    // shift; it must drop F_z as the collision does.
    let mut cfg = LbmConfig::d2q9(8, 8, 0.1);
    cfg.body_force = Some([1.0e-5, 0.0, 4.0e-5]);
    let mut sim = LbmSimulation::new(cfg);
    sim.step_n(20);
    let mass_before = total_mass(&sim);
    for _ in 0..40 {
        sim.step(0.5);
    }
    assert!(
        ((total_mass(&sim) - mass_before) / mass_before).abs() < 1e-13,
        "mass {mass_before} -> {}",
        total_mass(&sim)
    );
    let expected_physical = 1.0e-5 * (sim.elapsed_time + 0.5);
    for v in sim.get_velocity_field() {
        assert!(
            (v[0] / 0.5 - expected_physical).abs() < 1e-10,
            "u_phys = {} expected {expected_physical}",
            v[0] / 0.5
        );
        assert_eq!(v[2], 0.0);
    }
}

#[test]
fn forced_flow_without_force_component_stays_zero() {
    let mut cfg = LbmConfig::d3q19(5, 5, 5, 0.1);
    cfg.body_force = Some([2e-5, 0.0, 0.0]);
    let mut sim = LbmSimulation::new(cfg);
    sim.step_n(50);
    for v in sim.get_velocity_field() {
        assert!(v[0] > 0.0);
        assert!(v[1].abs() < 1e-14 && v[2].abs() < 1e-14, "{v:?}");
    }
}

// ---------------------------------------------------------------------------
// Smagorinsky model
// ---------------------------------------------------------------------------

/// Non-equilibrium populations that carry the second moment `pi` (a symmetric
/// 3x3 tensor) and nothing else: `f_neq_i = w_i/(2 cs^4) Pi_ab (c_ia c_ib - cs^2 delta_ab)`.
fn populations_with_stress(
    table: &LatticeTable,
    rho: f64,
    u: [f64; 3],
    pi: [[f64; 3]; 3],
) -> Vec<f64> {
    table
        .weights
        .iter()
        .zip(&table.velocities)
        .map(|(&w, &c)| {
            let cu = c[0] * u[0] + c[1] * u[1] + c[2] * u[2];
            let u_sq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
            let feq = w * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
            let mut contraction = 0.0;
            for a in 0..3 {
                for b in 0..3 {
                    let delta = if a == b { CS2 } else { 0.0 };
                    contraction += pi[a][b] * (c[a] * c[b] - delta);
                }
            }
            feq + w / (2.0 * CS2 * CS2) * contraction
        })
        .collect()
}

#[test]
fn smagorinsky_tau_solves_the_eddy_viscosity_closure() {
    // tau_eff = tau_0 + 3 Cs^2 |S| with S = -Pi / (2 rho cs^2 tau_eff),
    // |S| = sqrt(2 S:S).
    let (tau0, cs_smag, rho) = (0.62, 0.17, 1.03);
    let pi = [
        [1.2e-4, -3.0e-4, 0.5e-4],
        [-3.0e-4, -0.7e-4, 0.9e-4],
        [0.5e-4, 0.9e-4, -0.2e-4],
    ];
    let mut pi_sq: f64 = 0.0;
    for row in &pi {
        for value in row {
            pi_sq += value * value;
        }
    }
    let pi_norm = pi_sq.sqrt();
    let tau_eff = smagorinsky_tau(tau0, cs_smag, pi_norm, rho);
    assert!(tau_eff > tau0);
    let mut s_sq: f64 = 0.0;
    for row in &pi {
        for value in row {
            let s = -value / (2.0 * rho * CS2 * tau_eff);
            s_sq += s * s;
        }
    }
    let strain_magnitude = (2.0 * s_sq).sqrt();
    let eddy_viscosity = CS2 * (tau_eff - tau0);
    assert!(
        (eddy_viscosity - cs_smag * cs_smag * strain_magnitude).abs() < 1e-15,
        "nu_t = {eddy_viscosity:e} vs Cs^2 |S| = {:e}",
        cs_smag * cs_smag * strain_magnitude
    );
}

#[test]
fn smagorinsky_tau_limits() {
    assert!((smagorinsky_tau(0.7, 0.17, 0.0, 1.0) - 0.7).abs() < 1e-15);
    assert!((smagorinsky_tau(0.7, 0.0, 1e-3, 1.0) - 0.7).abs() < 1e-15);
    let weak = smagorinsky_tau(0.7, 0.17, 1e-4, 1.0);
    let strong = smagorinsky_tau(0.7, 0.17, 1e-2, 1.0);
    assert!(0.7 < weak && weak < strong);
}

#[test]
fn collision_kernel_applies_the_local_smagorinsky_frequency_on_every_lattice() {
    let (rho, cs_smag, omega0) = (1.01, 0.2, 1.0 / 0.6);
    for lattice_type in ALL_LATTICES {
        let table = table_for(lattice_type);
        let three_d = lattice_type != LatticeType::D2Q9;
        let pi = if three_d {
            [
                [2.0e-4, 1.5e-4, -0.8e-4],
                [1.5e-4, -1.0e-4, 0.6e-4],
                [-0.8e-4, 0.6e-4, -0.5e-4],
            ]
        } else {
            [
                [2.0e-4, 1.5e-4, 0.0],
                [1.5e-4, -2.0e-4, 0.0],
                [0.0, 0.0, 0.0],
            ]
        };
        let u = if three_d {
            [0.02, -0.01, 0.015]
        } else {
            [0.02, -0.01, 0.0]
        };
        let before = populations_with_stress(&table, rho, u, pi);
        let mut pop = before.clone();
        let setup = CollisionSetup {
            omega: omega0,
            smagorinsky_cs: Some(cs_smag),
            force: [0.0; 3],
        };
        collide_cell(&mut pop, &table, &setup);

        let mut pi_sq: f64 = 0.0;
        for row in &pi {
            for value in row {
                pi_sq += value * value;
            }
        }
        let omega_local = 1.0 / smagorinsky_tau(1.0 / omega0, cs_smag, pi_sq.sqrt(), rho);
        assert!(omega_local < omega0, "the model must add viscosity");
        for (i, (&w, &c)) in table.weights.iter().zip(&table.velocities).enumerate() {
            let cu = c[0] * u[0] + c[1] * u[1] + c[2] * u[2];
            let u_sq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
            let feq = w * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
            let expected = before[i] + omega_local * (feq - before[i]);
            assert!(
                (pop[i] - expected).abs() < 1e-14,
                "{lattice_type:?} direction {i}: {} vs {expected}",
                pop[i]
            );
        }
        let mass: f64 = pop.iter().sum();
        assert!((mass - rho).abs() < 1e-14);
    }
}

#[test]
fn smagorinsky_with_zero_constant_is_plain_bgk_in_3d() {
    let build = |cs: Option<f64>| {
        let mut cfg = LbmConfig::d3q19(6, 6, 6, 0.04);
        cfg.body_force = Some([1e-5, 0.0, 2e-5]);
        cfg.smagorinsky_cs = cs;
        let mut sim = LbmSimulation::new(cfg);
        if let SimGrid::Grid3D(g) = &mut sim.grid {
            for z in 0..6 {
                for y in 0..6 {
                    for x in 0..6 {
                        let uy = 0.03 * (2.0 * PI * x as f64 / 6.0).sin();
                        g.set_equilibrium(x, y, z, 1.0, 0.0, uy, 0.0);
                    }
                }
            }
        }
        sim.step_n(40);
        sim.get_velocity_field()
    };
    let plain = build(None);
    let zero = build(Some(0.0));
    for (a, b) in plain.iter().zip(&zero) {
        for d in 0..3 {
            assert!((a[d] - b[d]).abs() < 1e-14);
        }
    }
}

/// Amplitude of the `sin(2 pi x / nx)` component of `u_y`.
fn shear_amplitude(sim: &LbmSimulation) -> f64 {
    let nx = sim.config.nx;
    let vel = sim.get_velocity_field();
    let mut acc = 0.0;
    for (k, v) in vel.iter().enumerate() {
        let x = k % nx;
        acc += v[1] * (2.0 * PI * x as f64 / nx as f64).sin();
    }
    2.0 * acc / vel.len() as f64
}

fn init_shear_wave(sim: &mut LbmSimulation, amplitude: f64) {
    let (nx, ny, nz) = (sim.config.nx, sim.config.ny, sim.config.nz);
    match &mut sim.grid {
        SimGrid::Grid2D(g) => {
            for y in 0..ny {
                for x in 0..nx {
                    let uy = amplitude * (2.0 * PI * x as f64 / nx as f64).sin();
                    g.set_equilibrium(x, y, 1.0, 0.0, uy);
                }
            }
        }
        SimGrid::Grid3D(g) => {
            for z in 0..nz {
                for y in 0..ny {
                    for x in 0..nx {
                        let uy = amplitude * (2.0 * PI * x as f64 / nx as f64).sin();
                        g.set_equilibrium(x, y, z, 1.0, 0.0, uy, 0.0);
                    }
                }
            }
        }
    }
}

#[test]
fn smagorinsky_adds_dissipation_to_a_shear_wave_in_2d_and_3d() {
    for lattice_type in [LatticeType::D2Q9, LatticeType::D3Q19] {
        let run = |cs: Option<f64>| {
            let mut cfg = match lattice_type {
                LatticeType::D2Q9 => LbmConfig::d2q9(32, 4, 0.01),
                _ => LbmConfig::d3q19(32, 4, 4, 0.01),
            };
            cfg.smagorinsky_cs = cs;
            let mut sim = LbmSimulation::new(cfg);
            init_shear_wave(&mut sim, 0.08);
            sim.step_n(300);
            shear_amplitude(&sim)
        };
        let plain = run(None);
        let les = run(Some(0.3));
        assert!(plain > 0.0 && les > 0.0);
        assert!(
            les < plain * 0.995,
            "{lattice_type:?}: LES amplitude {les} must be clearly below plain BGK {plain}"
        );
        // The wave must still be a decaying wave, not noise.
        assert!(les > 0.2 * plain, "{lattice_type:?}: {les} vs {plain}");
    }
}

// ---------------------------------------------------------------------------
// 3-D walls
// ---------------------------------------------------------------------------

/// Viscosity for which BGK + half-way bounce-back reproduces the Poiseuille
/// parabola exactly: tau = 1/2 + sqrt(3/16), nu = cs^2 (tau - 1/2).
fn magic_viscosity() -> f64 {
    (3.0_f64 / 16.0).sqrt() / 3.0
}

fn poiseuille_profile_error(sim: &LbmSimulation, force: f64, viscosity: f64) -> f64 {
    let (nx, ny, nz) = (sim.config.nx, sim.config.ny, sim.config.nz);
    let vel = sim.get_velocity_field();
    let mut err_sq = 0.0;
    let mut ref_sq = 0.0;
    for z in 0..nz {
        for y in 1..ny - 1 {
            for x in 0..nx {
                let yf = y as f64;
                let analytic = force / (2.0 * viscosity) * (yf - 0.5) * ((ny as f64 - 1.5) - yf);
                let numeric = vel[z * ny * nx + y * nx + x][0];
                err_sq += (numeric - analytic).powi(2);
                ref_sq += analytic * analytic;
            }
        }
    }
    (err_sq / ref_sq).sqrt()
}

#[test]
fn body_force_poiseuille_flow_between_3d_walls_matches_the_parabola() {
    let force = 1.0e-5;
    for (lattice_type, viscosity, tolerance) in [
        (LatticeType::D3Q19, magic_viscosity(), 1e-4),
        (LatticeType::D3Q27, magic_viscosity(), 1e-4),
        // Plain tau = 1: the finding's 1 % criterion.
        (LatticeType::D3Q19, 1.0 / 6.0, 1e-2),
    ] {
        let (nx, ny, nz) = (4, 18, 4);
        let mut cfg = match lattice_type {
            LatticeType::D3Q27 => LbmConfig::d3q27(nx, ny, nz, viscosity),
            _ => LbmConfig::d3q19(nx, ny, nz, viscosity),
        };
        cfg.body_force = Some([force, 0.0, 0.0]);
        let mut sim = LbmSimulation::new(cfg);
        sim.set_solid_cells_3d(&channel_walls_3d(nx, ny, nz))
            .expect("channel walls fit the grid");
        sim.step_n(4000);
        let err = poiseuille_profile_error(&sim, force, viscosity);
        assert!(
            err < tolerance,
            "{lattice_type:?} nu={viscosity}: relative L2 error {err:e} exceeds {tolerance:e}"
        );
        // The flow is purely streamwise and the solid plates are at rest.
        for (k, v) in sim.get_velocity_field().iter().enumerate() {
            let y = (k / nx) % ny;
            assert!(v[1].abs() < 1e-10 && v[2].abs() < 1e-10, "{v:?}");
            if y == 0 || y == ny - 1 {
                assert!(v[0].abs() < 1e-15, "solid cell velocity {}", v[0]);
            }
        }
    }
}

#[test]
fn walls_conserve_the_fluid_mass_and_stop_the_flow() {
    // A stream normal to the plates must be stopped by them without any loss
    // of fluid.
    let (nx, ny, nz) = (4, 10, 4);
    let mut sim = LbmSimulation::new(LbmConfig::d3q19(nx, ny, nz, 1.0 / 6.0));
    sim.set_solid_cells_3d(&channel_walls_3d(nx, ny, nz))
        .expect("channel walls fit the grid");
    if let SimGrid::Grid3D(g) = &mut sim.grid {
        for z in 0..nz {
            for y in 1..ny - 1 {
                for x in 0..nx {
                    g.set_equilibrium(x, y, z, 1.0, 0.0, 0.05, 0.0);
                }
            }
        }
    }
    let fluid_mass = |sim: &LbmSimulation| -> f64 {
        sim.get_density_field()
            .iter()
            .enumerate()
            .filter(|(k, _)| {
                let y = (k / nx) % ny;
                y != 0 && y != ny - 1
            })
            .map(|(_, rho)| rho)
            .sum()
    };
    let mass_before = fluid_mass(&sim);
    sim.step_n(600);
    let mass_after = fluid_mass(&sim);
    assert!(
        ((mass_after - mass_before) / mass_before).abs() < 1e-12,
        "fluid mass must be conserved by bounce-back: {mass_before} -> {mass_after}"
    );
    let mean_uy: f64 = sim
        .get_velocity_field()
        .iter()
        .map(|v| v[1].abs())
        .sum::<f64>()
        / (nx * ny * nz) as f64;
    assert!(
        mean_uy < 1e-3,
        "the plates must stop the stream: <|uy|> = {mean_uy}"
    );
}

#[test]
fn solid_cell_api_validates_its_input() {
    let mut flat = LbmSimulation::new(LbmConfig::d2q9(8, 8, 0.1));
    let err = flat.set_solid_cells_3d(&[[1, 1, 0]]).expect_err("2D grid");
    assert!(
        matches!(err, Error::IncompatibleParameters { .. }),
        "{err:?}"
    );

    let mut sim = LbmSimulation::new(LbmConfig::d3q19(4, 5, 6, 0.1));
    sim.set_solid_cells_3d(&[[0, 0, 0], [3, 4, 5]])
        .expect("valid cells");
    assert_eq!(sim.solid_3d.iter().filter(|&&s| s).count(), 2);
    let err = sim
        .set_solid_cells_3d(&[[0, 0, 0], [4, 0, 0]])
        .expect_err("x out of range");
    assert_eq!(
        err,
        Error::IndexOutOfBounds3D {
            x: 4,
            y: 0,
            z: 0,
            nx: 4,
            ny: 5,
            nz: 6
        }
    );
    // The previous mask is kept after a rejected update.
    assert_eq!(sim.solid_3d.iter().filter(|&&s| s).count(), 2);
    // Replacing the mask replaces the walls.
    sim.set_solid_cells_3d(&[]).expect("empty mask");
    assert_eq!(sim.solid_3d.iter().filter(|&&s| s).count(), 0);
    assert!(sim.wall_links_3d.is_empty());
}

#[test]
fn channel_wall_cells_and_links_have_the_expected_counts() {
    let (nx, ny, nz) = (4, 6, 3);
    let cells = channel_walls_3d(nx, ny, nz);
    assert_eq!(cells.len(), 2 * nx * nz);
    assert!(cells.iter().all(|c| c[1] == 0 || c[1] == ny - 1));

    let mut mask = vec![false; nx * ny * nz];
    for &[x, y, z] in &cells {
        mask[z * ny * nx + y * nx + x] = true;
    }
    // Links per fluid cell next to a plate = directions pointing into it.
    for (lattice_type, into_plate) in [(LatticeType::D3Q19, 5), (LatticeType::D3Q27, 9)] {
        let table = table_for(lattice_type);
        let links = build_wall_links_3d(&mask, [nx, ny, nz], &table);
        assert_eq!(links.len(), 2 * nx * nz * into_plate, "{lattice_type:?}");
        for link in &links {
            assert!(!mask[link.fluid]);
            assert!(mask[link.solid]);
            assert!(table.velocities[link.dir][1] != 0.0);
        }
    }
    assert!(channel_walls_3d(4, 0, 4).is_empty());
    assert_eq!(channel_walls_3d(2, 1, 2).len(), 4);
}

// ---------------------------------------------------------------------------
// The dt argument of step()
// ---------------------------------------------------------------------------

#[test]
fn invalid_time_steps_are_rejected_without_touching_the_state() {
    let mut sim = LbmSimulation::new(LbmConfig::d2q9(6, 6, 0.1));
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let err = sim.try_step(dt).expect_err("invalid dt");
        assert!(matches!(err, Error::ParameterOutOfBounds { .. }), "{err:?}");
        sim.step(dt);
    }
    assert_eq!(sim.step_count, 0);
    assert_eq!(sim.elapsed_time, 0.0);
    sim.try_step(0.5).expect("valid dt");
    assert_eq!(sim.step_count, 1);
    assert!((sim.elapsed_time - 0.5).abs() < 1e-15);
}

#[test]
fn elapsed_time_accumulates_the_step_durations() {
    let mut sim = LbmSimulation::new(LbmConfig::d2q9(6, 6, 0.1));
    sim.step(1.0);
    sim.step(0.25);
    sim.step(0.25);
    sim.step_n(2);
    assert_eq!(sim.step_count, 5);
    assert!(
        (sim.elapsed_time - 3.5).abs() < 1e-15,
        "{}",
        sim.elapsed_time
    );
}

#[test]
fn dt_scales_the_lattice_force_with_its_square() {
    // Physical acceleration F: u(t) = F (t + dt/2) with the half-step offset
    // of the Guo velocity; the lattice velocity is dt times the physical one.
    let force = 1.0e-5;
    for dt in [1.0, 0.5, 0.25] {
        let mut cfg = LbmConfig::d2q9(8, 8, 0.1);
        cfg.body_force = Some([force, 0.0, 0.0]);
        let mut sim = LbmSimulation::new(cfg);
        let steps = (100.0 / dt) as usize;
        for _ in 0..steps {
            sim.step(dt);
        }
        assert!((sim.elapsed_time - 100.0).abs() < 1e-9);
        let expected_physical = force * (sim.elapsed_time + 0.5 * dt);
        for v in sim.get_velocity_field() {
            assert!(
                (v[0] / dt - expected_physical).abs() < 1e-12,
                "dt={dt}: u_phys = {} expected {expected_physical}",
                v[0] / dt
            );
        }
    }
}

#[test]
fn changing_dt_keeps_the_physical_velocity_continuous() {
    let force = 1.0e-5;
    let mut cfg = LbmConfig::d3q19(5, 5, 5, 0.08);
    cfg.body_force = Some([force, 0.0, 0.0]);
    let mut sim = LbmSimulation::new(cfg);
    for _ in 0..100 {
        sim.step(1.0);
    }
    let mass_before = total_mass(&sim);
    for _ in 0..200 {
        sim.step(0.5);
    }
    assert!(((total_mass(&sim) - mass_before) / mass_before).abs() < 1e-12);
    assert!((sim.elapsed_time - 200.0).abs() < 1e-9);
    // 100 steps of dt=1 give u = 100.5 F; each later half step adds 0.5 F.
    let expected_physical = force * (sim.elapsed_time + 0.5);
    for v in sim.get_velocity_field() {
        assert!(
            (v[0] / 0.5 - expected_physical).abs() < 1e-10,
            "u_phys = {} expected {expected_physical}",
            v[0] / 0.5
        );
    }
}

#[test]
fn shear_wave_decay_is_independent_of_the_time_step() {
    // The physical decay exp(-nu k^2 t) must not depend on how the time is
    // sliced.  Run A: 200 steps of dt = 1.  Run B: 100 steps of dt = 1 followed
    // by 200 steps of dt = 1/2 (also 200 time units).  The old code ignored dt
    // and therefore ran 300 steps of the same viscosity in run B.
    let nu = 1.0 / 6.0;
    let amplitude = 0.01;
    let make = || {
        let mut sim = LbmSimulation::new(LbmConfig::d2q9(64, 4, nu));
        init_shear_wave(&mut sim, amplitude);
        sim
    };
    let mut reference = make();
    for _ in 0..200 {
        reference.step(1.0);
    }
    let reference_amplitude = shear_amplitude(&reference);

    let mut sliced = make();
    for _ in 0..100 {
        sliced.step(1.0);
    }
    for _ in 0..200 {
        sliced.step(0.5);
    }
    let sliced_amplitude = shear_amplitude(&sliced) / 0.5;

    // Analytic decay for reference.
    let k = 2.0 * PI / 64.0;
    let analytic = amplitude * (-nu * k * k * 200.0).exp();
    assert!(
        (reference_amplitude - analytic).abs() / analytic < 0.01,
        "reference run {reference_amplitude} vs analytic {analytic}"
    );
    assert!(
        (sliced_amplitude - reference_amplitude).abs() / reference_amplitude < 0.01,
        "sliced run {sliced_amplitude} vs reference {reference_amplitude}"
    );
    assert!((sliced.elapsed_time - 200.0).abs() < 1e-9);

    // The lattice velocity of the sliced run is dt times smaller.
    let lattice_amplitude = shear_amplitude(&sliced);
    assert!(
        (lattice_amplitude / reference_amplitude - 0.5).abs() < 0.01,
        "lattice amplitude ratio {}",
        lattice_amplitude / reference_amplitude
    );
}

#[test]
fn rescaling_to_a_new_time_step_conserves_mass_and_scales_velocity() {
    let mut sim = LbmSimulation::new(LbmConfig::d2q9(32, 8, 0.1));
    init_shear_wave(&mut sim, 0.05);
    sim.step_n(10);
    let mass_before = total_mass(&sim);
    let amplitude_before = shear_amplitude(&sim);
    // A very short step: almost no evolution, but the lattice amplitude must
    // already be halved by the change of units.
    sim.step(0.5);
    assert!(((total_mass(&sim) - mass_before) / mass_before).abs() < 1e-13);
    let ratio = shear_amplitude(&sim) / amplitude_before;
    assert!((ratio - 0.5).abs() < 5e-3, "ratio = {ratio}");
}

#[test]
fn adaptive_dt_reduces_the_lattice_mach_number() {
    // A fast uniform stream (Ma > 0.5) must make the controller shrink dt and
    // the simulation must honour it: the lattice speed follows dt while the
    // physical speed is unchanged.
    let u0 = 0.3;
    let mut sim = LbmSimulation::new(LbmConfig::d2q9(16, 16, 0.1));
    if let SimGrid::Grid2D(g) = &mut sim.grid {
        for y in 0..16 {
            for x in 0..16 {
                g.set_equilibrium(x, y, 1.0, u0, 0.0);
            }
        }
    }
    let mut run = FullStepLoop::new(30, 1e-14, 5, 0);
    run.enable_adaptive_dt(0.1);
    run.run(&mut sim);

    let controller = run.adaptive_dt.as_ref().expect("controller enabled");
    assert!(controller.n_adjustments() >= 1);
    assert!(
        controller.dt_current < 0.25,
        "dt = {}",
        controller.dt_current
    );
    assert!(
        sim.elapsed_time < 0.5 * sim.step_count as f64,
        "most steps must be short: {} time units in {} steps",
        sim.elapsed_time,
        sim.step_count
    );
    // The uniform stream is a steady state, so the physical speed is conserved
    // exactly while the lattice speed shrinks with the time step.
    let last_dt = sim.last_dt;
    assert!(last_dt < 0.25);
    for v in sim.get_velocity_field() {
        assert!(
            (v[0] / last_dt - u0).abs() < 1e-9,
            "physical speed {} instead of {u0}",
            v[0] / last_dt
        );
        assert!(v[0] < 0.06, "lattice speed {} is still compressible", v[0]);
    }
}

// ---------------------------------------------------------------------------
// Ensemble perturbation
// ---------------------------------------------------------------------------

#[test]
fn ensemble_members_are_identical_without_perturbation() {
    let mut cfg = LbmConfig::d2q9(8, 8, 1.0 / 6.0);
    cfg.body_force = Some([1e-5, 0.0, 0.0]);
    let mut runner = EnsembleRunner::new(cfg, 3, 50);
    runner.perturbation_amplitude = 0.0;
    runner.run_all();
    assert!(runner.velocity_variance() < 1e-24);
    assert_eq!(
        runner.results[0].mean_velocity,
        runner.results[2].mean_velocity
    );
}

#[test]
fn ensemble_perturbation_survives_the_first_collision_in_2d_and_3d() {
    for cfg in [
        LbmConfig::d2q9(8, 8, 1.0 / 6.0),
        LbmConfig::d3q19(4, 4, 4, 1.0 / 6.0),
    ] {
        let mut cfg = cfg;
        cfg.body_force = Some([1e-5, 0.0, 0.0]);
        let mut runner = EnsembleRunner::new(cfg, 3, 50);
        runner.perturbation_amplitude = 1e-3;
        runner.run_all();
        let variance = runner.velocity_variance();
        assert!(
            variance > 1e-10,
            "members must differ, variance = {variance:e}"
        );
        assert!(runner.results[0].mean_velocity != runner.results[2].mean_velocity);
    }
}

#[test]
fn ensemble_perturbation_is_part_of_the_initial_state() {
    // After zero steps the perturbed member must carry its perturbation.
    let cfg = LbmConfig::d2q9(8, 8, 1.0 / 6.0);
    let mut runner = EnsembleRunner::new(cfg, 2, 0);
    runner.perturbation_amplitude = 2e-3;
    runner.run_all();
    // member 0: amplitude 2e-3, member 1: 4e-3 -> the second is perturbed
    // twice as much.
    let m0 = runner.results[0].mean_velocity;
    let m1 = runner.results[1].mean_velocity;
    assert!(m0 > 0.0 && m1 > 1.9 * m0 && m1 < 2.1 * m0, "{m0} vs {m1}");
}

#[test]
fn table_is_sized_for_the_largest_lattice() {
    for lattice_type in ALL_LATTICES {
        let table = table_for(lattice_type);
        assert!(table.q <= MAX_Q);
        let weight_sum: f64 = table.weights.iter().sum();
        assert!((weight_sum - 1.0).abs() < 1e-14);
    }
    assert_eq!(table_for(LatticeType::D3Q27).q, MAX_Q);
}
