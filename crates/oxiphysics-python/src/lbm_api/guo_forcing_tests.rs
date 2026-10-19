// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the Guo body-force scheme of the D2Q9 and D3Q19
//! simulations.
//!
//! The old collision shifted the equilibrium velocity by `a / (2 omega)`,
//! which injects only `rho a / 2` of momentum per step. The tests below pin the
//! exact momentum balance (`rho a` per step, for every `omega`), the velocity
//! read-back convention and the steady plane Poiseuille profile
//! `u(s) = a s (H - s) / (2 nu)` of a walled channel, whose maximum is
//! `a H^2 / (8 nu)`.

use super::forcing::{forced_velocity, guo_source, stored_velocity};
use super::{
    D3Q19_EX, D3Q19_EY, D3Q19_EZ, D3Q19_W, EX, EY, LbmBoundary, PyLbm3dConfig, PyLbm3dSimulation,
    PyLbmConfig, PyLbmSimulation, W, apply_moving_wall_3d, zou_he_velocity_inlet,
};

/// Number of fluid layers between the two walls of the Poiseuille channels
/// (odd, so that the centre line lies on a lattice node).
const CHANNEL_FLUID_LAYERS: usize = 21;
/// Centre-line speed the channel accelerations are chosen for.
const TARGET_CENTRE_SPEED: f64 = 0.02;
/// Relative tolerance of the discrete profile against the analytic one.
const PROFILE_TOLERANCE: f64 = 0.02;

/// Body acceleration that gives the centre speed `u_max` in a channel of
/// `layers` fluid layers: `a = 8 nu u_max / H^2`.
fn channel_acceleration(viscosity: f64, layers: usize) -> f64 {
    let width = layers as f64;
    8.0 * viscosity * TARGET_CENTRE_SPEED / (width * width)
}

/// Plane Poiseuille speed at distance `s` from a wall of a channel of width `width`.
fn poiseuille_speed(accel: f64, viscosity: f64, width: f64, s: f64) -> f64 {
    accel * s * (width - s) / (2.0 * viscosity)
}

// ---------------------------------------------------------------------------
// Source term
// ---------------------------------------------------------------------------

#[test]
fn test_guo_source_d2q9_conserves_mass_and_scales_the_force_by_one_minus_half_omega() {
    let u = [0.031, -0.017];
    let force = [1.3e-3, -2.1e-3];
    for omega in [0.6, 1.0, 1.25, 1.7] {
        let mut mass = 0.0;
        let mut momentum = [0.0; 2];
        for q in 0..9 {
            let e_q = [f64::from(EX[q]), f64::from(EY[q])];
            let source = guo_source(omega, W[q], e_q, u, force);
            mass += source;
            momentum[0] += e_q[0] * source;
            momentum[1] += e_q[1] * source;
        }
        let scale = 1.0 - 0.5 * omega;
        assert!(mass.abs() < 1e-18, "omega {omega}: sum S_q = {mass}");
        for k in 0..2 {
            assert!(
                (momentum[k] - scale * force[k]).abs() < 1e-17,
                "omega {omega}: momentum {k} is {}, expected {}",
                momentum[k],
                scale * force[k]
            );
        }
    }
}

#[test]
fn test_guo_source_d3q19_conserves_mass_and_scales_the_force_by_one_minus_half_omega() {
    let u = [0.031, -0.017, 0.022];
    let force = [1.3e-3, -2.1e-3, 0.7e-3];
    for omega in [0.6, 1.0, 1.25, 1.7] {
        let mut mass = 0.0;
        let mut momentum = [0.0; 3];
        for q in 0..19 {
            let e_q = [
                f64::from(D3Q19_EX[q]),
                f64::from(D3Q19_EY[q]),
                f64::from(D3Q19_EZ[q]),
            ];
            let source = guo_source(omega, D3Q19_W[q], e_q, u, force);
            mass += source;
            for k in 0..3 {
                momentum[k] += e_q[k] * source;
            }
        }
        let scale = 1.0 - 0.5 * omega;
        assert!(mass.abs() < 1e-18, "omega {omega}: sum S_q = {mass}");
        for k in 0..3 {
            assert!(
                (momentum[k] - scale * force[k]).abs() < 1e-17,
                "omega {omega}: momentum {k} is {}, expected {}",
                momentum[k],
                scale * force[k]
            );
        }
    }
}

#[test]
fn test_forced_and_stored_velocity_are_inverse() {
    let accel = [2.0e-4, -5.0e-5, 1.0e-4];
    let velocity = [0.03, 0.01, -0.02];
    let rho = 1.4;
    let stored = stored_velocity(velocity, accel);
    // The populations of a cell with density `rho` and stored velocity carry
    // the momentum rho * (u - a / 2).
    let momentum = stored.map(|v| rho * v);
    let read_back = forced_velocity(momentum, rho, accel);
    for k in 0..3 {
        assert!(
            (read_back[k] - velocity[k]).abs() < 1e-15,
            "component {k}: {} instead of {}",
            read_back[k],
            velocity[k]
        );
    }
}

#[test]
fn test_forced_velocity_of_a_vacuum_cell_is_zero() {
    assert_eq!(
        forced_velocity([1e-3, -2e-3], 0.0, [1e-4, 1e-4]),
        [0.0, 0.0]
    );
    assert_eq!(
        forced_velocity([1e-3, -2e-3], 1e-16, [1e-4, 1e-4]),
        [0.0, 0.0]
    );
}

// ---------------------------------------------------------------------------
// Momentum balance in a uniform periodic box
// ---------------------------------------------------------------------------

#[test]
fn test_d2q9_body_force_adds_exactly_its_acceleration_per_step() {
    let accel = [3.0e-5, -1.2e-5];
    // omega = 1.25, 1.0 and about 1.79: the gain per step must not depend on it.
    for viscosity in [0.1, 1.0 / 6.0, 0.02] {
        let mut cfg = PyLbmConfig::new(6, 5, viscosity);
        cfg.body_force = accel;
        cfg.init_velocity = [0.01, 0.02];
        let mut sim = PyLbmSimulation::new(&cfg);
        let mut previous = sim.velocity_at(2, 3);
        for (k, expected) in [0.01, 0.02].iter().enumerate() {
            assert!(
                (previous[k] - expected).abs() < 1e-14,
                "nu {viscosity}: initial velocity {k} is {}",
                previous[k]
            );
        }
        for step in 1..=10 {
            sim.step();
            let current = sim.velocity_at(2, 3);
            for k in 0..2 {
                assert!(
                    (current[k] - previous[k] - accel[k]).abs() < 1e-12,
                    "nu {viscosity}, step {step}: velocity {k} gained {} instead of {}",
                    current[k] - previous[k],
                    accel[k]
                );
            }
            previous = current;
        }
        // Every cell of the uniform periodic box agrees, and mass is conserved.
        let field = sim.get_velocity_field();
        for cell in 0..6 * 5 {
            assert!((field[2 * cell] - previous[0]).abs() < 1e-12);
            assert!((field[2 * cell + 1] - previous[1]).abs() < 1e-12);
        }
        for rho in sim.get_density_field() {
            assert!((rho - 1.0).abs() < 1e-12, "density drifted to {rho}");
        }
    }
}

#[test]
fn test_d3q19_body_force_adds_exactly_its_acceleration_per_step() {
    let accel = [3.0e-5, -1.2e-5, 2.0e-5];
    for viscosity in [0.1, 1.0 / 6.0, 0.02] {
        let mut cfg = PyLbm3dConfig::new(4, 3, 5, viscosity);
        cfg.body_force = accel;
        cfg.init_velocity = [0.01, 0.02, -0.015];
        let mut sim = PyLbm3dSimulation::new(&cfg);
        let mut previous = sim.velocity_at(1, 2, 3);
        for (k, expected) in [0.01, 0.02, -0.015].iter().enumerate() {
            assert!(
                (previous[k] - expected).abs() < 1e-14,
                "nu {viscosity}: initial velocity {k} is {}",
                previous[k]
            );
        }
        for step in 1..=10 {
            sim.step();
            let current = sim.velocity_at(1, 2, 3);
            for k in 0..3 {
                assert!(
                    (current[k] - previous[k] - accel[k]).abs() < 1e-12,
                    "nu {viscosity}, step {step}: velocity {k} gained {} instead of {}",
                    current[k] - previous[k],
                    accel[k]
                );
            }
            previous = current;
        }
        for rho in sim.density_field() {
            assert!((rho - 1.0).abs() < 1e-12, "density drifted to {rho}");
        }
    }
}

#[test]
fn test_d2q9_body_force_conserves_mass_around_a_wall_cell() {
    let mut cfg = PyLbmConfig::new(5, 5, 0.1);
    cfg.body_force = [1e-3, 0.0];
    let mut sim = PyLbmSimulation::new(&cfg);
    sim.set_boundary(2, 2, LbmBoundary::NoSlipWall);
    sim.run(20);
    assert_eq!(sim.velocity_at(2, 2), vec![0.0, 0.0]);
    // Mass is conserved with a wall in the box.
    let total: f64 = sim.get_density_field().iter().sum();
    assert!((total - 25.0).abs() < 1e-9, "total mass {total}");
}

// ---------------------------------------------------------------------------
// Velocity read-back conventions
// ---------------------------------------------------------------------------

#[test]
fn test_d2q9_reset_to_equilibrium_reads_the_velocity_back_under_body_force() {
    let mut cfg = PyLbmConfig::new(6, 6, 0.05);
    cfg.body_force = [4.0e-4, -2.0e-4];
    let mut sim = PyLbmSimulation::new(&cfg);
    sim.run(5);
    sim.reset_to_equilibrium(1.2, 0.03, -0.01);
    for (x, y) in [(0, 0), (3, 2), (5, 5)] {
        let v = sim.velocity_at(x, y);
        assert!((v[0] - 0.03).abs() < 1e-14, "ux = {}", v[0]);
        assert!((v[1] + 0.01).abs() < 1e-14, "uy = {}", v[1]);
        assert!((sim.density_at(x, y) - 1.2).abs() < 1e-14);
    }
}

#[test]
fn test_d3q19_moving_wall_reads_the_prescribed_velocity_back_under_body_force() {
    let mut cfg = PyLbm3dConfig::new(4, 4, 4, 0.05);
    cfg.body_force = [2.0e-4, 0.0, -1.0e-4];
    let mut sim = PyLbm3dSimulation::new(&cfg);
    apply_moving_wall_3d(&mut sim, 3, vec![0.08, 0.0, 0.0]);
    let v = sim.velocity_at(1, 2, 3);
    assert!((v[0] - 0.08).abs() < 1e-14, "ux = {}", v[0]);
    assert!(v[1].abs() < 1e-14 && v[2].abs() < 1e-14, "{v:?}");
    // Other layers keep the initial rest state.
    let rest = sim.velocity_at(1, 2, 1);
    assert!(rest.iter().all(|c| c.abs() < 1e-14), "{rest:?}");
}

#[test]
fn test_d2q9_inlet_reads_the_prescribed_velocity_back_under_body_force() {
    let mut cfg = PyLbmConfig::new(8, 6, 0.05);
    cfg.body_force = [3.0e-4, -1.0e-4];
    let mut sim = PyLbmSimulation::new(&cfg);
    zou_he_velocity_inlet(&mut sim, 0.05);
    for y in 1..5 {
        let v = sim.velocity_at(0, y);
        assert!((v[0] - 0.05).abs() < 1e-14, "inlet ux at y = {y}: {}", v[0]);
        assert!(v[1].abs() < 1e-14, "inlet uy at y = {y}: {}", v[1]);
    }
    // Cells behind the inlet are untouched.
    let inside = sim.velocity_at(3, 2);
    assert!(inside.iter().all(|c| c.abs() < 1e-14), "{inside:?}");
}

#[test]
fn test_d2q9_lid_reads_the_prescribed_velocity_back_under_body_force() {
    // A one-row grid: every cell is a lid cell and streaming keeps the row
    // uniform, so the lid velocity must be read back exactly after a step. The
    // old populations would read back `lid + a / 2`.
    let mut cfg = PyLbmConfig::new(4, 1, 0.1);
    cfg.body_force = [3.0e-3, -2.0e-3];
    let mut sim = PyLbmSimulation::new(&cfg);
    sim.set_lid_velocity(Some(vec![0.05, 0.01]));
    sim.step();
    for x in 0..4 {
        let v = sim.velocity_at(x, 0);
        assert!((v[0] - 0.05).abs() < 1e-14, "lid ux at x = {x}: {}", v[0]);
        assert!((v[1] - 0.01).abs() < 1e-14, "lid uy at x = {x}: {}", v[1]);
    }
}

#[test]
fn test_d3q19_body_force_does_not_shift_the_initial_velocity() {
    // The initial populations used to be seeded with `u0 + a`, which adds an
    // acceleration to a velocity.
    let mut cfg = PyLbm3dConfig::new(3, 3, 3, 0.1);
    cfg.body_force = [1.0e-2, 2.0e-2, -1.0e-2];
    let sim = PyLbm3dSimulation::new(&cfg);
    let v = sim.velocity_at(1, 1, 1);
    assert!(v.iter().all(|c| c.abs() < 1e-14), "{v:?}");
}

#[test]
fn test_channel_flow_preset_is_sized_for_a_low_mach_channel() {
    // With the exact forcing the steady centre speed of the walled preset is
    // a H^2 / (8 nu); it has to stay far below the lattice sound speed 0.577.
    let cfg = PyLbmConfig::channel_flow();
    let layers = (cfg.height - 2) as f64;
    let centre_speed = cfg.body_force[0] * layers * layers / (8.0 * cfg.viscosity);
    assert!(
        (0.01..0.1).contains(&centre_speed),
        "preset centre speed {centre_speed} is not a low-Mach flow"
    );
}

// ---------------------------------------------------------------------------
// Plane Poiseuille flow
// ---------------------------------------------------------------------------

/// Steps after which a channel of `CHANNEL_FLUID_LAYERS` layers is steady to
/// well below the profile tolerance: about ten times the decay time
/// `H^2 / (pi^2 nu)` of its slowest diffusive mode.
fn steady_state_steps(viscosity: f64) -> u64 {
    let width = CHANNEL_FLUID_LAYERS as f64;
    (10.0 * width * width / (std::f64::consts::PI.powi(2) * viscosity)).ceil() as u64
}

fn check_d2q9_channel(viscosity: f64) {
    let layers = CHANNEL_FLUID_LAYERS;
    let width = layers as f64;
    let accel = channel_acceleration(viscosity, layers);
    let mut cfg = PyLbmConfig::new(3, layers + 2, viscosity);
    cfg.body_force = [accel, 0.0];
    let mut sim = PyLbmSimulation::new(&cfg);
    sim.add_bottom_wall();
    sim.add_top_wall();
    assert_eq!(sim.fluid_cell_count(), 3 * layers);

    sim.run(steady_state_steps(viscosity));
    let centre = layers / 2 + 1;
    let before = sim.velocity_at(1, centre)[0];
    sim.run(100);
    let after = sim.velocity_at(1, centre)[0];
    assert!(
        (after - before).abs() < 1e-4 * TARGET_CENTRE_SPEED,
        "nu {viscosity}: the flow is not steady yet ({before} -> {after})"
    );

    // Full-way bounce-back puts the walls half a cell outside the fluid layers.
    let mut worst = 0.0f64;
    let mut flux = 0.0;
    for y in 1..=layers {
        let expected = poiseuille_speed(accel, viscosity, width, y as f64 - 0.5);
        for x in 0..3 {
            let v = sim.velocity_at(x, y);
            worst = worst.max((v[0] - expected).abs());
            assert!(
                v[1].abs() < 1e-9,
                "nu {viscosity}: transverse velocity {}",
                v[1]
            );
        }
        flux += sim.velocity_at(1, y)[0];
    }
    assert!(
        worst < PROFILE_TOLERANCE * TARGET_CENTRE_SPEED,
        "nu {viscosity}: the profile deviates by {worst} from the analytic one"
    );
    let u_max = sim.velocity_at(1, centre)[0];
    assert!(
        (u_max / TARGET_CENTRE_SPEED - 1.0).abs() < PROFILE_TOLERANCE,
        "nu {viscosity}: centre speed {u_max}, expected a H^2 / (8 nu) = {TARGET_CENTRE_SPEED}"
    );
    // Flow rate: a H^3 / (12 nu).
    let expected_flux = accel * width.powi(3) / (12.0 * viscosity);
    assert!(
        (flux / expected_flux - 1.0).abs() < PROFILE_TOLERANCE,
        "nu {viscosity}: flux {flux}, expected {expected_flux}"
    );
    // Walls are at rest, and mass is conserved.
    assert_eq!(sim.velocity_at(1, 0), vec![0.0, 0.0]);
    assert_eq!(sim.velocity_at(1, layers + 1), vec![0.0, 0.0]);
    let total: f64 = sim.get_density_field().iter().sum();
    let cells = (3 * (layers + 2)) as f64;
    assert!((total - cells).abs() < 1e-8 * cells, "total mass {total}");
}

#[test]
fn test_d2q9_channel_reproduces_the_poiseuille_profile_at_omega_one() {
    check_d2q9_channel(1.0 / 6.0);
}

#[test]
fn test_d2q9_channel_reproduces_the_poiseuille_profile_at_a_smaller_viscosity() {
    check_d2q9_channel(0.05);
}

fn check_d3q19_channel(viscosity: f64) {
    let layers = CHANNEL_FLUID_LAYERS;
    let width = layers as f64;
    let accel = channel_acceleration(viscosity, layers);
    // One cell in x and y: the flow is uniform in both, so a single column of
    // cells is the exact periodic reduction of the wide channel.
    let mut cfg = PyLbm3dConfig::new(1, 1, layers + 2, viscosity);
    cfg.body_force = [accel, 0.0, 0.0];
    let mut sim = PyLbm3dSimulation::new(&cfg);
    sim.set_boundary(0, 0, 0, LbmBoundary::NoSlipWall);
    sim.set_boundary(0, 0, layers + 1, LbmBoundary::NoSlipWall);

    sim.run(steady_state_steps(viscosity));
    let centre = layers / 2 + 1;
    let before = sim.velocity_at(0, 0, centre)[0];
    sim.run(100);
    let after = sim.velocity_at(0, 0, centre)[0];
    assert!(
        (after - before).abs() < 1e-4 * TARGET_CENTRE_SPEED,
        "nu {viscosity}: the flow is not steady yet ({before} -> {after})"
    );

    let mut worst = 0.0f64;
    for z in 1..=layers {
        let expected = poiseuille_speed(accel, viscosity, width, z as f64 - 0.5);
        let v = sim.velocity_at(0, 0, z);
        worst = worst.max((v[0] - expected).abs());
        assert!(
            v[1].abs() < 1e-9 && v[2].abs() < 1e-9,
            "nu {viscosity}: cross-flow velocity {v:?}"
        );
    }
    assert!(
        worst < PROFILE_TOLERANCE * TARGET_CENTRE_SPEED,
        "nu {viscosity}: the profile deviates by {worst} from the analytic one"
    );
    let u_max = sim.velocity_at(0, 0, centre)[0];
    assert!(
        (u_max / TARGET_CENTRE_SPEED - 1.0).abs() < PROFILE_TOLERANCE,
        "nu {viscosity}: centre speed {u_max}, expected a H^2 / (8 nu) = {TARGET_CENTRE_SPEED}"
    );
    assert_eq!(sim.velocity_at(0, 0, 0), vec![0.0; 3]);
    assert_eq!(sim.velocity_at(0, 0, layers + 1), vec![0.0; 3]);
    let total: f64 = sim.density_field().iter().sum();
    let cells = (layers + 2) as f64;
    assert!((total - cells).abs() < 1e-8 * cells, "total mass {total}");
}

#[test]
fn test_d3q19_channel_reproduces_the_poiseuille_profile_at_omega_one() {
    check_d3q19_channel(1.0 / 6.0);
}

#[test]
fn test_d3q19_channel_reproduces_the_poiseuille_profile_at_a_smaller_viscosity() {
    check_d3q19_channel(0.05);
}

#[test]
fn test_d3q19_channel_with_lateral_cells_matches_the_single_column() {
    // A 3 x 3 channel must give the same profile as the one-column reduction.
    let viscosity = 1.0 / 6.0;
    let layers = 9;
    let accel = 8.0 * viscosity * TARGET_CENTRE_SPEED / (layers * layers) as f64;
    let build = |nx: usize, ny: usize| {
        let mut cfg = PyLbm3dConfig::new(nx, ny, layers + 2, viscosity);
        cfg.body_force = [accel, 0.0, 0.0];
        let mut sim = PyLbm3dSimulation::new(&cfg);
        for y in 0..ny {
            for x in 0..nx {
                sim.set_boundary(x, y, 0, LbmBoundary::NoSlipWall);
                sim.set_boundary(x, y, layers + 1, LbmBoundary::NoSlipWall);
            }
        }
        sim.run(400);
        sim
    };
    let narrow = build(1, 1);
    let wide = build(3, 3);
    for z in 0..layers + 2 {
        let a = narrow.velocity_at(0, 0, z);
        let b = wide.velocity_at(1, 2, z);
        for k in 0..3 {
            assert!(
                (a[k] - b[k]).abs() < 1e-12,
                "layer {z}, component {k}: {} vs {}",
                a[k],
                b[k]
            );
        }
    }
}
