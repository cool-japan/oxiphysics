//! Regression tests for `MicrofluidicsLbm`.
//!
//! The solver used to collide but never stream, so no momentum ever propagated, and its
//! `ln(rho)` pseudo-potential vanished at the nominal carrier density, so the carrier never felt
//! the Shan-Chen force. These tests pin the streaming, the half-way bounce-back walls, the
//! interaction force, the conservation laws and the pressure-driven channel flow.

use super::functions::{NQ, OPP, W, feq};
use super::types::*;

fn make_sim(nx: usize, ny: usize, tau_c: f64, tau_d: f64, g_sc: f64) -> MicrofluidicsLbm {
    let hs = HeleShawParams {
        gap_height: 1e-4,
        width: 1e-3,
        viscosity: 1e-3,
    };
    let df = DeanFlowParams {
        radius: 1e-3,
        curvature_radius: 10e-3,
        mean_velocity: 0.1,
        kinematic_viscosity: 1e-6,
    };
    MicrofluidicsLbm::new(nx, ny, tau_c, tau_d, g_sc, hs, df)
}

/// Overwrite the populations of node `(ix, iy)` with equilibria of the given densities and velocity.
fn set_node(sim: &mut MicrofluidicsLbm, ix: usize, iy: usize, rho_c: f64, rho_d: f64, u: [f64; 2]) {
    let base = node_base(sim, ix, iy);
    let fc = feq(rho_c, u[0], u[1]);
    let fd = feq(rho_d, u[0], u[1]);
    sim.f[base..base + NQ].copy_from_slice(&fc);
    sim.g[base..base + NQ].copy_from_slice(&fd);
}

/// Offset of the populations of node `(ix, iy)`.
fn node_base(sim: &MicrofluidicsLbm, ix: usize, iy: usize) -> usize {
    (ix + iy * sim.nx) * NQ
}

fn fluid_mass(sim: &MicrofluidicsLbm) -> (f64, f64) {
    let mut mc = 0.0;
    let mut md = 0.0;
    for idx in 0..sim.nx * sim.ny {
        if !sim.is_wall[idx] {
            mc += sim.f[idx * NQ..(idx + 1) * NQ].iter().sum::<f64>();
            md += sim.g[idx * NQ..(idx + 1) * NQ].iter().sum::<f64>();
        }
    }
    (mc, md)
}

#[test]
fn opposite_table_reverses_every_direction() {
    use super::functions::{CX, CY};
    for q in 0..NQ {
        assert_eq!(CX[OPP[q]], -CX[q], "cx of the opposite of {q}");
        assert_eq!(CY[OPP[q]], -CY[q], "cy of the opposite of {q}");
        assert_eq!(OPP[OPP[q]], q);
    }
}

#[test]
fn uniform_fluid_stays_at_rest_and_feels_no_force() {
    for g_sc in [0.0, 0.4, -0.1] {
        let mut sim = make_sim(8, 10, 0.8, 0.9, g_sc);
        sim.run(60, 0.0);
        for iy in 1..sim.ny - 1 {
            for ix in 0..sim.nx {
                let idx = ix + iy * sim.nx;
                assert!(
                    (sim.rho_c[idx] - 1.0).abs() < 1e-13,
                    "rho_c = {}",
                    sim.rho_c[idx]
                );
                assert!(
                    (sim.rho_d[idx] - 0.1).abs() < 1e-13,
                    "rho_d = {}",
                    sim.rho_d[idx]
                );
                assert!(sim.ux[idx].abs() < 1e-14, "ux = {}", sim.ux[idx]);
                assert!(sim.uy[idx].abs() < 1e-14, "uy = {}", sim.uy[idx]);
                let (fc, fd) = sim.interaction_forces(ix, iy);
                assert!(
                    fc.iter().chain(fd.iter()).all(|v| v.abs() < 1e-15),
                    "{fc:?} {fd:?}"
                );
            }
        }
    }
}

#[test]
fn populations_are_streamed_with_periodic_wrap() {
    // tau = 1 and no forcing: the post-collision populations are the equilibria of the node.
    let mut sim = make_sim(8, 7, 1.0, 1.0, 0.0);
    let x0 = sim.nx - 1;
    let u = [0.05, 0.02];
    set_node(&mut sim, x0, 3, 1.0, 0.1, u);
    let pulse = feq(1.0, u[0], u[1]);
    let rest = feq(1.0, 0.0, 0.0);
    sim.step(0.0);

    // The +x population crossed the periodic boundary into x = 0, the +y one moved to y = 4.
    let east = node_base(&sim, 0, 3);
    assert!(
        (sim.f[east + 1] - pulse[1]).abs() < 1e-14,
        "f_1 at x = 0: {}",
        sim.f[east + 1]
    );
    let north = node_base(&sim, x0, 4);
    assert!(
        (sim.f[north + 2] - pulse[2]).abs() < 1e-14,
        "f_2 at y = 4: {}",
        sim.f[north + 2]
    );
    let north_east = node_base(&sim, 0, 4);
    assert!(
        (sim.f[north_east + 5] - pulse[5]).abs() < 1e-14,
        "f_5 at (0, 4)"
    );
    // The pulse node itself receives rest-state populations from its neighbours.
    let origin = node_base(&sim, x0, 3);
    for (q, &expected) in rest.iter().enumerate().skip(1) {
        assert!(
            (sim.f[origin + q] - expected).abs() < 1e-14,
            "refilled population {q}"
        );
    }
    // Droplet populations are streamed the same way.
    let pulse_d = feq(0.1, u[0], u[1]);
    assert!((sim.g[east + 1] - pulse_d[1]).abs() < 1e-14, "g_1 at x = 0");
}

#[test]
fn half_way_bounce_back_reflects_populations_at_the_wall() {
    let mut sim = make_sim(8, 8, 1.0, 1.0, 0.0);
    let u = [0.03, -0.04];
    set_node(&mut sim, 3, 1, 1.0, 0.1, u);
    let wall = node_base(&sim, 3, 0);
    let wall_before: Vec<f64> = sim.f[wall..wall + NQ].to_vec();
    let pulse = feq(1.0, u[0], u[1]);
    let pulse_d = feq(0.1, u[0], u[1]);
    sim.step(0.0);
    let base = node_base(&sim, 3, 1);
    // Populations that would arrive from the wall row come back reversed: q = 2, 5, 6 carry the
    // post-collision populations that left through the wall along 4, 7, 8.
    for q in [2usize, 5, 6] {
        assert!(
            (sim.f[base + q] - pulse[OPP[q]]).abs() < 1e-14,
            "carrier population {q}: {} vs {}",
            sim.f[base + q],
            pulse[OPP[q]]
        );
        assert!(
            (sim.g[base + q] - pulse_d[OPP[q]]).abs() < 1e-14,
            "droplet population {q}"
        );
    }
    // The wall node itself is not touched by collision or streaming.
    assert_eq!(&sim.f[wall..wall + NQ], wall_before.as_slice());
    // No flow has been created inside the wall row.
    assert_eq!(sim.ux[3], 0.0);
}

#[test]
fn carrier_feels_the_shan_chen_force_at_its_nominal_density() {
    // With psi = ln(rho) the carrier potential vanished at rho = 1 and the force was zero.
    let g_sc = 0.5;
    let mut sim = make_sim(10, 8, 0.8, 0.9, g_sc);
    set_node(&mut sim, 5, 4, 1.0, 0.8, [0.0, 0.0]);
    let psi = |rho: f64| 1.0 - (-rho).exp();
    let expected = -g_sc * psi(1.0) * (psi(0.8) - psi(0.1)) / 9.0;

    let (west_c, west_d) = sim.interaction_forces(4, 4);
    assert!(
        (west_c[0] - expected).abs() < 1e-15,
        "{} vs {expected}",
        west_c[0]
    );
    assert!(expected < 0.0 && west_c[1].abs() < 1e-15);
    assert!(
        west_d.iter().all(|v| v.abs() < 1e-15),
        "uniform carrier exerts no force on the droplet"
    );
    // Mirror node on the other side of the droplet-rich node: opposite force.
    let (east_c, _) = sim.interaction_forces(6, 4);
    assert!(
        (east_c[0] + expected).abs() < 1e-15,
        "{} vs {}",
        east_c[0],
        -expected
    );
    // The droplet-rich node itself has a symmetric neighbourhood.
    let (centre_c, centre_d) = sim.interaction_forces(5, 4);
    assert!(
        centre_c
            .iter()
            .chain(centre_d.iter())
            .all(|v| v.abs() < 1e-15)
    );
    // Solid nodes feel nothing.
    let (wall_c, wall_d) = sim.interaction_forces(5, 0);
    assert_eq!((wall_c, wall_d), ([0.0; 2], [0.0; 2]));
}

#[test]
fn repulsive_interaction_pushes_the_carrier_out_of_a_droplet() {
    let mut sim = make_sim(14, 10, 0.9, 0.9, 0.6);
    for iy in 4..7 {
        for ix in 5..8 {
            set_node(&mut sim, ix, iy, 1.0, 0.7, [0.0, 0.0]);
        }
    }
    sim.run(40, 0.0);
    let centre = 6 + 5 * sim.nx;
    let outside = 6 + 2 * sim.nx;
    assert!(
        sim.rho_c[centre] < sim.rho_c[outside] - 1e-4,
        "carrier must be displaced from the droplet: {} vs {}",
        sim.rho_c[centre],
        sim.rho_c[outside]
    );
}

#[test]
fn mass_of_each_phase_is_conserved_with_walls_obstacles_forcing_and_interaction() {
    let mut sim = make_sim(12, 10, 0.8, 0.9, 0.4);
    for ix in 5..8 {
        let k = ix + 5 * sim.nx;
        sim.is_wall[k] = true;
    }
    for iy in 2..4 {
        for ix in 1..4 {
            set_node(&mut sim, ix, iy, 1.0, 0.6, [0.0, 0.0]);
        }
    }
    let (mc0, md0) = fluid_mass(&sim);
    sim.run(250, 2e-5);
    let (mc1, md1) = fluid_mass(&sim);
    assert!(
        (mc1 - mc0).abs() < 1e-10 * mc0,
        "carrier mass {mc0} -> {mc1}"
    );
    assert!(
        (md1 - md0).abs() < 1e-10 * md0,
        "droplet mass {md0} -> {md1}"
    );
    // Populations never enter or leave the solid nodes.
    for k in 5 + 5 * sim.nx..8 + 5 * sim.nx {
        assert_eq!(sim.ux[k], 0.0);
        assert_eq!(sim.uy[k], 0.0);
    }
}

#[test]
fn pressure_gradient_drives_a_parabolic_profile_between_the_walls() {
    // Equal relaxation times and no interaction: the two phases form one fluid with
    // nu = (tau - 1/2) / 3 accelerated by G / (rho_c + rho_d); the no-slip walls sit half-way
    // between the first fluid row and the wall rows.
    let (nx, ny) = (4, 18);
    let tau = 0.8;
    let g = 1e-5;
    let mut sim = make_sim(nx, ny, tau, tau, 0.0);
    sim.run(3000, g);
    let nu = (tau - 0.5) / 3.0;
    let accel = g / 1.1;
    let (y_lo, y_hi) = (0.5, ny as f64 - 1.5);
    let analytic = |y: f64| accel / (2.0 * nu) * (y - y_lo) * (y_hi - y);
    let u_max = analytic(0.5 * (y_lo + y_hi));
    for iy in 1..ny - 1 {
        let idx = 2 + iy * nx;
        let expected = analytic(iy as f64);
        assert!(
            (sim.ux[idx] - expected).abs() < 0.03 * u_max,
            "row {iy}: ux = {} vs analytic {expected}",
            sim.ux[idx]
        );
        assert!(sim.uy[idx].abs() < 1e-12, "row {iy}: uy = {}", sim.uy[idx]);
        let mirror = 2 + (ny - 1 - iy) * nx;
        assert!(
            (sim.ux[idx] - sim.ux[mirror]).abs() < 1e-12,
            "profile symmetry at row {iy}"
        );
    }
}

#[test]
fn default_parameters_produce_a_centre_peaked_flow_that_vanishes_at_the_walls() {
    let (nx, ny) = (8, 18);
    let mut sim = make_sim(nx, ny, 0.8, 0.9, -0.1);
    sim.run(3000, 1e-5);
    let profile: Vec<f64> = (1..ny - 1).map(|iy| sim.ux[3 + iy * nx]).collect();
    assert!(
        profile.iter().all(|&u| u > 0.0),
        "all fluid rows move downstream: {profile:?}"
    );
    let (imax, umax) = profile
        .iter()
        .copied()
        .enumerate()
        .fold(
            (0, f64::MIN),
            |acc, (i, u)| if u > acc.1 { (i, u) } else { acc },
        );
    let mid = profile.len() / 2;
    assert!(
        imax == mid - 1 || imax == mid,
        "maximum at row {imax}, expected the centre"
    );
    assert!(
        profile[0] < 0.2 * umax && profile[profile.len() - 1] < 0.2 * umax,
        "{profile:?}"
    );
    assert!(
        sim.ux[3].abs() < 1e-15 && sim.ux[3 + (ny - 1) * nx].abs() < 1e-15,
        "wall rows are solid"
    );
}

#[test]
fn droplets_are_carried_by_the_flow_and_wrap_periodically() {
    let (nx, ny) = (10, 12);
    let mut sim = make_sim(nx, ny, 1.0, 1.0, 0.0);
    for iy in 1..ny - 1 {
        for ix in 0..nx {
            set_node(&mut sim, ix, iy, 1.0, 0.1, [0.1, 0.0]);
        }
    }
    sim.add_droplet(Droplet::new([nx as f64 - 0.05, 6.0], 1.5, 1.0, 0.001));
    sim.step(0.0);
    let d = &sim.droplets[0];
    assert!((d.vel[0] - 0.1).abs() < 1e-9, "droplet speed {}", d.vel[0]);
    assert!(
        (d.pos[0] - 0.05).abs() < 1e-9,
        "wrapped position {}",
        d.pos[0]
    );
    assert!((d.pos[1] - 6.0).abs() < 1e-9);
}

#[test]
fn empty_grid_steps_are_noops() {
    let mut sim = make_sim(0, 4, 0.8, 0.9, -0.1);
    sim.step(1e-4);
    assert!(sim.f.is_empty() && sim.rho_c.is_empty());
}

#[test]
fn equilibrium_weights_are_consistent_with_the_stored_initial_state() {
    let sim = make_sim(4, 6, 0.8, 0.9, -0.1);
    for (q, &weight) in W.iter().enumerate() {
        assert_eq!(sim.f[q], weight);
        assert!((sim.g[q] - 0.1 * weight).abs() < 1e-18);
    }
}
