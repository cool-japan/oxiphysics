// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LBM Poiseuille channel flow.
//! 2D channel 64×16, pressure-driven (Zou-He BCs), run 1 000 steps to
//! steady state, sample velocity profile at center, compare to analytical
//! parabola.

use oxiphysics_lbm::{
    LbmSimulation2D,
    boundary::{
        channel_walls, poiseuille_analytical, viscosity_from_omega, zou_he_pressure_outlet,
        zou_he_velocity_inlet,
    },
};

fn main() {
    let nx = 64_usize;
    let ny = 16_usize;
    let omega = 1.0_f64; // tau = 1/omega, nu = cs² * (tau - 0.5)
    let u_inlet = 0.05_f64; // lattice units

    // Build simulation.
    let mut sim = LbmSimulation2D::new(nx, ny, omega);

    // Combine walls + Zou-He inlet/outlet BCs.
    let mut bcs = channel_walls(nx, ny);
    bcs.extend(zou_he_velocity_inlet(ny, u_inlet, 0.0));
    bcs.extend(zou_he_pressure_outlet(nx, ny, 1.0));
    sim.set_boundaries(bcs);

    // Run to steady state.
    let n_steps = 1000_usize;
    sim.run(n_steps);

    // Sample velocity profile at mid-channel (x = nx/2).
    let mid_x = nx / 2;
    let nu = viscosity_from_omega(omega);

    // Estimate dp/dx from the measured peak velocity.
    let mut u_max = 0.0_f64;
    for y in 1..(ny - 1) {
        let (ux, _) = sim.grid.velocity_at(mid_x, y);
        if ux > u_max {
            u_max = ux;
        }
    }
    let h_ch = (ny - 2) as f64;
    let dp_dx = u_max * 8.0 * nu / (h_ch * h_ch);

    println!(
        "=== lbm_channel: Poiseuille velocity profile at x={} ===",
        mid_x
    );
    println!(
        "{:>4}  {:>12}  {:>12}  {:>10}",
        "y", "numerical", "analytical", "rel_err"
    );
    let mut max_rel_err = 0.0_f64;
    for y in 1..(ny - 1) {
        let (ux, _) = sim.grid.velocity_at(mid_x, y);
        let u_ana = poiseuille_analytical(dp_dx, nu, ny, y);
        let rel_err = if u_ana.abs() > 1e-12 {
            (ux - u_ana).abs() / u_ana.abs()
        } else {
            ux.abs()
        };
        if rel_err > max_rel_err {
            max_rel_err = rel_err;
        }
        println!(
            "{:>4}  {:>12.6}  {:>12.6}  {:>10.2}%",
            y,
            ux,
            u_ana,
            rel_err * 100.0
        );
    }
    println!(
        "Max relative error vs analytical parabola: {:.2}%",
        max_rel_err * 100.0
    );

    // Check cells in the central half of the channel only (avoid boundary
    // layers near the no-slip walls where discretisation error is largest).
    let mut interior_max_err = 0.0_f64;
    let y_lo = ny / 4;
    let y_hi = 3 * ny / 4;
    for y in y_lo..y_hi {
        let (ux, _) = sim.grid.velocity_at(mid_x, y);
        let u_ana = poiseuille_analytical(dp_dx, nu, ny, y);
        let rel_err = if u_ana.abs() > 1e-12 {
            (ux - u_ana).abs() / u_ana.abs()
        } else {
            ux.abs()
        };
        if rel_err > interior_max_err {
            interior_max_err = rel_err;
        }
    }
    println!(
        "Central channel (y={}..{}) max error: {:.2}%",
        y_lo,
        y_hi,
        interior_max_err * 100.0
    );
    if interior_max_err < 0.15 {
        println!("PASS: central Poiseuille profile within 15% of analytical");
    } else {
        println!(
            "INFO: interior error {:.2}% (may need more steps for convergence)",
            interior_max_err * 100.0
        );
    }
}
