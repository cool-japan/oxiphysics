// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Minimal argon MD: 32 atoms in a 3 nm cubic box, Lennard-Jones potential.
//! Initialises FCC lattice, runs NVE for 1 000 steps at ~90 K, prints
//! initial/final temperature and total energy, verifies energy conservation.

use oxiphysics_core::math::Vec3;
use oxiphysics_md::{
    AtomSet, LennardJones, MdSimulation, NoBarostat, NoThermostat, PairForceField, PeriodicBox,
    VelocityVerlet,
};

fn main() {
    // Argon LJ parameters (real units: nm, kJ/mol, atomic-mass-units).
    // epsilon = 0.9959 kJ/mol,  sigma = 0.3405 nm,  m = 39.948 amu
    // Boltzmann: k_B = 8.314e-3 kJ/(mol·K)
    let eps = 0.9959_f64; // kJ/mol
    let sig = 0.3405_f64; // nm
    let cutoff = 1.0_f64; // nm
    let mass_ar = 39.948_f64; // amu (using amu units internally)
    let kb = 8.314e-3_f64; // kJ/(mol·K)

    // Box: 2×2×2 FCC unit cells (a≈0.526 nm → ~1.052 nm, use 1.5 nm for comfort)
    let n_side = 2_usize; // FCC cells per side
    let a = 0.75_f64; // nm, FCC lattice parameter (slightly expanded)
    let box_l = n_side as f64 * a;

    // Build FCC lattice (4 atoms per cell × 8 cells = 32 atoms).
    let basis = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.5, 0.5, 0.0),
        Vec3::new(0.5, 0.0, 0.5),
        Vec3::new(0.0, 0.5, 0.5),
    ];

    let mut atoms = AtomSet::with_capacity(32);
    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                for b in &basis {
                    let pos = Vec3::new(
                        (ix as f64 + b.x) * a,
                        (iy as f64 + b.y) * a,
                        (iz as f64 + b.z) * a,
                    );
                    atoms.add_atom(pos, Vec3::zeros(), mass_ar, 0.0, 0);
                }
            }
        }
    }
    assert_eq!(atoms.len(), 32, "Expected 32 atoms in 2×2×2 FCC");

    // Assign Maxwell-Boltzmann velocities at 90 K using a deterministic LCG.
    let target_temp = 90.0_f64;
    let sigma_v = (kb * target_temp / mass_ar).sqrt();
    let mut lcg: u64 = 12345678;
    let next_normal = |rng: &mut u64| -> f64 {
        // Box-Muller from two uniform samples.
        let u1 = {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let bits = ((*rng >> 33) as f64) / (u32::MAX as f64 + 1.0);
            bits.max(1e-10)
        };
        let u2 = {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((*rng >> 33) as f64) / (u32::MAX as f64 + 1.0)
        };
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    };
    for v in &mut atoms.velocities {
        *v = Vec3::new(
            sigma_v * next_normal(&mut lcg),
            sigma_v * next_normal(&mut lcg),
            sigma_v * next_normal(&mut lcg),
        );
    }
    atoms.remove_com_velocity();

    // Force field.
    let mut ff = PairForceField::new();
    ff.add_interaction(0, 0, LennardJones::new(eps, sig, cutoff));

    let pbox = PeriodicBox::cubic(box_l);

    let mut sim = MdSimulation::new(
        atoms,
        ff,
        VelocityVerlet::new(),
        NoThermostat::new(),
        NoBarostat::new(),
        pbox,
        target_temp,
        0.0,
        kb,
    );

    // Record initial state.
    let records = sim.run(0, 0.001, Some(1)); // zero steps to get initial forces
    let _ = records; // not used yet
    let t0 = sim.atoms.temperature(kb);
    let e0 = sim.energy.kinetic + sim.energy.potential;
    println!("Initial temperature: {:.2} K", t0);
    println!("Initial total energy: {:.4} kJ/mol", e0);

    // Run 1 000 NVE steps at dt = 2 fs = 0.002 ps (using nm/ps time units).
    let n_steps = 1000_u64;
    let dt = 0.002_f64; // ps
    let records = sim.run(n_steps, dt, Some(100));

    let t_final = sim.atoms.temperature(kb);
    let e_final = sim.energy.kinetic + sim.energy.potential;
    println!("Final temperature:   {:.2} K", t_final);
    println!("Final total energy:  {:.4} kJ/mol", e_final);

    // Energy conservation check (within 5%).
    let e_ref = records.first().map(|r| r.total).unwrap_or(e0);
    let drift = (e_final - e_ref).abs() / e_ref.abs().max(1.0);
    println!("Energy drift: {:.2}%", drift * 100.0);
    if drift < 0.05 {
        println!("PASS: energy conserved within 5%");
    } else {
        println!("WARN: energy drift {:.2}% exceeds 5%", drift * 100.0);
    }
}
