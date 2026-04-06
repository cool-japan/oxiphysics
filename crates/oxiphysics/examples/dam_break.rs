// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH dam-break demo: 200 particles in a 4×4 column on the left side.
//! Gravity acts in the -Y direction; run 500 steps with dt=0.001.

use oxiphysics_core::math::Vec3;
use oxiphysics_sph::{
    kernel::CubicSplineKernel,
    neighbor::SpatialHash,
    particle::{ParticleSet, SphParticle},
    wcsph::{WcsphParams, step as wcsph_step},
};

fn main() {
    let dt = 0.001_f64;
    let n_steps = 500_usize;
    let gravity = Vec3::new(0.0, -9.81, 0.0);

    // Particle parameters (2D-ish column: 10×20 in XY, single layer in Z).
    let spacing = 0.05_f64;
    let mass = 1000.0 * spacing.powi(3); // rest_density * volume
    let h = 2.0 * spacing; // smoothing length

    let params = WcsphParams {
        rest_density: 1000.0,
        stiffness: 1000.0 * 9.81 * 4.0, // gh * rho
        gamma: 7.0,
        viscosity: 0.01,
        smoothing_length: h,
    };
    let kernel = CubicSplineKernel;

    // Build a 10×20 block positioned at x ∈ [0, 0.45], y ∈ [0, 0.95].
    let mut particles = ParticleSet::with_capacity(200);
    for ix in 0..10_usize {
        for iy in 0..20_usize {
            let pos = Vec3::new(
                ix as f64 * spacing + spacing * 0.5,
                iy as f64 * spacing + spacing * 0.5,
                0.0,
            );
            particles.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
        }
    }

    let initial_count = particles.len();
    let initial_cx: f64 =
        particles.positions.iter().map(|p| p.x).sum::<f64>() / initial_count as f64;
    let initial_cy: f64 =
        particles.positions.iter().map(|p| p.y).sum::<f64>() / initial_count as f64;
    println!(
        "Initial: {} particles, centroid = ({:.3}, {:.3})",
        initial_count, initial_cx, initial_cy
    );

    // Run simulation.
    for _ in 0..n_steps {
        let neighbors = SpatialHash::find_all_neighbors(&particles.positions, 2.0 * h);
        wcsph_step(&mut particles, &neighbors, &kernel, &params, dt, gravity);

        // Simple floor bounce (y >= 0).
        for (pos, vel) in particles
            .positions
            .iter_mut()
            .zip(particles.velocities.iter_mut())
        {
            if pos.y < 0.0 {
                pos.y = 0.0;
                if vel.y < 0.0 {
                    vel.y = -vel.y * 0.3;
                }
            }
        }
    }

    let final_count = particles.len();
    let final_cx: f64 = particles.positions.iter().map(|p| p.x).sum::<f64>() / final_count as f64;
    let final_cy: f64 = particles.positions.iter().map(|p| p.y).sum::<f64>() / final_count as f64;

    println!(
        "Final:   {} particles, centroid = ({:.3}, {:.3})",
        final_count, final_cx, final_cy
    );
    println!("Particle count: {}", final_count);

    if final_count == initial_count {
        println!("PASS: particle count conserved ({})", final_count);
    } else {
        println!(
            "FAIL: particle count changed {} -> {}",
            initial_count, final_count
        );
    }
}
