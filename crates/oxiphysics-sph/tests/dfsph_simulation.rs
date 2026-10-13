// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! End-to-end checks of the DFSPH pressure solver driven through
//! [`SphSimulation`] with `SolverType::Dfsph`.

use oxiphysics_core::math::Vec3;
use oxiphysics_sph::kernel::CubicSplineKernel;
use oxiphysics_sph::particle::SphParticle;
use oxiphysics_sph::simulation::{SolverType, SphSimulation, SphSimulationParams};

const RHO0: f64 = 1000.0;
const REST_SPACING: f64 = 0.05;

/// `n`³ block of particles with lattice spacing `spacing`; the particle mass
/// belongs to the rest spacing, so `spacing < REST_SPACING` is a compressed
/// block.  The kernel support radius is the usual `2 Δx`.
fn dfsph_block(n: usize, spacing: f64) -> SphSimulation {
    let params = SphSimulationParams {
        gravity: Vec3::zeros(),
        smoothing_length: REST_SPACING,
        rest_density: RHO0,
        solver_type: SolverType::Dfsph,
        ..Default::default()
    };
    let mass = RHO0 * REST_SPACING.powi(3);
    let mut sim = SphSimulation::new(params);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let pos = Vec3::new(i as f64 * spacing, j as f64 * spacing, k as f64 * spacing);
                sim.particles
                    .add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
            }
        }
    }
    sim
}

#[test]
fn dfsph_simulation_relieves_the_compression_of_a_dense_block() {
    let kernel = CubicSplineKernel;
    let n = 5;
    let centre = 2 * 25 + 2 * 5 + 2;
    // Spacing 0.0495 against the rest spacing 0.05: the block starts ~3 % too dense.
    let mut sim = dfsph_block(n, 0.0495);

    sim.step(1.0e-3, &kernel);
    let initial_max = sim.particles.densities.iter().cloned().fold(0.0, f64::max);
    assert!(
        initial_max > 1.02 * RHO0,
        "block must start compressed: {initial_max}"
    );

    sim.step(1.0e-3, &kernel);
    // The pressure solve of the first step has pushed the particles apart: no
    // particle is compressed any more at the start of the second step, and the
    // block has not been blown apart either.
    let max_density = sim.particles.densities.iter().cloned().fold(0.0, f64::max);
    assert!(
        max_density < 1.01 * RHO0,
        "compression must be relieved: max density {max_density}"
    );
    let centre_density = sim.particles.densities[centre];
    assert!(
        centre_density > 0.8 * RHO0,
        "centre density {centre_density} should stay near the rest density {RHO0}"
    );
    for v in &sim.particles.velocities {
        assert!(v.iter().all(|c| c.is_finite()));
    }
    for p in &sim.particles.positions {
        assert!(p.iter().all(|c| c.is_finite() && c.abs() < 1.0));
    }
}

#[test]
fn dfsph_simulation_conserves_momentum_and_leaves_a_relaxed_block_alone() {
    let kernel = CubicSplineKernel;

    // A block at its rest spacing has (almost) no compression to remove: the
    // solver returns immediately and the particles stay put.
    let mut rest = dfsph_block(4, REST_SPACING);
    let start: Vec<Vec3> = rest.particles.positions.clone();
    rest.step(1.0e-3, &kernel);
    let drift = rest
        .particles
        .positions
        .iter()
        .zip(start.iter())
        .map(|(p, q)| (p - q).norm())
        .fold(0.0_f64, f64::max);
    assert!(drift < 1.0e-9, "relaxed block moved by {drift}");

    // Pressure corrections are pairwise antisymmetric: they conserve momentum.
    let mut dense = dfsph_block(5, 0.0495);
    for _ in 0..3 {
        dense.step(1.0e-3, &kernel);
    }
    let mut momentum = Vec3::zeros();
    let mut scale = 0.0;
    for (m, v) in dense
        .particles
        .masses
        .iter()
        .zip(dense.particles.velocities.iter())
    {
        momentum += *m * v;
        scale += *m * v.norm();
    }
    assert!(scale > 0.0, "the dense block must have been set in motion");
    assert!(momentum.norm() < 1.0e-9 * scale, "momentum {momentum:?}");
}
