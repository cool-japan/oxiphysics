// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the ISPH free-surface detection and pressure solve on a
//! particle block whose particles are displaced randomly from their lattice
//! sites, which is the situation of every real simulation after a few steps.

use super::*;

/// Lattice spacing of the test blocks \[m\].
const SPACING: f64 = 0.1;
/// Rest density of the test fluid \[kg/m³\].
const REST_DENSITY: f64 = 1000.0;
/// Particles per edge of the test blocks.
const EDGE: usize = 9;
/// Chebyshev ring (`max(|ix|, |iy|, |iz|)` in lattice units) of the outer layer.
const OUTER_RING: usize = 4;
/// Displacement amplitude of the disorder, in units of the lattice spacing.
const JITTER: f64 = 0.1;

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit linear congruential
/// generator, top 53 bits).
struct UnitRandom(u64);

impl UnitRandom {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
}

/// `EDGE`³ lattice centred on the origin whose sites are displaced by independent
/// uniform offsets in `[−jitter·s, jitter·s]` per coordinate.
///
/// Every particle keeps the lattice volume `s³` (mass `ρ₀ s³`).  Returns the
/// particles and the Chebyshev ring of the lattice site of every particle
/// (`OUTER_RING` for the outermost layer).
fn disordered_lattice(jitter: f64, seed: u64) -> (Vec<IsphParticle>, Vec<usize>) {
    let half = (EDGE as f64 - 1.0) / 2.0;
    let mut random = UnitRandom::new(seed);
    let mut particles = Vec::with_capacity(EDGE * EDGE * EDGE);
    let mut rings = Vec::with_capacity(EDGE * EDGE * EDGE);
    for ix in 0..EDGE {
        for iy in 0..EDGE {
            for iz in 0..EDGE {
                let site = [ix as f64 - half, iy as f64 - half, iz as f64 - half];
                let mut position = [0.0_f64; 3];
                for (coordinate, site_coordinate) in position.iter_mut().zip(site.iter()) {
                    let offset = (2.0 * random.next_unit() - 1.0) * jitter;
                    *coordinate = (site_coordinate + offset) * SPACING;
                }
                let mass = REST_DENSITY * SPACING * SPACING * SPACING;
                particles.push(IsphParticle::new(position, mass, REST_DENSITY));
                let ring = site
                    .iter()
                    .map(|c| c.abs().round() as usize)
                    .max()
                    .unwrap_or(0);
                rings.push(ring);
            }
        }
    }
    (particles, rings)
}

/// Index of the particle closest to the origin.
fn most_central(particles: &[IsphParticle]) -> usize {
    particles
        .iter()
        .enumerate()
        .map(|(i, p)| (i, vec3_dot(p.position, p.position)))
        .fold((0, f64::INFINITY), |best, candidate| {
            if candidate.1 < best.1 {
                candidate
            } else {
                best
            }
        })
        .0
}

#[test]
fn test_detect_free_surface_is_robust_to_particle_disorder() {
    // A displacement of 10 % of the spacing changes the number density of the
    // interior by several percent: a kernel-sum criterion (complete support = 1)
    // flags some interior particles of every such block, whereas the position
    // divergence (complete support = 3) still separates exactly the outer layer.
    let solver = IsphPressureSolver::new(1.2 * SPACING, IsphPressureSolverParams::default());
    for seed in 0..6_u64 {
        let (particles, rings) = disordered_lattice(JITTER, 1000 + seed);
        let flags = solver.detect_free_surface(&particles);
        assert_eq!(
            flags.iter().filter(|&&flag| flag).count(),
            EDGE * EDGE * EDGE - (EDGE - 2) * (EDGE - 2) * (EDGE - 2),
            "seed {seed}"
        );
        for (i, (&flag, &ring)) in flags.iter().zip(rings.iter()).enumerate() {
            assert_eq!(
                flag,
                ring == OUTER_RING,
                "seed {seed}: particle {i} at {:?} on lattice ring {ring}",
                particles[i].position
            );
        }
    }
}

#[test]
fn test_detect_free_surface_flags_a_lone_particle_and_a_separated_pair() {
    let solver = IsphPressureSolver::new(1.2 * SPACING, IsphPressureSolverParams::default());
    let lone = vec![IsphParticle::new([0.0; 3], 1.0, REST_DENSITY)];
    assert_eq!(solver.detect_free_surface(&lone), vec![true]);
    let pair = vec![
        IsphParticle::new([0.0; 3], 1.0, REST_DENSITY),
        IsphParticle::new([SPACING, 0.0, 0.0], 1.0, REST_DENSITY),
    ];
    assert_eq!(solver.detect_free_surface(&pair), vec![true, true]);
}

#[test]
fn test_pressure_solver_compression_on_a_disordered_block() {
    // The compressed block of the ordered test, with every particle displaced by
    // up to 10 % of the spacing: the pressure at the centre still matches the
    // Poisson solution of the cube (0.0562 L² ρ₀/Δt, L = 8 s), is positive
    // everywhere inside and vanishes on the pinned outer layer.
    let dt = 1e-3;
    let h = 1.2 * SPACING;
    let (mut particles, rings) = disordered_lattice(JITTER, 1002);
    for p in particles.iter_mut() {
        p.velocity = [-p.position[0], 0.0, 0.0];
    }
    VelocityDivergence::new(h).compute(&mut particles);
    let params = IsphPressureSolverParams {
        max_iter: 3000,
        tolerance: 1e-12,
        relaxation: 1.0,
        rest_density: REST_DENSITY,
        clamp_negative_pressure: false,
    };
    let solver = IsphPressureSolver::new(h, params);
    let (iters, _) = solver.solve(&mut particles, dt);
    assert!(iters < 3000, "solver did not converge: {iters} sweeps");

    let centre = most_central(&particles);
    let expected = 0.0562 * (8.0 * SPACING).powi(2) * REST_DENSITY / dt;
    let p_centre = particles[centre].pressure;
    assert!(
        (p_centre - expected).abs() < 0.12 * expected,
        "p(centre) = {p_centre}, cube solution {expected}"
    );
    for (i, (p, &ring)) in particles.iter().zip(rings.iter()).enumerate() {
        if ring == OUTER_RING {
            assert_eq!(p.pressure, 0.0, "pinned particle {i}");
        } else {
            assert!(
                p.pressure > 0.0,
                "compression must raise the pressure of particle {i}: {}",
                p.pressure
            );
        }
    }
}
