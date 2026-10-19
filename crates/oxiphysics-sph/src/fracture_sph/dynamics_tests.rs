// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

// Regression tests for the dynamics of the SPH fracture simulation: the
// momentum and continuity integration (stress must move the particles and the
// density must evolve), the Jaumann co-rotational stress update, the single
// application of damage, and the Monaghan artificial viscosity.

use std::f64::consts::PI;

use super::functions::*;
use super::types::*;

const E_MOD: f64 = 70e9;
const RHO: f64 = 2500.0;

/// Damage model whose critical strain is so large that no damage ever starts.
fn elastic_model(poisson_ratio: f64) -> ContinuumDamageModel {
    ContinuumDamageModel::new_isotropic(1.0e3, 50.0, E_MOD, poisson_ratio, 0.01)
}

fn fracture_energy() -> FractureEnergy {
    FractureEnergy::new(50.0, E_MOD, 1e6, 1e-5, true)
}

fn particle(position: [f64; 3], mass: f64, h: f64) -> SphFractureParticle {
    SphFractureParticle::new(position, [0.0; 3], mass, h, RHO)
}

fn simulation(
    particles: Vec<SphFractureParticle>,
    poisson_ratio: f64,
    dt: f64,
    h: f64,
) -> SphFractureSimulation {
    SphFractureSimulation::new(
        particles,
        elastic_model(poisson_ratio),
        fracture_energy(),
        dt,
        h,
    )
}

/// `nx × ny × nz` lattice with spacing `dx`, mass `ρ dx³` and smoothing length `h`.
fn lattice(nx: usize, ny: usize, nz: usize, dx: f64, h: f64) -> Vec<SphFractureParticle> {
    let mass = RHO * dx * dx * dx;
    let mut particles = Vec::with_capacity(nx * ny * nz);
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                particles.push(particle(
                    [ix as f64 * dx, iy as f64 * dx, iz as f64 * dx],
                    mass,
                    h,
                ));
            }
        }
    }
    particles
}

/// Closed-form derivative `dW/dr` of the Wendland C2 kernel with support `2h`,
/// written independently of the crate's kernel code.
fn wendland_dw_dr(r: f64, h: f64) -> f64 {
    let q = r / (2.0 * h);
    if q >= 1.0 {
        return 0.0;
    }
    let factor = 21.0 / (16.0 * PI * h * h * h);
    factor * (-20.0 * q * (1.0 - q).powi(3)) / (2.0 * h)
}

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

fn total_momentum(sim: &SphFractureSimulation) -> [f64; 3] {
    let mut p = [0.0f64; 3];
    for q in &sim.particles {
        for (sum, &v) in p.iter_mut().zip(q.velocity.iter()) {
            *sum += q.mass * v;
        }
    }
    p
}

// ---------------------------------------------------------------------------
// Stress acts on the motion
// ---------------------------------------------------------------------------

/// Two particles with the stress `σ_xx = ±10⁸ Pa`: tension pulls them together,
/// compression pushes them apart, with the analytic SPH acceleration
/// `a₀ = m (σ₀/ρ₀² + σ₁/ρ₁²) ∂W/∂x`.
#[test]
fn test_stress_moves_the_particles_tension_attracts_compression_repels() {
    let (dx, h, dt) = (0.01, 0.015, 1e-7);
    let mass = RHO * dx * dx * dx;
    for stress in [1.0e8, -1.0e8] {
        let mut particles = vec![
            particle([0.0, 0.0, 0.0], mass, h),
            particle([dx, 0.0, 0.0], mass, h),
        ];
        for p in particles.iter_mut() {
            p.stress[0][0] = stress;
        }
        let mut sim = simulation(particles, 0.3, dt, h);
        sim.update_momentum_and_density();
        // ∇₀W = dW/dr (x₀ − x₁)/r points towards +x because dW/dr < 0.
        let grad_x = -wendland_dw_dr(dx, h);
        let expected = mass * (2.0 * stress / (RHO * RHO)) * grad_x * dt;
        let (v0, v1) = (sim.particles[0].velocity, sim.particles[1].velocity);
        assert!(
            (v0[0] - expected).abs() < 1e-9 * expected.abs(),
            "σ = {stress}: v0 = {}, expected {expected}",
            v0[0]
        );
        assert!((v0[0] + v1[0]).abs() < 1e-12 * expected.abs());
        if stress > 0.0 {
            assert!(v0[0] > 0.0 && v1[0] < 0.0, "tension must attract");
        } else {
            assert!(v0[0] < 0.0 && v1[0] > 0.0, "compression must repel");
        }
        assert!(v0[1].abs() < 1e-12 * expected.abs() && v0[2].abs() < 1e-12 * expected.abs());
    }
}

/// A full step moves the particles with the new velocities (symplectic
/// Euler) and the density follows the continuity equation: two particles
/// pulled together by tension get denser from the second step on.
#[test]
fn test_step_integrates_velocity_density_and_position() {
    let (dx, h, dt) = (0.01, 0.015, 1e-7);
    let mass = RHO * dx * dx * dx;
    let mut particles = vec![
        particle([0.0, 0.0, 0.0], mass, h),
        particle([dx, 0.0, 0.0], mass, h),
    ];
    for p in particles.iter_mut() {
        p.stress[0][0] = 1.0e8;
    }
    let mut sim = simulation(particles, 0.3, dt, h);
    sim.art_visc_alpha = 0.0;
    sim.art_visc_beta = 0.0;
    assert_eq!(sim.kinetic_energy(), 0.0);
    sim.step();
    // The stress has produced motion (previously the velocity never changed).
    let v_after_one = sim.particles[0].velocity[0];
    assert!(v_after_one > 0.0, "v = {v_after_one}");
    assert!(sim.kinetic_energy() > 0.0);
    // Symplectic Euler: positions use the updated velocities.
    assert!(
        (sim.particles[0].position[0] - v_after_one * dt).abs() < 1e-12 * (v_after_one * dt),
        "x = {}",
        sim.particles[0].position[0]
    );
    // The velocity was zero when the first continuity rate was formed.
    assert_eq!(sim.particles[0].density, RHO);
    sim.step();
    assert!(
        sim.particles[0].density > RHO && sim.particles[1].density > RHO,
        "approaching particles must get denser: {} {}",
        sim.particles[0].density,
        sim.particles[1].density
    );
}

/// A uniform velocity field does not change the density, a contraction
/// `v = −k (x − x_c)` raises it at the rate `−ρ ∇·v = 3 k ρ`, an expansion
/// lowers it.
#[test]
fn test_density_follows_the_continuity_equation() {
    let (dx, h, dt) = (0.01, 0.015, 1e-8);
    let n = 7;
    let centre = (3 * n + 3) * n + 3;
    let build = |k: f64, uniform: bool| {
        let mut particles = lattice(n, n, n, dx, h);
        for p in particles.iter_mut() {
            p.velocity = if uniform {
                [0.4, -0.2, 0.1]
            } else {
                [
                    -k * (p.position[0] - 3.0 * dx),
                    -k * (p.position[1] - 3.0 * dx),
                    -k * (p.position[2] - 3.0 * dx),
                ]
            };
        }
        let mut sim = simulation(particles, 0.3, dt, h);
        sim.art_visc_alpha = 0.0;
        sim.art_visc_beta = 0.0;
        sim.update_momentum_and_density();
        sim
    };
    let uniform = build(0.0, true);
    assert_eq!(uniform.particles[centre].density, RHO);
    let k = 50.0;
    let contraction = build(k, false);
    let rate = (contraction.particles[centre].density - RHO) / dt;
    assert!(
        (rate - 3.0 * k * RHO).abs() < 0.03 * 3.0 * k * RHO,
        "dρ/dt = {rate}, expected {}",
        3.0 * k * RHO
    );
    let expansion = build(-k, false);
    let rate = (expansion.particles[centre].density - RHO) / dt;
    assert!(
        (rate + 3.0 * k * RHO).abs() < 0.03 * 3.0 * k * RHO,
        "dρ/dt = {rate}"
    );
}

// ---------------------------------------------------------------------------
// Momentum conservation
// ---------------------------------------------------------------------------

/// Random stresses, damages, masses, densities and velocities (so the
/// artificial viscosity is active for approaching pairs): the total momentum
/// is conserved by the integration to rounding.
#[test]
fn test_momentum_integration_conserves_total_momentum() {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let (dx, h, dt) = (0.01, 0.015, 1e-7);
    let mut particles = lattice(4, 4, 4, dx, h);
    for p in particles.iter_mut() {
        for d in 0..3 {
            p.position[d] += 0.2 * dx * (lcg_unit(&mut state) - 0.5);
            p.velocity[d] = 20.0 * (lcg_unit(&mut state) - 0.5);
        }
        p.mass *= 0.5 + lcg_unit(&mut state);
        p.density *= 0.9 + 0.2 * lcg_unit(&mut state);
        p.damage = lcg_unit(&mut state);
        let mut entries = [0.0f64; 6];
        for entry in entries.iter_mut() {
            *entry = 2.0e7 * (lcg_unit(&mut state) - 0.5);
        }
        let [xx, yy, zz, xy, xz, yz] = entries;
        p.stress = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
    }
    let mut sim = simulation(particles, 0.3, dt, h);
    let before = total_momentum(&sim);
    let velocities_before: Vec<[f64; 3]> = sim.particles.iter().map(|p| p.velocity).collect();
    sim.update_momentum_and_density();
    let after = total_momentum(&sim);
    let mut exchanged = 0.0f64;
    for (p, v0) in sim.particles.iter().zip(&velocities_before) {
        for (&v, &v_before) in p.velocity.iter().zip(v0.iter()) {
            exchanged += (p.mass * (v - v_before)).abs();
        }
    }
    assert!(
        exchanged > 1e-6,
        "the particles must exchange momentum: {exchanged}"
    );
    for d in 0..3 {
        assert!(
            (after[d] - before[d]).abs() < 1e-12 * exchanged,
            "momentum changed by {} in direction {d} (exchanged {exchanged})",
            after[d] - before[d]
        );
    }
    // Self-check that the viscosity really took part: without it the result
    // differs.
    let mut inviscid = sim.clone();
    for (p, v0) in inviscid.particles.iter_mut().zip(&velocities_before) {
        p.velocity = *v0;
    }
    inviscid.art_visc_alpha = 0.0;
    inviscid.art_visc_beta = 0.0;
    inviscid.update_momentum_and_density();
    let differs = sim
        .particles
        .iter()
        .zip(&inviscid.particles)
        .any(|(a, b)| (a.velocity[0] - b.velocity[0]).abs() > 1e-9);
    assert!(differs, "the artificial viscosity had no effect");
}

// ---------------------------------------------------------------------------
// Damage enters exactly once
// ---------------------------------------------------------------------------

/// The nominal stress is `(1 − D) σ`: halving the stiffness halves the force,
/// and a fully damaged pair exerts no force at all, however large the stored
/// (undamaged) stress is.
#[test]
fn test_damage_scales_the_force_by_one_minus_d() {
    let (dx, h, dt) = (0.01, 0.015, 1e-7);
    let mass = RHO * dx * dx * dx;
    let accelerate = |damage: f64| {
        let mut particles = vec![
            particle([0.0, 0.0, 0.0], mass, h),
            particle([dx, 0.0, 0.0], mass, h),
        ];
        for p in particles.iter_mut() {
            p.stress[0][0] = 1.0e9;
            p.damage = damage;
        }
        let mut sim = simulation(particles, 0.3, dt, h);
        sim.update_momentum_and_density();
        sim.particles[0].velocity[0]
    };
    let intact = accelerate(0.0);
    let half = accelerate(0.5);
    let broken = accelerate(1.0);
    assert!(intact > 0.0);
    assert!(
        (half - 0.5 * intact).abs() < 1e-12 * intact,
        "D = 0.5: {half} vs intact {intact}"
    );
    assert_eq!(broken, 0.0);
}

/// The stored stress is the stress of the undamaged solid: `update_stress`
/// must not scale its increments by `1 − D` (the force already does).
#[test]
fn test_update_stress_does_not_apply_damage_a_second_time() {
    let (dx, h, dt) = (0.01, 0.015, 1e-6);
    let n = 7;
    let strain_rate = 10.0;
    let centre = (3 * n + 3) * n + 3;
    let run = |damage: f64| {
        let mut particles = lattice(n, n, n, dx, h);
        for p in particles.iter_mut() {
            p.velocity = [strain_rate * p.position[0], 0.0, 0.0];
            p.damage = damage;
        }
        let mut sim = simulation(particles, 0.3, dt, h);
        sim.update_stress();
        sim.particles[centre].stress
    };
    let intact = run(0.0);
    let damaged = run(0.5);
    for a in 0..3 {
        for b in 0..3 {
            assert!(
                (intact[a][b] - damaged[a][b]).abs() < 1e-12 * intact[0][0].abs(),
                "σ[{a}][{b}]: {} vs {}",
                intact[a][b],
                damaged[a][b]
            );
        }
    }
    assert!(
        intact[0][0] > 1.0e3,
        "the stress must have grown: {}",
        intact[0][0]
    );
}

// ---------------------------------------------------------------------------
// Jaumann co-rotation
// ---------------------------------------------------------------------------

/// A lattice that rotates rigidly about `z` with angular velocity `w` carries
/// its stress along: `σ = diag(s, 0, 0)` acquires `σ_xy = +s w dt` (the sign of
/// `R σ Rᵀ`), while the elastic response to the (zero) strain rate vanishes.
#[test]
fn test_update_stress_rotates_the_stress_of_a_rigidly_rotating_body() {
    let (dx, h) = (0.01, 0.015);
    let (w, dt, s0) = (100.0, 1e-6, 1.0e6);
    let n = 7;
    let centre = (3 * n + 3) * n + 3;
    let mut particles = lattice(n, n, n, dx, h);
    let c = 3.0 * dx;
    for p in particles.iter_mut() {
        p.velocity = [-w * (p.position[1] - c), w * (p.position[0] - c), 0.0];
        p.stress[0][0] = s0;
    }
    let mut sim = simulation(particles, 0.3, dt, h);
    sim.update_stress();
    let sigma = sim.particles[centre].stress;
    let expected_xy = s0 * w * dt;
    assert!(
        (sigma[0][1] - expected_xy).abs() < 0.03 * expected_xy,
        "σ_xy = {}, expected {expected_xy}",
        sigma[0][1]
    );
    assert!((sigma[1][0] - sigma[0][1]).abs() < 1e-9 * expected_xy);
    // Rigid motion has no strain rate: no elastic stress is generated.
    assert!((sigma[0][0] - s0).abs() < 1e-6 * expected_xy);
    assert!(sigma[1][1].abs() < 1e-6 * expected_xy);
    assert!(sigma[2][2].abs() < 1e-6 * expected_xy);
    assert!(sigma[0][2].abs() < 1e-6 * expected_xy && sigma[1][2].abs() < 1e-6 * expected_xy);
}

// ---------------------------------------------------------------------------
// Artificial viscosity
// ---------------------------------------------------------------------------

/// A stress-free pair approaching at 2 m/s is decelerated by the Monaghan
/// viscosity with the analytic `Π`, a separating pair is not touched, and
/// switching `α = β = 0` removes it.
#[test]
fn test_artificial_viscosity_decelerates_approaching_pairs_only() {
    let (dx, h, dt) = (0.01, 0.015, 1e-9);
    let mass = RHO * dx * dx * dx;
    let make = |v0: f64, v1: f64| {
        let mut particles = vec![
            particle([0.0, 0.0, 0.0], mass, h),
            particle([dx, 0.0, 0.0], mass, h),
        ];
        particles[0].velocity = [v0, 0.0, 0.0];
        particles[1].velocity = [v1, 0.0, 0.0];
        // ν = 0 makes the P-wave speed sqrt(E/ρ).
        simulation(particles, 0.0, dt, h)
    };
    // Approaching.
    let mut sim = make(1.0, -1.0);
    sim.update_momentum_and_density();
    let c = (E_MOD / RHO).sqrt();
    let eta_sq = 0.01 * h * h;
    let mu = h * (2.0 * -dx) / (dx * dx + eta_sq);
    let pi_ij = (-c * mu + 2.0 * mu * mu) / RHO;
    let grad_x = -wendland_dw_dr(dx, h);
    let expected = -mass * pi_ij * grad_x * dt;
    let dv0 = sim.particles[0].velocity[0] - 1.0;
    assert!(expected < 0.0);
    assert!(
        (dv0 - expected).abs() < 1e-9 * expected.abs(),
        "Δv0 = {dv0}, expected {expected}"
    );
    let dv1 = sim.particles[1].velocity[0] + 1.0;
    assert!((dv0 + dv1).abs() < 1e-9 * expected.abs());
    // Separating: Π = 0.
    let mut sim = make(-1.0, 1.0);
    sim.update_momentum_and_density();
    assert_eq!(sim.particles[0].velocity[0], -1.0);
    assert_eq!(sim.particles[1].velocity[0], 1.0);
    // Disabled.
    let mut sim = make(1.0, -1.0);
    sim.art_visc_alpha = 0.0;
    sim.art_visc_beta = 0.0;
    sim.update_momentum_and_density();
    assert_eq!(sim.particles[0].velocity[0], 1.0);
}

// ---------------------------------------------------------------------------
// Elastic wave speed
// ---------------------------------------------------------------------------

/// Momentum centroid `Σ mᵢ vᵢ xᵢ / Σ mᵢ vᵢ` of the particles `indices`.
///
/// The centroid of a pulse moves at the long-wavelength group velocity, so it
/// measures the wave speed without the lag that the dispersion of the SPH
/// lattice gives to the peak of a short pulse.
fn momentum_centroid(sim: &SphFractureSimulation, indices: &[usize]) -> f64 {
    let mut momentum = 0.0;
    let mut moment = 0.0;
    for &i in indices {
        let p = &sim.particles[i];
        momentum += p.mass * p.velocity[0];
        moment += p.mass * p.velocity[0] * p.position[0];
    }
    moment / momentum
}

/// A velocity pulse in an elastic bar splits into two pulses that travel at the
/// bar speed `√(E/ρ)`. The bar is a 3-D lattice (`m = ρ dx³`, `h = dx`, kernel
/// support `2 dx`) with `ν = 0`, so the uniaxial-strain modulus `λ + 2μ` equals
/// `E`; the speed is measured on the centre line, whose stencil is complete.
///
/// The SPH velocity gradient of the particles on the free lateral faces of the
/// 5 × 5 cross-section has a truncated stencil, which makes the discrete bar a
/// few per cent softer than the ideal one (about 6 % at `h = dx`, 15 % at
/// `h = 1.3 dx`, falling to 1.5 % for a 9 × 9 cross-section). The tolerance of
/// 10 % is far below the factor 2 or the missing motion of the old
/// implementation, which never moved the particles at all.
#[test]
fn test_velocity_pulse_in_a_lattice_bar_travels_at_the_bar_speed() {
    let dx = 0.01;
    let h = dx;
    let (nx, ny, nz) = (80usize, 5usize, 5usize);
    let c = (E_MOD / RHO).sqrt();
    let dt = 0.25 * h / c;
    let mut particles = lattice(nx, ny, nz, dx, h);
    let (x0, width, v0) = (0.4, 5.0 * dx, 0.01);
    for p in particles.iter_mut() {
        let d = p.position[0] - x0;
        p.velocity[0] = v0 * (-0.5 * d * d / (width * width)).exp();
    }
    let mut sim = simulation(particles, 0.0, dt, h);
    sim.art_visc_alpha = 0.0;
    sim.art_visc_beta = 0.0;
    let momentum_before = total_momentum(&sim);
    let steps = 60;
    sim.run(steps);
    let t = steps as f64 * dt;
    // Centre line of the bar: iy = iz = 2.
    let line: Vec<usize> = (0..nx).map(|ix| (ix * ny + 2) * nz + 2).collect();
    let travelled = c * t;
    let right: Vec<usize> = line
        .iter()
        .copied()
        .filter(|&i| sim.particles[i].position[0] > x0 + 0.25 * travelled)
        .collect();
    let left: Vec<usize> = line
        .iter()
        .copied()
        .filter(|&i| sim.particles[i].position[0] < x0 - 0.25 * travelled)
        .collect();
    let right_centroid = momentum_centroid(&sim, &right);
    let left_centroid = momentum_centroid(&sim, &left);
    let speed = (right_centroid - x0) / t;
    assert!(
        (speed - c).abs() < 0.10 * c,
        "pulse speed {speed} m/s, bar speed {c} m/s (centroid at {right_centroid})"
    );
    // The two half-pulses are (nearly) mirror images about the launch point.
    assert!(
        ((x0 - left_centroid) - (right_centroid - x0)).abs() < 0.02 * (right_centroid - x0),
        "left {left_centroid}, right {right_centroid}"
    );
    // The total momentum is conserved by the whole step sequence.
    let momentum_after = total_momentum(&sim);
    assert!(
        (momentum_after[0] - momentum_before[0]).abs() < 1e-12 * momentum_before[0].abs(),
        "momentum {} -> {}",
        momentum_before[0],
        momentum_after[0]
    );
}

/// The same measurement on a chain of particles along `x` whose volume is
/// chosen so that the 3-D kernel reproduces the 1-D gradient exactly
/// (`Σⱼ Vⱼ (xⱼ − xᵢ) ∂ₓW = 1`): the pulse speed then equals `√(E/ρ)`.
#[test]
fn test_velocity_pulse_in_a_consistent_chain_travels_at_the_bar_speed() {
    let dx = 0.01;
    let h = 1.3 * dx;
    let n = 200usize;
    // Volume per particle from the first-order consistency condition.
    let mut moment = 0.0;
    for k in 1..=3 {
        let r = k as f64 * dx;
        moment += -2.0 * r * wendland_c2_3d_grad(r, h);
    }
    let volume = 1.0 / moment;
    let mass = RHO * volume;
    let c = (E_MOD / RHO).sqrt();
    let dt = 0.25 * h / c;
    let (x0, width, v0) = (1.0, 5.0 * dx, 0.01);
    let mut particles = Vec::with_capacity(n);
    for i in 0..n {
        let x = i as f64 * dx;
        let mut p = particle([x, 0.0, 0.0], mass, h);
        let d = x - x0;
        p.velocity[0] = v0 * (-0.5 * d * d / (width * width)).exp();
        particles.push(p);
    }
    let mut sim = simulation(particles, 0.0, dt, h);
    sim.art_visc_alpha = 0.0;
    sim.art_visc_beta = 0.0;
    let steps = 200;
    sim.run(steps);
    let t = steps as f64 * dt;
    let travelled = c * t;
    let right: Vec<usize> = (0..n)
        .filter(|&i| sim.particles[i].position[0] > x0 + 0.25 * travelled)
        .collect();
    let centroid = momentum_centroid(&sim, &right);
    let speed = (centroid - x0) / t;
    assert!(
        (speed - c).abs() < 0.02 * c,
        "pulse speed {speed} m/s, bar speed {c} m/s"
    );
    // Half of the launched amplitude survives in each pulse (weak dispersion).
    let amplitude = sim
        .particles
        .iter()
        .map(|p| p.velocity[0])
        .fold(0.0f64, f64::max);
    assert!(
        amplitude > 0.45 * v0 && amplitude < 0.52 * v0,
        "pulse amplitude {amplitude}"
    );
}
