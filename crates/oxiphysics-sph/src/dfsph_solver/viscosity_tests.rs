// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

// Regression tests for the Morris (1997) viscosity of `VelocityAdvection`.
//
// The viscous term must (a) oppose the relative velocity of a pair, (b) carry
// the dynamic viscosity `μ = ρ ν` so that the result is a real acceleration in
// m/s², and (c) reduce to `ν ∇²v` for smooth velocity fields.

use super::*;

/// Smoothing length used throughout (kernel support `2h`).
const H: f64 = 0.1;

/// Closed-form derivative `dW/dr` of the Wendland C2 kernel with support `2h`,
/// written independently of the crate's kernel code:
/// `W = 21/(16π h³) (1 − q)⁴ (1 + 4q)` with `q = r/(2h)`.
fn wendland_c2_dw_dr(r: f64, h: f64) -> f64 {
    let q = r / (2.0 * h);
    if q >= 1.0 {
        return 0.0;
    }
    let factor = 21.0 / (16.0 * std::f64::consts::PI * h * h * h);
    factor * (-20.0 * q * (1.0 - q).powi(3)) / (2.0 * h)
}

/// Neighbour entries of particle `i`: all the other particles with their
/// position, velocity, mass and density.
fn neighbors_of(particles: &[DfsphParticle], i: usize) -> Vec<NeighborEntry> {
    particles
        .iter()
        .enumerate()
        .filter(|(j, _)| *j != i)
        .map(|(_, q)| (q.position, q.velocity, q.mass, q.density))
        .collect()
}

/// Neighbour data of every particle (see [`neighbors_of`]).
fn all_pairs_neighbor_data(particles: &[DfsphParticle]) -> Vec<Vec<NeighborEntry>> {
    (0..particles.len())
        .map(|i| neighbors_of(particles, i))
        .collect()
}

/// Two particles `0.5 h` apart in `y`, moving at `±1 m/s` along `x`.
fn shear_pair() -> Vec<DfsphParticle> {
    let mut lower = DfsphParticle::new([0.0, 0.0, 0.0], 0.125, 1000.0);
    lower.velocity = [1.0, 0.0, 0.0];
    let mut upper = DfsphParticle::new([0.0, 0.5 * H, 0.0], 0.125, 1000.0);
    upper.velocity = [-1.0, 0.0, 0.0];
    vec![lower, upper]
}

fn kinetic_energy(particles: &[DfsphParticle]) -> f64 {
    particles.iter().map(DfsphParticle::kinetic_energy).sum()
}

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

#[test]
fn test_viscosity_accel_opposes_the_relative_velocity_of_a_pair() {
    let adv = VelocityAdvection::new([0.0; 3], 1e-2, H, KernelType::WendlandC2);
    let particles = shear_pair();
    let data = all_pairs_neighbor_data(&particles);
    let a_lower = adv.viscosity_accel(
        particles[0].position,
        particles[0].velocity,
        particles[0].density,
        &data[0],
    );
    let a_upper = adv.viscosity_accel(
        particles[1].position,
        particles[1].velocity,
        particles[1].density,
        &data[1],
    );
    // The lower particle moves at +1 m/s: viscosity must slow it down.
    assert!(a_lower[0] < 0.0, "lower particle: a = {a_lower:?}");
    assert!(a_upper[0] > 0.0, "upper particle: a = {a_upper:?}");
    for a in [a_lower, a_upper] {
        assert!(a[1].abs() < 1e-12 && a[2].abs() < 1e-12, "a = {a:?}");
    }
    // Equal and opposite for equal masses.
    assert!((a_lower[0] + a_upper[0]).abs() < 1e-12 * a_lower[0].abs());
}

/// `a_i = ν (ρ_i + ρ_j) m_j / (ρ_i ρ_j) · (r dW/dr) / (r² + η²) · (v_i − v_j)`
/// evaluated with an independent closed-form kernel derivative; the missing
/// density factor of the old expression would be off by a factor ρ = 1000.
#[test]
fn test_viscosity_accel_has_the_morris_magnitude_in_metres_per_second_squared() {
    let nu = 1e-2;
    let adv = VelocityAdvection::new([0.0; 3], nu, H, KernelType::WendlandC2);
    let particles = shear_pair();
    let data = all_pairs_neighbor_data(&particles);
    let a = adv.viscosity_accel(
        particles[0].position,
        particles[0].velocity,
        particles[0].density,
        &data[0],
    );
    let (rho, mass, r) = (1000.0, 0.125, 0.5 * H);
    let eta_sq = 0.01 * H * H;
    let expected = nu * (rho + rho) * mass / (rho * rho) * wendland_c2_dw_dr(r, H) * r
        / (r * r + eta_sq)
        * (1.0 - (-1.0));
    assert!(expected < 0.0);
    assert!(
        (a[0] - expected).abs() < 1e-9 * expected.abs(),
        "a_x = {}, expected {expected}",
        a[0]
    );
    // About 0.42 m/s² for a 2 m/s relative velocity. The previous expression
    // lacked the factor (ρ_i + ρ_j) / 2 = 1000 kg/m³ (it gave ~4e-4 m/s²) and
    // had the opposite sign.
    assert!(a[0].abs() > 0.3 && a[0].abs() < 0.6, "a_x = {}", a[0]);
}

/// With zero gravity the pair loses kinetic energy in one advection step
/// (the old sign made it grow) and conserves momentum.
#[test]
fn test_advect_velocity_dissipates_the_kinetic_energy_of_a_shear_pair() {
    let adv = VelocityAdvection::new([0.0; 3], 1e-2, H, KernelType::WendlandC2);
    let mut particles = shear_pair();
    let data = all_pairs_neighbor_data(&particles);
    let before = kinetic_energy(&particles);
    adv.advect_velocity(&mut particles, &data, 1e-2);
    let after = kinetic_energy(&particles);
    assert!(after < before, "kinetic energy grew: {before} -> {after}");
    assert!(
        before - after > 1e-6 * before,
        "dissipation must be measurable: {before} -> {after}"
    );
    // Both particles are slowed down symmetrically; momentum stays zero.
    assert!(particles[0].velocity[0] < 1.0 && particles[0].velocity[0] > 0.0);
    assert!(particles[1].velocity[0] > -1.0 && particles[1].velocity[0] < 0.0);
    let momentum: f64 = particles.iter().map(|p| p.mass * p.velocity[0]).sum();
    assert!(momentum.abs() < 1e-14);
}

/// Over many steps the relative velocity decays monotonically towards zero and
/// never changes sign (explicit Euler with a stable step).
#[test]
fn test_advect_velocity_relaxes_the_relative_velocity_monotonically() {
    let adv = VelocityAdvection::new([0.0; 3], 1e-2, H, KernelType::WendlandC2);
    let mut particles = shear_pair();
    let mut previous = 2.0f64;
    for _ in 0..200 {
        let data = all_pairs_neighbor_data(&particles);
        adv.advect_velocity(&mut particles, &data, 1e-1);
        let relative = particles[0].velocity[0] - particles[1].velocity[0];
        assert!(relative > 0.0, "relative velocity changed sign: {relative}");
        assert!(relative < previous, "relative velocity must shrink");
        previous = relative;
    }
    assert!(previous < 1.9, "relaxation too weak: {previous}");
}

/// For a smooth shear flow `v_x = U + s y + c y²` the Morris sum is `ν ∇²v`,
/// i.e. `a_x = 2 c ν` at an interior lattice particle (the linear part cancels
/// by symmetry).
#[test]
fn test_viscosity_accel_of_a_smooth_flow_is_nu_times_the_laplacian() {
    let n = 7usize;
    let dx = 0.05;
    let h = 1.3 * dx;
    let nu = 1e-3;
    let (u0, slope, curvature) = (0.3, 2.0, 4.0);
    let mass = 1000.0 * dx * dx * dx;
    let mut particles = Vec::new();
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let pos = [ix as f64 * dx, iy as f64 * dx, iz as f64 * dx];
                let mut p = DfsphParticle::new(pos, mass, 1000.0);
                let y = pos[1] - 3.0 * dx;
                p.velocity = [u0 + slope * y + curvature * y * y, 0.0, 0.0];
                particles.push(p);
            }
        }
    }
    let centre = 3 * n * n + 3 * n + 3;
    let adv = VelocityAdvection::new([0.0; 3], nu, h, KernelType::WendlandC2);
    let neighbors = neighbors_of(&particles, centre);
    let a = adv.viscosity_accel(
        particles[centre].position,
        particles[centre].velocity,
        particles[centre].density,
        &neighbors,
    );
    // The lattice sum and the `η² = 0.01 h²` regularisation of the Brookshaw
    // form cost a few per cent against the continuum value; a missing density
    // factor (1000x) or a wrong pair factor (2x) is far outside this band.
    let expected = 2.0 * curvature * nu;
    assert!(
        (a[0] - expected).abs() < 0.05 * expected,
        "a_x = {}, expected {expected}",
        a[0]
    );
    assert!(a[1].abs() < 1e-9 * expected && a[2].abs() < 1e-9 * expected);
}

/// A uniform translation or a rigid-body-like constant velocity field is not
/// decelerated by viscosity.
#[test]
fn test_viscosity_accel_vanishes_for_a_uniform_velocity() {
    let adv = VelocityAdvection::new([0.0; 3], 1e-2, H, KernelType::WendlandC2);
    let mut particles = shear_pair();
    for p in particles.iter_mut() {
        p.velocity = [0.7, -0.2, 0.1];
    }
    let data = all_pairs_neighbor_data(&particles);
    let a = adv.viscosity_accel(
        particles[0].position,
        particles[0].velocity,
        particles[0].density,
        &data[0],
    );
    assert_eq!(a, [0.0; 3]);
}

/// Random cloud with unequal masses and densities: the pair exchange is
/// antisymmetric (momentum conserved) and dissipative (power never positive),
/// for both Wendland kernels.
#[test]
fn test_viscosity_conserves_momentum_and_only_dissipates_for_random_clouds() {
    for kernel in [KernelType::WendlandC2, KernelType::WendlandC4] {
        let mut state = 0xD1B5_4A32_D192_ED03_u64;
        let mut particles = Vec::new();
        for _ in 0..40 {
            let pos = [
                0.2 * lcg_unit(&mut state),
                0.2 * lcg_unit(&mut state),
                0.2 * lcg_unit(&mut state),
            ];
            let mass = 0.05 + 0.1 * lcg_unit(&mut state);
            let mut p = DfsphParticle::new(pos, mass, 800.0 + 400.0 * lcg_unit(&mut state));
            for d in 0..3 {
                p.velocity[d] = 2.0 * (lcg_unit(&mut state) - 0.5);
            }
            particles.push(p);
        }
        let adv = VelocityAdvection::new([0.0; 3], 5e-3, H, kernel.clone());
        let data = all_pairs_neighbor_data(&particles);
        let mut momentum_rate = [0.0f64; 3];
        let mut scale = 0.0f64;
        let mut power = 0.0f64;
        for (i, p) in particles.iter().enumerate() {
            let a = adv.viscosity_accel(p.position, p.velocity, p.density, &data[i]);
            for d in 0..3 {
                momentum_rate[d] += p.mass * a[d];
                scale += (p.mass * a[d]).abs();
                power += p.mass * p.velocity[d] * a[d];
            }
        }
        assert!(scale > 1e-3, "{kernel:?}: viscosity must act: {scale}");
        for (d, rate) in momentum_rate.iter().enumerate() {
            assert!(
                rate.abs() < 1e-12 * scale,
                "{kernel:?}: momentum rate {rate} in direction {d} (scale {scale})"
            );
        }
        assert!(
            power < 0.0,
            "{kernel:?}: viscous power must be negative: {power}"
        );
    }
}

/// A particle with a vanishing density has no defined viscous acceleration; the
/// guard returns zero instead of a division by zero.
#[test]
fn test_viscosity_accel_guards_against_zero_density() {
    let adv = VelocityAdvection::new([0.0; 3], 1e-2, H, KernelType::WendlandC2);
    let particles = shear_pair();
    let data = all_pairs_neighbor_data(&particles);
    let a = adv.viscosity_accel(particles[0].position, particles[0].velocity, 0.0, &data[0]);
    assert_eq!(a, [0.0; 3]);
}
