// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the elastic SPH force terms: the stress divergence,
//! the hourglass control and the Jaumann (co-rotational) stress update.

use std::f64::consts::PI;

use super::functions::{IDENTITY3, cubic_kernel, hourglass_control, stress_divergence};
use super::types::{ElasticParticle, ElasticSphParticle, SphElasticStress};

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

/// Two aluminium-like particles 0.04 m apart along `x` with smoothing length
/// 0.05 m, both carrying the isotropic Cauchy stress `value · I` (tension
/// positive).
fn isotropic_pair(value: f64) -> Vec<ElasticParticle> {
    let mut particles = vec![
        ElasticParticle::new([0.0, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
        ElasticParticle::new([0.04, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
    ];
    for p in particles.iter_mut() {
        for d in 0..3 {
            p.stress[d][d] = value;
        }
    }
    particles
}

/// Cubic lattice of `n³` particles of density 2700 kg/m³ with spacing `dx`,
/// smoothing length `1.3 dx` and mass `ρ dx³`.
fn lattice(n: usize, dx: f64) -> Vec<ElasticParticle> {
    let density = 2700.0;
    let mass = density * dx * dx * dx;
    let mut particles = Vec::with_capacity(n * n * n);
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                particles.push(ElasticParticle::new(
                    [ix as f64 * dx, iy as f64 * dx, iz as f64 * dx],
                    mass,
                    density,
                    1.3 * dx,
                    70e9,
                    0.33,
                ));
            }
        }
    }
    particles
}

/// `dW/dr` of the cubic-spline kernel with support `2h`, written out
/// independently of the crate's gradient routine.
fn cubic_spline_dw_dr(r: f64, h: f64) -> f64 {
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha / h * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        -alpha / h * 0.75 * (2.0 - q) * (2.0 - q)
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// stress_divergence
// ---------------------------------------------------------------------------

/// Compressed material (`σ = −p I`, `p > 0`) must push the particles apart.
#[test]
fn test_stress_divergence_compression_pushes_particles_apart() {
    let mut particles = isotropic_pair(-1.0e6);
    stress_divergence(&mut particles);
    assert!(
        particles[0].acceleration[0] < 0.0,
        "the left particle must be pushed to -x, a = {:?}",
        particles[0].acceleration
    );
    assert!(
        particles[1].acceleration[0] > 0.0,
        "the right particle must be pushed to +x, a = {:?}",
        particles[1].acceleration
    );
    for p in &particles {
        assert!(p.acceleration[1].abs() < 1e-9 && p.acceleration[2].abs() < 1e-9);
    }
}

/// Stretched material (`σ > 0`) must pull the particles together.
#[test]
fn test_stress_divergence_tension_pulls_particles_together() {
    let mut particles = isotropic_pair(1.0e6);
    stress_divergence(&mut particles);
    assert!(
        particles[0].acceleration[0] > 0.0,
        "the left particle must move towards +x, a = {:?}",
        particles[0].acceleration
    );
    assert!(
        particles[1].acceleration[0] < 0.0,
        "the right particle must move towards -x, a = {:?}",
        particles[1].acceleration
    );
}

/// The pair acceleration is `mⱼ (σᵢ/ρᵢ² + σⱼ/ρⱼ²) ∂W/∂x` with the analytic
/// cubic-spline gradient (`σ = −10⁶ Pa`).
#[test]
fn test_stress_divergence_pair_acceleration_has_the_analytic_magnitude() {
    let mut particles = isotropic_pair(-1.0e6);
    stress_divergence(&mut particles);
    let (h, r) = (0.05, 0.04);
    let dw_dr = cubic_spline_dw_dr(r, h);
    // ∇ᵢW for xᵢ − xⱼ = (−0.04, 0, 0) is dW/dr · (−1, 0, 0).
    let grad_x = -dw_dr;
    let expected = 1e-3 * (2.0 * -1.0e6 / (2700.0 * 2700.0)) * grad_x;
    let got = particles[0].acceleration[0];
    assert!(
        (got - expected).abs() < 1e-9 * expected.abs(),
        "a_0x = {got}, expected {expected}"
    );
    assert!(expected < 0.0, "sanity: the expected value points to -x");
}

/// The acceleration is built from the Cauchy stress; the first
/// Piola-Kirchhoff tensor is not consulted.
#[test]
fn test_stress_divergence_ignores_the_first_piola_kirchhoff_stress() {
    let mut particles = isotropic_pair(0.0);
    for p in particles.iter_mut() {
        p.pk1_stress = [[1.0e9, 0.0, 0.0], [0.0, 1.0e9, 0.0], [0.0, 0.0, 1.0e9]];
    }
    stress_divergence(&mut particles);
    for p in &particles {
        assert_eq!(p.acceleration, [0.0; 3]);
    }
}

/// A linear stress field `σ_xx = k x` has the divergence `∂σ_xx/∂x = k`, i.e. a
/// uniform acceleration `k/ρ` of an interior particle: this checks sign and
/// magnitude of the whole sum at once.
#[test]
fn test_stress_divergence_of_a_linear_stress_field_is_k_over_rho() {
    let dx = 0.01;
    let mut particles = lattice(7, dx);
    let k = 1.0e8; // Pa/m
    for p in particles.iter_mut() {
        p.stress[0][0] = k * p.position[0];
    }
    stress_divergence(&mut particles);
    let centre = 3 * 49 + 3 * 7 + 3;
    let expected = k / 2700.0;
    let a = particles[centre].acceleration;
    assert!(
        (a[0] - expected).abs() < 0.03 * expected,
        "a_x = {}, expected {expected}",
        a[0]
    );
    assert!(a[1].abs() < 1e-6 * expected && a[2].abs() < 1e-6 * expected);
}

/// Pressure `σ = −p I` with a linear pressure ramp accelerates the interior
/// along the *decreasing* pressure, `a = −∇p/ρ`.
#[test]
fn test_stress_divergence_pressure_ramp_accelerates_down_the_gradient() {
    let dx = 0.01;
    let mut particles = lattice(7, dx);
    let slope = 5.0e7; // Pa/m, pressure grows with y
    for p in particles.iter_mut() {
        let pressure = slope * p.position[1];
        for d in 0..3 {
            p.stress[d][d] = -pressure;
        }
    }
    stress_divergence(&mut particles);
    let centre = 3 * 49 + 3 * 7 + 3;
    let expected = -slope / 2700.0;
    let a = particles[centre].acceleration;
    assert!(
        (a[1] - expected).abs() < 0.03 * expected.abs(),
        "a_y = {}, expected {expected}",
        a[1]
    );
    assert!(a[0].abs() < 1e-6 * expected.abs() && a[2].abs() < 1e-6 * expected.abs());
}

/// Random masses, densities and symmetric stresses: the pair forces cancel, so
/// the total momentum rate `Σ mᵢ aᵢ` vanishes to rounding.
#[test]
fn test_stress_divergence_conserves_total_momentum() {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let dx = 0.01;
    let mut particles = lattice(4, dx);
    for p in particles.iter_mut() {
        p.mass *= 0.5 + lcg_unit(&mut state);
        p.density *= 0.9 + 0.2 * lcg_unit(&mut state);
        let mut entries = [0.0f64; 6];
        for entry in entries.iter_mut() {
            *entry = 2.0e6 * (lcg_unit(&mut state) - 0.5);
        }
        let [xx, yy, zz, xy, xz, yz] = entries;
        p.stress = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
    }
    stress_divergence(&mut particles);
    let mut total = [0.0f64; 3];
    let mut scale = 0.0f64;
    for p in &particles {
        for (sum, &a) in total.iter_mut().zip(p.acceleration.iter()) {
            *sum += p.mass * a;
            scale += (p.mass * a).abs();
        }
    }
    assert!(scale > 1.0e-3, "the forces must be non-trivial: {scale}");
    for (d, rate) in total.iter().enumerate() {
        assert!(
            rate.abs() < 1e-12 * scale,
            "momentum rate {rate} in direction {d} (scale {scale})"
        );
    }
}

// ---------------------------------------------------------------------------
// hourglass_control
// ---------------------------------------------------------------------------

/// Two particles approaching each other at ±1 m/s: after the hourglass
/// acceleration and one explicit Euler step their relative speed must shrink
/// by exactly `2 |a| dt` (the pair is damped, not accelerated).
#[test]
fn test_hourglass_control_damps_the_relative_velocity_of_a_pair() {
    let (alpha_hg, sound_speed, dt) = (0.05, 5000.0, 1.0e-6);
    let mut particles = vec![
        ElasticParticle::new([0.0, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
        ElasticParticle::new([0.04, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
    ];
    particles[0].velocity = [1.0, 0.0, 0.0];
    particles[1].velocity = [-1.0, 0.0, 0.0];
    let before = (particles[0].velocity[0] - particles[1].velocity[0]).abs();
    hourglass_control(&mut particles, alpha_hg, sound_speed);
    assert!(
        particles[0].acceleration[0] < 0.0,
        "the faster particle must be decelerated: {:?}",
        particles[0].acceleration
    );
    assert!(
        particles[1].acceleration[0] > 0.0,
        "the slower particle must be accelerated: {:?}",
        particles[1].acceleration
    );
    for p in particles.iter_mut() {
        for d in 0..3 {
            p.velocity[d] += p.acceleration[d] * dt;
        }
    }
    let after = (particles[0].velocity[0] - particles[1].velocity[0]).abs();
    assert!(after < before, "relative speed grew: {before} -> {after}");
    // a_0 = (α c / h) (m/ρ) W(r, h) (v_1 − v_0): the rate α c / h has the unit 1/s.
    let h = 0.05;
    let rate = alpha_hg * sound_speed / h * (1e-3 / 2700.0) * cubic_kernel(0.04, h);
    let expected_drop = 2.0 * (2.0 * rate) * dt;
    assert!(
        ((before - after) - expected_drop).abs() < 1e-9 * expected_drop,
        "relative speed dropped by {}, expected {expected_drop}",
        before - after
    );
    // The rate is a few 1/s, i.e. an acceleration of a few m/s² for 1 m/s
    // differences: a wrong unit (no 1/h) would give 5e-2 times that.
    assert!(
        (particles[0].acceleration[0] + 4.0).abs() < 0.1,
        "a_0x = {}",
        particles[0].acceleration[0]
    );
}

/// The kernel weight of the pair acceleration is dimensionless, so the
/// acceleration is linear in both the coefficient `α` and the characteristic
/// speed `c` (the prefactor `α c / h` is a rate).
#[test]
fn test_hourglass_control_scales_linearly_with_alpha_and_sound_speed() {
    let make = || {
        let mut particles = vec![
            ElasticParticle::new([0.0, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
            ElasticParticle::new([0.04, 0.0, 0.0], 1e-3, 2700.0, 0.05, 70e9, 0.33),
        ];
        particles[0].velocity = [2.0, 1.0, 0.5];
        particles
    };
    let mut base = make();
    hourglass_control(&mut base, 0.05, 5000.0);
    let mut double_alpha = make();
    hourglass_control(&mut double_alpha, 0.10, 5000.0);
    let mut double_speed = make();
    hourglass_control(&mut double_speed, 0.05, 10000.0);
    for d in 0..3 {
        let a = base[0].acceleration[d];
        assert!(a.abs() > 1e-3);
        assert!((double_alpha[0].acceleration[d] - 2.0 * a).abs() < 1e-12 * a.abs());
        assert!((double_speed[0].acceleration[d] - 2.0 * a).abs() < 1e-12 * a.abs());
    }
}

/// With unequal masses and densities the correction still conserves total
/// momentum, and it only ever removes kinetic energy.
#[test]
fn test_hourglass_control_conserves_momentum_and_never_adds_energy() {
    let mut state = 0xA076_1D64_78BD_642F_u64;
    let dx = 0.01;
    let mut particles = lattice(4, dx);
    for p in particles.iter_mut() {
        p.mass *= 0.5 + lcg_unit(&mut state);
        p.density *= 0.7 + 0.6 * lcg_unit(&mut state);
        for d in 0..3 {
            p.velocity[d] = 4.0 * (lcg_unit(&mut state) - 0.5);
        }
    }
    hourglass_control(&mut particles, 0.1, 5000.0);
    let mut momentum_rate = [0.0f64; 3];
    let mut scale = 0.0f64;
    let mut power = 0.0f64;
    for p in &particles {
        for (d, &a) in p.acceleration.iter().enumerate() {
            momentum_rate[d] += p.mass * a;
            scale += (p.mass * a).abs();
            power += p.mass * p.velocity[d] * a;
        }
    }
    assert!(scale > 1.0, "the correction must be non-trivial: {scale}");
    for (d, rate) in momentum_rate.iter().enumerate() {
        assert!(
            rate.abs() < 1e-12 * scale,
            "momentum rate {rate} in direction {d} (scale {scale})"
        );
    }
    assert!(power < 0.0, "kinetic energy must decrease, power = {power}");
}

/// A checkerboard velocity pattern (the classic zero-energy mode) is damped
/// at every particle: `aᵢ·vᵢ ≤ 0`, strictly for the interior ones.
#[test]
fn test_hourglass_control_damps_a_checkerboard_mode_everywhere() {
    let dx = 0.01;
    let n = 5;
    let mut particles = lattice(n, dx);
    for (idx, p) in particles.iter_mut().enumerate() {
        let (ix, iy, iz) = (idx / (n * n), (idx / n) % n, idx % n);
        let sign = if (ix + iy + iz) % 2 == 0 { 1.0 } else { -1.0 };
        p.velocity = [sign, 0.0, 0.0];
    }
    hourglass_control(&mut particles, 0.1, 5000.0);
    for (idx, p) in particles.iter().enumerate() {
        let along = p.acceleration[0] * p.velocity[0];
        assert!(along <= 0.0, "particle {idx}: a·v = {along}");
    }
    let centre = &particles[2 * 25 + 2 * 5 + 2];
    assert!(centre.acceleration[0] * centre.velocity[0] < -1.0);
}

// ---------------------------------------------------------------------------
// Jaumann stress rate
// ---------------------------------------------------------------------------

/// Velocity gradient `L_ab = ∂v_a/∂x_b` of a rigid rotation with angular
/// velocity `w` about `z` (`v = w ẑ × x = (−w y, w x, 0)`).
fn rotation_velocity_gradient(w: f64) -> [[f64; 3]; 3] {
    [[0.0, -w, 0.0], [w, 0.0, 0.0], [0.0, 0.0, 0.0]]
}

/// A rigid rotation about `z` turns `σ = diag(s, 0, 0)` into
/// `R σ R^T` whose off-diagonal entry grows as `+s w dt`.
#[test]
fn test_jaumann_increment_of_a_rigid_rotation_has_the_sign_of_r_sigma_rt() {
    let model = SphElasticStress::from_elastic_params(70e9, 0.3);
    let (s, w, dt) = (1.0e6, 3.0, 1.0e-4);
    let sigma = [[s, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let spin = rotation_velocity_gradient(w); // L is purely skew, Ω = L
    let ds = model.jaumann_stress_increment(sigma, [[0.0; 3]; 3], spin, dt);
    assert!(
        (ds[0][1] - s * w * dt).abs() < 1e-12 * s * w * dt,
        "ds_xy = {}, expected {}",
        ds[0][1],
        s * w * dt
    );
    assert!((ds[1][0] - s * w * dt).abs() < 1e-12 * s * w * dt);
    assert_eq!(ds[0][0], 0.0);
    assert_eq!(ds[1][1], 0.0);
    // Cross-check with the exact rotation R(w dt) σ R(w dt)^T to first order.
    let angle = w * dt;
    let exact_xy = s * angle.cos() * angle.sin();
    assert!((ds[0][1] - exact_xy).abs() < 1e-3 * exact_xy);
}

/// The same through the particle update: the stress of a body that rotates
/// rigidly through 0.5 rad follows `R σ₀ Rᵀ`.
#[test]
fn test_update_particle_stress_follows_a_finite_rigid_rotation() {
    let model = SphElasticStress::from_elastic_params(70e9, 0.3);
    let mut particle = ElasticSphParticle::new([0.0; 3], 1e-3, 2700.0, 0.05, 70e9, 0.3);
    let s0 = 1.0e6;
    particle.stress[0][0] = s0;
    let w = 10.0;
    let steps = 5000;
    let dt = 0.05 / steps as f64;
    for _ in 0..steps {
        model.update_particle_stress(&mut particle, rotation_velocity_gradient(w), dt);
    }
    let theta = w * dt * steps as f64;
    let (c, s) = (theta.cos(), theta.sin());
    let expected = [
        [s0 * c * c, s0 * c * s, 0.0],
        [s0 * c * s, s0 * s * s, 0.0],
        [0.0, 0.0, 0.0],
    ];
    for (a, (stress_row, expected_row)) in particle.stress.iter().zip(expected.iter()).enumerate() {
        for (b, (value, want)) in stress_row.iter().zip(expected_row.iter()).enumerate() {
            assert!(
                (value - want).abs() < 2.0e-3 * s0,
                "stress[{a}][{b}] = {value}, expected {want}"
            );
        }
    }
    // The off-diagonal component has the sign of the rotation, not its negative.
    assert!(particle.stress[0][1] > 0.4 * s0);
}

/// Pure rotation must not change the invariants of the stress (the mean
/// pressure and the von Mises equivalent stress).
#[test]
fn test_update_particle_stress_rotation_preserves_the_invariants() {
    let model = SphElasticStress::from_elastic_params(70e9, 0.3);
    let mut particle = ElasticSphParticle::new([0.0; 3], 1e-3, 2700.0, 0.05, 70e9, 0.3);
    particle.stress = [[2.0e6, 0.5e6, 0.0], [0.5e6, -1.0e6, 0.0], [0.0, 0.0, 0.3e6]];
    let trace_before: f64 = (0..3).map(|i| particle.stress[i][i]).sum();
    let steps = 2000;
    for _ in 0..steps {
        model.update_particle_stress(&mut particle, rotation_velocity_gradient(20.0), 2.5e-5);
    }
    let trace_after: f64 = (0..3).map(|i| particle.stress[i][i]).sum();
    assert!((trace_after - trace_before).abs() < 1e-9 * trace_before.abs());
    // Frobenius norm is invariant under rotation; explicit Euler inflates it by
    // ~ θ dθ / 2 (here 1 rad · 5e-4 rad / 2).
    let frobenius = |s: &[[f64; 3]; 3]| s.iter().flatten().map(|v| v * v).sum::<f64>().sqrt();
    let initial = frobenius(&[[2.0e6, 0.5e6, 0.0], [0.5e6, -1.0e6, 0.0], [0.0, 0.0, 0.3e6]]);
    assert!((frobenius(&particle.stress) - initial).abs() < 1.0e-3 * initial);
}

/// An isotropic stress `σ = −p I` is invariant under rotation, so the
/// co-rotation `Ω σ − σ Ω` must vanish for any spin.
#[test]
fn test_jaumann_increment_leaves_an_isotropic_stress_unchanged() {
    let model = SphElasticStress::from_elastic_params(70e9, 0.3);
    let mut sigma = IDENTITY3;
    for (d, row) in sigma.iter_mut().enumerate() {
        row[d] = -4.0e6;
    }
    let spin = [[0.0, -2.0, 1.0], [2.0, 0.0, -0.5], [-1.0, 0.5, 0.0]];
    let ds = model.jaumann_stress_increment(sigma, [[0.0; 3]; 3], spin, 1e-3);
    for row in &ds {
        for v in row {
            assert!(v.abs() < 1e-9, "isotropic stress changed under spin: {v}");
        }
    }
}
