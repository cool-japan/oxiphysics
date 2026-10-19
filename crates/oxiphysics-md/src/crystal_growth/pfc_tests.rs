//! Regression tests of the conserved, spectral phase-field-crystal dynamics.

use super::*;
use oxifft::rfft2d;
use oxiphysics_core::random::Xoshiro256;
use std::f64::consts::PI;

/// Field `phi_0 + amplitude * cos(2π (mx i / nx + my j / ny))` on an `nx × ny` grid.
fn single_mode(
    nx: usize,
    ny: usize,
    dx: f64,
    phi_0: f64,
    amplitude: f64,
    mode: (usize, usize),
) -> PhaseFieldCrystal {
    let mut pfc = PhaseFieldCrystal::new(nx, ny, dx, phi_0);
    for i in 0..nx {
        for j in 0..ny {
            let phase = 2.0
                * PI
                * (mode.0 as f64 * i as f64 / nx as f64 + mode.1 as f64 * j as f64 / ny as f64);
            pfc.phi[i * ny + j] = phi_0 + amplitude * phase.cos();
        }
    }
    pfc
}

/// Amplitude of the cosine mode `(mx, my)` of the field (`2/N Σ φ cos(...)`).
fn mode_amplitude(pfc: &PhaseFieldCrystal, mode: (usize, usize)) -> f64 {
    let (nx, ny) = (pfc.nx, pfc.ny);
    let mut sum = 0.0;
    for i in 0..nx {
        for j in 0..ny {
            let phase = 2.0
                * PI
                * (mode.0 as f64 * i as f64 / nx as f64 + mode.1 as f64 * j as f64 / ny as f64);
            sum += pfc.phi[i * ny + j] * phase.cos();
        }
    }
    2.0 * sum / (nx * ny) as f64
}

/// `phi_0` plus a reproducible uniform perturbation of the given amplitude.
fn perturbed(
    nx: usize,
    ny: usize,
    dx: f64,
    phi_0: f64,
    amplitude: f64,
    seed: u64,
) -> PhaseFieldCrystal {
    let mut rng = Xoshiro256::new(seed);
    let mut pfc = PhaseFieldCrystal::new(nx, ny, dx, phi_0);
    for value in pfc.phi.iter_mut() {
        *value = phi_0 + amplitude * (2.0 * rng.next_f64() - 1.0);
    }
    pfc
}

/// Squared wave number of the mode `(mx, my)` on the grid.
fn mode_k_squared(pfc: &PhaseFieldCrystal, mode: (usize, usize)) -> f64 {
    let kx = 2.0 * PI * mode.0 as f64 / (pfc.nx as f64 * pfc.dx);
    let ky = 2.0 * PI * mode.1 as f64 / (pfc.ny as f64 * pfc.dx);
    kx * kx + ky * ky
}

#[test]
fn test_pfc_defaults_and_builders() {
    let pfc = PhaseFieldCrystal::new(8, 6, 0.5, 0.1);
    assert_eq!(pfc.phi.len(), 48);
    assert_eq!(pfc.r, DEFAULT_PFC_R);
    assert_eq!(pfc.seed, DEFAULT_PFC_SEED);
    let tuned = pfc.with_r(-0.4).with_seed(99);
    assert_eq!((tuned.r, tuned.seed), (-0.4, 99));
}

#[test]
fn test_pfc_uniform_state_is_a_fixed_point() {
    let mut pfc = PhaseFieldCrystal::new(16, 12, PI / 4.0, 0.3);
    for _ in 0..5 {
        pfc.step(0.5, 1.0, 0.0);
    }
    for value in &pfc.phi {
        assert!((value - 0.3).abs() < 1e-14, "{value}");
    }
}

#[test]
fn test_pfc_linear_dispersion_matches_the_semi_implicit_factor() {
    // Tiny single modes (the cubic term is ~A³): one step multiplies the amplitude by
    // 1 / (1 + dt M k² (r + (1 - k²)²)).  dx = π/4 makes the mode (4, 0) of a 32-grid the
    // preferred wave number k = 1, where the symbol is just r = -0.25.
    let (dt, mobility, amplitude) = (0.5, 1.0, 1.0e-6);
    for (nx, ny, mode) in [
        (32, 32, (4, 0)),
        (32, 32, (8, 0)),
        (32, 32, (3, 5)),
        (32, 32, (0, 2)),
        (24, 40, (3, 0)),
        (24, 40, (0, 5)),
        (24, 40, (2, 7)),
    ] {
        let mut pfc = single_mode(nx, ny, PI / 4.0, 0.0, amplitude, mode);
        let k2 = mode_k_squared(&pfc, mode);
        let symbol = pfc.r + (1.0 - k2) * (1.0 - k2);
        let expected = amplitude / (1.0 + dt * mobility * k2 * symbol);
        pfc.step(dt, mobility, 0.0);
        let got = mode_amplitude(&pfc, mode);
        assert!(
            ((got - expected) / expected).abs() < 1e-9,
            "{nx}x{ny} mode {mode:?}: amplitude {got}, expected {expected} (k² = {k2})"
        );
    }
    // The k = 1 mode of the 32 x 32 grid grows by exactly 1/0.875 per step.
    let mut pfc = single_mode(32, 32, PI / 4.0, 0.0, amplitude, (4, 0));
    pfc.step(dt, mobility, 0.0);
    assert!((mode_amplitude(&pfc, (4, 0)) / amplitude - 1.0 / 0.875).abs() < 1e-9);
}

#[test]
fn test_pfc_mean_is_conserved_with_and_without_noise() {
    for noise in [0.0, 2.0e-3] {
        let mut pfc = perturbed(32, 24, PI / 4.0, 0.2, 0.05, 11).with_seed(5);
        let mean_before = pfc.mean_phi();
        for _ in 0..100 {
            pfc.step(0.2, 1.0, noise);
        }
        let drift = (pfc.mean_phi() - mean_before).abs();
        assert!(drift < 1e-12, "noise {noise}: mean drifted by {drift}");
        assert!(pfc.phi.iter().all(|v| v.is_finite()));
    }
}

#[test]
fn test_pfc_energy_of_a_uniform_state_and_a_single_mode_in_closed_form() {
    // F / (N dx²) = ½ (r + 1) φ0² + ¼ L(k) A² + ¼ (φ0⁴ + 3 φ0² A² + 3 A⁴ / 8) for
    // φ = φ0 + A cos(k·x), with L(k) = r + (1 - k²)².
    for (nx, ny, mode) in [
        (32, 32, (4, 0)),
        (24, 40, (3, 0)),
        (24, 40, (2, 7)),
        (15, 21, (2, 3)),
    ] {
        let (phi_0, amplitude) = (0.3, 0.2);
        let dx = PI / 4.0;
        let uniform = PhaseFieldCrystal::new(nx, ny, dx, phi_0);
        let cells = (nx * ny) as f64;
        let uniform_density = 0.5 * (uniform.r + 1.0) * phi_0 * phi_0 + 0.25 * phi_0.powi(4);
        assert!(
            (uniform.energy() - cells * dx * dx * uniform_density).abs()
                < 1e-12 * uniform.energy().abs(),
            "{nx}x{ny} uniform"
        );

        let pfc = single_mode(nx, ny, dx, phi_0, amplitude, mode);
        let k2 = mode_k_squared(&pfc, mode);
        let symbol = pfc.r + (1.0 - k2) * (1.0 - k2);
        let density = 0.5 * (pfc.r + 1.0) * phi_0 * phi_0
            + 0.25 * symbol * amplitude * amplitude
            + 0.25
                * (phi_0.powi(4)
                    + 3.0 * phi_0 * phi_0 * amplitude * amplitude
                    + 3.0 * amplitude.powi(4) / 8.0);
        let expected = cells * dx * dx * density;
        assert!(
            (pfc.energy() - expected).abs() < 1e-11 * expected.abs(),
            "{nx}x{ny} mode {mode:?}: energy {}, expected {expected}",
            pfc.energy()
        );
    }
}

#[test]
fn test_pfc_energy_is_non_increasing_without_noise() {
    let mut pfc = perturbed(32, 32, PI / 4.0, 0.28, 0.1, 2024);
    let initial = pfc.energy();
    let mut previous = initial;
    for step in 0..300 {
        pfc.step(0.1, 1.0, 0.0);
        let energy = pfc.energy();
        assert!(
            energy <= previous + 1e-9 * previous.abs(),
            "step {step}: energy rose from {previous} to {energy}"
        );
        previous = energy;
    }
    assert!(
        previous < initial - 0.05 * initial.abs(),
        "the energy must relax: {initial} -> {previous}"
    );
}

#[test]
fn test_pfc_dynamics_selects_the_preferred_wave_number() {
    // r = -0.25, mean 0.2: the band around k = 1 is linearly unstable and the saturated
    // pattern is periodic with that wave number.
    let dx = PI / 4.0;
    let mut pfc = perturbed(32, 32, dx, 0.2, 0.01, 7);
    for _ in 0..800 {
        pfc.step(0.25, 1.0, 0.0);
    }
    let amplitude = pfc.phi.iter().fold(0.0_f64, |m, v| m.max((v - 0.2).abs()));
    assert!(
        amplitude > 0.1,
        "a pattern must have formed, amplitude {amplitude}"
    );
    let spectrum = rfft2d(&pfc.phi, 32, 32);
    let ny_half = 17;
    let mut strongest = (0.0_f64, 0.0_f64);
    for i in 0..32 {
        for jh in 0..ny_half {
            if i == 0 && jh == 0 {
                continue;
            }
            let fx = if i <= 16 { i as f64 } else { i as f64 - 32.0 };
            let k = 2.0 * PI / (32.0 * dx) * (fx * fx + (jh * jh) as f64).sqrt();
            let power = spectrum[i * ny_half + jh].norm_sqr();
            if power > strongest.0 {
                strongest = (power, k);
            }
        }
    }
    assert!(
        (0.85..=1.15).contains(&strongest.1),
        "dominant wave number {} should be the preferred k = 1",
        strongest.1
    );
}

#[test]
fn test_pfc_noise_has_the_fluctuation_dissipation_amplitude() {
    // From φ = 0 one tiny step gives φ̂ = Δt (i k·ζ̂) / (1 + Δt M k² L), so
    // E|φ̂|² (1 + Δt M k² L)² / k² = 2 M T Δt N / dx² for white noise of variance
    // 2 M T / (dx² Δt).
    let (nx, ny, dx) = (64, 64, PI / 4.0);
    let (dt, mobility, temperature) = (1.0e-3, 1.0, 0.01);
    let mut pfc = PhaseFieldCrystal::new(nx, ny, dx, 0.0).with_seed(31337);
    pfc.step(dt, mobility, temperature);
    let spectrum = rfft2d(&pfc.phi, nx, ny);
    let ny_half = ny / 2 + 1;
    let mut sum = 0.0;
    let mut bins = 0;
    for i in 1..nx {
        if i == nx / 2 {
            continue;
        }
        for jh in 1..ny / 2 {
            let fx = if i <= nx / 2 {
                i as f64
            } else {
                i as f64 - nx as f64
            };
            let kx = 2.0 * PI * fx / (nx as f64 * dx);
            let ky = 2.0 * PI * jh as f64 / (ny as f64 * dx);
            let k2 = kx * kx + ky * ky;
            let denominator = 1.0 + dt * mobility * k2 * (pfc.r + (1.0 - k2) * (1.0 - k2));
            sum += spectrum[i * ny_half + jh].norm_sqr() * denominator * denominator / k2;
            bins += 1;
        }
    }
    let measured = sum / bins as f64;
    let expected = 2.0 * mobility * temperature * dt * (nx * ny) as f64 / (dx * dx);
    assert!(
        ((measured - expected) / expected).abs() < 0.1,
        "noise power {measured}, fluctuation-dissipation value {expected}"
    );
    assert!(pfc.mean_phi().abs() < 1e-14, "conserved noise has no mean");
}

/// Correlation coefficient of two fields.
fn correlation(a: &[f64], b: &[f64]) -> f64 {
    let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    let norm_a: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let norm_b: f64 = b.iter().map(|y| y * y).sum::<f64>().sqrt();
    dot / (norm_a * norm_b)
}

#[test]
fn test_pfc_noise_is_fresh_every_step_and_follows_the_seed() {
    // One tiny step from φ = 0 is dominated by the noise term.
    let single_step = |seed: u64| {
        let mut pfc = PhaseFieldCrystal::new(32, 32, PI / 4.0, 0.0).with_seed(seed);
        pfc.step(1.0e-3, 1.0, 0.01);
        (pfc.phi, pfc.seed)
    };
    let (first, next_seed) = single_step(1234);
    let (repeated, repeated_seed) = single_step(1234);
    assert_eq!(first, repeated, "the same seed reproduces the noise");
    assert_eq!(next_seed, repeated_seed);
    assert!(first.iter().any(|v| v.abs() > 1e-6), "the noise must act");
    assert_ne!(next_seed, 1234, "a noisy step advances the seed");

    // The seed a step leaves behind drives the next step: its noise is an independent
    // realisation.  The old per-cell pattern (idx · φ).fract() was identical in every
    // step (correlation 1).
    let (second, _) = single_step(next_seed);
    assert!(
        correlation(&first, &second).abs() < 0.1,
        "successive noise realisations are correlated: {}",
        correlation(&first, &second)
    );
    let (other_seed, _) = single_step(4321);
    assert!(correlation(&first, &other_seed).abs() < 0.1);

    // A multi-step run is reproducible from its initial seed and differs for another seed.
    let run = |seed: u64| {
        let mut pfc = PhaseFieldCrystal::new(32, 32, PI / 4.0, 0.0).with_seed(seed);
        for _ in 0..3 {
            pfc.step(1.0e-3, 1.0, 0.01);
        }
        pfc
    };
    assert_eq!(run(9).phi, run(9).phi);
    assert_ne!(run(9).phi, run(10).phi);

    // Without noise the seed is left alone.
    let mut quiet = PhaseFieldCrystal::new(8, 8, 1.0, 0.1).with_seed(77);
    quiet.step(0.1, 1.0, 0.0);
    assert_eq!(quiet.seed, 77);
}

#[test]
fn test_pfc_step_ignores_invalid_arguments_and_degenerate_grids() {
    let mut pfc = perturbed(8, 8, 1.0, 0.1, 0.05, 3);
    let before = pfc.phi.clone();
    for (dt, mobility) in [
        (0.0, 1.0),
        (-0.1, 1.0),
        (f64::NAN, 1.0),
        (0.1, -1.0),
        (0.1, f64::NAN),
        (f64::INFINITY, 1.0),
    ] {
        pfc.step(dt, mobility, 0.0);
        assert_eq!(pfc.phi, before, "dt = {dt}, M = {mobility}");
    }
    let mut empty = PhaseFieldCrystal::new(0, 5, 1.0, 0.0);
    empty.step(0.1, 1.0, 0.0);
    assert_eq!(empty.energy(), 0.0);
    assert_eq!(empty.mean_phi(), 0.0);
    let mut inconsistent = PhaseFieldCrystal::new(4, 4, 1.0, 0.0);
    inconsistent.phi.truncate(5);
    inconsistent.step(0.1, 1.0, 0.0);
    assert_eq!(inconsistent.phi.len(), 5);
    assert_eq!(inconsistent.energy(), 0.0);
    // Zero mobility freezes the field, even with noise requested.
    let mut frozen = perturbed(8, 8, 1.0, 0.1, 0.05, 3);
    frozen.step(0.1, 0.0, 0.5);
    assert_eq!(frozen.phi, before);
}

#[test]
fn test_pfc_odd_and_rectangular_grids_conserve_the_mean() {
    for (nx, ny) in [(15, 21), (12, 9), (9, 12)] {
        let mut pfc = perturbed(nx, ny, 0.8, 0.15, 0.05, 99).with_seed(8);
        let mean_before = pfc.mean_phi();
        for _ in 0..50 {
            pfc.step(0.2, 1.0, 1.0e-3);
        }
        assert!(
            (pfc.mean_phi() - mean_before).abs() < 1e-12,
            "{nx}x{ny}: mean drifted"
        );
        assert!(pfc.phi.iter().all(|v| v.is_finite()));
    }
}
