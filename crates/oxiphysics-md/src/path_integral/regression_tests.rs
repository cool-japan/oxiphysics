//! Regression tests of the centroid-virial estimators and of the adiabatic
//! centroid molecular dynamics (CMD) propagator against exact results for a
//! three-dimensional isotropic harmonic oscillator.
//!
//! For `V = ½ m ω² r²` the ring polymer with `P` beads is Gaussian in the
//! normal-mode coordinates `q_k` (orthogonal transform `C`), with the exact
//! path-integral variance (per Cartesian component)
//!
//! ```text
//! ⟨q_k²⟩ = P / (β m (ω_k² + ω²)),   ω_k = 2 ω_P sin(kπ/P),   ω_P = P k_B T / ħ,
//! ```
//!
//! so configurations can be drawn exactly and every estimator has a closed-form
//! expectation, independent of any dynamics.

use super::functions::{
    HBAR, KB, beads_to_normal_modes, normal_mode_matrix, normal_modes_to_beads,
    virial_kinetic_estimator,
};
use super::types::{CentroidMD, PileLThermostat, PimdAtom, PimdState, RingPolymer};
use crate::gaussian_noise::standard_normal;
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::f64::consts::PI;

/// Proton mass (kg).
const MASS: f64 = 1.67e-27;
/// Harmonic well frequency (rad s⁻¹); `βħω ≈ 2.5` at 300 K, a strongly quantum regime.
const OMEGA: f64 = 1.0e14;
/// Temperature (K).
const TEMPERATURE: f64 = 300.0;

fn beta(temperature: f64) -> f64 {
    1.0 / (KB * temperature)
}

/// Free ring-polymer frequency `ω_k = 2 ω_P sin(kπ/P)` of normal mode `k`.
fn mode_frequency(k: usize, p: usize, temperature: f64) -> f64 {
    let omega_p = p as f64 * KB * temperature / HBAR;
    2.0 * omega_p * (k as f64 * PI / p as f64).sin()
}

/// Exact variance `⟨q_k²⟩` of one Cartesian component of mode `k`.
fn mode_variance(k: usize, p: usize, mass: f64, omega: f64, temperature: f64) -> f64 {
    let w_k = mode_frequency(k, p, temperature);
    p as f64 / (beta(temperature) * mass * (w_k * w_k + omega * omega))
}

/// Draw the bead positions of a harmonic ring polymer from its exact distribution.
fn sample_ring(p: usize, c: &[Vec<f64>], rng: &mut StdRng) -> Vec<[f64; 3]> {
    let modes: Vec<[f64; 3]> = (0..p)
        .map(|k| {
            let sigma = mode_variance(k, p, MASS, OMEGA, TEMPERATURE).sqrt();
            [
                sigma * standard_normal(rng),
                sigma * standard_normal(rng),
                sigma * standard_normal(rng),
            ]
        })
        .collect();
    normal_modes_to_beads(&modes, c)
}

/// Physical force `-m ω² r` on every bead of a harmonic well.
fn harmonic_forces(positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
    let k = MASS * OMEGA * OMEGA;
    positions
        .iter()
        .map(|r| [-k * r[0], -k * r[1], -k * r[2]])
        .collect()
}

/// Exact expectation of the centroid-virial kinetic energy of a `P`-bead ring
/// polymer in the harmonic well: `(3/2β) [1 + Σ_{k=1}^{P-1} ω² / (ω_k² + ω²)]`.
fn expected_virial_kinetic_energy(p: usize) -> f64 {
    let sum: f64 = (1..p)
        .map(|k| {
            let w_k = mode_frequency(k, p, TEMPERATURE);
            OMEGA * OMEGA / (w_k * w_k + OMEGA * OMEGA)
        })
        .sum();
    1.5 / beta(TEMPERATURE) * (1.0 + sum)
}

/// Exact quantum kinetic energy `3 (ħω/4) coth(βħω/2)` of the oscillator.
fn exact_quantum_kinetic_energy() -> f64 {
    0.75 * HBAR * OMEGA / (0.5 * beta(TEMPERATURE) * HBAR * OMEGA).tanh()
}

// ---------------------------------------------------------------------------
// Centroid-virial kinetic-energy estimators
// ---------------------------------------------------------------------------

/// A ring of four beads on a circle of radius 10 pm around the origin.
fn spread_ring_positions() -> Vec<[f64; 3]> {
    vec![
        [1.0e-11, 0.0, 0.0],
        [0.0, 1.0e-11, 0.0],
        [-1.0e-11, 0.0, 0.0],
        [0.0, -1.0e-11, 0.0],
    ]
}

/// Evaluate the three estimator entry points on one configuration.
fn estimator_triplet(positions: &[[f64; 3]], forces: &[[f64; 3]]) -> [f64; 3] {
    let p = positions.len();
    let mut ring = RingPolymer::new(p, MASS, [0.0; 3]);
    ring.positions = positions.to_vec();
    let mut atom = PimdAtom::new(p, MASS);
    atom.replicas = positions.to_vec();
    atom.forces = forces.to_vec();
    let state = PimdState {
        atoms: vec![atom],
        temperature: TEMPERATURE,
        n_beads: p,
    };
    [
        ring.virial_kinetic_energy(TEMPERATURE, forces),
        state.quantum_kinetic_energy(),
        virial_kinetic_estimator(&[positions.to_vec()], &[forces.to_vec()], TEMPERATURE, p),
    ]
}

#[test]
fn virial_estimators_add_the_ring_spread_for_a_harmonic_well() {
    // F = -k r  =>  KE_vir = 3 k_B T / 2 + (k / 2P) Σ |r_i - r̄|²  >  classical.
    let positions = spread_ring_positions();
    let forces = harmonic_forces(&positions);
    let p = positions.len() as f64;
    let spread_sq: f64 = positions
        .iter()
        .map(|r| r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
        .sum();
    let classical = 1.5 * KB * TEMPERATURE;
    let expected = classical + MASS * OMEGA * OMEGA / (2.0 * p) * spread_sq;
    assert!(
        expected > 1.1 * classical,
        "the test geometry must matter: {expected} vs {classical}"
    );
    for (idx, ke) in estimator_triplet(&positions, &forces).iter().enumerate() {
        assert!(
            (ke - expected).abs() < 1e-12 * expected,
            "estimator {idx}: {ke} vs {expected}"
        );
        assert!(
            *ke > classical,
            "estimator {idx} must exceed the classical value"
        );
    }
}

#[test]
fn virial_estimators_are_odd_in_the_force_and_vanish_for_a_collapsed_ring() {
    let positions = spread_ring_positions();
    let forces = harmonic_forces(&positions);
    let flipped: Vec<[f64; 3]> = forces.iter().map(|f| [-f[0], -f[1], -f[2]]).collect();
    let classical = 1.5 * KB * TEMPERATURE;
    let plus = estimator_triplet(&positions, &forces);
    let minus = estimator_triplet(&positions, &flipped);
    for (a, b) in plus.iter().zip(minus.iter()) {
        assert!(
            ((a - classical) + (b - classical)).abs() < 1e-12 * classical,
            "correction must flip sign with the force: {a} {b}"
        );
    }
    // All beads at one point: every (r_i - r̄) vanishes, whatever the force.
    let collapsed = vec![[3.0e-11, -2.0e-11, 5.0e-12]; 4];
    let strong_forces = vec![[1.0e-9, -2.0e-9, 3.0e-9]; 4];
    for ke in estimator_triplet(&collapsed, &strong_forces) {
        assert!((ke - classical).abs() < 1e-12 * classical);
    }
}

#[test]
fn finite_p_virial_expectation_converges_to_the_exact_quantum_energy() {
    let exact = exact_quantum_kinetic_energy();
    let classical = 1.5 * KB * TEMPERATURE;
    assert!(exact > 1.4 * classical, "{exact} vs {classical}");
    let mut previous_error = f64::INFINITY;
    for p in [8_usize, 16, 32, 64, 128] {
        let error = (expected_virial_kinetic_energy(p) - exact).abs() / exact;
        assert!(
            error < previous_error,
            "P = {p}: error {error} did not shrink"
        );
        previous_error = error;
    }
    assert!(previous_error < 1e-4, "P = 128 error {previous_error}");
}

#[test]
fn sampled_virial_estimators_reproduce_the_harmonic_quantum_kinetic_energy() {
    // The estimators averaged over exact path-integral samples must give the
    // finite-P expectation, which is far above the classical 3 k_B T / 2 (the
    // sign-flipped estimator would land below it).
    let p = 16;
    let n_samples = 4000;
    let c = normal_mode_matrix(p);
    let mut rng = StdRng::seed_from_u64(2024);
    let mut sums = [0.0_f64; 3];
    for _ in 0..n_samples {
        let positions = sample_ring(p, &c, &mut rng);
        let forces = harmonic_forces(&positions);
        for (sum, ke) in sums.iter_mut().zip(estimator_triplet(&positions, &forces)) {
            *sum += ke;
        }
    }
    let expected = expected_virial_kinetic_energy(p);
    let exact = exact_quantum_kinetic_energy();
    let classical = 1.5 * KB * TEMPERATURE;
    for (idx, sum) in sums.iter().enumerate() {
        let mean = sum / n_samples as f64;
        assert!(
            (mean - expected).abs() / expected < 0.02,
            "estimator {idx}: <KE> = {mean}, finite-P expectation {expected}"
        );
        assert!(
            (mean - exact).abs() / exact < 0.03,
            "estimator {idx}: <KE> = {mean}, exact quantum value {exact}"
        );
        assert!(
            mean > 1.3 * classical,
            "estimator {idx}: {mean} vs {classical}"
        );
    }
    // The three entry points evaluate the same estimator.
    assert!((sums[0] - sums[1]).abs() < 1e-9 * sums[0]);
    assert!((sums[0] - sums[2]).abs() < 1e-9 * sums[0]);
}

// ---------------------------------------------------------------------------
// Bath helper used by CMD
// ---------------------------------------------------------------------------

#[test]
fn bath_temperature_override_sets_the_stationary_velocity_variance() {
    // A free particle with friction gamma0 = 1/dt relaxes to m <v²> = 3 k_B T_bath,
    // where T_bath is the explicit bath temperature, not the construction one.
    let (dt, bath) = (1.0e-15, 900.0);
    let th = PileLThermostat::new(1, TEMPERATURE, 1.0 / dt).expect("valid thermostat");
    let mut rng = StdRng::seed_from_u64(8);
    let mut mode_vel = [[0.0_f64; 3]];
    let n_steps = 100_000;
    let mut sum_v2 = 0.0;
    for _ in 0..n_steps {
        th.apply_at_bath_temperature(&mut mode_vel, MASS, bath, dt, &mut rng);
        sum_v2 += mode_vel[0].iter().map(|v| v * v).sum::<f64>();
    }
    let measured = MASS * sum_v2 / (3.0 * n_steps as f64 * KB);
    assert!(
        (measured - bath).abs() / bath < 0.02,
        "measured {measured} K, expected {bath} K"
    );
    // `apply` is the special case bath = construction temperature.
    let mut a = [[1.0e3_f64; 3]; 1];
    let mut b = a;
    let mut rng_a = StdRng::seed_from_u64(3);
    let mut rng_b = StdRng::seed_from_u64(3);
    th.apply(&mut a, MASS, dt, &mut rng_a);
    th.apply_at_bath_temperature(&mut b, MASS, TEMPERATURE, dt, &mut rng_b);
    assert_eq!(a, b);
}

// ---------------------------------------------------------------------------
// Adiabatic centroid MD
// ---------------------------------------------------------------------------

#[test]
fn cmd_construction_validates_its_parameters() {
    assert!(CentroidMD::new(0, 300.0, 0.5).is_err());
    assert!(CentroidMD::new(8, 0.0, 0.5).is_err());
    assert!(CentroidMD::new(8, -300.0, 0.5).is_err());
    assert!(CentroidMD::new(8, f64::NAN, 0.5).is_err());
    assert!(CentroidMD::new(8, 300.0, 0.0).is_err());
    assert!(CentroidMD::new(8, 300.0, 1.5).is_err());
    assert!(CentroidMD::new(8, 300.0, f64::NAN).is_err());
    assert!(CentroidMD::new(1, 300.0, 1.0).is_ok());
    let cmd = CentroidMD::new(8, 300.0, 0.25).expect("valid parameters");
    assert!((cmd.bath_temperature() - 2400.0).abs() < 1e-9);
    assert!((cmd.internal_mass(MASS) - MASS / 16.0).abs() < 1e-40);
    // Internal modes oscillate and relax 1/gamma times faster than the free ring.
    let free = CentroidMD::new(8, 300.0, 1.0).expect("valid parameters");
    assert_eq!(cmd.propagator.omega[0], 0.0);
    assert_eq!(cmd.thermostat.gamma[0], 0.0);
    for k in 1..8 {
        let ratio = cmd.propagator.omega[k] / free.propagator.omega[k];
        assert!((ratio - 4.0).abs() < 1e-12, "mode {k}: {ratio}");
        assert!((cmd.thermostat.gamma[k] - 2.0 * cmd.propagator.omega[k]).abs() < 1e-3);
    }
}

#[test]
fn cmd_step_rejects_inconsistent_inputs_and_leaves_the_state_untouched() {
    let cmd = CentroidMD::new(4, TEMPERATURE, 0.5).expect("valid parameters");
    let mut rng = StdRng::seed_from_u64(1);
    let mut ring = RingPolymer::new(4, MASS, [1.0e-11, 0.0, 0.0]);
    ring.velocities[2] = [10.0, 0.0, 0.0];
    let mut forces = vec![[0.0; 3]; 4];
    let before = (ring.clone(), forces.clone());
    let zero = |pos: &[[f64; 3]]| vec![[0.0; 3]; pos.len()];

    // Wrong number of forces handed in.
    let mut short = vec![[0.0; 3]; 3];
    assert!(
        cmd.step(&mut ring, &mut short, zero, 1e-16, &mut rng)
            .is_err()
    );
    // Closure returning the wrong number of forces.
    let bad_closure = |_: &[[f64; 3]]| vec![[0.0; 3]; 2];
    assert!(
        cmd.step(&mut ring, &mut forces, bad_closure, 1e-16, &mut rng)
            .is_err()
    );
    // Non-positive or non-finite time step.
    assert!(
        cmd.step(&mut ring, &mut forces, zero, 0.0, &mut rng)
            .is_err()
    );
    assert!(
        cmd.step(&mut ring, &mut forces, zero, f64::NAN, &mut rng)
            .is_err()
    );
    // Ring with a different bead count.
    let mut other = RingPolymer::new(3, MASS, [0.0; 3]);
    assert!(
        cmd.step(&mut other, &mut forces, zero, 1e-16, &mut rng)
            .is_err()
    );
    // Non-positive mass.
    let mut massless = RingPolymer::new(4, 0.0, [0.0; 3]);
    assert!(
        cmd.step(&mut massless, &mut forces, zero, 1e-16, &mut rng)
            .is_err()
    );

    assert_eq!(ring.positions, before.0.positions);
    assert_eq!(ring.velocities, before.0.velocities);
    assert_eq!(forces, before.1);
}

#[test]
fn cmd_keeps_internal_mode_momenta_instead_of_copying_the_centroid_velocity() {
    // Regression: the old step overwrote every bead velocity with the centroid
    // velocity.  Without the bath and without forces the free ring must conserve
    // the internal-mode energy Σ ½ m̃ v_k² + ½ m ω_k² q_k² exactly, keep the
    // centroid velocity constant and leave the beads with distinct velocities.
    let (p, gamma) = (6_usize, 0.5);
    let mut cmd = CentroidMD::new(p, TEMPERATURE, gamma).expect("valid parameters");
    for g in cmd.thermostat.gamma.iter_mut() {
        *g = 0.0;
    }
    let mut ring = RingPolymer::new(p, MASS, [0.0; 3]);
    for (i, (r, v)) in ring
        .positions
        .iter_mut()
        .zip(ring.velocities.iter_mut())
        .enumerate()
    {
        let x = i as f64 + 1.0;
        *r = [1.0e-12 * x, -2.0e-12 * x * x, 5.0e-13 * (3.0 - x)];
        *v = [200.0 * x, 150.0 * (x - 2.5), -80.0 * x * x];
    }
    let c = cmd.propagator.c_matrix.clone();
    let physical_omega: Vec<f64> = cmd.propagator.omega.iter().map(|w| w * gamma).collect();
    let internal_energy = |ring: &RingPolymer| -> f64 {
        let q = beads_to_normal_modes(&ring.positions, &c);
        let v = beads_to_normal_modes(&ring.velocities, &c);
        (1..p)
            .map(|k| {
                let kinetic: f64 = v[k].iter().map(|x| x * x).sum::<f64>();
                let potential: f64 = q[k].iter().map(|x| x * x).sum::<f64>();
                0.5 * cmd.internal_mass(MASS) * kinetic
                    + 0.5 * MASS * physical_omega[k] * physical_omega[k] * potential
            })
            .sum()
    };
    let e0 = internal_energy(&ring);
    let v_centroid = ring.centroid_velocity();
    let c0 = ring.centroid();
    let mut forces = vec![[0.0; 3]; p];
    let mut rng = StdRng::seed_from_u64(4);
    let dt = 2.0e-17;
    let steps = 300;
    for _ in 0..steps {
        cmd.step(
            &mut ring,
            &mut forces,
            |pos| vec![[0.0; 3]; pos.len()],
            dt,
            &mut rng,
        )
        .expect("consistent inputs");
    }
    let e1 = internal_energy(&ring);
    assert!(
        (e1 - e0).abs() < 1e-9 * e0,
        "internal energy drifted: {e0} -> {e1}"
    );
    let v1 = ring.centroid_velocity();
    for d in 0..3 {
        assert!((v1[d] - v_centroid[d]).abs() < 1e-9 * v_centroid[d].abs().max(1.0));
    }
    let c1 = ring.centroid();
    for d in 0..3 {
        let expected = c0[d] + v_centroid[d] * dt * steps as f64;
        assert!(
            (c1[d] - expected).abs() < 1e-9 * expected.abs().max(1.0e-12),
            "centroid {d}: {} vs {expected}",
            c1[d]
        );
    }
    let distinct = ring
        .velocities
        .windows(2)
        .any(|w| (w[0][0] - w[1][0]).abs() > 1.0);
    assert!(distinct, "internal-mode velocities must survive the step");
}

#[test]
fn cmd_never_thermostats_the_centroid_mode() {
    // Even a friction stored for mode 0 must not reach the centroid velocity.
    let mut cmd = CentroidMD::new(4, TEMPERATURE, 0.5).expect("valid parameters");
    cmd.thermostat.gamma[0] = 1.0e15;
    let mut ring = RingPolymer::new(4, MASS, [0.0; 3]);
    for v in ring.velocities.iter_mut() {
        *v = [500.0, -200.0, 50.0];
    }
    let mut forces = vec![[0.0; 3]; 4];
    let mut rng = StdRng::seed_from_u64(12);
    for _ in 0..100 {
        cmd.step(
            &mut ring,
            &mut forces,
            |pos| vec![[0.0; 3]; pos.len()],
            5.0e-17,
            &mut rng,
        )
        .expect("consistent inputs");
    }
    let v = ring.centroid_velocity();
    let expected = [500.0, -200.0, 50.0];
    for d in 0..3 {
        assert!(
            (v[d] - expected[d]).abs() < 1e-9 * expected[d].abs(),
            "centroid velocity {d}: {} vs {}",
            v[d],
            expected[d]
        );
    }
}

#[test]
fn cmd_centroid_oscillates_in_a_harmonic_well_at_the_physical_frequency() {
    // For V = ½ m ω² r² the mean force on the ring is -m ω² r̄: the centroid is
    // an exact classical oscillator whatever the (thermostatted) internal modes do.
    let (p, dt, x0) = (8_usize, 5.0e-17, 1.0e-11);
    let cmd = CentroidMD::new(p, TEMPERATURE, 0.4).expect("valid parameters");
    let mut ring = RingPolymer::new(p, MASS, [x0, 0.0, 0.0]);
    let force_fn = |pos: &[[f64; 3]]| harmonic_forces(pos);
    let mut forces = force_fn(&ring.positions);
    let mut rng = StdRng::seed_from_u64(77);
    let mut max_error = 0.0_f64;
    let mut crossed_zero = false;
    for step in 1..=2000 {
        cmd.step(&mut ring, &mut forces, force_fn, dt, &mut rng)
            .expect("consistent inputs");
        let t = step as f64 * dt;
        let expected = x0 * (OMEGA * t).cos();
        let got = ring.centroid()[0];
        max_error = max_error.max((got - expected).abs());
        crossed_zero |= got < 0.0;
        assert!(ring.centroid()[1].abs() < 1e-3 * x0 && ring.centroid()[2].abs() < 1e-3 * x0);
    }
    assert!(crossed_zero, "the centroid must swing through the origin");
    assert!(max_error < 1e-4 * x0, "max |x̄ - x0 cos ωt| = {max_error}");
}

/// Run CMD in the harmonic well and return the time averages of the bead
/// spread `(1/P) Σ_i |r_i - r̄|²` and of the internal-mode kinetic temperature.
fn harmonic_cmd_averages(p: usize, gamma: f64, dt: f64, steps: usize, seed: u64) -> (f64, f64) {
    let cmd = CentroidMD::new(p, TEMPERATURE, gamma).expect("valid parameters");
    let mut ring = RingPolymer::new(p, MASS, [0.0; 3]);
    let force_fn = |pos: &[[f64; 3]]| harmonic_forces(pos);
    let mut forces = force_fn(&ring.positions);
    let mut rng = StdRng::seed_from_u64(seed);
    let c = cmd.propagator.c_matrix.clone();
    let burn_in = 2_000;
    let (mut spread_sum, mut kinetic_sum, mut count) = (0.0, 0.0, 0.0);
    for step in 0..steps + burn_in {
        cmd.step(&mut ring, &mut forces, force_fn, dt, &mut rng)
            .expect("consistent inputs");
        if step < burn_in {
            continue;
        }
        let centroid = ring.centroid();
        let spread: f64 = ring
            .positions
            .iter()
            .map(|r| (0..3).map(|d| (r[d] - centroid[d]).powi(2)).sum::<f64>())
            .sum::<f64>()
            / p as f64;
        spread_sum += spread;
        let v = beads_to_normal_modes(&ring.velocities, &c);
        let internal_ke: f64 = v[1..]
            .iter()
            .map(|vk| vk.iter().map(|x| x * x).sum::<f64>())
            .sum::<f64>()
            * 0.5
            * cmd.internal_mass(MASS);
        // equipartition over 3 (P - 1) internal degrees of freedom
        kinetic_sum += 2.0 * internal_ke / (3.0 * (p - 1) as f64 * KB);
        count += 1.0;
    }
    (spread_sum / count, kinetic_sum / count)
}

/// Exact `⟨(1/P) Σ_i |r_i - r̄|²⟩ = Σ_{k≥1} 3 k_B T / (m (ω_k² + ω²))`.
fn exact_bead_spread(p: usize) -> f64 {
    (1..p)
        .map(|k| {
            let w_k = mode_frequency(k, p, TEMPERATURE);
            3.0 * KB * TEMPERATURE / (MASS * (w_k * w_k + OMEGA * OMEGA))
        })
        .sum()
}

#[test]
fn cmd_samples_the_path_integral_bead_spread_of_a_harmonic_well() {
    let p = 8;
    let (spread, kinetic_temperature) = harmonic_cmd_averages(p, 0.4, 5.0e-17, 200_000, 11);
    let expected = exact_bead_spread(p);
    assert!(
        (spread - expected).abs() / expected < 0.04,
        "CMD bead spread {spread} vs path-integral {expected}"
    );
    // The internal modes of the P-bead isomorphism sit at the bath temperature P T.
    let bath = cmd_bath(p);
    assert!(
        (kinetic_temperature - bath).abs() / bath < 0.03,
        "internal kinetic temperature {kinetic_temperature} K vs bath {bath} K"
    );
}

fn cmd_bath(p: usize) -> f64 {
    CentroidMD::new(p, TEMPERATURE, 0.5)
        .expect("valid parameters")
        .bath_temperature()
}

#[test]
fn cmd_samples_the_bead_spread_for_odd_bead_counts_and_gamma_one() {
    // gamma = 1 keeps the physical internal masses (plain thermostatted ring
    // polymer dynamics); an odd P exercises the sine/cosine mode pairing.
    let p = 5;
    let (spread, kinetic_temperature) = harmonic_cmd_averages(p, 1.0, 5.0e-17, 200_000, 29);
    let expected = exact_bead_spread(p);
    assert!(
        (spread - expected).abs() / expected < 0.04,
        "CMD bead spread {spread} vs path-integral {expected}"
    );
    let bath = cmd_bath(p);
    assert!((kinetic_temperature - bath).abs() / bath < 0.03);
}

#[test]
fn cmd_with_a_single_bead_is_classical_velocity_verlet() {
    let cmd = CentroidMD::new(1, TEMPERATURE, 1.0).expect("valid parameters");
    let (dt, x0) = (5.0e-17, 1.0e-11);
    let mut ring = RingPolymer::new(1, MASS, [x0, 0.0, 0.0]);
    let force_fn = |pos: &[[f64; 3]]| harmonic_forces(pos);
    let mut forces = force_fn(&ring.positions);
    let mut rng = StdRng::seed_from_u64(2);
    for step in 1..=1500 {
        cmd.step(&mut ring, &mut forces, force_fn, dt, &mut rng)
            .expect("consistent inputs");
        let expected = x0 * (OMEGA * step as f64 * dt).cos();
        assert!((ring.positions[0][0] - expected).abs() < 1e-4 * x0);
    }
    let energy = 0.5 * MASS * ring.velocities[0][0].powi(2)
        + 0.5 * MASS * OMEGA * OMEGA * ring.positions[0][0].powi(2);
    let e0 = 0.5 * MASS * OMEGA * OMEGA * x0 * x0;
    assert!((energy - e0).abs() < 1e-5 * e0, "{energy} vs {e0}");
}

#[test]
fn cmd_trajectory_is_reproducible_for_a_seeded_generator() {
    let run = |seed: u64| {
        let cmd = CentroidMD::new(4, TEMPERATURE, 0.5).expect("valid parameters");
        let mut ring = RingPolymer::new(4, MASS, [1.0e-11, 0.0, 0.0]);
        let force_fn = |pos: &[[f64; 3]]| harmonic_forces(pos);
        let mut forces = force_fn(&ring.positions);
        let mut rng = StdRng::seed_from_u64(seed);
        for _ in 0..50 {
            cmd.step(&mut ring, &mut forces, force_fn, 5.0e-17, &mut rng)
                .expect("consistent inputs");
        }
        ring.positions
    };
    assert_eq!(run(6), run(6));
    assert_ne!(run(6), run(7));
}
