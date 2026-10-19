// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the stochastic routines of the surface-chemistry models.
//!
//! * `Lcg::next_f64` returned `(state >> 33) / (2³² − 1)`, a 31-bit value over a
//!   32-bit range, so every "uniform" number lay in `[0, 0.5)`: Metropolis tests
//!   fired at twice their probability, every `u < 0.5` direction test was true and
//!   the Box–Muller radius was bounded below.
//! * `CatalysisMd::langevin_step` and `SurfaceReconstruction::mc_relax_step` used
//!   `k_B T` in joules next to energies in eV (a temperature of about `1e-20 K`
//!   in eV units), divided a force in eV/Å by a mass in amu without the
//!   `9648.533` conversion to Å/ps², and drew the thermostat noise from a
//!   generator that was reseeded on every call.
//!
//! The generators are seeded, so the statistical checks are deterministic; their
//! tolerances are several standard errors wide.

use super::*;

/// Thermal energy `k_B T` (eV) at `temperature` (K).
fn kt_ev(temperature: f64) -> f64 {
    KB_EV * temperature
}

/// Mean kinetic energy per adsorbate (eV) from the Å/ps velocities and amu masses.
fn mean_kinetic_energy_ev(cat: &CatalysisMd) -> f64 {
    let total: f64 = cat
        .adsorbate_velocities
        .iter()
        .zip(cat.adsorbate_masses.iter())
        .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / EV_TO_AMU_ANG2_PS2)
        .sum();
    total / cat.adsorbate_masses.len() as f64
}

/// Mean `½ m v_a²` (eV) of velocity component `axis`.
fn mean_component_energy_ev(cat: &CatalysisMd, axis: usize) -> f64 {
    let total: f64 = cat
        .adsorbate_velocities
        .iter()
        .zip(cat.adsorbate_masses.iter())
        .map(|(v, &m)| 0.5 * m * v[axis] * v[axis] / EV_TO_AMU_ANG2_PS2)
        .sum();
    total / cat.adsorbate_masses.len() as f64
}

/// `n` adsorbates of mass `mass` (amu) at rest at height `z` (Å).
fn adsorbates_at_rest(n: usize, mass: f64, z: f64) -> CatalysisMd {
    let mut cat = CatalysisMd::new(300.0, -1.5e-19, 1.0e-6, 100);
    for i in 0..n {
        cat.add_adsorbate([i as f64, 0.0, z], [0.0; 3], mass);
    }
    cat
}

// ---------------------------------------------------------------------------
// Lcg
// ---------------------------------------------------------------------------

#[test]
fn lcg_uniform_draws_cover_the_unit_interval() {
    let n = 100_000;
    let mut rng = Lcg::new(1);
    let draws: Vec<f64> = (0..n).map(|_| rng.next_f64()).collect();
    let mean = draws.iter().sum::<f64>() / n as f64;
    let maximum = draws.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let minimum = draws.iter().cloned().fold(f64::INFINITY, f64::min);
    // A 31-bit value divided by 2³² − 1 had mean 0.25 and maximum 0.5.
    assert!((mean - 0.5).abs() < 0.01, "mean = {mean}");
    assert!(maximum > 0.99, "max = {maximum}");
    assert!(minimum < 0.01, "min = {minimum}");
    assert!(draws.iter().all(|&u| (0.0..1.0).contains(&u)));
    // Ten equal bins hold 10 000 draws each (σ = 95).
    let mut bins = [0usize; 10];
    for &u in &draws {
        bins[(u * 10.0) as usize] += 1;
    }
    for (k, &count) in bins.iter().enumerate() {
        assert!(
            (count as f64 - 10_000.0).abs() < 500.0,
            "bin {k} holds {count} of {n} draws"
        );
    }
}

#[test]
fn lcg_range_covers_the_whole_interval() {
    let half_width = 0.3;
    let mut rng = Lcg::new(2);
    let draws: Vec<f64> = (0..100_000)
        .map(|_| rng.range(-half_width, half_width))
        .collect();
    let mean = draws.iter().sum::<f64>() / draws.len() as f64;
    let maximum = draws.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let minimum = draws.iter().cloned().fold(f64::INFINITY, f64::min);
    // [−a, a) was only ever sampled on [−a, 0).
    assert!(mean.abs() < 0.01 * half_width, "mean = {mean}");
    assert!(maximum > 0.99 * half_width, "max = {maximum}");
    assert!(minimum < -0.99 * half_width, "min = {minimum}");
}

#[test]
fn lcg_gaussian_is_a_standard_normal_variate() {
    let n = 100_000;
    let mut rng = Lcg::new(3);
    let draws: Vec<f64> = (0..n).map(|_| rng.next_gaussian()).collect();
    let mean = draws.iter().sum::<f64>() / n as f64;
    let variance = draws.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
    assert!(mean.abs() < 0.015, "mean = {mean}");
    assert!((variance - 1.0).abs() < 0.03, "variance = {variance}");
    // P(|ξ| < 1) = 0.6827 and P(|ξ| < 1.18) = 0.762; with u₁ confined to [0, 0.5) the
    // radius is at least 1.18 and the second fraction was zero.
    let within_one = draws.iter().filter(|x| x.abs() < 1.0).count() as f64 / n as f64;
    let within_1p18 = draws.iter().filter(|x| x.abs() < 1.18).count() as f64 / n as f64;
    assert!(
        (within_one - 0.6827).abs() < 0.01,
        "P(|ξ|<1) = {within_one}"
    );
    assert!(
        (within_1p18 - 0.762).abs() < 0.01,
        "P(|ξ|<1.18) = {within_1p18}"
    );
}

#[test]
fn lcg_streams_of_consecutive_seeds_are_independent() {
    // The step index of a loop is a typical seed: the first draws of the streams of
    // seeds 0, 1, 2, ... must look like independent uniform numbers.
    let first_draws: Vec<f64> = (0..20_000_u64)
        .map(|seed| Lcg::new(seed).next_f64())
        .collect();
    let n = first_draws.len() as f64;
    let mean = first_draws.iter().sum::<f64>() / n;
    assert!((mean - 0.5).abs() < 0.01, "mean = {mean}");
    let covariance = first_draws
        .windows(2)
        .map(|w| (w[0] - mean) * (w[1] - mean))
        .sum::<f64>()
        / (n - 1.0);
    let correlation = covariance / (1.0 / 12.0);
    assert!(
        correlation.abs() < 0.03,
        "lag-1 correlation = {correlation}"
    );
}

#[test]
fn a_repeated_seed_continues_the_stream_and_a_new_seed_restarts_it() {
    let mut slot = None;
    let first = stream_for_seed(&mut slot, 5).next_f64();
    let second = stream_for_seed(&mut slot, 5).next_f64();
    assert_ne!(first, second, "the same seed must continue the stream");
    // A different seed starts that seed's stream from its beginning ...
    let other = stream_for_seed(&mut slot, 6).next_f64();
    assert_eq!(other, Lcg::new(6).next_f64());
    // ... and so does returning to an earlier seed.
    let back = stream_for_seed(&mut slot, 5).next_f64();
    assert_eq!(back, first);
}

// ---------------------------------------------------------------------------
// CatalysisMd::langevin_step
// ---------------------------------------------------------------------------

/// Run free adsorbates (no surface spring) at `temperature` and return the
/// time-averaged mean kinetic energy per adsorbate and its three components.
fn thermalise_free_adsorbates(seed_of_step: impl Fn(usize) -> u64) -> (f64, [f64; 3]) {
    let temperature = 300.0;
    let mut cat = adsorbates_at_rest(1_000, 18.0, 5.0);
    cat.temperature = temperature;
    let (dt, gamma) = (0.01, 5.0);
    for step in 0..400 {
        cat.langevin_step(dt, 0.0, gamma, seed_of_step(step));
    }
    let samples = 200;
    let mut kinetic = 0.0;
    let mut components = [0.0; 3];
    for step in 400..400 + samples {
        cat.langevin_step(dt, 0.0, gamma, seed_of_step(step));
        kinetic += mean_kinetic_energy_ev(&cat);
        for (axis, component) in components.iter_mut().enumerate() {
            *component += mean_component_energy_ev(&cat, axis);
        }
    }
    let inverse = 1.0 / samples as f64;
    (
        kinetic * inverse,
        [
            components[0] * inverse,
            components[1] * inverse,
            components[2] * inverse,
        ],
    )
}

#[test]
fn free_adsorbates_reach_equipartition_with_a_constant_seed() {
    // ⟨½ m v²⟩ = 3/2 k_B T.  The former noise amplitude was sqrt(2 γ k_B T m / dt)
    // with k_B in J/K (about 1e-10 in these units): the adsorbates never moved.
    // A constant seed on every call must still draw fresh noise at every step.
    let (kinetic, components) = thermalise_free_adsorbates(|_| 12_345);
    let kt = kt_ev(300.0);
    assert!(
        (kinetic - 1.5 * kt).abs() < 0.05 * 1.5 * kt,
        "<1/2 m v^2> = {kinetic} eV, expected {} eV",
        1.5 * kt
    );
    for (axis, &component) in components.iter().enumerate() {
        assert!(
            (component - 0.5 * kt).abs() < 0.06 * 0.5 * kt,
            "axis {axis}: <1/2 m v_a^2> = {component} eV, expected {} eV",
            0.5 * kt
        );
    }
}

#[test]
fn free_adsorbates_reach_equipartition_with_step_index_seeds() {
    let (kinetic, _) = thermalise_free_adsorbates(|step| step as u64);
    let kt = kt_ev(300.0);
    assert!(
        (kinetic - 1.5 * kt).abs() < 0.05 * 1.5 * kt,
        "<1/2 m v^2> = {kinetic} eV, expected {} eV",
        1.5 * kt
    );
}

#[test]
fn free_adsorbate_temperature_follows_the_thermostat_setting() {
    // Doubling the temperature doubles the mean kinetic energy.
    let mut kinetic = Vec::new();
    for temperature in [150.0, 300.0, 600.0] {
        let mut cat = adsorbates_at_rest(1_000, 18.0, 5.0);
        cat.temperature = temperature;
        for step in 0..500 {
            cat.langevin_step(0.01, 0.0, 5.0, step);
        }
        let mut accumulated = 0.0;
        for step in 500..700 {
            cat.langevin_step(0.01, 0.0, 5.0, step);
            accumulated += mean_kinetic_energy_ev(&cat);
        }
        kinetic.push(accumulated / 200.0);
    }
    assert!((kinetic[1] / kinetic[0] - 2.0).abs() < 0.15, "{kinetic:?}");
    assert!((kinetic[2] / kinetic[1] - 2.0).abs() < 0.15, "{kinetic:?}");
}

#[test]
fn surface_spring_gives_equipartition_in_height_and_velocity() {
    // Harmonic well F_z = −k z with k = 1 eV/Å² and m = 18 amu: ω = 23 ps⁻¹.  Without
    // the 9648.533 eV/(Å amu) → Å/ps² conversion of the force the oscillator was
    // some 9648 times too soft and z wandered over hundreds of Å.
    let temperature = 300.0;
    let k_surf = 1.0;
    let mut cat = adsorbates_at_rest(1_000, 18.0, 0.0);
    cat.temperature = temperature;
    let (dt, gamma) = (0.004, 10.0);
    for _ in 0..400 {
        cat.langevin_step(dt, k_surf, gamma, 7);
    }
    let samples = 400;
    let (mut potential, mut kinetic_z) = (0.0, 0.0);
    for _ in 0..samples {
        cat.langevin_step(dt, k_surf, gamma, 7);
        potential += cat
            .adsorbate_positions
            .iter()
            .map(|p| 0.5 * k_surf * p[2] * p[2])
            .sum::<f64>()
            / cat.adsorbate_positions.len() as f64;
        kinetic_z += mean_component_energy_ev(&cat, 2);
    }
    let (potential, kinetic_z) = (potential / samples as f64, kinetic_z / samples as f64);
    let half_kt = 0.5 * kt_ev(temperature);
    assert!(
        (potential - half_kt).abs() < 0.06 * half_kt,
        "<1/2 k z^2> = {potential} eV, expected {half_kt} eV"
    );
    assert!(
        (kinetic_z - half_kt).abs() < 0.06 * half_kt,
        "<1/2 m v_z^2> = {kinetic_z} eV, expected {half_kt} eV"
    );
}

#[test]
fn frictionless_surface_oscillator_conserves_energy_and_keeps_its_frequency() {
    // γ = 0: the BAOAB step is velocity Verlet.  The energy ½ m v² / 9648.533 + ½ k z²
    // (eV) is only conserved if the force and the kinetic energy use the same unit
    // conversion, and z(t) = z₀ cos(ω t) with ω² = 9648.533 k / m.
    let (k_surf, mass, z0, dt) = (1.0, 18.0, 0.1, 0.002);
    let mut cat = CatalysisMd::new(300.0, -1.5e-19, 1.0e-6, 100);
    cat.add_adsorbate([0.0, 0.0, z0], [0.0; 3], mass);
    let omega = (EV_TO_AMU_ANG2_PS2 * k_surf / mass).sqrt();
    let energy = |cat: &CatalysisMd| {
        let v = cat.adsorbate_velocities[0];
        let z = cat.adsorbate_positions[0][2];
        0.5 * mass * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / EV_TO_AMU_ANG2_PS2
            + 0.5 * k_surf * z * z
    };
    let e0 = energy(&cat);
    let steps = 100;
    for step in 1..=steps {
        cat.langevin_step(dt, k_surf, 0.0, 1);
        let drift = (energy(&cat) - e0).abs() / e0;
        assert!(drift < 2.0e-3, "step {step}: relative energy drift {drift}");
    }
    let t = steps as f64 * dt;
    let z = cat.adsorbate_positions[0][2];
    assert!(
        (z - z0 * (omega * t).cos()).abs() < 1.0e-3 * z0,
        "z({t}) = {z}, expected {}",
        z0 * (omega * t).cos()
    );
    // No force acts along x and y, and the noise is off: no drift there.
    assert_eq!(cat.adsorbate_positions[0][0], 0.0);
    assert_eq!(cat.adsorbate_velocities[0][0], 0.0);
}

#[test]
fn langevin_step_is_reproducible_for_a_given_sequence_of_seeds() {
    let run = |seeds: &[u64]| -> Vec<[f64; 3]> {
        let mut cat = adsorbates_at_rest(20, 18.0, 1.0);
        for &seed in seeds {
            cat.langevin_step(0.01, 1.0, 5.0, seed);
        }
        cat.adsorbate_velocities
    };
    assert_eq!(run(&[1, 1, 2, 3]), run(&[1, 1, 2, 3]));
    assert_ne!(run(&[1, 1, 2, 3]), run(&[7, 7, 8, 9]));
}

#[test]
fn langevin_step_skips_massless_adsorbates() {
    let mut cat = CatalysisMd::new(300.0, -1.5e-19, 1.0e-6, 100);
    cat.add_adsorbate([0.0, 0.0, 1.0], [0.5, 0.0, 0.0], 0.0);
    cat.langevin_step(0.01, 1.0, 5.0, 1);
    assert_eq!(cat.adsorbate_positions[0], [0.0, 0.0, 1.0]);
    assert_eq!(cat.adsorbate_velocities[0], [0.5, 0.0, 0.0]);
}

// ---------------------------------------------------------------------------
// CatalysisMd::simulate_reactions
// ---------------------------------------------------------------------------

#[test]
fn reaction_events_fire_with_probability_rate_times_dt() {
    // Zero activation energy: k = A = 1e13 s⁻¹ and dt = 0.2e-13 s give a probability
    // of 0.2 per site and step; a generator on [0, 0.5) fired at 0.4.
    let mut cat = CatalysisMd::new(300.0, -1.5e-19, 1.0e-6, 1_000);
    cat.add_step(ReactionStep::new("barrierless", 0.0, 1.0e13, -1.0e-19));
    cat.simulate_reactions(100, 0.2e-13, 17);
    let trials: f64 = 1_000.0 * 100.0;
    let expected = 0.2 * trials;
    let sigma = (trials * 0.2 * 0.8).sqrt();
    assert!(
        (cat.total_turnovers as f64 - expected).abs() < 5.0 * sigma,
        "turnovers {}, expected {expected} (σ = {sigma})",
        cat.total_turnovers
    );
}

// ---------------------------------------------------------------------------
// SurfaceReconstruction::mc_relax_step
// ---------------------------------------------------------------------------

/// Surface atoms at their equilibrium positions with the spring constant
/// `spring_constant` (eV/Å²) at `temperature` (K).
fn reconstruction(n: usize, spring_constant: f64, temperature: f64) -> SurfaceReconstruction {
    let atoms: Vec<SurfaceAtom> = (0..n)
        .map(|i| SurfaceAtom::new([i as f64, 0.0, 2.0], -3.5))
        .collect();
    SurfaceReconstruction::new(atoms, 0.5, spring_constant, temperature)
}

/// Mean harmonic energy per atom (eV) after relaxing the surface by Metropolis
/// sweeps; `seed_of_sweep` maps the sweep index to the seed passed to the routine.
fn sampled_energy_per_atom(
    surface: &mut SurfaceReconstruction,
    amplitude: f64,
    seed_of_sweep: impl Fn(usize) -> u64,
) -> f64 {
    let n = surface.atoms.len() as f64;
    for sweep in 0..150 {
        surface.mc_relax_step(amplitude, seed_of_sweep(sweep));
    }
    let samples = 300;
    let mut accumulated = 0.0;
    for sweep in 150..150 + samples {
        surface.mc_relax_step(amplitude, seed_of_sweep(sweep));
        accumulated += surface.potential_energy() / n;
    }
    accumulated / samples as f64
}

#[test]
fn mc_relax_samples_the_boltzmann_distribution_of_the_harmonic_well() {
    // ⟨½ k d²⟩ = 3/2 k_B T for a 3-D harmonic well.  With k_B T in joules the
    // Metropolis factor exp(−ΔE / k_B T) vanished for every uphill move and the
    // surface was quenched to zero temperature.
    for temperature in [150.0, 300.0, 600.0] {
        let mut surface = reconstruction(1_000, 10.0, temperature);
        // Proposal width of about two thermal standard deviations sqrt(k_B T / k).
        let amplitude = 2.0 * (kt_ev(temperature) / 10.0).sqrt();
        let mean = sampled_energy_per_atom(&mut surface, amplitude, |sweep| sweep as u64);
        let expected = 1.5 * kt_ev(temperature);
        assert!(
            (mean - expected).abs() < 0.05 * expected,
            "T = {temperature} K: <1/2 k d^2> = {mean} eV, expected {expected} eV"
        );
    }
}

#[test]
fn mc_relax_with_a_constant_seed_still_samples_the_boltzmann_distribution() {
    let mut surface = reconstruction(1_000, 10.0, 300.0);
    let amplitude = 2.0 * (kt_ev(300.0) / 10.0).sqrt();
    let mean = sampled_energy_per_atom(&mut surface, amplitude, |_| 99);
    let expected = 1.5 * kt_ev(300.0);
    assert!(
        (mean - expected).abs() < 0.05 * expected,
        "<1/2 k d^2> = {mean} eV, expected {expected} eV"
    );
}

#[test]
fn mc_relax_accepts_downhill_moves_and_rejects_far_uphill_moves() {
    // A very stiff well at 1 K: any proposal of 0.1 Å costs far more than k_B T
    // and is rejected, so every atom stays where it is.
    let mut frozen = reconstruction(200, 1.0e6, 1.0);
    let accepted = frozen.mc_relax_step(0.1, 5);
    assert_eq!(accepted, 0);
    assert!(frozen.atoms.iter().all(|a| a.displacement() == 0.0));
    // Atoms displaced by 1 Å from the minimum of a 10 eV/Å² well at 1 K: a symmetric
    // proposal of amplitude a lowers the energy with probability 1/2 − a/4 (the
    // x-step has to beat the quadratic term) and every uphill step, costing at least
    // meV, is rejected.
    let amplitude = 0.05;
    let mut displaced = reconstruction(2_000, 10.0, 1.0);
    for atom in displaced.atoms.iter_mut() {
        atom.pos[0] += 1.0;
    }
    let accepted = displaced.mc_relax_step(amplitude, 6);
    let fraction = accepted as f64 / 2_000.0;
    let expected = 0.5 - amplitude / 4.0;
    assert!(
        (fraction - expected).abs() < 0.05,
        "accepted fraction {fraction}, expected {expected}"
    );
}

// ---------------------------------------------------------------------------
// SurfaceDiffusion
// ---------------------------------------------------------------------------

/// Flat-barrier diffusion model whose adatoms hop on every call (zero barrier,
/// `k dt = 10 > 1` is clamped to a hop probability of one).
fn always_hopping(n_adatoms: usize) -> SurfaceDiffusion {
    SurfaceDiffusion::new(
        vec![500.0; n_adatoms],
        vec![0.0; 1000],
        300.0,
        1.0e13,
        3.0e-10,
    )
}

/// Mean displacement and mean square displacement of an adatom ensemble after
/// `steps` hops, drawing the seed of step `i` as `seed_of_step(i)`.
fn walk(n_adatoms: usize, steps: usize, seed_of_step: impl Fn(usize) -> u64) -> (f64, f64) {
    let mut model = always_hopping(n_adatoms);
    let initial = model.positions.clone();
    for step in 0..steps {
        let hops = model.monte_carlo_step(1.0e-12, seed_of_step(step));
        assert_eq!(hops, n_adatoms);
    }
    let mean = model
        .positions
        .iter()
        .zip(initial.iter())
        .map(|(x, x0)| x - x0)
        .sum::<f64>()
        / n_adatoms as f64;
    (mean, model.mean_square_displacement(&initial))
}

#[test]
fn hops_are_equally_likely_in_both_directions() {
    // Unit hops in a random direction on every call: a symmetric random walk with
    // mean 0 and MSD = number of steps.  A generator on [0, 0.5) always chose +1.
    let (n_adatoms, steps) = (4_000, 500);
    let (mean, msd) = walk(n_adatoms, steps, |step| step as u64);
    let standard_error = (msd / n_adatoms as f64).sqrt();
    assert!(
        mean.abs() < 4.0 * standard_error,
        "drift {mean} exceeds four standard errors ({standard_error})"
    );
    assert!(
        (msd - steps as f64).abs() < 0.1 * steps as f64,
        "MSD = {msd}, expected {steps}"
    );
}

#[test]
fn hops_are_diffusive_with_a_constant_seed() {
    // A generator reseeded on every call replays the same hop directions: every
    // adatom then moves ballistically, MSD = steps² instead of steps.
    let (n_adatoms, steps) = (4_000, 500);
    let (mean, msd) = walk(n_adatoms, steps, |_| 42);
    let standard_error = (msd / n_adatoms as f64).sqrt();
    assert!(
        mean.abs() < 4.0 * standard_error,
        "drift {mean} exceeds four standard errors ({standard_error})"
    );
    assert!(
        (msd - steps as f64).abs() < 0.1 * steps as f64,
        "MSD = {msd}, expected {steps}"
    );
}

#[test]
fn flat_barriers_give_no_drift_over_a_thousand_steps() {
    // The configuration of the original finding: 1000 adatoms, 1000 steps.
    let (n_adatoms, steps) = (1_000, 1_000);
    let (mean, msd) = walk(n_adatoms, steps, |step| step as u64);
    assert!(
        mean.abs() < 3.0 * (msd / n_adatoms as f64).sqrt(),
        "drift {mean}, MSD {msd}"
    );
}

#[test]
fn hop_probability_equals_rate_times_time_step() {
    // Zero barrier: k = ν₀ = 1e13 / s and dt = 0.3e-13 s give a hop probability of 0.3;
    // a generator on [0, 0.5) hopped with probability 0.6.
    let n_adatoms = 20_000;
    let mut model = always_hopping(n_adatoms);
    let hops = model.monte_carlo_step(0.3e-13, 7);
    let fraction = hops as f64 / n_adatoms as f64;
    let sigma = (0.3_f64 * 0.7 / n_adatoms as f64).sqrt();
    assert!(
        (fraction - 0.3).abs() < 5.0 * sigma,
        "hopping fraction {fraction}, expected 0.3 (σ = {sigma})"
    );
}

#[test]
fn precursors_chemisorb_with_probability_rate_times_dt() {
    // e_chem = 0: k = ν₀ = 1e13 / s, dt = 0.25e-13 s → probability 0.25.  The draws
    // that decide the diffusion of the precursor come first, so this also checks that
    // the chemisorption test is fed uniform numbers.
    let n_precursors = 20_000;
    let mut model = SurfaceDiffusion::new(Vec::new(), vec![0.0], 300.0, 1.0e13, 3.0e-10);
    for _ in 0..n_precursors {
        model.add_precursor(0.0, -5.0e-21);
    }
    let chemisorbed = model.precursor_step(0.0, 0.25e-13, 11);
    let fraction = chemisorbed as f64 / n_precursors as f64;
    let sigma = (0.25_f64 * 0.75 / n_precursors as f64).sqrt();
    assert!(
        (fraction - 0.25).abs() < 5.0 * sigma,
        "chemisorbed fraction {fraction}, expected 0.25 (σ = {sigma})"
    );
    assert!((model.chemisorption_fraction() - fraction).abs() < 1e-12);
}

#[test]
fn diffusion_runs_are_reproducible_for_a_given_sequence_of_seeds() {
    let run = |seeds: &[u64]| -> Vec<f64> {
        let mut model = always_hopping(100);
        for &seed in seeds {
            model.monte_carlo_step(1.0e-12, seed);
        }
        model.positions
    };
    assert_eq!(run(&[1, 1, 2, 3]), run(&[1, 1, 2, 3]));
    assert_ne!(run(&[1, 1, 2, 3]), run(&[7, 7, 8, 9]));
}

// ---------------------------------------------------------------------------
// Tribochemistry::step
// ---------------------------------------------------------------------------

/// Tribochemistry whose bond breaking is stress independent (`E_a = 0`, no
/// activation volume): `k = ν₀ = 1e13 / s`.
fn barrierless_tribo() -> Tribochemistry {
    Tribochemistry::new(
        10.0,
        0.1,
        0.3,
        300.0,
        0.0,
        0.0,
        1.0e13,
        1.0e9,
        TriboFilm::new(0.0, 0.0, 0.0, 0.0, 0.0),
    )
}

#[test]
fn tribochemistry_breaks_bonds_and_sheds_particles_with_the_modelled_probabilities() {
    // dt = 0.3e-13 s: a bond breaks with probability k dt = 0.3 per step and a wear
    // particle is shed with probability 0.1.  The same seed is passed at every step, so
    // a generator that was reseeded on every call would draw one and the same number
    // for ever: bonds_broken would be 0 or 20 000.
    let steps = 20_000;
    let mut tribo = barrierless_tribo();
    for step in 0..steps {
        tribo.step(0.3e-13, 1.0e-6, 1.0e-4, step, 99);
    }
    let n = steps as f64;
    let sigma_bonds = (n * 0.3 * 0.7).sqrt();
    assert!(
        (tribo.bonds_broken as f64 - 0.3 * n).abs() < 5.0 * sigma_bonds,
        "{} bonds broke, expected {}",
        tribo.bonds_broken,
        0.3 * n
    );
    let sigma_particles = (n * 0.1 * 0.9).sqrt();
    assert!(
        (tribo.n_wear_particles() as f64 - 0.1 * n).abs() < 5.0 * sigma_particles,
        "{} particles, expected {}",
        tribo.n_wear_particles(),
        0.1 * n
    );
    // Particle volumes are uniform on [0.5, 1.5] of the volume per step, mean 1.
    let volume_per_step = 1.0e-4 * 10.0 / 1.0e9 * 0.1 * 0.3e-13;
    let ratio = tribo.mean_particle_volume() / volume_per_step;
    assert!(
        (ratio - 1.0).abs() < 0.05,
        "mean volume / nominal = {ratio}"
    );
    assert!(
        (tribo.wear_volume - tribo.wear_particles.iter().map(|p| p.volume).sum::<f64>()).abs()
            < 1.0e-12 * tribo.wear_volume
    );
}
