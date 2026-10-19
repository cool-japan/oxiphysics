// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit-system regression tests of the self-contained MD drivers.
//!
//! The drivers work in Å, ps, amu and kJ mol⁻¹, where
//! `1 kJ mol⁻¹ = 100 amu Å² ps⁻²`.  Every expectation below is derived from SI
//! constants or from closed-form physics (springs, constant forces, the
//! equipartition theorem), independently of the conversion factors of the
//! module, so that a missing or doubled factor of 100 is caught.

use super::stochastic::SimRng;
use super::tests_support::{ARGON_MASS, ShiftedForceLj};
use super::*;

/// Atomic mass unit (kg), CODATA 2018.
const AMU_KG: f64 = 1.660_539_066_60e-27;
/// Avogadro constant (mol⁻¹), exact SI value.
const AVOGADRO: f64 = 6.022_140_76e23;
/// Boltzmann constant (J K⁻¹), exact SI value.
const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// Single-particle MdSim in free space.
fn single_particle_sim(
    mass: f64,
    position: [f64; 3],
    velocity: [f64; 3],
    force: [f64; 3],
    dt: f64,
) -> MdSim {
    let config = MdConfig {
        dt,
        n_steps: 0,
        temperature: 0.0,
        pressure: 0.0,
        box_lengths: [1.0e6; 3],
        pbc: false,
        ensemble: Ensemble::NVE,
    };
    let mut state = MdState::new(vec![position], vec![velocity], vec![mass]);
    state.forces = vec![force];
    MdSim::new(config, state)
}

#[test]
fn test_unit_factors_match_the_si_definition_of_the_amu() {
    // One amu Å² ps⁻² per particle is amu · (1e-10 m)² / (1e-12 s)² joule; per
    // mole and in kilojoule that is 0.01 (CODATA 2018 amu: 3.5e-10 relative).
    let kj_per_mol = AMU_KG * 1.0e-20 / 1.0e-24 * AVOGADRO / 1.0e3;
    assert!(
        ((kj_per_mol - KJMOL_PER_AMU_A2_PS2) / KJMOL_PER_AMU_A2_PS2).abs() < 1e-8,
        "1 amu Å² ps⁻² = {kj_per_mol} kJ/mol, constant = {KJMOL_PER_AMU_A2_PS2}"
    );
    assert!((ACCEL_PER_FORCE_OVER_MASS * KJMOL_PER_AMU_A2_PS2 - 1.0).abs() < 1e-15);
    assert_eq!(ACCEL_PER_FORCE_OVER_MASS, 100.0);
}

#[test]
fn test_boltzmann_constant_is_the_molar_gas_constant_in_kj() {
    let gas_constant = BOLTZMANN_J_PER_K * AVOGADRO / 1.0e3;
    assert!(
        ((KB_REDUCED - gas_constant) / gas_constant).abs() < 1e-9,
        "KB_REDUCED = {KB_REDUCED}, R = {gas_constant} kJ/(mol K)"
    );
}

#[test]
fn test_thermal_velocity_sigma_of_argon_matches_si_theory() {
    // sigma = sqrt(k_B T / m) in m/s; one Å/ps is 100 m/s.
    let temperature = 300.0;
    let sigma_si = (BOLTZMANN_J_PER_K * temperature / (ARGON_MASS * AMU_KG)).sqrt();
    let expected = sigma_si / 100.0;
    let sigma = thermal_velocity_sigma(temperature, ARGON_MASS);
    assert!(
        ((sigma - expected) / expected).abs() < 1e-8,
        "sigma = {sigma} Å/ps, SI theory {expected} Å/ps"
    );
    assert!((sigma - 2.4987).abs() < 2e-4, "sigma = {sigma} Å/ps");
}

#[test]
fn test_maxwell_boltzmann_velocities_have_the_physical_thermal_speed() {
    let n = 6000_usize;
    let temperature = 300.0;
    let masses = vec![ARGON_MASS; n];
    let velocities = maxwell_boltzmann_velocities(&masses, temperature, 2024);
    let mean_square: f64 = velocities
        .iter()
        .map(|v| v.iter().map(|c| c * c).sum::<f64>())
        .sum::<f64>()
        / (3.0 * n as f64);
    let sigma = thermal_velocity_sigma(temperature, ARGON_MASS);
    // 18 000 Gaussian components: relative standard error of <v²> is 1.05 %.
    assert!(
        (mean_square.sqrt() / sigma - 1.0).abs() < 0.02,
        "rms component {} Å/ps, sigma {sigma} Å/ps",
        mean_square.sqrt()
    );
    // The temperature recovered from the kinetic energy agrees with the target.
    let sim = MdSim::new(
        MdConfig::default(),
        MdState::new(vec![[0.0; 3]; n], velocities, masses),
    );
    let t = sim.compute_temperature();
    assert!(((t - temperature) / temperature).abs() < 0.03, "T = {t} K");
}

#[test]
fn test_spring_oscillates_with_the_physical_period() {
    // k = 100 kJ mol⁻¹ Å⁻², m = 1 amu: omega = sqrt(100 k / m) = 100 ps⁻¹ and
    // T = 0.2 pi sqrt(m/k) = 0.062 832 ps.  (Without the factor 100 the period
    // would be ten times longer.)
    let well = HarmonicWell::new(100.0);
    let dt = 2.0e-4;
    let mut sim = single_particle_sim(1.0, [0.1, 0.0, 0.0], [0.0; 3], [0.0; 3], dt);
    sim.state.forces = well.compute_all(&sim.state.positions).1;
    let force_fn = |positions: &[[f64; 3]], _box: &[f64; 3]| well.compute_all(positions).1;

    // Times of the downward zero crossings of x(t), by linear interpolation.
    let mut crossings: Vec<f64> = Vec::new();
    let mut previous_x = sim.state.positions[0][0];
    for step in 0..1200_usize {
        sim.velocity_verlet_step(force_fn);
        let x = sim.state.positions[0][0];
        if previous_x > 0.0 && x <= 0.0 {
            let fraction = previous_x / (previous_x - x);
            crossings.push((step as f64 + fraction) * dt);
        }
        previous_x = x;
    }
    assert!(crossings.len() >= 3, "crossings: {crossings:?}");
    let expected = 0.2 * std::f64::consts::PI;
    let expected_period = expected * (1.0_f64 / 100.0).sqrt();
    for pair in crossings.windows(2) {
        let period = pair[1] - pair[0];
        assert!(
            ((period - expected_period) / expected_period).abs() < 1e-3,
            "period {period} ps, expected {expected_period} ps"
        );
    }
    // The analytic helper agrees with the simulation.
    assert!((well.period(1.0) - expected_period).abs() < 1e-12);
}

#[test]
fn test_kinetic_energy_at_the_centre_equals_the_spring_energy_at_the_turning_point() {
    // E = ½ k A² at the turning point must reappear as KE = ½ m v² · 0.01 at the
    // centre of the well, with v = A omega.
    let k = 37.0;
    let mass = 12.0;
    let amplitude = 0.2;
    let well = HarmonicWell::new(k);
    let period = well.period(mass);
    let steps = 4000_usize;
    let dt = period / steps as f64;
    let mut sim = single_particle_sim(mass, [amplitude, 0.0, 0.0], [0.0; 3], [0.0; 3], dt);
    sim.state.forces = well.compute_all(&sim.state.positions).1;
    let force_fn = |positions: &[[f64; 3]], _box: &[f64; 3]| well.compute_all(positions).1;
    let turning_point_energy = well.energy([amplitude, 0.0, 0.0]);
    assert!((turning_point_energy - 0.5 * k * amplitude * amplitude).abs() < 1e-15);

    // A quarter period later the particle is at the centre with maximal speed.
    for _ in 0..steps / 4 {
        sim.velocity_verlet_step(force_fn);
    }
    assert!(sim.state.positions[0][0].abs() < 1e-3 * amplitude);
    let kinetic = sim.compute_kinetic_energy();
    assert!(
        ((kinetic - turning_point_energy) / turning_point_energy).abs() < 1e-4,
        "KE = {kinetic} kJ/mol, spring energy = {turning_point_energy} kJ/mol"
    );
    // The total energy stays constant over the whole period.
    for _ in 0..3 * steps / 4 {
        sim.velocity_verlet_step(force_fn);
        let total = sim.compute_kinetic_energy() + well.energy(sim.state.positions[0]);
        assert!(
            ((total - turning_point_energy) / turning_point_energy).abs() < 1e-4,
            "total energy {total} drifted from {turning_point_energy}"
        );
    }
}

#[test]
fn test_constant_force_acceleration_is_one_hundred_force_over_mass() {
    // F = 1.2 kJ mol⁻¹ Å⁻¹ on 12 amu: a = 100 · 1.2 / 12 = 10 Å ps⁻².
    // After t = 0.1 ps: v = 1 Å/ps and x = ½ a t² = 0.05 Å (velocity Verlet is
    // exact for a constant force).
    let dt = 1.0e-3;
    let force = [1.2, 0.0, 0.0];
    let mut sim = single_particle_sim(12.0, [0.0; 3], [0.0; 3], force, dt);
    let constant = move |positions: &[[f64; 3]], _box: &[f64; 3]| vec![force; positions.len()];
    for _ in 0..100 {
        sim.velocity_verlet_step(constant);
    }
    assert!(
        (sim.state.velocities[0][0] - 1.0).abs() < 1e-12,
        "v = {}",
        sim.state.velocities[0][0]
    );
    assert!(
        (sim.state.positions[0][0] - 0.05).abs() < 1e-12,
        "x = {}",
        sim.state.positions[0][0]
    );

    // The same law in the MdDriver and the plain driver.
    let state = MdState::new(vec![[0.0; 3]], vec![[0.0; 3]], vec![12.0]);
    let mut driver = MdDriver::new(state, dt);
    let fill = move |s: &mut MdState| {
        for f in s.forces.iter_mut() {
            *f = force;
        }
    };
    fill(&mut driver.state);
    for _ in 0..100 {
        driver.velocity_verlet_step(fill);
    }
    assert!((driver.state.velocities[0][0] - 1.0).abs() < 1e-12);
    assert!((driver.state.positions[0][0] - 0.05).abs() < 1e-12);

    let mut atom = Atom::new([0.0; 3], [0.0; 3], 12.0, 0.0);
    atom.force = force;
    let mut plain = MdSimulationPlain::new(vec![atom], [1.0e6; 3]);
    for _ in 0..100 {
        plain.step(dt).expect("frozen-force step");
    }
    assert!((plain.atoms[0].velocity[0] - 1.0).abs() < 1e-12);
    assert!((plain.atoms[0].position[0] - 0.05).abs() < 1e-12);
}

#[test]
fn test_adaptive_timestep_applies_the_same_acceleration() {
    let force = [1.2, 0.0, 0.0];
    let dt = 1.0e-3;
    let mut sim = single_particle_sim(12.0, [0.0; 3], [0.0; 3], force, dt);
    let constant = move |positions: &[[f64; 3]], _box: &[f64; 3]| vec![force; positions.len()];
    // |F| = 1.2 lies between threshold/4 and threshold, so dt is kept.
    let mut elapsed = 0.0;
    for _ in 0..100 {
        elapsed += sim.adaptive_timestep(constant, 2.0, 1.0e-3);
    }
    assert!((elapsed - 0.1).abs() < 1e-12, "elapsed = {elapsed}");
    assert!((sim.state.velocities[0][0] - 1.0).abs() < 1e-12);
    assert!((sim.state.positions[0][0] - 0.05).abs() < 1e-12);
}

#[test]
fn test_andersen_thermostat_samples_the_target_temperature() {
    let n = 4000_usize;
    let target = 200.0;
    let mut state = MdState::new(vec![[0.0; 3]; n], vec![[0.0; 3]; n], vec![12.0; n]);
    // nu · dt = 1: every atom is redrawn from the Maxwell-Boltzmann distribution.
    AndersenThermostat::new(target, 1.0e3).apply(&mut state, 1.0e-3);
    let sim = MdSim::new(MdConfig::default(), state);
    let t = sim.compute_temperature();
    // Relative standard error of T for 12 000 degrees of freedom: 1.3 %.
    assert!(((t - target) / target).abs() < 0.07, "T = {t} K");
}

#[test]
fn test_langevin_o_step_equilibrates_to_the_target_temperature() {
    let n = 4000_usize;
    let target = 150.0;
    let mut state = MdState::new(vec![[0.0; 3]; n], vec![[0.0; 3]; n], vec![39.948; n]);
    let thermostat = LangevinThermostat::new(target, 20.0);
    for _ in 0..60 {
        thermostat.apply_o_step(&mut state, 0.01);
    }
    let sim = MdSim::new(MdConfig::default(), state);
    let t = sim.compute_temperature();
    assert!(((t - target) / target).abs() < 0.07, "T = {t} K");
    // The per-atom target kinetic energy is 3/2 k_B T in kJ/mol.
    let per_atom = sim.compute_kinetic_energy() / n as f64;
    let expected = thermostat.target_kinetic_energy_per_atom();
    assert!(
        ((per_atom - expected) / expected).abs() < 0.07,
        "{per_atom} vs {expected}"
    );
}

#[test]
fn test_pressure_tensor_kinetic_part_is_twice_the_kinetic_energy_density() {
    // Tr Σ m v v / V · 0.01 = 2 KE / V with KE in kJ/mol.
    let velocities = [[1.0, 2.0, -0.5], [-1.5, 0.25, 0.75], [0.3, -0.2, 1.1]];
    let masses = [12.0, 39.948, 1.008];
    let volume = 512.0;
    let tensor = PressureTensor::kinetic_contribution(&velocities, &masses, volume);
    let trace = tensor[0][0] + tensor[1][1] + tensor[2][2];
    let kinetic: f64 = velocities
        .iter()
        .zip(masses.iter())
        .map(|(v, &m)| atom_kinetic_energy(v, m))
        .sum();
    assert!((trace * volume - 2.0 * kinetic).abs() < 1e-14, "{trace}");
    // Symmetric and consistent with the scalar pressure of an ideal gas.
    assert!((tensor[0][1] - tensor[1][0]).abs() < 1e-18);
    let pressure = PressureTensor::from_kinetic_virial(
        &[
            [2.0 * kinetic / 3.0, 0.0, 0.0],
            [0.0, 2.0 * kinetic / 3.0, 0.0],
            [0.0, 0.0, 2.0 * kinetic / 3.0],
        ],
        &[[0.0; 3]; 3],
        volume,
    )
    .scalar_pressure();
    let expected = 2.0 * kinetic / (3.0 * volume) * KJ_MOL_PER_ANG3_TO_BAR;
    assert!(((pressure - expected) / expected).abs() < 1e-12);
}

#[test]
fn test_isothermal_compressibility_is_returned_in_inverse_bar() {
    // Isothermal-isobaric ideal gas: V ~ Gamma(N + 1, k_B T / P), whose
    // fluctuations give kappa_T = Var(V) / (k_B T <V>) = 1 / P exactly.
    let n_atoms = 100.0;
    let temperature = 300.0;
    let pressure_bar = 1000.0;
    let kt = KB_REDUCED * temperature;
    let pressure_internal = pressure_bar / KJ_MOL_PER_ANG3_TO_BAR;
    let mut rng = SimRng::new(77);
    let mut monitor = VolumeFluctuationMonitor::new();
    for _ in 0..300_000 {
        monitor.push(rng.gamma(n_atoms + 1.0) * kt / pressure_internal);
    }
    let kappa = monitor.isothermal_compressibility(temperature);
    let expected = 1.0 / pressure_bar;
    assert!(
        ((kappa - expected) / expected).abs() < 0.03,
        "kappa_T = {kappa} bar⁻¹, expected {expected} bar⁻¹"
    );
    let mean_volume = monitor.mean_volume();
    let expected_volume = (n_atoms + 1.0) * kt / pressure_internal;
    assert!(((mean_volume - expected_volume) / expected_volume).abs() < 0.01);
}

#[test]
fn test_lj_dimer_velocity_verlet_conserves_the_total_energy_in_kj_per_mol() {
    // Two argon atoms oscillating in the LJ well: KE (kJ/mol) + PE (kJ/mol) is
    // constant only if the kick conversion and the energy conversion agree.
    let lj = ShiftedForceLj::argon(2.5);
    let r0 = 3.6;
    let positions = vec![[0.0, 0.0, 0.0], [r0, 0.0, 0.0]];
    let state = MdState::new(positions, vec![[0.0; 3]; 2], vec![ARGON_MASS; 2]);
    let config = MdConfig {
        dt: 0.002,
        n_steps: 0,
        temperature: 0.0,
        pressure: 0.0,
        box_lengths: [200.0; 3],
        pbc: false,
        ensemble: Ensemble::NVE,
    };
    let mut sim = MdSim::new(config, state);
    let virial_fn = lj.virial_fn();
    sim.state.forces = virial_fn(&sim.state.positions, &sim.config.box_lengths).0;
    let energy = |sim: &MdSim| {
        let potential = lj
            .evaluate(&sim.state.positions, &sim.config.box_lengths)
            .potential_energy;
        sim.compute_kinetic_energy() + potential
    };
    let e0 = energy(&sim);
    let mut min_distance = r0;
    let mut max_distance = r0;
    for _ in 0..4000 {
        sim.velocity_verlet_step_virial(&virial_fn);
        let d = sim.state.positions[1][0] - sim.state.positions[0][0];
        min_distance = min_distance.min(d);
        max_distance = max_distance.max(d);
        assert!(
            ((energy(&sim) - e0) / e0.abs()).abs() < 1e-4,
            "energy {} drifted from {e0}",
            energy(&sim)
        );
    }
    // The dimer really moved through the well: it starts at rest at the inner
    // turning point 3.6 Å and swings out to the outer one at 4.19 Å, past the
    // minimum r_min = 2^(1/6) sigma = 3.82 Å (period about 1.3 ps).
    assert!(min_distance > 3.59, "inner turning point {min_distance}");
    assert!(
        max_distance > 4.1 && max_distance < 4.3,
        "outer turning point {max_distance}"
    );
}
