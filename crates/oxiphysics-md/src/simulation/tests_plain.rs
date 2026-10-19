// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of [`MdSimulationPlain`]: velocity Verlet with the forces
//! recomputed at the new positions, the canonical thermostats (CSVR and the
//! Nosé-Hoover chain), the barostats that act on the pair virial, and the error
//! paths of invalid configurations.
//!
//! Every expectation is a closed form -- the harmonic oscillator, the ideal gas
//! `P V = N k_B T`, the gamma law of the canonical kinetic energy, the
//! isothermal-isobaric volume distribution of an ideal gas -- so the tests fail
//! for the earlier behaviour: stale forces in both half kicks, Berendsen
//! rescaling under every thermostat name, and a barostat field that was never
//! read.

use super::stochastic::SimRng;
use super::tests_support::{ARGON_MASS, ShiftedForceLj, argon_fcc_fluid, reduced_density};
use super::*;

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

/// Force routine of non-interacting atoms: zero forces, zero virial, zero
/// potential energy.
fn free_gas(atoms: &[Atom], _box_lengths: &[f64; 3]) -> ForceEvaluation {
    ForceEvaluation::new(vec![[0.0; 3]; atoms.len()])
        .with_virial(0.0)
        .with_potential_energy(0.0)
}

/// `n` free argon atoms at uniformly random positions in a cube of edge `edge`,
/// with exactly thermal velocities at `temperature` (K) scaled by `heat`.
fn ideal_gas_plain(
    n: usize,
    temperature: f64,
    heat: f64,
    edge: f64,
    seed: u64,
) -> MdSimulationPlain {
    let masses = vec![ARGON_MASS; n];
    let mut velocities = maxwell_boltzmann_velocities(&masses, temperature, seed);
    rescale_to_temperature(&mut velocities, &masses, heat * temperature);
    let mut rng = SimRng::new(seed ^ 0x5A5A_5A5A);
    let atoms = velocities
        .into_iter()
        .map(|velocity| {
            let position = [
                rng.uniform() * edge,
                rng.uniform() * edge,
                rng.uniform() * edge,
            ];
            Atom::new(position, velocity, ARGON_MASS, 0.0)
        })
        .collect();
    let mut sim = MdSimulationPlain::new(atoms, [edge; 3]).with_seed(seed);
    sim.target_temp = temperature;
    sim
}

/// Centre of the oscillator well; the driver wraps coordinates into `[0, L)`, so
/// the well sits in the middle of the box instead of at the origin.
const WELL_CENTRE: f64 = 50.0;

/// One atom of mass 1 amu at rest at the displacement `x0` from the centre of the
/// isotropic well `k = 100 kJ mol⁻¹ Å⁻²`, i.e. `ω = 100 ps⁻¹` (period 0.0628 ps).
fn oscillator(x0: f64) -> (MdSimulationPlain, HarmonicWell) {
    let atom = Atom::new(
        [WELL_CENTRE + x0, WELL_CENTRE, WELL_CENTRE],
        [0.0; 3],
        1.0,
        0.0,
    );
    (
        MdSimulationPlain::new(vec![atom], [2.0 * WELL_CENTRE; 3]),
        HarmonicWell::with_center(100.0, [WELL_CENTRE; 3]),
    )
}

/// Displacement of the oscillator atom from the centre of its well (Å).
fn displacement(sim: &MdSimulationPlain) -> f64 {
    sim.atoms[0].position[0] - WELL_CENTRE
}

/// Force routine of the external well (forces and potential energy; a virial
/// does not exist for an external field).
fn well_forces(well: HarmonicWell) -> impl FnMut(&[Atom], &[f64; 3]) -> ForceEvaluation {
    move |atoms, _box_lengths| {
        let positions: Vec<[f64; 3]> = atoms.iter().map(|atom| atom.position).collect();
        let (energy, forces) = well.compute_all(&positions);
        ForceEvaluation::new(forces).with_potential_energy(energy)
    }
}

/// A small argon cluster in free space (box of 1000 Å, centred), at rest apart
/// from thermal velocities at `temperature` (K); the Lennard-Jones interaction
/// is the shifted-force 12-6 potential with the cutoff `2.5 σ`.
fn argon_cluster(temperature: f64, seed: u64) -> (MdSimulationPlain, ShiftedForceLj) {
    let centre = 500.0;
    let edge = 3.9;
    let offsets = [
        [0.0, 0.0, 0.0],
        [edge, 0.0, 0.0],
        [0.5 * edge, 0.866_025_403_784_438_6 * edge, 0.0],
        [
            0.5 * edge,
            0.288_675_134_594_812_9 * edge,
            0.816_496_580_927_726 * edge,
        ],
    ];
    let masses = vec![ARGON_MASS; offsets.len()];
    let mut velocities = maxwell_boltzmann_velocities(&masses, temperature, seed);
    rescale_to_temperature(&mut velocities, &masses, temperature);
    let atoms = offsets
        .iter()
        .zip(velocities)
        .map(|(offset, velocity)| {
            let position = [centre + offset[0], centre + offset[1], centre + offset[2]];
            Atom::new(position, velocity, ARGON_MASS, 0.0)
        })
        .collect();
    (
        MdSimulationPlain::new(atoms, [1000.0; 3]),
        ShiftedForceLj::argon(2.5),
    )
}

/// Mean and (population) variance of the kinetic energy over `steps` steps of a
/// non-interacting gas, after `warm_up` steps of equilibration.
fn gas_kinetic_statistics(
    sim: &mut MdSimulationPlain,
    dt: f64,
    warm_up: u64,
    steps: u64,
) -> (f64, f64) {
    sim.run_with(warm_up, dt, free_gas).expect("warm-up");
    let (mut sum, mut sum_sq) = (0.0, 0.0);
    for _ in 0..steps {
        sim.step_with(dt, free_gas).expect("sampling step");
        let kinetic = sim.kinetic_energy();
        sum += kinetic;
        sum_sq += kinetic * kinetic;
    }
    let mean = sum / steps as f64;
    (mean, sum_sq / steps as f64 - mean * mean)
}

/// Argon-like Lennard-Jones fluid on a 3×3×3 FCC lattice (108 atoms).
fn lj_fluid_plain(
    density: f64,
    temperature: f64,
    seed: u64,
) -> (MdSimulationPlain, ShiftedForceLj) {
    let (atoms, box_lengths) = argon_fcc_fluid(3, density, temperature, seed);
    let mut sim = MdSimulationPlain::new(atoms, box_lengths).with_seed(seed);
    sim.target_temp = temperature;
    (sim, ShiftedForceLj::argon(2.0))
}

// ---------------------------------------------------------------------------
// Velocity Verlet with recomputed forces
// ---------------------------------------------------------------------------

#[test]
fn test_step_with_is_second_order_accurate_velocity_verlet() {
    // x(t) = x0 cos(omega t) with omega = 100 ps⁻¹.  Velocity Verlet has a global
    // error proportional to dt²: halving dt divides it by four.  Two half kicks with
    // the same stale force are only first order (ratio two).
    let t_end = 0.05;
    let x0 = 0.1;
    let exact = x0 * (100.0_f64 * t_end).cos();
    let error_for = |steps: u64| {
        let (mut sim, well) = oscillator(x0);
        sim.run_with(steps, t_end / steps as f64, well_forces(well))
            .expect("valid run");
        (displacement(&sim) - exact).abs()
    };
    let coarse = error_for(50);
    let fine = error_for(100);
    let finest = error_for(200);
    assert!(
        (coarse / fine - 4.0).abs() < 0.1,
        "error ratio {} for dt -> dt/2, expected 4",
        coarse / fine
    );
    assert!(
        (fine / finest - 4.0).abs() < 0.1,
        "error ratio {} for dt -> dt/2, expected 4",
        fine / finest
    );
    assert!(finest < 2.0e-5, "absolute error {finest} Å");
}

#[test]
fn test_step_with_oscillator_has_the_physical_period_and_conserves_energy() {
    let dt = 2.0e-4;
    let (mut sim, well) = oscillator(0.1);
    let mut evaluate = well_forces(well);
    let e0 = well.energy(sim.atoms[0].position);
    let mut crossings: Vec<f64> = Vec::new();
    let mut previous_x = displacement(&sim);
    for step in 0..1200_usize {
        sim.step_with(dt, &mut evaluate).expect("valid step");
        let x = displacement(&sim);
        if previous_x > 0.0 && x <= 0.0 {
            crossings.push((step as f64 + previous_x / (previous_x - x)) * dt);
        }
        previous_x = x;
        let total = sim.total_energy().expect("potential energy supplied");
        assert!(
            ((total - e0) / e0).abs() < 1e-3,
            "total energy {total} drifted from {e0}"
        );
    }
    assert!(crossings.len() >= 3, "crossings: {crossings:?}");
    let expected_period = well.period(1.0);
    for pair in crossings.windows(2) {
        let period = pair[1] - pair[0];
        assert!(
            ((period - expected_period) / expected_period).abs() < 1e-3,
            "period {period} ps, expected {expected_period} ps"
        );
    }
}

#[test]
fn test_step_with_conserves_the_total_energy_of_an_argon_cluster() {
    let (mut sim, lj) = argon_cluster(40.0, 17);
    let mut evaluate = lj.plain_fn();
    sim.refresh_forces(&mut evaluate).expect("force evaluation");
    let e0 = sim.total_energy().expect("potential energy supplied");
    assert!(e0 < 0.0, "a bound cluster has negative energy: {e0}");
    let mut worst: f64 = 0.0;
    for _ in 0..4000 {
        sim.step_with(0.002, &mut evaluate).expect("valid step");
        let e = sim.total_energy().expect("potential energy supplied");
        worst = worst.max(((e - e0) / e0).abs());
    }
    assert!(worst < 1e-4, "relative energy drift {worst}");
    // The cluster really moved and exchanged kinetic and potential energy.
    assert!(sim.kinetic_energy() > 0.0);
}

#[test]
fn test_step_with_is_time_reversible() {
    // Velocity Verlet with forces recomputed between the kicks is time reversible:
    // propagating, flipping the velocities and propagating again returns to the
    // start up to round-off.  A scheme that kicks twice with the same stale force
    // is not.
    let (mut sim, lj) = argon_cluster(60.0, 3);
    let mut evaluate = lj.plain_fn();
    let start: Vec<[f64; 3]> = sim.atoms.iter().map(|atom| atom.position).collect();
    sim.run_with(1500, 0.002, &mut evaluate).expect("forward");
    let moved: f64 = sim
        .atoms
        .iter()
        .zip(&start)
        .map(|(atom, p)| {
            (0..3)
                .map(|a| (atom.position[a] - p[a]).powi(2))
                .sum::<f64>()
        })
        .sum::<f64>()
        .sqrt();
    assert!(moved > 0.3, "the cluster must have moved: {moved} Å");
    for atom in &mut sim.atoms {
        for v in atom.velocity.iter_mut() {
            *v = -*v;
        }
    }
    sim.run_with(1500, 0.002, &mut evaluate).expect("backward");
    for (atom, p) in sim.atoms.iter().zip(&start) {
        for (a, (now, before)) in atom.position.iter().zip(p.iter()).enumerate() {
            assert!(
                (now - before).abs() < 1e-7,
                "axis {a}: {now} vs start {before}"
            );
        }
    }
}

#[test]
fn test_invalidated_forces_are_recomputed_before_the_next_step() {
    // The first half kick of a step uses the force at the *current* position.  After
    // the position was changed by hand and the stored forces were invalidated, the
    // step must start from a fresh evaluation.
    let (mut sim, well) = oscillator(0.1);
    let mut evaluate = well_forces(well);
    sim.step_with(1.0e-4, &mut evaluate).expect("first step");
    sim.atoms[0].position = [WELL_CENTRE + 0.3, WELL_CENTRE, WELL_CENTRE];
    sim.atoms[0].velocity = [0.0; 3];
    sim.invalidate_forces();
    let dt = 1.0e-4;
    sim.step_with(dt, &mut evaluate).expect("second step");
    // Two half kicks with the force at 0.3 (and a negligible change over the step):
    // v = 100 * F/m * dt = 100 * (-100 * 0.3) * 1e-4 = -0.3 Å/ps.
    let v = sim.atoms[0].velocity[0];
    assert!((v + 0.3).abs() < 2e-3, "v = {v}, expected about -0.3");
    // Without the invalidation the stale force at 0.1 would give -0.1.
    assert!(v < -0.25);
}

// ---------------------------------------------------------------------------
// Canonical thermostats
// ---------------------------------------------------------------------------

#[test]
fn test_csvr_reproduces_the_canonical_kinetic_energy_distribution() {
    // The canonical kinetic energy of N_f = 3N degrees of freedom is gamma
    // distributed: <K> = N_f kT / 2 and var K = (2 / N_f) <K>².  A non-interacting
    // gas exchanges no energy between its degrees of freedom, so only a stochastic
    // thermostat can produce the fluctuations.
    let n = 20_usize;
    let temperature = 250.0;
    let mut sim = ideal_gas_plain(n, temperature, 1.0, 40.0, 3);
    sim.thermostat = Some(ThermostatType::Csvr);
    sim.tau_t = 0.1;
    let n_dof = 3.0 * n as f64;
    let target = 0.5 * n_dof * KB_REDUCED * temperature;
    let (mean, variance) = gas_kinetic_statistics(&mut sim, 0.01, 2000, 400_000);
    let canonical_variance = 2.0 * target * target / n_dof;
    assert!(
        ((mean - target) / target).abs() < 0.02,
        "<K> = {mean} kJ/mol, expected {target}"
    );
    assert!(
        ((variance - canonical_variance) / canonical_variance).abs() < 0.10,
        "var K = {variance}, expected (2/3N) <K>² = {canonical_variance}"
    );
}

#[test]
fn test_csvr_heats_a_cold_gas_to_the_target_temperature() {
    let temperature = 300.0;
    let mut sim = ideal_gas_plain(30, temperature, 0.1, 40.0, 8);
    sim.thermostat = Some(ThermostatType::Csvr);
    sim.tau_t = 0.05;
    assert!(sim.temperature() < 0.11 * temperature);
    let (mean, _) = gas_kinetic_statistics(&mut sim, 0.01, 1000, 100_000);
    let mean_temperature = 2.0 * mean / (sim.degrees_of_freedom() as f64 * KB_REDUCED);
    assert!(
        ((mean_temperature - temperature) / temperature).abs() < 0.03,
        "T = {mean_temperature} K"
    );
}

#[test]
fn test_berendsen_thermostat_suppresses_the_kinetic_energy_fluctuations() {
    // Weak-coupling rescaling drives K deterministically to the target; the
    // canonical fluctuations that CSVR reproduces are absent.  (This is the
    // behaviour that used to be shared by the Csvr and NoseHoover settings.)
    let n = 20_usize;
    let temperature = 250.0;
    let mut sim = ideal_gas_plain(n, temperature, 0.5, 40.0, 3);
    sim.thermostat = Some(ThermostatType::Berendsen);
    sim.tau_t = 0.1;
    let n_dof = 3.0 * n as f64;
    let target = 0.5 * n_dof * KB_REDUCED * temperature;
    let (mean, variance) = gas_kinetic_statistics(&mut sim, 0.01, 3000, 20_000);
    let canonical_variance = 2.0 * target * target / n_dof;
    assert!(((mean - target) / target).abs() < 1e-6, "<K> = {mean}");
    assert!(
        variance < 1e-6 * canonical_variance,
        "var K = {variance} vs canonical {canonical_variance}"
    );
}

#[test]
fn test_csvr_runs_are_reproducible_and_depend_on_the_seed() {
    let run = |seed: u64| {
        let mut sim = ideal_gas_plain(10, 200.0, 1.0, 30.0, 5);
        sim.thermostat = Some(ThermostatType::Csvr);
        sim.tau_t = 0.05;
        sim.set_seed(seed);
        sim.run_with(200, 0.01, free_gas).expect("valid run");
        sim.kinetic_energy()
    };
    assert_eq!(run(11), run(11));
    assert_ne!(run(11), run(12));
}

#[test]
fn test_nose_hoover_chain_conserves_the_extended_energy() {
    // H' = K + U + H_bath is the invariant of the Nosé-Hoover chain equations of
    // motion; it is conserved by a correct reversible integrator (up to the
    // O(dt²) error of velocity Verlet) while K and U separately fluctuate.
    let temperature = 120.0;
    let (mut sim, lj) = lj_fluid_plain(0.75, temperature, 23);
    sim.thermostat = Some(ThermostatType::NoseHoover);
    sim.tau_t = 0.3;
    let mut evaluate = lj.plain_fn();
    sim.refresh_forces(&mut evaluate).expect("force evaluation");
    let h0 = sim.extended_energy().expect("potential energy supplied");
    let mut worst: f64 = 0.0;
    let mut kinetic_min = f64::INFINITY;
    let mut kinetic_max: f64 = 0.0;
    for _ in 0..3000 {
        sim.step_with(0.004, &mut evaluate).expect("valid step");
        let h = sim.extended_energy().expect("potential energy supplied");
        worst = worst.max((h - h0).abs());
        let kinetic = sim.kinetic_energy();
        kinetic_min = kinetic_min.min(kinetic);
        kinetic_max = kinetic_max.max(kinetic);
    }
    let kinetic_scale = 0.5 * sim.degrees_of_freedom() as f64 * KB_REDUCED * temperature;
    assert!(
        worst < 2.0e-3 * kinetic_scale,
        "extended energy drifted by {worst} kJ/mol (kinetic scale {kinetic_scale})"
    );
    // The thermostat really acts: the kinetic energy moves by far more than H'.
    assert!(
        kinetic_max - kinetic_min > 20.0 * worst,
        "K varied by {} kJ/mol, H' by {worst}",
        kinetic_max - kinetic_min
    );
    // The chain holds state between steps.
    assert!(sim.nose_hoover.eta_dot.iter().any(|&v| v != 0.0));
}

#[test]
fn test_nose_hoover_chain_equilibrates_to_the_target_temperature() {
    let temperature = 150.0;
    let (mut sim, lj) = lj_fluid_plain(0.7, 60.0, 29);
    sim.thermostat = Some(ThermostatType::NoseHoover);
    sim.target_temp = temperature;
    sim.tau_t = 0.2;
    let mut evaluate = lj.plain_fn();
    sim.run_with(2500, 0.004, &mut evaluate).expect("warm-up");
    let records = sim.run_with(5000, 0.004, &mut evaluate).expect("sampling");
    let mean_t = records.iter().map(|r| r.1).sum::<f64>() / records.len() as f64;
    assert!(
        ((mean_t - temperature) / temperature).abs() < 0.03,
        "<T> = {mean_t} K, target {temperature} K"
    );
}

#[test]
fn test_nose_hoover_chain_kinetic_energy_fluctuations_are_canonical() {
    // In the canonical ensemble the kinetic energy of N_f degrees of freedom has the
    // mean N_f kT / 2 and the variance (2 / N_f) <K>², whatever the interaction.  The
    // chain is deterministic, so the fluctuations come from the energy exchange with
    // the interacting fluid and test the integrator and the chain masses together.
    let temperature = 150.0;
    let (mut sim, lj) = lj_fluid_plain(0.7, 100.0, 37);
    sim.thermostat = Some(ThermostatType::NoseHoover);
    sim.target_temp = temperature;
    sim.tau_t = 0.2;
    let mut evaluate = lj.plain_fn();
    sim.run_with(3000, 0.004, &mut evaluate).expect("warm-up");
    let steps = 20_000_u64;
    let (mut sum, mut sum_sq) = (0.0, 0.0);
    for _ in 0..steps {
        sim.step_with(0.004, &mut evaluate).expect("sampling step");
        let kinetic = sim.kinetic_energy();
        sum += kinetic;
        sum_sq += kinetic * kinetic;
    }
    let mean = sum / steps as f64;
    let variance = sum_sq / steps as f64 - mean * mean;
    let n_dof = sim.degrees_of_freedom() as f64;
    let target = 0.5 * n_dof * KB_REDUCED * temperature;
    let canonical_variance = 2.0 * target * target / n_dof;
    assert!(
        ((mean - target) / target).abs() < 0.03,
        "<K> = {mean} kJ/mol, expected {target}"
    );
    assert!(
        ((variance - canonical_variance) / canonical_variance).abs() < 0.30,
        "var K = {variance}, expected (2/N_f) <K>² = {canonical_variance}"
    );
}

// ---------------------------------------------------------------------------
// Pressure coupling
// ---------------------------------------------------------------------------

/// Volume at which `n` ideal-gas atoms at `temperature` (K) exert `pressure_bar`.
fn ideal_gas_volume(n: usize, temperature: f64, pressure_bar: f64) -> f64 {
    n as f64 * KB_REDUCED * temperature / (pressure_bar / KJ_MOL_PER_ANG3_TO_BAR)
}

#[test]
fn test_pressure_from_the_stored_virial_follows_the_virial_theorem() {
    let (mut sim, lj) = lj_fluid_plain(0.8, 120.0, 5);
    assert_eq!(sim.pressure(), None, "nothing evaluated yet");
    assert_eq!(sim.total_energy(), None);
    sim.refresh_forces(lj.plain_fn()).expect("force evaluation");
    let positions: Vec<[f64; 3]> = sim.atoms.iter().map(|atom| atom.position).collect();
    let result = lj.evaluate(&positions, &sim.box_lengths);
    let expected = (2.0 * sim.kinetic_energy() + result.virial) / (3.0 * sim.volume())
        * KJ_MOL_PER_ANG3_TO_BAR;
    let pressure = sim.pressure().expect("virial attached");
    assert!(
        (pressure - expected).abs() < 1e-9 * expected.abs().max(1.0),
        "P = {pressure} bar, expected {expected} bar"
    );
    assert_eq!(sim.potential_energy, Some(result.potential_energy));
    assert_eq!(sim.virial, Some(result.virial));
}

#[test]
fn test_berendsen_barostat_sees_the_pair_virial() {
    // Atoms at rest: P = W / (3 V) · 16605.4 bar.  The virial of the force routine
    // alone dilates the box by mu = (1 + beta dt / tau_p (P - P0))^(1/3).
    let edge = 10.0;
    let virial = 5000.0;
    let dt = 0.002;
    let atoms = vec![
        Atom::new([2.0, 3.0, 4.0], [0.0; 3], ARGON_MASS, 0.0),
        Atom::new([6.0, 5.0, 7.0], [0.0; 3], ARGON_MASS, 0.0),
    ];
    let mut sim = MdSimulationPlain::new(atoms, [edge; 3]);
    sim.barostat = Some(BarostatType::Berendsen);
    sim.tau_p = 0.1;
    sim.target_pressure = 1.0;
    let routine = move |atoms: &[Atom], _box: &[f64; 3]| {
        ForceEvaluation::new(vec![[0.0; 3]; atoms.len()]).with_virial(virial)
    };
    sim.step_with(dt, routine).expect("valid step");
    let pressure = virial / (3.0 * edge.powi(3)) * KJ_MOL_PER_ANG3_TO_BAR;
    let mu = (1.0 + sim.compressibility * dt / sim.tau_p * (pressure - 1.0)).cbrt();
    assert!(mu > 1.0 + 1e-3, "visible dilation needed: {mu}");
    for edge_now in sim.box_lengths {
        assert!(
            (edge_now - edge * mu).abs() < 1e-12,
            "edge {edge_now}, expected {}",
            edge * mu
        );
    }
    // The coordinates are dilated with the box.
    assert!((sim.atoms[0].position[0] - 2.0 * mu).abs() < 1e-12);
}

#[test]
fn test_parrinello_rahman_box_velocity_is_persistent_state_of_the_simulation() {
    let edge = 10.0;
    let virial = 5000.0;
    let dt = 0.002;
    let atoms = vec![
        Atom::new([2.0, 3.0, 4.0], [0.0; 3], ARGON_MASS, 0.0),
        Atom::new([6.0, 5.0, 7.0], [0.0; 3], ARGON_MASS, 0.0),
    ];
    let mut sim = MdSimulationPlain::new(atoms, [edge; 3]);
    sim.barostat = Some(BarostatType::Parrinello);
    sim.tau_p = 0.5;
    sim.target_pressure = 1.0;
    let pushing = move |atoms: &[Atom], _box: &[f64; 3]| {
        ForceEvaluation::new(vec![[0.0; 3]; atoms.len()]).with_virial(virial)
    };
    sim.step_with(dt, pushing).expect("valid step");

    let pressure = virial / (3.0 * edge.powi(3)) * KJ_MOL_PER_ANG3_TO_BAR;
    let reference = ParrinelloRahmanBarostat::new(1.0, 0.5, sim.compressibility);
    let acceleration = reference.box_acceleration(pressure, &[edge; 3]);
    for (axis, &a) in acceleration.iter().enumerate() {
        assert!(a > 0.0);
        assert!(
            (sim.box_velocity[axis] - dt * a).abs() < 1e-15,
            "box velocity {}, expected {}",
            sim.box_velocity[axis],
            dt * a
        );
        assert!(
            (sim.box_lengths[axis] - (edge + dt * dt * a)).abs() < 1e-14,
            "edge {}",
            sim.box_lengths[axis]
        );
    }
    // From now on the force routine reports the target pressure (the atoms stay at
    // rest, so W = 3 V P0 gives P = P0).  The step that follows is still driven by
    // the virial stored by the previous evaluation; the one after it feels no force
    // and the box coasts with its velocity.
    let at_target = |atoms: &[Atom], box_lengths: &[f64; 3]| {
        let volume = box_lengths[0] * box_lengths[1] * box_lengths[2];
        let virial = 3.0 * volume / KJ_MOL_PER_ANG3_TO_BAR;
        ForceEvaluation::new(vec![[0.0; 3]; atoms.len()]).with_virial(virial)
    };
    sim.step_with(dt, at_target).expect("valid step");
    let velocity = sim.box_velocity;
    let before = sim.box_lengths;
    sim.step_with(dt, at_target).expect("valid step");
    for axis in 0..3 {
        assert!(
            (sim.box_velocity[axis] - velocity[axis]).abs() < 1e-12,
            "no force at P = P0: velocity persists ({} vs {})",
            sim.box_velocity[axis],
            velocity[axis]
        );
        assert!(
            (sim.box_lengths[axis] - (before[axis] + dt * velocity[axis])).abs() < 1e-12,
            "the box coasts with its velocity"
        );
        assert!(sim.box_lengths[axis] > before[axis]);
    }
}

#[test]
fn test_ideal_gas_volume_under_berendsen_barostat_equals_nkt_over_p() {
    let n = 64_usize;
    let temperature = 300.0;
    let target_bar = 1000.0;
    let equilibrium = ideal_gas_volume(n, temperature, target_bar);
    let edge = (2.2 * equilibrium).cbrt();
    let mut sim = ideal_gas_plain(n, temperature, 1.0, edge, 41);
    sim.thermostat = Some(ThermostatType::Berendsen);
    sim.tau_t = 0.2;
    sim.barostat = Some(BarostatType::Berendsen);
    sim.tau_p = 2.0;
    sim.target_pressure = target_bar;
    sim.compressibility = 1.0 / target_bar;
    sim.run_with(15_000, 0.004, free_gas).expect("valid run");
    let volume = sim.volume();
    // The thermostat holds K at its target and the weak coupling relaxes V to the
    // fixed point P = N k_B T / V of the ideal gas.
    assert!(
        ((volume - equilibrium) / equilibrium).abs() < 0.01,
        "V = {volume} Å³, N k_B T / P = {equilibrium} Å³"
    );
}

#[test]
fn test_ideal_gas_volume_under_parrinello_rahman_equals_nkt_over_p() {
    let n = 64_usize;
    let temperature = 300.0;
    let target_bar = 1000.0;
    let equilibrium = ideal_gas_volume(n, temperature, target_bar);
    let edge = (2.2 * equilibrium).cbrt();
    let mut sim = ideal_gas_plain(n, temperature, 1.0, edge, 41);
    sim.thermostat = Some(ThermostatType::Berendsen);
    sim.tau_t = 0.2;
    sim.barostat = Some(BarostatType::Parrinello);
    sim.tau_p = 1.0;
    sim.target_pressure = target_bar;
    sim.compressibility = 1.0 / target_bar;
    let steps = 25_000_u64;
    let mut volume_sum = 0.0;
    for step in 0..steps {
        sim.step_with(0.004, free_gas).expect("valid step");
        if step >= steps / 2 {
            volume_sum += sim.volume();
        }
    }
    let mean_volume = volume_sum / (steps - steps / 2) as f64;
    assert!(
        ((mean_volume - equilibrium) / equilibrium).abs() < 0.01,
        "mean V = {mean_volume} Å³, N k_B T / P = {equilibrium} Å³"
    );
}

/// Chunked NPT trajectory of the 108-atom Lennard-Jones fluid; returns the
/// `(reduced density, pressure in bar, temperature)` after every chunk.
fn lj_npt_series(
    start_density: f64,
    barostat: BarostatType,
    chunks: usize,
) -> Vec<(f64, f64, f64)> {
    let temperature = 150.0;
    let target_bar = 500.0;
    let (mut sim, lj) = lj_fluid_plain(start_density, temperature, 31);
    sim.thermostat = Some(ThermostatType::Berendsen);
    sim.tau_t = 0.5;
    sim.barostat = Some(barostat);
    sim.tau_p = if barostat == BarostatType::Parrinello {
        1.0
    } else {
        2.0
    };
    sim.target_pressure = target_bar;
    sim.compressibility = 1.0e-4;
    let mut evaluate = lj.plain_fn();
    let mut series = Vec::with_capacity(chunks);
    for _ in 0..chunks {
        sim.run_with(100, 0.004, &mut evaluate).expect("valid run");
        series.push((
            reduced_density(sim.n_atoms(), &sim.box_lengths),
            sim.pressure().expect("virial supplied"),
            sim.temperature(),
        ));
    }
    series
}

/// Run the LJ fluid at 150 K and 500 bar under `barostat` from a compressed
/// (`ρ* = 0.85`) and an expanded (`ρ* = 0.55`) start; both must relax to the target
/// pressure within `pressure_tolerance` and to the same density.
fn check_lj_fluid_relaxes_to_the_target_pressure(barostat: BarostatType, pressure_tolerance: f64) {
    let mut finals = Vec::new();
    for start in [0.85, 0.55] {
        let series = lj_npt_series(start, barostat, 200);
        let tail = &series[series.len() / 2..];
        let mean = |pick: fn(&(f64, f64, f64)) -> f64| {
            tail.iter().map(pick).sum::<f64>() / tail.len() as f64
        };
        let mean_p = mean(|s| s.1);
        let mean_t = mean(|s| s.2);
        assert!(
            ((mean_p - 500.0) / 500.0).abs() < pressure_tolerance,
            "{barostat:?}, start {start}: mean P = {mean_p} bar, target 500 bar"
        );
        assert!(
            ((mean_t - 150.0) / 150.0).abs() < 0.03,
            "{barostat:?}, start {start}: T = {mean_t} K"
        );
        finals.push(mean(|s| s.0));
    }
    assert!(
        ((finals[0] - finals[1]) / finals[1]).abs() < 0.05,
        "{barostat:?}: the two starting densities must relax to the same state: {finals:?}"
    );
}

#[test]
fn test_lj_fluid_reaches_the_target_pressure_from_both_sides_under_berendsen_barostat() {
    // Weak coupling leaves a residual pressure offset (finite relaxation time).
    check_lj_fluid_relaxes_to_the_target_pressure(BarostatType::Berendsen, 0.10);
}

#[test]
fn test_lj_fluid_reaches_the_target_pressure_from_both_sides_under_parrinello_rahman() {
    // The Parrinello-Rahman box integrates the pressure error and averages to the
    // target.
    check_lj_fluid_relaxes_to_the_target_pressure(BarostatType::Parrinello, 0.05);
}

#[test]
fn test_monte_carlo_barostat_samples_the_isothermal_isobaric_volume_distribution() {
    // For an ideal gas the NPT volume density is V^N exp(-P V / kT), i.e.
    // V ~ Gamma(N + 1, kT / P): <V> = (N + 1) kT / P and var V = (N + 1) (kT / P)².
    // Both the P dV term and the N ln(V'/V) term of the acceptance rule enter.
    let n = 20_usize;
    let temperature = 300.0;
    let pressure_bar = 500.0;
    let kt_over_p = KB_REDUCED * temperature / (pressure_bar / KJ_MOL_PER_ANG3_TO_BAR);
    let mut sim = ideal_gas_plain(n, temperature, 1.0, 12.0, 77);
    sim.barostat = Some(BarostatType::MonteCarlo);
    sim.target_pressure = pressure_bar;
    sim.monte_carlo = McVolumeMoves {
        max_relative_change: 0.3,
        interval: 1,
        attempted: 0,
        accepted: 0,
    };
    let moves = 300_000_u64;
    sim.run_with(2000, 0.002, free_gas).expect("warm-up");
    let (mut sum, mut sum_sq) = (0.0, 0.0);
    for _ in 0..moves {
        sim.step_with(0.002, free_gas).expect("valid step");
        let volume = sim.volume();
        sum += volume;
        sum_sq += volume * volume;
    }
    let mean = sum / moves as f64;
    let variance = sum_sq / moves as f64 - mean * mean;
    let shape = n as f64 + 1.0;
    assert!(
        ((mean - shape * kt_over_p) / (shape * kt_over_p)).abs() < 0.02,
        "<V> = {mean} Å³, expected {}",
        shape * kt_over_p
    );
    let expected_variance = shape * kt_over_p * kt_over_p;
    assert!(
        ((variance - expected_variance) / expected_variance).abs() < 0.08,
        "var V = {variance}, expected {expected_variance}"
    );
    let rate = sim.monte_carlo.acceptance_rate();
    assert!(rate > 0.3 && rate < 0.95, "acceptance rate {rate}");
    assert_eq!(sim.monte_carlo.attempted, moves + 2000);
}

#[test]
fn test_monte_carlo_barostat_uses_the_potential_energy_change() {
    // A potential energy U = c V (forces zero) acts like the extra pressure
    // c: V ~ Gamma(N + 1, kT / (P + c)).
    let n = 20_usize;
    let temperature = 300.0;
    let pressure_bar = 500.0;
    let c = 0.02; // kJ mol⁻¹ Å⁻³
    let kt = KB_REDUCED * temperature;
    let effective_pressure = pressure_bar / KJ_MOL_PER_ANG3_TO_BAR + c;
    let mut sim = ideal_gas_plain(n, temperature, 1.0, 12.0, 78);
    sim.barostat = Some(BarostatType::MonteCarlo);
    sim.target_pressure = pressure_bar;
    sim.monte_carlo = McVolumeMoves {
        max_relative_change: 0.3,
        interval: 1,
        attempted: 0,
        accepted: 0,
    };
    let linear = move |atoms: &[Atom], box_lengths: &[f64; 3]| {
        let volume = box_lengths[0] * box_lengths[1] * box_lengths[2];
        ForceEvaluation::new(vec![[0.0; 3]; atoms.len()])
            .with_virial(0.0)
            .with_potential_energy(c * volume)
    };
    sim.run_with(2000, 0.002, linear).expect("warm-up");
    let moves = 300_000_u64;
    let mut sum = 0.0;
    for _ in 0..moves {
        sim.step_with(0.002, linear).expect("valid step");
        sum += sim.volume();
    }
    let mean = sum / moves as f64;
    let expected = (n as f64 + 1.0) * kt / effective_pressure;
    assert!(
        ((mean - expected) / expected).abs() < 0.02,
        "<V> = {mean} Å³, expected {expected} Å³"
    );
}

#[test]
fn test_monte_carlo_barostat_keeps_the_forces_of_the_accepted_configuration() {
    // After a rejected move the box, the coordinates and the stored forces refer to
    // the old configuration; after an accepted one to the trial configuration.
    let (mut sim, lj) = lj_fluid_plain(0.7, 100.0, 9);
    sim.barostat = Some(BarostatType::MonteCarlo);
    sim.target_pressure = 500.0;
    sim.monte_carlo.interval = 1;
    sim.monte_carlo.max_relative_change = 0.005;
    let mut evaluate = lj.plain_fn();
    for _ in 0..400 {
        sim.step_with(0.002, &mut evaluate).expect("valid step");
        let positions: Vec<[f64; 3]> = sim.atoms.iter().map(|atom| atom.position).collect();
        let fresh = lj.evaluate(&positions, &sim.box_lengths);
        for (atom, force) in sim.atoms.iter().zip(&fresh.forces) {
            for (stored, current) in atom.force.iter().zip(force.iter()) {
                assert!(
                    (stored - current).abs() < 1e-9 * current.abs().max(1.0),
                    "stored force {stored} differs from the force of the current configuration {current}"
                );
            }
        }
        assert_eq!(sim.virial, Some(fresh.virial));
        assert_eq!(sim.potential_energy, Some(fresh.potential_energy));
    }
    assert!(sim.monte_carlo.accepted > 0, "no move was accepted");
    assert!(
        sim.monte_carlo.accepted < sim.monte_carlo.attempted,
        "every move was accepted"
    );
}

// ---------------------------------------------------------------------------
// Configuration errors
// ---------------------------------------------------------------------------

/// Message of the error returned by `result`, panicking if it is `Ok`.
fn error_message(result: crate::error::Result<()>) -> String {
    match result {
        Ok(()) => panic!("an error was expected"),
        Err(error) => error.to_string(),
    }
}

fn two_atoms() -> MdSimulationPlain {
    let atoms = vec![
        Atom::new([1.0, 1.0, 1.0], [0.5, 0.0, 0.0], 12.0, 0.0),
        Atom::new([4.0, 4.0, 4.0], [0.0, 0.5, 0.0], 12.0, 0.0),
    ];
    MdSimulationPlain::new(atoms, [10.0; 3])
}

#[test]
fn test_invalid_time_steps_are_rejected_without_changing_the_state() {
    for dt in [0.0, -0.001, f64::NAN, f64::INFINITY] {
        let mut sim = two_atoms();
        let message = error_message(sim.step(dt));
        assert!(message.contains("time step"), "{message}");
        assert_eq!(sim.step_count, 0);
        assert_eq!(sim.atoms[0].position, [1.0, 1.0, 1.0]);
        assert!(sim.step_with(dt, free_gas).is_err());
    }
}

#[test]
fn test_invalid_masses_and_boxes_are_rejected() {
    let mut sim = two_atoms();
    sim.atoms[1].mass = 0.0;
    assert!(error_message(sim.step(0.001)).contains("atom 1"));
    sim.atoms[1].mass = f64::NAN;
    assert!(sim.step(0.001).is_err());

    let mut sim = two_atoms();
    sim.box_lengths = [10.0, -1.0, 10.0];
    assert!(error_message(sim.step(0.001)).contains("box lengths"));
    sim.box_lengths = [10.0, f64::INFINITY, 10.0];
    assert!(sim.step(0.001).is_err());
}

#[test]
fn test_invalid_thermostat_settings_are_rejected() {
    for thermostat in [
        ThermostatType::Berendsen,
        ThermostatType::Csvr,
        ThermostatType::NoseHoover,
    ] {
        let mut sim = two_atoms();
        sim.thermostat = Some(thermostat);
        sim.tau_t = 0.0;
        assert!(
            error_message(sim.step(0.001)).contains("tau_t"),
            "{thermostat:?}"
        );
        sim.tau_t = 0.1;
        sim.target_temp = -5.0;
        assert!(sim.step(0.001).is_err(), "{thermostat:?}");
        sim.target_temp = f64::NAN;
        assert!(sim.step(0.001).is_err(), "{thermostat:?}");
    }
    // A chain cannot hold a system at absolute zero.
    let mut sim = two_atoms();
    sim.thermostat = Some(ThermostatType::NoseHoover);
    sim.target_temp = 0.0;
    assert!(sim.step(0.001).is_err());
    // Rescaling thermostats can.
    sim.thermostat = Some(ThermostatType::Csvr);
    assert!(sim.step(0.001).is_ok());
}

#[test]
fn test_invalid_barostat_settings_are_rejected() {
    let mut sim = two_atoms();
    sim.barostat = Some(BarostatType::Berendsen);
    sim.tau_p = -1.0;
    assert!(error_message(sim.step(0.001)).contains("tau_p"));
    sim.tau_p = 1.0;
    sim.compressibility = 0.0;
    assert!(error_message(sim.step(0.001)).contains("compressibility"));
    sim.compressibility = 4.5e-5;
    sim.target_pressure = f64::NAN;
    assert!(error_message(sim.step(0.001)).contains("target pressure"));

    let mut sim = two_atoms();
    sim.barostat = Some(BarostatType::MonteCarlo);
    for change in [0.0, 1.0, -0.1, f64::NAN] {
        sim.monte_carlo.max_relative_change = change;
        assert!(
            error_message(sim.step(0.001)).contains("max_relative_change"),
            "{change}"
        );
    }
    sim.monte_carlo.max_relative_change = 0.01;
    sim.target_temp = 0.0;
    assert!(error_message(sim.step(0.001)).contains("temperature"));
}

#[test]
fn test_pressure_coupling_without_a_virial_is_an_error() {
    let forces_only =
        |atoms: &[Atom], _box: &[f64; 3]| ForceEvaluation::new(vec![[0.0; 3]; atoms.len()]);
    for barostat in [BarostatType::Berendsen, BarostatType::Parrinello] {
        let mut sim = two_atoms();
        sim.barostat = Some(barostat);
        let message = error_message(sim.step_with(0.001, forces_only));
        assert!(message.contains("virial"), "{barostat:?}: {message}");
    }
    // The frozen-force step knows its virial (zero), so it is accepted.
    let mut sim = two_atoms();
    sim.barostat = Some(BarostatType::Berendsen);
    assert!(sim.step(0.001).is_ok());
}

#[test]
fn test_monte_carlo_barostat_without_potential_energy_is_an_error() {
    let no_energy = |atoms: &[Atom], _box: &[f64; 3]| {
        ForceEvaluation::new(vec![[0.0; 3]; atoms.len()]).with_virial(0.0)
    };
    let mut sim = two_atoms();
    sim.barostat = Some(BarostatType::MonteCarlo);
    let message = error_message(sim.step_with(0.001, no_energy));
    assert!(message.contains("potential energy"), "{message}");
}

#[test]
fn test_wrong_number_of_forces_is_an_error() {
    let too_few = |_atoms: &[Atom], _box: &[f64; 3]| ForceEvaluation::new(vec![[0.0; 3]]);
    let mut sim = two_atoms();
    let message = error_message(sim.step_with(0.001, too_few));
    assert!(message.contains("1 forces for 2 atoms"), "{message}");
    assert!(sim.refresh_forces(too_few).is_err());
}

// ---------------------------------------------------------------------------
// Small API pieces
// ---------------------------------------------------------------------------

#[test]
fn test_force_evaluation_builders_and_counters() {
    let evaluation = ForceEvaluation::new(vec![[1.0, 2.0, 3.0]]);
    assert_eq!(evaluation.virial, None);
    assert_eq!(evaluation.potential_energy, None);
    let full = evaluation
        .clone()
        .with_virial(-2.5)
        .with_potential_energy(7.0);
    assert_eq!(full.virial, Some(-2.5));
    assert_eq!(full.potential_energy, Some(7.0));
    assert_eq!(full.forces, vec![[1.0, 2.0, 3.0]]);
    assert_eq!(ForceEvaluation::default().forces.len(), 0);

    let mut moves = McVolumeMoves::default();
    assert_eq!(moves.acceptance_rate(), 0.0);
    moves.attempted = 8;
    moves.accepted = 2;
    assert_eq!(moves.acceptance_rate(), 0.25);

    let sim = two_atoms();
    assert_eq!(sim.degrees_of_freedom(), 6);
    assert_eq!(sim.volume(), 1000.0);
    assert_eq!(sim.extended_energy(), None);
}

#[test]
fn test_extended_energy_equals_the_total_energy_without_a_chain() {
    let (mut sim, lj) = lj_fluid_plain(0.7, 100.0, 4);
    sim.thermostat = Some(ThermostatType::Csvr);
    sim.refresh_forces(lj.plain_fn()).expect("force evaluation");
    assert_eq!(sim.extended_energy(), sim.total_energy());
    sim.thermostat = Some(ThermostatType::NoseHoover);
    // The chain is still at rest before the first step.
    assert_eq!(sim.extended_energy(), sim.total_energy());
}
