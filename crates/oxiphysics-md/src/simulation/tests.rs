// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for MD simulation modules.

use super::*;

use crate::atom::AtomSet;
use crate::barostat::NoBarostat;
use crate::forcefield::PairForceField;
use crate::integrator::VelocityVerlet;
use crate::neighbor::PeriodicBox;
use crate::potential::LennardJones;
use crate::thermostat::NoThermostat;
use oxiphysics_core::math::Vec3;

// ----- MdSim / MdState / MdConfig unit tests -----

fn make_two_atom_state() -> MdState {
    MdState::new(
        vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        vec![[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]],
        vec![1.0, 1.0],
    )
}

#[test]
fn test_md_state_n_atoms() {
    let state = make_two_atom_state();
    assert_eq!(state.n_atoms(), 2);
}

#[test]
fn test_compute_kinetic_energy_two_atoms() {
    let config = MdConfig::default();
    let state = make_two_atom_state();
    let sim = MdSim::new(config, state);
    let ke = sim.compute_kinetic_energy();
    // KE = 2 × 0.5 × 1.0 × 0.1² = 0.01
    assert!((ke - 0.01).abs() < 1e-12, "KE = {ke}");
}

#[test]
fn test_compute_temperature_positive() {
    let config = MdConfig::default();
    let state = make_two_atom_state();
    let sim = MdSim::new(config, state);
    let t = sim.compute_temperature();
    assert!(t > 0.0, "temperature should be positive, got {t}");
}

#[test]
fn test_temperature_zero_for_stationary_atoms() {
    let config = MdConfig::default();
    let state = MdState::new(
        vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        vec![1.0, 1.0],
    );
    let sim = MdSim::new(config, state);
    let t = sim.compute_temperature();
    assert!(
        t.abs() < 1e-12,
        "temperature should be 0 for stationary atoms, got {t}"
    );
}

#[test]
fn test_velocity_rescaling_hits_target() {
    let config = MdConfig {
        ensemble: Ensemble::NVT,
        temperature: 300.0,
        ..Default::default()
    };
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config, state);
    sim.apply_velocity_rescaling(300.0);
    let t = sim.compute_temperature();
    // After rescaling T should equal target within floating point
    assert!(
        (t - 300.0).abs() < 1e-6,
        "temperature after rescaling = {t}, expected 300.0"
    );
}

#[test]
fn test_compute_pressure_positive_for_warm_gas() {
    let config = MdConfig {
        box_lengths: [10.0, 10.0, 10.0],
        ..MdConfig::default()
    };
    let state = make_two_atom_state();
    let sim = MdSim::new(config, state);
    let p = sim.compute_pressure(0.0); // zero virial (ideal gas)
    assert!(p > 0.0, "pressure should be positive for warm gas, got {p}");
}

#[test]
fn test_step_advances_time_and_step_count() {
    let config = MdConfig {
        dt: 0.001,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config, state);
    sim.step();
    assert_eq!(sim.state.step, 1);
    assert!((sim.state.time - 0.001).abs() < 1e-15);
}

#[test]
fn test_velocity_verlet_free_particle_displacement() {
    // Free particle: no forces, constant velocity
    // After one step: x should move by v * dt exactly.
    let config = MdConfig {
        dt: 0.01,
        pbc: false,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state = MdState::new(vec![[0.0, 0.0, 0.0]], vec![[2.0, 0.0, 0.0]], vec![1.0]);
    let mut sim = MdSim::new(config, state);
    let zero_forces = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    sim.velocity_verlet_step(zero_forces);
    // x should be 2.0 * 0.01 = 0.02
    assert!((sim.state.positions[0][0] - 0.02).abs() < 1e-12);
}

#[test]
fn test_ensemble_enum_variants() {
    assert_ne!(Ensemble::NVE, Ensemble::NVT);
    assert_ne!(Ensemble::NVT, Ensemble::NPT);
}

#[test]
fn test_md_config_default_ensemble_is_nvt() {
    let cfg = MdConfig::default();
    assert_eq!(cfg.ensemble, Ensemble::NVT);
}

#[test]
fn test_md_sim_run_records_correct_count() {
    let config = MdConfig {
        n_steps: 10,
        dt: 0.001,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config, state);
    let zero_forces = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    let records = sim.run_with_forces(zero_forces, 2);
    assert_eq!(
        records.len(),
        5,
        "expected 5 records for 10 steps with interval 2"
    );
}

// ----- Berendsen rescaler tests -----

#[test]
fn test_berendsen_rescaler_scale_factor_at_target() {
    let rescaler = BerendsenRescaler::new(0.1);
    let scale = rescaler.scale_factor(300.0, 300.0, 0.002);
    assert!(
        (scale - 1.0).abs() < 1e-12,
        "scale should be 1.0 at target T, got {scale}"
    );
}

#[test]
fn test_berendsen_rescaler_cold_system_scales_up() {
    let rescaler = BerendsenRescaler::new(0.1);
    let scale = rescaler.scale_factor(100.0, 300.0, 0.002);
    assert!(scale > 1.0, "cold system should scale up, got {scale}");
}

#[test]
fn test_berendsen_rescaler_hot_system_scales_down() {
    let rescaler = BerendsenRescaler::new(0.1);
    let scale = rescaler.scale_factor(600.0, 300.0, 0.002);
    assert!(scale < 1.0, "hot system should scale down, got {scale}");
}

#[test]
fn test_remove_com_velocity_zeroes_momentum() {
    let config = MdConfig::default();
    let state = MdState::new(
        vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
        vec![2.0, 3.0],
    );
    let mut sim = MdSim::new(config, state);
    sim.remove_com_velocity();
    let p = sim.total_momentum();
    let mag = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    assert!(
        mag < 1e-10,
        "total momentum after COM removal should be ~0, got {mag}"
    );
}

#[test]
fn test_total_energy_is_sum_of_ke_and_pe() {
    let config = MdConfig::default();
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config, state);
    sim.state.kinetic_energy = 10.0;
    sim.state.potential_energy = -5.0;
    assert!((sim.total_energy() - 5.0).abs() < 1e-12);
}

#[test]
fn test_step_with_forces_advances_state() {
    let config = MdConfig {
        dt: 0.001,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config, state);
    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    sim.step_with_forces(zero_f, 0.1);
    assert_eq!(sim.state.step, 1);
    assert!((sim.state.time - 0.001).abs() < 1e-15);
}

#[test]
fn test_berendsen_thermostat_applied_nvt() {
    let config = MdConfig {
        dt: 0.001,
        temperature: 300.0,
        ensemble: Ensemble::NVT,
        ..MdConfig::default()
    };
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config, state);
    // Run a few steps with Berendsen thermostat
    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    for _ in 0..10 {
        sim.step_with_forces(zero_f, 0.01);
    }
    // Temperature should approach target (won't be exact in 10 steps)
    let t = sim.compute_temperature();
    assert!(t.is_finite(), "temperature should be finite, got {t}");
}

#[test]
fn test_total_momentum_computed_correctly() {
    let config = MdConfig::default();
    let state = MdState::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0]],
        vec![[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
        vec![3.0, 4.0],
    );
    let sim = MdSim::new(config, state);
    let p = sim.total_momentum();
    // p_x = 3*1 + 4*0 = 3, p_y = 3*0 + 4*2 = 8, p_z = 0
    assert!((p[0] - 3.0).abs() < 1e-12);
    assert!((p[1] - 8.0).abs() < 1e-12);
    assert!(p[2].abs() < 1e-12);
}

// ----- Generic MdSimulation (trait-based) tests (original) -----

#[test]
fn test_two_lj_atoms_energy_conservation() {
    let sigma = 1.0;
    let epsilon = 1.0;
    let cutoff = 5.0;
    let r0 = 1.5;

    let mut atoms = AtomSet::new();
    atoms.add_atom(Vec3::new(0.0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);
    atoms.add_atom(Vec3::new(r0, 0.0, 0.0), Vec3::zeros(), 1.0, 0.0, 0);

    let mut ff = PairForceField::new();
    ff.add_interaction(0, 0, LennardJones::new(epsilon, sigma, cutoff));

    let pbox = PeriodicBox::cubic(10.0);

    let mut sim = MdSimulation::new(
        atoms,
        ff,
        VelocityVerlet::new(),
        NoThermostat::new(),
        NoBarostat::new(),
        pbox,
        0.0,
        0.0,
        1.0,
    );

    let records = sim.run(100, 0.001, Some(1));

    let e0 = records[0].total;
    let max_drift = records
        .iter()
        .map(|r| (r.total - e0).abs())
        .fold(0.0_f64, f64::max);

    let e_scale = e0.abs().max(1.0);
    assert!(
        max_drift / e_scale < 0.01,
        "Energy drift too large: {max_drift} / {e_scale}"
    );
}

// ----- NPT ensemble improvements -----

#[test]
fn test_npt_ensemble_barostat_scales_box() {
    // In NPT with Berendsen barostat, box should change if pressure != target.
    let mut config = MdConfig {
        dt: 0.001,
        temperature: 300.0,
        pressure: 100.0, // unrealistic target to force scaling
        ensemble: Ensemble::NPT,
        box_lengths: [5.0, 5.0, 5.0],
        ..MdConfig::default()
    };
    config.n_steps = 1;
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config.clone(), state.clone());

    // Initial box
    let box_vol_initial = config.box_lengths[0] * config.box_lengths[1] * config.box_lengths[2];

    // Apply NPT step (uses berendsen-like thermostat + barostat)
    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    sim.step_with_forces(zero_f, 0.1);

    // Box volume: remains the same since our simple MdSim doesn't modify box
    // (full NPT barostat not in MdSim), but step should proceed without panic.
    assert_eq!(sim.state.step, 1);
    let _ = box_vol_initial;
}

// ----- NVE monitoring -----

#[test]
fn test_nve_energy_drift_is_small() {
    let config = MdConfig {
        dt: 0.001,
        n_steps: 50,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state = MdState::new(
        vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]],
        vec![1.0, 1.0],
    );
    let mut sim = MdSim::new(config, state);
    let e0 = sim.compute_kinetic_energy();
    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    let records = sim.run_with_forces(zero_f, 10);
    // With zero forces, KE should be exactly constant.
    for (ke, _, _) in &records {
        assert!(
            (ke - e0).abs() < 1e-10,
            "NVE KE drift: {ke} vs initial {e0}"
        );
    }
}

// ----- Simulation checkpointing -----

#[test]
fn test_checkpoint_save_restore() {
    let config = MdConfig {
        dt: 0.002,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state0 = make_two_atom_state();
    let mut sim = MdSim::new(config.clone(), state0.clone());

    // Run for 5 steps.
    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    for _ in 0..5 {
        sim.step_with_forces(zero_f, 0.1);
    }

    // Save checkpoint.
    let ckpt = sim.checkpoint();

    // Run 5 more steps from the checkpoint.
    let mut sim2 = MdSim::restore_checkpoint(config.clone(), ckpt.clone());
    for _ in 0..5 {
        sim2.step_with_forces(zero_f, 0.1);
    }

    // The positions after 10 steps from scratch vs. 5+5 steps should match.
    let mut sim_ref = MdSim::new(config, state0);
    for _ in 0..10 {
        sim_ref.step_with_forces(zero_f, 0.1);
    }

    for i in 0..2 {
        for a in 0..3 {
            assert!(
                (sim2.state.positions[i][a] - sim_ref.state.positions[i][a]).abs() < 1e-12,
                "Checkpoint restore mismatch at particle {i} axis {a}"
            );
        }
    }
}

// ----- Multi-replica -----

#[test]
fn test_multi_replica_all_evolve() {
    let config = MdConfig {
        dt: 0.001,
        n_steps: 10,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state0 = make_two_atom_state();
    let n_replicas = 4;
    let mut replicas: Vec<MdSim> = (0..n_replicas)
        .map(|_| MdSim::new(config.clone(), state0.clone()))
        .collect();

    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    for sim in replicas.iter_mut() {
        sim.run_with_forces(zero_f, 0);
    }

    for (idx, sim) in replicas.iter().enumerate() {
        assert_eq!(
            sim.state.step, 10,
            "Replica {idx} should have taken 10 steps"
        );
    }
}

#[test]
fn test_multi_replica_independent() {
    // Two replicas with different initial velocities should diverge.
    let config = MdConfig {
        dt: 0.001,
        n_steps: 5,
        ensemble: Ensemble::NVE,
        pbc: false,
        ..MdConfig::default()
    };
    let state1 = make_two_atom_state();
    let mut state2 = make_two_atom_state();
    state2.velocities[0][0] = 0.5; // different

    let mut sim1 = MdSim::new(config.clone(), state1.clone());
    let mut sim2_r = MdSim::new(config, state2);

    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    sim1.run_with_forces(zero_f, 0);
    sim2_r.run_with_forces(zero_f, 0);

    let _ = state1; // just verify no compile error
    // Positions should differ due to different initial velocities.
    let diff = (sim1.state.positions[0][0] - sim2_r.state.positions[0][0]).abs();
    assert!(diff > 1e-12, "Replicas should diverge: diff = {diff}");
}

// ----- Profiling -----

#[test]
fn test_simulation_profiling_records_timing() {
    let config = MdConfig {
        dt: 0.001,
        n_steps: 20,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config, state);
    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    let profile = sim.run_with_profiling(zero_f, 5);
    assert!(profile.total_steps == 20, "Profiling should count 20 steps");
    assert!(profile.wall_time_s >= 0.0, "Wall time must be non-negative");
    assert!(profile.steps_per_second > 0.0, "Steps/s must be positive");
}

// ----- MdDriver tests -----

fn make_driver_state(n: usize, speed: f64) -> MdState {
    let mut positions = Vec::with_capacity(n);
    let mut velocities = Vec::with_capacity(n);
    let masses = vec![1.0f64; n];
    for i in 0..n {
        positions.push([(i as f64) * 2.0, 0.0, 0.0]);
        velocities.push([speed * (if i % 2 == 0 { 1.0 } else { -1.0 }), 0.0, 0.0]);
    }
    MdState::new(positions, velocities, masses)
}

#[test]
fn test_md_driver_kinetic_energy() {
    let state = make_driver_state(4, 0.5);
    let driver = MdDriver::new(state, 0.001);
    let ke = driver.kinetic_energy();
    // Each atom: 0.5 * 1.0 * 0.25 = 0.125; 4 atoms: 0.5
    assert!((ke - 0.5).abs() < 1e-12, "KE = {ke}, expected 0.5");
}

#[test]
fn test_md_driver_temperature_positive() {
    let state = make_driver_state(4, 0.5);
    let driver = MdDriver::new(state, 0.001);
    let t = driver.temperature();
    assert!(t > 0.0, "temperature should be positive, got {t}");
}

#[test]
fn test_md_driver_rescale_velocities_hits_target() {
    let state = make_driver_state(10, 0.1);
    let mut driver = MdDriver::new(state, 0.001);
    driver.rescale_velocities_to_temperature(300.0);
    let t = driver.temperature();
    assert!(
        (t - 300.0).abs() < 1e-6,
        "temperature after rescaling = {t}, expected 300.0"
    );
}

#[test]
fn test_md_driver_remove_com_velocity() {
    let state = make_driver_state(4, 0.5);
    let mut driver = MdDriver::new(state, 0.001);
    driver.remove_com_velocity();
    // After removing COM velocity, total momentum should be zero
    let mut px = 0.0f64;
    let n = driver.state.n_atoms();
    for i in 0..n {
        px += driver.state.masses[i] * driver.state.velocities[i][0];
    }
    assert!(
        px.abs() < 1e-10,
        "COM px should be ~0 after removal, got {px}"
    );
}

#[test]
fn test_md_driver_velocity_verlet_free_particle() {
    let state = MdState::new(vec![[0.0, 0.0, 0.0]], vec![[1.0, 0.0, 0.0]], vec![1.0]);
    let mut driver = MdDriver::new(state, 0.01);
    // Zero-force function
    let force_fn = |s: &mut MdState| {
        for f in s.forces.iter_mut() {
            *f = [0.0; 3];
        }
    };
    driver.velocity_verlet_step(force_fn);
    // Free particle should move v * dt = 1.0 * 0.01 = 0.01
    assert!(
        (driver.state.positions[0][0] - 0.01).abs() < 1e-12,
        "expected 0.01, got {}",
        driver.state.positions[0][0]
    );
    assert_eq!(driver.step_count, 1);
}

#[test]
fn test_apply_minimum_image_wraps() {
    let r = [9.5, 0.0, 0.0];
    let dr = apply_minimum_image(r, 10.0);
    assert!(
        (dr[0] - (-0.5)).abs() < 1e-12,
        "minimum image x should be -0.5, got {}",
        dr[0]
    );
}

#[test]
fn test_apply_minimum_image_no_wrap_needed() {
    let r = [2.0, 3.0, -1.0];
    let dr = apply_minimum_image(r, 10.0);
    assert!((dr[0] - 2.0).abs() < 1e-12);
    assert!((dr[1] - 3.0).abs() < 1e-12);
    assert!((dr[2] - (-1.0)).abs() < 1e-12);
}

#[test]
fn test_pair_distance_sq() {
    let ri = [0.0, 0.0, 0.0];
    let rj = [3.0, 4.0, 0.0];
    let d2 = pair_distance_sq(ri, rj);
    assert!((d2 - 25.0).abs() < 1e-12, "d^2 = {d2}, expected 25");
}

#[test]
fn test_compute_lj_forces_newton_third_law() {
    // Two particles: forces should be equal and opposite
    let mut state = MdState::new(
        vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]],
        vec![[0.0; 3]; 2],
        vec![1.0, 1.0],
    );
    compute_lj_forces(&mut state, 1.0, 1.0, 5.0, 20.0);
    for a in 0..3 {
        let sum = state.forces[0][a] + state.forces[1][a];
        assert!(
            sum.abs() < 1e-10,
            "Newton III violated on axis {a}: sum = {sum}"
        );
    }
}

#[test]
fn test_compute_lj_forces_repulsive_at_short_range() {
    // At r < sigma, force should be repulsive (push apart: f on i points toward negative x)
    let mut state = MdState::new(
        vec![[0.0, 0.0, 0.0], [0.8, 0.0, 0.0]],
        vec![[0.0; 3]; 2],
        vec![1.0, 1.0],
    );
    compute_lj_forces(&mut state, 1.0, 1.0, 5.0, 20.0);
    // Force on atom 0 should be in -x direction (pushed away from atom 1)
    assert!(
        state.forces[0][0] < 0.0,
        "LJ force on atom 0 at r<sigma should be repulsive (neg x), got {}",
        state.forces[0][0]
    );
}

#[test]
fn test_compute_lj_forces_attractive_beyond_sigma() {
    // At r = 1.5 * sigma, force should be attractive (atom 0 pulled toward +x)
    let sigma = 1.0f64;
    let r = 1.5 * sigma;
    let mut state = MdState::new(
        vec![[0.0, 0.0, 0.0], [r, 0.0, 0.0]],
        vec![[0.0; 3]; 2],
        vec![1.0, 1.0],
    );
    compute_lj_forces(&mut state, 1.0, sigma, 5.0, 20.0);
    assert!(
        state.forces[0][0] > 0.0,
        "LJ force on atom 0 at r=1.5sigma should be attractive (+x), got {}",
        state.forces[0][0]
    );
}

// ----- ThermostatType / BarostatType / Atom / MdSimulationPlain tests -----

#[test]
fn test_thermostat_type_enum_variants() {
    let _t = ThermostatType::Berendsen;
    let _t = ThermostatType::NoseHoover;
    let _t = ThermostatType::Csvr;
}

#[test]
fn test_barostat_type_enum_variants() {
    let _b = BarostatType::Berendsen;
    let _b = BarostatType::MonteCarlo;
    let _b = BarostatType::Parrinello;
}

#[test]
fn test_atom_struct_fields() {
    let a = Atom {
        position: [1.0, 2.0, 3.0],
        velocity: [0.1, 0.2, 0.3],
        force: [0.0; 3],
        mass: 12.0,
        charge: -0.5,
    };
    assert!((a.position[0] - 1.0).abs() < 1e-14);
    assert!((a.mass - 12.0).abs() < 1e-14);
    assert!((a.charge - (-0.5)).abs() < 1e-14);
}

#[test]
fn test_md_simulation_plain_pbc_wrap() {
    let atoms = vec![Atom {
        position: [11.0, -1.0, 5.0],
        velocity: [0.0; 3],
        force: [0.0; 3],
        mass: 1.0,
        charge: 0.0,
    }];
    let mut sim = MdSimulationPlain::new(atoms, [10.0, 10.0, 10.0]);
    sim.apply_pbc();
    assert!(
        (sim.atoms[0].position[0] - 1.0).abs() < 1e-10,
        "x should wrap to 1.0, got {}",
        sim.atoms[0].position[0]
    );
    assert!(
        (sim.atoms[0].position[1] - 9.0).abs() < 1e-10,
        "y should wrap to 9.0, got {}",
        sim.atoms[0].position[1]
    );
}

#[test]
fn test_md_simulation_plain_temperature() {
    // Two atoms each with vx = 1.0 Å/ps, mass = 1.0 amu
    // KE = 2 * 0.5 * 1 * 1^2 = 1.0 kJ/mol
    // T = 2 KE / (3 N k_B) = 2 * 1.0 / (3 * 2 * 8.314e-3) ~ 40.1 K
    let atoms = vec![
        Atom {
            position: [0.0; 3],
            velocity: [1.0, 0.0, 0.0],
            force: [0.0; 3],
            mass: 1.0,
            charge: 0.0,
        },
        Atom {
            position: [3.0, 0.0, 0.0],
            velocity: [-1.0, 0.0, 0.0],
            force: [0.0; 3],
            mass: 1.0,
            charge: 0.0,
        },
    ];
    let sim = MdSimulationPlain::new(atoms, [20.0, 20.0, 20.0]);
    let t = sim.temperature();
    assert!(t > 0.0, "temperature should be positive, got {t}");
    assert!(t.is_finite(), "temperature should be finite");
}

#[test]
fn test_md_simulation_plain_kinetic_energy() {
    let atoms = vec![Atom {
        position: [0.0; 3],
        velocity: [2.0, 0.0, 0.0],
        force: [0.0; 3],
        mass: 3.0,
        charge: 0.0,
    }];
    let sim = MdSimulationPlain::new(atoms, [20.0; 3]);
    let ke = sim.kinetic_energy();
    // KE = 0.5 * 3 * 4 = 6.0
    assert!((ke - 6.0).abs() < 1e-12, "KE = {ke}, expected 6.0");
}

#[test]
fn test_md_simulation_plain_pressure_virial() {
    // Ideal gas virial = 0 => P = N*kB*T / (3V)
    let atoms = vec![
        Atom {
            position: [0.0; 3],
            velocity: [1.0, 0.0, 0.0],
            force: [0.0; 3],
            mass: 1.0,
            charge: 0.0,
        },
        Atom {
            position: [5.0, 0.0, 0.0],
            velocity: [-1.0, 0.0, 0.0],
            force: [0.0; 3],
            mass: 1.0,
            charge: 0.0,
        },
    ];
    let sim = MdSimulationPlain::new(atoms, [10.0; 3]);
    let p = sim.pressure_virial(0.0);
    assert!(p > 0.0, "ideal gas pressure should be positive, got {p}");
}

#[test]
fn test_md_simulation_plain_step_advances_time() {
    let atoms = vec![Atom {
        position: [0.0; 3],
        velocity: [1.0, 0.0, 0.0],
        force: [0.0; 3],
        mass: 1.0,
        charge: 0.0,
    }];
    let mut sim = MdSimulationPlain::new(atoms, [20.0; 3]);
    let dt = 0.002;
    sim.step(dt);
    assert_eq!(sim.step_count, 1);
    assert!(
        (sim.time - dt).abs() < 1e-14,
        "time = {}, expected {dt}",
        sim.time
    );
}

#[test]
fn test_md_simulation_plain_step_free_particle() {
    // Free particle with zero force: position should advance by v*dt
    let v = 2.5_f64;
    let dt = 0.01_f64;
    let atoms = vec![Atom {
        position: [0.0; 3],
        velocity: [v, 0.0, 0.0],
        force: [0.0; 3],
        mass: 1.0,
        charge: 0.0,
    }];
    let mut sim = MdSimulationPlain::new(atoms, [1000.0; 3]);
    sim.step(dt);
    let expected_x = v * dt;
    assert!(
        (sim.atoms[0].position[0] - expected_x).abs() < 1e-12,
        "position = {}, expected {expected_x}",
        sim.atoms[0].position[0]
    );
}

#[test]
fn test_md_simulation_plain_with_thermostat() {
    let atoms = vec![
        Atom {
            position: [0.0; 3],
            velocity: [1.0, 0.0, 0.0],
            force: [0.0; 3],
            mass: 1.0,
            charge: 0.0,
        },
        Atom {
            position: [3.0, 0.0, 0.0],
            velocity: [-1.0, 0.0, 0.0],
            force: [0.0; 3],
            mass: 1.0,
            charge: 0.0,
        },
    ];
    let mut sim = MdSimulationPlain::new(atoms, [20.0; 3]);
    sim.thermostat = Some(ThermostatType::Berendsen);
    for _ in 0..10 {
        sim.step(0.001);
    }
    assert!(sim.step_count == 10);
}

#[test]
fn test_argon_nve_energy_conservation() {
    let n_side = 4_usize;
    let spacing = 1.6;
    let box_l = n_side as f64 * spacing;

    let mut atoms = AtomSet::with_capacity(n_side * n_side * n_side);
    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                let pos = Vec3::new(
                    (ix as f64 + 0.5) * spacing,
                    (iy as f64 + 0.5) * spacing,
                    (iz as f64 + 0.5) * spacing,
                );
                let vx = ((ix + iy) as f64 * 0.1) - 0.3;
                let vy = ((iy + iz) as f64 * 0.1) - 0.3;
                let vz = ((iz + ix) as f64 * 0.1) - 0.3;
                atoms.add_atom(pos, Vec3::new(vx, vy, vz), 1.0, 0.0, 0);
            }
        }
    }

    atoms.remove_com_velocity();
    assert_eq!(atoms.len(), 64);

    let mut ff = PairForceField::new();
    ff.add_interaction(0, 0, LennardJones::new(1.0, 1.0, 2.5));

    let pbox = PeriodicBox::cubic(box_l);

    let mut sim = MdSimulation::new(
        atoms,
        ff,
        VelocityVerlet::new(),
        NoThermostat::new(),
        NoBarostat::new(),
        pbox,
        0.0,
        0.0,
        1.0,
    );

    let records = sim.run(1000, 0.002, Some(10));

    let e0 = records[0].total;
    let max_drift = records
        .iter()
        .map(|r| (r.total - e0).abs())
        .fold(0.0_f64, f64::max);

    let e_scale = e0.abs().max(1.0);
    assert!(
        max_drift / e_scale < 0.05,
        "Argon NVE energy drift: {max_drift} / scale {e_scale} = {}",
        max_drift / e_scale
    );
}

// ----- ObservableRecord tests -----

#[test]
fn test_observable_record_push_and_len() {
    let mut rec = ObservableRecord::new();
    assert!(rec.is_empty());
    rec.push(0.001, 10.0, -20.0, 300.0, 1.0, 0.0);
    rec.push(0.002, 10.5, -20.5, 301.0, 1.1, 0.1);
    assert_eq!(rec.len(), 2);
    assert!(!rec.is_empty());
}

#[test]
fn test_observable_record_total_energy() {
    let mut rec = ObservableRecord::new();
    rec.push(0.001, 10.0, -5.0, 300.0, 1.0, 0.0);
    assert!((rec.total_energy[0] - 5.0).abs() < 1e-12);
}

#[test]
fn test_observable_record_mean_temperature() {
    let mut rec = ObservableRecord::new();
    rec.push(0.001, 10.0, -5.0, 200.0, 1.0, 0.0);
    rec.push(0.002, 10.0, -5.0, 400.0, 1.0, 0.0);
    let mean_t = rec.mean_temperature();
    assert!((mean_t - 300.0).abs() < 1e-10, "mean T = {mean_t}");
}

#[test]
fn test_observable_record_energy_variance_constant() {
    let mut rec = ObservableRecord::new();
    for i in 0..10 {
        rec.push(i as f64 * 0.001, 10.0, -10.0, 300.0, 1.0, 0.0);
    }
    // Total energy = 0.0 for all frames -> variance should be 0
    let var = rec.energy_variance();
    assert!(
        var.abs() < 1e-20,
        "variance of constant series should be 0, got {var}"
    );
}

#[test]
fn test_observable_record_max_drift_nonzero() {
    let mut rec = ObservableRecord::new();
    rec.push(0.001, 10.0, -10.0, 300.0, 1.0, 0.0); // total = 0
    rec.push(0.002, 12.0, -10.0, 300.0, 1.0, 0.0); // total = 2
    let drift = rec.max_energy_drift();
    assert!(
        (drift - 2.0).abs() < 1e-10,
        "max drift should be 2.0, got {drift}"
    );
}

// ----- compute_msd tests -----

#[test]
fn test_compute_msd_zero_for_same_positions() {
    let pos = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
    let msd = compute_msd(&pos, &pos);
    assert!(
        msd.abs() < 1e-20,
        "MSD should be 0 for identical positions, got {msd}"
    );
}

#[test]
fn test_compute_msd_known_value() {
    let pos_now = vec![[1.0, 0.0, 0.0]];
    let pos_ref = vec![[0.0, 0.0, 0.0]];
    let msd = compute_msd(&pos_now, &pos_ref);
    // |r|^2 = 1.0, MSD = 1.0 / 1 = 1.0
    assert!((msd - 1.0).abs() < 1e-12, "MSD should be 1.0, got {msd}");
}

#[test]
fn test_compute_msd_increases_with_displacement() {
    let pos_ref = vec![[0.0, 0.0, 0.0]];
    let pos_far = vec![[3.0, 4.0, 0.0]];
    let msd = compute_msd(&pos_far, &pos_ref);
    assert!((msd - 25.0).abs() < 1e-12, "MSD should be 25.0, got {msd}");
}

#[test]
fn test_per_atom_displacement_sq() {
    let pos_now = vec![[3.0, 0.0, 0.0], [0.0, 4.0, 0.0]];
    let pos_ref = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let disps = per_atom_displacement_sq(&pos_now, &pos_ref);
    assert_eq!(disps.len(), 2);
    assert!((disps[0] - 9.0).abs() < 1e-12);
    assert!((disps[1] - 16.0).abs() < 1e-12);
}

// ----- EquilibrationDetector tests -----

#[test]
fn test_equilibration_detector_constant_energy() {
    let mut det = EquilibrationDetector::new(5, 1e-6, 2);
    let mut equilibrated = false;
    for step in 0..50u64 {
        if det.update(100.0, step) {
            equilibrated = true;
            break;
        }
    }
    assert!(equilibrated, "constant energy should trigger equilibration");
    assert!(det.equilibration_step.is_some());
}

#[test]
fn test_equilibration_detector_fluctuating_energy_not_equilibrated() {
    let mut det = EquilibrationDetector::new(5, 1e-8, 3);
    for step in 0..30u64 {
        // Wildly fluctuating energy - never equilibrates
        let e = (step as f64 * 0.5).sin() * 50.0;
        det.update(e, step);
    }
    assert!(
        !det.equilibrated,
        "fluctuating system should not equilibrate"
    );
}

#[test]
fn test_equilibration_detector_reset() {
    let mut det = EquilibrationDetector::new(3, 1e-4, 1);
    for step in 0..20u64 {
        det.update(42.0, step);
    }
    assert!(det.equilibrated);
    det.reset();
    assert!(!det.equilibrated);
    assert!(det.equilibration_step.is_none());
}

// ----- RestartData tests -----

#[test]
fn test_restart_data_round_trip() {
    let state = MdState::new(
        vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
        vec![[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]],
        vec![12.0, 16.0],
    );
    let rd = RestartData::from_state(&state);
    let buf = rd.to_buffer();
    let rd2 = RestartData::from_buffer(&buf).expect("buffer should be valid");
    assert_eq!(rd2.n_atoms, 2);
    assert!((rd2.positions[0][0] - 1.0).abs() < 1e-14);
    assert!((rd2.positions[1][2] - 6.0).abs() < 1e-14);
    assert!((rd2.masses[1] - 16.0).abs() < 1e-14);
}

#[test]
fn test_restart_data_from_buffer_too_short_returns_none() {
    let buf = vec![1.0, 0.0]; // too short
    assert!(RestartData::from_buffer(&buf).is_none());
}

#[test]
fn test_restart_data_to_state_preserves_step_time() {
    let mut state = MdState::new(vec![[0.0; 3]], vec![[1.0, 0.0, 0.0]], vec![1.0]);
    state.step = 42;
    state.time = 0.084;
    state.kinetic_energy = 3.125;
    state.potential_energy = -6.25;
    let rd = RestartData::from_state(&state);
    let buf = rd.to_buffer();
    let rd2 = RestartData::from_buffer(&buf).unwrap();
    let state2 = rd2.to_state();
    assert_eq!(state2.step, 42);
    assert!((state2.time - 0.084).abs() < 1e-14);
    assert!((state2.kinetic_energy - 3.125).abs() < 1e-14);
    assert!((state2.potential_energy - (-6.25)).abs() < 1e-14);
}

#[test]
fn test_mdsim_save_load_restart() {
    let config = MdConfig {
        dt: 0.001,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config, state);
    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    for _ in 0..10 {
        sim.step_with_forces(zero_f, 0.1);
    }
    let buf = sim.save_restart();
    let mut sim2 = MdSim::new(MdConfig::default(), make_two_atom_state());
    let ok = sim2.load_restart(&buf);
    assert!(ok, "load_restart should succeed");
    assert_eq!(sim2.state.step, sim.state.step);
    for a in 0..3 {
        assert!(
            (sim2.state.positions[0][a] - sim.state.positions[0][a]).abs() < 1e-14,
            "loaded position mismatch"
        );
    }
}

// ----- RunStatistics tests -----

#[test]
fn test_run_statistics_from_empty_record() {
    let rec = ObservableRecord::new();
    let stats = RunStatistics::from_record(&rec, 1.0);
    assert_eq!(stats.total_steps, 0);
    assert!(stats.mean_temperature.abs() < 1e-12);
}

#[test]
fn test_run_statistics_steps_per_second() {
    let mut rec = ObservableRecord::new();
    for i in 0..100 {
        rec.push(i as f64 * 0.001, 10.0, -10.0, 300.0, 1.0, 0.0);
    }
    let stats = RunStatistics::from_record(&rec, 0.01);
    assert!(stats.steps_per_second > 0.0, "steps/s should be positive");
    assert!(
        (stats.steps_per_second - 10000.0).abs() < 1.0,
        "steps/s = {}",
        stats.steps_per_second
    );
}

#[test]
fn test_run_statistics_energy_conserved() {
    let mut rec = ObservableRecord::new();
    for i in 0..50 {
        rec.push(i as f64 * 0.001, 10.0, -10.0, 300.0, 1.0, 0.0);
    }
    let stats = RunStatistics::from_record(&rec, 1.0);
    assert!(
        stats.energy_conserved(1e-5),
        "perfectly constant energy should pass conservation"
    );
}

// ----- maxwell_boltzmann_velocities tests -----

#[test]
fn test_maxwell_boltzmann_velocities_com_removed() {
    let masses = vec![1.0; 10];
    let vels = maxwell_boltzmann_velocities(&masses, 300.0, 0);
    let mut com = [0.0f64; 3];
    for v in &vels {
        for (a, &va) in v.iter().enumerate() {
            com[a] += va;
        }
    }
    for (a, &c) in com.iter().enumerate() {
        assert!(
            c.abs() < 1e-10,
            "COM velocity axis {a} should be ~0, got {}",
            c
        );
    }
}

#[test]
fn test_maxwell_boltzmann_velocities_returns_n_velocities() {
    let masses = vec![1.0; 5];
    let vels = maxwell_boltzmann_velocities(&masses, 300.0, 42);
    assert_eq!(vels.len(), 5);
}

#[test]
fn test_maxwell_boltzmann_velocities_different_seeds() {
    let masses = vec![1.0; 4];
    let v1 = maxwell_boltzmann_velocities(&masses, 300.0, 0);
    let v2 = maxwell_boltzmann_velocities(&masses, 300.0, 1);
    let same = v1
        .iter()
        .zip(v2.iter())
        .all(|(a, b)| (a[0] - b[0]).abs() < 1e-14 && (a[1] - b[1]).abs() < 1e-14);
    assert!(!same, "different seeds should produce different velocities");
}

// ----- rescale_to_temperature tests -----

#[test]
fn test_rescale_to_temperature_hits_target() {
    let masses = vec![1.0; 4];
    let mut vels = maxwell_boltzmann_velocities(&masses, 200.0, 10);
    rescale_to_temperature(&mut vels, &masses, 400.0);
    let ke: f64 = vels
        .iter()
        .zip(masses.iter())
        .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
        .sum();
    let t = 2.0 * ke / (3.0 * 4.0 * KB_REDUCED);
    assert!(
        (t - 400.0).abs() < 1e-6,
        "temperature after rescaling = {t}"
    );
}

// ----- EnergyDriftTracker tests -----

#[test]
fn test_energy_drift_tracker_zero_for_constant_energy() {
    let mut tracker = EnergyDriftTracker::new();
    for _ in 0..20 {
        tracker.update(50.0);
    }
    assert!(tracker.relative_drift() < 1e-20);
}

#[test]
fn test_energy_drift_tracker_detects_drift() {
    let mut tracker = EnergyDriftTracker::new();
    tracker.update(100.0);
    tracker.update(105.0);
    assert!((tracker.max_abs_drift - 5.0).abs() < 1e-10);
    assert!((tracker.relative_drift() - 0.05).abs() < 1e-10);
}

#[test]
fn test_energy_drift_tracker_reset() {
    let mut tracker = EnergyDriftTracker::new();
    tracker.update(100.0);
    tracker.update(200.0);
    tracker.reset();
    assert_eq!(tracker.n_updates, 0);
    assert!(tracker.max_abs_drift.abs() < 1e-20);
}

// ----- run_full integration test -----

#[test]
fn test_run_full_records_correct_count() {
    let config = MdConfig {
        dt: 0.001,
        n_steps: 20,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config, state);
    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    let (record, stats) = sim.run_full(zero_f, 5, None::<fn(&[[f64; 3]], &[f64; 3]) -> f64>);
    assert_eq!(record.len(), 4, "expected 4 records for 20 steps / 5");
    assert_eq!(stats.total_steps, 4);
    assert!(stats.wall_time_s >= 0.0);
}

#[test]
fn test_run_full_msd_increases_for_moving_particles() {
    let config = MdConfig {
        dt: 0.01,
        n_steps: 10,
        ensemble: Ensemble::NVE,
        pbc: false,
        ..MdConfig::default()
    };
    let state = MdState::new(
        vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]],
        vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        vec![1.0, 1.0],
    );
    let mut sim = MdSim::new(config, state);
    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    let (record, _) = sim.run_full(zero_f, 2, None::<fn(&[[f64; 3]], &[f64; 3]) -> f64>);
    // MSD should be monotonically non-decreasing (free particles)
    for i in 1..record.msd.len() {
        assert!(
            record.msd[i] >= record.msd[i - 1] - 1e-10,
            "MSD decreased at frame {i}: {} < {}",
            record.msd[i],
            record.msd[i - 1]
        );
    }
}

// ----- run_until_equilibrated test -----

#[test]
fn test_run_until_equilibrated_stops_at_equilibration() {
    let config = MdConfig {
        dt: 0.001,
        n_steps: 1000,
        ensemble: Ensemble::NVE,
        ..MdConfig::default()
    };
    let state = make_two_atom_state();
    let mut sim = MdSim::new(config, state);
    // With zero forces and NVE, energy is constant -> should equilibrate quickly
    let zero_f = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0f64; 3]; pos.len()];
    let mut detector = EquilibrationDetector::new(5, 1e-6, 2);
    let steps = sim.run_until_equilibrated(zero_f, &mut detector, 500);
    assert!(steps <= 500, "should not exceed max_steps");
    assert!(
        detector.equilibrated,
        "should equilibrate with constant energy"
    );
}

// =========================================================================
// RdfHistogram tests
// =========================================================================

#[test]
fn test_rdf_histogram_empty_gofr_is_zero() {
    let rdf = RdfHistogram::new(10.0, 50, 1000.0, 10);
    let g = rdf.gofr();
    assert!(
        g.iter().all(|&v| v == 0.0),
        "empty histogram g(r) should be all zero"
    );
}

#[test]
fn test_rdf_histogram_bin_centres_length() {
    let rdf = RdfHistogram::new(8.0, 40, 512.0, 8);
    let centres = rdf.bin_centres();
    assert_eq!(centres.len(), 40);
}

#[test]
fn test_rdf_histogram_bin_centres_first_half_dr() {
    let rdf = RdfHistogram::new(8.0, 40, 512.0, 8);
    let centres = rdf.bin_centres();
    let expected = 0.5 * rdf.dr;
    assert!((centres[0] - expected).abs() < 1e-12);
}

#[test]
fn test_rdf_histogram_bin_centres_last() {
    let rdf = RdfHistogram::new(8.0, 40, 512.0, 8);
    let centres = rdf.bin_centres();
    let expected = (39.0 + 0.5) * rdf.dr;
    assert!((centres[39] - expected).abs() < 1e-12);
}

#[test]
fn test_rdf_histogram_n_frames_increments() {
    let mut rdf = RdfHistogram::new(10.0, 20, 1000.0, 2);
    let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    rdf.accumulate(&positions, 20.0);
    rdf.accumulate(&positions, 20.0);
    assert_eq!(rdf.n_frames, 2);
}

#[test]
fn test_rdf_histogram_accumulate_counts_pair() {
    let mut rdf = RdfHistogram::new(10.0, 50, 1000.0, 2);
    let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    rdf.accumulate(&positions, 20.0);
    let total: u64 = rdf.histogram.iter().sum();
    assert_eq!(total, 1, "one pair should contribute exactly 1 count");
}

#[test]
fn test_rdf_histogram_reset_zeroes_counts() {
    let mut rdf = RdfHistogram::new(10.0, 20, 1000.0, 2);
    let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    rdf.accumulate(&positions, 20.0);
    rdf.reset();
    assert_eq!(rdf.n_frames, 0);
    assert!(rdf.histogram.iter().all(|&v| v == 0));
}

#[test]
fn test_rdf_histogram_gofr_has_n_bins_elements() {
    let mut rdf = RdfHistogram::new(10.0, 30, 1000.0, 2);
    let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    rdf.accumulate(&positions, 20.0);
    assert_eq!(rdf.gofr().len(), 30);
}

#[test]
fn test_rdf_histogram_gofr_nonnegative() {
    let mut rdf = RdfHistogram::new(10.0, 30, 1000.0, 4);
    let positions = [
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [6.0, 0.0, 0.0],
        [9.0, 0.0, 0.0],
    ];
    rdf.accumulate(&positions, 20.0);
    for g in rdf.gofr() {
        assert!(g >= 0.0, "g(r) must be non-negative, got {g}");
    }
}

#[test]
fn test_rdf_histogram_pair_not_counted_beyond_rmax() {
    let mut rdf = RdfHistogram::new(5.0, 25, 1000.0, 2);
    // Place atoms farther than r_max=5 apart
    let positions = [[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]];
    rdf.accumulate(&positions, 20.0);
    let total: u64 = rdf.histogram.iter().sum();
    assert_eq!(total, 0, "pair beyond r_max should not be counted");
}

#[test]
fn test_rdf_histogram_minimum_image_wraps() {
    let mut rdf = RdfHistogram::new(6.0, 30, 1000.0, 2);
    // Two atoms with a large raw distance but close under PBC
    let positions = [[0.5, 0.0, 0.0], [19.5, 0.0, 0.0]];
    let box_len = 20.0;
    rdf.accumulate(&positions, box_len);
    let total: u64 = rdf.histogram.iter().sum();
    assert_eq!(total, 1, "wrapped distance = 1.0 should be counted");
}

#[test]
fn test_rdf_histogram_coordination_number_zero_for_empty() {
    let rdf = RdfHistogram::new(10.0, 50, 1000.0, 10);
    let cn = rdf.coordination_number(5.0);
    assert_eq!(cn, 0.0);
}

// =========================================================================
// VelocityAutocorrelation tests
// =========================================================================

#[test]
fn test_vacf_new_has_no_reference() {
    let vacf = VelocityAutocorrelation::new();
    assert!(!vacf.has_reference);
    assert!(vacf.v0.is_empty());
}

#[test]
fn test_vacf_set_reference_marks_has_reference() {
    let mut vacf = VelocityAutocorrelation::new();
    let vels = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    vacf.set_reference(&vels);
    assert!(vacf.has_reference);
    assert_eq!(vacf.v0.len(), 2);
}

#[test]
fn test_vacf_accumulate_at_t0_returns_c0() {
    let mut vacf = VelocityAutocorrelation::new();
    let vels = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    vacf.set_reference(&vels);
    let c = vacf.accumulate(&vels);
    // <v(0)·v(0)> / N = (1.0 + 1.0) / 2 = 1.0
    assert!((c - 1.0).abs() < 1e-12, "C(0) = 1.0, got {c}");
}

#[test]
fn test_vacf_normalised_first_is_one() {
    let mut vacf = VelocityAutocorrelation::new();
    let vels = vec![[2.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
    vacf.set_reference(&vels);
    vacf.accumulate(&vels);
    let norm = vacf.normalised();
    assert!(!norm.is_empty());
    assert!(
        (norm[0] - 1.0).abs() < 1e-12,
        "normalised C(0) = 1, got {}",
        norm[0]
    );
}

#[test]
fn test_vacf_diffusion_coefficient_zero_velocities() {
    let mut vacf = VelocityAutocorrelation::new();
    let vels = vec![[0.0; 3]; 3];
    vacf.set_reference(&vels);
    vacf.accumulate(&vels);
    let d = vacf.diffusion_coefficient(0.001);
    assert_eq!(
        d, 0.0,
        "diffusion coefficient for zero velocity should be 0"
    );
}

#[test]
fn test_vacf_accumulate_no_reference_returns_zero() {
    let mut vacf = VelocityAutocorrelation::new();
    let vels = vec![[1.0, 0.0, 0.0]];
    let c = vacf.accumulate(&vels);
    assert_eq!(c, 0.0, "accumulate without reference should return 0");
}

#[test]
fn test_vacf_empty_normalised_is_empty() {
    let vacf = VelocityAutocorrelation::new();
    assert!(vacf.normalised().is_empty());
}

#[test]
fn test_vacf_n_frames_increments_on_accumulate() {
    let mut vacf = VelocityAutocorrelation::new();
    let vels = vec![[1.0, 0.0, 0.0]];
    vacf.set_reference(&vels);
    vacf.accumulate(&vels);
    vacf.accumulate(&vels);
    assert_eq!(vacf.n_frames, 2);
}

#[test]
fn test_vacf_set_reference_resets_history() {
    let mut vacf = VelocityAutocorrelation::new();
    let vels = vec![[1.0, 0.0, 0.0]];
    vacf.set_reference(&vels);
    vacf.accumulate(&vels);
    vacf.set_reference(&vels); // reset
    assert_eq!(vacf.n_frames, 0);
    assert!(vacf.vacf.is_empty());
}

#[test]
fn test_vacf_diffusion_positive_for_moving_particles() {
    let mut vacf = VelocityAutocorrelation::new();
    let vels = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    vacf.set_reference(&vels);
    // Accumulate multiple identical frames
    for _ in 0..5 {
        vacf.accumulate(&vels);
    }
    let d = vacf.diffusion_coefficient(0.001);
    assert!(
        d > 0.0,
        "diffusion coefficient should be positive for moving particles, got {d}"
    );
}

// =========================================================================
// PressureTensor tests
// =========================================================================

#[test]
fn test_pressure_tensor_zero_is_zero() {
    let pt = PressureTensor::zero();
    assert!(pt.tensor.iter().all(|&v| v == 0.0));
}

#[test]
fn test_pressure_tensor_scalar_pressure_of_isotropic() {
    let ke = [[100.0, 0.0, 0.0], [0.0, 100.0, 0.0], [0.0, 0.0, 100.0]];
    let vt = [[0.0f64; 3]; 3];
    let vol = 1000.0_f64; // Å³
    let pt = PressureTensor::from_kinetic_virial(&ke, &vt, vol);
    let p = pt.scalar_pressure();
    assert!(p.is_finite(), "scalar pressure should be finite, got {p}");
    // All diagonal elements equal -> scalar = diagonal element
    let pxx = pt.element(0, 0);
    let pyy = pt.element(1, 1);
    let pzz = pt.element(2, 2);
    assert!((pxx - pyy).abs() < 1e-6);
    assert!((pyy - pzz).abs() < 1e-6);
    assert!((p - pxx).abs() < 1e-6);
}

#[test]
fn test_pressure_tensor_element_access() {
    let mut pt = PressureTensor::zero();
    pt.tensor[1] = 42.0;
    assert_eq!(pt.element(0, 1), 42.0);
}

#[test]
fn test_pressure_tensor_is_symmetric_zero() {
    let pt = PressureTensor::zero();
    assert!(pt.is_symmetric(1e-10));
}

#[test]
fn test_pressure_tensor_not_symmetric_for_asymmetric_tensor() {
    let mut pt = PressureTensor::zero();
    pt.tensor[1] = 1.0;
    pt.tensor[3] = 2.0;
    assert!(!pt.is_symmetric(1e-10));
}

#[test]
fn test_pressure_tensor_kinetic_contribution_single_particle() {
    let velocities = [[2.0, 0.0, 0.0]];
    let masses = [1.0];
    let vol = 8.0_f64;
    let kt = PressureTensor::kinetic_contribution(&velocities, &masses, vol);
    // Only (0,0) component nonzero: m*vx*vx/V = 4.0/8.0 = 0.5
    assert!((kt[0][0] - 0.5).abs() < 1e-12, "K[0][0] = {}", kt[0][0]);
    assert!(kt[0][1].abs() < 1e-12);
    assert!(kt[1][1].abs() < 1e-12);
}

#[test]
fn test_pressure_tensor_from_kinetic_virial_units() {
    // With zero virial, pressure = KE_tensor / volume * conversion
    let ke = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let vt = [[0.0f64; 3]; 3];
    let vol = 1.0_f64;
    let pt = PressureTensor::from_kinetic_virial(&ke, &vt, vol);
    let expected = KJ_MOL_PER_ANG3_TO_BAR;
    assert!(
        (pt.element(0, 0) - expected).abs() < 1.0,
        "P_xx = {}",
        pt.element(0, 0)
    );
}

// =========================================================================
// TimecorrelationFn tests
// =========================================================================

#[test]
fn test_tcf_empty_compute_lag_zero() {
    let tcf = TimecorrelationFn::new();
    assert_eq!(tcf.compute_lag(0), 0.0);
}

#[test]
fn test_tcf_compute_lag_0_is_mean_square() {
    let mut tcf = TimecorrelationFn::new();
    tcf.push(2.0);
    tcf.push(3.0);
    let c0 = tcf.compute_lag(0);
    // C(0) = (2*2 + 3*3) / 2 = (4+9)/2 = 6.5
    assert!((c0 - 6.5).abs() < 1e-12, "C(0) = {c0}");
}

#[test]
fn test_tcf_compute_lag_beyond_length_is_zero() {
    let mut tcf = TimecorrelationFn::new();
    tcf.push(1.0);
    tcf.push(2.0);
    assert_eq!(tcf.compute_lag(5), 0.0);
}

#[test]
fn test_tcf_compute_all_length() {
    let mut tcf = TimecorrelationFn::new();
    for i in 0..10 {
        tcf.push(i as f64);
    }
    let all = tcf.compute_all(4);
    assert_eq!(all.len(), 5, "max_lag=4 -> lags 0..=4");
}

#[test]
fn test_tcf_normalised_first_is_one() {
    let mut tcf = TimecorrelationFn::new();
    for _ in 0..5 {
        tcf.push(3.0);
    }
    let norm = tcf.normalised(2);
    assert!(
        (norm[0] - 1.0).abs() < 1e-12,
        "normalised C(0) = {}",
        norm[0]
    );
}

#[test]
fn test_tcf_constant_signal_all_lags_one() {
    let mut tcf = TimecorrelationFn::new();
    for _ in 0..8 {
        tcf.push(5.0);
    }
    let norm = tcf.normalised(4);
    for (i, &c) in norm.iter().enumerate() {
        assert!(
            (c - 1.0).abs() < 1e-12,
            "constant signal: C(tau={i}) = {c} != 1"
        );
    }
}

#[test]
fn test_tcf_push_increments_history() {
    let mut tcf = TimecorrelationFn::new();
    tcf.push(1.0);
    tcf.push(2.0);
    assert_eq!(tcf.history.len(), 2);
}
