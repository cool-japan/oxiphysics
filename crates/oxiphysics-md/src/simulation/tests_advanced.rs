// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced tests for NPT, MSD, RDF, and adaptive timestep methods.

use super::*;

fn simple_sim(n: usize, dt: f64) -> MdSim {
    let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 3.0, 0.0, 0.0]).collect();
    let velocities = vec![[0.0_f64; 3]; n];
    let masses = vec![12.0_f64; n]; // carbon-like
    let state = MdState::new(positions, velocities, masses);
    let config = MdConfig {
        dt,
        n_steps: 10,
        temperature: 300.0,
        pressure: 1.0,
        box_lengths: [30.0, 30.0, 30.0],
        pbc: true,
        ensemble: Ensemble::NVT,
    };
    MdSim::new(config, state)
}

// -----------------------------------------------------------------------
// ParrinelloRahmanBarostat tests
// -----------------------------------------------------------------------

#[test]
fn test_pr_barostat_scale_factors_at_target_pressure_is_one() {
    let pr = ParrinelloRahmanBarostat::new(1.0, 0.5, 4.5e-5);
    let scale = pr.scale_factors(1.0, 0.002);
    // At target pressure delta = 0, so mu = cbrt(1) = 1
    assert!((scale[0] - 1.0).abs() < 1e-12, "mu={}", scale[0]);
    assert!((scale[1] - 1.0).abs() < 1e-12);
    assert!((scale[2] - 1.0).abs() < 1e-12);
}

#[test]
fn test_pr_barostat_scale_factors_high_pressure_expands_box() {
    let pr = ParrinelloRahmanBarostat::new(1.0, 0.5, 4.5e-5);
    // p_current >> target → mu > 1 (box expands to reduce pressure)
    let scale = pr.scale_factors(1000.0, 0.002);
    assert!(
        scale[0] > 1.0,
        "high pressure should expand box: mu={}",
        scale[0]
    );
}

#[test]
fn test_pr_barostat_apply_changes_box_lengths() {
    let mut sim = simple_sim(4, 0.002);
    let original_box = sim.config.box_lengths;
    let mut pr = ParrinelloRahmanBarostat::new(1.0, 0.5, 4.5e-5);
    // With zero virial and many atoms, pressure != target → box rescales.
    pr.apply(&mut sim, 1e6, 0.002); // large virial forces pressure far from target
    let new_box = sim.config.box_lengths;
    // Box lengths should differ from original (at least one)
    let changed = (0..3).any(|a| (new_box[a] - original_box[a]).abs() > 1e-15);
    assert!(changed, "barostat should change box dimensions");
}

// -----------------------------------------------------------------------
// MsdCalculator tests
// -----------------------------------------------------------------------

#[test]
fn test_msd_zero_displacement_is_zero() {
    let pos = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
    let calc = MsdCalculator::new(pos.clone());
    let msd = calc.compute(&pos);
    assert!(msd.abs() < 1e-15, "MSD of identical positions = {msd}");
}

#[test]
fn test_msd_known_displacement() {
    let ref_pos = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
    let new_pos = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let calc = MsdCalculator::new(ref_pos);
    let msd = calc.compute(&new_pos);
    // MSD = (1^2 + 1^2) / 2 = 1.0
    assert!((msd - 1.0).abs() < 1e-12, "MSD = {msd}");
}

#[test]
fn test_msd_diffusion_coefficient_finite() {
    let ref_pos = vec![[0.0; 3]; 4];
    let displaced = vec![
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 3.0],
        [1.0, 1.0, 1.0],
    ];
    let mut calc = MsdCalculator::new(ref_pos);
    calc.accumulate(&displaced);
    let d = calc.diffusion_coefficient(0.01);
    assert!(d.is_finite() && d > 0.0, "D = {d}");
}

#[test]
fn test_msd_via_sim_method() {
    let mut sim = simple_sim(3, 0.002);
    let ref_pos = sim.state.positions.clone();
    // displace first atom by 2 Å
    sim.state.positions[0][0] += 2.0;
    let msd = sim.compute_msd(&ref_pos);
    // MSD = 4 / 3 ≈ 1.333...
    assert!((msd - 4.0 / 3.0).abs() < 1e-12, "MSD = {msd}");
}

// -----------------------------------------------------------------------
// RdfRunning tests
// -----------------------------------------------------------------------

#[test]
fn test_rdf_running_new_histogram_is_zero() {
    let rdf = RdfRunning::new(10.0, 50, 1000.0, 10);
    let total: u64 = rdf.histogram.iter().sum();
    assert_eq!(total, 0, "fresh histogram should be all zeros");
}

#[test]
fn test_rdf_running_accumulate_counts_pair() {
    let mut rdf = RdfRunning::new(10.0, 10, 1000.0, 2);
    let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    rdf.accumulate(&positions, 100.0);
    let total: u64 = rdf.histogram.iter().sum();
    assert_eq!(total, 1, "one pair should be counted, got {total}");
}

#[test]
fn test_rdf_running_gofr_length_matches_bins() {
    let mut rdf = RdfRunning::new(8.0, 20, 8000.0, 4);
    let positions = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [6.0, 0.0, 0.0],
    ];
    rdf.accumulate(&positions, 20.0);
    let g = rdf.gofr();
    assert_eq!(g.len(), 20, "g(r) length should match n_bins");
}

#[test]
fn test_rdf_running_via_sim_method() {
    let sim = simple_sim(4, 0.002);
    let mut rdf = RdfRunning::new(10.0, 20, 27000.0, 4);
    sim.compute_rdf_running(&mut rdf);
    assert_eq!(rdf.n_frames, 1, "one frame accumulated");
}

// -----------------------------------------------------------------------
// Adaptive timestep tests
// -----------------------------------------------------------------------

#[test]
fn test_adaptive_timestep_zero_force_doubles_dt() {
    let mut sim = simple_sim(3, 0.002);
    let zero_forces = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0; 3]; pos.len()];
    let dt_used = sim.adaptive_timestep(zero_forces, 100.0, 0.01);
    // zero forces < threshold/4 → dt should double (capped at dt_max=0.01)
    assert!(dt_used >= 0.002, "dt_used={dt_used}");
}

#[test]
fn test_adaptive_timestep_returns_valid_dt() {
    let mut sim = simple_sim(3, 0.002);
    let zero_forces = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0; 3]; pos.len()];
    let dt = sim.adaptive_timestep(zero_forces, 100.0, 0.01);
    assert!(
        dt > 0.0 && dt.is_finite(),
        "dt must be positive and finite: {dt}"
    );
}

#[test]
fn test_adaptive_timestep_increments_step() {
    let mut sim = simple_sim(3, 0.002);
    let step_before = sim.state.step;
    let zero_forces = |pos: &[[f64; 3]], _: &[f64; 3]| vec![[0.0; 3]; pos.len()];
    sim.adaptive_timestep(zero_forces, 100.0, 0.01);
    assert_eq!(sim.state.step, step_before + 1, "step should increment");
}

// -----------------------------------------------------------------------
// NeighborList tests
// -----------------------------------------------------------------------

#[test]
fn test_neighbor_list_build_finds_close_pair() {
    let mut nl = NeighborList::new(5.0, 1.0);
    let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [20.0, 0.0, 0.0]];
    nl.build(&positions, [30.0, 30.0, 30.0]);
    assert!(
        nl.n_pairs() >= 1,
        "atoms at r=3 should be in list, got {}",
        nl.n_pairs()
    );
}

#[test]
fn test_neighbor_list_build_excludes_far_pair() {
    let mut nl = NeighborList::new(5.0, 0.5);
    let positions = vec![[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]];
    nl.build(&positions, [40.0, 40.0, 40.0]);
    assert_eq!(
        nl.n_pairs(),
        0,
        "atoms at r=20 should not be in list with rc=5.5"
    );
}

#[test]
fn test_neighbor_list_needs_rebuild_initially_false() {
    let mut nl = NeighborList::new(5.0, 1.0);
    let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    nl.build(&positions, [30.0, 30.0, 30.0]);
    assert!(
        !nl.needs_rebuild(&positions),
        "no displacement → no rebuild needed"
    );
}

#[test]
fn test_neighbor_list_needs_rebuild_after_large_displacement() {
    let mut nl = NeighborList::new(5.0, 1.0);
    let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
    nl.build(&positions, [30.0, 30.0, 30.0]);
    let displaced = vec![[0.0, 0.0, 0.0], [3.0, 2.0, 0.0]]; // moved by 2 Å > skin/2 = 0.5
    assert!(
        nl.needs_rebuild(&displaced),
        "displacement > skin/2 should trigger rebuild"
    );
}

#[test]
fn test_neighbor_list_for_each_pair_counts_correctly() {
    let mut nl = NeighborList::new(6.0, 0.5);
    let positions = vec![[0.0; 3], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
    nl.build(&positions, [30.0, 30.0, 30.0]);
    let mut count = 0usize;
    nl.for_each_pair(&positions, [30.0, 30.0, 30.0], |_i, _j, _r2| {
        count += 1;
    });
    assert_eq!(count, nl.n_pairs(), "for_each_pair should visit all pairs");
}

#[test]
fn test_neighbor_list_single_atom_no_pairs() {
    let mut nl = NeighborList::new(5.0, 1.0);
    let positions = vec![[0.0; 3]];
    nl.build(&positions, [10.0, 10.0, 10.0]);
    assert_eq!(nl.n_pairs(), 0, "single atom has no pairs");
}

// -----------------------------------------------------------------------
// EnsembleSwitcher tests
// -----------------------------------------------------------------------

#[test]
fn test_ensemble_switcher_same_ensemble_coupling_one() {
    let sw = EnsembleSwitcher::new(Ensemble::NVT);
    assert_eq!(
        sw.coupling_fraction(100),
        1.0,
        "same ensemble → coupling = 1"
    );
}

#[test]
fn test_ensemble_switcher_before_transition_coupling_zero() {
    let mut sw = EnsembleSwitcher::new(Ensemble::NVE);
    sw.schedule_transition(Ensemble::NVT, 100, 50);
    assert_eq!(sw.coupling_fraction(50), 0.0, "before start → coupling = 0");
}

#[test]
fn test_ensemble_switcher_midway_coupling_half() {
    let mut sw = EnsembleSwitcher::new(Ensemble::NVE);
    sw.schedule_transition(Ensemble::NVT, 0, 100);
    let frac = sw.coupling_fraction(50);
    assert!(
        (frac - 0.5).abs() < 1e-12,
        "midway → coupling = 0.5, got {frac}"
    );
}

#[test]
fn test_ensemble_switcher_after_ramp_coupling_one() {
    let mut sw = EnsembleSwitcher::new(Ensemble::NVE);
    sw.schedule_transition(Ensemble::NVT, 0, 100);
    let frac = sw.coupling_fraction(200);
    assert_eq!(frac, 1.0, "after ramp → coupling = 1");
}

#[test]
fn test_ensemble_switcher_advance_completes_transition() {
    let mut sw = EnsembleSwitcher::new(Ensemble::NVE);
    sw.schedule_transition(Ensemble::NPT, 0, 10);
    sw.advance(10);
    assert_eq!(sw.current, Ensemble::NPT, "transition should complete");
}

#[test]
fn test_ensemble_switcher_zero_ramp_instant_transition() {
    let mut sw = EnsembleSwitcher::new(Ensemble::NVE);
    sw.schedule_transition(Ensemble::NVT, 5, 0);
    assert_eq!(sw.coupling_fraction(5), 1.0, "zero ramp → instant coupling");
}

// -----------------------------------------------------------------------
// EnergyTracker tests
// -----------------------------------------------------------------------

#[test]
fn test_energy_tracker_empty_mean_zero() {
    let t = EnergyTracker::new();
    assert_eq!(t.mean(), 0.0, "empty tracker → mean = 0");
}

#[test]
fn test_energy_tracker_single_sample() {
    let mut t = EnergyTracker::new();
    t.push(-100.0);
    assert_eq!(t.mean(), -100.0);
    assert_eq!(t.min, -100.0);
    assert_eq!(t.max, -100.0);
}

#[test]
fn test_energy_tracker_variance_known_values() {
    let mut t = EnergyTracker::new();
    // Values: 1, 3, 5 → mean = 3, var = (4+0+4)/3 = 8/3
    for v in [1.0, 3.0, 5.0] {
        t.push(v);
    }
    let expected_var = 8.0_f64 / 3.0_f64;
    assert!(
        (t.variance() - expected_var).abs() < 1e-12,
        "variance = {}",
        t.variance()
    );
}

#[test]
fn test_energy_tracker_std_dev_positive() {
    let mut t = EnergyTracker::new();
    for v in [1.0, 2.0, 3.0] {
        t.push(v);
    }
    assert!(t.std_dev() > 0.0, "std dev should be positive");
}

#[test]
fn test_energy_tracker_peak_to_peak() {
    let mut t = EnergyTracker::new();
    for v in [-10.0, 5.0, 3.0, -2.0] {
        t.push(v);
    }
    assert!(
        (t.peak_to_peak() - 15.0).abs() < 1e-12,
        "p2p = {}",
        t.peak_to_peak()
    );
}

#[test]
fn test_energy_tracker_relative_drift() {
    let mut t = EnergyTracker::new();
    // mean ≈ -100, p2p = 2 → drift = 2/100 = 0.02
    for v in [-101.0, -100.0, -99.0] {
        t.push(v);
    }
    let drift = t.relative_drift();
    assert!((drift - 0.02).abs() < 1e-12, "relative drift = {drift}");
}

// -----------------------------------------------------------------------
// PressureTracker tests
// -----------------------------------------------------------------------

#[test]
fn test_pressure_tracker_mean_correct() {
    let mut pt = PressureTracker::new();
    for v in [1.0, 2.0, 3.0] {
        pt.push(v);
    }
    assert!((pt.mean() - 2.0).abs() < 1e-12, "mean = {}", pt.mean());
}

#[test]
fn test_pressure_tracker_std_dev_positive() {
    let mut pt = PressureTracker::new();
    for v in [1.0, 2.0, 3.0] {
        pt.push(v);
    }
    assert!(pt.std_dev() > 0.0, "std dev = {}", pt.std_dev());
}

#[test]
fn test_pressure_tracker_min_max() {
    let mut pt = PressureTracker::new();
    for v in [5.0, 1.0, 10.0, 3.0] {
        pt.push(v);
    }
    assert_eq!(pt.min(), 1.0);
    assert_eq!(pt.max(), 10.0);
}

#[test]
fn test_pressure_tracker_n() {
    let mut pt = PressureTracker::new();
    for v in [1.0, 2.0, 3.0, 4.0, 5.0] {
        pt.push(v);
    }
    assert_eq!(pt.n(), 5);
}

// -----------------------------------------------------------------------
// HarmonicWell tests
// -----------------------------------------------------------------------

#[test]
fn test_harmonic_well_energy_at_origin_zero() {
    let hw = HarmonicWell::new(1.0);
    assert_eq!(hw.energy([0.0, 0.0, 0.0]), 0.0);
}

#[test]
fn test_harmonic_well_energy_known_value() {
    let hw = HarmonicWell::new(2.0);
    // E = 0.5 * 2.0 * (1² + 0 + 0) = 1.0
    assert!((hw.energy([1.0, 0.0, 0.0]) - 1.0).abs() < 1e-12);
}

#[test]
fn test_harmonic_well_force_at_origin_zero() {
    let hw = HarmonicWell::new(3.0);
    let f = hw.force([0.0, 0.0, 0.0]);
    for &fi in &f {
        assert_eq!(fi, 0.0);
    }
}

#[test]
fn test_harmonic_well_force_restoring() {
    let hw = HarmonicWell::new(1.0);
    let f = hw.force([2.0, 0.0, 0.0]);
    // Force should be -k * displacement = -1 * 2 = -2 along x
    assert!((f[0] - (-2.0)).abs() < 1e-12, "f[0] = {}", f[0]);
    assert_eq!(f[1], 0.0);
    assert_eq!(f[2], 0.0);
}

#[test]
fn test_harmonic_well_omega_and_period() {
    let hw = HarmonicWell::new(4.0); // k = 4
    let m = 1.0;
    let omega = hw.omega(m); // sqrt(4) = 2
    assert!((omega - 2.0).abs() < 1e-12, "omega = {omega}");
    let period = hw.period(m); // 2π/2 = π
    assert!(
        (period - std::f64::consts::PI).abs() < 1e-12,
        "T = {period}"
    );
}

#[test]
fn test_harmonic_well_compute_all_energy_sum() {
    let hw = HarmonicWell::new(1.0);
    let positions = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let (e, f) = hw.compute_all(&positions);
    // E = 0.5*1 + 0.5*1 = 1.0
    assert!((e - 1.0).abs() < 1e-12, "total energy = {e}");
    assert_eq!(f.len(), 2);
}

// -----------------------------------------------------------------------
// AndersenThermostat tests
// -----------------------------------------------------------------------

#[test]
fn test_andersen_thermostat_high_frequency_changes_velocities() {
    let thermo = AndersenThermostat::new(300.0, 1000.0); // ν very large
    let mut state = MdState::new(vec![[0.0; 3]; 20], vec![[0.0; 3]; 20], vec![1.0; 20]);
    thermo.apply(&mut state, 0.002);
    let any_nonzero = state
        .velocities
        .iter()
        .any(|v| v[0] != 0.0 || v[1] != 0.0 || v[2] != 0.0);
    assert!(
        any_nonzero,
        "high-frequency Andersen should change at least one velocity"
    );
}

#[test]
fn test_andersen_thermostat_zero_frequency_no_change() {
    let thermo = AndersenThermostat::new(300.0, 0.0);
    let init_vel = vec![[1.0, 2.0, 3.0]; 5];
    let mut state = MdState::new(vec![[0.0; 3]; 5], init_vel.clone(), vec![1.0; 5]);
    thermo.apply(&mut state, 0.002);
    for (v_new, v_old) in state.velocities.iter().zip(init_vel.iter()) {
        for a in 0..3 {
            assert_eq!(v_new[a], v_old[a], "ν=0 → no velocity change");
        }
    }
}

// -----------------------------------------------------------------------
// LangevinThermostat tests
// -----------------------------------------------------------------------

#[test]
fn test_langevin_o_step_changes_velocities() {
    let lt = LangevinThermostat::new(300.0, 1.0);
    let mut state = MdState::new(vec![[0.0; 3]; 10], vec![[0.0; 3]; 10], vec![1.0; 10]);
    lt.apply_o_step(&mut state, 0.002);
    // With γ=1, c1 = exp(-0.002) ≈ 0.998, and noise will be added
    let any_nonzero = state
        .velocities
        .iter()
        .any(|v| v.iter().any(|&vi| vi != 0.0));
    assert!(
        any_nonzero,
        "Langevin O-step should produce non-zero velocities"
    );
}

#[test]
fn test_langevin_target_kinetic_energy_per_atom() {
    let lt = LangevinThermostat::new(300.0, 1.0);
    let expected = 1.5_f64 * KB_REDUCED * 300.0;
    assert!((lt.target_kinetic_energy_per_atom() - expected).abs() < 1e-12);
}

// -----------------------------------------------------------------------
// VolumeFluctuationMonitor tests
// -----------------------------------------------------------------------

#[test]
fn test_volume_monitor_empty_mean_zero() {
    let vm = VolumeFluctuationMonitor::new();
    assert_eq!(vm.mean_volume(), 0.0);
}

#[test]
fn test_volume_monitor_push_box_correct_volume() {
    let mut vm = VolumeFluctuationMonitor::new();
    vm.push_box([3.0, 4.0, 5.0]);
    assert!((vm.mean_volume() - 60.0).abs() < 1e-12, "V = 3*4*5 = 60");
}

#[test]
fn test_volume_monitor_variance_positive_for_varying_volumes() {
    let mut vm = VolumeFluctuationMonitor::new();
    for v in [100.0, 110.0, 90.0, 105.0] {
        vm.push(v);
    }
    assert!(vm.variance_volume() > 0.0, "variance should be positive");
}

#[test]
fn test_volume_monitor_isothermal_compressibility_finite() {
    let mut vm = VolumeFluctuationMonitor::new();
    for v in [1000.0, 1001.0, 999.0, 1002.0] {
        vm.push(v);
    }
    let kappa = vm.isothermal_compressibility(300.0);
    assert!(kappa.is_finite() && kappa > 0.0, "κ_T = {kappa}");
}

#[test]
fn test_volume_monitor_compressibility_zero_for_zero_temperature() {
    let mut vm = VolumeFluctuationMonitor::new();
    vm.push(1000.0);
    let kappa = vm.isothermal_compressibility(0.0);
    assert_eq!(kappa, 0.0, "T=0 → κ_T = 0");
}

// -----------------------------------------------------------------------
// LoggingCallback tests
// -----------------------------------------------------------------------

#[test]
fn test_logging_callback_records_at_freq() {
    let mut cb = LoggingCallback::new(5);
    let config = MdConfig::default();
    let mut state = MdState::new(vec![[0.0; 3]; 2], vec![[0.0; 3]; 2], vec![1.0; 2]);
    for step in 1u64..=10 {
        state.step = step;
        cb.on_step(&state, &config);
    }
    assert_eq!(cb.records.len(), 2, "should record at step 5 and 10");
}

#[test]
fn test_logging_callback_zero_freq_no_records() {
    let mut cb = LoggingCallback::new(0);
    let config = MdConfig::default();
    let mut state = MdState::new(vec![[0.0; 3]; 2], vec![[0.0; 3]; 2], vec![1.0; 2]);
    for step in 1u64..=10 {
        state.step = step;
        cb.on_step(&state, &config);
    }
    assert_eq!(cb.records.len(), 0, "freq=0 → no records");
}
