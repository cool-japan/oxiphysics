// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for [`super`] (bead and system containers, normal modes, staging
//! coordinates, thermostat coefficients, CMD/RPMD state, VACF, tunnelling and
//! diagnostics).  The estimator and integrator regression tests live in
//! `pimd_tests.rs`.

use super::*;

// ── QuantumBead ──────────────────────────────────────────────────────────

#[test]
fn test_bead_new() {
    let b = QuantumBead::new(4, 1.0, 1.0);
    assert_eq!(b.n_beads(), 4);
    assert_eq!(b.positions.len(), 4);
    assert_eq!(b.velocities.len(), 4);
}

#[test]
fn test_bead_positions_zero() {
    let b = QuantumBead::new(3, 2.0, 1.0);
    for pos in &b.positions {
        assert_eq!(*pos, [0.0, 0.0, 0.0]);
    }
}

#[test]
fn test_bead_mass() {
    let b = QuantumBead::new(2, 3.5, 1.0);
    assert!((b.mass - 3.5).abs() < 1e-12);
}

#[test]
fn test_bead_ring_spread_zero_for_identical_positions() {
    let mut b = QuantumBead::new(4, 1.0, 1.0);
    for p in &mut b.positions {
        *p = [1.0, 2.0, 3.0];
    }
    assert!(b.ring_spread_sq() < 1e-12);
}

#[test]
fn test_bead_ring_spread_positive_when_dispersed() {
    let mut b = QuantumBead::new(2, 1.0, 1.0);
    b.positions[0] = [0.0, 0.0, 0.0];
    b.positions[1] = [2.0, 0.0, 0.0];
    assert!(b.ring_spread_sq() > 0.0);
}

// ── QuantumSystem ────────────────────────────────────────────────────────

#[test]
fn test_system_new() {
    let s = QuantumSystem::new(3, 4, 1.0, 1.0, 1.0);
    assert_eq!(s.n_atoms, 3);
    assert_eq!(s.n_beads, 4);
    assert_eq!(s.beads.len(), 3);
}

#[test]
fn test_system_beta() {
    let s = QuantumSystem::new(1, 2, 2.5, 1.0, 1.0);
    assert!((s.beta - 2.5).abs() < 1e-12);
}

#[test]
fn test_system_temperature() {
    let s = QuantumSystem::new(1, 4, 2.0, 1.0, 1.0);
    assert!((s.temperature() - 0.5).abs() < 1e-12);
}

#[test]
fn test_system_centroid_ke_zero_for_zero_velocity() {
    let s = QuantumSystem::new(2, 4, 1.0, 1.0, 1.0);
    assert!(s.centroid_kinetic_energy() < 1e-12);
}

// ── ring_polymer_spring_constant ─────────────────────────────────────────

#[test]
fn test_spring_constant_positive() {
    let k = ring_polymer_spring_constant(1.0, 4, 1.0, 1.0);
    assert!(k > 0.0);
}

#[test]
fn test_spring_constant_scales_with_p_squared() {
    let k4 = ring_polymer_spring_constant(1.0, 4, 1.0, 1.0);
    let k8 = ring_polymer_spring_constant(1.0, 8, 1.0, 1.0);
    assert!((k8 - 4.0 * k4).abs() < 1e-10, "k4={k4}, k8={k8}");
}

#[test]
fn test_spring_constant_scales_with_mass() {
    let k1 = ring_polymer_spring_constant(1.0, 4, 1.0, 1.0);
    let k2 = ring_polymer_spring_constant(2.0, 4, 1.0, 1.0);
    assert!((k2 - 2.0 * k1).abs() < 1e-10);
}

#[test]
fn test_spring_constant_scales_with_inverse_hbar_squared() {
    let k1 = ring_polymer_spring_constant(1.5, 6, 0.8, 1.0);
    let k2 = ring_polymer_spring_constant(1.5, 6, 0.8, 2.0);
    assert!((k2 - 0.25 * k1).abs() < 1e-12 * k1, "k1={k1}, k2={k2}");
    // m (P / (β ℏ))² evaluated directly.
    assert!((k1 - 1.5 * (6.0_f64 / 0.8).powi(2)).abs() < 1e-12 * k1);
}

// ── ring_polymer_spring_force ────────────────────────────────────────────

#[test]
fn test_spring_force_uniform_ring_zero() {
    let beads = vec![[1.0_f64, 2.0, 3.0]; 4];
    let forces = ring_polymer_spring_force(&beads, 1.0);
    for f in forces {
        assert!(f[0].abs() < 1e-12);
        assert!(f[1].abs() < 1e-12);
        assert!(f[2].abs() < 1e-12);
    }
}

#[test]
fn test_spring_force_length() {
    let beads = vec![[0.0_f64; 3]; 5];
    let forces = ring_polymer_spring_force(&beads, 1.0);
    assert_eq!(forces.len(), 5);
}

#[test]
fn test_spring_force_empty() {
    let forces = ring_polymer_spring_force(&[], 1.0);
    assert!(forces.is_empty());
}

#[test]
fn test_spring_force_sum_zero() {
    let beads = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let forces = ring_polymer_spring_force(&beads, 2.0);
    let sum_x: f64 = forces.iter().map(|f| f[0]).sum();
    let sum_y: f64 = forces.iter().map(|f| f[1]).sum();
    let sum_z: f64 = forces.iter().map(|f| f[2]).sum();
    assert!(sum_x.abs() < 1e-10, "sum_x={sum_x}");
    assert!(sum_y.abs() < 1e-10, "sum_y={sum_y}");
    assert!(sum_z.abs() < 1e-10, "sum_z={sum_z}");
}

// ── centroid_position ────────────────────────────────────────────────────

#[test]
fn test_centroid_empty() {
    assert_eq!(centroid_position(&[]), [0.0, 0.0, 0.0]);
}

#[test]
fn test_centroid_single() {
    let c = centroid_position(&[[1.0, 2.0, 3.0]]);
    assert_eq!(c, [1.0, 2.0, 3.0]);
}

#[test]
fn test_centroid_two_beads() {
    let beads = vec![[0.0_f64, 0.0, 0.0], [2.0, 4.0, 6.0]];
    let c = centroid_position(&beads);
    assert!((c[0] - 1.0).abs() < 1e-12);
    assert!((c[1] - 2.0).abs() < 1e-12);
    assert!((c[2] - 3.0).abs() < 1e-12);
}

// ── centroid_velocity ────────────────────────────────────────────────────

#[test]
fn test_centroid_velocity() {
    let vels = vec![[1.0_f64, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    let cv = centroid_velocity(&vels);
    assert!(cv[0].abs() < 1e-12);
}

// ── quantum_kinetic_energy ───────────────────────────────────────────────

#[test]
fn test_qke_positive_high_temp() {
    let sys = QuantumSystem::new(2, 4, 0.01, 1.0, 1.0);
    let forces = vec![vec![[0.0; 3]; 4]; 2];
    let ke = quantum_kinetic_energy(&sys, &forces).expect("forces match the system");
    assert!(ke > 0.0);
    // Without forces the virial correction vanishes: pure classical 3N/(2β).
    assert!((ke - 3.0 * 2.0 / (2.0 * 0.01)).abs() < 1e-9);
}

#[test]
fn test_qke_single_bead_is_classical() {
    let sys = QuantumSystem::new(1, 1, 1.0, 1.0, 1.0);
    let forces = vec![vec![[0.0; 3]; 1]];
    let ke = quantum_kinetic_energy(&sys, &forces).expect("forces match the system");
    let expected = 3.0 * 1.0 * 1.0 / (2.0 * 1.0);
    assert!((ke - expected).abs() < 1e-10);
}

// ── primitive_kinetic_energy ─────────────────────────────────────────────

#[test]
fn test_primitive_ke_positive_high_temp() {
    let sys = QuantumSystem::new(1, 4, 0.01, 1.0, 1.0);
    let ke = primitive_kinetic_energy(&sys);
    assert!(ke > 0.0);
}

#[test]
fn test_primitive_ke_single_bead_matches_virial() {
    let sys = QuantumSystem::new(1, 1, 1.0, 1.0, 1.0);
    let prim = primitive_kinetic_energy(&sys);
    let forces = vec![vec![[0.0; 3]; 1]];
    let virial = quantum_kinetic_energy(&sys, &forces).expect("forces match the system");
    // For P=1 and all positions at origin, both estimators should agree
    assert!((prim - virial).abs() < 1e-10);
}

// ── thermal_de_broglie ───────────────────────────────────────────────────

#[test]
fn test_de_broglie_positive() {
    assert!(thermal_de_broglie(1.0, 1.0, 1.0) > 0.0);
}

#[test]
fn test_de_broglie_heavier_mass_smaller() {
    let l1 = thermal_de_broglie(1.0, 1.0, 1.0);
    let l2 = thermal_de_broglie(4.0, 1.0, 1.0);
    assert!(l1 > l2);
}

#[test]
fn test_de_broglie_formula() {
    let m = 2.0;
    let beta = 1.0;
    let hbar = 1.0;
    let expected = hbar * (2.0 * PI * beta / m).sqrt();
    assert!((thermal_de_broglie(m, beta, hbar) - expected).abs() < 1e-12);
}

#[test]
fn test_de_broglie_scales_with_hbar() {
    let l1 = thermal_de_broglie(1.0, 1.0, 1.0);
    let l2 = thermal_de_broglie(1.0, 1.0, 2.0);
    assert!((l2 - 2.0 * l1).abs() < 1e-12);
}

// ── quantum_correction_factor ────────────────────────────────────────────

#[test]
fn test_correction_equal() {
    assert!((quantum_correction_factor(1.0, 1.0) - 1.0).abs() < 1e-12);
}

#[test]
fn test_correction_zero_classical() {
    assert!((quantum_correction_factor(0.0, 5.0) - 1.0).abs() < 1e-12);
}

#[test]
fn test_correction_ratio() {
    let qc = quantum_correction_factor(2.0, 3.0);
    assert!((qc - 1.5).abs() < 1e-12);
}

// ── step_pimd ────────────────────────────────────────────────────────────

#[test]
fn test_step_pimd_zero_force_constant_velocity() {
    let mut sys = QuantumSystem::new(1, 2, 1.0, 1.0, 1.0);
    sys.beads[0].velocities[0] = [1.0, 0.0, 0.0];
    sys.beads[0].velocities[1] = [1.0, 0.0, 0.0];
    let mut forces: Vec<Vec<[f64; 3]>> = vec![vec![[0.0; 3]; 2]];
    let dt = 0.01;
    let pos_before = sys.beads[0].positions[0];
    step_pimd(&mut sys, &mut forces, dt, |_, f| {
        for bead_forces in f.iter_mut() {
            bead_forces.fill([0.0; 3]);
        }
    })
    .expect("force shape matches the system");
    let pos_after = sys.beads[0].positions[0];
    assert!((pos_after[0] - pos_before[0]).abs() > 1e-12);
    // Both beads move rigidly with the common velocity: a pure centroid drift.
    assert!((pos_after[0] - dt).abs() < 1e-14);
}

#[test]
fn test_step_pimd_momentum_change_with_force() {
    let mut sys = QuantumSystem::new(1, 1, 1.0, 1.0, 1.0);
    let mut forces: Vec<Vec<[f64; 3]>> = vec![vec![[1.0, 0.0, 0.0]]];
    let dt = 0.1;
    step_pimd(&mut sys, &mut forces, dt, |_, f| f[0][0] = [1.0, 0.0, 0.0])
        .expect("force shape matches the system");
    // Constant force F = 1, m = 1: v = F dt exactly.
    assert!((sys.beads[0].velocities[0][0] - 0.1).abs() < 1e-14);
}

// ── NormalModeTransform ──────────────────────────────────────────────────

#[test]
fn test_normal_mode_forward_inverse_roundtrip() {
    let p = 4;
    let nm = NormalModeTransform::new(p);
    let original: Vec<[f64; 3]> = (0..p)
        .map(|i| [i as f64, (i + 1) as f64 * 0.5, 0.0])
        .collect();
    let modes = nm.forward(&original);
    let recovered = nm.inverse(&modes);
    for k in 0..p {
        for d in 0..3 {
            assert!(
                (recovered[k][d] - original[k][d]).abs() < 1e-10,
                "k={k} d={d}: got {} expected {}",
                recovered[k][d],
                original[k][d]
            );
        }
    }
}

#[test]
fn test_normal_mode_centroid_is_mode_zero() {
    // The zero-th normal mode should equal the centroid (scaled by sqrt(P))
    let p = 4;
    let nm = NormalModeTransform::new(p);
    let beads: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; p];
    let modes = nm.forward(&beads);
    let centroid = centroid_position(&beads);
    // mode[0] = centroid * sqrt(P) (from the 1/sqrt(P) prefactor in C)
    let expected = v3_scale((p as f64).sqrt(), centroid);
    for d in 0..3 {
        assert!(
            (modes[0][d] - expected[d]).abs() < 1e-10,
            "d={d}: got {} expected {}",
            modes[0][d],
            expected[d]
        );
    }
}

// ── normal_mode_frequencies ──────────────────────────────────────────────

#[test]
fn test_normal_mode_freq_zero_for_centroid() {
    let freqs = normal_mode_frequencies(4, 1.0, 1.0);
    assert!(freqs[0].abs() < 1e-12, "centroid frequency should be zero");
}

#[test]
fn test_normal_mode_freq_positive_for_non_zero_modes() {
    let freqs = normal_mode_frequencies(4, 1.0, 1.0);
    for &f in &freqs[1..] {
        assert!(f >= 0.0);
    }
}

#[test]
fn test_normal_mode_freq_count() {
    let freqs = normal_mode_frequencies(8, 1.0, 1.0);
    assert_eq!(freqs.len(), 8);
}

// ── StagingCoords ────────────────────────────────────────────────────────

#[test]
fn test_staging_mass_bead_1_is_physical_mass() {
    let sc = StagingCoords::new(4, 2.0);
    assert!((sc.staging_mass(1) - 2.0).abs() < 1e-12);
}

#[test]
fn test_staging_mass_decreases_with_level() {
    // staging_mass(j) = j/(j-1) * m; the ratio decreases towards 1 as j grows
    let sc = StagingCoords::new(8, 1.0);
    let m2 = sc.staging_mass(2); // 2/1 = 2.0
    let m4 = sc.staging_mass(4); // 4/3 ≈ 1.333
    assert!(m2 > m4, "m2={m2} should exceed m4={m4}");
}

#[test]
fn test_staging_mass_formula() {
    let sc = StagingCoords::new(4, 1.0);
    // m_s(4) = 4/3 * m
    let expected = 4.0 / 3.0;
    assert!((sc.staging_mass(4) - expected).abs() < 1e-12);
}

// ── PileThermostat ───────────────────────────────────────────────────────

#[test]
fn test_pile_c1_less_than_one() {
    let pile = PileThermostat::new(4, 1.0, 1.0, 1.0, 0.01);
    for &c in &pile.c1 {
        assert!(c <= 1.0 && c > 0.0, "c1 should be in (0,1], got {c}");
    }
}

#[test]
fn test_pile_c2_non_negative() {
    let pile = PileThermostat::new(4, 1.0, 1.0, 1.0, 0.01);
    for &c in &pile.c2 {
        assert!(c >= 0.0, "c2 should be non-negative, got {c}");
    }
}

#[test]
fn test_pile_apply_half_step_preserves_order() {
    let pile = PileThermostat::new(4, 1.0, 1.0, 0.5, 0.001);
    // With xi=0 the new velocity is just c1 * v_old (damped)
    let v_new = pile.apply_half_step(0, 1.0, 1.0, 0.0);
    assert!(v_new < 1.0 && v_new > 0.0, "Damped velocity: {v_new}");
}

#[test]
fn test_pile_gamma_modes_count() {
    let pile = PileThermostat::new(6, 1.0, 1.0, 1.0, 0.01);
    assert_eq!(pile.gamma_modes.len(), 6);
}

// ── CentroidMdState ──────────────────────────────────────────────────────

#[test]
fn test_cmd_state_creation() {
    let state = CentroidMdState::new(3, 4, 1.0, 1.0);
    assert_eq!(state.centroid_positions.len(), 3);
    assert_eq!(state.centroid_velocities.len(), 3);
}

#[test]
fn test_cmd_temperature_zero_for_zero_velocities() {
    let state = CentroidMdState::new(4, 4, 1.0, 1.0);
    assert!(state.temperature() < 1e-12);
}

#[test]
fn test_cmd_step_moves_position() {
    let mut state = CentroidMdState::new(1, 4, 1.0, 1.0);
    state.centroid_velocities[0] = [1.0, 0.0, 0.0];
    state.step(0.1);
    assert!((state.centroid_positions[0][0] - 0.1_f64).abs() < 1e-10);
}

#[test]
fn test_cmd_kinetic_energy_formula() {
    let mut state = CentroidMdState::new(1, 4, 1.0, 2.0);
    state.centroid_velocities[0] = [1.0, 0.0, 0.0];
    // KE = 0.5 * m * v^2 = 0.5 * 2 * 1 = 1.0
    assert!((state.kinetic_energy() - 1.0).abs() < 1e-12);
}

// ── RingPolymerMdState ───────────────────────────────────────────────────

#[test]
fn test_rpmd_creation() {
    let state = RingPolymerMdState::new(2, 4, 1.0, 1.0, 1.0);
    assert_eq!(state.system.n_atoms, 2);
    assert_eq!(state.system.n_beads, 4);
    assert!((state.time - 0.0).abs() < 1e-12);
}

#[test]
fn test_rpmd_step_advances_time() {
    let mut state = RingPolymerMdState::new(1, 2, 1.0, 1.0, 1.0);
    state
        .step(0.05, |_, f| {
            for bead_forces in f.iter_mut() {
                bead_forces.fill([0.0; 3]);
            }
        })
        .expect("force shape matches the system");
    assert!((state.time - 0.05).abs() < 1e-12);
}

#[test]
fn test_rpmd_spring_energy_zero_for_uniform_ring() {
    let state = RingPolymerMdState::new(1, 4, 1.0, 1.0, 1.0);
    // All beads at origin → spring energy = 0
    assert!(state.spring_energy() < 1e-12);
}

// ── VacfAccumulator ──────────────────────────────────────────────────────

#[test]
fn test_vacf_tau_zero_is_one() {
    let mut acc = VacfAccumulator::new(10);
    acc.push(vec![[1.0, 0.0, 0.0]]);
    acc.push(vec![[0.5, 0.0, 0.0]]);
    let c0 = acc.vacf(0).unwrap();
    assert!((c0 - 1.0).abs() < 1e-10, "VACF(0) should be 1, got {c0}");
}

#[test]
fn test_vacf_insufficient_data_returns_none() {
    let acc = VacfAccumulator::new(10);
    assert!(acc.vacf(5).is_none());
}

#[test]
fn test_vacf_unnorm_zero_for_empty() {
    let acc = VacfAccumulator::new(5);
    assert!((acc.vacf_unnorm(0) - 0.0).abs() < 1e-12);
}

#[test]
fn test_vacf_max_snapshots_respected() {
    let mut acc = VacfAccumulator::new(3);
    for i in 0..5 {
        acc.push(vec![[i as f64, 0.0, 0.0]]);
    }
    assert_eq!(acc.snapshots.len(), 3);
}

// ── quantum_diffusion_coefficient ────────────────────────────────────────

#[test]
fn test_diffusion_zero_for_zero_vacf() {
    let vacf = vec![0.0_f64; 10];
    let d = quantum_diffusion_coefficient(&vacf, 0.01);
    assert!(d.abs() < 1e-12);
}

#[test]
fn test_diffusion_positive_for_positive_vacf() {
    let vacf = vec![1.0_f64; 20];
    let d = quantum_diffusion_coefficient(&vacf, 0.01);
    assert!(d > 0.0);
}

#[test]
fn test_diffusion_empty_vacf() {
    let d = quantum_diffusion_coefficient(&[], 0.01);
    assert!((d - 0.0).abs() < 1e-12);
}

// ── TunnelingRate ────────────────────────────────────────────────────────

#[test]
fn test_wkb_transmission_above_barrier_is_one() {
    let tr = TunnelingRate::new(1.0, 1.0, 1.0, 1.0);
    assert!((tr.wkb_transmission(2.0) - 1.0).abs() < 1e-12);
}

#[test]
fn test_wkb_transmission_below_barrier_less_than_one() {
    let tr = TunnelingRate::new(1.0, 1.0, 1.0, 1.0);
    assert!(tr.wkb_transmission(0.5) < 1.0);
}

#[test]
fn test_wkb_transmission_positive() {
    let tr = TunnelingRate::new(5.0, 1.0, 1.0, 1.0);
    assert!(tr.wkb_transmission(1.0) > 0.0);
}

#[test]
fn test_instanton_rate_positive() {
    let tr = TunnelingRate::new(1.0, 0.5, 1.0, 1.0);
    assert!(tr.instanton_rate() > 0.0);
}

#[test]
fn test_wkb_tunneling_rate_function_positive() {
    let rate = wkb_tunneling_rate(1.0, 1.0, 1.0, 1.0);
    assert!(rate > 0.0);
}

// ── zero_point_energy ────────────────────────────────────────────────────

#[test]
fn test_zpe_formula() {
    let zpe = zero_point_energy(2.0, 1.0);
    assert!((zpe - 1.0).abs() < 1e-12);
}

#[test]
fn test_zpe_positive() {
    assert!(zero_point_energy(1.0, 1.0) > 0.0);
}

// ── path_integral_partition_ratio ────────────────────────────────────────

#[test]
fn test_partition_ratio_approaches_one_large_p() {
    // For large P, quantum discretisation error → 0 so ratio → 1
    let r = path_integral_partition_ratio(1.0, 1.0, 1.0, 64);
    assert!((r - 1.0).abs() < 0.05, "ratio={r}");
}

#[test]
fn test_partition_ratio_positive() {
    let r = path_integral_partition_ratio(1.0, 2.0, 1.0, 8);
    assert!(r > 0.0);
}

// ── ring_polymer_momentum_spread ─────────────────────────────────────────

#[test]
fn test_momentum_spread_positive() {
    let sigma = ring_polymer_momentum_spread(1.0, 4, 1.0, 1.0);
    assert!(sigma > 0.0);
}

#[test]
fn test_momentum_spread_scales_with_p() {
    let s4 = ring_polymer_momentum_spread(1.0, 4, 1.0, 1.0);
    let s16 = ring_polymer_momentum_spread(1.0, 16, 1.0, 1.0);
    assert!(s16 > s4, "More beads → larger momentum spread");
}

// ── harmonic_kubo_correlation ────────────────────────────────────────────

#[test]
fn test_kubo_correlation_at_t_zero_positive() {
    let c0 = harmonic_kubo_correlation(1.0, 0.0, 1.0, 1.0, 1.0);
    assert!(c0 > 0.0);
}

#[test]
fn test_kubo_correlation_finite() {
    let c = harmonic_kubo_correlation(1.0, 0.5, 1.0, 1.0, 1.0);
    assert!(c.is_finite());
}

// ── QuantumDiagnostics ───────────────────────────────────────────────────

#[test]
fn test_diagnostics_summary_not_empty() {
    let sys = QuantumSystem::new(2, 4, 1.0, 1.0, 1.0);
    let forces = vec![vec![[0.0; 3]; 4]; 2];
    let diag = QuantumDiagnostics::from_system(&sys, &forces).expect("forces match");
    assert!(!diag.summary().is_empty());
}

#[test]
fn test_diagnostics_de_broglie_positive() {
    let sys = QuantumSystem::new(1, 4, 1.0, 1.0, 1.0);
    let forces = vec![vec![[0.0; 3]; 4]];
    let diag = QuantumDiagnostics::from_system(&sys, &forces).expect("forces match");
    assert!(diag.de_broglie > 0.0);
}

#[test]
fn test_diagnostics_temperature_matches_beta() {
    let sys = QuantumSystem::new(1, 4, 2.0, 1.0, 1.0);
    let forces = vec![vec![[0.0; 3]; 4]];
    let diag = QuantumDiagnostics::from_system(&sys, &forces).expect("forces match");
    assert!((diag.temperature - 0.5).abs() < 1e-12);
}

// ── normal_mode_spring_energy ────────────────────────────────────────────

#[test]
fn test_normal_mode_spring_energy_zero_for_zero_amplitudes() {
    let amps = vec![[0.0_f64; 3]; 4];
    let freqs = vec![0.0, 1.0, 2.0, 1.0];
    let e = normal_mode_spring_energy(&amps, &freqs, 1.0);
    assert!(e.abs() < 1e-12);
}

#[test]
fn test_normal_mode_spring_energy_positive_for_nonzero_modes() {
    let amps = vec![[0.0_f64; 3], [1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]];
    let freqs = vec![0.0, 2.0, 3.0, 2.0];
    let e = normal_mode_spring_energy(&amps, &freqs, 1.0);
    // E = 0.5 * 1.0 * 4.0 * 1.0 = 2.0
    assert!((e - 2.0).abs() < 1e-12, "energy={e}");
}
