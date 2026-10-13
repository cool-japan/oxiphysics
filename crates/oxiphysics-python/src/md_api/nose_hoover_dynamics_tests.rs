// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the Nose-Hoover NVT dynamics of
//! [`PyNvtSimulation`](super::PyNvtSimulation).

use super::*;

/// Deterministic generator of uniform numbers in `[0, 1)` (64-bit LCG).
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    /// Approximately standard-normal number (Irwin-Hall sum of twelve uniforms).
    fn next_gaussian(&mut self) -> f64 {
        (0..12).map(|_| self.next_unit()).sum::<f64>() - 6.0
    }
}

/// 64 Lennard-Jones atoms on a jittered 4x4x4 lattice in a box of edge 6
/// with Maxwell-like velocities at `initial_temperature` and zero drift.
fn lj_gas(initial_temperature: f64, target: f64, tau: f64) -> PyNvtSimulation {
    let cfg = PyMdConfig {
        box_size: 6.0,
        target_temperature: Some(target),
        ..PyMdConfig::default()
    };
    let mut sim = PyNvtSimulation::new(cfg, tau);
    let mut rng = Lcg(2024);
    let spacing = 1.5;
    for ix in 0..4 {
        for iy in 0..4 {
            for iz in 0..4 {
                let jitter = |rng: &mut Lcg| 0.05 * (2.0 * rng.next_unit() - 1.0);
                let p = vec![
                    (f64::from(ix) + 0.5) * spacing + jitter(&mut rng),
                    (f64::from(iy) + 0.5) * spacing + jitter(&mut rng),
                    (f64::from(iz) + 0.5) * spacing + jitter(&mut rng),
                ];
                sim.add_atom(p, 0);
            }
        }
    }
    let n = sim.atom_count();
    let mut velocities: Vec<[f64; 3]> = (0..n)
        .map(|_| {
            [
                rng.next_gaussian(),
                rng.next_gaussian(),
                rng.next_gaussian(),
            ]
        })
        .collect();
    let mut drift = [0.0; 3];
    for v in &velocities {
        for k in 0..3 {
            drift[k] += v[k] / n as f64;
        }
    }
    for v in &mut velocities {
        for k in 0..3 {
            v[k] -= drift[k];
        }
    }
    let kinetic: f64 = velocities
        .iter()
        .map(|v| 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
        .sum();
    let scale = (initial_temperature / (2.0 * kinetic / (3.0 * n as f64))).sqrt();
    for (i, v) in velocities.iter().enumerate() {
        sim.md
            .set_velocity(i, vec![v[0] * scale, v[1] * scale, v[2] * scale]);
    }
    sim
}

#[test]
fn test_nvt_new_switches_the_inner_rescaling_thermostat_off() {
    let sim = PyNvtSimulation::new(PyMdConfig::default(), 0.1);
    assert!(
        !sim.md.thermostat_active(),
        "Nose-Hoover must be the only thermostat"
    );
}

#[test]
fn test_nvt_thermostat_mass_follows_tau_and_atom_count() {
    let mut sim = PyNvtSimulation::new(PyMdConfig::default(), 0.5);
    for k in 0..10 {
        sim.add_atom(vec![f64::from(k), 0.0, 0.0], 0);
    }
    // Q = 3 N T tau^2 = 30 * 1.2 * 0.25
    assert_eq!(sim.thermostat.n_dof, 30);
    assert!((sim.thermostat.thermostat_mass - 9.0).abs() < 1e-12);
    let mut slow = PyNvtSimulation::new(PyMdConfig::default(), 1.0);
    for k in 0..10 {
        slow.add_atom(vec![f64::from(k), 0.0, 0.0], 0);
    }
    assert!((slow.thermostat.thermostat_mass - 36.0).abs() < 1e-12);
}

#[test]
fn test_thermostat_energy_and_kick() {
    let mut nh = NoseHooverThermostat::new(2.0, 2, 0.5);
    nh.thermostat_mass = 2.0;
    nh.xi = 3.0;
    nh.eta = 0.5;
    // 1/2 Q xi^2 + g T eta = 9 + 6 * 2 * 0.5
    assert!((nh.thermostat_energy() - 15.0).abs() < 1e-12);

    let mut a = NoseHooverThermostat::new(1.0, 3, 0.2);
    let mut b = a.clone();
    let scale = a.half_step_xi(7.0, 0.01);
    b.kick_xi(7.0, 0.01);
    assert!((a.xi - b.xi).abs() < 1e-15);
    assert!((scale - b.damping_factor(0.01)).abs() < 1e-15);
    // KE above (g T / 2) heats the bath up: xi grows, velocities are damped.
    assert!(b.xi > 0.0 && scale < 1.0);
}

#[test]
fn test_nvt_step_applies_the_second_half_of_the_thermostat() {
    // With xi already non-zero a single step damps the velocities twice,
    // before and after the Verlet step. The two atoms are far apart and feel
    // no force, and the huge thermostat mass freezes xi, so the speed after
    // one step is exactly v * exp(-xi dt / 2)^2 = v * exp(-xi dt).
    let mut sim = PyNvtSimulation::new(PyMdConfig::default(), 0.1);
    sim.add_atom(vec![1.0, 1.0, 1.0], 0);
    sim.add_atom(vec![6.0, 6.0, 6.0], 0);
    sim.md.set_velocity(0, vec![1.0, 0.0, 0.0]);
    sim.md.set_velocity(1, vec![-1.0, 0.0, 0.0]);
    sim.thermostat.thermostat_mass = 1.0e12;
    sim.thermostat.xi = 0.4;
    let dt = 0.01;
    sim.step(dt);
    let expected = (-0.4_f64 * dt).exp();
    let v0 = sim.md.velocity(0).expect("atom 0")[0];
    assert!(
        (v0 - expected).abs() < 1e-9,
        "v = {v0}, expected exp(-xi dt) = {expected}"
    );
}

#[test]
fn test_nvt_samples_the_target_temperature_and_conserves_h() {
    let target = 1.2;
    let mut sim = lj_gas(2.0, target, 0.1);
    let dt = 0.002;
    let h0 = sim.conserved_energy();
    let mut temperatures = Vec::new();
    let mut max_drift = 0.0f64;
    let mut xi_extremes = (f64::INFINITY, f64::NEG_INFINITY);
    for step in 0..3000 {
        sim.step(dt);
        max_drift = max_drift.max(((sim.conserved_energy() - h0) / h0).abs());
        xi_extremes = (
            xi_extremes.0.min(sim.thermostat.xi),
            xi_extremes.1.max(sim.thermostat.xi),
        );
        if step >= 1000 {
            temperatures.push(sim.temperature());
        }
    }
    let count = temperatures.len() as f64;
    let mean = temperatures.iter().sum::<f64>() / count;
    let variance = temperatures.iter().map(|t| (t - mean).powi(2)).sum::<f64>() / count;
    let relative_spread = variance.sqrt() / mean;
    assert!(
        (mean - target).abs() < 0.05 * target,
        "mean temperature {mean} should be within 5 % of {target}"
    );
    assert!(
        relative_spread > 0.02,
        "a canonical thermostat lets T fluctuate (relative spread {relative_spread}); \
         a hard rescale would pin it"
    );
    assert!(
        xi_extremes.0 < 0.0 && xi_extremes.1 > 0.0,
        "the bath must both heat and cool: xi in [{}, {}]",
        xi_extremes.0,
        xi_extremes.1
    );
    assert!(
        max_drift < 1e-3,
        "conserved quantity drifted by {max_drift} (relative)"
    );
}

#[test]
fn test_conserved_energy_includes_the_bath_terms() {
    let mut sim = lj_gas(1.2, 1.2, 0.1);
    let h0 = sim.conserved_energy();
    let q = sim.thermostat.thermostat_mass;
    sim.thermostat.xi = 1.0;
    let with_momentum = sim.conserved_energy();
    assert!(
        (with_momentum - h0 - 0.5 * q).abs() < 1e-9,
        "the bath momentum adds 1/2 Q xi^2"
    );
    sim.thermostat.eta = 0.25;
    let g_t = sim.thermostat.n_dof as f64 * sim.thermostat.target_temperature;
    assert!(
        (sim.conserved_energy() - with_momentum - 0.25 * g_t).abs() < 1e-9,
        "the bath coordinate adds g k_B T eta"
    );
}
