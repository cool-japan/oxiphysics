// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the constant-pressure drivers of [`MdSim`]
//! ([`MdSim::run_npt`] with the Berendsen barostat and
//! [`MdSim::run_npt_parrinello_rahman`]).
//!
//! The barostats must (1) see the *pair virial* of the interacting system, not a
//! virial fixed at zero, (2) obey the virial-theorem pressure
//! `P = (2·KE + W) / (3V)`, and (3) for Parrinello-Rahman carry a persistent box
//! velocity that produces the oscillation period `tau_p`.  The checks use closed
//! forms (ideal gas `P V = N k_B T`, linear box oscillator, derivative of the
//! potential with respect to the volume) as well as Lennard-Jones fluids that
//! must converge to the same density from different starting densities.

use super::stochastic::SimRng;
use super::tests_support::{
    ARGON_MASS, ShiftedForceLj, argon_fcc_fluid, fcc_lattice, reduced_density,
};
use super::*;

/// Zero forces and zero virial: the ideal gas.
fn ideal_gas_forces(positions: &[[f64; 3]], _box_lengths: &[f64; 3]) -> (Vec<[f64; 3]>, f64) {
    (vec![[0.0; 3]; positions.len()], 0.0)
}

/// An MdSim with `n` argon atoms at `temperature` (K), random positions and
/// exactly thermal velocities in a cube of edge `edge`.
fn ideal_gas_sim(n: usize, temperature: f64, edge: f64, dt: f64, steps: u64) -> MdSim {
    let masses = vec![ARGON_MASS; n];
    let mut velocities = maxwell_boltzmann_velocities(&masses, temperature, 99);
    rescale_to_temperature(&mut velocities, &masses, temperature);
    let mut rng = SimRng::new(5);
    let positions: Vec<[f64; 3]> = (0..n)
        .map(|_| {
            [
                rng.uniform() * edge,
                rng.uniform() * edge,
                rng.uniform() * edge,
            ]
        })
        .collect();
    let config = MdConfig {
        dt,
        n_steps: steps,
        temperature,
        pressure: 1000.0,
        box_lengths: [edge; 3],
        pbc: true,
        ensemble: Ensemble::NPT,
    };
    MdSim::new(config, MdState::new(positions, velocities, masses))
}

/// Volume at which `n` ideal-gas atoms at `temperature` (K) exert `pressure_bar`.
fn ideal_gas_volume(n: usize, temperature: f64, pressure_bar: f64) -> f64 {
    n as f64 * KB_REDUCED * temperature / (pressure_bar / KJ_MOL_PER_ANG3_TO_BAR)
}

#[test]
fn test_lj_helper_virial_is_minus_three_v_du_dv() {
    // The pair virial W = Σ r_ij · F_ij satisfies W / (3V) = −dU/dV for a uniform
    // dilation of box and coordinates.  This ties the virial convention of the
    // drivers (repulsive pair: +|F| r) to the thermodynamic pressure.
    let lj = ShiftedForceLj::argon(2.0);
    let (mut positions, box_lengths) = fcc_lattice(3, 1.65 * ARGON_SIGMA_LOCAL);
    let mut rng = SimRng::new(11);
    for p in positions.iter_mut() {
        for x in p.iter_mut() {
            *x += 0.4 * (rng.uniform() - 0.5);
        }
    }
    let result = lj.evaluate(&positions, &box_lengths);
    let h = 1.0e-5;
    let energy_at = |scale: f64| {
        let scaled: Vec<[f64; 3]> = positions
            .iter()
            .map(|p| [p[0] * scale, p[1] * scale, p[2] * scale])
            .collect();
        let scaled_box = [
            box_lengths[0] * scale,
            box_lengths[1] * scale,
            box_lengths[2] * scale,
        ];
        lj.evaluate(&scaled, &scaled_box).potential_energy
    };
    let du_dscale = (energy_at(1.0 + h) - energy_at(1.0 - h)) / (2.0 * h);
    // U(λ) with r → λ r: dU/dλ = −W.
    assert!(
        ((du_dscale + result.virial) / result.virial).abs() < 1e-6,
        "dU/dλ = {du_dscale}, −W = {}",
        -result.virial
    );
    // Forces sum to zero (Newton's third law) and match −dU/dx numerically.
    let total: [f64; 3] = result.forces.iter().fold([0.0; 3], |mut acc, f| {
        for a in 0..3 {
            acc[a] += f[a];
        }
        acc
    });
    assert!(total.iter().all(|c| c.abs() < 1e-9), "net force {total:?}");
    let displaced = |sign: f64| {
        let mut moved = positions.clone();
        moved[7][1] += sign * h;
        lj.evaluate(&moved, &box_lengths).potential_energy
    };
    let numeric = -(displaced(1.0) - displaced(-1.0)) / (2.0 * h);
    assert!(
        (numeric - result.forces[7][1]).abs() < 1e-6 * result.forces[7][1].abs().max(1.0),
        "numeric force {numeric}, analytic {}",
        result.forces[7][1]
    );
}

/// Alias keeping the line lengths of the tests short.
const ARGON_SIGMA_LOCAL: f64 = super::tests_support::ARGON_SIGMA;

#[test]
fn test_run_npt_passes_the_virial_to_the_berendsen_barostat() {
    // Atoms at rest (KE = 0), so the pressure of the first step is exactly
    // P = W / (3V) · 16605.4 bar.  The box is dilated by
    // mu = (1 + beta dt / tau_p (P − P0))^(1/3) in that very step.
    let box_edge = 10.0;
    let virial = 5000.0;
    let positions = vec![[2.0, 3.0, 4.0], [6.0, 5.0, 7.0]];
    let state = MdState::new(positions, vec![[0.0; 3]; 2], vec![ARGON_MASS; 2]);
    let config = MdConfig {
        dt: 0.002,
        n_steps: 1,
        temperature: 300.0,
        pressure: 1.0,
        box_lengths: [box_edge; 3],
        pbc: true,
        ensemble: Ensemble::NPT,
    };
    let mut sim = MdSim::new(config, state);
    let barostat = NptBarostat::new(1.0, 0.1, 4.5e-5);
    let forces =
        move |positions: &[[f64; 3]], _box: &[f64; 3]| (vec![[0.0; 3]; positions.len()], virial);
    sim.run_npt(forces, &barostat, 0.1, 0);

    let pressure = virial / (3.0 * box_edge.powi(3)) * KJ_MOL_PER_ANG3_TO_BAR;
    let mu = (1.0 + 4.5e-5 * 0.002 / 0.1 * (pressure - 1.0)).cbrt();
    assert!(
        mu > 1.0 + 1e-3,
        "the test needs a clearly visible dilation: {mu}"
    );
    for edge in sim.config.box_lengths {
        assert!(
            (edge - box_edge * mu).abs() < 1e-12,
            "edge {edge}, expected {}",
            box_edge * mu
        );
    }
    // Without the virial the same atoms would not move the box at all by more
    // than the (zero) kinetic pressure: P = 0 → mu < 1.
    let mut ideal = MdSim::new(
        MdConfig {
            dt: 0.002,
            n_steps: 1,
            temperature: 300.0,
            pressure: 1.0,
            box_lengths: [box_edge; 3],
            pbc: true,
            ensemble: Ensemble::NPT,
        },
        MdState::new(
            vec![[2.0, 3.0, 4.0], [6.0, 5.0, 7.0]],
            vec![[0.0; 3]; 2],
            vec![ARGON_MASS; 2],
        ),
    );
    ideal.run_npt(ideal_gas_forces, &barostat, 0.1, 0);
    assert!(
        ideal.config.box_lengths[0] < box_edge,
        "an ideal gas at rest contracts"
    );
}

#[test]
fn test_run_npt_parrinello_rahman_passes_the_virial_to_the_box_equation() {
    // Atoms at rest: P = W / (3V) · 16605.4 bar at the start of the step, the box
    // velocity becomes dt · a and the edge L + dt² a with the acceleration a of
    // the PR equation W_box L̈ = (V/L)(P − P0).
    let box_edge = 10.0;
    let virial = 5000.0;
    let dt = 0.002;
    let state = MdState::new(
        vec![[2.0, 3.0, 4.0], [6.0, 5.0, 7.0]],
        vec![[0.0; 3]; 2],
        vec![ARGON_MASS; 2],
    );
    let config = MdConfig {
        dt,
        n_steps: 1,
        temperature: 300.0,
        pressure: 1.0,
        box_lengths: [box_edge; 3],
        pbc: true,
        ensemble: Ensemble::NPT,
    };
    let mut sim = MdSim::new(config, state);
    let mut barostat = ParrinelloRahmanBarostat::new(1.0, 0.5, 4.5e-5);
    let forces =
        move |positions: &[[f64; 3]], _box: &[f64; 3]| (vec![[0.0; 3]; positions.len()], virial);
    sim.run_npt_parrinello_rahman(forces, &mut barostat, 0.1, 0);

    let pressure = virial / (3.0 * box_edge.powi(3)) * KJ_MOL_PER_ANG3_TO_BAR;
    let reference = ParrinelloRahmanBarostat::new(1.0, 0.5, 4.5e-5);
    let acceleration = reference.box_acceleration(pressure, &[box_edge; 3]);
    for (axis, &a) in acceleration.iter().enumerate() {
        assert!(a > 0.0, "high pressure pushes the box out");
        assert!(
            (barostat.box_velocity[axis] - dt * a).abs() < 1e-15,
            "box velocity {}, expected {}",
            barostat.box_velocity[axis],
            dt * a
        );
        let expected_edge = box_edge + dt * dt * a;
        assert!(
            (sim.config.box_lengths[axis] - expected_edge).abs() < 1e-14,
            "edge {}, expected {expected_edge}",
            sim.config.box_lengths[axis]
        );
    }
}

#[test]
fn test_parrinello_rahman_box_keeps_its_velocity_between_steps() {
    // A pressure pulse during one step sets the box in motion; afterwards, at the
    // target pressure, the box coasts with constant velocity (extended-Lagrangian
    // dynamics).  Berendsen scaling would stop immediately.
    let mut sim = MdSim::new(
        MdConfig {
            dt: 0.002,
            box_lengths: [20.0; 3],
            ..MdConfig::default()
        },
        MdState::new(
            vec![[1.0, 2.0, 3.0]],
            vec![[0.1, 0.0, 0.0]],
            vec![ARGON_MASS],
        ),
    );
    let mut pr = ParrinelloRahmanBarostat::new(1.0, 0.5, 4.5e-5);
    pr.advance(&mut sim, 1500.0, 0.002);
    let velocity = pr.box_velocity;
    assert!(
        velocity.iter().all(|&v| v > 0.0),
        "pulse must start the box: {velocity:?}"
    );
    let edge_after_pulse = sim.config.box_lengths[0];
    for _ in 0..10 {
        pr.advance(&mut sim, 1.0, 0.002);
    }
    assert_eq!(
        pr.box_velocity, velocity,
        "no force at P = P0: velocity persists"
    );
    let expected = edge_after_pulse + 10.0 * 0.002 * velocity[0];
    assert!(
        (sim.config.box_lengths[0] - expected).abs() < 1e-12,
        "edge {}, expected {expected}",
        sim.config.box_lengths[0]
    );
}

#[test]
fn test_parrinello_rahman_box_oscillates_with_period_tau_p() {
    // A medium of compressibility beta obeys P(V) = P0 − (V − V0) / (beta V0).  The
    // barostat mass W = 3 tau_p² L / (4 pi² beta) makes the box oscillate with the
    // period tau_p; this validates the mass formula and the unit conversions of
    // the box equation of motion.
    let beta_per_bar = 1.0e-4;
    let tau_p = 2.0;
    let p0 = 1.0;
    let l0 = 30.0_f64;
    let v0 = l0.powi(3);
    let dt = 2.0e-3;
    let mut sim = MdSim::new(
        MdConfig {
            dt,
            box_lengths: [1.004 * l0; 3],
            ..MdConfig::default()
        },
        MdState::new(vec![], vec![], vec![]),
    );
    let mut pr = ParrinelloRahmanBarostat::new(p0, tau_p, beta_per_bar);

    let mut crossings: Vec<f64> = Vec::new();
    let volume = |sim: &MdSim| sim.config.box_lengths.iter().product::<f64>();
    let mut previous = volume(&sim) - v0;
    let steps = 6000_usize;
    for step in 0..steps {
        let v = volume(&sim);
        let pressure = p0 - (v - v0) / (beta_per_bar * v0);
        pr.advance(&mut sim, pressure, dt);
        let now = volume(&sim) - v0;
        if previous > 0.0 && now <= 0.0 {
            crossings.push((step as f64 + previous / (previous - now)) * dt);
        }
        previous = now;
    }
    assert!(crossings.len() >= 3, "crossings {crossings:?}");
    for pair in crossings.windows(2) {
        let period = pair[1] - pair[0];
        assert!(
            ((period - tau_p) / tau_p).abs() < 0.02,
            "box period {period} ps, expected {tau_p} ps"
        );
    }
    // A cubic box stays cubic under isotropic coupling.
    let edges = sim.config.box_lengths;
    assert!((edges[0] - edges[1]).abs() < 1e-12 && (edges[1] - edges[2]).abs() < 1e-12);
}

#[test]
fn test_berendsen_barostat_scale_factor_follows_the_weak_coupling_law() {
    let barostat = NptBarostat::new(1.0, 0.5, 4.5e-5);
    let dt = 0.002;
    // At the target pressure nothing happens.
    assert!((barostat.scale_factor(1.0, dt) - 1.0).abs() < 1e-15);
    // μ³ = 1 − β dt / τ (P0 − P)
    let p = 20_000.0;
    let expected = (1.0 - 4.5e-5 * dt / 0.5 * (1.0 - p)).cbrt();
    assert!((barostat.scale_factor(p, dt) - expected).abs() < 1e-15);
    assert!(expected > 1.0, "excess pressure expands the box");
    assert!(
        barostat.scale_factor(-20_000.0, dt) < 1.0,
        "tension contracts it"
    );
    // A runaway pressure cannot mirror the box: μ³ is bounded below.
    let floor = barostat.scale_factor(-1.0e12, dt);
    assert!(floor > 0.0 && floor < 1.0, "floor {floor}");
}

#[test]
fn test_ideal_gas_volume_under_parrinello_rahman_equals_nkt_over_p() {
    let n = 64_usize;
    let temperature = 300.0;
    let target_bar = 1000.0;
    let equilibrium = ideal_gas_volume(n, temperature, target_bar);
    let start_edge = (2.2 * equilibrium).cbrt();
    let mut sim = ideal_gas_sim(n, temperature, start_edge, 0.004, 25_000);
    // An ideal gas has the compressibility 1/P.
    let mut pr = ParrinelloRahmanBarostat::new(target_bar, 1.0, 1.0 / target_bar);
    let records = sim.run_npt_parrinello_rahman(ideal_gas_forces, &mut pr, 0.2, 10);
    let tail = &records[records.len() / 2..];
    let mean_volume = tail.iter().map(|r| r.3.powi(3)).sum::<f64>() / tail.len() as f64;
    let mean_t = tail.iter().map(|r| r.2).sum::<f64>() / tail.len() as f64;
    assert!(
        ((mean_volume - equilibrium) / equilibrium).abs() < 0.05,
        "mean V = {mean_volume} Å³, N k_B T / P = {equilibrium} Å³"
    );
    assert!(
        ((mean_t - temperature) / temperature).abs() < 0.03,
        "T = {mean_t} K"
    );
}

#[test]
fn test_ideal_gas_volume_under_berendsen_barostat_equals_nkt_over_p() {
    let n = 64_usize;
    let temperature = 300.0;
    let target_bar = 1000.0;
    let equilibrium = ideal_gas_volume(n, temperature, target_bar);
    let start_edge = (2.2 * equilibrium).cbrt();
    let mut sim = ideal_gas_sim(n, temperature, start_edge, 0.004, 15_000);
    let barostat = NptBarostat::new(target_bar, 2.0, 1.0 / target_bar);
    let records = sim.run_npt(ideal_gas_forces, &barostat, 0.2, 10);
    let _ = records;
    let volume = sim.config.box_lengths.iter().product::<f64>();
    // Equilibrium of the relaxation dV/dt ∝ (P − P0): the instantaneous pressure
    // fluctuates by sqrt(2/(3N)) = 10 %.
    assert!(
        ((volume - equilibrium) / equilibrium).abs() < 0.15,
        "V = {volume} Å³, N k_B T / P = {equilibrium} Å³"
    );
}

/// Run an NPT trajectory of the 108-atom LJ fluid in chunks and return the time
/// series of `(density, pressure, temperature)` sampled after every chunk.
fn lj_npt_series(
    start_density: f64,
    parrinello_rahman: bool,
    chunks: usize,
) -> Vec<(f64, f64, f64)> {
    let temperature = 150.0;
    let target_bar = 500.0;
    let lj = ShiftedForceLj::argon(2.0);
    let (atoms, box_lengths) = argon_fcc_fluid(3, start_density, temperature, 31);
    let positions: Vec<[f64; 3]> = atoms.iter().map(|a| a.position).collect();
    let velocities: Vec<[f64; 3]> = atoms.iter().map(|a| a.velocity).collect();
    let masses = vec![ARGON_MASS; atoms.len()];
    let config = MdConfig {
        dt: 0.004,
        n_steps: 100,
        temperature,
        pressure: target_bar,
        box_lengths,
        pbc: true,
        ensemble: Ensemble::NPT,
    };
    let mut sim = MdSim::new(config, MdState::new(positions, velocities, masses));
    let forces = lj.virial_fn();
    let mut pr = ParrinelloRahmanBarostat::new(target_bar, 1.0, 1.0e-4);
    let berendsen = NptBarostat::new(target_bar, 2.0, 1.0e-4);
    let mut series = Vec::with_capacity(chunks);
    for _ in 0..chunks {
        if parrinello_rahman {
            sim.run_npt_parrinello_rahman(&forces, &mut pr, 0.5, 0);
        } else {
            sim.run_npt(&forces, &berendsen, 0.5, 0);
        }
        let (_, virial) = forces(&sim.state.positions, &sim.config.box_lengths);
        series.push((
            reduced_density(sim.state.n_atoms(), &sim.config.box_lengths),
            sim.compute_pressure(virial),
            sim.compute_temperature(),
        ));
    }
    series
}

fn mean_of(series: &[(f64, f64, f64)], pick: impl Fn(&(f64, f64, f64)) -> f64) -> f64 {
    series.iter().map(&pick).sum::<f64>() / series.len() as f64
}

#[test]
fn test_lj_fluid_under_parrinello_rahman_reaches_the_target_pressure_from_both_sides() {
    let ideal_density = 500.0 / KJ_MOL_PER_ANG3_TO_BAR / (150.0 * KB_REDUCED)
        * super::tests_support::ARGON_SIGMA.powi(3);
    let mut finals = Vec::new();
    for start in [0.85, 0.55] {
        let series = lj_npt_series(start, true, 200);
        let tail = &series[series.len() / 2..];
        let mean_p = mean_of(tail, |s| s.1);
        let mean_t = mean_of(tail, |s| s.2);
        let mean_rho = mean_of(tail, |s| s.0);
        assert!(
            ((mean_p - 500.0) / 500.0).abs() < 0.05,
            "start {start}: mean P = {mean_p} bar, target 500 bar"
        );
        assert!(
            ((mean_t - 150.0) / 150.0).abs() < 0.03,
            "start {start}: T = {mean_t} K"
        );
        finals.push(mean_rho);
        // The attractive pair virial lowers the density far below the value an
        // ideal gas (zero virial) would reach at this pressure and temperature.
        assert!(
            mean_rho < 0.75 * ideal_density,
            "start {start}: rho* = {mean_rho}, ideal-gas density {ideal_density}"
        );
    }
    assert!(
        ((finals[0] - finals[1]) / finals[1]).abs() < 0.05,
        "the two starting densities must relax to the same state: {finals:?}"
    );
}

#[test]
fn test_lj_fluid_under_berendsen_barostat_reaches_the_target_pressure_from_both_sides() {
    let mut finals = Vec::new();
    for start in [0.85, 0.55] {
        let series = lj_npt_series(start, false, 200);
        let tail = &series[series.len() / 2..];
        let mean_p = mean_of(tail, |s| s.1);
        let mean_t = mean_of(tail, |s| s.2);
        assert!(
            ((mean_p - 500.0) / 500.0).abs() < 0.10,
            "start {start}: mean P = {mean_p} bar, target 500 bar"
        );
        assert!(
            ((mean_t - 150.0) / 150.0).abs() < 0.03,
            "start {start}: T = {mean_t} K"
        );
        finals.push(mean_of(tail, |s| s.0));
    }
    assert!(
        ((finals[0] - finals[1]) / finals[1]).abs() < 0.06,
        "the two starting densities must relax to the same state: {finals:?}"
    );
}
