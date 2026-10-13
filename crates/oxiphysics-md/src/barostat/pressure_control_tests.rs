//! Regression tests for pressure control.
//!
//! Covers the translation-invariant pair virial under periodic boundaries, the
//! inertia of the Parrinello-Rahman box, the explicit target-pressure argument
//! of the barostats and the Monte Carlo volume acceptance criterion.

use super::functions::{Barostat, compute_virial, virial_pressure};
use super::types::{BerendsenBarostat, MonteCarloBarostat, ParrinelloRahmanBarostat};
use crate::atom::AtomSet;
use crate::forcefield::{ForceField, PairForceField};
use crate::integrator::VelocityVerlet;
use crate::neighbor::PeriodicBox;
use crate::potential::{LennardJones, Potential};
use crate::simulation::MdSimulation;
use crate::thermostat::NoThermostat;
use oxiphysics_core::math::Vec3;

/// Edge length of the periodic test cell.
const BOX_EDGE: f64 = 10.0;

fn lj_field() -> PairForceField {
    let mut ff = PairForceField::new();
    ff.add_interaction(0, 0, LennardJones::new(1.0, 1.0, 2.5));
    ff
}

/// Two atoms on the x axis, mass 1, at rest, with the given x coordinates.
fn atom_pair(x_a: f64, x_b: f64) -> AtomSet {
    let mut atoms = AtomSet::new();
    atoms.add_atom(Vec3::new(x_a, 5.0, 5.0), Vec3::zeros(), 1.0, 0.0, 0);
    atoms.add_atom(Vec3::new(x_b, 5.0, 5.0), Vec3::zeros(), 1.0, 0.0, 0);
    atoms
}

/// Pair deep inside the cell, separated by 0.8.
fn interior_pair() -> AtomSet {
    atom_pair(4.6, 5.4)
}

/// The same pair translated so that it straddles the periodic boundary
/// (minimum-image separation again 0.8).
fn straddling_pair() -> AtomSet {
    atom_pair(9.6, 0.4)
}

/// Small deterministic generator of uniform numbers in `[0, 1)`.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

// ---------------------------------------------------------------------------
// Pair virial under periodic boundaries (F3314)
// ---------------------------------------------------------------------------

#[test]
fn pair_virial_equals_f_times_r_for_interior_pair() {
    let ff = lj_field();
    let pbox = PeriodicBox::cubic(BOX_EDGE);
    let mut atoms = interior_pair();
    let (energy, virial) = ff.compute_forces_virial(&mut atoms, &pbox);
    let lj = LennardJones::new(1.0, 1.0, 2.5);
    let expected = lj.force(0.8) * 0.8;
    assert!(expected > 100.0, "pair must be strongly repulsive");
    assert!(
        (virial - expected).abs() <= 1e-12 * expected,
        "virial {virial} vs f*r {expected}"
    );
    assert!((energy - lj.energy(0.8)).abs() < 1e-12);
}

#[test]
fn pair_virial_is_translation_invariant_across_the_boundary() {
    let ff = lj_field();
    let pbox = PeriodicBox::cubic(BOX_EDGE);
    let mut interior = interior_pair();
    let mut straddling = straddling_pair();
    let (e_in, w_in) = ff.compute_forces_virial(&mut interior, &pbox);
    let (e_st, w_st) = ff.compute_forces_virial(&mut straddling, &pbox);
    assert!((e_in - e_st).abs() < 1e-12, "energies {e_in} vs {e_st}");
    assert!(
        (w_in - w_st).abs() <= 1e-12 * w_in.abs(),
        "pair virial must not depend on wrapping: {w_in} vs {w_st}"
    );
}

#[test]
fn atomic_virial_of_wrapped_coordinates_is_not_translation_invariant() {
    // Documents why `Σ r_i·F_i` cannot be used under PBC: the straddling pair
    // picks up the spurious term L·F_ij = BOX_EDGE · f.
    let ff = lj_field();
    let pbox = PeriodicBox::cubic(BOX_EDGE);
    let mut interior = interior_pair();
    let mut straddling = straddling_pair();
    let (_, w_pair) = ff.compute_forces_virial(&mut interior, &pbox);
    let _ = ff.compute_forces_virial(&mut straddling, &pbox);
    let atomic = |atoms: &AtomSet| {
        let positions: Vec<[f64; 3]> = atoms.positions.iter().map(|p| [p.x, p.y, p.z]).collect();
        let forces: Vec<[f64; 3]> = atoms.forces.iter().map(|f| [f.x, f.y, f.z]).collect();
        compute_virial(&positions, &forces)
    };
    let f_mag = LennardJones::new(1.0, 1.0, 2.5).force(0.8);
    let atomic_interior = atomic(&interior);
    let atomic_straddling = atomic(&straddling);
    assert!((atomic_interior - w_pair).abs() <= 1e-9 * w_pair.abs());
    assert!(
        ((atomic_interior - atomic_straddling) - BOX_EDGE * f_mag).abs() <= 1e-9 * f_mag,
        "spurious boundary term: {atomic_interior} vs {atomic_straddling}"
    );
}

#[test]
fn pressure_of_boundary_pair_equals_pressure_of_interior_pair() {
    let ff = lj_field();
    let pbox = PeriodicBox::cubic(BOX_EDGE);
    let mut interior = interior_pair();
    let mut straddling = straddling_pair();
    let (_, w_in) = ff.compute_forces_virial(&mut interior, &pbox);
    let (_, w_st) = ff.compute_forces_virial(&mut straddling, &pbox);
    let volume = pbox.volume();
    let p_in = virial_pressure(volume, interior.kinetic_energy(), w_in);
    let p_st = virial_pressure(volume, straddling.kinetic_energy(), w_st);
    assert!(p_in > 0.0);
    assert!((p_in - p_st).abs() <= 1e-12 * p_in, "{p_in} vs {p_st}");
}

#[test]
fn berendsen_scaling_is_identical_for_interior_and_boundary_pair() {
    let ff = lj_field();
    let target = 0.0;
    let mut results = Vec::new();
    for mut atoms in [interior_pair(), straddling_pair()] {
        let mut pbox = PeriodicBox::cubic(BOX_EDGE);
        let (_, virial) = ff.compute_forces_virial(&mut atoms, &pbox);
        let mut baro = BerendsenBarostat::new(target, 1.0, 1e-3);
        Barostat::apply(&mut baro, &mut atoms, &mut pbox, target, virial, 0.01);
        results.push(pbox.dims.x);
    }
    assert!(
        results[0] > BOX_EDGE,
        "an over-pressured system must expand, got {}",
        results[0]
    );
    assert!(
        (results[0] - results[1]).abs() < 1e-12,
        "box edge {} vs {}",
        results[0],
        results[1]
    );
}

#[test]
fn md_simulation_reports_the_pair_virial_and_pressure() {
    let lj = LennardJones::new(1.0, 1.0, 2.5);
    let expected_virial = lj.force(0.8) * 0.8;
    let mut pressures = Vec::new();
    for atoms in [interior_pair(), straddling_pair()] {
        let mut sim = MdSimulation::new(
            atoms,
            lj_field(),
            VelocityVerlet::new(),
            NoThermostat::new(),
            BerendsenBarostat::new(1.0, 1.0, 1e-3),
            PeriodicBox::cubic(BOX_EDGE),
            0.0,
            1.0,
            1.0,
        );
        // dt = 0 evaluates forces and the virial without moving anything.
        sim.step(0.0);
        assert!(
            (sim.virial - expected_virial).abs() <= 1e-9 * expected_virial,
            "stored virial {} vs {expected_virial}",
            sim.virial
        );
        pressures.push(sim.pressure());
    }
    let expected_pressure = expected_virial / (3.0 * BOX_EDGE.powi(3));
    for p in &pressures {
        assert!(
            (p - expected_pressure).abs() <= 1e-9 * expected_pressure,
            "pressure {p} vs {expected_pressure}"
        );
    }
}

#[test]
fn md_simulation_run_primes_the_virial_before_the_first_step() {
    let mut sim = MdSimulation::new(
        straddling_pair(),
        lj_field(),
        VelocityVerlet::new(),
        NoThermostat::new(),
        crate::barostat::NoBarostat::new(),
        PeriodicBox::cubic(BOX_EDGE),
        0.0,
        0.0,
        1.0,
    );
    assert_eq!(sim.virial, 0.0);
    // Zero steps: only the initial force evaluation of `run` takes place.
    let records = sim.run(0, 0.01, None);
    assert!(records.is_empty());
    let expected = LennardJones::new(1.0, 1.0, 2.5).force(0.8) * 0.8;
    assert!(
        (sim.virial - expected).abs() <= 1e-9 * expected,
        "virial primed by run: {} vs {expected}",
        sim.virial
    );
}

// ---------------------------------------------------------------------------
// Explicit target pressure (F3314 / F3315)
// ---------------------------------------------------------------------------

#[test]
fn berendsen_uses_the_target_pressure_argument() {
    let mut atoms = interior_pair();
    let volume = BOX_EDGE.powi(3);
    let virial = 3.0 * volume * 2.0; // P = 2.0 for a gas at rest
    // The barostat's own target (5.0) differs from the argument (= the
    // current pressure), so only the argument leaves the box unchanged.
    let mut baro = BerendsenBarostat::new(5.0, 1.0, 1e-2);
    let mut pbox = PeriodicBox::cubic(BOX_EDGE);
    Barostat::apply(&mut baro, &mut atoms, &mut pbox, 2.0, virial, 0.1);
    assert!((pbox.dims.x - BOX_EDGE).abs() < 1e-12, "{}", pbox.dims.x);

    // A non-finite target falls back to the barostat's own setting (5.0 > P).
    let mut pbox = PeriodicBox::cubic(BOX_EDGE);
    Barostat::apply(&mut baro, &mut atoms, &mut pbox, f64::NAN, virial, 0.1);
    let mu = baro.scale_factor(2.0, 0.1);
    assert!(mu < 1.0);
    assert!((pbox.dims.x - BOX_EDGE * mu).abs() < 1e-12);
}

#[test]
fn berendsen_scale_factor_to_matches_scale_factor_for_own_target() {
    let baro = BerendsenBarostat::new(1.5, 0.5, 4.5e-5);
    for p in [0.0, 1.0, 1.5, 3.0] {
        assert_eq!(
            baro.scale_factor(p, 0.002),
            baro.scale_factor_to(p, 1.5, 0.002)
        );
    }
}

// ---------------------------------------------------------------------------
// Parrinello-Rahman box inertia (F3315)
// ---------------------------------------------------------------------------

/// Pressure excess `P - P0` used by the constant-mismatch PR tests.
const DELTA_P: f64 = 0.01;
const PR_MASS: f64 = 1.0e4;
const PR_DT: f64 = 0.1;

/// Virial that makes a gas at rest in a volume-1000 box exert `P0 + DELTA_P`.
fn pr_virial() -> f64 {
    3.0 * BOX_EDGE.powi(3) * (1.0 + DELTA_P)
}

#[test]
fn parrinello_rahman_box_velocity_persists_between_calls() {
    let mut atoms = interior_pair();
    let mut pbox = PeriodicBox::cubic(BOX_EDGE);
    let mut baro = ParrinelloRahmanBarostat::new(1.0, 1.0, PR_MASS);
    assert_eq!(baro.box_velocity, 0.0);

    // Box acceleration a = 3 V (P - P0) / W for a constant pressure mismatch.
    let accel = 3.0 * BOX_EDGE.powi(3) * DELTA_P / PR_MASS;
    baro.apply(&mut atoms, &mut pbox, 1.0, pr_virial(), PR_DT);
    let l1 = pbox.dims.x - BOX_EDGE;
    assert!(
        (l1 - accel * PR_DT * PR_DT).abs() <= 1e-3 * accel * PR_DT * PR_DT,
        "first call moves L by a dt^2: {l1}"
    );
    baro.apply(&mut atoms, &mut pbox, 1.0, pr_virial(), PR_DT);
    let l2 = pbox.dims.x - BOX_EDGE;
    // With inertia the displacement after two calls is a dt^2 (1 + 2) = 3 a dt^2;
    // the old memoryless rescaler gave 2 a dt^2.
    let inertial = 3.0 * accel * PR_DT * PR_DT;
    let memoryless = 2.0 * accel * PR_DT * PR_DT;
    assert!(
        (l2 - inertial).abs() <= 1e-2 * inertial,
        "two calls move L by 3 a dt^2 = {inertial}, got {l2}"
    );
    assert!((l2 - memoryless).abs() > 0.2 * memoryless);
    assert!(
        (baro.box_velocity - 2.0 * accel * PR_DT).abs() <= 1e-2 * 2.0 * accel * PR_DT,
        "box velocity after two calls: {}",
        baro.box_velocity
    );
}

#[test]
fn parrinello_rahman_box_keeps_moving_after_the_mismatch_vanishes() {
    let mut atoms = interior_pair();
    let mut pbox = PeriodicBox::cubic(BOX_EDGE);
    let mut baro = ParrinelloRahmanBarostat::new(1.0, 1.0, PR_MASS);
    baro.apply(&mut atoms, &mut pbox, 1.0, pr_virial(), PR_DT);
    let l_after_kick = pbox.dims.x;
    let v_after_kick = baro.box_velocity;
    assert!(v_after_kick > 0.0);
    // Feed exactly the target pressure: the force vanishes, but the box must
    // continue to coast with its stored velocity.
    let volume = pbox.volume();
    let neutral_virial = 3.0 * volume * 1.0;
    baro.apply(&mut atoms, &mut pbox, 1.0, neutral_virial, PR_DT);
    assert!((baro.box_velocity - v_after_kick).abs() < 1e-15);
    assert!(
        (pbox.dims.x - (l_after_kick + v_after_kick * PR_DT)).abs() < 1e-12,
        "coasting displacement"
    );
}

#[test]
fn parrinello_rahman_scales_positions_and_velocities_with_the_box() {
    let mut atoms = AtomSet::new();
    atoms.add_atom(
        Vec3::new(2.0, 4.0, 6.0),
        Vec3::new(1.0, -2.0, 0.5),
        2.0,
        0.0,
        0,
    );
    atoms.add_atom(
        Vec3::new(8.0, 1.0, 3.0),
        Vec3::new(-0.5, 0.25, 3.0),
        1.0,
        0.0,
        0,
    );
    let positions_before = atoms.positions.clone();
    let velocities_before = atoms.velocities.clone();
    let ke_before = atoms.kinetic_energy();
    let mut pbox = PeriodicBox::cubic(BOX_EDGE);
    let mut baro = ParrinelloRahmanBarostat::new(0.0, 1.0, 50.0);
    // Large virial: strong over-pressure, the box expands.
    baro.apply(&mut atoms, &mut pbox, 0.0, 3.0e4, 0.05);
    let scale = pbox.dims.x / BOX_EDGE;
    assert!(
        scale > 1.0,
        "over-pressure must expand the box, scale {scale}"
    );
    assert!((pbox.dims.y - pbox.dims.x).abs() < 1e-12 && (pbox.dims.z - pbox.dims.x).abs() < 1e-12);
    for (after, before) in atoms.positions.iter().zip(&positions_before) {
        assert!((after - before * scale).norm() < 1e-12);
    }
    for (after, before) in atoms.velocities.iter().zip(&velocities_before) {
        assert!((after - before / scale).norm() < 1e-12);
    }
    let ke_after = atoms.kinetic_energy();
    assert!(
        (ke_after - ke_before / (scale * scale)).abs() <= 1e-12 * ke_before,
        "expansion cools the system: {ke_before} -> {ke_after}"
    );
}

#[test]
fn parrinello_rahman_uses_the_target_pressure_argument() {
    let virial = pr_virial(); // P = 1.01
    // Own target equals the pressure, argument asks for lower pressure.
    let mut atoms = interior_pair();
    let mut pbox = PeriodicBox::cubic(BOX_EDGE);
    let mut baro = ParrinelloRahmanBarostat::new(1.01, 1.0, PR_MASS);
    baro.apply(&mut atoms, &mut pbox, 1.01, virial, PR_DT);
    assert!(baro.box_velocity.abs() < 1e-15, "no force at P = P0");
    baro.apply(&mut atoms, &mut pbox, 0.5, virial, PR_DT);
    assert!(baro.box_velocity > 0.0, "argument target below P expands");

    // NaN falls back to the barostat's own target pressure.
    let mut atoms = interior_pair();
    let mut pbox = PeriodicBox::cubic(BOX_EDGE);
    let mut baro = ParrinelloRahmanBarostat::new(0.5, 1.0, PR_MASS);
    baro.apply(&mut atoms, &mut pbox, f64::NAN, virial, PR_DT);
    assert!(baro.box_velocity > 0.0);
}

#[test]
fn parrinello_rahman_mass_is_derived_from_tau_p_when_not_positive() {
    let tau = 0.5_f64;
    let baro = ParrinelloRahmanBarostat::new(1.0, tau, 0.0);
    let ke = 10.0;
    let edge = 10.0;
    let expected = 3.0 * ke * tau * tau / (2.0 * std::f64::consts::PI.powi(2) * edge);
    let derived = baro
        .effective_mass(ke, edge)
        .expect("mass must be derivable");
    assert!((derived - expected).abs() < 1e-15 * expected);
    // A configured positive mass wins over tau_p.
    let fixed = ParrinelloRahmanBarostat::new(1.0, tau, 7.5);
    assert_eq!(fixed.effective_mass(ke, edge), Some(7.5));
    // No kinetic energy or no tau_p: nothing to derive the mass from.
    assert_eq!(baro.effective_mass(0.0, edge), None);
    assert_eq!(
        ParrinelloRahmanBarostat::new(1.0, 0.0, -1.0).effective_mass(ke, edge),
        None
    );
    let from_nan_mass = ParrinelloRahmanBarostat::new(1.0, tau, f64::NAN)
        .effective_mass(ke, edge)
        .expect("a NaN mass selects the derived mass");
    assert!((from_nan_mass - expected).abs() < 1e-12 * expected);
}

#[test]
fn parrinello_rahman_leaves_the_box_alone_without_a_usable_mass() {
    // Atoms at rest and mass <= 0: W cannot be derived from tau_p.
    let mut atoms = interior_pair();
    let mut pbox = PeriodicBox::cubic(BOX_EDGE);
    let mut baro = ParrinelloRahmanBarostat::new(1.0, 1.0, 0.0);
    baro.apply(&mut atoms, &mut pbox, 1.0, pr_virial(), PR_DT);
    assert_eq!(pbox.dims.x, BOX_EDGE);
    assert_eq!(baro.box_velocity, 0.0);
}

#[test]
fn parrinello_rahman_never_collapses_the_box() {
    let mut atoms = interior_pair();
    let mut pbox = PeriodicBox::cubic(BOX_EDGE);
    // A huge negative-pressure kick would drive the edge through zero.
    let mut baro = ParrinelloRahmanBarostat::new(1.0e9, 1.0, 1.0e-3);
    baro.apply(&mut atoms, &mut pbox, 1.0e9, 0.0, 1.0);
    assert_eq!(pbox.dims.x, BOX_EDGE, "collapsing step must be rejected");
    assert_eq!(baro.box_velocity, 0.0);
    assert!(atoms.positions[0].x.is_finite());
}

#[test]
fn parrinello_rahman_box_velocity_update_uses_three_v_over_w() {
    let baro = ParrinelloRahmanBarostat::new(1.0, 1.0, 100.0);
    let v = baro.box_velocity_update(0.0, 2.0, 1000.0, 0.001);
    // 3 * 1000 * (2 - 1) / 100 * 0.001
    assert!((v - 0.03).abs() < 1e-14, "{v}");
    let v2 = baro.box_velocity_update(v, 2.0, 1000.0, 0.001);
    assert!((v2 - 0.06).abs() < 1e-14, "{v2}");
}

#[test]
fn parrinello_rahman_box_oscillates_instead_of_creeping() {
    // Adiabatic ideal gas: P(L) = P(L0) (L0/L)^5 because the peculiar
    // velocities scale as 1/L.  With target pressure equal to the pressure at
    // L_eq = 12 the box (started at 10) overshoots L_eq and rings; the old
    // memoryless rescaler approached L_eq monotonically from below.
    let n_atoms = 10;
    let mut atoms = AtomSet::new();
    for i in 0..n_atoms {
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        atoms.add_atom(
            Vec3::new(1.0 + i as f64, 5.0, 5.0),
            Vec3::new(sign * 1.5, 0.0, 0.0),
            1.0,
            0.0,
            0,
        );
    }
    let ke0 = atoms.kinetic_energy();
    let mut pbox = PeriodicBox::cubic(BOX_EDGE);
    let p_initial = 2.0 * ke0 / (3.0 * pbox.volume());
    let l_eq = 12.0;
    let p_target = p_initial * (BOX_EDGE / l_eq).powi(5);
    let mut baro = ParrinelloRahmanBarostat::new(p_target, 1.0, 0.7);
    let mut max_edge = BOX_EDGE;
    let mut crossings = 0;
    let mut previous_side = false;
    for step in 0..4000 {
        baro.apply(&mut atoms, &mut pbox, p_target, 0.0, 0.01);
        let edge = pbox.dims.x;
        max_edge = max_edge.max(edge);
        let side = edge > l_eq;
        if step > 0 && side != previous_side {
            crossings += 1;
        }
        previous_side = side;
    }
    assert!(
        max_edge > 1.1 * l_eq,
        "box must overshoot, max edge {max_edge}"
    );
    assert!(
        crossings >= 4,
        "box must keep oscillating, crossings {crossings}"
    );
}

// ---------------------------------------------------------------------------
// Monte Carlo barostat acceptance (F3319)
// ---------------------------------------------------------------------------

#[test]
fn mc_acceptance_log_prob_matches_hand_value() {
    // Delta E = 0, V 1000 -> 1010, N = 100, P0 = 1, k_B T = 1:
    // -(1 * 10) + 101 ln(1.01) = -8.995016584
    let baro = MonteCarloBarostat::new(1.0);
    let log_p = baro.acceptance_log_prob(0.0, 1000.0, 1010.0, 100, 1.0, 1.0);
    assert!((log_p + 8.995_016_584).abs() < 1e-8, "log_p = {log_p}");
}

#[test]
fn mc_acceptance_jacobian_is_dimensionless_and_uses_n_plus_one() {
    // With P0 = 0 and Delta E = 0 only the Jacobian remains: (N + 1) ln(V'/V).
    let baro = MonteCarloBarostat::new(0.0);
    let log_p = baro.acceptance_log_prob(0.0, 500.0, 1000.0, 49, 300.0, 1.380649e-23);
    assert!(
        (log_p - 50.0 * 2.0_f64.ln()).abs() < 1e-12,
        "log_p = {log_p}"
    );
    // The old code returned N * Delta V, which would be 49 * 500 here.
    assert!(log_p < 40.0);
}

#[test]
fn mc_acceptance_rejects_non_positive_volumes() {
    let baro = MonteCarloBarostat::new(1.0);
    for (v_old, v_new) in [(0.0, 10.0), (10.0, 0.0), (-1.0, 10.0), (10.0, f64::NAN)] {
        assert_eq!(
            baro.acceptance_log_prob(0.0, v_old, v_new, 10, 1.0, 1.0),
            f64::NEG_INFINITY
        );
    }
}

#[test]
fn mc_volume_change_energy_is_consistent_with_the_acceptance_log_prob() {
    let baro = MonteCarloBarostat::new(0.7);
    let (kb, temp) = (1.0, 2.0);
    let beta = 1.0 / (kb * temp);
    for (de, v_old, v_new, n) in [(0.0, 1000.0, 1010.0, 100), (-3.0, 500.0, 480.0, 8)] {
        let du = baro.compute_volume_change_energy(de, v_old, v_new, n, temp, kb);
        let log_p = baro.acceptance_log_prob(de, v_old, v_new, n, temp, kb);
        assert!(
            (log_p + beta * du).abs() < 1e-12,
            "ln(acc) = {log_p} must equal -beta*dU = {}",
            -beta * du
        );
    }
}

#[test]
fn mc_barostat_samples_the_ideal_gas_volume_distribution() {
    // Ideal gas (Delta E = 0): p(V) ~ V^N exp(-beta P V), <V> = (N + 1) k_B T / P.
    // The N-instead-of-(N+1) Jacobian would give 20 / 0.05 = 400 here.
    let n_atoms = 20_usize;
    let (kb, temp, pressure) = (1.0, 1.0, 0.05);
    let mut baro = MonteCarloBarostat::new(pressure);
    baro.move_size = 0.4;
    let mut rng = Lcg(3);
    let mut volume = 400.0;
    let mut sum = 0.0;
    let mut count = 0_u32;
    for step in 0..300_000_u32 {
        let proposed = baro.propose_volume_change(volume, rng.next_unit());
        let log_acc = baro.acceptance_log_prob(0.0, volume, proposed, n_atoms, temp, kb);
        if rng.next_unit().ln() < log_acc {
            volume = proposed;
        }
        if step >= 5_000 {
            sum += volume;
            count += 1;
        }
    }
    let mean = sum / f64::from(count);
    let expected = (n_atoms as f64 + 1.0) * kb * temp / pressure;
    assert!(
        (mean - expected).abs() <= 0.02 * expected,
        "<V> = {mean}, expected {expected}"
    );
}

#[test]
fn mc_rng_barostat_samples_the_ideal_gas_volume_distribution() {
    use super::types::McBarostatRng;
    let n_atoms = 20_usize;
    let (kb, temp, pressure) = (1.0, 1.0, 0.05);
    let mut mc = McBarostatRng::with_seed(pressure, 0.4, 7);
    let mut volume = 400.0;
    let mut sum = 0.0;
    let mut count = 0_u32;
    for step in 0..300_000_u32 {
        if let Some(accepted) = mc.try_volume_move(volume, 0.0, n_atoms, temp, kb) {
            volume = accepted;
        }
        if step >= 5_000 {
            sum += volume;
            count += 1;
        }
    }
    let mean = sum / f64::from(count);
    let expected = (n_atoms as f64 + 1.0) * kb * temp / pressure;
    assert!(
        (mean - expected).abs() <= 0.02 * expected,
        "<V> = {mean}, expected {expected}"
    );
}
