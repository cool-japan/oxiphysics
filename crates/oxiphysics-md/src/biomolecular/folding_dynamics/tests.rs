//! Tests of the Gō-model Langevin dynamics of [`ProteinFolding`].

use super::*;
use crate::biomolecular::types::ProteinTopology;

/// Cα trace of an ideal α-helix (radius 2.3 Å, rise 1.5 Å, 100° per residue): the
/// consecutive beads are 3.83 Å apart and every bead has native contacts with the
/// residues `i ± 2 … 4` inside an 8 Å cutoff.
fn helix_topology(n_residues: usize) -> ProteinTopology {
    let mut topology = ProteinTopology::new();
    for i in 0..n_residues {
        let angle = i as f64 * 100.0_f64.to_radians();
        topology
            .positions
            .push([2.3 * angle.cos(), 2.3 * angle.sin(), 1.5 * i as f64]);
        topology.residues.push("ALA".to_string());
        topology.atom_names.push("CA".to_string());
        topology.atom_residue.push(i);
        if i > 0 {
            topology.bonds.push([i - 1, i]);
        }
    }
    topology.n_atoms = n_residues;
    topology
}

fn helix_folding(n_residues: usize, temperature: f64, dt: f64) -> ProteinFolding {
    ProteinFolding::new(helix_topology(n_residues), temperature, 8.0, dt)
}

/// Native contacts that are not bonds (the 10-12 contacts of the model).
fn non_bonded_native_contacts(fold: &ProteinFolding) -> usize {
    fold.native_contacts
        .iter()
        .filter(|&&[i, j]| !fold.topology.bonds.contains(&[i.min(j), i.max(j)]))
        .count()
}

/// Mean of the last `window` entries of the Q history.
fn tail_mean(history: &[f64], window: usize) -> f64 {
    let tail = &history[history.len() - window..];
    tail.iter().sum::<f64>() / window as f64
}

#[test]
fn test_native_structure_is_a_force_free_minimum() {
    let fold = helix_folding(14, 0.0, 0.005);
    // |i - j| = 2, 3, 4 contacts of a 14-residue helix: 12 + 11 + 10.
    assert_eq!(fold.native_contacts.len(), 33);
    assert_eq!(non_bonded_native_contacts(&fold), 33);
    let forces = fold.forces();
    let largest = forces
        .iter()
        .flat_map(|f| f.iter())
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    // Only the (4 Å / 8.7 Å)^12 repulsion of the non-native pairs is left.
    assert!(largest < 5e-3, "native forces up to {largest} kcal/mol/Å");
    let energy = fold.potential_energy();
    assert!(
        (energy + 33.0).abs() < 0.05,
        "E = {energy}, expected -33 ε from the 33 contacts"
    );
    // The native structure is the minimum: any displacement raises the energy.
    let mut displaced = fold.clone();
    displaced.topology.positions[7][0] += 0.3;
    assert!(displaced.potential_energy() > energy + 0.1);
}

#[test]
fn test_forces_are_the_negative_gradient_of_the_potential() {
    let mut fold = helix_folding(10, 300.0, 0.005);
    // A distorted conformation: stretched bonds, broken contacts and close non-native pairs.
    for (index, position) in fold.topology.positions.iter_mut().enumerate() {
        let t = index as f64;
        position[0] += 0.9 * (1.7 * t).sin();
        position[1] += 0.8 * (2.3 * t + 0.4).cos();
        position[2] += 0.7 * (0.9 * t + 1.1).sin();
    }
    let forces = fold.forces();
    let scale = forces
        .iter()
        .flat_map(|f| f.iter())
        .fold(1.0_f64, |m, v| m.max(v.abs()));
    assert!(scale > 1.0, "the distorted chain must be under strain");
    let step = 1.0e-6;
    for (atom, force) in forces.iter().enumerate() {
        for (axis, &analytic) in force.iter().enumerate() {
            let mut plus = fold.clone();
            plus.topology.positions[atom][axis] += step;
            let mut minus = fold.clone();
            minus.topology.positions[atom][axis] -= step;
            let numeric = -(plus.potential_energy() - minus.potential_energy()) / (2.0 * step);
            assert!(
                (analytic - numeric).abs() <= 1e-6 * scale,
                "atom {atom} axis {axis}: F = {analytic}, -dV/dx = {numeric}"
            );
        }
    }
    // Internal forces: no net force.
    for axis in 0..3 {
        let net: f64 = forces.iter().map(|f| f[axis]).sum();
        assert!(net.abs() < 1e-9 * scale, "net force {net}");
    }
}

#[test]
fn test_bond_oscillation_period_pins_the_unit_conversion() {
    // Two beads joined by one bond: omega = sqrt(k_b · 418.4 / μ) with k_b = 200 ε/Å²
    // = 200 kcal/mol/Å², μ = 55 amu; the period is 2π/ω = 0.1611 ps.
    let mut topology = ProteinTopology::new();
    topology.positions = vec![[0.0, 0.0, 0.0], [3.8, 0.0, 0.0]];
    topology.bonds.push([0, 1]);
    topology.n_atoms = 2;
    let mut fold = ProteinFolding::new(topology, 0.0, 8.0, 0.0005);
    fold.gamma = 0.0;
    fold.topology.positions[1][0] = 3.9;
    let omega =
        (fold.bond_force_constant() * KCAL_MOL_TO_AMU_A2_PS2 / (0.5 * DEFAULT_BEAD_MASS)).sqrt();
    let period = 2.0 * std::f64::consts::PI / omega;
    assert!((period - 0.1611).abs() < 1e-3, "analytic period {period}");

    let stretch =
        |f: &ProteinFolding| f.topology.positions[1][0] - f.topology.positions[0][0] - 3.8;
    let initial_energy = fold.potential_energy() + fold.kinetic_energy();
    let mut crossings: Vec<f64> = Vec::new();
    let mut previous = stretch(&fold);
    for step in 1..=1000 {
        fold.step();
        let current = stretch(&fold);
        if previous.signum() != current.signum() {
            // Linear interpolation of the zero crossing inside the step.
            let fraction = previous / (previous - current);
            crossings.push((step as f64 - 1.0 + fraction) * fold.dt);
        }
        previous = current;
        let total = fold.potential_energy() + fold.kinetic_energy();
        assert!(
            (total - initial_energy).abs() < 2e-3 * initial_energy,
            "step {step}: energy drift {total} vs {initial_energy}"
        );
    }
    assert!(crossings.len() >= 4, "{crossings:?}");
    // Consecutive zero crossings are half a period apart.
    let measured = crossings[3] - crossings[1];
    assert!(
        (measured - period).abs() < 0.01 * period,
        "measured period {measured} ps vs {period} ps"
    );
    // The centre of mass does not drift for this internal force.
    let com = 0.5 * (fold.topology.positions[0][0] + fold.topology.positions[1][0]);
    assert!((com - 1.95).abs() < 1e-9, "centre of mass moved to {com}");
}

#[test]
fn test_native_state_is_stable_at_zero_temperature() {
    let mut fold = helix_folding(14, 0.0, 0.005);
    fold.run(1000);
    assert_eq!(fold.step_count, 1000);
    assert_eq!(fold.q_history.len(), 1000);
    assert!(fold.rmsd() < 0.02, "rmsd {}", fold.rmsd());
    assert!(fold.q_history.iter().all(|&q| (q - 1.0).abs() < 1e-12));
}

#[test]
fn test_helix_stays_folded_at_50_k_and_unfolds_at_1000_k() {
    let mut cold = helix_folding(14, 50.0, 0.005).with_seed(1);
    cold.thermalize_velocities();
    cold.run(5000);
    let q_cold = tail_mean(&cold.q_history, 1000);
    assert!(q_cold > 0.9, "Q at 50 K = {q_cold}");

    let mut hot = helix_folding(14, 1000.0, 0.005).with_seed(1);
    hot.thermalize_velocities();
    hot.run(5000);
    let q_hot = tail_mean(&hot.q_history, 1000);
    assert!(q_hot < 0.5, "Q at 1000 K = {q_hot}");
    assert!(hot.rmsd() > 3.0, "rmsd at 1000 K = {}", hot.rmsd());
    // Q really moves: the history is not a constant.
    let (min, max) = hot
        .q_history
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &q| {
            (lo.min(q), hi.max(q))
        });
    assert!(max - min > 0.3, "Q history spans {min}..{max}");
}

#[test]
fn test_kinetic_temperature_follows_the_langevin_thermostat() {
    let mut fold = helix_folding(14, 300.0, 0.005).with_seed(7);
    fold.thermalize_velocities();
    fold.run(2000);
    let mut sum = 0.0;
    let samples = 6000;
    for _ in 0..samples {
        fold.step();
        sum += fold.kinetic_temperature();
    }
    let mean = sum / samples as f64;
    assert!(
        (mean - 300.0).abs() < 45.0,
        "mean kinetic temperature {mean} K at a 300 K bath"
    );
}

#[test]
fn test_langevin_noise_is_reproducible_from_the_seed() {
    let run = |seed: u64| {
        let mut fold = helix_folding(8, 400.0, 0.005).with_seed(seed);
        fold.run(200);
        fold
    };
    let first = run(5);
    let second = run(5);
    let other = run(6);
    assert_eq!(first.topology.positions, second.topology.positions);
    assert_eq!(first.velocities, second.velocities);
    assert_eq!(first.q_history, second.q_history);
    assert_ne!(first.topology.positions, other.topology.positions);
    assert_ne!(first.rng_state, [0; 4]);
}

#[test]
fn test_quench_lowers_the_thermostat_temperature() {
    let mut fold = helix_folding(14, 800.0, 0.005).with_seed(3);
    fold.thermalize_velocities();
    fold.run(2000);
    let hot = fold.kinetic_temperature();
    for _ in 0..40 {
        fold.quench(50.0, 0.2);
    }
    assert!(
        (fold.temperature - 50.0).abs() < 1.0,
        "{}",
        fold.temperature
    );
    fold.run(3000);
    let cold = fold.kinetic_temperature();
    assert!(cold < 0.3 * hot, "kinetic temperature {hot} K -> {cold} K");
}

/// Three beads on a line with the given spacing (one native contact, beads 0 and 2).
fn three_bead_line(spacing: f64) -> ProteinFolding {
    let mut topology = ProteinTopology::new();
    topology.positions = vec![[0.0; 3], [spacing, 0.0, 0.0], [2.0 * spacing, 0.0, 0.0]];
    topology.bonds = vec![[0, 1], [1, 2]];
    topology.n_atoms = 3;
    ProteinFolding::new(topology, 300.0, 8.0, 0.002)
}

#[test]
fn test_contact_fraction_uses_the_native_distance_of_every_contact() {
    // One contact with native separation r0 is formed up to 1.2 r0.  The old test against a
    // fixed 1.2 * 8 Å = 9.6 Å counted a 3 Å contact as formed out to 9.6 Å.
    for (spacing, inside, outside) in [(1.5, 3.59, 3.61), (3.8, 9.1, 9.14)] {
        let mut fold = three_bead_line(spacing);
        assert_eq!(fold.native_contacts, vec![[0, 2]]);
        assert_eq!(fold.contact_fraction_q(), 1.0);
        fold.topology.positions[2][0] = inside;
        assert_eq!(
            fold.contact_fraction_q(),
            1.0,
            "d = {inside} Å is inside 1.2 r0"
        );
        fold.topology.positions[2][0] = outside;
        assert_eq!(
            fold.contact_fraction_q(),
            0.0,
            "d = {outside} Å is beyond 1.2 r0"
        );
    }

    // Helix of six: 4 + 3 + 2 = 9 contacts; sending the last bead away breaks the three
    // contacts it takes part in and no other.
    let mut fold = helix_folding(6, 300.0, 0.005);
    assert_eq!(fold.native_contacts.len(), 9);
    assert_eq!(fold.contact_fraction_q(), 1.0);
    fold.topology.positions[5][0] += 50.0;
    assert!((fold.contact_fraction_q() - 6.0 / 9.0).abs() < 1e-12);
    // A fully separated chain has lost every contact.
    for (i, position) in fold.topology.positions.iter_mut().enumerate() {
        *position = [30.0 * i as f64, 0.0, 0.0];
    }
    assert_eq!(fold.contact_fraction_q(), 0.0);
    // Without native contacts Q is 1 by convention.
    let mut lone = ProteinTopology::new();
    lone.positions = vec![[0.0; 3]];
    lone.n_atoms = 1;
    assert_eq!(
        ProteinFolding::new(lone, 300.0, 8.0, 0.002).contact_fraction_q(),
        1.0
    );
}

#[test]
fn test_step_with_invalid_parameters_only_records_q() {
    for (dt, mass) in [
        (0.0, 110.0),
        (-0.005, 110.0),
        (f64::NAN, 110.0),
        (0.005, 0.0),
    ] {
        let mut fold = helix_folding(6, 500.0, dt);
        fold.bead_mass = mass;
        let before = fold.topology.positions.clone();
        fold.run(3);
        assert_eq!(fold.topology.positions, before, "dt = {dt}, mass = {mass}");
        assert_eq!(fold.step_count, 3);
        assert_eq!(fold.q_history, vec![1.0; 3]);
    }
    // A velocity list of the wrong length is reset instead of panicking.
    let mut fold = helix_folding(6, 0.0, 0.005);
    fold.velocities.clear();
    fold.step();
    assert_eq!(fold.velocities.len(), 6);
    // An all-zero generator state is replaced rather than producing silent zeros.
    let mut stuck = helix_folding(6, 600.0, 0.005);
    stuck.rng_state = [0; 4];
    stuck.run(10);
    assert_ne!(stuck.rng_state, [0; 4]);
    assert!(stuck.velocities.iter().flatten().any(|v| *v != 0.0));
}

#[test]
fn test_bonded_and_adjacent_pairs_carry_no_non_native_repulsion() {
    // The residue backbone N, CA, C, O: the O ... N(next) pair is two bonds apart and
    // only 2.2 Å away; the (4 Å / r)^12 repulsion would be ~1000 ε if it were not excluded.
    let mut topology = ProteinTopology::new();
    topology.add_residue(
        "ALA",
        [0.0, 0.0, 0.0],
        [1.46, 0.0, 0.0],
        [2.0, 1.4, 0.0],
        [1.2, 2.4, 0.0],
    );
    topology.add_residue(
        "GLY",
        [3.3, 1.5, 0.2],
        [4.2, 2.6, 0.6],
        [5.6, 2.3, 0.9],
        [6.0, 1.2, 1.1],
    );
    let fold = ProteinFolding::new(topology, 0.0, 10.0, 0.002);
    let classes = fold.pair_classes();
    // O of residue 0 (atom 3) and N of residue 1 (atom 4): consecutive atoms, 1-3 via C.
    assert!(classes.no_repulsion.contains(3, 4));
    // The 1-3 pair N(1) ... C(1) shares CA; the 1-3 pair CA(0) ... N(1) shares C(0).
    assert!(classes.no_repulsion.contains(4, 6));
    assert!(classes.no_repulsion.contains(1, 4));
    assert!(classes.bonded.contains(2, 4), "C(0)-N(1) peptide bond");
    assert!(!classes.bonded.contains(3, 4));
    let native_energy = fold.potential_energy();
    // Every non-bonded pair two or more atoms apart is a native contact: 21 pairs, one of
    // which (the peptide bond) is a bond and not counted again.
    assert!((native_energy + 20.0).abs() < 1e-9, "E = {native_energy}");
    let largest = fold
        .forces()
        .iter()
        .flat_map(|f| f.iter())
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    assert!(largest < 1e-9, "native backbone forces up to {largest}");
}
