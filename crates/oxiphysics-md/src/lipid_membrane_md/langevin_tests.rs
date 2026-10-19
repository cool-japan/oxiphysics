// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the Langevin dynamics, the periodic bonded forces and
//! the molecule wrapping of [`BilayerSimulation`].
//!
//! The previous integrator drew its "noise" from a hash of the step counter
//! (the same number for every lipid, uniformly distributed, with variance 1/12
//! instead of 1) and evaluated bonds on raw coordinate differences although the
//! beads were wrapped one by one.  The tests below check the statistics of the
//! new BAOAB integrator against exact equipartition results.

use super::bonded::{BOND_FORCE_CONSTANT, BOND_R0_NM};
use super::*;
use rand::RngExt;
use rand::rngs::StdRng;

const TEMPERATURE: f64 = 310.0;

fn thermal_energy() -> f64 {
    KB_KJ_PER_MOL_K * TEMPERATURE
}

/// A lipid made of unbonded single beads.
fn free_bead(position: [f64; 3]) -> CgLipid {
    CgLipid {
        species: LipidSpecies::Dppc,
        beads: vec![CgBead::new(BeadType::ApolarTail, position, 72.0)],
        bonds: Vec::new(),
        angles: Vec::new(),
        leaflet: 0,
        flipping: false,
    }
}

/// A chain lipid whose bonds and angles are removed: `n_beads` free beads.
fn unbonded_chain(n_beads: usize) -> CgLipid {
    let mut lipid = CgLipid::new(LipidSpecies::Dppc, [1.0, 1.0, 1.0], 0);
    lipid.beads.truncate(n_beads);
    lipid.bonds.clear();
    lipid.angles.clear();
    lipid
}

/// Two beads joined by one harmonic bond at its rest length.
fn dimer(position: [f64; 3], direction: [f64; 3]) -> CgLipid {
    let norm = (direction[0].powi(2) + direction[1].powi(2) + direction[2].powi(2)).sqrt();
    let second = [
        position[0] + BOND_R0_NM * direction[0] / norm,
        position[1] + BOND_R0_NM * direction[1] / norm,
        position[2] + BOND_R0_NM * direction[2] / norm,
    ];
    CgLipid {
        species: LipidSpecies::Dppc,
        beads: vec![
            CgBead::new(BeadType::ApolarTail, position, 72.0),
            CgBead::new(BeadType::ApolarTail, second, 72.0),
        ],
        bonds: vec![(0, 1)],
        angles: Vec::new(),
        leaflet: 0,
        flipping: false,
    }
}

fn random_direction(rng: &mut StdRng) -> [f64; 3] {
    loop {
        let v = [
            rng.random::<f64>() - 0.5,
            rng.random::<f64>() - 0.5,
            rng.random::<f64>() - 0.5,
        ];
        let norm2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        if norm2 > 1.0e-4 && norm2 < 0.25 {
            return v;
        }
    }
}

fn assert_close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{what}: {actual} vs {expected} (tolerance {tolerance})"
    );
}

#[test]
fn single_step_noise_has_the_fluctuation_dissipation_variance() {
    // From rest and without forces one step leaves v = sqrt((1 - c1^2) kT / m) xi
    // with c1 = exp(-gamma dt), i.e. m <v_x^2> = (1 - c1^2) kT.
    let n = 1000;
    let lipids = (0..n).map(|_| free_bead([1.0, 1.0, 1.0])).collect();
    let mut sim = BilayerSimulation::from_lipids(lipids, SimBox::cubic(5.0), TEMPERATURE, 5);
    sim.step_langevin();
    let c1 = (-sim.friction * sim.dt).exp();
    let expected = (1.0 - c1 * c1) * thermal_energy();
    let mean_m_v2 = sim
        .lipids
        .iter()
        .flat_map(|lipid| lipid.beads.iter())
        .flat_map(|bead| bead.velocity.iter().map(|v| bead.mass * v * v))
        .sum::<f64>()
        / (3 * n) as f64;
    assert_close(mean_m_v2 / expected, 1.0, 0.1, "single-step m<v^2> ratio");
}

#[test]
fn free_beads_equilibrate_to_equipartition_with_gaussian_velocities() {
    let n = 1000;
    let lipids = (0..n).map(|_| free_bead([1.0, 1.0, 1.0])).collect();
    let mut sim = BilayerSimulation::from_lipids(lipids, SimBox::cubic(5.0), TEMPERATURE, 17);
    // Ten velocity relaxation times 1/gamma = 1 ps = 50 steps.
    for _ in 0..500 {
        sim.step_langevin();
    }
    let production_steps = 20_000;
    let mut kinetic_sum = 0.0;
    let mut velocity_samples: Vec<f64> = Vec::new();
    for step in 0..production_steps {
        sim.step_langevin();
        kinetic_sum += sim.kinetic_energy();
        if step % 100 == 0 {
            velocity_samples.extend(
                sim.lipids
                    .iter()
                    .flat_map(|lipid| lipid.beads.iter().map(|bead| bead.velocity[0])),
            );
        }
    }
    let mean_kinetic_per_bead = kinetic_sum / production_steps as f64 / n as f64;
    assert_close(
        mean_kinetic_per_bead / (1.5 * thermal_energy()),
        1.0,
        0.03,
        "<1/2 m v^2> / (3/2 kT)",
    );
    let samples = velocity_samples.len() as f64;
    let mean = velocity_samples.iter().sum::<f64>() / samples;
    let variance = velocity_samples
        .iter()
        .map(|v| (v - mean).powi(2))
        .sum::<f64>()
        / samples;
    let fourth = velocity_samples
        .iter()
        .map(|v| (v - mean).powi(4))
        .sum::<f64>()
        / samples;
    // Maxwell-Boltzmann: Var(v_x) = kT / m and kurtosis 3.  Uniform noise would
    // give a kurtosis of 1.8.
    assert_close(
        variance * 72.0 / thermal_energy(),
        1.0,
        0.03,
        "m Var(v_x) / kT",
    );
    assert_close(fourth / (variance * variance), 3.0, 0.15, "kurtosis of v_x");
}

#[test]
fn noise_is_independent_between_lipids_and_beads() {
    // Two identical unbonded 12-bead lipids at the same state: with the old
    // hash noise bead k of both lipids received exactly the same kick.
    let lipids = vec![unbonded_chain(12), unbonded_chain(12)];
    let mut sim = BilayerSimulation::from_lipids(lipids, SimBox::cubic(5.0), TEMPERATURE, 23);
    for _ in 0..200 {
        sim.step_langevin();
    }
    let (mut sum_ab, mut sum_aa, mut sum_bb) = (0.0, 0.0, 0.0);
    let mut sum_within = 0.0;
    let mut sum_within_norm = 0.0;
    for _ in 0..10_000 {
        sim.step_langevin();
        for bead in 0..12 {
            for dim in 0..3 {
                let a = sim.lipids[0].beads[bead].velocity[dim];
                let b = sim.lipids[1].beads[bead].velocity[dim];
                sum_ab += a * b;
                sum_aa += a * a;
                sum_bb += b * b;
            }
            // Neighbouring beads of one lipid must not share their noise either.
            if bead + 1 < 12 {
                let here = sim.lipids[0].beads[bead].velocity[0];
                let next = sim.lipids[0].beads[bead + 1].velocity[0];
                sum_within += here * next;
                sum_within_norm += here * here;
            }
        }
    }
    let between_lipids = sum_ab / (sum_aa * sum_bb).sqrt();
    let within_lipid = sum_within / sum_within_norm;
    assert!(
        between_lipids.abs() < 0.1,
        "velocity correlation between lipids: {between_lipids}"
    );
    assert!(
        within_lipid.abs() < 0.1,
        "velocity correlation between neighbouring beads: {within_lipid}"
    );
}

#[test]
fn seeded_runs_are_reproducible_and_seed_dependent() {
    let make = |seed: u64| {
        BilayerSimulation::new_flat_bilayer_with_seed(LipidSpecies::Dppc, 4, TEMPERATURE, seed)
    };
    let positions = |sim: &BilayerSimulation| -> Vec<[f64; 3]> {
        sim.lipids
            .iter()
            .flat_map(|lipid| lipid.beads.iter().map(|bead| bead.position))
            .collect()
    };
    let advance = |sim: &mut BilayerSimulation| {
        for _ in 0..40 {
            sim.step_langevin();
        }
    };
    let (mut a, mut b, mut c) = (make(7), make(7), make(8));
    advance(&mut a);
    advance(&mut b);
    advance(&mut c);
    assert_eq!(positions(&a), positions(&b), "same seed, same trajectory");
    assert_ne!(
        positions(&a),
        positions(&c),
        "different seed, different noise"
    );
    // `new_flat_bilayer` is the default seed.
    let mut default_run = BilayerSimulation::new_flat_bilayer(LipidSpecies::Dppc, 4, TEMPERATURE);
    let mut explicit = make_with_default_seed(LipidSpecies::Dppc, 4, TEMPERATURE);
    advance(&mut default_run);
    advance(&mut explicit);
    assert_eq!(positions(&default_run), positions(&explicit));
    // Reseeding restarts the noise stream.
    let mut fresh = make(7);
    let mut restarted = make(99);
    restarted.reseed(7);
    advance(&mut fresh);
    advance(&mut restarted);
    assert_eq!(positions(&fresh), positions(&restarted));
}

fn make_with_default_seed(
    species: LipidSpecies,
    n_per_leaflet: usize,
    temperature: f64,
) -> BilayerSimulation {
    BilayerSimulation::new_flat_bilayer_with_seed(
        species,
        n_per_leaflet,
        temperature,
        BilayerSimulation::DEFAULT_SEED,
    )
}

#[test]
fn bonded_dimers_sample_the_canonical_bond_and_velocity_distribution() {
    let mut rng = StdRng::seed_from_u64(3);
    let box_ = SimBox::cubic(30.0);
    let lipids: Vec<CgLipid> = (0..60)
        .map(|_| {
            let position = [
                10.0 + 10.0 * rng.random::<f64>(),
                10.0 + 10.0 * rng.random::<f64>(),
                10.0 + 10.0 * rng.random::<f64>(),
            ];
            dimer(position, random_direction(&mut rng))
        })
        .collect();
    let mut sim = BilayerSimulation::from_lipids(lipids, box_, TEMPERATURE, 41);
    for _ in 0..1000 {
        sim.step_langevin();
    }
    let production_steps = 8000;
    let (mut bond_energy_sum, mut temperature_sum) = (0.0, 0.0);
    for _ in 0..production_steps {
        sim.step_langevin();
        bond_energy_sum += sim.total_intramolecular_energy();
        temperature_sum += sim.kinetic_temperature();
    }
    // One harmonic degree of freedom per dimer: <1/2 k (r - r0)^2> = kT / 2.
    let mean_bond_energy = bond_energy_sum / production_steps as f64 / 60.0;
    assert_close(
        mean_bond_energy / (0.5 * thermal_energy()),
        1.0,
        0.08,
        "bond energy per dimer / (kT/2)",
    );
    let mean_temperature = temperature_sum / production_steps as f64;
    assert_close(
        mean_temperature / TEMPERATURE,
        1.0,
        0.03,
        "kinetic temperature",
    );
}

#[test]
fn bonded_forces_are_the_negative_gradient_of_the_box_energy() {
    let box_ = SimBox::new(3.0, 3.5, 6.0);
    let mut rng = StdRng::seed_from_u64(77);
    // Chain with angles, randomly distorted so that bonds and angles are away
    // from their minima; the second copy is split across the x boundary.
    let mut lipid = CgLipid::new(LipidSpecies::Dppc, [2.95, 1.0, 3.0], 0);
    for bead in &mut lipid.beads {
        for coordinate in &mut bead.position {
            *coordinate += 0.2 * (rng.random::<f64>() - 0.5);
        }
    }
    let mut split = lipid.clone();
    for bead in &mut split.beads {
        bead.position[0] = bead.position[0].rem_euclid(box_.lx);
        bead.position[1] = bead.position[1].rem_euclid(box_.ly);
    }
    assert!(
        split.beads.iter().any(|b| b.position[0] < 1.0)
            && split.beads.iter().any(|b| b.position[0] > 2.0),
        "the chain straddles the periodic boundary"
    );
    for molecule in [&lipid, &split] {
        let forces = molecule.intramolecular_forces(&box_);
        let step = 1.0e-6;
        for (bead, bead_force) in forces.iter().enumerate() {
            for (dim, &force) in bead_force.iter().enumerate() {
                let mut plus = molecule.clone();
                let mut minus = molecule.clone();
                plus.beads[bead].position[dim] += step;
                minus.beads[bead].position[dim] -= step;
                let gradient = (plus.intramolecular_energy_in_box(&box_)
                    - minus.intramolecular_energy_in_box(&box_))
                    / (2.0 * step);
                assert!(
                    (force + gradient).abs() <= 1.0e-5 * force.abs().max(1.0),
                    "bead {bead} axis {dim}: force {force} vs -dE/dx {}",
                    -gradient
                );
            }
        }
        // No external field: the molecule feels no net force.
        for dim in 0..3 {
            let net: f64 = forces.iter().map(|f| f[dim]).sum();
            assert!(net.abs() < 1.0e-9, "net force {net} along axis {dim}");
        }
    }
    // Splitting the molecule does not change its energy.
    assert_close(
        split.intramolecular_energy_in_box(&box_),
        lipid.intramolecular_energy_in_box(&box_),
        1.0e-9,
        "energy of the split molecule",
    );
}

#[test]
fn bond_across_the_boundary_has_minimum_image_energy_and_force() {
    let box_ = SimBox::new(5.0, 5.0, 5.0);
    // Bead 1 sits 0.47 nm "behind" bead 0 through the x boundary.
    let mut molecule = dimer([4.9, 2.0, 2.0], [1.0, 0.0, 0.0]);
    molecule.beads[1].position[0] = (4.9 + BOND_R0_NM).rem_euclid(box_.lx);
    assert!(molecule.beads[1].position[0] < 0.5);
    assert!(
        molecule.bond_energy() > 1.0e4,
        "raw coordinates see a bond of length 4.4 nm: {}",
        molecule.bond_energy()
    );
    assert!(molecule.bond_energy_in_box(&box_) < 1.0e-12);
    for force in molecule.intramolecular_forces(&box_) {
        assert!(force.iter().all(|f| f.abs() < 1.0e-8), "{force:?}");
    }
    // Stretch it by 0.03 nm: the restoring force is k * 0.03 along the bond.
    molecule.beads[1].position[0] += 0.03;
    let forces = molecule.intramolecular_forces(&box_);
    assert_close(
        forces[0][0],
        BOND_FORCE_CONSTANT * 0.03,
        1.0e-6,
        "pull on bead 0",
    );
    assert_close(
        forces[1][0],
        -BOND_FORCE_CONSTANT * 0.03,
        1.0e-6,
        "pull on bead 1",
    );
}

#[test]
fn wrap_into_box_makes_split_molecules_whole() {
    let box_ = SimBox::new(5.0, 4.0, 8.0);
    // A 12-bead chain lying along x that crosses the boundary, then wrapped bead
    // by bead the way the old integrator did.
    let mut chain = CgLipid::new(LipidSpecies::Dppc, [4.7, 1.0, 4.0], 0);
    for (index, bead) in chain.beads.iter_mut().enumerate() {
        bead.position = [4.7 + BOND_R0_NM * index as f64, 1.0, 4.0];
    }
    let whole_energy = chain.bond_energy();
    assert!(whole_energy < 1.0e-9);
    for bead in &mut chain.beads {
        bead.position[0] = bead.position[0].rem_euclid(box_.lx);
    }
    assert!(
        chain.bond_energy() > 1.0e3,
        "split molecule: raw bonds are huge"
    );
    chain.wrap_into_box(&box_);
    assert!(chain.bond_energy() < 1.0e-9, "{}", chain.bond_energy());
    assert_close(
        chain.beads[0].position[0],
        4.7,
        1.0e-12,
        "head stays inside the box",
    );
    for (index, bead) in chain.beads.iter().enumerate() {
        assert_close(
            bead.position[0],
            4.7 + BOND_R0_NM * index as f64,
            1.0e-9,
            "bead position along the unwrapped chain",
        );
    }
    // A molecule whose head left the box is translated back as a whole.
    let mut outside = chain.clone();
    for bead in &mut outside.beads {
        bead.position[0] += 5.0;
        bead.position[1] -= 4.0;
    }
    outside.wrap_into_box(&box_);
    assert_close(outside.beads[0].position[0], 4.7, 1.0e-12, "head x");
    assert_close(outside.beads[0].position[1], 1.0, 1.0e-12, "head y");
    assert!(outside.bond_energy() < 1.0e-9);
    // Single free beads are simply wrapped.
    let mut bead = free_bead([-0.25, 9.0, 1.0]);
    bead.wrap_into_box(&box_);
    assert_close(bead.beads[0].position[0], 4.75, 1.0e-12, "free bead x");
    assert_close(bead.beads[0].position[1], 1.0, 1.0e-12, "free bead y");
}

#[test]
fn langevin_dynamics_keeps_molecules_whole_across_the_boundary() {
    let mut rng = StdRng::seed_from_u64(9);
    let box_ = SimBox::new(4.0, 4.0, 4.0);
    // Dimers straddling the x = 0 / x = 4 boundary, plus some in the bulk.
    let lipids: Vec<CgLipid> = (0..40)
        .map(|i| {
            let x = if i % 2 == 0 {
                -0.2 + 0.4 * rng.random::<f64>()
            } else {
                1.0 + 2.0 * rng.random::<f64>()
            };
            let position = [x.rem_euclid(box_.lx), 2.0 * rng.random::<f64>() + 1.0, 2.0];
            dimer(position, random_direction(&mut rng))
        })
        .collect();
    let mut sim = BilayerSimulation::from_lipids(lipids, box_, TEMPERATURE, 31);
    sim.wrap_molecules();
    let mut crossed = false;
    for _ in 0..3000 {
        sim.step_langevin();
        for lipid in &sim.lipids {
            crossed |= lipid
                .beads
                .iter()
                .any(|b| b.position[0] < 0.0 || b.position[0] > 4.0);
            let raw = lipid.bond_energy();
            let minimum_image = lipid.bond_energy_in_box(&sim.box_);
            assert!(
                (raw - minimum_image).abs() <= 1.0e-9 * (1.0 + raw),
                "molecule split: raw {raw} vs minimum image {minimum_image}"
            );
            let head = lipid.beads[0].position;
            assert!((0.0..sim.box_.lx).contains(&head[0]) && (0.0..sim.box_.ly).contains(&head[1]));
        }
    }
    assert!(crossed, "some molecules must have crossed the boundary");
    // Whole molecules have the thermal bond energy, not a box-sized artefact.
    let mean_energy = sim.total_intramolecular_energy() / 40.0;
    assert!(
        mean_energy < 5.0 * thermal_energy(),
        "mean bond energy {mean_energy}"
    );
}

#[test]
fn zero_friction_conserves_the_total_energy() {
    // gamma = 0 reduces BAOAB to velocity Verlet.
    let mut rng = StdRng::seed_from_u64(13);
    let box_ = SimBox::new(6.0, 6.0, 8.0);
    let lipids: Vec<CgLipid> = (0..6)
        .map(|i| {
            let mut lipid = CgLipid::new(
                LipidSpecies::Dppc,
                [1.0 + i as f64, 1.0 + 0.5 * i as f64, 2.0],
                0,
            );
            for bead in &mut lipid.beads {
                for coordinate in &mut bead.position {
                    *coordinate += 0.05 * (rng.random::<f64>() - 0.5);
                }
                for velocity in &mut bead.velocity {
                    *velocity = 0.3 * (rng.random::<f64>() - 0.5);
                }
            }
            lipid
        })
        .collect();
    let mut sim = BilayerSimulation::from_lipids(lipids, box_, TEMPERATURE, 1);
    sim.friction = 0.0;
    sim.dt = 0.005;
    let initial = sim.kinetic_energy() + sim.total_intramolecular_energy();
    let mut worst = 0.0_f64;
    for _ in 0..4000 {
        sim.step_langevin();
        let total = sim.kinetic_energy() + sim.total_intramolecular_energy();
        worst = worst.max((total - initial).abs() / initial);
    }
    assert!(initial > 1.0, "initial energy {initial}");
    assert!(worst < 0.01, "energy drift {worst}");
}

#[test]
fn kinetic_temperature_of_an_empty_system_is_zero() {
    let sim = BilayerSimulation::from_lipids(Vec::new(), SimBox::cubic(2.0), TEMPERATURE, 1);
    assert_eq!(sim.kinetic_energy(), 0.0);
    assert_eq!(sim.kinetic_temperature(), 0.0);
}

#[test]
fn beads_without_a_positive_mass_stay_fixed() {
    let mut lipids = vec![free_bead([1.0, 1.0, 1.0]), free_bead([2.0, 1.0, 1.0])];
    lipids[0].beads[0].mass = 0.0;
    lipids[1].beads[0].mass = f64::NAN;
    let mut sim = BilayerSimulation::from_lipids(lipids, SimBox::cubic(4.0), TEMPERATURE, 2);
    for _ in 0..20 {
        sim.step_langevin();
    }
    for lipid in &sim.lipids {
        assert_eq!(lipid.beads[0].velocity, [0.0; 3]);
    }
    assert_eq!(sim.lipids[0].beads[0].position, [1.0, 1.0, 1.0]);
    assert_eq!(sim.lipids[1].beads[0].position, [2.0, 1.0, 1.0]);
}
