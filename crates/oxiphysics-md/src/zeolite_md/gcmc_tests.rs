// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the Grand Canonical Monte Carlo sampler: genuine
//! Metropolis acceptance (no deterministic thresholds), the thermal wavelength
//! in the activity, the framework in the energy, panic-free indexing and
//! reproducible random numbers.

use super::*;

/// Cubic periodic framework of edge `edge` (Å) with the given T-sites
/// (fractional coordinates).
fn cubic_framework(edge: f64, t_sites: Vec<[f64; 3]>) -> ZeoliteFramework {
    ZeoliteFramework {
        framework_type: FrameworkType::Custom,
        cell_lengths: [edge; 3],
        cell_angles: [90.0; 3],
        pore_limiting_diameter: edge,
        largest_included_sphere: edge,
        t_sites,
        o_sites: Vec::new(),
        si_al_ratio: f64::INFINITY,
    }
}

/// Mean and population variance of the loading series.
fn mean_and_variance(series: &[f64]) -> (f64, f64) {
    let n = series.len() as f64;
    let mean = series.iter().sum::<f64>() / n;
    let var = series.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    (mean, var)
}

/// Chemical potential (J) at which non-interacting guests have `⟨N⟩ = n_mean`.
fn ideal_gas_mu(sim: &GcmcSimulation, n_mean: f64) -> f64 {
    K_B * sim.temperature * (n_mean * sim.thermal_wavelength().powi(3) / sim.volume_m3()).ln()
}

#[test]
fn thermal_wavelength_matches_the_textbook_value() {
    let ff = MoleculeZeoliteFF::new(0.0, 3.5, 0.0, 0.0, 12.0);
    let sim = GcmcSimulation::new(300.0, 0.0, ff, cubic_framework(10.0, Vec::new()));
    // Methane at 300 K: Lambda = h / sqrt(2 pi m k_B T) = 0.2517 Angstrom.
    let lambda_angstrom = sim.thermal_wavelength() / ANGSTROM;
    assert!(
        (lambda_angstrom - 0.2517).abs() < 1.0e-3,
        "Lambda = {lambda_angstrom} Å"
    );
    // Lambda ~ T^(-1/2) and ~ m^(-1/2).
    let mut hot = sim.clone();
    hot.temperature = 1200.0;
    assert!((hot.thermal_wavelength() / sim.thermal_wavelength() - 0.5).abs() < 1.0e-12);
    let heavy = sim.clone().with_mass(4.0 * sim.mass);
    assert!((heavy.thermal_wavelength() / sim.thermal_wavelength() - 0.5).abs() < 1.0e-12);
    assert!((sim.volume_m3() - 1.0e-27).abs() < 1.0e-39);
}

#[test]
fn ideal_gas_loading_follows_the_poisson_mean_of_the_grand_canonical_ensemble() {
    // No framework, no interactions: P(N) is Poisson with mean V exp(beta mu) / Lambda^3.
    let ff = MoleculeZeoliteFF::new(0.0, 3.5, 0.0, 0.0, 12.0);
    let mut sim =
        GcmcSimulation::new(300.0, 0.0, ff, cubic_framework(10.0, Vec::new())).with_seed(2024);
    let target = 20.0;
    sim.chemical_potential = ideal_gas_mu(&sim, target);
    sim.run(20_000, 1.0);
    let mut series = Vec::with_capacity(400_000);
    for _ in 0..400_000 {
        sim.mc_step(1.0);
        series.push(sim.state.n_molecules() as f64);
    }
    let (mean, var) = mean_and_variance(&series);
    assert!(
        (mean / target - 1.0).abs() < 0.05,
        "<N> = {mean}, grand-canonical ideal-gas value {target}"
    );
    assert!(
        (var / mean - 1.0).abs() < 0.2,
        "Var(N)/<N> = {} (Poisson statistics expected)",
        var / mean
    );
    assert!(sim.n_accept_ins > 0 && sim.n_accept_del > 0 && sim.n_accept_trans > 0);
    assert_eq!(
        sim.state.energy, 0.0,
        "non-interacting guests carry no energy"
    );
}

#[test]
fn the_framework_changes_the_acceptance_probability() {
    // One T-site in the middle of a 20 Å box.
    let site = vec![[0.5, 0.5, 0.5]];
    let epsilon = 2.5e-19; // 60 k_B T at 300 K
    let r_min = 2.0_f64.powf(1.0 / 6.0) * 3.5;
    let ff = MoleculeZeoliteFF::new(epsilon, 3.5, 0.0, 0.0, 9.0);
    let mu = -60.0 * K_B * 300.0;
    let make = |t_sites: Vec<[f64; 3]>| {
        GcmcSimulation::new(300.0, mu, ff.clone(), cubic_framework(20.0, t_sites)).with_seed(5)
    };

    // 1 Å from the T-site the Lennard-Jones core is ~3e6 k_B T high: the
    // insertion is impossible, however favourable the chemical potential.
    let mut sim = make(site.clone());
    let mut open = make(Vec::new());
    open.chemical_potential = 0.0;
    sim.chemical_potential = 0.0;
    let overlap = [10.0, 10.0, 11.0];
    let accepted_with_framework = (0..500).filter(|_| sim.try_insert(overlap)).count();
    assert_eq!(accepted_with_framework, 0, "overlapping insertion accepted");
    assert_eq!(sim.state.n_molecules(), 0);
    assert_eq!(sim.state.energy, 0.0);
    assert_eq!(sim.n_attempts, 500);
    assert!(
        open.try_insert(overlap),
        "without a framework the site is free"
    );

    // At the Lennard-Jones minimum the gain of 60 k_B T beats a chemical
    // potential that makes the bare insertion hopeless (probability e^-47).
    let attractive = [10.0, 10.0, 10.0 + r_min];
    let mut bare = make(Vec::new());
    let accepted_bare = (0..200).filter(|_| bare.try_insert(attractive)).count();
    assert_eq!(accepted_bare, 0, "bare insertion should be hopeless");
    let mut sim = make(site);
    let e_site = sim.framework_energy(attractive);
    assert!(
        (e_site + epsilon).abs() < 1.0e-9 * epsilon,
        "U at the LJ minimum = {e_site}, expected -{epsilon}"
    );
    assert!(
        sim.try_insert(attractive),
        "attractive site must be accepted"
    );
    assert_eq!(sim.state.n_molecules(), 1);
    assert!((sim.state.energy - e_site).abs() < 1.0e-12 * epsilon);
    assert!((sim.compute_energy() - e_site).abs() < 1.0e-12 * epsilon);
}

#[test]
fn guest_guest_energy_is_lennard_jones_and_host_guest_energy_includes_coulomb() {
    let ff = MoleculeZeoliteFF::new(3.0e-21, 3.4, 0.3, -0.5, 9.0);
    let fw = cubic_framework(20.0, vec![[0.5, 0.5, 0.5]]);
    let mut sim = GcmcSimulation::new(300.0, 0.0, ff.clone(), fw);
    // Two guests 4 Å apart along x, 6 Å from the T-site at (10, 10, 10) ...
    let (a, b) = ([4.0, 10.0, 10.0], [8.0, 10.0, 10.0]);
    sim.state.positions.push(a);
    sim.state.positions.push(b);
    let expected = ff.lj_energy(4.0) + ff.total_energy(6.0) + ff.total_energy(2.0);
    let total = sim.compute_energy();
    assert!(
        (total - expected).abs() < 1.0e-12 * expected.abs(),
        "E = {total}, expected {expected}"
    );
    // The guest charge couples to the framework charge, not to other guests.
    let mut free = GcmcSimulation::new(300.0, 0.0, ff.clone(), cubic_framework(20.0, Vec::new()));
    free.state.positions.push(a);
    free.state.positions.push(b);
    assert_eq!(free.compute_energy(), ff.lj_energy(4.0));
    assert!(ff.coulomb_energy(4.0).abs() > 0.0);
}

#[test]
fn out_of_range_indices_are_not_trials_and_do_not_panic() {
    let ff = MoleculeZeoliteFF::new(0.0, 3.5, 0.0, 0.0, 9.0);
    let mut sim = GcmcSimulation::new(300.0, 0.0, ff, cubic_framework(12.0, Vec::new()));
    assert!(!sim.try_delete(0), "empty cell");
    for pos in [[1.0, 1.0, 1.0], [5.0, 5.0, 5.0], [9.0, 9.0, 9.0]] {
        sim.state.positions.push(pos);
    }
    assert!(!sim.try_delete(3));
    assert!(!sim.try_delete(usize::MAX));
    assert!(!sim.try_translate(7, 1.0));
    assert!(!sim.try_move_to(3, [1.0, 2.0, 3.0]));
    assert_eq!(sim.n_attempts, 0);
    assert_eq!(sim.state.n_molecules(), 3);
}

#[test]
fn translation_is_downhill_certain_and_uphill_forbidden_at_overlap() {
    let ff = MoleculeZeoliteFF::new(1.0e-21, 3.5, 0.0, 0.0, 9.0);
    let fw = cubic_framework(20.0, vec![[0.5, 0.5, 0.5]]);
    let mut sim = GcmcSimulation::new(300.0, 0.0, ff, fw);
    // A guest 2 Å from the T-site sits on the repulsive wall (~770 k_B T).
    let repulsive = [10.0, 10.0, 12.0];
    let attractive = [10.0, 10.0, 10.0 + 2.0_f64.powf(1.0 / 6.0) * 3.5];
    sim.state.positions.push(repulsive);
    sim.state.energy = sim.compute_energy();
    assert!(sim.state.energy > 100.0 * K_B * 300.0);

    assert!(
        sim.try_move_to(0, attractive),
        "downhill move must be accepted"
    );
    assert_eq!(sim.state.positions[0], attractive);
    assert_eq!(sim.n_accept_trans, 1);
    assert!(sim.state.energy < 0.0);
    assert!((sim.state.energy - sim.compute_energy()).abs() < 1.0e-9 * 1.0e-21);

    for _ in 0..50 {
        assert!(!sim.try_move_to(0, repulsive), "uphill move onto the wall");
    }
    assert_eq!(sim.state.positions[0], attractive);
    assert_eq!(sim.n_accept_trans, 1);
}

#[test]
fn metropolis_accept_follows_the_boltzmann_factor() {
    let ff = MoleculeZeoliteFF::new(0.0, 3.5, 0.0, 0.0, 9.0);
    let mut sim =
        GcmcSimulation::new(300.0, 0.0, ff, cubic_framework(10.0, Vec::new())).with_seed(99);
    let kt = K_B * 300.0;
    assert!(sim.accept(0.0));
    assert!(sim.accept(-5.0 * kt));
    assert!(!sim.accept(1.0e-15), "enormous uphill energy");
    assert!(!sim.accept(f64::NAN));
    let trials = 40_000;
    let accepted = (0..trials)
        .filter(|_| sim.accept(kt * 2.0_f64.ln()))
        .count();
    let fraction = accepted as f64 / f64::from(trials);
    assert!(
        (fraction - 0.5).abs() < 0.02,
        "exp(-ln 2) = 0.5, observed acceptance {fraction}"
    );
}

#[test]
fn runs_are_reproducible_for_a_fixed_seed() {
    let ff = MoleculeZeoliteFF::new(2.0e-21, 3.4, 0.0, 0.0, 6.0);
    let build = |seed: u64| {
        let mut sim =
            GcmcSimulation::new(300.0, 0.0, ff.clone(), cubic_framework(14.0, Vec::new()))
                .with_seed(seed);
        sim.chemical_potential = ideal_gas_mu(&sim, 5.0);
        sim
    };
    let (mut a, mut b, mut c) = (build(7), build(7), build(8));
    for sim in [&mut a, &mut b, &mut c] {
        sim.run(3_000, 1.0);
    }
    assert!(
        a.state.n_molecules() > 0,
        "the run should have adsorbed guests"
    );
    assert_eq!(a.state.positions, b.state.positions);
    assert_eq!(a.state.energy, b.state.energy);
    assert_eq!(a.n_accept_ins, b.n_accept_ins);
    assert_ne!(
        a.state.positions, c.state.positions,
        "different seed, same run"
    );
}

#[test]
fn incremental_energy_matches_a_full_recomputation() {
    // Two T-sites with charges, interacting guests, all three move types.
    let ff = MoleculeZeoliteFF::new(4.0e-21, 3.4, 0.2, -0.4, 7.0);
    let fw = cubic_framework(15.0, vec![[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]);
    let mut sim = GcmcSimulation::new(300.0, 0.0, ff.clone(), fw).with_seed(31);
    sim.chemical_potential = ideal_gas_mu(&sim, 6.0);
    for _ in 0..20 {
        sim.run(1_000, 1.5);
        let direct = sim.compute_energy();
        assert!(
            (sim.state.energy - direct).abs() < 1.0e-9 * direct.abs().max(ff.epsilon),
            "bookkeeping drifted: incremental {} vs direct {direct}",
            sim.state.energy
        );
    }
    assert!(sim.n_accept_ins > 0, "no insertion accepted");
    assert!(sim.n_accept_del > 0, "no deletion accepted");
    assert!(sim.n_accept_trans > 0, "no translation accepted");
    assert!(sim.state.n_molecules() > 0);
    assert!(sim.insertion_rate() > 0.0 && sim.insertion_rate() < 1.0);
}

/// Boltzmann-factor average `W = <exp(-beta U)>` over the box for one
/// Lennard-Jones site at its centre, by radial Simpson quadrature
/// `W = 1 + (4 pi / V) ∫ (exp(-beta U(r)) - 1) r² dr`.
fn single_site_boltzmann_average(sim: &GcmcSimulation) -> f64 {
    let beta = 1.0 / (K_B * sim.temperature);
    let volume = sim.state.pbox.volume();
    let (r_start, r_end) = (1.0_f64, sim.ff.cutoff);
    let intervals = 20_000_usize;
    let h = (r_end - r_start) / intervals as f64;
    let integrand = |r: f64| ((-beta * sim.ff.total_energy(r)).exp() - 1.0) * r * r;
    let mut sum = integrand(r_start) + integrand(r_end);
    for k in 1..intervals {
        let weight = if k % 2 == 1 { 4.0 } else { 2.0 };
        sum += weight * integrand(r_start + k as f64 * h);
    }
    // Below 1 Å the Lennard-Jones core excludes the guest: exp(-beta U) = 0.
    let core = -1.0 / 3.0;
    1.0 + 4.0 * PI / volume * (sum * h / 3.0 + core)
}

#[test]
fn dilute_loading_is_the_ideal_gas_value_times_the_framework_boltzmann_factor() {
    // A single attractive site (5 k_B T deep) in a 30 Å box: at low density the
    // guests are independent and <N> = (V e^{beta mu} / Lambda^3) * <e^{-beta U}>.
    let epsilon = 5.0 * K_B * 300.0;
    let ff = MoleculeZeoliteFF::new(epsilon, 3.5, 0.0, 0.0, 12.0);
    let mut sim = GcmcSimulation::new(300.0, 0.0, ff, cubic_framework(30.0, vec![[0.5, 0.5, 0.5]]))
        .with_seed(12_345);
    let ideal_mean = 0.01;
    sim.chemical_potential = ideal_gas_mu(&sim, ideal_mean);
    let w = single_site_boltzmann_average(&sim);
    assert!(w > 1.3, "the test needs a clearly attractive site, W = {w}");

    sim.run(100_000, 2.0);
    let moves = 2_000_000;
    let mut total = 0.0;
    for _ in 0..moves {
        sim.mc_step(2.0);
        total += sim.state.n_molecules() as f64;
    }
    let mean = total / f64::from(moves);
    let expected = ideal_mean * w;
    assert!(
        (mean / expected - 1.0).abs() < 0.10,
        "<N> = {mean}, expected {expected} (ideal value {ideal_mean} x W = {w}); \
         a sampler that ignores the framework would give {ideal_mean}"
    );
}
