//! Regression tests for the virial stress and the MD temperature estimate of the
//! continuum-atomistic coupling.

use super::*;
use oxiphysics_core::random::{Xoshiro256, box_muller};

/// Boltzmann constant in eV/K (CODATA 2018).
const K_B_EV_PER_K: f64 = 8.617_333_262e-5;

/// Lennard-Jones energy 4 eps ((sigma/r)^12 - (sigma/r)^6).
fn lj_energy(r: f64, epsilon: f64, sigma: f64) -> f64 {
    let sr6 = (sigma / r).powi(6);
    4.0 * epsilon * (sr6 * sr6 - sr6)
}

/// Lennard-Jones force magnitude -dU/dr (positive when repulsive).
fn lj_force_magnitude(r: f64, epsilon: f64, sigma: f64) -> f64 {
    let sr6 = (sigma / r).powi(6);
    24.0 * epsilon * (2.0 * sr6 * sr6 - sr6) / r
}

/// Minimum-image separation in a cubic box of edge `length`.
fn min_image(d: [f64; 3], length: f64) -> [f64; 3] {
    std::array::from_fn(|k| d[k] - length * (d[k] / length).round())
}

/// All `(i, j, f_ij)` Lennard-Jones pairs of `positions` in a cubic periodic box.
fn lj_pairs(
    positions: &[[f64; 3]],
    length: f64,
    epsilon: f64,
    sigma: f64,
) -> Vec<(usize, usize, [f64; 3])> {
    let mut pairs = Vec::new();
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            let r_vec = min_image(
                std::array::from_fn(|k| positions[i][k] - positions[j][k]),
                length,
            );
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            let magnitude = lj_force_magnitude(r, epsilon, sigma);
            pairs.push((i, j, std::array::from_fn(|k| magnitude * r_vec[k] / r)));
        }
    }
    pairs
}

/// Total Lennard-Jones energy with every position scaled by `scale` (box scaled too).
fn lj_energy_scaled(
    positions: &[[f64; 3]],
    length: f64,
    scale: f64,
    epsilon: f64,
    sigma: f64,
) -> f64 {
    let scaled_length = scale * length;
    let mut energy = 0.0;
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            let r_vec = min_image(
                std::array::from_fn(|k| scale * (positions[i][k] - positions[j][k])),
                scaled_length,
            );
            let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
            energy += lj_energy(r, epsilon, sigma);
        }
    }
    energy
}

fn maxwell_boltzmann_velocities(n: usize, mass: f64, temperature: f64, seed: u64) -> Vec<[f64; 3]> {
    let sigma_v = (K_B_AMU_A2_PS2 * temperature / mass).sqrt();
    let mut rng = Xoshiro256::new(seed);
    (0..n)
        .map(|_| {
            let (gx, gy) = box_muller(&mut rng);
            let (gz, _) = box_muller(&mut rng);
            [sigma_v * gx, sigma_v * gy, sigma_v * gz]
        })
        .collect()
}

fn fluctuating_hydrodynamics() -> FluctuatingHydrodynamics {
    FluctuatingHydrodynamics::new(0.0, 1e-3, 0.0, 0.6, 1e-27, 1e-15, 0)
}

// ── pair virial stress ────────────────────────────────────────────────────

#[test]
fn repulsive_pair_is_compressive_and_raises_the_pressure() {
    // Two atoms 2 Å apart in a 10 Å cube; atom 0 is pushed towards -x by 0.5 eV/Å.
    let positions = [[4.0, 5.0, 5.0], [6.0, 5.0, 5.0]];
    let pairs = [(0, 1, [-0.5, 0.0, 0.0])];
    let mut extractor = AtomisticToContin::new(10.0);
    let s = extractor
        .pair_virial_stress(&pairs, &positions, [10.0; 3])
        .expect("valid pair list");
    // r_01 = r_0 - r_1 = (-2, 0, 0), r_x f_x = +1 eV, V = 1000 Å^3.
    assert!((s[0] + 1.0e-3).abs() < 1e-15, "sigma_xx = {}", s[0]);
    assert!(s[1..].iter().all(|v| v.abs() < 1e-18), "{s:?}");
    assert_eq!(extractor.stress_voigt, s);
    let pressure = extractor.pressure();
    assert!((pressure - 1.0 / 3000.0).abs() < 1e-15, "P = {pressure}");

    // Reversing the force makes the pair attractive: tension and negative pressure.
    let attractive = [(0, 1, [0.5, 0.0, 0.0])];
    let s = extractor
        .pair_virial_stress(&attractive, &positions, [10.0; 3])
        .expect("valid pair list");
    assert!(s[0] > 0.0 && extractor.pressure() < 0.0);
}

#[test]
fn per_atom_form_equals_the_pair_form_without_periodicity() {
    // Non-periodic cluster: sum_i r_i (x) f_i == sum_{i<j} r_ij (x) f_ij for pair
    // forces that obey Newton's third law.
    let positions = [
        [1.0, 2.0, 3.0],
        [2.5, 2.2, 3.4],
        [1.7, 3.9, 2.6],
        [3.1, 3.3, 4.2],
    ];
    let strengths = [0.30, -0.12, 0.05, 0.21, -0.40, 0.09];
    let mut pairs = Vec::new();
    let mut forces = vec![[0.0_f64; 3]; positions.len()];
    let mut k = 0;
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            let f_ij: [f64; 3] =
                std::array::from_fn(|c| strengths[k] * (positions[i][c] - positions[j][c]));
            k += 1;
            for c in 0..3 {
                forces[i][c] += f_ij[c];
                forces[j][c] -= f_ij[c];
            }
            pairs.push((i, j, f_ij));
        }
    }
    let volume = 1.0e9;
    let mut per_atom = AtomisticToContin::new(10.0);
    let a = per_atom.potential_stress(&positions, &forces, volume);
    let mut periodic = AtomisticToContin::new(10.0);
    let b = periodic
        .pair_virial_stress(&pairs, &positions, [1.0e3; 3])
        .expect("valid pair list");
    for c in 0..6 {
        assert!(
            (a[c] - b[c]).abs() <= 1e-12 * a[c].abs().max(1e-12),
            "component {c}: per-atom {} vs pair {}",
            a[c],
            b[c]
        );
    }
}

#[test]
fn pair_split_across_the_boundary_gives_the_unwrapped_stress() {
    let length = 10.0;
    let f_ab = [0.4, 0.1, -0.2];
    let atom_a = [0.5, 5.0, 5.0];
    let atom_b_wrapped = [9.5, 5.0, 5.0];
    let atom_b_unwrapped = [-0.5, 5.0, 5.0];
    let pairs = [(0, 1, f_ab)];
    let mut extractor = AtomisticToContin::new(10.0);
    let wrapped = extractor
        .pair_virial_stress(&pairs, &[atom_a, atom_b_wrapped], [length; 3])
        .expect("valid pair list");
    let unwrapped = extractor
        .pair_virial_stress(&pairs, &[atom_a, atom_b_unwrapped], [length; 3])
        .expect("valid pair list");
    for c in 0..6 {
        assert!((wrapped[c] - unwrapped[c]).abs() < 1e-15, "component {c}");
    }
    // r_AB = (1, 0, 0) after the minimum-image wrap, V = 1000.
    assert!((wrapped[0] + 0.4e-3).abs() < 1e-15);
    assert!((wrapped[3] + 5.0e-5).abs() < 1e-15);
    assert!((wrapped[5] - 1.0e-4).abs() < 1e-15);
    assert!(wrapped[1].abs() < 1e-18 && wrapped[2].abs() < 1e-18 && wrapped[4].abs() < 1e-18);

    // The per-atom form uses the absolute positions and gets the wrapped pair wrong.
    let naive = extractor.potential_stress(
        &[atom_a, atom_b_wrapped],
        &[f_ab, [-0.4, -0.1, 0.2]],
        1000.0,
    );
    assert!((naive[0] - wrapped[0]).abs() > 1e-3, "naive = {}", naive[0]);
}

#[test]
fn pair_stress_is_invariant_under_translation_and_lattice_shifts() {
    let length = 12.0;
    let positions = [
        [1.0, 2.0, 3.0],
        [4.2, 2.9, 3.4],
        [7.9, 8.1, 2.6],
        [3.1, 10.8, 9.2],
        [11.3, 5.5, 6.1],
    ];
    let pairs = lj_pairs(&positions, length, 0.01, 2.5);
    let mut extractor = AtomisticToContin::new(10.0);
    let reference = extractor
        .pair_virial_stress(&pairs, &positions, [length; 3])
        .expect("valid pair list");
    assert!(reference.iter().any(|v| v.abs() > 1e-6));

    let shift = [3.7, -2.2, 8.1];
    let lattice = [[1, 0, -2], [0, 3, 0], [-1, -1, 1], [2, 0, 0], [0, 0, -3]];
    let moved: Vec<[f64; 3]> = positions
        .iter()
        .zip(lattice.iter())
        .map(|(p, n)| std::array::from_fn(|c| p[c] + shift[c] + f64::from(n[c]) * length))
        .collect();
    let shifted = extractor
        .pair_virial_stress(&pairs, &moved, [length; 3])
        .expect("valid pair list");
    for c in 0..6 {
        assert!(
            (shifted[c] - reference[c]).abs() <= 1e-12 * reference[c].abs().max(1e-9),
            "component {c}: {} vs {}",
            shifted[c],
            reference[c]
        );
    }
}

#[test]
fn virial_pressure_matches_the_negative_energy_volume_derivative() {
    // P_virial = -tr(sigma_pot) / 3 must equal -dU/dV under a uniform dilation.
    let (length, epsilon, sigma) = (12.0, 0.01, 2.5);
    let positions = [
        [1.0, 2.0, 3.0],
        [3.6, 2.4, 3.3],
        [2.1, 5.2, 5.9],
        [5.4, 4.6, 2.2],
    ];
    let pairs = lj_pairs(&positions, length, epsilon, sigma);
    let mut extractor = AtomisticToContin::new(10.0);
    let s = extractor
        .pair_virial_stress(&pairs, &positions, [length; 3])
        .expect("valid pair list");
    let pressure = extractor.pressure();
    assert!(s.iter().all(|v| v.is_finite()));

    let delta = 1.0e-5;
    let u_plus = lj_energy_scaled(&positions, length, 1.0 + delta, epsilon, sigma);
    let u_minus = lj_energy_scaled(&positions, length, 1.0 - delta, epsilon, sigma);
    let v_plus = ((1.0 + delta) * length).powi(3);
    let v_minus = ((1.0 - delta) * length).powi(3);
    let numeric = -(u_plus - u_minus) / (v_plus - v_minus);
    assert!(
        (pressure - numeric).abs() <= 1e-6 * numeric.abs(),
        "virial pressure {pressure} vs -dU/dV = {numeric}"
    );
}

#[test]
fn ideal_gas_pressure_is_n_k_t_over_v() {
    // A third of the atoms moves along each axis with |v|² = 3 k_B T / m, so Σ m v_x² =
    // Σ m v_y² = Σ m v_z² = N k_B T and the shear components vanish.
    let n = 1200;
    let (mass, temperature) = (39.948, 300.0);
    let v_th = (K_B_AMU_A2_PS2 * temperature / mass).sqrt();
    let speed = 3.0_f64.sqrt() * v_th;
    let velocities: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let mut v = [0.0; 3];
            v[i % 3] = speed;
            v
        })
        .collect();
    let masses = vec![mass; n];
    let positions = vec![[1.0, 2.0, 3.0]; n];
    let box_lengths = [20.0, 20.0, 20.0];
    let volume = 8000.0;
    let mut extractor = AtomisticToContin::new(10.0);
    let s = extractor
        .total_pair_stress(&[], &positions, box_lengths, &velocities, &masses)
        .expect("valid box");
    let expected = n as f64 * K_B_EV_PER_K * temperature / volume;
    let pressure = extractor.pressure();
    assert!(
        ((pressure - expected) / expected).abs() < 1e-9,
        "P = {pressure} eV/Å^3, N k_B T / V = {expected}"
    );
    assert!(s[3..].iter().all(|v| v.abs() < 1e-18), "{s:?}");

    // With interactions the pressure is the ideal-gas term plus the virial term
    // (1 / 3V) sum r_ij . f_ij.
    let two = [[4.0, 5.0, 5.0], [6.0, 5.0, 5.0]];
    let pairs = [(0, 1, [-0.5, 0.0, 0.0])];
    let mut extractor = AtomisticToContin::new(10.0);
    extractor
        .total_pair_stress(&pairs, &two, [10.0; 3], &[[v_th; 3]; 2], &[mass; 2])
        .expect("valid box");
    let ideal = 2.0 * K_B_EV_PER_K * temperature / 1000.0;
    assert!(
        (extractor.pressure() - (ideal + 1.0 / 3000.0)).abs() < 1e-9 * ideal,
        "P = {}",
        extractor.pressure()
    );
}

#[test]
fn pair_virial_rejects_invalid_boxes_and_pairs() {
    let positions = [[0.0; 3], [1.0, 0.0, 0.0]];
    let pairs = [(0, 1, [0.1, 0.0, 0.0])];
    let mut extractor = AtomisticToContin::new(10.0);
    for bad_box in [
        [0.0, 5.0, 5.0],
        [5.0, -5.0, 5.0],
        [5.0, 5.0, f64::NAN],
        [f64::INFINITY, 5.0, 5.0],
    ] {
        assert!(
            extractor
                .pair_virial_stress(&pairs, &positions, bad_box)
                .is_err(),
            "box {bad_box:?}"
        );
        assert!(
            extractor
                .total_pair_stress(&pairs, &positions, bad_box, &[], &[])
                .is_err()
        );
    }
    assert!(
        extractor
            .pair_virial_stress(&[(0, 2, [0.0; 3])], &positions, [5.0; 3])
            .is_err()
    );
    assert!(
        extractor
            .pair_virial_stress(&[(1, 1, [0.0; 3])], &positions, [5.0; 3])
            .is_err()
    );
}

// ── temperature from the MD kinetic energy ───────────────────────────────

#[test]
fn temperature_of_a_maxwell_boltzmann_argon_gas_is_recovered() {
    // 1e4 argon atoms sampled at 300 K: sigma_v = sqrt(k_B T / m) in Å/ps.  With the old
    // k_B = 8.3e-3 the estimate came out 100 times too high.
    let n = 10_000;
    let mass = 39.948;
    let velocities = maxwell_boltzmann_velocities(n, mass, 300.0, 20_241_019);
    let masses = vec![mass; n];
    let mut fh = fluctuating_hydrodynamics();
    fh.update_temperature_from_md(&velocities, &masses);
    assert!(
        (fh.temperature - 300.0).abs() < 6.0,
        "T = {} K, expected 300 K within 2%",
        fh.temperature
    );
}

#[test]
fn temperature_is_exact_for_a_hand_computed_velocity() {
    // One atom of 39.948 amu with v = (1, 2, 2) Å/ps: KE = 0.5 * 39.948 * 9 amu Å^2/ps^2,
    // T = 2 KE / (3 k_B) with k_B = 0.8314462618 amu Å^2/ps^2/K.
    let mut fh = fluctuating_hydrodynamics();
    fh.update_temperature_from_md(&[[1.0, 2.0, 2.0]], &[39.948]);
    let expected = 39.948 * 9.0 / (3.0 * 0.831_446_261_815_324);
    assert!(
        (fh.temperature - expected).abs() < 1e-9,
        "T = {}",
        fh.temperature
    );
    // 39.948 * 9 / (3 * 0.8314) = 144.1 K
    assert!((fh.temperature - 144.14).abs() < 0.01);
}

#[test]
fn com_momentum_removal_uses_three_n_minus_three_degrees_of_freedom() {
    let n = 500;
    let mass = 39.948;
    let mut velocities = maxwell_boltzmann_velocities(n, mass, 300.0, 7);
    let mean: [f64; 3] =
        std::array::from_fn(|c| velocities.iter().map(|v| v[c]).sum::<f64>() / n as f64);
    for v in velocities.iter_mut() {
        for c in 0..3 {
            v[c] -= mean[c];
        }
    }
    let masses = vec![mass; n];

    let mut all_dof = fluctuating_hydrodynamics();
    all_dof.update_temperature_from_md(&velocities, &masses);
    let mut without_com = fluctuating_hydrodynamics().with_com_momentum_removed(true);
    assert!(without_com.com_momentum_removed);
    without_com.update_temperature_from_md(&velocities, &masses);

    let ratio = without_com.temperature / all_dof.temperature;
    let expected_ratio = (3 * n) as f64 / (3 * n - 3) as f64;
    assert!(
        (ratio - expected_ratio).abs() < 1e-12,
        "ratio {ratio} vs {expected_ratio}"
    );
}

#[test]
fn temperature_update_handles_degenerate_inputs() {
    // Empty system: the temperature is left untouched.
    let mut fh = FluctuatingHydrodynamics::new(123.0, 1e-3, 0.0, 0.6, 1e-27, 1e-15, 0);
    fh.update_temperature_from_md(&[], &[]);
    assert_eq!(fh.temperature, 123.0);

    // A lone atom has no centre-of-mass-free degrees of freedom: fall back to 3N.
    let mut single = fluctuating_hydrodynamics().with_com_momentum_removed(true);
    single.update_temperature_from_md(&[[1.0, 0.0, 0.0]], &[2.0]);
    let expected = 2.0 * 0.5 * 2.0 / (3.0 * K_B_AMU_A2_PS2);
    assert!((single.temperature - expected).abs() < 1e-12);

    // Mismatched lengths use the atoms that have both a velocity and a mass.
    let mut partial = fluctuating_hydrodynamics();
    partial.update_temperature_from_md(
        &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [5.0, 5.0, 5.0]],
        &[2.0, 2.0],
    );
    let expected = 2.0 * (0.5 * 2.0 + 0.5 * 2.0) / (6.0 * K_B_AMU_A2_PS2);
    assert!((partial.temperature - expected).abs() < 1e-12);
}
