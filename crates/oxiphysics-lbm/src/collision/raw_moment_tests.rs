//! Regression tests for the raw-moment collision operator.
//!
//! The operator used to relax the three normal-stress moments `m_200, m_020, m_002` with the
//! bulk rate (so the trace-free normal stresses, i.e. extensional shear, relaxed at the wrong
//! rate) and to rebuild the populations from `m_000 … m_011` only, discarding every ghost
//! moment. It now relaxes the trace with `omega_b`, the trace-free stresses with `omega_v`, the
//! ghosts with `omega_g`, and inverts the full 19×19 raw-moment matrix exactly.

use super::test_support::{
    assert_conserved, perturbed_populations, shear_wave_amplitude_ratio, shear_wave_decay,
};
use super::*;
use crate::lattice::D3Q19_VELOCITIES;

/// The 19×19 raw-moment matrix: row `k` holds `c_x^a c_y^b c_z^c` of the exponents of moment `k`
/// at every velocity.
fn moment_matrix() -> [[f64; 19]; 19] {
    std::array::from_fn(|k| {
        let e = MOMENT_EXPONENTS[k];
        std::array::from_fn(|i| {
            let c = D3Q19_VELOCITIES[i];
            f64::from(c[0]).powi(e[0] as i32)
                * f64::from(c[1]).powi(e[1] as i32)
                * f64::from(c[2]).powi(e[2] as i32)
        })
    })
}

/// Solve `a x = b` by Gaussian elimination with partial pivoting.
fn solve(mut a: [[f64; 19]; 19], mut b: [f64; 19]) -> [f64; 19] {
    for col in 0..19 {
        let pivot = (col..19)
            .max_by(|&r1, &r2| a[r1][col].abs().total_cmp(&a[r2][col].abs()))
            .unwrap_or(col);
        a.swap(col, pivot);
        b.swap(col, pivot);
        let diagonal = a[col][col];
        for row in col + 1..19 {
            let factor = a[row][col] / diagonal;
            let pivot_row = a[col];
            for (entry, &p) in a[row].iter_mut().zip(pivot_row.iter()).skip(col) {
                *entry -= factor * p;
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = [0.0_f64; 19];
    for row in (0..19).rev() {
        let mut acc = b[row];
        for (k, &xk) in x.iter().enumerate().skip(row + 1) {
            acc -= a[row][k] * xk;
        }
        x[row] = acc / a[row][row];
    }
    x
}

/// Deterministic pseudo-random vector with entries in `[-1, 1)`.
fn random_vector(seed: u64) -> [f64; 19] {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    std::array::from_fn(|_| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        2.0 * ((state >> 11) as f64 / (1u64 << 53) as f64) - 1.0
    })
}

fn max_abs_diff(a: &[f64; 19], b: &[f64; 19]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

#[test]
fn moment_matrix_reproduces_the_raw_moments() {
    let matrix = moment_matrix();
    for seed in 1..=4 {
        let f = perturbed_populations(seed, 0.3);
        let m = compute_raw_moments_d3q19(&f);
        for (k, row) in matrix.iter().enumerate() {
            let direct: f64 = row.iter().zip(f.iter()).map(|(a, b)| a * b).sum();
            assert!((direct - m[k]).abs() < 1e-14, "moment {k}");
        }
    }
}

#[test]
fn inverse_agrees_with_gaussian_elimination_of_the_moment_matrix() {
    let matrix = moment_matrix();
    for seed in 1..=8 {
        let m = random_vector(seed);
        let reference = solve(matrix, m);
        let closed_form = raw_moments_to_f_d3q19(&m);
        assert!(
            max_abs_diff(&reference, &closed_form) < 1e-12,
            "seed {seed}: {reference:?} vs {closed_form:?}"
        );
    }
}

#[test]
fn raw_moment_transform_round_trips_to_round_off() {
    for seed in 1..=10 {
        let f = perturbed_populations(seed, 0.3);
        let back = raw_moments_to_f_d3q19(&compute_raw_moments_d3q19(&f));
        assert!(max_abs_diff(&f, &back) < 1e-13, "seed {seed}");
    }
    for k in 0..19 {
        // Unit populations (every column of M) and unit moments (every column of M^-1).
        let mut f = [0.0_f64; 19];
        f[k] = 1.0;
        let back = raw_moments_to_f_d3q19(&compute_raw_moments_d3q19(&f));
        assert!(max_abs_diff(&f, &back) < 1e-13, "unit population {k}");
        let mut m = [0.0_f64; 19];
        m[k] = 1.0;
        let back = compute_raw_moments_d3q19(&raw_moments_to_f_d3q19(&m));
        assert!(max_abs_diff(&m, &back) < 1e-13, "unit moment {k}");
    }
}

#[test]
fn collision_conserves_mass_and_momentum_for_every_rate_choice() {
    for omega_v in [0.6, 1.25, 1.9] {
        for omega_b in [0.5, 1.0, 1.9] {
            for omega_g in [0.4, 1.0, 1.7] {
                let coll = RawMomentCollision {
                    omega_v,
                    omega_b,
                    omega_g,
                };
                for seed in 1..=3 {
                    let f = perturbed_populations(seed, 0.2);
                    let (rho, u) = compute_macroscopic(&f);
                    assert_conserved(
                        &f,
                        &coll.collide(&f, rho, u),
                        1e-14,
                        &format!("rates ({omega_v}, {omega_b}, {omega_g}), seed {seed}"),
                    );
                }
            }
        }
    }
}

#[test]
fn equilibrium_is_a_fixed_point() {
    let (rho, u) = (0.8, [0.05, -0.02, 0.03]);
    let feq = compute_equilibrium_d3q19(rho, u);
    let coll = RawMomentCollision::from_viscosity(0.05, 0.02);
    let out = coll.collide(&feq, rho, u);
    assert!(max_abs_diff(&out, &feq) < 1e-14);
}

#[test]
fn equal_rates_reproduce_bgk() {
    // A single rate on every non-conserved moment is BGK; this fails whenever a group of
    // moments is dropped or mapped to a different rate.
    for omega in [0.7, 1.3, 1.95] {
        let coll = RawMomentCollision {
            omega_v: omega,
            omega_b: omega,
            omega_g: omega,
        };
        for seed in 1..=5 {
            let f = perturbed_populations(seed, 0.2);
            let (rho, u) = compute_macroscopic(&f);
            let bgk = collide_bgk(&f, rho, u, omega);
            let out = coll.collide(&f, rho, u);
            assert!(
                max_abs_diff(&out, &bgk) < 1e-13,
                "omega = {omega}, seed {seed}"
            );
        }
    }
}

/// Raw moments of the equilibrium and of the post-collision state of `equilibrium + perturbation`,
/// where `perturbation` lists `(moment index, change of that raw moment)` pairs.
fn relax_moment_perturbation(
    coll: &RawMomentCollision,
    perturbation: &[(usize, f64)],
) -> ([f64; 19], [f64; 19]) {
    let (rho, u) = (1.0, [0.02, -0.01, 0.015]);
    let feq = compute_equilibrium_d3q19(rho, u);
    let meq = compute_raw_moments_d3q19(&feq);
    let mut m = meq;
    for &(k, delta) in perturbation {
        m[k] += delta;
    }
    let f = raw_moments_to_f_d3q19(&m);
    // The perturbation leaves density and momentum untouched.
    let (rho_f, u_f) = compute_macroscopic(&f);
    assert!((rho_f - rho).abs() < 1e-14);
    assert!(u_f.iter().zip(u.iter()).all(|(a, b)| (a - b).abs() < 1e-14));
    (meq, compute_raw_moments_d3q19(&coll.collide(&f, rho, u)))
}

const RATES: RawMomentCollision = RawMomentCollision {
    omega_v: 1.25,
    omega_b: 0.7,
    omega_g: 0.6,
};

#[test]
fn extensional_stress_relaxes_with_the_shear_rate() {
    // m_200 - m_020 is a trace-free normal-stress difference, i.e. shear: it must follow omega_v,
    // not omega_b (the old rate map used omega_b for m_200, m_020, m_002).
    let a = 1e-3;
    let (meq, out) = relax_moment_perturbation(&RATES, &[(4, a), (5, -a)]);
    let factor = 1.0 - RATES.omega_v;
    assert!(((out[4] - meq[4]) - factor * a).abs() < 1e-14, "m_200");
    assert!(((out[5] - meq[5]) + factor * a).abs() < 1e-14, "m_020");
    assert!(
        (out[6] - meq[6]).abs() < 1e-14,
        "m_002 must stay at equilibrium"
    );
    // The extensional decay rate differs clearly from the bulk one.
    assert!((factor - (1.0 - RATES.omega_b)).abs() > 0.5);
}

#[test]
fn trace_of_the_stress_relaxes_with_the_bulk_rate() {
    let a = 1e-3;
    let (meq, out) = relax_moment_perturbation(&RATES, &[(4, a), (5, a), (6, a)]);
    let factor = 1.0 - RATES.omega_b;
    for k in 4..7 {
        assert!(((out[k] - meq[k]) - factor * a).abs() < 1e-14, "m_{k}");
    }
}

#[test]
fn trace_free_part_and_trace_relax_independently_in_a_mixed_perturbation() {
    // dm = (3a, a, -2a): trace change 2a, deviator changes (2a - 2a/3, a - 2a/3, -2a - 2a/3).
    let a = 1e-3;
    let (meq, out) = relax_moment_perturbation(&RATES, &[(4, 3.0 * a), (5, a), (6, -2.0 * a)]);
    let trace_change = 2.0 * a;
    let deviators = [
        3.0 * a - trace_change / 3.0,
        a - trace_change / 3.0,
        -2.0 * a - trace_change / 3.0,
    ];
    for (k, deviator) in deviators.iter().enumerate() {
        let expected =
            (1.0 - RATES.omega_b) * trace_change / 3.0 + (1.0 - RATES.omega_v) * deviator;
        assert!(
            ((out[4 + k] - meq[4 + k]) - expected).abs() < 1e-14,
            "m_{}",
            4 + k
        );
    }
}

#[test]
fn off_diagonal_stress_relaxes_with_the_shear_rate() {
    let a = 1e-3;
    let (meq, out) = relax_moment_perturbation(&RATES, &[(7, a), (8, -2.0 * a), (9, 0.5 * a)]);
    let factor = 1.0 - RATES.omega_v;
    assert!(((out[7] - meq[7]) - factor * a).abs() < 1e-14);
    assert!(((out[8] - meq[8]) + 2.0 * factor * a).abs() < 1e-14);
    assert!(((out[9] - meq[9]) - 0.5 * factor * a).abs() < 1e-14);
}

#[test]
fn ghost_moments_relax_with_the_ghost_rate() {
    // Third- and fourth-order raw moments were discarded by the old reconstruction, so omega_g
    // had no effect at all.
    let a = 1e-3;
    let perturbation: Vec<(usize, f64)> = (10..19).map(|k| (k, a * (k as f64 - 9.0))).collect();
    let (meq, out) = relax_moment_perturbation(&RATES, &perturbation);
    let factor = 1.0 - RATES.omega_g;
    for &(k, delta) in &perturbation {
        assert!(
            ((out[k] - meq[k]) - factor * delta).abs() < 1e-14,
            "ghost moment {k}: {} vs {}",
            out[k] - meq[k],
            factor * delta
        );
    }
    // The stress moments are not touched by a pure ghost perturbation.
    for k in 4..10 {
        assert!((out[k] - meq[k]).abs() < 1e-14, "stress moment {k}");
    }
}

#[test]
fn viscosities_map_to_the_relaxation_rates() {
    let coll = RawMomentCollision::from_viscosity(0.1, 0.05);
    assert!((coll.shear_viscosity() - 0.1).abs() < 1e-14);
    assert!((coll.bulk_viscosity() - 0.05).abs() < 1e-14);
    assert!((coll.omega_v - 1.0 / 0.8).abs() < 1e-14);
    assert!((coll.omega_b - 1.0 / 0.725).abs() < 1e-14);
    assert_eq!(coll.omega_g, 1.0);
    // zeta = 0 puts the trace mode on the edge of stability, omega_b = 2.
    assert!((RawMomentCollision::from_viscosity(0.1, 0.0).omega_b - 2.0).abs() < 1e-14);
}

#[test]
fn shear_wave_decays_with_the_requested_viscosity() {
    let (nu, ny, steps, amplitude) = (0.02, 32, 400, 0.02);
    let coll = RawMomentCollision::from_viscosity(nu, 0.05);
    let ratio =
        shear_wave_amplitude_ratio(|f, rho, u| coll.collide(f, rho, u), ny, amplitude, steps);
    let expected = shear_wave_decay(nu, ny, steps);
    assert!(
        (ratio - expected).abs() < 0.01,
        "decay {ratio} vs exp(-nu k^2 t) = {expected}"
    );
}
