//! Regression tests for the D3Q27 MRT operator.
//!
//! The previous implementation ignored its moment matrix, relaxed every population with the
//! single stress rate `s[9]` (a BGK collision in disguise), discarded the bulk rate, and its
//! "inverse" was a weighted transpose of a rank-deficient matrix. These tests pin the exact
//! basis, the exact inverse, the per-moment relaxation and the conservation laws.

use super::*;

/// Deterministic pseudo-random populations around a moving equilibrium.
fn perturbed_cell(seed: u64, amplitude: f64) -> [f64; Q] {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let feq = D3Q27Lattice::equilibrium(1.1, 0.03, -0.02, 0.025);
    std::array::from_fn(|i| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let unit = (state >> 11) as f64 / (1u64 << 53) as f64;
        feq[i] * (1.0 + amplitude * (2.0 * unit - 1.0))
    })
}

/// Density and velocity of a cell.
fn macroscopic(cell: &[f64; Q]) -> (f64, f64, f64, f64) {
    D3Q27Lattice::new(1, 1, 1, 1.0).compute_macroscopic(cell)
}

/// Mass and momentum `[rho, jx, jy, jz]`.
fn mass_and_momentum(cell: &[f64; Q]) -> [f64; 4] {
    let mut m = [0.0_f64; 4];
    for (&fi, c) in cell.iter().zip(VELOCITIES.iter()) {
        m[0] += fi;
        m[1] += fi * f64::from(c[0]);
        m[2] += fi * f64::from(c[1]);
        m[3] += fi * f64::from(c[2]);
    }
    m
}

#[test]
fn moment_matrix_rows_are_orthogonal_under_the_lattice_weights() {
    let m = mrt_moment_matrix();
    for k in 0..Q {
        for l in 0..Q {
            let inner: f64 = (0..Q)
                .map(|j| WEIGHTS[j] * m[k * Q + j] * m[l * Q + j])
                .sum();
            if k == l {
                assert!(inner > 0.0, "row {k} has non-positive norm {inner}");
            } else {
                assert!(
                    inner.abs() < 1e-14,
                    "rows {k} and {l} are not orthogonal: {inner}"
                );
            }
        }
    }
}

#[test]
fn moment_matrix_inverse_is_exact() {
    let m = mrt_moment_matrix();
    let inv = mrt_inverse_moment_matrix();
    for i in 0..Q {
        for j in 0..Q {
            let left: f64 = (0..Q).map(|k| inv[i * Q + k] * m[k * Q + j]).sum();
            let right: f64 = (0..Q).map(|k| m[i * Q + k] * inv[k * Q + j]).sum();
            let identity = if i == j { 1.0 } else { 0.0 };
            assert!(
                (left - identity).abs() < 1e-14,
                "(M^-1 M)[{i}][{j}] = {left}"
            );
            assert!(
                (right - identity).abs() < 1e-14,
                "(M M^-1)[{i}][{j}] = {right}"
            );
        }
    }
}

/// Invert a 27x27 matrix with Gauss-Jordan elimination and partial pivoting (test oracle).
fn gauss_jordan_inverse(m: &[f64; Q * Q]) -> Vec<f64> {
    let mut a = vec![0.0_f64; Q * 2 * Q];
    for i in 0..Q {
        for j in 0..Q {
            a[i * 2 * Q + j] = m[i * Q + j];
        }
        a[i * 2 * Q + Q + i] = 1.0;
    }
    for col in 0..Q {
        let pivot_row = (col..Q)
            .max_by(|&r1, &r2| {
                a[r1 * 2 * Q + col]
                    .abs()
                    .total_cmp(&a[r2 * 2 * Q + col].abs())
            })
            .expect("non-empty pivot range");
        assert!(a[pivot_row * 2 * Q + col].abs() > 1e-12, "singular matrix");
        for j in 0..2 * Q {
            a.swap(col * 2 * Q + j, pivot_row * 2 * Q + j);
        }
        let pivot = a[col * 2 * Q + col];
        for j in 0..2 * Q {
            a[col * 2 * Q + j] /= pivot;
        }
        for row in 0..Q {
            if row != col {
                let factor = a[row * 2 * Q + col];
                for j in 0..2 * Q {
                    let delta = factor * a[col * 2 * Q + j];
                    a[row * 2 * Q + j] -= delta;
                }
            }
        }
    }
    let mut inv = vec![0.0_f64; Q * Q];
    for i in 0..Q {
        for j in 0..Q {
            inv[i * Q + j] = a[i * 2 * Q + Q + j];
        }
    }
    inv
}

#[test]
fn moment_matrix_inverse_matches_gauss_jordan() {
    let reference = gauss_jordan_inverse(&mrt_moment_matrix());
    let inv = mrt_inverse_moment_matrix();
    for (idx, (&a, &b)) in inv.iter().zip(reference.iter()).enumerate() {
        assert!((a - b).abs() < 1e-12, "entry {idx}: {a} vs {b}");
    }
}

#[test]
fn moment_transform_round_trips_exactly() {
    let m = mrt_moment_matrix();
    let inv = mrt_inverse_moment_matrix();
    for seed in 1..=6 {
        let f = perturbed_cell(seed, 0.4);
        let back = transform_from_moments(&inv, &transform_to_moments(&m, &f));
        for i in 0..Q {
            assert!(
                (back[i] - f[i]).abs() < 1e-14,
                "seed {seed}, population {i}"
            );
        }
    }
}

#[test]
fn hydrodynamic_moments_have_their_physical_meaning() {
    let m = mrt_moment_matrix();
    let f = perturbed_cell(3, 0.3);
    let moments = transform_to_moments(&m, &f);
    let sum = |weight: &dyn Fn(f64, f64, f64) -> f64| -> f64 {
        f.iter()
            .zip(VELOCITIES.iter())
            .map(|(&fi, c)| fi * weight(f64::from(c[0]), f64::from(c[1]), f64::from(c[2])))
            .sum()
    };
    let checks: [(usize, f64); 9] = [
        (0, sum(&|_, _, _| 1.0)),
        (1, sum(&|x, y, z| 3.0 * (x * x + y * y + z * z) - 3.0)),
        (3, sum(&|x, _, _| x)),
        (5, sum(&|_, y, _| y)),
        (7, sum(&|_, _, z| z)),
        (9, sum(&|x, y, _| 3.0 * (x * x - y * y))),
        (11, sum(&|x, y, _| x * y)),
        (12, sum(&|_, y, z| y * z)),
        (13, sum(&|x, _, z| x * z)),
    ];
    for (row, expected) in checks {
        assert!(
            (moments[row] - expected).abs() < 1e-13,
            "moment {row}: {} vs {expected}",
            moments[row]
        );
    }
}

#[test]
fn uniform_rates_reproduce_bgk() {
    for seed in 1..=6 {
        let f = perturbed_cell(seed, 0.2);
        let (rho, ux, uy, uz) = macroscopic(&f);
        let feq = D3Q27Lattice::equilibrium(rho, ux, uy, uz);
        for omega in [0.7, 1.0, 1.3, 1.85] {
            let mut mrt = f;
            collide_mrt_full(
                &mut mrt,
                rho,
                ux,
                uy,
                uz,
                &MrtRelaxationRates::uniform(omega),
            );
            for i in 0..Q {
                let bgk = f[i] - omega * (f[i] - feq[i]);
                assert!(
                    (mrt[i] - bgk).abs() < 1e-12,
                    "seed {seed}, omega {omega}, population {i}: {} vs BGK {bgk}",
                    mrt[i]
                );
            }
        }
    }
}

#[test]
fn non_uniform_rates_conserve_mass_and_momentum() {
    for seed in 1..=8 {
        let f = perturbed_cell(seed, 0.25);
        let (rho, ux, uy, uz) = macroscopic(&f);
        for (nu, xi) in [(0.02, 0.01), (1.0 / 6.0, 0.0), (0.3, 0.2)] {
            let mut out = f;
            let rates = MrtRelaxationRates::for_viscosity(nu, xi);
            collide_mrt_full(&mut out, rho, ux, uy, uz, &rates);
            let before = mass_and_momentum(&f);
            let after = mass_and_momentum(&out);
            for (c, name) in ["rho", "jx", "jy", "jz"].iter().enumerate() {
                assert!(
                    (before[c] - after[c]).abs() < 1e-13,
                    "seed {seed}, nu {nu}, xi {xi}: {name} {} -> {}",
                    before[c],
                    after[c]
                );
            }
        }
    }
}

#[test]
fn rates_of_conserved_moments_are_ignored() {
    let f = perturbed_cell(5, 0.2);
    let (rho, ux, uy, uz) = macroscopic(&f);
    let mut reference = f;
    collide_mrt_full(
        &mut reference,
        rho,
        ux,
        uy,
        uz,
        &MrtRelaxationRates::for_viscosity(0.1, 0.05),
    );
    let mut rates = MrtRelaxationRates::for_viscosity(0.1, 0.05);
    for mode in [0usize, 3, 5, 7] {
        rates.s[mode] = 1.7;
    }
    let mut out = f;
    collide_mrt_full(&mut out, rho, ux, uy, uz, &rates);
    assert_eq!(reference, out);
}

#[test]
fn every_moment_relaxes_with_its_own_rate() {
    let f = perturbed_cell(11, 0.3);
    let (rho, ux, uy, uz) = macroscopic(&f);
    let (nu, xi) = (0.08, 0.3);
    let rates = MrtRelaxationRates::for_viscosity(nu, xi);
    let mut out = f;
    collide_mrt_full(&mut out, rho, ux, uy, uz, &rates);

    let m = mrt_moment_matrix();
    let feq = D3Q27Lattice::equilibrium(rho, ux, uy, uz);
    let neq_before = transform_to_moments(&m, &std::array::from_fn::<f64, Q, _>(|i| f[i] - feq[i]));
    let neq_after =
        transform_to_moments(&m, &std::array::from_fn::<f64, Q, _>(|i| out[i] - feq[i]));
    for k in 0..Q {
        let s_k = if [0, 3, 5, 7].contains(&k) {
            0.0
        } else {
            rates.s[k]
        };
        assert!(
            (neq_after[k] - (1.0 - s_k) * neq_before[k]).abs() < 1e-13,
            "moment {k}: {} vs {}",
            neq_after[k],
            (1.0 - s_k) * neq_before[k]
        );
    }
    // The rates are genuinely different: shear, bulk (trace) and ghost modes.
    assert!((rates.s[9] - 1.0 / (0.5 + 3.0 * nu)).abs() < 1e-14);
    assert!((rates.s[1] - 1.0 / (0.5 + 4.5 * xi)).abs() < 1e-14);
    assert!((rates.s[1] - rates.s[9]).abs() > 1e-3);
    assert!((rates.s[17] - rates.s[9]).abs() > 1e-3);

    // Physical meaning: shear stress relaxes with s_nu, the trace with s_bulk.
    let before = non_equilibrium_stress_tensor(&f, rho, ux, uy, uz);
    let after = non_equilibrium_stress_tensor(&out, rho, ux, uy, uz);
    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
        assert!(
            (after[a][b] - (1.0 - rates.s[9]) * before[a][b]).abs() < 1e-13,
            "Pi[{a}][{b}]"
        );
    }
    let trace = |pi: &[[f64; 3]; 3]| pi[0][0] + pi[1][1] + pi[2][2];
    assert!(
        (trace(&after) - (1.0 - rates.s[1]) * trace(&before)).abs() < 1e-13,
        "trace of the non-equilibrium stress"
    );
}

#[test]
fn bulk_rate_is_the_trace_mode_rate() {
    // xi = 2 nu / 3 is the BGK-equivalent bulk viscosity: s_bulk == s_nu.
    let nu = 0.1;
    let rates = MrtRelaxationRates::for_viscosity(nu, 2.0 * nu / 3.0);
    assert!((rates.s[1] - rates.s[9]).abs() < 1e-14);
    // The default constructor no longer discards the bulk viscosity.
    let a = MrtRelaxationRates::for_viscosity(nu, 0.1);
    let b = MrtRelaxationRates::for_viscosity(nu, 0.4);
    assert!((a.s[1] - b.s[1]).abs() > 0.1);
}

#[test]
fn collide_all_mrt_full_conserves_mass_and_momentum_over_many_steps() {
    let rates = MrtRelaxationRates::for_viscosity(0.05, 0.05);
    let mut lat = D3Q27Lattice::new(4, 4, 4, 1.0);
    let mut state = 12345_u64;
    for cell in lat.f.iter_mut() {
        for fi in cell.iter_mut() {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let unit = (state >> 11) as f64 / (1u64 << 53) as f64;
            *fi *= 1.0 + 0.1 * (2.0 * unit - 1.0);
        }
    }
    let total = |lat: &D3Q27Lattice| {
        let mut m = [0.0_f64; 4];
        for cell in &lat.f {
            let cm = mass_and_momentum(cell);
            for c in 0..4 {
                m[c] += cm[c];
            }
        }
        m
    };
    let before = total(&lat);
    for _ in 0..30 {
        collide_all_mrt_full(&mut lat, &rates);
        lat.stream_periodic();
    }
    let after = total(&lat);
    for c in 0..4 {
        assert!(
            (before[c] - after[c]).abs() < 1e-10,
            "moment {c}: {} -> {}",
            before[c],
            after[c]
        );
    }
}

#[test]
fn mrt_shear_wave_decays_with_the_prescribed_viscosity() {
    // A sinusoidal shear wave u_x(y) = U sin(k y) must decay as exp(-nu k^2 t).
    let (nx, ny, nz) = (4, 16, 4);
    let nu = 0.05;
    let rates = MrtRelaxationRates::for_viscosity(nu, 2.0 * nu / 3.0);
    let mut lat = D3Q27Lattice::new(nx, ny, nz, 1.0);
    let u0 = 0.02;
    let k = std::f64::consts::TAU / ny as f64;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                lat.set_equilibrium(x, y, z, 1.0, u0 * (k * y as f64).sin(), 0.0, 0.0);
            }
        }
    }
    let amplitude = |lat: &D3Q27Lattice| -> f64 {
        let mut a = 0.0;
        for y in 0..ny {
            let idx = lat.idx(0, y, 0);
            let (_, ux, _, _) = lat.compute_macroscopic(&lat.f[idx]);
            a += 2.0 / ny as f64 * ux * (k * y as f64).sin();
        }
        a
    };
    let steps = 200;
    for _ in 0..steps {
        collide_all_mrt_full(&mut lat, &rates);
        lat.stream_periodic();
    }
    let expected = u0 * (-nu * k * k * steps as f64).exp();
    let measured = amplitude(&lat);
    assert!(
        (measured - expected).abs() < 0.03 * expected,
        "shear wave amplitude {measured} vs analytic {expected}"
    );
}
