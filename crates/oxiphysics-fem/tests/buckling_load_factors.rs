// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`BucklingProblem::solve_buckling_load_factors`].
//!
//! The old solver located higher modes by shifted inverse iteration with the
//! shift `previous·1.01`; inverse iteration converges to the eigenvalue nearest
//! the shift, which is the *first* mode again (`K = diag(1, 2, 3)`, `Kg = I`
//! returned `[1, 1, 1]`), and a singular solve returned the fabricated
//! `shift + 1`.  The load factors now come from the exact symmetric
//! eigen-decomposition of `L⁻¹·Kg·L⁻ᵀ` (`K = L·Lᵀ`).

use oxiphysics_fem::analysis::FemError;
use oxiphysics_fem::buckling::{
    BucklingProblem, cantilever_fem_buckling, column_fem_buckling, geometric_stiffness_beam,
};
use std::f64::consts::PI;

/// A problem with the given row-major `K` and `Kg`.
fn problem(k: &[f64], kg: &[f64], n: usize) -> BucklingProblem {
    let mut p = BucklingProblem::new(n);
    for i in 0..n {
        for j in 0..n {
            p.set_k(i, j, k[i * n + j]);
            p.set_kg(i, j, kg[i * n + j]);
        }
    }
    p
}

fn identity(n: usize) -> Vec<f64> {
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        m[i * n + i] = 1.0;
    }
    m
}

fn diagonal(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    let mut m = vec![0.0; n * n];
    for (i, &v) in values.iter().enumerate() {
        m[i * n + i] = v;
    }
    m
}

fn matvec(a: &[f64], x: &[f64]) -> Vec<f64> {
    let n = x.len();
    (0..n)
        .map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum())
        .collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[test]
fn higher_modes_are_not_copies_of_the_first_mode() {
    let p = problem(&diagonal(&[1.0, 2.0, 3.0]), &identity(3), 3);
    let factors = p.solve_buckling_load_factors(3).expect("K is SPD");
    assert_eq!(factors.len(), 3);
    for (got, want) in factors.iter().zip([1.0, 2.0, 3.0]) {
        assert!((got - want).abs() < 1e-12, "{factors:?}");
    }
}

#[test]
fn coupled_system_has_the_closed_form_spectrum() {
    // K = tridiag(−1, 2, −1), Kg = I: λ = 2 − √2, 2, 2 + √2.
    let k = [2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0];
    let p = problem(&k, &identity(3), 3);
    let factors = p.solve_buckling_load_factors(3).expect("K is SPD");
    let sqrt2 = std::f64::consts::SQRT_2;
    for (got, want) in factors.iter().zip([2.0 - sqrt2, 2.0, 2.0 + sqrt2]) {
        assert!((got - want).abs() < 1e-12, "{factors:?}");
    }
    // Requesting fewer modes returns the lowest ones.
    let two = p.solve_buckling_load_factors(2).expect("K is SPD");
    assert_eq!(two.len(), 2);
    assert!((two[0] - (2.0 - sqrt2)).abs() < 1e-12);
    assert!((two[1] - 2.0).abs() < 1e-12);
}

#[test]
fn physically_scaled_matrices_give_the_same_load_factors() {
    // K ~ 1e9, Kg ~ 1e3: the ratio, not the magnitude, decides the factors.
    let k = diagonal(&[1.0e9, 2.0e9, 4.0e9]);
    let kg = diagonal(&[1.0e3, 1.0e3, 2.0e3]);
    let factors = problem(&k, &kg, 3)
        .solve_buckling_load_factors(3)
        .expect("K is SPD");
    for (got, want) in factors.iter().zip([1.0e6, 2.0e6, 2.0e6]) {
        assert!((got - want).abs() <= 1e-9 * want, "{factors:?}");
    }
}

#[test]
fn indefinite_geometric_stiffness_reports_only_positive_load_factors() {
    // μ = 1/λ ∈ {1, −2, 1/2}: the negative one is a buckling load of the
    // *reversed* reference load and is not reported.
    let p = problem(&identity(3), &diagonal(&[1.0, -2.0, 0.5]), 3);
    let factors = p.solve_buckling_load_factors(3).expect("K is SPD");
    assert_eq!(factors.len(), 2, "{factors:?}");
    assert!((factors[0] - 1.0).abs() < 1e-12);
    assert!((factors[1] - 2.0).abs() < 1e-12);
}

#[test]
fn a_geometric_stiffness_with_a_null_space_gives_fewer_factors_not_spurious_ones() {
    let p = problem(&identity(3), &diagonal(&[1.0, 0.0, 0.0]), 3);
    let factors = p.solve_buckling_load_factors(3).expect("K is SPD");
    assert_eq!(factors, vec![1.0]);
}

#[test]
fn a_structure_in_tension_does_not_buckle() {
    let p = problem(&identity(2), &diagonal(&[-1.0, -3.0]), 2);
    assert!(
        p.solve_buckling_load_factors(2)
            .expect("K is SPD")
            .is_empty()
    );
    let zero_load = problem(&identity(2), &[0.0; 4], 2);
    assert!(
        zero_load
            .solve_buckling_load_factors(2)
            .expect("K is SPD")
            .is_empty()
    );
}

#[test]
fn zero_modes_and_empty_problems_are_empty_results() {
    let p = problem(&identity(2), &identity(2), 2);
    assert!(
        p.solve_buckling_load_factors(0)
            .expect("K is SPD")
            .is_empty()
    );
    let empty = BucklingProblem::new(0);
    assert!(
        empty
            .solve_buckling_load_factors(3)
            .expect("empty")
            .is_empty()
    );
}

#[test]
fn mode_shapes_solve_the_pencil_and_are_k_orthogonal() {
    let k = [4.0, -1.0, 0.0, -1.0, 3.0, -1.0, 0.0, -1.0, 5.0];
    let kg = [1.0, 0.2, 0.0, 0.2, 2.0, 0.1, 0.0, 0.1, 0.5];
    let p = problem(&k, &kg, 3);
    let modes = p.solve_buckling_modes(3).expect("K is SPD");
    assert_eq!(modes.len(), 3);
    assert!(
        modes
            .windows(2)
            .all(|w| w[0].load_factor <= w[1].load_factor)
    );
    for (i, mode_i) in modes.iter().enumerate() {
        let phi = &mode_i.mode_shape;
        let kphi = matvec(&k, phi);
        let kgphi = matvec(&kg, phi);
        let residual: f64 = kphi
            .iter()
            .zip(&kgphi)
            .map(|(a, b)| (a - mode_i.load_factor * b).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(residual < 1e-12 * dot(&kphi, &kphi).sqrt(), "mode {i}");
        assert!((dot(phi, phi) - 1.0).abs() < 1e-14, "unit norm, mode {i}");
        // Sign convention: the first component of (nearly) largest magnitude is positive.
        let largest = phi.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let reference = phi
            .iter()
            .copied()
            .find(|v| v.abs() >= (1.0 - 1.0e-9) * largest)
            .unwrap_or(0.0);
        assert!(reference > 0.0, "sign convention, mode {i}");
        for (j, mode_j) in modes.iter().enumerate() {
            if i != j {
                let k_inner = dot(&mode_j.mode_shape, &kphi);
                assert!(k_inner.abs() < 1e-12, "φ{j}ᵀKφ{i} = {k_inner}");
            }
        }
    }
}

#[test]
fn k_must_be_symmetric_positive_definite_and_inputs_symmetric() {
    // A rigid-body mode (singular K) is not a valid buckling problem.
    let singular = problem(&[1.0, -1.0, -1.0, 1.0], &identity(2), 2);
    assert!(matches!(
        singular.solve_buckling_load_factors(1),
        Err(FemError::InvalidInput(_))
    ));
    let indefinite = problem(&diagonal(&[1.0, -1.0]), &identity(2), 2);
    assert!(matches!(
        indefinite.solve_buckling_load_factors(1),
        Err(FemError::InvalidInput(_))
    ));
    let asymmetric_kg = problem(&identity(2), &[1.0, 0.5, 0.0, 1.0], 2);
    assert!(matches!(
        asymmetric_kg.solve_buckling_load_factors(1),
        Err(FemError::InvalidInput(_))
    ));
    let non_finite = problem(&[f64::NAN, 0.0, 0.0, 1.0], &identity(2), 2);
    assert!(matches!(
        non_finite.solve_buckling_load_factors(1),
        Err(FemError::InvalidInput(_))
    ));
}

/// Beam K and Kg of a column of `n_elements` Euler-Bernoulli elements with
/// DOFs `[w, θ]` per node, reduced by the supports in `fixed`.
fn beam_column(
    n_elements: usize,
    ei: f64,
    length: f64,
    load: f64,
    fixed: &[usize],
) -> BucklingProblem {
    let le = length / n_elements as f64;
    let total = 2 * (n_elements + 1);
    let mut k_full = vec![0.0; total * total];
    let mut kg_full = vec![0.0; total * total];
    let c = ei / (le * le * le);
    let ke = [
        [12.0 * c, 6.0 * c * le, -12.0 * c, 6.0 * c * le],
        [
            6.0 * c * le,
            4.0 * c * le * le,
            -6.0 * c * le,
            2.0 * c * le * le,
        ],
        [-12.0 * c, -6.0 * c * le, 12.0 * c, -6.0 * c * le],
        [
            6.0 * c * le,
            2.0 * c * le * le,
            -6.0 * c * le,
            4.0 * c * le * le,
        ],
    ];
    let kge = geometric_stiffness_beam(load, le);
    for element in 0..n_elements {
        let dofs = [
            2 * element,
            2 * element + 1,
            2 * element + 2,
            2 * element + 3,
        ];
        for a in 0..4 {
            for b in 0..4 {
                k_full[dofs[a] * total + dofs[b]] += ke[a][b];
                kg_full[dofs[a] * total + dofs[b]] += kge[a][b];
            }
        }
    }
    let free: Vec<usize> = (0..total).filter(|d| !fixed.contains(d)).collect();
    let mut p = BucklingProblem::new(free.len());
    for (ia, &ga) in free.iter().enumerate() {
        for (ib, &gb) in free.iter().enumerate() {
            p.set_k(ia, ib, k_full[ga * total + gb]);
            p.set_kg(ia, ib, kg_full[ga * total + gb]);
        }
    }
    p
}

/// A pinned-pinned column buckles at `P_k = k²·π²·EI/L²`.  The old solver
/// returned the first Euler load for every requested mode.
#[test]
fn pinned_pinned_column_modes_follow_the_euler_series() {
    let (ei, length, load) = (200.0e9 * 1.0e-6, 2.0, 1.0);
    let n_elements = 10;
    let p = beam_column(n_elements, ei, length, load, &[0, 2 * n_elements]);
    let factors = p.solve_buckling_load_factors(4).expect("K is SPD");
    assert_eq!(factors.len(), 4);
    let tolerances = [1.0e-4, 1.0e-3, 3.0e-3, 8.0e-3];
    for (index, (&factor, &tol)) in factors.iter().zip(&tolerances).enumerate() {
        let order = (index + 1) as f64;
        let euler = order * order * PI * PI * ei / (length * length);
        let relative = (factor * load - euler).abs() / euler;
        assert!(
            relative < tol,
            "mode {}: P = {}, Euler {euler}, relative error {relative:e}",
            index + 1,
            factor * load
        );
    }
    assert!(factors.windows(2).all(|w| w[0] < w[1]), "{factors:?}");
}

#[test]
fn cantilever_column_modes_follow_the_euler_series() {
    let (ei, length, load) = (200.0e9 * 1.0e-6, 2.0, 1.0);
    let p = beam_column(10, ei, length, load, &[0, 1]);
    let factors = p.solve_buckling_load_factors(3).expect("K is SPD");
    for (index, &factor) in factors.iter().enumerate() {
        let order = (2 * index + 1) as f64;
        let euler = (order * PI / 2.0).powi(2) * ei / (length * length);
        assert!(
            (factor * load - euler).abs() / euler < 2.0e-3,
            "mode {}: {}, Euler {euler}",
            index + 1,
            factor * load
        );
    }
}

#[test]
fn convenience_wrappers_agree_with_the_solver() {
    let (e, i_mom, length, load) = (200.0e9, 1.0e-6, 2.0, 1.0);
    let euler_pp = PI * PI * e * i_mom / (length * length);
    let pinned = column_fem_buckling(10, e, i_mom, length, load);
    assert!(
        (pinned * load - euler_pp).abs() / euler_pp < 1.0e-4,
        "{pinned}"
    );
    let euler_cf = PI * PI * e * i_mom / (4.0 * length * length);
    let cantilever = cantilever_fem_buckling(10, e, i_mom, length, load);
    assert!(
        (cantilever * load - euler_cf).abs() / euler_cf < 1.0e-5,
        "{cantilever}"
    );
    // Compression is positive: reversing the load leaves no positive factor.
    let reversed = column_fem_buckling(10, e, i_mom, length, -load);
    assert!(reversed.is_nan(), "{reversed}");
}
