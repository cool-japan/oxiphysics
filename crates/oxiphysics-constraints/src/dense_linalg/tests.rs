//! Tests of the dense linear-algebra kernels against reference values
//! computed with LAPACK (numpy) and against algebraic invariants.

use super::*;

fn assert_close(actual: f64, expected: f64, tol: f64, what: &str) {
    let scale = expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tol * scale,
        "{what}: got {actual}, expected {expected}"
    );
}

/// Check that `actual` and `expected` are the same multiset of eigenvalues
/// (order independent: eigenvalues with equal real parts may be sorted either way).
fn assert_eigenvalues_match(actual: &[(f64, f64)], expected: &[(f64, f64)], tol: f64) {
    assert_eq!(actual.len(), expected.len(), "eigenvalue count");
    let mut unused: Vec<(f64, f64)> = actual.to_vec();
    for want in expected {
        let nearest = unused
            .iter()
            .enumerate()
            .map(|(k, got)| (k, (got.0 - want.0).hypot(got.1 - want.1)))
            .min_by(|a, b| a.1.total_cmp(&b.1));
        let Some((k, dist)) = nearest else {
            panic!("no eigenvalue left to match {want:?}");
        };
        let scale = want.0.hypot(want.1).max(1.0);
        assert!(
            dist <= tol * scale,
            "expected eigenvalue {want:?} but the closest computed one is {:?} (distance {dist})",
            unused[k]
        );
        unused.swap_remove(k);
    }
}

fn mat_mul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            for j in 0..n {
                c[i * n + j] += a[i * n + k] * b[k * n + j];
            }
        }
    }
    c
}

// ── Cholesky ────────────────────────────────────────────────────────────

const SPD3: [f64; 9] = [4.0, 2.0, 0.6, 2.0, 5.0, 1.0, 0.6, 1.0, 3.0];

#[test]
fn cholesky_matches_lapack_factor() {
    let l = cholesky_factor(&SPD3, 3).expect("matrix is positive definite");
    let expected = [2.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.3, 0.35, 1.6695807857064];
    for (k, (got, want)) in l.iter().zip(expected.iter()).enumerate() {
        assert_close(*got, *want, 1e-13, &format!("L[{k}]"));
    }
}

#[test]
fn cholesky_rejects_indefinite_and_semidefinite_matrices() {
    // Eigenvalues 3 and -1.
    assert!(cholesky_factor(&[1.0, 2.0, 2.0, 1.0], 2).is_none());
    // Rank one: the second pivot is exactly zero.
    assert!(cholesky_factor(&[1.0, 1.0, 1.0, 1.0], 2).is_none());
    // Zero diagonal entry.
    assert!(cholesky_factor(&[0.0, 0.0, 0.0, 1.0], 2).is_none());
    // Non-finite entries.
    assert!(cholesky_factor(&[f64::NAN, 0.0, 0.0, 1.0], 2).is_none());
}

#[test]
fn cholesky_solve_reproduces_the_right_hand_side() {
    let l = cholesky_factor(&SPD3, 3).expect("matrix is positive definite");
    let b = [1.0, -2.0, 0.5, 0.0, 3.0, -1.0];
    let x = cholesky_solve(&l, 3, &b, 2);
    for c in 0..2 {
        for i in 0..3 {
            let ax: f64 = (0..3).map(|k| SPD3[i * 3 + k] * x[k * 2 + c]).sum();
            assert_close(ax, b[i * 2 + c], 1e-13, "A x = b");
        }
    }
}

// ── LU ──────────────────────────────────────────────────────────────────

#[test]
fn lu_solves_the_classic_three_by_three_system() {
    let a = [2.0, 1.0, -1.0, -3.0, -1.0, 2.0, -2.0, 1.0, 2.0];
    let b = [8.0, -11.0, -3.0];
    let x = lu_solve(&a, 3, &b, 1).expect("system is regular");
    for (got, want) in x.iter().zip([2.0, 3.0, -1.0]) {
        assert_close(*got, want, 1e-13, "solution");
    }
}

#[test]
fn lu_pivots_around_a_zero_leading_entry() {
    let x = lu_solve(&[0.0, 1.0, 1.0, 0.0], 2, &[3.0, 5.0], 1).expect("permutation matrix");
    assert_close(x[0], 5.0, 1e-14, "x0");
    assert_close(x[1], 3.0, 1e-14, "x1");
}

#[test]
fn lu_handles_several_right_hand_sides() {
    let a = [4.0, 1.0, 2.0, 3.0];
    // Identity right-hand side yields the inverse.
    let inv = lu_solve(&a, 2, &[1.0, 0.0, 0.0, 1.0], 2).expect("regular");
    let prod = mat_mul(&a, &inv, 2);
    for (k, want) in [1.0, 0.0, 0.0, 1.0].iter().enumerate() {
        assert_close(prod[k], *want, 1e-14, "A A^-1");
    }
}

#[test]
fn lu_rejects_singular_and_non_finite_matrices() {
    assert!(lu_solve(&[1.0, 2.0, 2.0, 4.0], 2, &[1.0, 1.0], 1).is_none());
    assert!(lu_solve(&[0.0, 0.0, 0.0, 0.0], 2, &[1.0, 1.0], 1).is_none());
    assert!(lu_solve(&[f64::INFINITY, 0.0, 0.0, 1.0], 2, &[1.0, 1.0], 1).is_none());
    assert!(lu_solve(&[f64::NAN, 0.0, 0.0, 1.0], 2, &[1.0, 1.0], 1).is_none());
}

// ── Singular values ─────────────────────────────────────────────────────

const TALL_6X3: [f64; 18] = [
    1.0, 2.0, 0.5, 0.3, -1.2, 2.2, 2.5, 0.7, -0.4, -1.1, 0.9, 1.6, 0.2, 0.4, -2.0, 1.7, -0.6, 0.8,
];

#[test]
fn singular_values_match_lapack() {
    let s = singular_values(&TALL_6X3, 6, 3);
    let expected = [3.6834643277816155, 3.3467446246933465, 2.5340463616608386];
    assert_eq!(s.len(), 3);
    for (k, (got, want)) in s.iter().zip(expected.iter()).enumerate() {
        assert_close(*got, *want, 1e-13, &format!("sigma_{k}"));
    }
}

#[test]
fn singular_values_of_a_wide_matrix_equal_those_of_its_transpose() {
    let mut wide = vec![0.0; 18];
    for r in 0..6 {
        for c in 0..3 {
            wide[c * 6 + r] = TALL_6X3[r * 3 + c];
        }
    }
    let s = singular_values(&wide, 3, 6);
    let expected = [3.6834643277816155, 3.3467446246933465, 2.5340463616608386];
    assert_eq!(s.len(), 3);
    for (got, want) in s.iter().zip(expected.iter()) {
        assert_close(*got, *want, 1e-13, "wide sigma");
    }
}

#[test]
fn singular_values_resolve_the_null_space_of_a_rank_deficient_matrix() {
    // Second row is twice the first: rank 2.  A Gram-matrix approach would
    // return ~1e-8 for the third value, the one-sided Jacobi method does not.
    let s = singular_values(&[1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 1.0, 0.0, 1.0], 3, 3);
    assert_close(s[0], 8.43544851578705, 1e-13, "sigma_0");
    assert_close(s[1], 0.9182637624920782, 1e-13, "sigma_1");
    assert!(s[2] < 1e-15, "null singular value {} should vanish", s[2]);
}

#[test]
fn singular_values_of_degenerate_shapes() {
    assert!(singular_values(&[], 0, 3).is_empty());
    assert!(singular_values(&[], 3, 0).is_empty());
    assert_eq!(singular_values(&[-3.0], 1, 1), vec![3.0]);
    assert_eq!(singular_values(&[0.0; 6], 2, 3), vec![0.0; 2]);
    // A diagonal matrix returns its absolute entries in descending order.
    let s = singular_values(&[1.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 2.0], 3, 3);
    assert_eq!(s, vec![4.0, 2.0, 1.0]);
}

#[test]
fn singular_value_product_equals_the_absolute_determinant() {
    let a = [2.0, -1.0, 0.5, 0.0, 3.0, 1.0, 1.0, 1.0, -2.0];
    // Laplace expansion along the first row: 2 M00 - (-1) M01 + 0.5 M02.
    let det: f64 =
        2.0 * (3.0 * -2.0 - 1.0 * 1.0) + (0.0 * -2.0 - 1.0 * 1.0) + 0.5 * (0.0 * 1.0 - 3.0 * 1.0);
    let prod: f64 = singular_values(&a, 3, 3).iter().product();
    assert_close(prod, det.abs(), 1e-13, "|det|");
}

// ── Eigenvalues ─────────────────────────────────────────────────────────

#[test]
fn eigenvalues_of_a_companion_matrix_are_the_polynomial_roots() {
    // (s+1)(s+2)(s+3)(s^2+2s+5): roots -1, -2, -3, -1 +- 2i.
    let c = [
        -8.0, -28.0, -58.0, -67.0, -30.0, //
        1.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, 0.0, //
        0.0, 0.0, 1.0, 0.0, 0.0, //
        0.0, 0.0, 0.0, 1.0, 0.0,
    ];
    let ev = real_eigenvalues(&c, 5).expect("converges");
    let expected = [
        (-3.0, 0.0),
        (-2.0, 0.0),
        (-1.0, -2.0),
        (-1.0, 0.0),
        (-1.0, 2.0),
    ];
    // Sorted by real part: -3, -2, then the three eigenvalues with real part -1.
    assert_eigenvalues_match(&ev, &expected, 1e-9);
}

#[test]
fn eigenvalues_match_lapack_for_a_dense_six_by_six() {
    let m = [
        2.041, -2.556, 0.418, -0.568, -0.453, -0.216, //
        -2.02, -0.232, -0.865, 3.323, 0.226, -0.353, //
        -0.281, -0.668, -1.055, -0.391, 0.482, -0.239, //
        0.958, -0.2, 0.024, 1.546, 0.545, -0.505, //
        -0.183, 0.541, 1.935, -0.27, -0.244, 1.002, //
        -0.886, -0.292, 0.883, 0.58, 0.092, 0.67,
    ];
    let expected = [
        (-1.8178848454730492, -0.17773452514573565),
        (-1.8178848454730492, 0.17773452514573565),
        (0.00833379288741885, 0.0),
        (1.0593474522864352, 0.0),
        (2.647044222886124, -1.2028396012234617),
        (2.647044222886124, 1.2028396012234617),
    ];
    let ev = real_eigenvalues(&m, 6).expect("converges");
    assert_eigenvalues_match(&ev, &expected, 1e-10);
}

#[test]
fn eigenvalues_match_lapack_for_a_dense_ten_by_ten() {
    let m = [
        -2.828, 1.021, -0.96, -1.669, 0.276, 0.701, -0.445, -1.076, 0.026, -0.053, //
        1.406, 0.747, 0.194, 1.112, -0.206, -0.926, 0.584, 0.583, -0.215, -0.783, //
        0.229, -2.494, 0.69, 0.491, -1.639, 0.061, -0.964, 0.757, -2.034, -0.914, //
        0.71, 1.156, -2.158, -0.498, 0.328, -0.609, 1.591, -1.191, 0.355, -1.048, //
        1.406, -0.022, -0.372, -1.718, 1.682, 0.753, 0.754, 1.138, 0.349, -0.639, //
        -0.8, -0.8, 1.37, -1.46, -0.596, -0.321, 0.225, 0.575, -1.249, -1.73, //
        -0.004, 1.214, 0.757, 0.216, -0.317, 0.293, -0.243, 0.817, -0.794, 0.134, //
        -0.111, 0.543, 0.225, 2.55, 1.499, 1.497, -2.04, -0.34, -0.609, 0.533, //
        -2.279, 1.174, 1.067, -1.302, -0.979, -0.801, 0.043, 0.641, 2.048, -0.197, //
        0.768, 0.155, 1.76, 0.742, 1.369, -1.078, -0.192, -0.814, 1.505, 0.658,
    ];
    let expected = [
        (-3.1202963538099784, 0.0),
        (-1.1050898058774803, -1.5686418703578116),
        (-1.1050898058774803, 1.5686418703578116),
        (-0.5826699295290272, 0.0),
        (0.11837035495251863, -2.9646705297750735),
        (0.11837035495251863, 2.9646705297750735),
        (0.7930193187681567, -0.8224046239573449),
        (0.7930193187681567, 0.8224046239573449),
        (2.8426832738263066, -0.1937836766679309),
        (2.8426832738263066, 0.1937836766679309),
    ];
    let ev = real_eigenvalues(&m, 10).expect("converges");
    assert_eigenvalues_match(&ev, &expected, 1e-10);
}

#[test]
fn balancing_recovers_the_eigenvalues_of_a_badly_scaled_matrix() {
    // Eigenvalues (3 +- sqrt(5)) / 2 although the off-diagonal entries differ
    // by twelve orders of magnitude.
    let ev = real_eigenvalues(&[1.0, 1e6, 1e-6, 2.0], 2).expect("converges");
    assert_close(ev[0].0, 0.3819660112501051, 1e-12, "small root");
    assert_close(ev[1].0, 2.618033988749895, 1e-12, "large root");
    assert_eq!(ev[0].1, 0.0);
    assert_eq!(ev[1].1, 0.0);
}

#[test]
fn eigenvalues_of_a_jordan_block_a_rotation_and_trivial_matrices() {
    let jordan = real_eigenvalues(&[0.0, 1.0, 0.0, 0.0], 2).expect("converges");
    assert_eq!(jordan, vec![(0.0, 0.0), (0.0, 0.0)]);
    let rotation = real_eigenvalues(&[0.0, -1.0, 1.0, 0.0], 2).expect("converges");
    assert_eigenvalues_match(&rotation, &[(0.0, -1.0), (0.0, 1.0)], 1e-14);
    assert_eq!(real_eigenvalues(&[-7.5], 1), Some(vec![(-7.5, 0.0)]));
    assert_eq!(real_eigenvalues(&[], 0), Some(Vec::new()));
    let zeros = real_eigenvalues(&[0.0; 9], 3).expect("converges");
    assert_eq!(zeros, vec![(0.0, 0.0); 3]);
    // A cyclic permutation matrix has the three cube roots of unity.
    let cyc =
        real_eigenvalues(&[0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0], 3).expect("converges");
    let h = 3.0_f64.sqrt() / 2.0;
    assert_eigenvalues_match(&cyc, &[(-0.5, -h), (-0.5, h), (1.0, 0.0)], 1e-12);
}

#[test]
fn eigenvalues_reject_malformed_input() {
    assert!(real_eigenvalues(&[1.0, f64::NAN, 0.0, 1.0], 2).is_none());
    assert!(real_eigenvalues(&[1.0, f64::INFINITY, 0.0, 1.0], 2).is_none());
    assert!(real_eigenvalues(&[1.0, 2.0, 3.0], 2).is_none());
}

#[test]
fn eigenvalues_preserve_trace_and_determinant() {
    let a = [
        0.5, -1.0, 2.0, 0.0, //
        1.5, 0.25, -0.5, 1.0, //
        -2.0, 0.75, 1.0, 0.3, //
        0.1, -0.4, 0.9, -1.5,
    ];
    let ev = real_eigenvalues(&a, 4).expect("converges");
    let trace: f64 = (0..4).map(|i| a[i * 4 + i]).sum();
    assert_close(ev.iter().map(|e| e.0).sum::<f64>(), trace, 1e-12, "trace");
    let imag_sum: f64 = ev.iter().map(|e| e.1).sum();
    assert!(imag_sum.abs() < 1e-12, "conjugate pairs cancel: {imag_sum}");
    // det = product of the eigenvalues (complex product, real result).
    let (mut re, mut im) = (1.0, 0.0);
    for (r, i) in &ev {
        let next_re = re * r - im * i;
        im = re * i + im * r;
        re = next_re;
    }
    let sv_product: f64 = singular_values(&a, 4, 4).iter().product();
    assert_close(
        (re * re + im * im).sqrt(),
        sv_product,
        1e-11,
        "|det| from eigenvalues vs singular values",
    );
}

#[test]
fn hamiltonian_spectrum_is_symmetric_about_the_imaginary_axis() {
    // H = [[A, B B'], [-C' C, -A']] for a lightly damped oscillator.
    let h = [
        0.0, 1.0, 0.0, 0.0, //
        -1.0, -0.1, 0.0, 1.0, //
        -1.0, 0.0, 0.0, 1.0, //
        0.0, 0.0, -1.0, 0.1,
    ];
    let ev = real_eigenvalues(&h, 4).expect("converges");
    for (re, im) in &ev {
        let mirrored = ev
            .iter()
            .any(|(r2, i2)| (r2 + re).abs() < 1e-9 && (i2 - im).abs() < 1e-9);
        assert!(
            mirrored,
            "eigenvalue {re}+{im}i has no mirror image in {ev:?}"
        );
    }
}
